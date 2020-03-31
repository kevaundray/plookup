use crate::kzg10;
use crate::multiset::MultiSet;
use crate::proof::{Commitments, Evaluations, MultiSetEqualityProof};
use crate::quotient_poly;
use crate::transcript::TranscriptProtocol;
use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use ff_fft::{DensePolynomial as Polynomial, EvaluationDomain};
use num_traits::identities::One;
use poly_commit::kzg10::Powers;
/// Computes the multisets h_1 and h_2
pub fn compute_h1_h2(f: &MultiSet, t: &MultiSet) -> (MultiSet, MultiSet) {
    //
    // 1. Compute s
    // XXX: check if s is (f,t) sorted by t? (Tests will fail anyways according to the proof, so may be better to skip)
    let s = f.concatenate(&t).sort();

    //2 . Compute h_1 and h_2
    let (h_1, h_2) = s.halve();
    // assert that the last element of h_1 is equal to the first element of h_2
    assert_eq!(h_1.0.last().unwrap(), &h_2.0[0]);
    (h_1, h_2)
}

// Computes the i+1'th term of F(beta, gamma)
fn compute_f_i(i: usize, f: &MultiSet, t: &MultiSet, beta: Fr, gamma: Fr) -> Fr {
    let gamma_beta_one = gamma * (beta + Fr::one());
    // (gamma + f_i) * [gamma *(1 + beta) + t_i + beta * t_{i+1}]
    (gamma + f.0[i]) * (gamma_beta_one + t.0[i] + (beta * t.0[i + 1]))
}

// Computes the i+1'th term of F(beta, gamma)
fn compute_g_i(i: usize, h_1: &MultiSet, h_2: &MultiSet, beta: Fr, gamma: Fr) -> Fr {
    let gamma_one_b = gamma * (Fr::one() + beta);

    // gamma * (1 + beta) + s_j + beta * s_{j+1}
    let d = gamma_one_b + h_1.0[i] + (beta * h_1.0[i + 1]);
    // gamma * (1 + beta) + s_{n+j} + beta * s_{n+j+1}
    let e = gamma_one_b + h_2.0[i] + (beta * h_2.0[i + 1]);

    d * e
}

/// Computes the values for Z(X)
pub fn compute_accumulator_values(
    f: &MultiSet,
    t: &MultiSet,
    h_1: &MultiSet,
    h_2: &MultiSet,
    beta: Fr,
    gamma: Fr,
) -> Vec<Fr> {
    let n = f.len();

    // F(beta, gamma)
    let mut numerator: Vec<Fr> = Vec::with_capacity(n + 1);
    // G(beta, gamma)
    let mut denominator: Vec<Fr> = Vec::with_capacity(n + 1);

    // Z evaluated at the first root of unity is 1
    numerator.push(Fr::one());
    denominator.push(Fr::one());

    let beta_one = Fr::one() + beta;

    // Compute values for Z(X)
    for i in 0..n {
        let f_i = beta_one * compute_f_i(i, f, t, beta, gamma);
        let g_i = compute_g_i(i, h_1, h_2, beta, gamma);

        let last_numerator = *numerator.last().unwrap();
        let last_denominator = *denominator.last().unwrap();

        numerator.push(f_i * last_numerator);
        denominator.push(g_i * last_denominator);
    }

    // Check that Z(g^{n+1}) = 1
    let last_numerator = *numerator.last().unwrap();
    let last_denominator = *denominator.last().unwrap();

    assert_eq!(last_numerator / last_denominator, Fr::one());

    // Combine numerator and denominator
    assert_eq!(numerator.len(), denominator.len());
    assert_eq!(numerator.len(), n + 1);
    let mut evaluations = Vec::with_capacity(numerator.len());
    for (n, d) in numerator.into_iter().zip(denominator) {
        evaluations.push(n / d)
    }
    evaluations
}

pub fn prove(
    f: MultiSet,
    t: MultiSet,
    proving_key: &Powers<Bls12_381>,
    transcript: &mut dyn TranscriptProtocol,
) -> MultiSetEqualityProof {
    let domain: EvaluationDomain<Fr> = EvaluationDomain::new(f.len()).unwrap();

    // Convert witness and table to polynomials
    let f_poly = f.to_polynomial(&domain);
    let f_commit = kzg10::commit(proving_key, &f_poly);

    let t_poly = t.to_polynomial(&domain);
    let t_commit = kzg10::commit(proving_key, &t_poly);

    // Compute h_1 and h_2
    let (h_1, h_2) = compute_h1_h2(&f, &t);

    // Convert h_1 and h_2 to polynomials
    let h_1_poly = h_1.to_polynomial(&domain);
    let h_2_poly = h_2.to_polynomial(&domain);

    // Commit to h_1(X) and h_2(X)
    let h_1_commit = kzg10::commit(proving_key, &h_1_poly);
    let h_2_commit = kzg10::commit(proving_key, &h_2_poly);

    // Add commitments to transcript
    transcript.append_commitment(b"h_1_poly", &h_1_commit);
    transcript.append_commitment(b"h_2_poly", &h_2_commit);

    let beta = transcript.challenge_scalar(b"beta");
    let gamma = transcript.challenge_scalar(b"gamma");

    // Compute Z(X)
    let z_evaluations = compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
    let z_poly = Polynomial::from_coefficients_vec(domain.ifft(&z_evaluations));

    // Commit to Z(X)
    let z_commit = kzg10::commit(proving_key, &z_poly);
    transcript.append_commitment(b"accumulator_poly", &z_commit);

    // Compute quotient polynomial
    let (quotient_poly, _) = quotient_poly::compute(
        &domain, &z_poly, &f_poly, &t_poly, &h_1_poly, &h_2_poly, beta, gamma,
    );

    // Commit to quotient polynomial
    let q_commit = kzg10::commit(proving_key, &quotient_poly);
    transcript.append_commitment(b"quotient_poly", &q_commit);

    // Compute the Witness that f was a subset of t
    //
    let evaluation_challenge = transcript.challenge_scalar(b"evaluation_challenge");
    transcript.append_scalar(b"evaluation_challenge", &evaluation_challenge);
    let evaluation_omega = evaluation_challenge * domain.group_gen;

    // Compute evaluations at `z`
    let f_eval = f_poly.evaluate(evaluation_challenge);
    let t_eval = t_poly.evaluate(evaluation_challenge);
    let h_1_eval = h_1_poly.evaluate(evaluation_challenge);
    let h_2_eval = h_2_poly.evaluate(evaluation_challenge);
    let z_eval = z_poly.evaluate(evaluation_challenge);
    let q_eval = quotient_poly.evaluate(evaluation_challenge);

    // Compute evaluations at `z * omega`
    let t_omega_eval = t_poly.evaluate(evaluation_omega);
    let h_1_omega_eval = h_1_poly.evaluate(evaluation_omega);
    let h_2_omega_eval = h_2_poly.evaluate(evaluation_omega);
    let z_omega_eval = z_poly.evaluate(evaluation_omega);

    transcript.append_scalar(b"f_eval", &f_eval);
    transcript.append_scalar(b"t_eval", &t_eval);
    transcript.append_scalar(b"h_1_eval", &h_1_eval);
    transcript.append_scalar(b"h_2_eval", &h_2_eval);
    transcript.append_scalar(b"z_eval", &z_eval);
    transcript.append_scalar(b"q_eval", &q_eval);
    transcript.append_scalar(b"t_omega_eval", &t_omega_eval);
    transcript.append_scalar(b"h_1_omega_eval", &h_1_omega_eval);
    transcript.append_scalar(b"h_2_omega_eval", &h_2_omega_eval);
    transcript.append_scalar(b"z_omega_eval", &z_omega_eval);

    let aggregation_challenge = transcript.challenge_scalar(b"witness_aggregation");

    // Compute opening proof for f(X) evaluated at `z`
    let agg_witness = kzg10::compute_aggregate_witness(
        vec![
            &f_poly,
            &t_poly,
            &h_1_poly,
            &h_2_poly,
            &z_poly,
            &quotient_poly,
        ],
        evaluation_challenge,
        aggregation_challenge,
    );
    let agg_witness_comm = kzg10::commit(proving_key, &agg_witness);

    // Compute opening proofs for f(X) evaluated at `z * omega`
    let shifted_agg_witness = kzg10::compute_aggregate_witness(
        vec![&t_poly, &h_1_poly, &h_2_poly, &z_poly],
        evaluation_omega,
        aggregation_challenge,
    );
    let shifted_agg_witness_comm = kzg10::commit(proving_key, &shifted_agg_witness);

    MultiSetEqualityProof {
        n: domain.size(),
        evaluations: Evaluations {
            f: f_eval,
            t: t_eval,
            t_omega: t_omega_eval,
            h_1: h_1_eval,
            h_1_omega: h_1_omega_eval,
            h_2: h_2_eval,
            h_2_omega: h_2_omega_eval,
            z: z_eval,
            z_omega: z_omega_eval,
        },
        commitments: Commitments {
            f: f_commit,
            q: q_commit,
            t: t_commit,
            h_1: h_1_commit,
            h_2: h_2_commit,
            z: z_commit,
        },
        aggregate_witness_comm: agg_witness_comm,
        shifted_aggregate_witness_comm: shifted_agg_witness_comm,
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ff_fft::{DensePolynomial as Polynomial, EvaluationDomain};

    #[test]
    fn test_manually_compute_z() {
        // This test manually computes the values of the accumulator Z(x)
        // Notice that the test will fail if:
        // - You add a value to 'f' that is not in 't'
        // - 't' is unordered
        // Now notice that the test will pass if:
        // - (1) len(t) != len(f) + 1
        // - (2) len(t) is not a power of two
        // The reason why the tests pass for these values is :
        // (1) We made this restriction, so that it would be easier to check that h_1 and h_2 are continuous. This should not affect the outcome of Z(X) if h_1 and h_2 when merged
        // indeed form 's'
        // (2) This is a restriction that is placed upon the protocol due to using roots of unity. It would not affect the computation of Z(X)

        let mut f = MultiSet::new();
        f.push(Fr::from(1u8));
        f.push(Fr::from(1u8));
        f.push(Fr::from(2u8));
        f.push(Fr::from(2u8));
        f.push(Fr::from(3u8));
        f.push(Fr::from(3u8));
        f.push(Fr::from(3u8));

        // Table of values
        let mut t = MultiSet::new();
        t.push(Fr::from(1u8));
        t.push(Fr::from(1u8));
        t.push(Fr::from(2u8));
        t.push(Fr::from(2u8));
        t.push(Fr::from(3u8));
        t.push(Fr::from(3u8));
        t.push(Fr::from(4u8));
        t.push(Fr::from(5u8));

        let beta = Fr::from(8u8);
        let gamma = Fr::from(10u8);

        let (h_1, h_2) = compute_h1_h2(&f, &t);

        let beta_one = Fr::one() + beta;

        // Manually compute the accumulator values
        //
        // First value of z_0 is 1 / 1
        let z_0_numerator = Fr::one();
        let z_0_denominator = Fr::one();
        let z_0 = z_0_numerator / z_0_denominator;
        //
        // Next value z_1 is (1+beta) * (z_0_numerator * f_0) / (z_0_denominator * g_0)
        let f_0 = compute_f_i(0, &f, &t, beta, gamma);
        let g_0 = compute_g_i(0, &h_1, &h_2, beta, gamma);
        let z_1_numerator = beta_one * z_0_numerator * f_0;
        let z_1_denominator = z_0_denominator * g_0;
        let z_1 = z_1_numerator / z_1_denominator;
        //
        // Next value z_2 is (1+beta)^2 * (z_1_numerator * f_1) / (z_1_denominator * g_1)
        let f_1 = compute_f_i(1, &f, &t, beta, gamma);
        let g_1 = compute_g_i(1, &h_1, &h_2, beta, gamma);
        let z_2_numerator = beta_one * z_1_numerator * f_1;
        let z_2_denominator = z_1_denominator * g_1;
        let z_2 = z_2_numerator / z_2_denominator;
        //
        // Next value z_3 is (1+beta)^3 * (z_2_numerator * f_2) / (z_2_denominator * g_2)
        let f_2 = compute_f_i(2, &f, &t, beta, gamma);
        let g_2 = compute_g_i(2, &h_1, &h_2, beta, gamma);
        let z_3_numerator = beta_one * z_2_numerator * f_2;
        let z_3_denominator = z_2_denominator * g_2;
        let z_3 = z_3_numerator / z_3_denominator;
        //
        // Next value z_4 is (1+beta)^4 * (z_3_numerator * f_3) / (z_3_denominator * g_3)
        let f_3 = compute_f_i(3, &f, &t, beta, gamma);
        let g_3 = compute_g_i(3, &h_1, &h_2, beta, gamma);
        let z_4_numerator = beta_one * z_3_numerator * f_3;
        let z_4_denominator = z_3_denominator * g_3;
        let z_4 = z_4_numerator / z_4_denominator;
        //
        // Next value z_5 is (1+beta)^5 * (z_4_numerator * f_4) / (z_4_denominator * g_4)
        let f_4 = compute_f_i(4, &f, &t, beta, gamma);
        let g_4 = compute_g_i(4, &h_1, &h_2, beta, gamma);
        let z_5_numerator = beta_one * z_4_numerator * f_4;
        let z_5_denominator = z_4_denominator * g_4;
        let z_5 = z_5_numerator / z_5_denominator;
        //
        // Next value z_6 is (1+beta)^6 * (z_5_numerator * f_5) / (z_5_denominator * g_5)
        let f_5 = compute_f_i(5, &f, &t, beta, gamma);
        let g_5 = compute_g_i(5, &h_1, &h_2, beta, gamma);
        let z_6_numerator = beta_one * z_5_numerator * f_5;
        let z_6_denominator = z_5_denominator * g_5;
        let z_6 = z_6_numerator / z_6_denominator;
        //
        // Last value z_7 is (1+beta)^7 * (z_6_numerator * f_6) / (z_6_denominator * g_6)
        // For an honest prover, this should be 1
        let f_6 = compute_f_i(6, &f, &t, beta, gamma);
        let g_6 = compute_g_i(6, &h_1, &h_2, beta, gamma);
        let z_7_numerator = beta_one * z_6_numerator * f_6;
        let z_7_denominator = z_6_denominator * g_6;
        let z_7 = z_7_numerator / z_7_denominator;

        // Check that the next value in z can be computed using the previously accumulated values multiplied by the next term
        // ie z_{n+1} = (1+beta) * z_n * (f_n / g_n)
        // Except for the last element which should be equal to 1
        assert_eq!(z_1, beta_one * z_0 * (f_0 / g_0));
        assert_eq!(z_2, beta_one * z_1 * (f_1 / g_1));
        assert_eq!(z_3, beta_one * z_2 * (f_2 / g_2));
        assert_eq!(z_4, beta_one * z_3 * (f_3 / g_3));
        assert_eq!(z_5, beta_one * z_4 * (f_4 / g_4));
        assert_eq!(z_6, beta_one * z_5 * (f_5 / g_5));
        assert_eq!(z_7, Fr::one());

        // Now check if we get the same values when computed by our function
        let expected_z_evaluations = vec![z_0, z_1, z_2, z_3, z_4, z_5, z_6, z_7];
        let z_evaluations = compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
        assert_eq!(expected_z_evaluations.len(), z_evaluations.len());
        for (should_be, got) in expected_z_evaluations.iter().zip(z_evaluations.iter()) {
            assert_eq!(should_be, got)
        }
    }
    #[test]
    fn test_h1_h2() {
        // Checks whether h_1 and h_2 are well formed(continuous) in s
        // This is done by checking that the last element of h_1 is the first element of h_2

        let (f, t, _) = setup_correct_test();

        // Compute h_1(x) and h_2(x) from f and t
        let (h_1, h_2) = compute_h1_h2(&f, &t);
        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(h_1.len()).unwrap();
        let h_1_poly = h_1.to_polynomial(&domain);
        let h_2_poly = h_2.to_polynomial(&domain);

        // compute the last and first element in the domain
        let last_element = domain.elements().last().unwrap();
        let first_element = Fr::one();

        let h_1_last = h_1_poly.evaluate(last_element);
        let h_2_first = h_2_poly.evaluate(first_element);

        assert_eq!(h_1_last, h_2_first);
    }

    #[test]
    fn test_term_check() {
        // Checks whether z_{n+1} = z_n * (f_n / g_n) holds for every value in the domain except for the last element
        // Then checks that Z(X) evaluated at the last element in the domain is 1

        let (f, t, domain) = setup_correct_test();
        let f_poly = f.to_polynomial(&domain);
        let t_poly = t.to_polynomial(&domain);

        let beta = Fr::from(5u8);
        let gamma = Fr::from(6u8);

        let (h_1, h_2) = compute_h1_h2(&f, &t);
        let h_1_poly = h_1.to_polynomial(&domain);
        let h_2_poly = h_2.to_polynomial(&domain);

        let z_evaluations = compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
        let z_poly = Polynomial::from_coefficients_vec(domain.ifft(&z_evaluations));

        let beta_one = Fr::one() + beta;
        let last_element = domain.elements().last().unwrap();

        for (_, element) in domain.elements().enumerate() {
            // Evaluate polynomials
            // evaluate z(X)
            let z_x = z_poly.evaluate(element);
            // evaluate z(Xg)
            let z_x_next = z_poly.evaluate(element * domain.group_gen);
            // evaluate f(X)
            let f_x = f_poly.evaluate(element);
            // evaluate t(X)
            let t_x = t_poly.evaluate(element);
            // evaluate t(Xg)
            let t_x_next = t_poly.evaluate(element * domain.group_gen);
            // evaluate h_1(X)
            let h_1_x = h_1_poly.evaluate(element);
            // evaluate h_1(Xg)
            let h_1_x_next = h_1_poly.evaluate(element * domain.group_gen);
            // evaluate h_2(X)
            let h_2_x = h_2_poly.evaluate(element);
            // evaluate h_2(Xg)
            let h_2_x_next = h_2_poly.evaluate(element * domain.group_gen);

            // LHS = (x - g^n)[z(x) * (1+b) * (gamma + f(x))] ( gamma(1 + beta) + t(x) + beta * t(Xg))

            // x - g^n
            let a = element - last_element;

            // z(x) * (1+b) * (gamma + f(x))
            let b = z_x * beta_one * (gamma + f_x);

            // gamma(1 + beta) + t(x) + beta * t(Xg)
            let c = (gamma * beta_one) + t_x + (beta * t_x_next);

            let lhs = a * b * c;

            // RHS = z(Xg)(x - g^n) [gamma *(1 + beta) + h_1(X) + beta*h_1(Xg)] [gamma * (1 + beta)] + h_2(X) + (beta * h_2(Xg))

            // z(Xg)(x - g^n)
            let a = z_x_next * (element - last_element);

            // [gamma *(1 + beta) + h_1(X) + beta*h_1(Xg)]
            let b = (gamma * beta_one) + h_1_x + (beta * h_1_x_next);

            // [gamma * (1 + beta)] + h_2(X) + (beta * h_2(Xg))
            let c = (gamma * beta_one) + h_2_x + (beta * h_2_x_next);

            let rhs = a * b * c;

            assert_eq!(lhs, rhs,);
        }

        // Now check that the last element is equal to 1
        assert_eq!(z_poly.evaluate(last_element), Fr::one())
    }

    // This is just a helper function to setup tests with values that work
    fn setup_correct_test() -> (MultiSet, MultiSet, EvaluationDomain<Fr>) {
        // n is the amount of witness values
        // d is the amount of table values
        // We need d = n+1

        let mut f = MultiSet::new();
        f.push(Fr::from(1u8));
        f.push(Fr::from(2u8));
        f.push(Fr::from(3u8));
        f.push(Fr::from(4u8));
        f.push(Fr::from(5u8));
        f.push(Fr::from(6u8));
        f.push(Fr::from(7u8));

        // Table of values
        let mut t = MultiSet::new();
        t.push(Fr::from(1u8));
        t.push(Fr::from(2u8));
        t.push(Fr::from(3u8));
        t.push(Fr::from(4u8));
        t.push(Fr::from(5u8));
        t.push(Fr::from(6u8));
        t.push(Fr::from(7u8));
        t.push(Fr::from(7u8));

        assert_eq!(t.len(), f.len() + 1);
        assert_eq!(t.len().next_power_of_two(), t.len());

        // The domain will be n+1
        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(f.len() + 1).unwrap();

        (f, t, domain)
    }
}
