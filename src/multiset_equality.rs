use crate::multiset::MultiSet;
use algebra::bls12_381::Fr;
use ff_fft::EvaluationDomain;
use num_traits::identities::One;

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

// Computes the i'th term of F(beta, gamma)
fn compute_f_i(i: usize, f: &MultiSet, t: &MultiSet, beta: Fr, gamma: Fr) -> Fr {
    let gamma_beta_one = gamma * (beta + Fr::one());
    // (gamma + f_i) * [gamma *(1 + beta) + t_i + beta * t_{i+1}]
    (gamma + f.0[i]) * (gamma_beta_one + t.0[i] + (beta * t.0[i + 1]))
}

// Computes the i'th term of F(beta, gamma)
fn compute_g_i(i: usize, h_1: &MultiSet, h_2: &MultiSet, beta: Fr, gamma: Fr) -> Fr {
    let gamma_one_b = gamma * (Fr::one() + beta);

    // gamma * (1 + beta) + s_j + beta * s_{j+1}
    let d = gamma_one_b + h_1.0[i] + (beta * h_1.0[i + 1]);
    // gamma * (1 + beta) + s_{n+j} + beta * s_{n+j+1}
    let e = gamma_one_b + h_2.0[i] + (beta * h_2.0[i + 1]);

    d * e
}

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

    // (1+\beta)
    let beta_one = Fr::one() + beta;

    // First value of z_0 is 1 / 1
    let z_0_numerator = Fr::one();
    let z_0_denominator = Fr::one();
    let z_0 = z_0_numerator / z_0_denominator;
    // Next value z_1 is (1+beta)^n * (z_0_numerator * f_0) / (z_0_denominator * g_0)
    let f_0 = compute_f_i(0, &f, &t, beta, gamma);
    let g_0 = compute_g_i(0, &h_1, &h_2, beta, gamma);
    let z_1_numerator = beta_one * z_0_numerator * f_0;
    let z_1_denominator = z_0_denominator * g_0;
    let z_1 = z_1_numerator / z_1_denominator;
    // Next value z_2 is (1+beta)^n * (z_1_numerator * f_1) / (z_1_denominator * g_1)
    let f_1 = compute_f_i(1, &f, &t, beta, gamma);
    let g_1 = compute_g_i(1, &h_1, &h_2, beta, gamma);
    let z_2_numerator = beta_one * z_1_numerator * f_1;
    let z_2_denominator = z_1_denominator * g_1;
    let z_2 = z_2_numerator / z_2_denominator;
    // Next value z_3 is (2+beta)^n * (z_2_numerator * f_2) / (z_2_denominator * g_2)
    let f_2 = compute_f_i(2, &f, &t, beta, gamma);
    let g_2 = compute_g_i(2, &h_1, &h_2, beta, gamma);
    let z_3_numerator = beta_one * z_2_numerator * f_2;
    let z_3_denominator = z_2_denominator * g_2;
    let z_3 = z_3_numerator / z_3_denominator;
    // Next value z_4 is (1+beta)^n * (z_3_numerator * f_3) / (z_3_denominator * g_3)
    let f_3 = compute_f_i(3, &f, &t, beta, gamma);
    let g_3 = compute_g_i(3, &h_1, &h_2, beta, gamma);
    let z_4_numerator = beta_one * z_3_numerator * f_3;
    let z_4_denominator = z_3_denominator * g_3;
    let z_4 = z_4_numerator / z_4_denominator;
    // Next value z_5 is (1+beta)^n * (z_4_numerator * f_4) / (z_4_denominator * g_4)
    let f_4 = compute_f_i(4, &f, &t, beta, gamma);
    let g_4 = compute_g_i(4, &h_1, &h_2, beta, gamma);
    let z_5_numerator = beta_one * z_4_numerator * f_4;
    let z_5_denominator = z_4_denominator * g_4;
    let z_5 = z_5_numerator / z_5_denominator;
    // Next value z_6 is (1+beta)^n * (z_5_numerator * f_5) / (z_5_denominator * g_5)
    let f_5 = compute_f_i(5, &f, &t, beta, gamma);
    let g_5 = compute_g_i(5, &h_1, &h_2, beta, gamma);
    let z_6_numerator = beta_one * z_5_numerator * f_5;
    let z_6_denominator = z_5_denominator * g_5;
    let z_6 = z_6_numerator / z_6_denominator;
    // Last value z_7 is (1+beta)^n * (z_6_numerator * f_6) / (z_6_denominator * g_6)
    // For an honest prover, this should be 1
    let f_6 = compute_f_i(6, &f, &t, beta, gamma);
    let g_6 = compute_g_i(6, &h_1, &h_2, beta, gamma);
    let z_7_numerator = beta_one * z_6_numerator * f_6;
    let z_7_denominator = z_6_denominator * g_6;
    let z_7 = z_7_numerator / z_7_denominator;

    // Check that the next value in z can be computed using the previously accumulated values multiplied by the next term
    // ie z_n = z_n-1 * (f_n / g_n)
    // Except for the last element which should be equal to 1
    assert_eq!(z_1, beta_one * z_0 * (f_0 / g_0));
    assert_eq!(z_2, beta_one * z_1 * (f_1 / g_1));
    assert_eq!(z_3, beta_one * z_2 * (f_2 / g_2));
    assert_eq!(z_4, beta_one * z_3 * (f_3 / g_3));
    assert_eq!(z_5, beta_one * z_4 * (f_4 / g_4));
    assert_eq!(z_6, beta_one * z_5 * (f_5 / g_5));
    assert_eq!(z_7, Fr::one());
}
#[test]
fn test_h1_h2() {
    // Checks whether h_1 and h_2 are well formed(continuous) in s
    // This is done by checking that the last element of h_1 is the first element of h_2

    let (f, t, _) = setup_test();

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

// This is just a helper function to setup tests with values that work
fn setup_test() -> (MultiSet, MultiSet, EvaluationDomain<Fr>) {
    // n is the amount of witness values
    // d is the amount of table values
    // We need d = n+1

    let mut f = MultiSet::new();
    f.push(Fr::from(1u8));
    f.push(Fr::from(1u8));
    f.push(Fr::from(1u8));
    f.push(Fr::from(1u8));
    f.push(Fr::from(1u8));
    f.push(Fr::from(1u8));
    f.push(Fr::from(1u8));

    // Table of values
    let mut t = MultiSet::new();
    t.push(Fr::from(1u8));
    t.push(Fr::from(2u8));
    t.push(Fr::from(1u8));
    t.push(Fr::from(1u8));
    t.push(Fr::from(1u8));
    t.push(Fr::from(1u8));
    t.push(Fr::from(1u8));
    t.push(Fr::from(1u8));

    assert_eq!(t.len(), f.len() + 1);
    assert_eq!(t.len().next_power_of_two(), t.len());

    // The domain will be n+1
    let domain: EvaluationDomain<Fr> = EvaluationDomain::new(f.len() + 1).unwrap();

    (f, t, domain)
}
