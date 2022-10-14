use crate::{
    kzg10,
    multiset::{multiset_equality, quotient_poly, MultiSet},
    transcript::TranscriptProtocol,
};
use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::{
    polynomial::univariate::DensePolynomial as Polynomial, EvaluationDomain, Polynomial as Poly,
    Radix2EvaluationDomain, UVPolynomial,
};
use ark_poly_commit::kzg10::{Commitment, Powers, VerifierKey};
// Evaluations store the evaluations of different polynomial.
// `t` denotes that the polynomial was evaluated at t(z) for some random evaluation challenge `z`
// `t_omega` denotes the polynomial was evaluated at t(z * omega) where omega is the group generator
// In the FFT context, the normal terminology is that t(z*omega) means to evaluate a polynomial at the next root of unity from `z`.
pub struct Evaluations {
    pub f: Fr,
    pub t: Fr,
    pub t_omega: Fr,
    pub h_1: Fr,
    pub h_1_omega: Fr,
    pub h_2: Fr,
    pub h_2_omega: Fr,
    pub z: Fr,
    pub z_omega: Fr,
}
// Commitments of different polynomials
pub struct Commitments {
    pub f: Commitment<Bls12_381>,
    pub q: Commitment<Bls12_381>,
    pub h_1: Commitment<Bls12_381>,
    pub h_2: Commitment<Bls12_381>,
    pub z: Commitment<Bls12_381>,
}

// In the best case, this protocol requires 4 extra G1 elements (Commitment)
// These are: h_1_commit,h_2_commit, f_commit,t_commit
//
// The commitment to the accumulator and the quotient polynomial would ideally be joined into the existing ones in PLONK
//
// We would require 7 extra Scalar elements (Evaluations)
// These are: h_1_eval, h_1_omega_eval, h_2_eval, h_2_omega_eval, f_eval, t_eval, t_omega_eval
//
// We would ideally be able to combine the accumulator, Z(X) with the permutation accumulator in plonk, and the quotient polynomial with the quotient polynomial in PLONK
// Which would save us 2 evaluation points: z_eval and z_omega_eval
// q_eval which is the quotient evaluation is usually created from the prover messages
//
// Lastly, the Witness commitments can also be batched with the PLONK opening Proof.
pub struct EqualityProof {
    pub aggregate_witness_comm: Commitment<Bls12_381>,
    pub shifted_aggregate_witness_comm: Commitment<Bls12_381>,

    pub evaluations: Evaluations,

    pub commitments: Commitments,
}

impl EqualityProof {
    pub fn prove(
        f: MultiSet,
        t: MultiSet,
        proving_key: &Powers<Bls12_381>,
        transcript: &mut dyn TranscriptProtocol,
    ) -> EqualityProof {
        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(t.len()).unwrap();
        // Convert witness and table to polynomials
        let f_poly = f.to_polynomial(&domain);
        let f_commit = kzg10::commit(proving_key, &f_poly);
        let t_poly = t.to_polynomial(&domain);
        // Compute h_1 and h_2
        let (h_1, h_2) = multiset_equality::compute_h1_h2(&f, &t);
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
        let z_evaluations =
            multiset_equality::compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
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
        let f_eval = f_poly.evaluate(&evaluation_challenge);
        let t_eval = t_poly.evaluate(&evaluation_challenge);
        let h_1_eval = h_1_poly.evaluate(&evaluation_challenge);
        let h_2_eval = h_2_poly.evaluate(&evaluation_challenge);
        let z_eval = z_poly.evaluate(&evaluation_challenge);
        let q_eval = quotient_poly.evaluate(&evaluation_challenge);
        // Compute evaluations at `z * omega`
        let t_omega_eval = t_poly.evaluate(&evaluation_omega);
        let h_1_omega_eval = h_1_poly.evaluate(&evaluation_omega);
        let h_2_omega_eval = h_2_poly.evaluate(&evaluation_omega);
        let z_omega_eval = z_poly.evaluate(&evaluation_omega);
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
        EqualityProof {
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
                h_1: h_1_commit,
                h_2: h_2_commit,
                z: z_commit,
            },
            aggregate_witness_comm: agg_witness_comm,
            shifted_aggregate_witness_comm: shifted_agg_witness_comm,
        }
    }

    pub fn verify(
        &self,
        n: usize,
        verification_key: &VerifierKey<Bls12_381>,
        commitment_to_t: Commitment<Bls12_381>,
        transcript: &mut dyn TranscriptProtocol,
    ) -> bool {
        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(n).unwrap();

        transcript.append_commitment(b"h_1_poly", &self.commitments.h_1);
        transcript.append_commitment(b"h_2_poly", &self.commitments.h_2);
        let beta = transcript.challenge_scalar(b"beta");
        let gamma = transcript.challenge_scalar(b"gamma");
        transcript.append_commitment(b"accumulator_poly", &self.commitments.z);
        transcript.append_commitment(b"quotient_poly", &self.commitments.q);
        let evaluation_challenge = transcript.challenge_scalar(b"evaluation_challenge");
        transcript.append_scalar(b"evaluation_challenge", &evaluation_challenge);
        let evaluation_omega = evaluation_challenge * domain.group_gen;

        // Compute quotient evaluation (Q(z)) from the provers messages
        let q_eval =
            self.compute_quotient_evaluation(&beta, &gamma, &evaluation_challenge, &domain);

        transcript.append_scalar(b"f_eval", &self.evaluations.f);
        transcript.append_scalar(b"t_eval", &self.evaluations.t);
        transcript.append_scalar(b"h_1_eval", &self.evaluations.h_1);
        transcript.append_scalar(b"h_2_eval", &self.evaluations.h_2);
        transcript.append_scalar(b"z_eval", &self.evaluations.z);
        transcript.append_scalar(b"q_eval", &q_eval);
        transcript.append_scalar(b"t_omega_eval", &self.evaluations.t_omega);
        transcript.append_scalar(b"h_1_omega_eval", &self.evaluations.h_1_omega);
        transcript.append_scalar(b"h_2_omega_eval", &self.evaluations.h_2_omega);
        transcript.append_scalar(b"z_omega_eval", &self.evaluations.z_omega);

        let aggregation_challenge = transcript.challenge_scalar(b"witness_aggregation");

        // Create aggregate opening proof for all polynomials evaluated at the evaluation challenge `z`
        let agg_commitment = kzg10::aggregate_commitments(
            vec![
                &self.commitments.f,
                &commitment_to_t,
                &self.commitments.h_1,
                &self.commitments.h_2,
                &self.commitments.z,
                &self.commitments.q,
            ],
            aggregation_challenge,
        );
        let agg_value = kzg10::aggregate_values(
            vec![
                &self.evaluations.f,
                &self.evaluations.t,
                &self.evaluations.h_1,
                &self.evaluations.h_2,
                &self.evaluations.z,
                &q_eval,
            ],
            aggregation_challenge,
        );

        // Create aggregate opening proof for all polynomials evaluated at the shifted evaluation challenge `z * omega`
        let shifted_agg_commitment = kzg10::aggregate_commitments(
            vec![
                &commitment_to_t,
                &self.commitments.h_1,
                &self.commitments.h_2,
                &self.commitments.z,
            ],
            aggregation_challenge,
        );
        let shifted_agg_value = kzg10::aggregate_values(
            vec![
                &self.evaluations.t_omega,
                &self.evaluations.h_1_omega,
                &self.evaluations.h_2_omega,
                &self.evaluations.z_omega,
            ],
            aggregation_challenge,
        );

        // Batch Verify both opening proofs
        let ok = kzg10::batch_verify(
            &verification_key,
            vec![agg_commitment, shifted_agg_commitment],
            vec![
                self.aggregate_witness_comm,
                self.shifted_aggregate_witness_comm,
            ],
            vec![evaluation_challenge, evaluation_omega],
            vec![agg_value, shifted_agg_value],
        );

        ok
    }
    /// Computes the quotient evaluation from the prover messages
    fn compute_quotient_evaluation<E: EvaluationDomain<Fr>>(
        &self,
        beta: &Fr,
        gamma: &Fr,
        evaluation_challenge: &Fr,
        domain: &E,
    ) -> Fr {
        // g^{n+1}
        let last_element = domain.elements().last().unwrap();

        let lagrange_evaluations = domain.evaluate_all_lagrange_coefficients(*evaluation_challenge);
        // L_1(Z);
        let l1_z = lagrange_evaluations[0];
        // L_{n+1}(Z);
        let ln_plus_1_z = lagrange_evaluations[domain.size() - 1];

        // Z_H(Z)
        let v_h = domain.evaluate_vanishing_polynomial(*evaluation_challenge);

        let beta_one = Fr::from(1u8) + beta;
        let gamma_beta_one = (Fr::from(1u8) + beta) * gamma;

        // L_1(X) [ Z(X) -1]
        let a = { (self.evaluations.z - Fr::from(1u8)) * l1_z };

        // x-g^{n+1} * Z(X)(1+beta) * (gamma + f(x)) (gamma(1+beta) + t(x) + beta * t(Xg))
        let b = {
            let b_0 = *evaluation_challenge - last_element;
            let b_1 = self.evaluations.z * beta_one;
            let b_2 = self.evaluations.f + gamma;
            let b_3 = gamma_beta_one + self.evaluations.t + (self.evaluations.t_omega * beta);
            b_0 * b_1 * b_2 * b_3
        };
        // x-g^{n+1} * Z(Xg)[(gamma(1+beta) + h_1(X) + beta * h_1(Xg)][(gamma(1+beta) + h_2(X) + beta * h_2(Xg)]
        let c = {
            let c_0 = (*evaluation_challenge - last_element) * self.evaluations.z_omega;

            let c_1 = gamma_beta_one + self.evaluations.h_1 + (self.evaluations.h_1_omega * beta);

            let c_2 = gamma_beta_one + self.evaluations.h_2 + (self.evaluations.h_2_omega * beta);

            c_0 * c_1 * c_2
        };

        // L_{n+1}(X)[h_1(X) - h_2(Xg)]
        let d = ln_plus_1_z * (self.evaluations.h_1 - self.evaluations.h_2_omega);
        // L_{n+1}(X)[Z(X) - 1]
        let e = (self.evaluations.z - Fr::from(1u8)) * ln_plus_1_z;

        (a + b - c + d + e) / v_h
    }
}
