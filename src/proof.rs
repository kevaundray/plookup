use crate::kzg10;
use crate::lookup_table::{LookUpTable, XOR4BitTable};
use crate::multiset;
use crate::transcript::TranscriptProtocol;
use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use ff_fft::EvaluationDomain;
use poly_commit::kzg10::Commitment;
use poly_commit::kzg10::VerifierKey;
// Opening Proof consists of three components:
// - The commitment to the polynomial, which may be `None` if we have already included the commitment
// - The Witness for an opening `z`
// - The result of evaluating the committed polynomial at `z`
pub struct OpeningProof {
    pub poly_commitment: Option<Commitment<Bls12_381>>,
    pub witness_commitment: Commitment<Bls12_381>,
    pub evaluation: Fr,
}

impl OpeningProof {
    pub fn new(proof: (Option<Commitment<Bls12_381>>, Commitment<Bls12_381>, Fr)) -> OpeningProof {
        OpeningProof {
            poly_commitment: proof.0,
            witness_commitment: proof.1,
            evaluation: proof.2,
        }
    }
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
pub struct MultiSetEqualityProof {
    //Size of the domain
    // XXX: Verifier should have this value
    pub n: usize,
    // Opening proofs for polynomials evaluated at `z`
    pub h_1_proof: OpeningProof,
    pub h_2_proof: OpeningProof,
    pub z_proof: OpeningProof,
    pub t_proof: OpeningProof,
    pub f_proof: OpeningProof,
    pub q_proof: OpeningProof,
    // Opening proofs for polynomials evaluated at `z * omega`
    pub t_omega_proof: OpeningProof,
    pub h_1_omega_proof: OpeningProof,
    pub h_2_omega_proof: OpeningProof,
    pub z_omega_proof: OpeningProof,
    // XXX: We include also the evaluations for t since the verifier cannot know pre-prover, what it is.
    // Surely, we can take advantage of the fact that the verifier has the table in the preprocessing stage? Yes. The paper states this, so will need to figure this out
}

impl MultiSetEqualityProof {
    // XXX: For debugging purposes, verify either returns true or it panics. Will remove this once we aggregate all commitments
    pub fn verify(
        &self,
        verification_key: &VerifierKey<Bls12_381>,
        transcript: &mut dyn TranscriptProtocol,
    ) -> bool {
        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(self.n).unwrap();

        // Generate alpha and fold the table using alpha
        let alpha = transcript.challenge_scalar(b"alpha");

        transcript.append_commitment(b"h_1_poly", &self.h_1_proof.poly_commitment.unwrap());
        transcript.append_commitment(b"h_2_poly", &self.h_2_proof.poly_commitment.unwrap());

        let beta = transcript.challenge_scalar(b"beta");
        let gamma = transcript.challenge_scalar(b"gamma");
        transcript.append_commitment(b"accumulator_poly", &self.z_proof.poly_commitment.unwrap());
        transcript.append_commitment(b"quotient_poly", &self.q_proof.poly_commitment.unwrap());

        let evaluation_challenge = transcript.challenge_scalar(b"evaluation_challenge");
        let evaluation_omega = evaluation_challenge * domain.group_gen;

        transcript.append_scalar(b"f_eval", &self.f_proof.evaluation);
        transcript.append_scalar(b"t_eval", &self.t_proof.evaluation);
        transcript.append_scalar(b"h_1_eval", &self.h_1_proof.evaluation);
        transcript.append_scalar(b"h_2_eval", &self.h_2_proof.evaluation);
        transcript.append_scalar(b"z_eval", &self.z_proof.evaluation);
        transcript.append_scalar(b"q_eval", &self.q_proof.evaluation);
        transcript.append_scalar(b"t_omega_eval", &self.t_omega_proof.evaluation);
        transcript.append_scalar(b"h_1_omega_eval", &self.h_1_omega_proof.evaluation);
        transcript.append_scalar(b"h_2_omega_eval", &self.h_2_omega_proof.evaluation);
        transcript.append_scalar(b"z_omega_eval", &self.z_omega_proof.evaluation);

        // Check each commitment individually -- XXX: Can aggregate all of these into one using an aggregation challenge
        // h_1(X)
        let ok = kzg10::verify(
            &verification_key,
            &self.h_1_proof.poly_commitment.unwrap(),
            &self.h_1_proof.witness_commitment,
            evaluation_challenge,
            self.h_1_proof.evaluation,
        );
        assert!(ok);
        let ok = kzg10::verify(
            &verification_key,
            &self.h_1_proof.poly_commitment.unwrap(),
            &self.h_1_omega_proof.witness_commitment,
            evaluation_omega,
            self.h_1_omega_proof.evaluation,
        );
        assert!(ok);

        // h_2(X)
        let ok = kzg10::verify(
            &verification_key,
            &self.h_2_proof.poly_commitment.unwrap(),
            &self.h_2_proof.witness_commitment,
            evaluation_challenge,
            self.h_2_proof.evaluation,
        );
        assert!(ok);
        let ok = kzg10::verify(
            &verification_key,
            &self.h_2_proof.poly_commitment.unwrap(),
            &self.h_2_omega_proof.witness_commitment,
            evaluation_omega,
            self.h_2_omega_proof.evaluation,
        );
        assert!(ok);

        // Z(X)
        let ok = kzg10::verify(
            &verification_key,
            &self.z_proof.poly_commitment.unwrap(),
            &self.z_proof.witness_commitment,
            evaluation_challenge,
            self.z_proof.evaluation,
        );
        assert!(ok);
        let ok = kzg10::verify(
            &verification_key,
            &self.z_proof.poly_commitment.unwrap(),
            &self.z_omega_proof.witness_commitment,
            evaluation_omega,
            self.z_omega_proof.evaluation,
        );
        assert!(ok);

        // t(X)
        let ok = kzg10::verify(
            &verification_key,
            &self.t_proof.poly_commitment.unwrap(),
            &self.t_proof.witness_commitment,
            evaluation_challenge,
            self.t_proof.evaluation,
        );
        assert!(ok);
        let ok = kzg10::verify(
            &verification_key,
            &self.t_proof.poly_commitment.unwrap(),
            &self.t_omega_proof.witness_commitment,
            evaluation_omega,
            self.t_omega_proof.evaluation,
        );
        assert!(ok);

        // f(X)
        let ok = kzg10::verify(
            &verification_key,
            &self.f_proof.poly_commitment.unwrap(),
            &self.f_proof.witness_commitment,
            evaluation_challenge,
            self.f_proof.evaluation,
        );
        assert!(ok);

        // Q(X)
        // Compute quotient evaluation from the provers messages
        let expected_q =
            self.compute_quotient_evaluation(&beta, &gamma, &evaluation_challenge, &domain);

        assert_eq!(expected_q, self.q_proof.evaluation);

        let ok = kzg10::verify(
            &verification_key,
            &self.q_proof.poly_commitment.unwrap(),
            &self.q_proof.witness_commitment,
            evaluation_challenge,
            self.q_proof.evaluation,
        );
        assert!(ok);

        ok
    }
    /// Computes the quotient evaluation from the prover messages
    fn compute_quotient_evaluation(
        &self,
        beta: &Fr,
        gamma: &Fr,
        evaluation_challenge: &Fr,
        domain: &EvaluationDomain<Fr>,
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
        let a = { (self.z_proof.evaluation - Fr::from(1u8)) * l1_z };

        // x-g^{n+1} * Z(X)(1+beta) * (gamma + f(x)) (gamma(1+beta) + t(x) + beta * t(Xg))
        let b = {
            let b_0 = *evaluation_challenge - last_element;
            let b_1 = self.z_proof.evaluation * beta_one;
            let b_2 = self.f_proof.evaluation + gamma;
            let b_3 =
                gamma_beta_one + self.t_proof.evaluation + (self.t_omega_proof.evaluation * beta);
            b_0 * b_1 * b_2 * b_3
        };
        // x-g^{n+1} * Z(Xg)[(gamma(1+beta) + h_1(X) + beta * h_1(Xg)][(gamma(1+beta) + h_2(X) + beta * h_2(Xg)]
        let c = {
            let c_0 = (*evaluation_challenge - last_element) * self.z_omega_proof.evaluation;

            let c_1 = gamma_beta_one
                + self.h_1_proof.evaluation
                + (self.h_1_omega_proof.evaluation * beta);

            let c_2 = gamma_beta_one
                + self.h_2_proof.evaluation
                + (self.h_2_omega_proof.evaluation * beta);

            c_0 * c_1 * c_2
        };

        // L_{n+1}(X)[h_1(X) - h_2(Xg)]
        let d = ln_plus_1_z * (self.h_1_proof.evaluation - self.h_2_omega_proof.evaluation);
        // L_{n+1}(X)[Z(X) - 1]
        let e = (self.z_proof.evaluation - Fr::from(1u8)) * ln_plus_1_z;

        (a + b - c + d + e) / v_h
    }
}
