use crate::kzg10;
use crate::lookup_table::{LookUpTable, XOR4BitTable};
use crate::multiset;
use crate::transcript::TranscriptProtocol;
use algebra::bls12_381::Fr;
use algebra::Bls12_381;
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
    // Surely, we can take advantage of the fact that the verifier has the table in the preprocessing stage?
}
