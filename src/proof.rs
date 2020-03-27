use crate::kzg10;
use crate::lookup_table::{LookUpTable, XOR4BitTable};
use crate::multiset;
use crate::transcript::TranscriptProtocol;
use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use poly_commit::kzg10::Commitment;
use poly_commit::kzg10::VerifierKey;
pub struct Evaluations {
    pub h_1_eval: Fr,
    pub h_2_eval: Fr,
    pub z_eval: Fr,
    pub q_eval: Fr,
}

// Opening Proof consists of three components:
// - The commitment to the polynomial, which may be `None` if we have already included the commitment
// - The Witness for an opening `z`
// - The result of evaluating the committed polynomial at `z`
pub struct OpeningProof(
    pub Option<Commitment<Bls12_381>>,
    pub Commitment<Bls12_381>,
    pub Fr,
);

// This protocol requires 3 extra G1 elements (Commitment)
// and 7 extra Scalar elements (Evaluations)
// - The commitment to the quotient polynomial can be batched with the PLONK quotient polynomial
// - The Witness commitments can also be batched with the PLONK opening Proof.
pub struct Proof {
    // Opening proofs for polynomials evaluated at `z`
    pub h_1_proof: OpeningProof,
    pub h_2_proof: OpeningProof,
    pub z_proof: OpeningProof,
    pub q_proof: OpeningProof,
    // Opening proofs for polynomials evaluated at `z * omega`
    pub h_1_omega_proof: OpeningProof,
    pub h_2_omega_proof: OpeningProof,
    pub z_omega_proof: OpeningProof,
}
