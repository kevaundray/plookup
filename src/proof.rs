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
// This protocol requires 3 extra G1 elements (Commitment)
// and 4 extra Scalar elements (Evaluations)
// - The commitment to the quotient polynomial can be batched with the PLONK quotient polynomial
// - The Opening Proof can also be batched with the PLONK opening Proof.
//
// Also note this is a Proof of Concept, which is why we have a separate proof struct
// In reality, we will adds these to the PLONK proof data structure
pub struct Proof {
    // Two commitments to h_1 and h_2
    pub h_1_comm: Commitment<Bls12_381>,
    pub h_1_witness_comm: Commitment<Bls12_381>,
    pub h_2_comm: Commitment<Bls12_381>,
    pub h_2_witness_comm: Commitment<Bls12_381>,
    // Commitment to Z(X) ; the accumulator polynomial
    pub z_comm: Commitment<Bls12_381>,
    pub z_witness_comm: Commitment<Bls12_381>,
    // Commitment to Q(X); the quotient polynomial
    pub q_comm: Commitment<Bls12_381>,
    pub q_witness_comm: Commitment<Bls12_381>,
    // Commitment to the witness polynomial
    pub evaluations: Evaluations,
}
