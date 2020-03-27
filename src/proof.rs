use algebra::Bls12_381;
use poly_commit::kzg10::Commitment;
// This protocol requires 3 extra G1 elements (Commitment)
// This is a Proof of Concept, which is why we have a separate proof struct
// In reality, we will adds these to the PLONK proof data structure
pub struct Proof {
    // Two commitments to h_1 and h_2
    pub h_1_comm: Commitment<Bls12_381>,
    pub h_2_comm: Commitment<Bls12_381>,
    // Commitment to Z
    pub z_comm: Commitment<Bls12_381>,
    // Commitment to the quotient polynomial
    pub q_comm: Commitment<Bls12_381>,
}
