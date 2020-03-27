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

// This protocol requires 3 extra G1 elements (Commitment)
// and 7 extra Scalar elements (Evaluations)
// - The commitment to the quotient polynomial can be batched with the PLONK quotient polynomial
// - The Witness commitments can also be batched with the PLONK opening Proof.
pub struct MultiSetEqualityProof {
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
