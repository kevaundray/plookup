use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{to_bytes, UniformRand};
use ark_poly_commit::kzg10::Commitment;
use ark_std::test_rng;
use merlin::Transcript;

pub trait TranscriptProtocol {
    /// Append a `commitment` with the given `label`.
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<Bls12_381>);

    /// Append a `Scalar` with the given `label`.
    fn append_scalar(&mut self, label: &'static [u8], s: &Fr);

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> Fr;
}

impl TranscriptProtocol for Transcript {
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<Bls12_381>) {
        self.append_message(label, &to_bytes![comm].unwrap());
    }

    fn append_scalar(&mut self, label: &'static [u8], s: &Fr) {
        self.append_message(label, &to_bytes![s].unwrap())
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> Fr {
        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);

        let mut rng = test_rng();
        Fr::rand(&mut rng)
    }
}
