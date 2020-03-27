use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use algebra::{to_bytes, ToBytes};
use merlin::Transcript;
use poly_commit::kzg10::Commitment;

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
        use algebra::UniformRand;
        use rand_chacha::ChaChaRng;
        use rand_core::SeedableRng;

        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);

        let mut rng = &mut self.build_rng().finalize(&mut ChaChaRng::from_seed(buf));
        Fr::rand(&mut rng)
    }
}
