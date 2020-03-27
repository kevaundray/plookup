use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use ff_fft::DensePolynomial as Polynomial;
use poly_commit::kzg10::{Commitment, Powers, Proof, UniversalParams, VerifierKey, KZG10};
use rand_chacha::ChaChaRng;
use rand_core::SeedableRng;
// Modification of https://github.com/scipr-lab/poly-commit/blob/master/src/kzg10/mod.rs
type KZG_Bls12_381 = KZG10<Bls12_381>;

pub fn trusted_setup(max_deg: usize, seed: &[u8]) -> UniversalParams<Bls12_381> {
    let mut rng = ChaChaRng::from_seed(to_32_bytes(seed));
    KZG_Bls12_381::setup(max_deg, false, &mut rng).unwrap()
}

fn to_32_bytes(bytes: &[u8]) -> [u8; 32] {
    let mut array: [u8; 32] = [0; 32];
    for (a, b) in bytes.iter().zip(array.iter_mut()) {
        *b = *a
    }
    array
}

pub fn trim<'a>(
    pp: &UniversalParams<Bls12_381>,
    mut supported_degree: usize,
) -> (Powers<'a, Bls12_381>, VerifierKey<Bls12_381>) {
    if supported_degree == 1 {
        supported_degree += 1;
    }
    let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
    let powers_of_gamma_g = pp.powers_of_gamma_g[..=supported_degree].to_vec();

    let powers = Powers {
        powers_of_g: std::borrow::Cow::Owned(powers_of_g),
        powers_of_gamma_g: std::borrow::Cow::Owned(powers_of_gamma_g),
    };
    let vk = VerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    (powers, vk)
}

pub fn commit(powers: &Powers<Bls12_381>, p: &Polynomial<Fr>) -> Commitment<Bls12_381> {
    let hiding_bound = None;
    let (comm, _) = KZG10::commit(&powers, &p, hiding_bound, None).unwrap();
    comm
}

// Computes the witness for a set of polynomials evaluated at the same point
// W(X) = f(x) - f(z) / x-z
// However, the quotient is invariant under `f(z)`,
// So we can compute the witness as f(x) / x-z
pub fn compute_witness(polynomial: &Polynomial<Fr>, point: Fr) -> Polynomial<Fr> {
    let divisor = Polynomial::from_coefficients_vec(vec![-point, Fr::from(1u8)]);

    polynomial / &divisor
}

pub fn verify(
    vk: &VerifierKey<Bls12_381>,
    commitment_to_poly: &Commitment<Bls12_381>,
    commitment_to_witness: &Commitment<Bls12_381>,
    evaluation_point: Fr,
    value: Fr,
) -> bool {
    let proof = Proof {
        w: commitment_to_witness.0,
        random_v: Fr::from(0u8),
    };

    KZG10::check(vk, commitment_to_poly, evaluation_point, value, &proof).unwrap()
}
