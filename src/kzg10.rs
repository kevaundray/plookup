use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use ff_fft::DensePolynomial as Polynomial;
use poly_commit::kzg10::{Commitment, Powers, UniversalParams, VerifierKey, KZG10};
use rand_core::RngCore;

// Modification of https://github.com/scipr-lab/poly-commit/blob/master/src/kzg10/mod.rs
type KZG_Bls12_381 = KZG10<Bls12_381>;

pub fn trusted_setup<R: RngCore>(max_deg: usize, rng: &mut R) -> UniversalParams<Bls12_381> {
    KZG_Bls12_381::setup(max_deg, false, rng).unwrap()
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
