use ark_bls12_381::{Bls12_381, Fr, G1Projective};
use ark_ec::AffineCurve;
use ark_poly::{polynomial::univariate::DensePolynomial as Polynomial, UVPolynomial};
use ark_poly_commit::kzg10::{Commitment, Powers, Proof, UniversalParams, VerifierKey, KZG10};
use ark_std::test_rng;
use num_traits::identities::Zero;

type KzgBls12_381 = KZG10<Bls12_381, Polynomial<Fr>>;

pub fn trusted_setup<'a>(max_deg: usize) -> (Powers<'a, Bls12_381>, VerifierKey<Bls12_381>) {
    let mut rng = test_rng();
    let pp = KzgBls12_381::setup(max_deg, false, &mut rng).unwrap();
    trim(&pp, max_deg)
}

fn trim<'a>(
    pp: &UniversalParams<Bls12_381>,
    mut supported_degree: usize,
) -> (Powers<'a, Bls12_381>, VerifierKey<Bls12_381>) {
    if supported_degree == 1 {
        supported_degree += 1;
    }
    let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
    let powers_of_gamma_g = (0..=supported_degree)
        .map(|i| pp.powers_of_gamma_g[&i])
        .collect();

    let powers = Powers {
        powers_of_g: std::borrow::Cow::Owned(powers_of_g),
        powers_of_gamma_g: std::borrow::Cow::Owned(powers_of_gamma_g),
    };
    let vk = VerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[&0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    (powers, vk)
}

pub fn commit(powers: &Powers<Bls12_381>, p: &Polynomial<Fr>) -> Commitment<Bls12_381> {
    let hiding_bound = None;
    let (comm, _) = KZG10::commit(&powers, p, hiding_bound, None).unwrap();
    comm
}

pub fn commit_vec(powers: &Powers<Bls12_381>, p_vec: &Vec<Fr>) -> Commitment<Bls12_381> {
    let p = Polynomial::from_coefficients_slice(p_vec);
    commit(powers, &p)
}

// Computes the witness for a set of polynomials evaluated at the same point
// W(X) = f(x) - f(z) / x-z
// However, the quotient is invariant under `f(z)`,
// So we can compute the witness as f(x) / x-z
pub fn compute_witness(polynomial: &Polynomial<Fr>, point: Fr) -> Polynomial<Fr> {
    let divisor = Polynomial::from_coefficients_vec(vec![-point, Fr::from(1u8)]);
    polynomial / &divisor
}
// For some challenge v, a list of polynomials p_i and a point z
// We compute the aggregate witness as (v^0 * p_0 + v^1 * p_1 + ...+ v^n * p_n ) / x-z
pub fn compute_aggregate_witness(
    polynomials: Vec<&Polynomial<Fr>>,
    point: Fr,
    aggregation_challenge: Fr,
) -> Polynomial<Fr> {
    let mut powers = Fr::from(1u8);
    let mut result = Polynomial::zero();

    for polynomial in polynomials {
        let intermediate_poly = polynomial * &Polynomial::from_coefficients_vec(vec![powers]);
        result += &intermediate_poly;
        powers = powers * aggregation_challenge;
    }

    let divisor = Polynomial::from_coefficients_vec(vec![-point, Fr::from(1u8)]);

    &result / &divisor
}

pub fn aggregate_commitments(
    commitments: Vec<&Commitment<Bls12_381>>,
    aggregation_challenge: Fr,
) -> Commitment<Bls12_381> {
    let mut powers = Fr::from(1u8);
    let mut result = G1Projective::zero();

    for commitment in commitments {
        let intermediate_comm = commitment.0.mul(powers);
        result += &intermediate_comm;
        powers = powers * aggregation_challenge;
    }

    Commitment(result.into())
}
pub fn aggregate_values(values: Vec<&Fr>, aggregation_challenge: Fr) -> Fr {
    let mut powers = Fr::from(1u8);
    let mut result = Fr::zero();

    for value in values {
        let intermediate_value = *value * powers;
        result += &intermediate_value;
        powers = powers * aggregation_challenge;
    }

    result
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
        random_v: Some(Fr::zero()),
    };

    KzgBls12_381::check(vk, commitment_to_poly, evaluation_point, value, &proof).unwrap()
}

pub fn batch_verify(
    vk: &VerifierKey<Bls12_381>,
    commitment_to_polynomials: Vec<Commitment<Bls12_381>>,
    commitment_to_witnesses: Vec<Commitment<Bls12_381>>,
    evaluation_points: Vec<Fr>,
    values: Vec<Fr>,
) -> bool {
    let mut proofs: Vec<Proof<Bls12_381>> = Vec::new();
    for witness in commitment_to_witnesses {
        let proof = Proof {
            w: witness.0,
            random_v: Some(Fr::zero()),
        };
        proofs.push(proof);
    }

    KzgBls12_381::batch_check(
        vk,
        commitment_to_polynomials.as_slice(),
        &evaluation_points,
        &values,
        &proofs,
        &mut test_rng(),
    )
    .unwrap()
}
