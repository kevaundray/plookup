use algebra::bls12_381::Fr;
use ff_fft::{DensePolynomial as Polynomial, EvaluationDomain};
use num_traits::identities::{One, Zero};
// The quotient polynomial will encode the four checks for the multiset equality argument
// These checks are:
// 1) Z(X) evaluated at the first root of unity is 1
// 2) Z(X) is correct accumulated. z_(Xg) * g(X) = (1+beta)^n * z(X) * f(X)
// 3) The last element of h_1(x) is equal to the first element of h_2(x)
// 4) Z(x) evaluated at the last root of unity is 1
//
// We can denote check 1 and check 4 as point checks because they are checking the evaluation of Z(x) at a specific point
// We can denote check 3 as an interval check because it checks whether h_1 and h_2 combined form 's' without any gaps. See paper for more details on 's'
// We can denote check 2 as the term check
//
// Notice that the term check equation will determine the degree of the quotient polynomial
// We can compute it by adding the degrees of Z(x), f(x) and t(x).
// deg(Z(x)) = n because it has n + 1 elements
// deg(f(x)) = n  although it has n elements, we must zero pad to ensure that f(x) evaluated on the n+1'th element is zero
// deg(t(x)) = n because we define it to have n + 1 elements.
// Summing the degrees gives us n + n + n = 3n
// However, similar to [GWC19](PLONK) we must divide by the vanishing polynomial
// So the degree of the quotient polynomial Q(x) is 3n - n = 2n
// Significance: Adding this protocol into PLONK will not "blow up" the degree of the quotient polynomial
// Where "blow up" denotes increasing the overall degree past 4n for standard plonk
pub fn compute(
    domain: &EvaluationDomain<Fr>,
    z_poly: &Polynomial<Fr>,
    f_poly: &Polynomial<Fr>,
    t_poly: &Polynomial<Fr>,
    h_1_poly: &Polynomial<Fr>,
    h_2_poly: &Polynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> (Polynomial<Fr>, Polynomial<Fr>) {
    // 1. Compute Point check polynomial
    let point_check = compute_point_checks(z_poly, domain);
    //2. Compute interval check polynomial
    let interval_check = compute_interval_check(h_1_poly, h_2_poly, domain);
    //3. Compute term check polynomial
    let term_check = compute_term_check(
        domain, z_poly, f_poly, t_poly, h_1_poly, h_2_poly, beta, gamma,
    );
    // Compute quotient polynomial
    let sum = &(&interval_check + &point_check) + &term_check;

    sum.divide_by_vanishing_poly(*domain).unwrap()
}

fn compute_point_checks(z_poly: &Polynomial<Fr>, domain: &EvaluationDomain<Fr>) -> Polynomial<Fr> {
    // Compute lagrange polynomials
    let l1_poly = compute_n_lagrange_poly(domain, 0);
    let ln_poly = compute_n_lagrange_poly(domain, domain.size() - 1);

    // Compute Z'(X) = Z(x) - 1
    let z_prime_poly = z_poly - &Polynomial::from_coefficients_vec(vec![Fr::one()]);

    // We can batch the two point checks into one with the following: (Z(X)-1)[L_1(x) + L_n(x)]
    //
    let l_poly = &l1_poly + &ln_poly;
    &z_prime_poly * &l_poly
}

fn compute_interval_check(
    h_1_poly: &Polynomial<Fr>,
    h_2_poly: &Polynomial<Fr>,
    domain: &EvaluationDomain<Fr>,
) -> Polynomial<Fr> {
    // Increase domain size by two
    let domain_2n: EvaluationDomain<Fr> = EvaluationDomain::new(2 * domain.size()).unwrap();

    // Compute last lagrange polynomial in evaluation form
    let ln_evals = compute_n_lagrange_evaluations(domain_2n.size(), domain_2n.size() - 1);

    // Convert h_1 and h_2 to evaluation form
    let h_1_evals = domain_2n.fft(&h_1_poly);
    let mut h_2_evals = domain_2n.fft(&h_2_poly);
    // We need h_2(x * g) so push 2 extra elements into the domain
    h_2_evals.push(h_2_evals[0]);
    h_2_evals.push(h_2_evals[1]);

    // Compute [L_n(x)](h_1(x) - h_2(x * g))
    let i_evals: Vec<_> = (0..domain_2n.size())
        .into_iter()
        .map(|i| {
            let ln_i = ln_evals[i];
            let h_1_i = h_1_evals[i];
            let h_2_i_next = h_2_evals[i + 2];

            ln_i * (h_1_i - h_2_i_next)
        })
        .collect();

    // Convert the evaluations for our point check to coefficient form
    let i_poly = Polynomial::from_coefficients_vec(domain_2n.ifft(&i_evals));
    i_poly
}

fn compute_term_check(
    domain: &EvaluationDomain<Fr>,
    z_poly: &Polynomial<Fr>,
    f_poly: &Polynomial<Fr>,
    t_poly: &Polynomial<Fr>,
    h_1_poly: &Polynomial<Fr>,
    h_2_poly: &Polynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> Polynomial<Fr> {
    // The equation for this is quite big. Similar to PLONK, we can split the point check into two.
    // The first part will compute the grand product Z(X) term
    // The second part will compute the grand product Z(Xg) term

    // First Part
    let part_a = compute_term_check_a(domain, z_poly, f_poly, t_poly, beta, gamma);
    // Second part
    let part_b = compute_term_check_b(domain, z_poly, h_1_poly, h_2_poly, beta, gamma);

    &part_a - &part_b
}
// This computes the grand product term for Z(X) or F(\beta, \gamma)
pub fn compute_term_check_a(
    domain: &EvaluationDomain<Fr>,
    z_poly: &Polynomial<Fr>,
    f_poly: &Polynomial<Fr>,
    t_poly: &Polynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> Polynomial<Fr> {
    // Increase the domain size by 4
    let domain_4n: &EvaluationDomain<Fr> = &EvaluationDomain::new(4 * domain.size()).unwrap();

    // Convert all polynomials into evaluation form
    let z_evals = domain_4n.fft(&z_poly);
    let f_evals = domain_4n.fft(f_poly);
    let mut t_evals = domain_4n.fft(t_poly);
    // Add four terms to the t(x) evaluations as we need to compute t(Xg)
    t_evals.push(t_evals[0]);
    t_evals.push(t_evals[1]);
    t_evals.push(t_evals[2]);
    t_evals.push(t_evals[3]);

    let beta_one = Fr::one() + beta;

    // Compute the last element in the domain
    let g_n = domain_4n.elements().last().unwrap();

    let i_evals: Vec<_> = (0..domain_4n.size())
        .into_iter()
        .zip(domain_4n.elements())
        .map(|(i, root_i)| {
            let z_i = z_evals[i];
            let f_i = f_evals[i];
            let t_i = t_evals[i];
            let t_i_next = t_evals[i + 4];

            // Compute X - g^n
            let a = root_i - g_n;

            // Compute Z(X)(1+beta)
            let b = z_i * beta_one;

            // Compute gamma + f(X)
            let c = gamma + f_i;

            // Compute gamma(1+beta) +t(x) + beta * t(Xg)
            let d = (gamma * beta_one) + t_i + (beta * t_i_next);

            a * b * c * d
        })
        .collect();

    // Convert the evaluations for our term check to coefficient form
    let i_poly = Polynomial::from_coefficients_vec(domain_4n.ifft(&i_evals));

    assert_eq!(
        i_poly.evaluate(domain_4n.elements().last().unwrap()),
        Fr::zero()
    );

    i_poly
}
// This computes the grand product term for Z(Xg) or G(\beta, \gamma)
fn compute_term_check_b(
    domain: &EvaluationDomain<Fr>,
    z_poly: &Polynomial<Fr>,
    h_1_poly: &Polynomial<Fr>,
    h_2_poly: &Polynomial<Fr>,
    beta: Fr,
    gamma: Fr,
) -> Polynomial<Fr> {
    // Increase the domain size by 4
    let domain_4n: &EvaluationDomain<Fr> = &EvaluationDomain::new(4 * domain.size()).unwrap();

    // Convert all polynomials into evaluation form, then add four terms to each evaluation as we need to compute their evaluations at the next root of unity
    let mut z_evals = domain_4n.fft(z_poly);
    z_evals.push(z_evals[0]);
    z_evals.push(z_evals[1]);
    z_evals.push(z_evals[2]);
    z_evals.push(z_evals[3]);
    let mut h_1_evals = domain_4n.fft(h_1_poly);
    h_1_evals.push(h_1_evals[0]);
    h_1_evals.push(h_1_evals[1]);
    h_1_evals.push(h_1_evals[2]);
    h_1_evals.push(h_1_evals[3]);

    let mut h_2_evals = domain_4n.fft(h_2_poly);
    h_2_evals.push(h_2_evals[0]);
    h_2_evals.push(h_2_evals[1]);
    h_2_evals.push(h_2_evals[2]);
    h_2_evals.push(h_2_evals[3]);

    // Compute (1 + beta)
    let beta_one = Fr::one() + beta;
    // Compute the last element in the domain
    let g_n = domain_4n.elements().last().unwrap();

    let i_evals: Vec<_> = (0..domain_4n.size())
        .into_iter()
        .zip(domain_4n.elements())
        .map(|(i, root_i)| {
            let z_i_next = z_evals[i + 4];
            let h_1_i = h_1_evals[i];
            let h_1_i_next = h_1_evals[i + 4];
            let h_2_i = h_2_evals[i];
            let h_2_i_next = h_2_evals[i + 4];

            // Compute (X - g^n) Z(Xg)
            let a = (root_i - g_n) * z_i_next;

            // Compute [gamma * (1+beta)] + h_1(x) + beta * h_1(Xg)
            let b = (gamma * beta_one) + h_1_i + (beta * h_1_i_next);

            // Compute [gamma * (1+beta)] + h_2(x) + beta * h_2(Xg)
            let c = (gamma * beta_one) + h_2_i + (beta * h_2_i_next);

            a * b * c
        })
        .collect();

    // Convert the evaluations for our term check to coefficient form
    let i_poly = Polynomial::from_coefficients_vec(domain_4n.ifft(&i_evals));

    assert_eq!(
        i_poly.evaluate(domain_4n.elements().last().unwrap()),
        Fr::zero()
    );

    i_poly
}

// Computes the n'th lagrange poly for a particular domain
// Easiest way is to compute the evaluation points, which will be zero at every position except for n
// Then IFFT to get the coefficient form
// Note: n=0 is the first lagrange polynomial and n = domain.size() -1 is the last lagrange polynomial
fn compute_n_lagrange_poly(domain: &EvaluationDomain<Fr>, n: usize) -> Polynomial<Fr> {
    assert!(n <= domain.size() - 1);
    let mut evaluations = compute_n_lagrange_evaluations(domain.size(), n);
    domain.ifft_in_place(&mut evaluations);
    Polynomial::from_coefficients_vec(evaluations)
}
fn compute_n_lagrange_evaluations(domain_size: usize, n: usize) -> Vec<Fr> {
    let mut lagrange_evaluations = vec![Fr::zero(); domain_size];
    lagrange_evaluations[n] = Fr::one();
    lagrange_evaluations
}

mod test {
    use super::*;
    use crate::multiset::MultiSet;
    use crate::multiset_equality::*;
    #[test]
    #[ignore]
    fn test_quotient_poly() {
        // Compute f
        let mut f = MultiSet::new();
        f.push(Fr::from(2u8));
        f.push(Fr::from(3u8));
        f.push(Fr::from(4u8));

        // Compute t
        let mut t = MultiSet::new();
        t.push(Fr::from(2u8));
        t.push(Fr::from(3u8));
        t.push(Fr::from(4u8));
        t.push(Fr::from(5u8));

        // Setup domain
        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(f.len()).unwrap();
        let domain_2n: EvaluationDomain<Fr> = EvaluationDomain::new(2 * f.len()).unwrap();
        let beta = Fr::from(10u8);
        let gamma = Fr::from(11u8);

        // Compute h_1 and h_2
        let (h_1, h_2) = compute_h1_h2(&f, &t);

        // Convert h_1 and h_2 to polynomials
        let h_1_poly = h_1.to_polynomial(&domain);
        let h_2_poly = h_2.to_polynomial(&domain);

        // Compute f(x)
        let f_poly = f.to_polynomial(&domain);
        assert_eq!(f_poly.degree(), f.len());
        // Compute t(x)
        let t_poly = t.to_polynomial(&domain_2n);
        assert_eq!(t_poly.degree(), t.len());

        // Compute Z(x) poly
        let z_evaluations = compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
        let z_poly = Polynomial::from_coefficients_vec(domain_2n.ifft(&z_evaluations));

        let (q, remainder) = compute(
            &domain, &z_poly, &f_poly, &t_poly, &h_1_poly, &h_2_poly, beta, gamma,
        );
        // XXX: Last thing to do, the term check validation seems to produce a non-zero remainder.
        assert!(remainder.is_zero());
    }
}
