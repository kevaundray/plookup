use crate::{kzg10, multiset::MultiSet};
use ark_bls12_381::{Bls12_381, Fr};
use ark_poly::{
    polynomial::univariate::DensePolynomial as Polynomial, EvaluationDomain,
    Radix2EvaluationDomain, UVPolynomial,
};
use ark_poly_commit::kzg10::{Commitment, Powers};
use std::collections::HashMap;

pub mod four_bits;
pub mod generic;
pub use generic::Generic;

pub struct PreProcessedTable {
    pub n: usize,
    pub t_1: (MultiSet, Commitment<Bls12_381>, Polynomial<Fr>),
    pub t_2: (MultiSet, Commitment<Bls12_381>, Polynomial<Fr>),
    pub t_3: (MultiSet, Commitment<Bls12_381>, Polynomial<Fr>),
}

pub trait LookUpTable {
    /// Returns the number of entries in the lookup table
    fn len(&self) -> usize {
        self.borrow_map().keys().len()
    }

    /// We represent the lookup table as a map
    /// Returns an immutable copy of the map
    fn borrow_map(&self) -> &HashMap<(Fr, Fr), Fr>;

    /// Fetches the appropriate lookup table value
    /// Given its input
    fn read(&self, key: &(Fr, Fr)) -> Option<&Fr> {
        self.borrow_map().get(key)
    }

    /// Given a lookup table where each row contains three entries (a,b,c)
    /// Create three multisets of the form
    /// a = {a_0, a_1, a_2, a_3,...,a_n}
    /// b = {b_0, b_1, b_2, b_3,...,b_n}
    /// c = {c_0, c_1, c_2, c_3,...,c_n}
    fn to_multiset(&self) -> (MultiSet, MultiSet, MultiSet) {
        let mut table_multiset_left = MultiSet::new();
        let mut table_multiset_right = MultiSet::new();
        let mut table_multiset_out = MultiSet::new();

        for (key, value) in self.borrow_map().iter() {
            let input_0 = key.0;
            let input_1 = key.1;
            let output = *value;

            table_multiset_left.push(input_0);
            table_multiset_right.push(input_1);
            table_multiset_out.push(output);
        }

        (
            table_multiset_left,
            table_multiset_right,
            table_multiset_out,
        )
    }
    /// Pre-process a table by padding it to a size `n` commitment to each column in the table
    /// `n` will usually be equal to the size of your circuit, when padded.
    fn preprocess(&self, commit_key: &Powers<Bls12_381>, n: usize) -> PreProcessedTable {
        assert!(n.is_power_of_two());

        let (mut t_1, mut t_2, mut t_3) = self.to_multiset();

        let k = t_1.len();
        assert_eq!(t_1.len(), k);
        assert_eq!(t_2.len(), k);
        assert_eq!(t_3.len(), k);

        let domain: Radix2EvaluationDomain<Fr> = EvaluationDomain::new(n).unwrap();

        // Pad
        let pad_by = n - t_1.len();
        t_1.extend(pad_by, t_1.last());
        t_2.extend(pad_by, t_2.last());
        t_3.extend(pad_by, t_3.last());

        let t_1_poly = Polynomial::from_coefficients_vec(domain.ifft(&t_1.0));
        let t_2_poly = Polynomial::from_coefficients_vec(domain.ifft(&t_2.0));
        let t_3_poly = Polynomial::from_coefficients_vec(domain.ifft(&t_3.0));

        let t_1_commit = kzg10::commit(commit_key, &t_1_poly);
        let t_2_commit = kzg10::commit(commit_key, &t_2_poly);
        let t_3_commit = kzg10::commit(commit_key, &t_3_poly);

        PreProcessedTable {
            n: n,
            t_1: (t_1, t_1_commit, t_1_poly),
            t_2: (t_2, t_2_commit, t_2_poly),
            t_3: (t_3, t_3_commit, t_3_poly),
        }
    }
}
