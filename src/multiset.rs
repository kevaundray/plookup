use algebra::bls12_381::Fr;
use ff_fft::{DensePolynomial as Polynomial, EvaluationDomain};

/// A MultiSet is a variation of a set, where we allow duplicate members
/// This can be emulated in Rust by using vectors
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct MultiSet(pub Vec<Fr>);

impl MultiSet {
    // Creates an empty Multiset
    pub fn new() -> MultiSet {
        MultiSet(vec![])
    }
    /// Pushes a value onto the end of the set
    pub fn push(&mut self, value: Fr) {
        self.0.push(value)
    }
    fn from_slice(slice: &[Fr]) -> MultiSet {
        MultiSet(slice.to_vec())
    }
    /// Returns the cardinality of the multiset
    pub fn len(&self) -> usize {
        self.0.len()
    }
    /// Sorts an multiset in ascending order
    pub fn sort(&self) -> MultiSet {
        let mut cloned = self.0.clone();
        cloned.sort();
        MultiSet(cloned)
    }
    /// Concatenates two sets together
    /// Does not sort the concatenated multisets
    pub fn concatenate(&self, other: &MultiSet) -> MultiSet {
        let mut result: Vec<Fr> = Vec::with_capacity(self.0.len() + other.0.len());
        result.extend(&self.0);
        result.extend(&other.0);
        MultiSet(result)
    }
    /// SortedBy checks whether every value in self appears in the same order as t
    /// Example: self = [1,2,2] t = [1,2,3] returns true
    /// Example : self = [2,1] t = [1,2] returns false
    pub fn sorted_by(&self, _t: &MultiSet) -> bool {
        todo!()
    }
    /// Splits a multiset into halves as specified by the paper
    /// If s = [1,2,3,4,5,6,7], we can deduce n using |s| = 2 * n + 1 = 7
    /// n is therefore 3
    /// We split s into two MultiSets of size n+1 each  
    /// s_0 = [1,2,3,4] ,|s_0| = n+1 = 4
    /// s_1 = [4,5,6,7] , |s_1| = n+1 = 4
    /// Notice that the last element of the first half equals the first element in the second half
    /// This is specified in the paper
    pub fn halve(&self) -> (MultiSet, MultiSet) {
        let length = self.0.len();

        let first_half = MultiSet::from_slice(&self.0[0..=length / 2]);
        let second_half = MultiSet::from_slice(&self.0[length / 2..]);

        (first_half, second_half)
    }
    /// Treats each element in the multiset as evaluation points
    /// Computes IFFT of the set of evaluation points
    /// and returns the coefficients as a Polynomial data structure
    pub fn to_polynomial(&self, domain: &EvaluationDomain<Fr>) -> Polynomial<Fr> {
        Polynomial::from_coefficients_vec(domain.ifft(&self.0))
    }
}
#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_sort() {
        let unsorted_set = MultiSet(vec![
            Fr::from(50u64),
            Fr::from(20u64),
            Fr::from(30u64),
            Fr::from(40u64),
        ]);
        let expected_sorted_multiset = MultiSet(vec![
            Fr::from(20u64),
            Fr::from(30u64),
            Fr::from(40u64),
            Fr::from(50u64),
        ]);

        let sorted_set = unsorted_set.sort();
        assert_eq!(sorted_set, expected_sorted_multiset);
        assert_ne!(sorted_set, unsorted_set);
    }

    #[test]
    fn test_concat() {
        let mut a = MultiSet::new();
        a.push(Fr::from(1u64));
        a.push(Fr::from(2u64));
        a.push(Fr::from(3u64));
        let mut b = MultiSet::new();
        b.push(Fr::from(4u64));
        b.push(Fr::from(5u64));
        b.push(Fr::from(6u64));

        let c = a.concatenate(&b);

        let expected_set = MultiSet(vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(5u64),
            Fr::from(6u64),
        ]);
        assert_eq!(expected_set, c);
    }

    #[test]
    fn test_concat_sort() {
        let mut a = MultiSet::new();
        a.push(Fr::from(2u64));
        a.push(Fr::from(2u64));
        a.push(Fr::from(3u64));
        a.push(Fr::from(1u64));
        let mut b = MultiSet::new();
        b.push(Fr::from(6u64));
        b.push(Fr::from(4u64));
        b.push(Fr::from(4u64));
        b.push(Fr::from(5u64));
        let c = a.concatenate(&b).sort();

        let expected_set = MultiSet(vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(4u64),
            Fr::from(5u64),
            Fr::from(6u64),
        ]);
        assert_eq!(expected_set, c);
    }

    #[test]
    fn test_halve() {
        let mut a = MultiSet::new();
        a.push(Fr::from(1u64));
        a.push(Fr::from(2u64));
        a.push(Fr::from(3u64));
        a.push(Fr::from(4u64));
        a.push(Fr::from(5u64));
        a.push(Fr::from(6u64));
        a.push(Fr::from(7u64));

        let (h_1, h_2) = a.halve();
        assert_eq!(h_1.len(), 4);
        assert_eq!(h_2.len(), 4);

        assert_eq!(
            MultiSet(vec![
                Fr::from(1u64),
                Fr::from(2u64),
                Fr::from(3u64),
                Fr::from(4u64)
            ]),
            h_1
        );

        assert_eq!(
            MultiSet(vec![
                Fr::from(4u64),
                Fr::from(5u64),
                Fr::from(6u64),
                Fr::from(7u64)
            ]),
            h_2
        );

        // Last element in the first half should equal first element in the second half
        assert_eq!(h_1.0.last().unwrap(), &h_2.0[0])
    }
}
