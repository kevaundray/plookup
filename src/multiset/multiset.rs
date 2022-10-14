use ark_bls12_381::Fr;
use ark_ff::{One, Zero};
use ark_poly::{polynomial::univariate::DensePolynomial as Polynomial, EvaluationDomain, UVPolynomial};
use std::ops::{Add, Mul};
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
    /// Pushes 'n' elements into the multiset
    pub fn extend(&mut self, n: usize, value: Fr) {
        let elements = vec![value; n];
        self.0.extend(elements);
    }
    /// Fetches last element in multiset
    /// Panics if there are no elements
    pub fn last(&self) -> Fr {
        *self.0.last().unwrap()
    }
    fn from_slice(slice: &[Fr]) -> MultiSet {
        MultiSet(slice.to_vec())
    }
    /// Returns the cardinality of the multiset
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Concatenates two sets together
    /// Does not sort the concatenated multisets
    pub fn concatenate(&self, other: &MultiSet) -> MultiSet {
        let mut result: Vec<Fr> = Vec::with_capacity(self.0.len() + other.0.len());
        result.extend(&self.0);
        result.extend(&other.0);
        MultiSet(result)
    }

    /// Returns the position of the element in the Multiset
    /// Panics if element is not in the Multiset
    fn position(&self, element: &Fr) -> usize {
        let index = self.0.iter().position(|&x| x == *element).unwrap();
        index
    }

    /// Performs a element-wise insertion into the second multiset
    /// Example: f {1,2,3,1} t : {3,1,2,3}
    /// We now take each element from f and find the element in `t` then insert the element from `f` into right next to it's duplicate
    /// We are assuming that `f` is contained in `t`
    pub fn concatenate_and_sort(&self, t: &MultiSet) -> MultiSet {
        assert!(self.is_subset_of(t));
        let mut result = t.clone();

        for element in self.0.iter() {
            let index = result.position(element);
            result.0.insert(index, *element);
        }

        result
    }

    /// Checks whether self is a subset of other
    pub fn is_subset_of(&self, other: &MultiSet) -> bool {
        let mut is_subset = true;

        for x in self.0.iter() {
            is_subset = other.contains(x);
            if is_subset == false {
                break;
            }
        }
        is_subset
    }
    /// Checks if an element is in the MultiSet
    pub fn contains(&self, element: &Fr) -> bool {
        self.0.contains(element)
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
    pub fn to_polynomial<E :EvaluationDomain<Fr>>(&self, domain: &E) -> Polynomial<Fr> {
        Polynomial::from_coefficients_vec(domain.ifft(&self.0))
    }
    /// Aggregates multisets together using a random challenge
    /// Eg. for three sets A,B,C and a random challenge `k`
    /// The aggregate is k^0 *A + k^1 * B + k^2 * C
    pub fn aggregate(sets: Vec<&MultiSet>, challenge: Fr) -> MultiSet {
        // First find the set with the most elements
        let mut max = 0usize;
        for set in sets.iter() {
            if set.len() > max {
                max = set.len()
            }
        }

        let mut result = MultiSet(vec![Fr::zero(); max]);
        let mut powers = Fr::one();

        for set in sets {
            let intermediate_set = set * powers;

            result = result + intermediate_set;

            powers = powers * challenge;
        }

        result
    }
}

impl Add for MultiSet {
    type Output = MultiSet;
    fn add(self, other: MultiSet) -> Self::Output {
        let result = self
            .0
            .into_iter()
            .zip(other.0.iter())
            .map(|(x, y)| x + y)
            .collect();

        MultiSet(result)
    }
}
impl Mul<Fr> for MultiSet {
    type Output = MultiSet;
    fn mul(self, other: Fr) -> Self::Output {
        let result = self.0.into_iter().map(|x| x * other).collect();
        MultiSet(result)
    }
}
impl Mul<Fr> for &MultiSet {
    type Output = MultiSet;
    fn mul(self, other: Fr) -> Self::Output {
        let result = self.0.iter().map(|x| other * x).collect();
        MultiSet(result)
    }
}
#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_concatenate() {
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

    #[test]
    fn test_to_polynomial() {

        let mut a = MultiSet::new();
        a.push(Fr::from(1u8));
        a.push(Fr::from(2u8));
        a.push(Fr::from(3u8));
        a.push(Fr::from(4u8));
        a.push(Fr::from(5u8));
        a.push(Fr::from(6u8));
        a.push(Fr::from(7u8));

        let domain = EvaluationDomain::new(a.len() + 1).unwrap();
        let a_poly = a.to_polynomial(&domain);

        assert_eq!(a_poly.degree(), 7)
    }
    #[test]
    fn test_is_subset() {
        let mut a = MultiSet::new();
        a.push(Fr::from(1u8));
        a.push(Fr::from(2u8));
        a.push(Fr::from(3u8));
        a.push(Fr::from(4u8));
        a.push(Fr::from(5u8));
        a.push(Fr::from(6u8));
        a.push(Fr::from(7u8));
        let mut b = MultiSet::new();
        b.push(Fr::from(1u8));
        b.push(Fr::from(2u8));
        let mut c = MultiSet::new();
        c.push(Fr::from(100u8));

        assert!(b.is_subset_of(&a));
        assert!(!c.is_subset_of(&a));
    }
    #[test]
    fn test_sort_by() {
        let mut f = MultiSet::new();
        f.push(Fr::from(2u8));
        f.push(Fr::from(1u8));
        f.push(Fr::from(2u8));
        f.push(Fr::from(4u8));
        f.push(Fr::from(3u8));

        let mut t = MultiSet::new();
        t.push(Fr::from(3u8));
        t.push(Fr::from(1u8));
        t.push(Fr::from(2u8));
        t.push(Fr::from(4u8));

        let sorted_s = f.concatenate_and_sort(&t);

        let mut expected_sorted_s = MultiSet::new();
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(1u8));
        expected_sorted_s.push(Fr::from(1u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(4u8));
        expected_sorted_s.push(Fr::from(4u8));

        assert_eq!(expected_sorted_s, sorted_s);
    }

    #[test]
    fn test_concate_sort() {
        let mut f = MultiSet::new();
        f.push(Fr::from(1u8));
        f.push(Fr::from(2u8));

        // Table of values
        let mut t = MultiSet::new();
        t.push(Fr::from(2u8));
        t.push(Fr::from(1u8));
        t.push(Fr::from(2u8));

        let mut expected_s = MultiSet::new();
        expected_s.push(Fr::from(2u8));
        expected_s.push(Fr::from(2u8));
        expected_s.push(Fr::from(1u8));
        expected_s.push(Fr::from(1u8));
        expected_s.push(Fr::from(2u8));

        let sorted_s_1 = f.concatenate_and_sort(&t);
        assert_eq!(sorted_s_1, expected_s)
    }
    #[test]
    fn test_concate_sort2() {
        let mut f = MultiSet::new();
        f.push(Fr::from(1u8));
        f.push(Fr::from(1u8));
        f.push(Fr::from(2u8));
        f.push(Fr::from(2u8));
        f.push(Fr::from(3u8));
        f.push(Fr::from(3u8));
        f.push(Fr::from(3u8));

        // Table of values
        let mut t = MultiSet::new();
        t.push(Fr::from(1u8));
        t.push(Fr::from(1u8));
        t.push(Fr::from(2u8));
        t.push(Fr::from(2u8));
        t.push(Fr::from(3u8));
        t.push(Fr::from(3u8));
        t.push(Fr::from(4u8));
        t.push(Fr::from(5u8));

        let mut expected_sorted_s = MultiSet::new();
        expected_sorted_s.push(Fr::from(1u8));
        expected_sorted_s.push(Fr::from(1u8));
        expected_sorted_s.push(Fr::from(1u8));
        expected_sorted_s.push(Fr::from(1u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(2u8));
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(3u8));
        expected_sorted_s.push(Fr::from(4u8));
        expected_sorted_s.push(Fr::from(5u8));

        let sorted_s = f.concatenate_and_sort(&t);
        assert_eq!(sorted_s.len(), f.len() + t.len());
        assert_eq!(sorted_s.len(), expected_sorted_s.len());
        assert_eq!(sorted_s, expected_sorted_s)
    }
}
