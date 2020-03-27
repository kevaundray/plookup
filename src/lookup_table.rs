use crate::multiset::MultiSet;
use algebra::bls12_381::Fr;
use std::collections::HashMap;

pub trait LookUpTable {
    /// Creates a new lookup table with its entries populated
    fn new() -> Self;

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
}

pub struct XOR4BitTable(HashMap<(Fr, Fr), Fr>);

impl LookUpTable for XOR4BitTable {
    // Initialise all 4 bit combinations of XOR
    fn new() -> Self {
        let mut table = XOR4BitTable(HashMap::new());

        for i in 0..=15 {
            for k in 0..=15 {
                let result = i ^ k;
                table.0.insert(
                    (Fr::from(i as u8), Fr::from(k as u8)),
                    Fr::from(result as u8),
                );
            }
        }
        table
    }

    fn borrow_map(&self) -> &HashMap<(Fr, Fr), Fr> {
        &self.0
    }
}

#[test]
fn test_size_bit_table() {
    let four_bit_table = XOR4BitTable::new();
    assert_eq!(four_bit_table.0.len(), 2usize.pow(8))
}
