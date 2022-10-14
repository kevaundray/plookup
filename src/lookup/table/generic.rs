use crate::lookup::table::LookUpTable;
use ark_bls12_381::Fr;

use std::collections::HashMap;

/// Construct a Generic lookup table over a bi-variate function
pub struct Generic(HashMap<(Fr, Fr), Fr>);

impl Generic {
    // Initialise a table using a bi-variate function over some bit-range
    pub fn with_fn<F>(f: F, n: usize) -> Self
    where
        F: Fn(usize, usize) -> Fr,
    {
        let mut table = Generic(HashMap::new());

        for i in 0..n {
            let i_fr = Fr::from(i as u8);

            for k in 0..n {
                let k_fr = Fr::from(k as u8);

                let result = f(i, k);
                table.0.insert((i_fr, k_fr), result);
            }
        }
        table
    }
    // Initialise a table by passing all of its entries to the table
    pub fn with_hashmap(map: HashMap<(Fr, Fr), Fr>) -> Self {
        Generic(map)
    }
}

impl LookUpTable for Generic {
    fn borrow_map(&self) -> &HashMap<(Fr, Fr), Fr> {
        &self.0
    }
}
