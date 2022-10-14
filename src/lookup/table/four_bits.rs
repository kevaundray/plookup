use ark_bls12_381::Fr;

use super::{Generic, LookUpTable};
use std::collections::HashMap;

const BIT_RANGE: usize = 16;

/// Constructs a 4-bit Add table
pub struct Add4Bit(Generic);

impl LookUpTable for Add4Bit {
    fn borrow_map(&self) -> &HashMap<(Fr, Fr), Fr> {
        &self.0.borrow_map()
    }
}

impl Add4Bit {
    // Populate table with all 4 bit combinations of Add
    pub fn new() -> Self {
        // Define the bi-variate function for constructing the table using XOR
        let func = |a: usize, b: usize| -> Fr {
            let result = a + b;
            Fr::from(result as u128)
        };

        let add_bit_table = Generic::with_fn(func, BIT_RANGE);
        Add4Bit(add_bit_table)
    }
}

/// Constructs a 4-bit XOR table
pub struct XOR4Bit(Generic);

impl XOR4Bit {
    // Populate table with all 4 bit combinations of XOR
    pub fn new() -> Self {
        // Define the bi-variate function for constructing the table using XOR
        let func = |a: usize, b: usize| -> Fr {
            let result = a ^ b;
            Fr::from(result as u128)
        };

        let xor_bit_table = Generic::with_fn(func, BIT_RANGE);
        XOR4Bit(xor_bit_table)
    }
}
impl LookUpTable for XOR4Bit {
    fn borrow_map(&self) -> &HashMap<(Fr, Fr), Fr> {
        &self.0.borrow_map()
    }
}
