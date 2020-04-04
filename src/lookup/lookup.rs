use super::{
    proof::LookUpProof,
    table::{LookUpTable, PreProcessedTable},
};
use crate::{multiset::MultiSet, transcript::TranscriptProtocol};
use algebra::{bls12_381::Fr, Bls12_381};
use poly_commit::kzg10::Powers;

pub struct LookUp<T: LookUpTable> {
    table: T,
    // This is the set of values which we want to prove is a subset of the
    // table values. This may or may not be equal to the whole witness.
    left_wires: MultiSet,
    right_wires: MultiSet,
    output_wires: MultiSet,
}

impl<T: LookUpTable> LookUp<T> {
    pub fn new(table: T) -> LookUp<T> {
        LookUp {
            table: table,
            left_wires: MultiSet::new(),
            right_wires: MultiSet::new(),
            output_wires: MultiSet::new(),
        }
    }
    // First reads a value from the underlying table
    // Then we add the key and value to their respective multisets
    // Returns true if the value existed in the table
    pub fn read(&mut self, key: &(Fr, Fr)) -> bool {
        let option_output = self.table.read(key);
        if option_output.is_none() {
            return false;
        }
        let output = *option_output.unwrap();

        // Add (input, output) combination into the corresponding multisets
        self.left_wires.push(key.0);
        self.right_wires.push(key.1);
        self.output_wires.push(output);

        return true;
    }

    /// Creates a proof that the values (f_1, f_2, f_3)  are within the table of values (t_1, t_2,t_3)
    pub fn prove(
        &mut self,
        proving_key: &Powers<Bls12_381>,
        preprocessed_table: &PreProcessedTable,
        transcript: &mut dyn TranscriptProtocol,
    ) -> LookUpProof {
        LookUpProof::prove(
            &self.left_wires,
            &self.right_wires,
            &self.output_wires,
            proving_key,
            preprocessed_table,
            transcript,
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::kzg10;
    use crate::lookup::table::four_bits::XOR4Bit;
    use merlin::Transcript;

    #[test]
    fn test_proof() {
        // Setup SRS
        let (proving_key, verifier_key) = kzg10::trusted_setup(2usize.pow(12), b"insecure_seed");

        // Setup Lookup with a 4 bit table
        let table = XOR4Bit::new();
        let preprocessed_table = table.preprocess(&proving_key, 2usize.pow(8));

        let mut lookup = LookUp::new(table);

        // Adds 1 XOR 2
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
        // Adds 2 XOR 4
        lookup.read(&(Fr::from(2u8), Fr::from(4u8)));
        // Adds 3 XOR 5
        lookup.read(&(Fr::from(3u8), Fr::from(5u8)));

        let mut prover_transcript = Transcript::new(b"lookup");
        let proof = lookup.prove(&proving_key, &preprocessed_table, &mut prover_transcript);

        let mut verifier_transcript = Transcript::new(b"lookup");
        let ok = proof.verify(&verifier_key, &preprocessed_table, &mut verifier_transcript);
        assert!(ok);
    }
}
