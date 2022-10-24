extern crate plookup;

use ark_bls12_381::Fr;
use merlin::Transcript;
use plookup::kzg10::trusted_setup;
use plookup::lookup::{
    lookup::LookUp, table::four_bits::XOR4Bit, table::Generic, table::LookUpTable,
};
use std::collections::HashMap;

#[test]
fn test_xor_four_bit_lookup() {
    let (prover_key, verifier_key) = trusted_setup(2usize.pow(10));

    // Setup Lookup with a 4 bit table
    let table = XOR4Bit::new();
    let preprocessed_table = table.preprocess(&prover_key, 2usize.pow(8));

    let mut lookup = LookUp::new(table);

    // Adds 1 XOR 2
    lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
    // Adds 2 XOR 4
    lookup.read(&(Fr::from(2u8), Fr::from(4u8)));
    // Adds 3 XOR 5
    lookup.read(&(Fr::from(3u8), Fr::from(5u8)));

    let mut prover_transcript = Transcript::new(b"lookup");
    let proof = lookup.prove(&prover_key, &preprocessed_table, &mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"lookup");
    let ok = proof.verify(&verifier_key, &preprocessed_table, &mut verifier_transcript);
    assert!(ok);
}

#[test]
fn test_inefficient_range_lookups() {
    let (prover_key, verifier_key) = trusted_setup(2usize.pow(10));

    // Setup a 8-bit range protocol using Lookup API
    let mut hashmap: HashMap<(Fr, Fr), Fr> = HashMap::new();
    for i in 0..2usize.pow(8) {
        hashmap.insert(
            (Fr::from(i as u128), Fr::from(i as u128)),
            Fr::from(i as u128),
        );
    }
    let table = Generic::with_hashmap(hashmap);
    let preprocessed_table = table.preprocess(&prover_key, 2usize.pow(8));

    let mut lookup = LookUp::new(table);

    // Adds 1 to the rangeproof
    lookup.read(&(Fr::from(1u8), Fr::from(1u8)));
    // Adds 2 to the rangeproof
    lookup.read(&(Fr::from(2u8), Fr::from(2u8)));
    // Adds 10 to the rangeproof
    lookup.read(&(Fr::from(10u8), Fr::from(10u8)));
    // Adds 200 to the rangeproof
    lookup.read(&(Fr::from(200u8), Fr::from(200u8)));
    // Adds 30000 to the rangeproof
    lookup.read(&(Fr::from(30000u128), Fr::from(30000u128)));

    let mut prover_transcript = Transcript::new(b"range_lookup");
    let proof = lookup.prove(&prover_key, &preprocessed_table, &mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"range_lookup");
    let ok = proof.verify(&verifier_key, &preprocessed_table, &mut verifier_transcript);
    assert!(ok);
}

#[test]
fn test_mul_four_bit_lookup() {
    let (prover_key, verifier_key) = trusted_setup(2usize.pow(10));

    // 4-bit multiplication lookup table
    let table = Generic::with_fn(
        |a: usize, b: usize| -> Fr {
            let a_fr = Fr::from(a as u128);
            let b_fr = Fr::from(b as u128);
            a_fr * b_fr
        },
        16,
    );
    let preprocessed_table = table.preprocess(&prover_key, 2usize.pow(8));

    let mut lookup = LookUp::new(table);

    // Mul 1 X 2
    lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
    // Adds 2 X 4
    lookup.read(&(Fr::from(2u8), Fr::from(4u8)));
    // Adds 3 X 5
    lookup.read(&(Fr::from(3u8), Fr::from(5u8)));
    // Adds 10 X 3
    lookup.read(&(Fr::from(10u8), Fr::from(3u8)));

    let mut prover_transcript = Transcript::new(b"lookup");
    let proof = lookup.prove(&prover_key, &preprocessed_table, &mut prover_transcript);

    let mut verifier_transcript = Transcript::new(b"lookup");
    let ok = proof.verify(&verifier_key, &preprocessed_table, &mut verifier_transcript);
    assert!(ok);
}
