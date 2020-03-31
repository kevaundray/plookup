use crate::kzg10;
use crate::lookup_table::{LookUpTable, PreProcessedTable};
use crate::multiset::MultiSet;
use crate::multiset_equality;
use crate::proof::{Commitments, Evaluations, MultiSetEqualityProof};
use crate::quotient_poly;
use crate::transcript::TranscriptProtocol;
use algebra::bls12_381::Fr;
use algebra::Bls12_381;
use ff_fft::{DensePolynomial as Polynomial, EvaluationDomain};
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

    /// Aggregates the table and witness values into one multiset
    /// sorts, and pads the witness and or table to be the correct size
    pub fn to_multiset(
        &self,
        preprocessed_table: &PreProcessedTable,
        alpha: Fr,
    ) -> (MultiSet, MultiSet) {
        // Now we need to aggregate our witness values into one multiset
        let mut merged_witness = MultiSet::aggregate(
            vec![&self.left_wires, &self.right_wires, &self.output_wires],
            alpha,
        );

        // Now we need to aggregate our table values into one multiset
        let mut merged_table = MultiSet::aggregate(
            vec![
                &preprocessed_table.t_1.0,
                &preprocessed_table.t_2.0,
                &preprocessed_table.t_3.0,
            ],
            alpha,
        );
        // Sort merged table values
        merged_table = merged_table.sort();

        // Pad witness to be one less than `n`
        assert!(merged_witness.len() < merged_table.len()); // XXX: We could incorporate this in the API by counting the number of reads
        let pad_by = preprocessed_table.n - 1 - merged_witness.len();
        merged_witness.extend(pad_by, merged_witness.last());

        (merged_witness, merged_table)
    }

    /// Creates a proof that the multiset is within the table
    pub fn prove(
        &self,
        proving_key: &Powers<Bls12_381>,
        preprocessed_table: &PreProcessedTable,
        transcript: &mut dyn TranscriptProtocol,
    ) -> MultiSetEqualityProof {
        // Generate alpha challenge
        let alpha = transcript.challenge_scalar(b"alpha");
        transcript.append_scalar(b"alpha", &alpha);

        // Aggregate witness and table values using a random challenge
        let (f, t) = self.to_multiset(preprocessed_table, alpha);
        assert_eq!(f.len() + 1, t.len());

        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(f.len()).unwrap();

        // Convert witness and table to polynomials
        let f_poly = f.to_polynomial(&domain);
        let f_commit = kzg10::commit(proving_key, &f_poly);

        let t_poly = t.to_polynomial(&domain);
        let t_commit = kzg10::commit(proving_key, &t_poly);

        // Compute h_1 and h_2
        let (h_1, h_2) = multiset_equality::compute_h1_h2(&f, &t);

        // Convert h_1 and h_2 to polynomials
        let h_1_poly = h_1.to_polynomial(&domain);
        let h_2_poly = h_2.to_polynomial(&domain);

        // Commit to h_1(X) and h_2(X)
        let h_1_commit = kzg10::commit(proving_key, &h_1_poly);
        let h_2_commit = kzg10::commit(proving_key, &h_2_poly);

        // Add commitments to transcript
        transcript.append_commitment(b"h_1_poly", &h_1_commit);
        transcript.append_commitment(b"h_2_poly", &h_2_commit);

        let beta = transcript.challenge_scalar(b"beta");
        let gamma = transcript.challenge_scalar(b"gamma");

        // Compute Z(X)
        let z_evaluations =
            multiset_equality::compute_accumulator_values(&f, &t, &h_1, &h_2, beta, gamma);
        let z_poly = Polynomial::from_coefficients_vec(domain.ifft(&z_evaluations));

        // Commit to Z(X)
        let z_commit = kzg10::commit(proving_key, &z_poly);
        transcript.append_commitment(b"accumulator_poly", &z_commit);

        // Compute quotient polynomial
        let (quotient_poly, _) = quotient_poly::compute(
            &domain, &z_poly, &f_poly, &t_poly, &h_1_poly, &h_2_poly, beta, gamma,
        );

        // Commit to quotient polynomial
        let q_commit = kzg10::commit(proving_key, &quotient_poly);
        transcript.append_commitment(b"quotient_poly", &q_commit);

        // Compute the Witness that f was a subset of t
        //
        let evaluation_challenge = transcript.challenge_scalar(b"evaluation_challenge");
        transcript.append_scalar(b"evaluation_challenge", &evaluation_challenge);
        let evaluation_omega = evaluation_challenge * domain.group_gen;

        // Compute evaluations at `z`
        let f_eval = f_poly.evaluate(evaluation_challenge);
        let t_eval = t_poly.evaluate(evaluation_challenge);
        let h_1_eval = h_1_poly.evaluate(evaluation_challenge);
        let h_2_eval = h_2_poly.evaluate(evaluation_challenge);
        let z_eval = z_poly.evaluate(evaluation_challenge);
        let q_eval = quotient_poly.evaluate(evaluation_challenge);

        // Compute evaluations at `z * omega`
        let t_omega_eval = t_poly.evaluate(evaluation_omega);
        let h_1_omega_eval = h_1_poly.evaluate(evaluation_omega);
        let h_2_omega_eval = h_2_poly.evaluate(evaluation_omega);
        let z_omega_eval = z_poly.evaluate(evaluation_omega);

        transcript.append_scalar(b"f_eval", &f_eval);
        transcript.append_scalar(b"t_eval", &t_eval);
        transcript.append_scalar(b"h_1_eval", &h_1_eval);
        transcript.append_scalar(b"h_2_eval", &h_2_eval);
        transcript.append_scalar(b"z_eval", &z_eval);
        transcript.append_scalar(b"q_eval", &q_eval);
        transcript.append_scalar(b"t_omega_eval", &t_omega_eval);
        transcript.append_scalar(b"h_1_omega_eval", &h_1_omega_eval);
        transcript.append_scalar(b"h_2_omega_eval", &h_2_omega_eval);
        transcript.append_scalar(b"z_omega_eval", &z_omega_eval);

        let aggregation_challenge = transcript.challenge_scalar(b"witness_aggregation");

        // Compute opening proof for f(X) evaluated at `z`
        let agg_witness = kzg10::compute_aggregate_witness(
            vec![
                &f_poly,
                &t_poly,
                &h_1_poly,
                &h_2_poly,
                &z_poly,
                &quotient_poly,
            ],
            evaluation_challenge,
            aggregation_challenge,
        );
        let agg_witness_comm = kzg10::commit(proving_key, &agg_witness);

        // Compute opening proofs for f(X) evaluated at `z * omega`
        let shifted_agg_witness = kzg10::compute_aggregate_witness(
            vec![&t_poly, &h_1_poly, &h_2_poly, &z_poly],
            evaluation_omega,
            aggregation_challenge,
        );
        let shifted_agg_witness_comm = kzg10::commit(proving_key, &shifted_agg_witness);

        MultiSetEqualityProof {
            n: preprocessed_table.n,
            evaluations: Evaluations {
                f: f_eval,
                t: t_eval,
                t_omega: t_omega_eval,
                h_1: h_1_eval,
                h_1_omega: h_1_omega_eval,
                h_2: h_2_eval,
                h_2_omega: h_2_omega_eval,
                z: z_eval,
                z_omega: z_omega_eval,
            },
            commitments: Commitments {
                f: f_commit,
                q: q_commit,
                t: t_commit,
                h_1: h_1_commit,
                h_2: h_2_commit,
                z: z_commit,
            },
            aggregate_witness_comm: agg_witness_comm,
            shifted_aggregate_witness_comm: shifted_agg_witness_comm,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::lookup_table::XOR4BitTable;
    use merlin::Transcript;

    #[test]
    fn test_pad_correct() {
        // Setup SRS
        let universal_parameters = kzg10::trusted_setup(2usize.pow(12), b"insecure_seed");
        let (proving_key, _) = kzg10::trim(&universal_parameters, 2usize.pow(12));

        let table = XOR4BitTable::new();
        let preprocessed_table = table.preprocess(&proving_key, 2usize.pow(8));

        // Setup lookup and add 3 XOR reads into it
        let mut lookup = LookUp::new(table);

        // Add 1 XOR 2
        lookup.read(&(Fr::from(2u8), Fr::from(2u8)));
        // Add 2 XOR 4
        lookup.read(&(Fr::from(3u8), Fr::from(2u8)));
        // Add 3 XOR 5
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));

        let (f, t) = lookup.to_multiset(&preprocessed_table, Fr::from(5u8));
        assert_eq!(f.len() + 1, t.len());

        assert!(t.len().is_power_of_two());
    }

    #[test]
    fn test_inclusion() {
        // Setup SRS
        let universal_parameters = kzg10::trusted_setup(2usize.pow(12), b"insecure_seed");
        let (proving_key, _) = kzg10::trim(&universal_parameters, 2usize.pow(12));

        let table = XOR4BitTable::new();
        let preprocessed_table = table.preprocess(&proving_key, 2usize.pow(8));

        let mut lookup = LookUp::new(table);

        // Add 2 XOR 2
        lookup.read(&(Fr::from(2u8), Fr::from(2u8)));
        // Add 1 XOR 2
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
        // Add 3 XOR 5
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));

        let (f, t) = lookup.to_multiset(&preprocessed_table, Fr::from(5u8));
        assert!(f.is_subset_of(&t));
    }
    #[test]
    fn test_len() {
        // Check that the correct values are being added to the witness
        // If the value is not in the XOR4BitTable, it is not added to the witness
        // For a 4-bit XOR table the range is [0,15]

        // Setup SRS
        let universal_parameters = kzg10::trusted_setup(2usize.pow(12), b"insecure_seed");
        let (proving_key, _) = kzg10::trim(&universal_parameters, 2usize.pow(12));

        let table = XOR4BitTable::new();
        let preprocessed_table = table.preprocess(&proving_key, 2usize.pow(8));

        let mut lookup = LookUp::new(table);

        let added = lookup.read(&(Fr::from(16u8), Fr::from(6u8)));
        assert!(!added);

        let added = lookup.read(&(Fr::from(8u8), Fr::from(17u8)));
        assert!(!added);
        let added = lookup.read(&(Fr::from(15u8), Fr::from(13u8)));
        assert!(added);

        assert_eq!(lookup.left_wires.len(), 1);
        assert_eq!(lookup.right_wires.len(), 1);
        assert_eq!(lookup.output_wires.len(), 1);

        let (f, t) = lookup.to_multiset(&preprocessed_table, Fr::from(5u8));
        assert!(f.is_subset_of(&t));
    }
    #[test]
    fn test_proof() {
        // Setup SRS
        let universal_parameters = kzg10::trusted_setup(2usize.pow(12), b"insecure_seed");
        let (proving_key, verifier_key) = kzg10::trim(&universal_parameters, 2usize.pow(12));

        // Setup Lookup with a 4 bit table
        let table = XOR4BitTable::new();
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
