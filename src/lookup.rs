use crate::kzg10;
use crate::lookup_table::LookUpTable;
use crate::multiset::MultiSet;
use crate::multiset_equality;
use crate::proof::{MultiSetEqualityProof, OpeningProof};
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

    // Pads the witness or table, so that len(table) = len(witness) + 1
    // Due to FFTs, we apply additional padding so that table size is a power of two
    fn pad(&self, witness: &mut MultiSet, table: &mut MultiSet) {
        if witness.len() < table.len() {
            // First pad the table multiset to next power of two
            let pad_to_power_of_two = table.len().next_power_of_two() - table.len();
            table.extend(pad_to_power_of_two, table.last());

            // Then pad the witness to be one less than this power of two
            let pad_amount = table.len() - witness.len() - 1;
            witness.extend(pad_amount, witness.last())
        } else {
            // First pad the witness to be one less than a power of two
            let pad_to_power_of_two = witness.len().next_power_of_two() - witness.len();
            witness.extend(pad_to_power_of_two - 1, table.last());

            // Then pad the table to be one more than the witness
            let pad_amount = witness.len() + 1;
            table.extend(pad_amount, table.last())
        }
    }

    /// Aggregates the table and witness values into one multiset
    /// sorts, and pads the witness and or table to be the correct size
    pub fn to_multiset(&self, alpha: Fr) -> (MultiSet, MultiSet) {
        // Compute challenge alpha^0, alpha^1, alpha^2
        let one = Fr::from(1u8);
        let alpha_sq = alpha * alpha;

        // First get the witness as multisets
        let left = &self.left_wires;
        let right = &self.right_wires;
        let output = &self.output_wires;

        // Now lets get the table values as multisets
        let (t_left, t_right, t_output) = self.table.to_multiset();

        // Now we need to aggregate our witness values into one multiset
        let mut merged_witness = MultiSet::aggregate((left, right, output), (one, alpha, alpha_sq));
        // Now we need to aggregate our table values into one multiset
        let mut merged_table =
            MultiSet::aggregate((&t_left, &t_right, &t_output), (one, alpha, alpha_sq));
        // Sort merged table values
        merged_table = merged_table.sort();

        // Pad values
        self.pad(&mut merged_witness, &mut merged_table);
        (merged_witness, merged_table)
    }

    /// Creates a proof that the multiset is within the table
    fn prove(
        &self,
        proving_key: &Powers<Bls12_381>,
        transcript: &mut dyn TranscriptProtocol,
    ) -> MultiSetEqualityProof {
        // First we convert the table to a multiset and apply appropriate padding
        let alpha = transcript.challenge_scalar(b"alpha");
        let (f, t) = self.to_multiset(alpha);

        assert_eq!(f.len() + 1, t.len());

        let domain: EvaluationDomain<Fr> = EvaluationDomain::new(f.len()).unwrap();

        // Convert witness and table to polynomials
        let f_poly = f.to_polynomial(&domain);
        let t_poly = t.to_polynomial(&domain);

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

        // Compute the witness
        //
        let evaluation_challenge = transcript.challenge_scalar(b"evaluation_challenge");
        let evaluation_omega = evaluation_challenge * domain.group_gen;

        let h_1_eval = h_1_poly.evaluate(evaluation_challenge);
        let h_1_omega_eval = h_1_poly.evaluate(evaluation_omega);
        let h_2_eval = h_2_poly.evaluate(evaluation_challenge);
        let h_2_omega_eval = h_2_poly.evaluate(evaluation_omega);
        let z_eval = z_poly.evaluate(evaluation_challenge);
        let z_omega_eval = z_poly.evaluate(evaluation_omega);
        let q_eval = quotient_poly.evaluate(evaluation_challenge);

        transcript.append_scalar(b"h_1_eval", &h_1_eval);
        transcript.append_scalar(b"h_1_omega_eval", &h_1_omega_eval);
        transcript.append_scalar(b"h_2_eval", &h_2_eval);
        transcript.append_scalar(b"h_2_omega_eval", &h_2_omega_eval);
        transcript.append_scalar(b"z_eval", &z_eval);
        transcript.append_scalar(b"z_omega_eval", &z_eval);
        transcript.append_scalar(b"q_eval", &q_eval);

        // Compute opening proofs for f(X) evaluated at `z`
        let h_1_witness = kzg10::compute_witness(&h_1_poly, evaluation_challenge);
        let h_1_witness_comm = kzg10::commit(proving_key, &h_1_witness);

        let h_2_witness = kzg10::compute_witness(&h_2_poly, evaluation_challenge);
        let h_2_witness_comm = kzg10::commit(proving_key, &h_2_witness);

        let z_witness = kzg10::compute_witness(&z_poly, evaluation_challenge);
        let z_witness_comm = kzg10::commit(proving_key, &z_witness);

        let q_witness = kzg10::compute_witness(&quotient_poly, evaluation_challenge);
        let q_witness_comm = kzg10::commit(proving_key, &q_witness);

        // Compute opening proofs for f(X) evaluated at `z * omega`
        let h_1_omega_witness = kzg10::compute_witness(&h_1_poly, evaluation_omega);
        let h_1_omega_witness_comm = kzg10::commit(proving_key, &h_1_omega_witness);

        let h_2_omega_witness = kzg10::compute_witness(&h_2_poly, evaluation_omega);
        let h_2_omega_witness_comm = kzg10::commit(proving_key, &h_2_omega_witness);

        let z_omega_witness = kzg10::compute_witness(&z_poly, evaluation_omega);
        let z_omega_witness_comm = kzg10::commit(proving_key, &z_omega_witness);

        MultiSetEqualityProof {
            h_1_proof: OpeningProof::new((Some(h_1_commit), h_1_witness_comm, h_1_eval)),
            h_2_proof: OpeningProof::new((Some(h_2_commit), h_2_witness_comm, h_2_eval)),
            z_proof: OpeningProof::new((Some(z_commit), z_witness_comm, z_eval)),
            q_proof: OpeningProof::new((Some(q_commit), q_witness_comm, q_eval)),
            h_1_omega_proof: OpeningProof::new((None, h_1_omega_witness_comm, h_1_omega_eval)),
            h_2_omega_proof: OpeningProof::new((None, h_2_omega_witness_comm, h_2_omega_eval)),
            z_omega_proof: OpeningProof::new((None, z_omega_witness_comm, z_omega_eval)),
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
        let table = XOR4BitTable::new();

        // Setup lookup and add 3 XOR reads into it
        let mut lookup = LookUp::new(table);

        // Add 1 XOR 2
        lookup.read(&(Fr::from(2u8), Fr::from(2u8)));
        // Add 2 XOR 4
        lookup.read(&(Fr::from(3u8), Fr::from(2u8)));
        // Add 3 XOR 5
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));

        let (f, t) = lookup.to_multiset(Fr::from(5u8));
        assert_eq!(f.len() + 1, t.len());

        assert!(t.len().is_power_of_two());
    }

    #[test]
    fn test_inclusion() {
        let table = XOR4BitTable::new();

        let mut lookup = LookUp::new(table);

        // Add 2 XOR 2
        lookup.read(&(Fr::from(2u8), Fr::from(2u8)));
        // Add 1 XOR 2
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
        // Add 3 XOR 5
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
        let (f, t) = lookup.to_multiset(Fr::from(5u8));
        assert!(f.is_subset_of(&t));
    }
    #[test]
    fn test_len() {
        // Check that the correct values are being added to the witness
        // If the value is not in the XOR4BitTable, it is not added to the witness
        // For a 4-bit XOR table the range is [0,15]

        let table = XOR4BitTable::new();
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

        let (f, t) = lookup.to_multiset(Fr::from(5u8));
        assert!(f.is_subset_of(&t));
    }
    #[test]
    fn test_proof() {
        // Setup SRS
        let universal_parameters = kzg10::trusted_setup(2usize.pow(12), b"insecure_seed");
        let (proving_key, verifier_key) = kzg10::trim(&universal_parameters, 2usize.pow(12));

        // Setup Lookup with a 4 bit table
        let table = XOR4BitTable::new();
        let mut lookup = LookUp::new(table);

        // Adds 1 XOR 2
        lookup.read(&(Fr::from(1u8), Fr::from(2u8)));
        // Adds 2 XOR 4
        lookup.read(&(Fr::from(2u8), Fr::from(4u8)));
        // Adds 3 XOR 5
        lookup.read(&(Fr::from(3u8), Fr::from(5u8)));

        let mut prover_transcript = Transcript::new(b"lookup");
        let proof = lookup.prove(&proving_key, &mut prover_transcript);
    }
}
