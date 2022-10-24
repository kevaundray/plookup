use ark_bls12_381::Bls12_381;
use ark_poly_commit::kzg10::{Powers, VerifierKey};

use super::table::PreProcessedTable;
use crate::{
    kzg10,
    multiset::{EqualityProof, MultiSet},
    transcript::TranscriptProtocol,
};

pub struct LookUpProof {
    pub multiset_equality_proof: EqualityProof,
}

impl LookUpProof {
    pub fn prove(
        f_1: &MultiSet,
        f_2: &MultiSet,
        f_3: &MultiSet,
        proving_key: &Powers<Bls12_381>,
        preprocessed_table: &PreProcessedTable,
        transcript: &mut dyn TranscriptProtocol,
    ) -> LookUpProof {
        // Generate alpha challenge
        let alpha = transcript.challenge_scalar(b"alpha");
        transcript.append_scalar(b"alpha", &alpha);

        // Aggregates the table and witness values into one multiset
        // and pads the witness to be the correct size
        //
        // Aggregate our table values into one multiset
        let merged_table = MultiSet::aggregate(
            vec![
                &preprocessed_table.t_1.0,
                &preprocessed_table.t_2.0,
                &preprocessed_table.t_3.0,
            ],
            alpha,
        );

        // Aggregate witness values into one multiset
        let mut merged_witness = MultiSet::aggregate(vec![f_1, f_2, f_3], alpha);

        // Pad merged Witness to be one less than `n`
        assert!(merged_witness.len() < preprocessed_table.n);
        let pad_by = preprocessed_table.n - 1 - merged_witness.len();
        merged_witness.extend(pad_by, merged_witness.last());

        // Create a Multi-set equality proof
        let multiset_equality_proof =
            EqualityProof::prove(merged_witness, merged_table, proving_key, transcript);

        LookUpProof {
            multiset_equality_proof,
        }
    }
    pub fn verify(
        &self,
        verification_key: &VerifierKey<Bls12_381>,
        preprocessed_table: &PreProcessedTable,
        transcript: &mut dyn TranscriptProtocol,
    ) -> bool {
        // Merge preprocessed commitments to table using `alpha` challenge
        let alpha = transcript.challenge_scalar(b"alpha");
        let merged_table_commit = kzg10::aggregate_commitments(
            vec![
                &preprocessed_table.t_1.1,
                &preprocessed_table.t_2.1,
                &preprocessed_table.t_3.1,
            ],
            alpha,
        );
        transcript.append_scalar(b"alpha", &alpha);

        // Call Multiset Equality Proof as a sub-routine
        self.multiset_equality_proof.verify(
            preprocessed_table.n,
            verification_key,
            merged_table_commit,
            transcript,
        )
    }
}
