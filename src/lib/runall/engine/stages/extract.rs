//! Extract stage: `TemplateBatch` → serialized unmapped BAM records with
//! UMI/barcode tags.
//!
//! Parallel pool stage. Wraps the existing
//! [`crate::extract::make_raw_records`] function. Each worker has its
//! own [`ExtractStage`] instance; [`ExtractParams`] and the read
//! structures are `Arc`-shared (read-only), so cloning is cheap.
//!
//! ## Input
//!
//! [`TemplateBatch`] — a batch of FASTQ templates, each with one
//! `OwnedFastqRecord` per input stream and a stable `ordinal` assigned
//! by the upstream pairing/merge stage.
//!
//! ## Output
//!
//! [`SerializedBatch`] — concatenated length-prefixed unmapped BAM
//! record bytes (4-byte little-endian `block_size` + payload per
//! record) and the record count, with `ordinal` preserved unchanged
//! from input. Ready for `CorrectStage`, `AlignAndMerge`, or
//! `BgzfCompress` → `BamFileWrite`.
//!
//! ## Ordering guarantees
//!
//! Record order within the batch follows template order in the input.
//! Cross-batch ordering is the upstream `ordinal` sequence, preserved
//! verbatim.
//!
//! ## Memory model
//!
//! 1:1 with input: one output batch per input batch. Per-batch `Vec`
//! grows with the number of records; no persistent scratch buffers.
//!
//! ## Determinism
//!
//! Byte-identical across runs for a given input: UMI/barcode extraction
//! is a pure function of `(template bytes, ReadStructure, ExtractParams,
//! QualityEncoding)`.

use std::sync::Arc;

use anyhow::Result;
use read_structure::ReadStructure;

use crate::extract::{ExtractParams, QualityEncoding, make_raw_records};
use crate::fastq::FastqSet;
use crate::runall::engine::fastq_types::{FastqTemplateV2, TemplateBatch};
use crate::runall::engine::output_types::{RawBytes, SerializedBatch};
use crate::runall::engine::stage::{Parallelism, Stage};

/// Parallel stage that extracts UMIs from FASTQ templates and emits raw
/// unmapped BAM record bytes with UMI/barcode tags applied.
pub struct ExtractStage {
    params: Arc<ExtractParams>,
    read_structures: Arc<Vec<ReadStructure>>,
    encoding: QualityEncoding,
}

impl ExtractStage {
    /// Construct a new [`ExtractStage`].
    ///
    /// `read_structures` must have one entry per input FASTQ stream, in the
    /// same order as the per-stream records carried in each [`FastqTemplateV2`].
    #[must_use]
    pub fn new(
        params: Arc<ExtractParams>,
        read_structures: Arc<Vec<ReadStructure>>,
        encoding: QualityEncoding,
    ) -> Self {
        Self { params, read_structures, encoding }
    }
}

impl Clone for ExtractStage {
    fn clone(&self) -> Self {
        Self {
            params: self.params.clone(),
            read_structures: self.read_structures.clone(),
            encoding: self.encoding,
        }
    }
}

impl Stage for ExtractStage {
    type Input = TemplateBatch;
    type Output = SerializedBatch;

    #[tracing::instrument(name = "extract", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let TemplateBatch { templates, ordinal } = input;
        let mut out_data: Vec<u8> = Vec::new();
        let mut count: u64 = 0;
        for template in &templates {
            let read_set = template_to_fastq_set(template, &self.read_structures)?;
            let extracted = make_raw_records(&read_set, self.encoding, &self.params)?;
            out_data.extend_from_slice(&extracted.data);
            count += extracted.num_records;
        }
        out(SerializedBatch {
            primary: RawBytes { data: out_data, record_count: count },
            secondary: None,
            ordinal,
        });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, out: &Self::Output) -> usize {
        out.primary.data.capacity() + out.secondary.as_ref().map_or(0, |s| s.data.capacity())
    }

    fn name(&self) -> &'static str {
        "ExtractStage"
    }
}

/// Adapter: convert a [`FastqTemplateV2`] (one raw [`crate::runall::engine::fastq_types::OwnedFastqRecord`] per input stream)
/// into the [`FastqSet`] expected by [`make_raw_records`].
///
/// For each record in the template, applies the corresponding [`ReadStructure`]
/// via [`FastqSet::from_record_with_structure`] to split the record into
/// segments (template / UMI / barcode), then merges all per-stream sets into
/// one via [`FastqSet::combine_readsets`]. This mirrors the v1 extract flow in
/// `src/lib/commands/runall/run_extract.rs` and `src/lib/commands/extract.rs`.
fn template_to_fastq_set(
    template: &FastqTemplateV2,
    read_structures: &[ReadStructure],
) -> Result<FastqSet> {
    anyhow::ensure!(
        template.records.len() == read_structures.len(),
        "ExtractStage: template has {} records but {} read structures were configured",
        template.records.len(),
        read_structures.len(),
    );

    let mut fastq_sets: Vec<FastqSet> = Vec::with_capacity(template.records.len());
    for (record, rs) in template.records.iter().zip(read_structures.iter()) {
        let fastq_set = FastqSet::from_record_with_structure(
            &record.name,
            &record.sequence,
            &record.quality,
            rs,
            &[], // no skip reasons — mismatched lengths are hard errors here
        )?;
        fastq_sets.push(fastq_set);
    }

    Ok(FastqSet::combine_readsets(fastq_sets))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    use crate::runall::engine::fastq_types::{OwnedFastqRecord, TemplateBatch};

    /// Minimal `ExtractParams` suitable for unit tests. Mirrors the `test_params`
    /// helper in `src/lib/extract/mod.rs`.
    fn minimal_extract_params() -> ExtractParams {
        ExtractParams {
            read_group_id: "A".to_string(),
            umi_tag: *b"RX",
            cell_tag: *b"CB",
            umi_qual_tag: None,
            cell_qual_tag: None,
            single_tag: None,
            annotate_read_names: false,
            extract_umis_from_read_names: false,
            store_sample_barcode_qualities: false,
        }
    }

    fn minimal_read_structures() -> Vec<ReadStructure> {
        vec![
            ReadStructure::from_str("4M+T").expect("valid read structure"),
            ReadStructure::from_str("+T").expect("valid read structure"),
        ]
    }

    #[test]
    fn test_extract_stage_cloning_shares_params_and_structures() {
        let params = Arc::new(minimal_extract_params());
        let structures = Arc::new(minimal_read_structures());
        let stage =
            ExtractStage::new(params.clone(), structures.clone(), QualityEncoding::Standard);
        let cloned = stage.clone();

        // Arc-shared state must remain shared across clones (cheap, read-only).
        assert!(Arc::ptr_eq(&stage.params, &cloned.params));
        assert!(Arc::ptr_eq(&stage.read_structures, &cloned.read_structures));
        assert_eq!(stage.encoding, cloned.encoding);
        assert_eq!(stage.name(), "ExtractStage");
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
    }

    #[test]
    fn test_extract_stage_processes_paired_template() {
        let params = Arc::new(minimal_extract_params());
        let structures = Arc::new(minimal_read_structures());
        let mut stage = ExtractStage::new(params, structures, QualityEncoding::Standard);

        // R1 has 4-base UMI (ACGT) + 8-base template; R2 is 8-base template only.
        let r1 = OwnedFastqRecord {
            name: b"read1".to_vec(),
            sequence: b"ACGTAAAAAAAA".to_vec(),
            quality: b"IIIIIIIIIIII".to_vec(),
        };
        let r2 = OwnedFastqRecord {
            name: b"read1".to_vec(),
            sequence: b"CCCCCCCC".to_vec(),
            quality: b"IIIIIIII".to_vec(),
        };
        let template = FastqTemplateV2 { records: vec![r1, r2], name: b"read1".to_vec() };
        let batch = TemplateBatch { templates: vec![template], ordinal: 7 };

        let mut captured = None;
        stage
            .process(batch, &mut |v| {
                captured = Some(v);
            })
            .expect("extract stage should succeed");
        let out = captured.expect("stage must emit");

        // Two template segments (one per read) produce two BAM records.
        assert_eq!(out.primary.record_count, 2, "expected two BAM records for paired template");
        assert!(!out.primary.data.is_empty(), "BAM payload must not be empty");
        assert!(out.secondary.is_none());
        assert_eq!(out.ordinal, 7, "ordinal must be preserved through ExtractStage");
        assert!(stage.output_memory_estimate(&out) >= out.primary.data.len());
    }

    #[test]
    fn test_extract_stage_rejects_template_stream_count_mismatch() {
        let params = Arc::new(minimal_extract_params());
        let structures = Arc::new(minimal_read_structures()); // expects 2 streams
        let mut stage = ExtractStage::new(params, structures, QualityEncoding::Standard);

        // Only one record per template: deliberate mismatch.
        let r1 = OwnedFastqRecord {
            name: b"readX".to_vec(),
            sequence: b"ACGTAAAAAAAA".to_vec(),
            quality: b"IIIIIIIIIIII".to_vec(),
        };
        let template = FastqTemplateV2 { records: vec![r1], name: b"readX".to_vec() };
        let batch = TemplateBatch { templates: vec![template], ordinal: 0 };

        let err = stage
            .process(batch, &mut |_| panic!("no output on error"))
            .expect_err("mismatch must error");
        let msg = err.to_string();
        assert!(
            msg.contains("template has 1 records but 2 read structures"),
            "unexpected error message: {msg}"
        );
    }
}
