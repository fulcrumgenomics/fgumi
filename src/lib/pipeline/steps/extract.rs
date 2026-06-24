//! `ExtractStep` — parallel FASTQ-to-BAM conversion step.
//!
//! Transforms a [`FastqTemplateBatch`] (N `FastqTemplate`s) into a
//! [`BamTemplateBatch`] (N `Template`s) by:
//!
//! 1. Applying read structures to each `FastqRecord` in the template,
//!    yielding a combined [`FastqSet`](crate::fastq::FastqSet).
//! 2. Calling [`make_raw_records_from_fastq_set`] to produce
//!    `Vec<RawRecord>`.
//! 3. Building a [`Template`](crate::template::Template) via
//!    [`Template::from_records`].
//!
//! The step is `Parallel` (closure-driven via [`ProcessOrdered`]) and
//! preserves batch ordering via `BamTemplateBatch::batch_serial`.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use read_structure::ReadStructure;

use crate::commands::extract::{ExtractOptions, make_raw_records_from_fastq_set};
use crate::fastq::FastqSet;
use crate::pipeline::steps::process::ProcessOrdered;
use crate::pipeline::steps::source::zip_fastq::FastqTemplateBatch;
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::template::Template;

/// Build a `Parallel` step that converts [`FastqTemplateBatch`] into
/// [`BamTemplateBatch`].
///
/// Each input batch's `FastqTemplate`s are independently converted: read
/// structures are applied, UMI/barcode tags are extracted via
/// [`make_raw_records_from_fastq_set`], and the resulting `RawRecord`s are
/// assembled into `Template`s.
///
/// The returned step preserves batch ordering (output ordinal ==
/// `batch_serial` from the input).
///
/// # Arguments
///
/// * `read_structures` — one per FASTQ input file, shared across workers.
/// * `extract_opts` — tag-output and name-annotation options, shared across
///   workers.
/// * `records_emitted` — running counter incremented with the number of BAM
///   records (reads) emitted per batch.
/// * `output_byte_limit` — byte-bounded queue limit for the output branch.
pub fn build_extract_step(
    read_structures: Arc<Vec<ReadStructure>>,
    extract_opts: Arc<ExtractOptions>,
    records_emitted: Arc<AtomicU64>,
    output_byte_limit: u64,
) -> ProcessOrdered<
    FastqTemplateBatch,
    BamTemplateBatch,
    impl Fn(FastqTemplateBatch) -> io::Result<BamTemplateBatch> + Send + Sync + 'static,
> {
    ProcessOrdered::new(
        "ExtractStep",
        output_byte_limit,
        move |batch: FastqTemplateBatch| -> io::Result<BamTemplateBatch> {
            extract_batch(&read_structures, &extract_opts, &records_emitted, &batch)
        },
    )
}

/// Convert one [`FastqTemplateBatch`] into a [`BamTemplateBatch`], applying read
/// structures, extracting tags, and assembling `Template`s. This is the body of
/// the [`build_extract_step`] closure, factored out so it can be unit-tested
/// directly (the closure inside `ProcessOrdered` is otherwise unreachable).
///
/// `records_emitted` is incremented by the number of BAM records (reads)
/// emitted — summed across templates, not the template count.
///
/// # Errors
///
/// Returns `io::ErrorKind::InvalidData` if a template's record count does not
/// equal the read-structure count (an internal-invariant violation — see the
/// per-template guard below), or if read-structure application / tag extraction
/// / template assembly fails for any record.
pub(crate) fn extract_batch(
    read_structures: &[ReadStructure],
    extract_opts: &ExtractOptions,
    records_emitted: &AtomicU64,
    batch: &FastqTemplateBatch,
) -> io::Result<BamTemplateBatch> {
    let serial = batch.batch_serial;
    let mut templates = Vec::with_capacity(batch.templates.len());

    for fq_template in &batch.templates {
        // Convert each FastqRecord into a FastqSet via its read structure. The
        // per-template record count must equal the read-structure count: this
        // is an internal invariant enforced upstream (CLI parse validates
        // `inputs.len() == read_structures.len()`, and the zipper rejects
        // per-chunk stream desync), so a mismatch here indicates a pipeline
        // logic regression, NOT malformed user input. Fail loud rather than
        // letting the `.zip()` below silently truncate sequencing reads. A
        // `debug_assert_eq!` surfaces it loudly in test/debug builds; the
        // release `Err` makes it a safe abort instead of silent data loss (do
        // not rely on the debug assert alone — that re-introduces the
        // silent-truncation gap in release builds).
        debug_assert_eq!(
            fq_template.records.len(),
            read_structures.len(),
            "ExtractStep: template record count != read-structure count (internal invariant \
             violated)"
        );
        if fq_template.records.len() != read_structures.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "ExtractStep internal invariant violated in batch_serial {serial}: template \
                     {:?} has {} record(s) but {} read structure(s) were provided; counts are \
                     validated equal at CLI parse and enforced per-chunk by the zipper, so this \
                     indicates a pipeline logic error, not malformed input",
                    String::from_utf8_lossy(&fq_template.name),
                    fq_template.records.len(),
                    read_structures.len(),
                ),
            ));
        }

        let mut fastq_sets: Vec<FastqSet> = Vec::with_capacity(fq_template.records.len());
        for (record, rs) in fq_template.records.iter().zip(read_structures.iter()) {
            let fastq_set = FastqSet::from_record_with_structure(
                record.name(),
                record.sequence(),
                record.quality(),
                rs,
                &[], // No skip reasons
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            fastq_sets.push(fastq_set);
        }

        let combined = FastqSet::combine_readsets(fastq_sets);

        let raw_records = make_raw_records_from_fastq_set(&combined, extract_opts)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let template = Template::from_records(raw_records)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        templates.push(template);
    }

    // Count emitted BAM records, not templates: each template carries one
    // record per read (R1, R2, …), and the finalize hook reports "records
    // emitted". Summing `templates.len()` would undercount on paired-end input.
    let count: u64 = templates.iter().map(|t| t.read_count() as u64).sum();
    records_emitted.fetch_add(count, Ordering::Relaxed);
    Ok(BamTemplateBatch::new(serial, templates))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::extract::{ExtractOptions, QualityEncoding};
    use crate::fastq::FastqSet;
    use crate::sam::SamTag;
    use fgumi_raw_bam::AsTagBytes;
    use fgumi_raw_bam::fields::RawRecordView;

    /// Build a minimal `ExtractOptions` with standard quality encoding.
    fn default_extract_opts() -> ExtractOptions {
        ExtractOptions {
            sample: "S1".to_string(),
            library: "L1".to_string(),
            platform: None,
            platform_unit: None,
            read_group_id: "A".to_string(),
            comments: Vec::new(),
            barcode: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            run_date: None,
            quality_encoding: QualityEncoding::Standard,
            store_umi_quals: false,
            store_cell_quals: false,
            single_tag: None,
            annotate_read_names: false,
            extract_umis_from_read_names: false,
            store_sample_barcode_qualities: false,
            async_reader: false,
        }
    }

    #[test]
    fn make_raw_records_single_end_template() {
        // Construct a FastqSet for a single-end read with a +T read structure.
        let rs = "10T".parse::<ReadStructure>().unwrap();
        let fq_set =
            FastqSet::from_record_with_structure(b"read1", b"ACGTACGTAC", b"IIIIIIIIII", &rs, &[])
                .unwrap();

        let opts = default_extract_opts();
        let records = make_raw_records_from_fastq_set(&fq_set, &opts).unwrap();

        assert_eq!(records.len(), 1, "single-end should produce one record");
        let view = RawRecordView::new(&records[0]);
        assert_eq!(view.read_name(), b"read1");
        // Single-end: should NOT have PAIRED flag
        assert_eq!(
            view.flags() & fgumi_raw_bam::fields::flags::PAIRED,
            0,
            "single-end should not be paired"
        );
        // Should be UNMAPPED
        assert_ne!(view.flags() & fgumi_raw_bam::fields::flags::UNMAPPED, 0, "should be unmapped");
        // Verify RG tag is present
        let aux = fgumi_raw_bam::fields::aux_data_slice(&records[0]);
        let rg_tag = SamTag::RG.as_tag_bytes();
        assert!(aux.windows(2).any(|w| w == rg_tag), "RG tag should be present in the record");
    }

    #[test]
    fn make_raw_records_paired_end_template() {
        let rs = "5T".parse::<ReadStructure>().unwrap();
        let fq_set1 =
            FastqSet::from_record_with_structure(b"read1", b"ACGTG", b"IIIII", &rs, &[]).unwrap();
        let fq_set2 =
            FastqSet::from_record_with_structure(b"read1", b"TGCAA", b"IIIII", &rs, &[]).unwrap();

        let combined = FastqSet::combine_readsets(vec![fq_set1, fq_set2]);
        let opts = default_extract_opts();
        let records = make_raw_records_from_fastq_set(&combined, &opts).unwrap();

        assert_eq!(records.len(), 2, "paired-end should produce two records");

        let v0 = RawRecordView::new(&records[0]);
        let v1 = RawRecordView::new(&records[1]);

        // Both should be PAIRED + UNMAPPED + MATE_UNMAPPED
        for v in [&v0, &v1] {
            assert_ne!(v.flags() & fgumi_raw_bam::fields::flags::PAIRED, 0);
            assert_ne!(v.flags() & fgumi_raw_bam::fields::flags::UNMAPPED, 0);
            assert_ne!(v.flags() & fgumi_raw_bam::fields::flags::MATE_UNMAPPED, 0);
        }

        // R1 = FIRST_SEGMENT, R2 = LAST_SEGMENT
        assert_ne!(v0.flags() & fgumi_raw_bam::fields::flags::FIRST_SEGMENT, 0);
        assert_ne!(v1.flags() & fgumi_raw_bam::fields::flags::LAST_SEGMENT, 0);
    }

    #[test]
    fn make_raw_records_with_umi() {
        // Read structure: 3M7T (3bp molecular barcode, 7bp template)
        let rs = "3M7T".parse::<ReadStructure>().unwrap();
        let fq_set =
            FastqSet::from_record_with_structure(b"read1", b"AAACCCCCCC", b"IIIIIIIIII", &rs, &[])
                .unwrap();

        let opts = default_extract_opts();
        let records = make_raw_records_from_fastq_set(&fq_set, &opts).unwrap();

        assert_eq!(records.len(), 1);
        // Verify the RX tag contains the UMI "AAA"
        let aux = fgumi_raw_bam::fields::aux_data_slice(&records[0]);
        let rx_tag = SamTag::RX.as_tag_bytes();
        let pos = aux.windows(2).position(|w| w == rx_tag);
        assert!(pos.is_some(), "RX tag should be present when UMI segments exist");
    }

    #[test]
    fn make_raw_records_builds_valid_template() {
        // End-to-end: build RawRecords, then verify Template::from_records succeeds.
        let rs = "5T".parse::<ReadStructure>().unwrap();
        let fq_set1 =
            FastqSet::from_record_with_structure(b"qname", b"ACGTG", b"IIIII", &rs, &[]).unwrap();
        let fq_set2 =
            FastqSet::from_record_with_structure(b"qname", b"TGCAA", b"IIIII", &rs, &[]).unwrap();

        let combined = FastqSet::combine_readsets(vec![fq_set1, fq_set2]);
        let opts = default_extract_opts();
        let records = make_raw_records_from_fastq_set(&combined, &opts).unwrap();

        let template = Template::from_records(records).unwrap();
        assert_eq!(template.name, b"qname");
        assert_eq!(template.read_count(), 2);
        assert!(template.r1.is_some());
        assert!(template.r2.is_some());
    }

    // ─────────────────────────────────────────────────────────────────────────
    // S5a2-002: end-to-end coverage of the extract batch conversion, including
    // the S5a2-001 record-count / read-structure mismatch hard-error.
    // ─────────────────────────────────────────────────────────────────────────

    use crate::fastq_parse::FastqRecord;
    use crate::grouper::FastqTemplate;
    use crate::pipeline::core::item::Ordered;

    /// Build a `FastqRecord` from a name/sequence (quality = all `I`).
    fn fq_record(name: &str, seq: &str) -> FastqRecord {
        let qual: String = std::iter::repeat_n('I', seq.len()).collect();
        FastqRecord::from_slice(format!("@{name}\n{seq}\n+\n{qual}\n").as_bytes()).unwrap()
    }

    /// A paired-end `FastqTemplateBatch`: each template has two records (R1+R2).
    fn paired_template_batch(batch_serial: u64, n_templates: usize) -> FastqTemplateBatch {
        let templates = (0..n_templates)
            .map(|i| FastqTemplate {
                name: format!("read{i}").into_bytes(),
                records: vec![
                    fq_record(&format!("read{i}"), "ACGTG"),
                    fq_record(&format!("read{i}"), "TGCAA"),
                ],
            })
            .collect();
        FastqTemplateBatch::new(batch_serial, templates)
    }

    /// `extract_batch` propagates the input `batch_serial` to the output and
    /// bumps `records_emitted` by the number of BAM records (reads) emitted,
    /// not the number of templates: 3 paired templates × 2 reads each = 6.
    #[test]
    fn extract_batch_propagates_serial_and_counts_records() {
        let read_structures = vec!["5T".parse::<ReadStructure>().unwrap(); 2];
        let opts = default_extract_opts();
        let emitted = AtomicU64::new(0);

        let batch = paired_template_batch(7, 3);
        let out = extract_batch(&read_structures, &opts, &emitted, &batch).unwrap();

        assert_eq!(out.ordinal(), 7, "output batch_serial must equal the input batch_serial");
        assert_eq!(out.templates.len(), 3, "all templates converted");
        let total_records: usize = out.templates.iter().map(Template::read_count).sum();
        assert_eq!(total_records, 6, "3 paired templates yield 6 records");
        assert_eq!(
            emitted.load(Ordering::Relaxed),
            6,
            "records_emitted bumped by emitted record count, not template count"
        );
    }

    /// A template whose record count differs from `read_structures.len()` is a
    /// structural invariant violation: `extract_batch` must hard-error with
    /// `InvalidData` rather than silently truncating reads via `zip`. (S5a2-001)
    #[test]
    fn extract_batch_errors_on_record_count_mismatch() {
        // Two read structures but a template carrying only ONE record.
        let read_structures = vec!["5T".parse::<ReadStructure>().unwrap(); 2];
        let opts = default_extract_opts();
        let emitted = AtomicU64::new(0);

        let short_template =
            FastqTemplate { name: b"read0".to_vec(), records: vec![fq_record("read0", "ACGTG")] };
        let batch = FastqTemplateBatch::new(0, vec![short_template]);

        // In debug builds the `debug_assert_eq!` fires first; catch the panic so
        // the test asserts the loud-failure contract in both build profiles. In
        // release the `Err(InvalidData)` path is taken.
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            extract_batch(&read_structures, &opts, &emitted, &batch)
        }));

        match result {
            Err(_panic) => { /* debug_assert_eq! fired — loud failure, contract met */ }
            Ok(Ok(_)) => panic!("expected a mismatch error or assertion, but conversion succeeded"),
            Ok(Err(err)) => {
                assert_eq!(
                    err.kind(),
                    io::ErrorKind::InvalidData,
                    "mismatch must surface as InvalidData"
                );
                assert!(
                    err.to_string().contains("internal invariant violated"),
                    "error must be framed as an internal invariant: {err}"
                );
            }
        }
    }
}
