//! `ExtractStep` — parallel FASTQ-to-BAM conversion step.
//!
//! Transforms a [`FastqTemplateBatch`] (N `FastqTemplate`s) into a
//! [`BamTemplateBatch`] (N `Template`s) by:
//!
//! 1. Applying read structures to each `FastqRecord` in the template,
//!    segmenting them zero-copy into a combined `FastqSetView`.
//! 2. Calling `make_raw_records_from_view` to produce
//!    `Vec<RawRecord>`.
//! 3. Building a [`Template`] via
//!    [`Template::from_records`].
//!
//! The step is `Parallel` (closure-driven via [`ProcessOrdered`]) and
//! preserves batch ordering via `BamTemplateBatch::batch_serial`.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::read_structure::ReadStructure;

use crate::commands::extract::{ExtractOptions, make_raw_records_from_view};
use crate::fastq::{FastqSegmentView, FastqSetView};
use crate::pipeline::steps::process::ProcessOrdered;
use crate::pipeline::steps::source::fastq_zip::FastqTemplateBatch;
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::template::Template;

/// Build a `Parallel` step that converts [`FastqTemplateBatch`] into
/// [`BamTemplateBatch`].
///
/// Each input batch's `FastqTemplate`s are independently converted: read
/// structures are applied (segmenting the source records zero-copy into a
/// `FastqSetView`), UMI/barcode tags are extracted via
/// `make_raw_records_from_view`, and the resulting `RawRecord`s are
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
        // logic regression, NOT malformed user input. Fail loud with a
        // release-safe `Err` rather than letting the `.zip()` below silently
        // truncate sequencing reads. This is deliberately NOT a `debug_assert!`:
        // a debug-only assert both downgrades the abort to a panic AND shadows
        // this branch from CI's debug-build tests, leaving the silent-truncation
        // guard the comment warns about unexercised.
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

        // Zero-copy segmentation: build one borrowed `FastqSetView` over the
        // source records' bases/qualities instead of allocating an owned
        // `FastqSet` per stream and then `combine_readsets`'ing them. The
        // combined header is the FIRST record's name (matching
        // `FastqSet::combine_readsets`, which keeps the first set's header), and
        // segments are appended in input order (R1 then R2, …) — byte-identical
        // to the owned path, but with no per-segment `to_vec` and no intermediate
        // per-stream `FastqSet` allocation.
        let mut segments: Vec<FastqSegmentView<'_>> = Vec::new();
        for (record, rs) in fq_template.records.iter().zip(read_structures.iter()) {
            FastqSetView::segment_into(
                record.name(),
                record.sequence(),
                record.quality(),
                rs,
                &[], // No skip reasons
                &mut segments,
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }
        // `.first()` (with a descriptive `expect`) rather than `records[0]`: the
        // count-equality guard above plus the non-empty `read_structures`
        // (validated at CLI parse / `validate_template_count`) make this
        // unreachable, but a labeled panic preserves the diagnostic the old
        // `FastqSet::combine_readsets` empty-vec assert carried.
        let header = fq_template
            .records
            .first()
            .expect("ExtractStep: template has no records despite the count-equality check above")
            .name();
        let combined = FastqSetView { header, segments };

        let raw_records = make_raw_records_from_view(&combined, extract_opts)
            // Render the full anyhow context chain into the io::Error string (the
            // alternate `{:#}` form: "context: cause"). A bare
            // `io::Error::new(_, e)` Displays only the top context, dropping the
            // root cause — e.g. it would keep "could not write the record for read
            // …" but lose the "read name too long" underneath, so the chain path's
            // error would diverge from the serial oracle's. Preserving the chain
            // keeps the two paths' messages in parity.
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{e:#}")))?;

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
    use crate::commands::extract::{
        ExtractOptions, QualityEncoding, make_raw_records_from_fastq_set,
    };
    use crate::fastq::FastqSet;
    use crate::sam::SamTag;
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
            check_crc: false,
            no_check_crc: false,
        }
    }

    /// Options exercising every optional tag branch, so the parity test below
    /// covers CB/CY, BC/QT, RX/QX, `--single-tag`, and `--annotate-read-names`.
    fn all_tags_opts() -> ExtractOptions {
        ExtractOptions {
            store_umi_quals: true,
            store_cell_quals: true,
            store_sample_barcode_qualities: true,
            single_tag: Some("ZU".parse().expect("valid tag")),
            annotate_read_names: true,
            ..default_extract_opts()
        }
    }

    /// The zero-copy view path (`make_raw_records_from_view`, used by the
    /// pipeline `extract_batch`) must emit BYTE-IDENTICAL `RawRecord`s to the
    /// owned path (`make_raw_records_from_fastq_set`, via `FastqSet` +
    /// `combine_readsets`) — the refactor removes per-segment `to_vec` copies but
    /// must not change a single output byte. Exercised across single-end,
    /// paired-end, and a full-structure paired read with UMI + cell + sample
    /// barcodes under the all-tags options.
    #[rstest::rstest]
    // (r1_structure, r1_seq, r1_qual, r2 (structure,seq,qual) or None, opts_all_tags)
    #[case::single_end_template("10T", b"ACGTACGTAC".as_slice(), b"IIIIIIIIII".as_slice(), None, false)]
    #[case::paired_template("5T", b"ACGTG".as_slice(), b"IIIII".as_slice(), Some(("5T", b"TGCAA".as_slice(), b"JJJJJ".as_slice())), false)]
    #[case::umi_only("3M7T", b"AAACCCCCCC".as_slice(), b"IIIIIIIIII".as_slice(), None, true)]
    #[case::full_structure_all_tags(
        "3C2B4M8T", b"CCCBBUUUUTTTTTTTT".as_slice(), b"IIIIIIIIIIIIIIIII".as_slice(),
        Some(("8T", b"GGGGGGGG".as_slice(), b"JJJJJJJJ".as_slice())), true
    )]
    fn view_path_matches_owned_path_byte_for_byte(
        #[case] r1_rs: &str,
        #[case] r1_seq: &[u8],
        #[case] r1_qual: &[u8],
        #[case] r2: Option<(&str, &[u8], &[u8])>,
        #[case] all_tags: bool,
    ) {
        let opts = if all_tags { all_tags_opts() } else { default_extract_opts() };
        // R1 and R2 carry DISTINCT names (both 8-field so UMI-from-name is inert
        // unless enabled). The combined header must come from the FIRST record
        // (R1): `combine_readsets` keeps the first set's header and `extract_batch`
        // uses `records.first()`. Distinct names let the checks below detect a
        // regression that sourced the header from the wrong record.
        let r1_name = b"inst:1:fc:1:1:1:1:ACGT".as_slice();
        let r2_name = b"inst:1:fc:1:1:1:1:TTTT".as_slice();

        // ---- Owned path: FastqSet(s) + combine_readsets ----
        let rs1 = r1_rs.parse::<ReadStructure>().unwrap();
        let mut owned_sets = vec![
            FastqSet::from_record_with_structure(r1_name, r1_seq, r1_qual, &rs1, &[]).unwrap(),
        ];
        if let Some((r2_rs, r2_seq, r2_qual)) = r2 {
            let rs2 = r2_rs.parse::<ReadStructure>().unwrap();
            owned_sets.push(
                FastqSet::from_record_with_structure(r2_name, r2_seq, r2_qual, &rs2, &[]).unwrap(),
            );
        }
        let owned = FastqSet::combine_readsets(owned_sets);
        let owned_records = make_raw_records_from_fastq_set(&owned, &opts).unwrap();

        // ---- View path: FastqSetView::segment_into over borrowed slices ----
        // Header is the FIRST record's name, mirroring `extract_batch`.
        let mut segments = Vec::new();
        FastqSetView::segment_into(r1_name, r1_seq, r1_qual, &rs1, &[], &mut segments).unwrap();
        if let Some((r2_rs, r2_seq, r2_qual)) = r2 {
            let rs2 = r2_rs.parse::<ReadStructure>().unwrap();
            FastqSetView::segment_into(r2_name, r2_seq, r2_qual, &rs2, &[], &mut segments).unwrap();
        }
        let view = FastqSetView { header: r1_name, segments };
        let view_records = make_raw_records_from_view(&view, &opts).unwrap();

        // Guard against a vacuous pass: a regression emitting zero records from
        // both paths would otherwise satisfy the length + zip assertions trivially.
        assert!(!owned_records.is_empty(), "expected at least one record");

        // Byte-for-byte identical record bytes (and count).
        assert_eq!(owned_records.len(), view_records.len(), "record count differs");
        for (i, (o, v)) in owned_records.iter().zip(view_records.iter()).enumerate() {
            assert_eq!(o.as_ref(), v.as_ref(), "record {i} bytes differ (view vs owned)");
        }

        // Pin that the header is sourced from R1, not R2: a view built with R2's
        // (distinct) name must produce DIFFERENT bytes, proving the header
        // participates in the output and that first-record-wins is what both
        // paths do.
        if r2.is_some() {
            let mut alt_segments = Vec::new();
            FastqSetView::segment_into(r2_name, r1_seq, r1_qual, &rs1, &[], &mut alt_segments)
                .unwrap();
            if let Some((r2_rs, r2_seq, r2_qual)) = r2 {
                let rs2 = r2_rs.parse::<ReadStructure>().unwrap();
                FastqSetView::segment_into(r2_name, r2_seq, r2_qual, &rs2, &[], &mut alt_segments)
                    .unwrap();
            }
            let alt_view = FastqSetView { header: r2_name, segments: alt_segments };
            let alt_records = make_raw_records_from_view(&alt_view, &opts).unwrap();
            assert!(!alt_records.is_empty(), "expected at least one record");
            assert_ne!(
                alt_records[0].as_ref(),
                owned_records[0].as_ref(),
                "header did not affect output — first-record-wins is not pinned"
            );
        }
    }

    /// How a segmentation guard classified an input.
    #[derive(Debug, PartialEq, Clone, Copy)]
    enum GuardOutcome {
        /// Rejected with an error (length mismatch, over-long fixed, or too-few
        /// bases without a `TooFewBases` skip).
        Rejected,
        /// Accepted-but-skipped: too few bases with `TooFewBases` in skip reasons
        /// (owned path returns a `skip_reason`-tagged empty set; view path leaves
        /// `out` unchanged).
        Skipped,
        /// Segmented normally into one or more segments.
        Segmented,
    }

    /// `FastqSetView::segment_into` re-implements the same three input guards as
    /// its owned sibling `FastqSet::from_record_with_structure` (length mismatch,
    /// over-long fixed, and too-few-bases bail-vs-skip). The byte-parity test
    /// above only feeds valid inputs, so this pins that the two siblings CLASSIFY
    /// every error/skip input identically — a future edit to one guard cannot
    /// silently diverge them (guard-set parity, the repo's key check for these
    /// paired parsing entry points).
    #[rstest::rstest]
    // (structure, seq, qual, allow_too_few_skip, expected)
    #[case::len_mismatch("5T", b"ACGTG".as_slice(), b"IIII".as_slice(), false, GuardOutcome::Rejected)]
    #[case::over_long_fixed("5T", b"ACGTGT".as_slice(), b"IIIIII".as_slice(), false, GuardOutcome::Rejected)]
    #[case::over_long_fixed_even_with_skip("5T", b"ACGTGT".as_slice(), b"IIIIII".as_slice(), true, GuardOutcome::Rejected)]
    #[case::too_few_no_skip("10T", b"ACGT".as_slice(), b"IIII".as_slice(), false, GuardOutcome::Rejected)]
    #[case::too_few_with_skip("10T", b"ACGT".as_slice(), b"IIII".as_slice(), true, GuardOutcome::Skipped)]
    #[case::ok_exact("5T", b"ACGTG".as_slice(), b"IIIII".as_slice(), false, GuardOutcome::Segmented)]
    fn segment_into_matches_owned_guard_classification(
        #[case] structure: &str,
        #[case] seq: &[u8],
        #[case] qual: &[u8],
        #[case] allow_too_few_skip: bool,
        #[case] expected: GuardOutcome,
    ) {
        use crate::fastq::SkipReason;
        let name = b"read1".as_slice();
        let rs = structure.parse::<ReadStructure>().unwrap();
        let skip: &[SkipReason] = if allow_too_few_skip { &[SkipReason::TooFewBases] } else { &[] };

        // Owned sibling classification.
        let owned_outcome = match FastqSet::from_record_with_structure(name, seq, qual, &rs, skip) {
            Err(_) => GuardOutcome::Rejected,
            Ok(set) if set.skip_reason.is_some() => GuardOutcome::Skipped,
            Ok(_) => GuardOutcome::Segmented,
        };

        // Borrowed sibling classification. A valid structure always yields >=1
        // segment, so an empty `out` after `Ok(())` is unambiguously a skip.
        let mut out = Vec::new();
        let view_outcome = match FastqSetView::segment_into(name, seq, qual, &rs, skip, &mut out) {
            Err(_) => GuardOutcome::Rejected,
            Ok(()) if out.is_empty() => GuardOutcome::Skipped,
            Ok(()) => GuardOutcome::Segmented,
        };

        assert_eq!(owned_outcome, expected, "owned path misclassified");
        assert_eq!(view_outcome, expected, "view path misclassified");
        assert_eq!(owned_outcome, view_outcome, "owned and view guards diverged");
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
        // Verify the RG tag carries the read-group id value "A" (not just that
        // the two bytes [R,G] appear somewhere in the aux blob).
        let aux = fgumi_raw_bam::fields::aux_data_slice(&records[0]);
        assert_eq!(
            fgumi_raw_bam::tags::find_string_tag(aux, SamTag::RG),
            Some(&b"A"[..]),
            "RG must carry the read-group id 'A'"
        );
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
        // Verify the RX tag carries the 3bp UMI VALUE "AAA" (the leading 3M of
        // the "3M7T" read structure), not merely that an RX tag is present.
        let aux = fgumi_raw_bam::fields::aux_data_slice(&records[0]);
        assert_eq!(
            fgumi_raw_bam::tags::find_string_tag(aux, SamTag::RX),
            Some(&b"AAA"[..]),
            "RX must carry the UMI value 'AAA'"
        );
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
        assert_eq!(template.name(), b"qname");
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
        assert_eq!(out.templates().len(), 3, "all templates converted");
        let total_records: usize = out.templates().iter().map(Template::read_count).sum();
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

        // The mismatch must surface as a release-safe `Err(InvalidData)` in ALL
        // build profiles (there is no `debug_assert` to panic first), so CI's
        // debug build exercises exactly the silent-truncation guard the code
        // comment warns about — not merely a debug-only panic.
        let err = extract_batch(&read_structures, &opts, &emitted, &batch)
            .expect_err("record-count mismatch must return Err, not truncate");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "mismatch must surface as InvalidData");
        assert!(
            err.to_string().contains("internal invariant violated"),
            "error must be framed as an internal invariant: {err}"
        );
    }

    #[test]
    fn extract_batch_errors_on_record_surplus() {
        // Mirror of the shortfall case: a template carrying MORE records than
        // there are read structures. The guard rejects both directions, so a
        // future rewrite to a one-sided check (`records.len() < structures.len()`)
        // must fail here — the surplus direction is tested independently.
        let read_structures = vec!["5T".parse::<ReadStructure>().unwrap(); 2];
        let opts = default_extract_opts();
        let emitted = AtomicU64::new(0);

        let long_template = FastqTemplate {
            name: b"read0".to_vec(),
            records: vec![
                fq_record("read0", "ACGTG"),
                fq_record("read0", "ACGTG"),
                fq_record("read0", "ACGTG"),
            ],
        };
        let batch = FastqTemplateBatch::new(0, vec![long_template]);

        let err = extract_batch(&read_structures, &opts, &emitted, &batch)
            .expect_err("record-count surplus must return Err, not ignore the extra records");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "surplus must surface as InvalidData");
        assert!(
            err.to_string().contains("internal invariant violated"),
            "error must be framed as an internal invariant: {err}"
        );
    }
}
