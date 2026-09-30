//! `MergeAlignedStep` — the shared, `Parallel` consumer of an align backend's
//! `ZipperBatch` stream (the subprocess backend's; the
//! in-process backend calls `merge_zipper_batch` in its own pair/emit step).
//!
//! Each `ZipperBatch` carries both halves of a batch (aligner-emitted `mapped`
//! templates positionally paired with the original `unmapped` templates) plus a
//! dense serial, so the merge is per-template and needs no cross-batch state:
//! `merge_one_template_with` (`merge_raw_with` plus optional bisulfite restore)
//! is run for each pair, folding record-count and heap-size accounting into the
//! same pass. That independence is why this step is `Parallel` with an
//! `OrderedBytesSingle<BamTemplateBatch>` + `BranchOrdering::ByItemOrdinal`
//! output: workers merge batches concurrently and the framework restores input
//! order from each `BamTemplateBatch`'s serial.
//!
//! This was previously the `merge_zipper_batch` call inlined on the (Serial)
//! `AlignAndMergeStep`'s dispatching worker; promoting it to its own `Parallel`
//! step is the split the module doc on the former `align_and_merge.rs`
//! anticipated. The merge logic itself is unchanged.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use noodles::sam::Header;

use crate::commands::zipper::{ZipperTags, merge_one_template_with};
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::align::ZipperBatch;
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::reference::ReferenceReader;
use crate::template::Template;
use crate::umi::TagInfo;

/// Immutable, Arc-friendly configuration for [`MergeAlignedStep`].
///
/// One `Arc` clone is shared by every per-worker copy of the step; the cost is
/// negligible compared to the data flowing through.
pub(crate) struct MergeConfig {
    /// Tag-merge rules (remove / reverse / revcomp) plumbed to `merge_raw`.
    pub(crate) tag_info: Arc<TagInfo>,

    /// Whether to skip TC (template-coordinate) tag handling. Mirrors
    /// `ZipperMergeConfig::skip_tc_tags`.
    pub(crate) skip_tc_tags: bool,

    /// Optional reference reader used by
    /// `restore_unconverted_bases_in_raw_template` (bisulfite path).
    /// `None` for normal alignment.
    pub(crate) reference: Option<Arc<ReferenceReader>>,

    /// Partial output header (dict-derived `@SQ` + unmapped-derived
    /// `@HD`/`@CO`/`@RG`/`@PG` + fgumi `@PG`) used by `merge_one_template_with`.
    pub(crate) partial_output_header: Arc<Header>,

    /// Counter for records emitted downstream. Exposed back to the caller after
    /// `Pipeline::run` returns so summary logging can report a real throughput
    /// number.
    pub(crate) records_emitted: Arc<AtomicU64>,

    /// Byte limit for the downstream output queue (`OrderedBytesSingle`
    /// `ByteBounded`) that merged `BamTemplateBatch`es are pushed onto.
    pub(crate) output_byte_limit: u64,
}

/// `Parallel + ByItemOrdinal` step that merges [`ZipperBatch`]es into
/// `BamTemplateBatch`es via [`merge_zipper_batch`].
pub(crate) struct MergeAlignedStep {
    cfg: Arc<MergeConfig>,
    /// Precomputed tag-merge bitsets, built once from `cfg.tag_info` and reused
    /// for every template across every batch — `cfg.tag_info` is immutable for
    /// the step's whole lifetime, so there is no reason to rebuild `ZipperTags`
    /// (three `TagBitset` allocations) per batch, let alone per template.
    tags: Arc<ZipperTags>,
    held: HeldSlot<Unpushed<BamTemplateBatch>>,
}

impl MergeAlignedStep {
    /// Build the shared merge step from `cfg`, precomputing the `ZipperTags`
    /// bitsets once.
    pub(crate) fn from_shared(cfg: Arc<MergeConfig>) -> Self {
        let tags = Arc::new(ZipperTags::from_tag_info(&cfg.tag_info));
        Self { cfg, tags, held: HeldSlot::new() }
    }
}

impl Step for MergeAlignedStep {
    type Input = ZipperBatch;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "MergeAligned",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.cfg.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Drain the held output slot first, before popping more input.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    // `Contention` (not `NoProgress`) — returning `NoProgress`
                    // from a Parallel step on a held slot risks the framework
                    // marking this worker `Skip` if upstream is also drained,
                    // silently dropping the held item. Mirrors
                    // `TemplatesToRecordBatch` / `BgzfDecompress`.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(zb) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output flushed above) and this step will never
            // push again — report Finished. For a Parallel step only the last
            // clone to finish closes the shared output (gated by the driver's
            // StepDrainCounter).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let merged = merge_zipper_batch(zb, &self.cfg, &self.tags)?;
        match ctx.outputs.push(merged) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        Self { cfg: Arc::clone(&self.cfg), tags: Arc::clone(&self.tags), held: HeldSlot::new() }
    }
}

/// Merge one `ZipperBatch` into a `BamTemplateBatch`. Takes the caller's
/// precomputed `tags` (built once for the whole step, since `cfg.tag_info` is
/// immutable) and, per template, runs `merge_one_template_with`
/// (`merge_raw_with` plus optional bisulfite restore); folds record-count and
/// heap-size accounting into the same single pass so the resulting
/// `BamTemplateBatch` doesn't re-walk the templates to compute `total_bytes`.
pub(crate) fn merge_zipper_batch(
    zb: ZipperBatch,
    cfg: &MergeConfig,
    tags: &ZipperTags,
) -> io::Result<BamTemplateBatch> {
    let ZipperBatch { serial, mapped, unmapped } = zb;
    // Release-safe guard (not a debug_assert): a violated length invariant would
    // otherwise make the `zip` below silently truncate to the shorter side,
    // dropping templates with no error. Fail loudly instead — mirroring
    // `extract_batch`'s record-count guard. (A bare debug_assert would also
    // shadow this Err path from CI's debug-build tests.)
    if mapped.len() != unmapped.templates().len() {
        return Err(io::Error::other(format!(
            "align-and-merge: ZipperBatch invariant violated — {} mapped vs {} unmapped \
             templates; refusing to zip mismatched halves",
            mapped.len(),
            unmapped.templates().len(),
        )));
    }

    let mut merged: Vec<Template> = Vec::with_capacity(mapped.len());
    let mut total_records: u64 = 0;
    let mut total_bytes: usize = 0;
    // One aux-rebuild scratch buffer reused across every template in this batch,
    // mirroring the standalone `Zipper::run` path (its allocation is reused, not
    // re-allocated per template) — see fgumi #971.
    let mut aux_scratch: Vec<u8> = Vec::new();
    for (mut mapped_template, unmapped_template) in mapped.into_iter().zip(unmapped.templates()) {
        merge_one_template_with(
            unmapped_template,
            &mut mapped_template,
            tags,
            cfg.skip_tc_tags,
            cfg.reference.as_deref(),
            &cfg.partial_output_header,
            &mut aux_scratch,
        )
        .map_err(|e| io::Error::other(format!("align-and-merge: {e:#}")))?;

        total_records += mapped_template.records.len() as u64;
        total_bytes += mapped_template.heap_size();
        merged.push(mapped_template);
    }

    cfg.records_emitted.fetch_add(total_records, Ordering::Relaxed);
    Ok(BamTemplateBatch::from_parts(serial, merged, total_bytes))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::item::Ordered;

    fn make_test_cfg() -> MergeConfig {
        MergeConfig {
            tag_info: Arc::new(TagInfo::new(vec![], vec![], vec![])),
            skip_tc_tags: true,
            reference: None,
            partial_output_header: Arc::new(Header::default()),
            records_emitted: Arc::new(AtomicU64::new(0)),
            output_byte_limit: 1024 * 1024,
        }
    }

    fn make_record(qname: &[u8], flags: u16) -> fgumi_raw_bam::RawRecord {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(qname).flags(flags).sequence(b"ACGT").qualities(b"IIII");
        b.build()
    }

    /// Build a single-record `RawRecord` carrying one string tag.
    fn make_record_with_string_tag(
        qname: &[u8],
        flags: u16,
        tag: crate::sam::SamTag,
        value: &[u8],
    ) -> fgumi_raw_bam::RawRecord {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(qname)
            .flags(flags)
            .sequence(b"ACGT")
            .qualities(b"IIII")
            .add_string_tag(tag, value);
        b.build()
    }

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let step = MergeAlignedStep::from_shared(Arc::new(make_test_cfg()));
        let p = step.profile();
        assert_eq!(p.name, "MergeAligned");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky, "Parallel so any worker can dispatch");
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        assert_eq!(p.output_queues.len(), 1);
    }

    #[test]
    fn new_worker_copy_shares_config_and_resets_held_slot() {
        let step = MergeAlignedStep::from_shared(Arc::new(make_test_cfg()));
        let copy = step.new_worker_copy();
        assert!(Arc::ptr_eq(&step.cfg, &copy.cfg), "config Arc is shared, not rebuilt");
        assert!(Arc::ptr_eq(&step.tags, &copy.tags), "tags Arc is shared, not rebuilt");
        assert!(!copy.held.is_held(), "worker copy gets a fresh held slot");
    }

    /// `merge_zipper_batch` transfers the unmapped half's tags onto the paired
    /// mapped template, preserves the batch serial, and bumps `records_emitted`
    /// by the emitted record count — the core merge contract.
    #[test]
    fn merge_zipper_batch_transfers_unmapped_tags_and_counts_records() {
        let cfg = make_test_cfg();
        let tags = ZipperTags::from_tag_info(&cfg.tag_info);
        // Unmapped half carries RX; the mapped half (aligner output) does not.
        let unmapped_rec = make_record_with_string_tag(
            b"readA",
            fgumi_raw_bam::flags::UNMAPPED,
            crate::sam::SamTag::RX,
            b"ACGT",
        );
        let unmapped = Template::from_records(vec![unmapped_rec]).expect("unmapped template");
        let mapped =
            Template::from_records(vec![make_record(b"readA", 0)]).expect("mapped template");

        let zb = ZipperBatch {
            serial: 5,
            mapped: vec![mapped],
            unmapped: BamTemplateBatch::new(5, vec![unmapped]),
        };
        let out = merge_zipper_batch(zb, &cfg, &tags).expect("merge ok");

        assert_eq!(out.ordinal(), 5, "batch serial is preserved through the merge");
        assert_eq!(out.templates().len(), 1, "one merged template out");
        assert_eq!(
            cfg.records_emitted.load(Ordering::Relaxed),
            1,
            "records_emitted bumped by merged record count, not template count"
        );
        // The RX tag from the unmapped half must land on the merged mapped record.
        let merged_rec = &out.templates()[0].records[0];
        let aux = fgumi_raw_bam::fields::aux_data_slice(merged_rec);
        assert_eq!(
            fgumi_raw_bam::tags::find_string_tag(aux, crate::sam::SamTag::RX),
            Some(&b"ACGT"[..]),
            "the unmapped RX tag must transfer onto the correct mapped record"
        );
    }

    /// A pair split across a mid-pair `-K` cut zips each half with the aligner's
    /// unpaired record for that read. Both halves (and the second read's
    /// supplementary) must receive their own unmapped read's tags, and the
    /// second read's QC-fail flag must transfer. Before the split halves'
    /// pairing bits were cleared, the second half (still `PAIRED |
    /// LAST_SEGMENT`) looked for a mapped R2, found none, and silently copied
    /// nothing.
    #[test]
    fn split_pair_halves_receive_their_own_tags_and_qc_flag() {
        use fgumi_raw_bam::flags::{
            FIRST_SEGMENT, LAST_SEGMENT, MATE_UNMAPPED, PAIRED, QC_FAIL, REVERSE, SUPPLEMENTARY,
            UNMAPPED,
        };
        let cfg = make_test_cfg();
        let tags = ZipperTags::from_tag_info(&cfg.tag_info);
        let r1 = make_record_with_string_tag(
            b"pe10",
            PAIRED | FIRST_SEGMENT | UNMAPPED | MATE_UNMAPPED,
            crate::sam::SamTag::RX,
            b"AAAA",
        );
        let r2 = make_record_with_string_tag(
            b"pe10",
            PAIRED | LAST_SEGMENT | UNMAPPED | MATE_UNMAPPED | QC_FAIL,
            crate::sam::SamTag::RX,
            b"CCCC",
        );
        let pair = Template::from_records(vec![r1, r2]).expect("unmapped pair");
        let (first, second) =
            crate::pipeline::steps::align::split_pair_into_singles(pair).expect("split");
        let flags_of = |t: &Template| t.records()[0].flags();
        assert_eq!(flags_of(&first), UNMAPPED, "first half keeps only non-pairing bits");
        assert_eq!(flags_of(&second), UNMAPPED | QC_FAIL, "second half keeps QC_FAIL");

        // The aligner emitted both reads unpaired; the second has a supplementary.
        let mapped_first =
            Template::from_records(vec![make_record(b"pe10", 0)]).expect("mapped first half");
        let mapped_second = Template::from_records(vec![
            make_record(b"pe10", REVERSE),
            make_record(b"pe10", SUPPLEMENTARY),
        ])
        .expect("mapped second half");

        let zb = ZipperBatch {
            serial: 0,
            mapped: vec![mapped_first, mapped_second],
            unmapped: BamTemplateBatch::new(0, vec![first, second]),
        };
        let out = merge_zipper_batch(zb, &cfg, &tags).expect("merge ok");

        let rx_of = |rec: &fgumi_raw_bam::RawRecord| {
            fgumi_raw_bam::tags::find_string_tag(
                fgumi_raw_bam::fields::aux_data_slice(rec),
                crate::sam::SamTag::RX,
            )
            .map(<[u8]>::to_vec)
        };
        let first_out = &out.templates()[0].records;
        let second_out = &out.templates()[1].records;
        assert_eq!(rx_of(&first_out[0]), Some(b"AAAA".to_vec()), "first half gets R1's RX");
        assert_eq!(second_out.len(), 2, "second half keeps its supplementary");
        for rec in second_out {
            assert_eq!(rx_of(rec), Some(b"CCCC".to_vec()), "second half gets R2's RX");
            assert_ne!(rec.flags() & QC_FAIL, 0, "second half gets R2's QC-fail flag");
        }
        assert_eq!(first_out[0].flags() & QC_FAIL, 0, "first half's QC-pass status is unchanged");
    }

    /// A `ZipperBatch` whose mapped and unmapped halves differ in length is a
    /// structural invariant violation: `merge_zipper_batch` must hard-error
    /// (release-safe) rather than let `zip` silently truncate to the shorter
    /// side and drop templates.
    #[test]
    fn merge_zipper_batch_errors_on_length_mismatch() {
        let cfg = make_test_cfg();
        let tags = ZipperTags::from_tag_info(&cfg.tag_info);
        let mapped = Template::from_records(vec![make_record(b"readA", 0)]).expect("mapped");
        // One mapped template, zero unmapped -> lengths differ.
        let zb = ZipperBatch {
            serial: 0,
            mapped: vec![mapped],
            unmapped: BamTemplateBatch::new(0, Vec::new()),
        };
        let err = merge_zipper_batch(zb, &cfg, &tags)
            .expect_err("length mismatch must error, not truncate");
        assert!(
            err.to_string().contains("ZipperBatch invariant violated"),
            "error must name the invariant: {err}"
        );
    }

    /// Regression guard for hoisting the `ZipperTags` bitset build out of the
    /// per-template merge loop: `merge_zipper_batch` takes the tags as a
    /// caller-supplied `&ZipperTags` (built once for the whole step) and reused
    /// across every batch. A non-trivial `TagInfo` (one remove + one reverse +
    /// one revcomp tag) is applied across THREE templates in a single batch;
    /// every template's negative-strand read must get the same
    /// remove/reverse/revcomp treatment, not just the first one merged.
    #[test]
    fn merge_zipper_batch_applies_transforms_to_every_template() {
        let mut cfg = make_test_cfg();
        cfg.tag_info = Arc::new(TagInfo::new(
            vec!["XA".to_string()],
            vec!["XV".to_string()],
            vec!["XC".to_string()],
        ));
        let tags = ZipperTags::from_tag_info(&cfg.tag_info);

        let names: [&[u8]; 3] = [b"readA", b"readB", b"readC"];
        let mut mapped_templates = Vec::new();
        let mut unmapped_templates = Vec::new();
        for name in names {
            // Mapped (aligner) record: negative strand, single unpaired read,
            // carrying a stale XA tag that must be removed on merge.
            let mut mb = fgumi_raw_bam::SamBuilder::new();
            mb.read_name(name)
                .flags(fgumi_raw_bam::flags::REVERSE)
                .sequence(b"ACGT")
                .qualities(b"IIII")
                .add_string_tag(*b"XA", b"stale");
            let mapped_rec = mb.build();
            mapped_templates
                .push(Template::from_records(vec![mapped_rec]).expect("mapped template"));

            // Unmapped record carries the tags to remove/reverse/revcomp.
            let mut ub = fgumi_raw_bam::SamBuilder::new();
            ub.read_name(name)
                .flags(fgumi_raw_bam::flags::UNMAPPED)
                .sequence(b"ACGT")
                .qualities(b"IIII")
                .add_string_tag(*b"XV", b"abcde")
                .add_string_tag(*b"XC", b"AGAGG")
                .add_string_tag(*b"XA", b"drop-me");
            let unmapped_rec = ub.build();
            unmapped_templates
                .push(Template::from_records(vec![unmapped_rec]).expect("unmapped template"));
        }

        let zb = ZipperBatch {
            serial: 0,
            mapped: mapped_templates,
            unmapped: BamTemplateBatch::new(0, unmapped_templates),
        };
        let out = merge_zipper_batch(zb, &cfg, &tags).expect("merge ok");

        assert_eq!(out.templates().len(), 3, "all three templates survive the merge");
        for (i, template) in out.templates().iter().enumerate() {
            let rec = &template.records[0];
            let aux = fgumi_raw_bam::fields::aux_data_slice(rec);

            assert_eq!(
                fgumi_raw_bam::tags::find_string_tag(aux, *b"XV"),
                Some(&b"edcba"[..]),
                "template {i}: XV must be reversed on the negative-strand read"
            );
            assert_eq!(
                fgumi_raw_bam::tags::find_string_tag(aux, *b"XC"),
                Some(&b"CCTCT"[..]),
                "template {i}: XC must be reverse-complemented on the negative-strand read"
            );
            assert!(
                fgumi_raw_bam::tags::find_string_tag(aux, *b"XA").is_none(),
                "template {i}: XA must be removed (stale mapped copy + skipped on tag-copy)"
            );
        }
    }
}
