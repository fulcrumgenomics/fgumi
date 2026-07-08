//! Chain-builder support for [`crate::pipeline::chains::Stage::Fastq`] — the
//! terminal BAM → FASTQ encode.
//!
//! The encode step consumes the [`DecodedRecordBatch`] tail produced by the
//! BAM source preamble and emits a [`DecompressedBlock`] of FASTQ text per
//! batch, reusing [`write_fastq_record`] (the same encoder the standalone
//! command used). It is a `Parallel` step, so the work-stealing pool spreads
//! the encode across N workers; `add_sink` then routes the byte tail to a
//! (optionally BGZF-compressing) raw-file/stdout writer.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::commands::fastq::{FastqOptions, Segment, classify_segment, write_fastq_record};
use crate::logging::OperationTimer;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::process::{
    Process3Output, Process3WithWorkerState, ProcessWithWorkerState, process_with_worker_state,
    process3_with_worker_state,
};
use crate::pipeline::steps::types::{DecodedRecordBatch, DecompressedBlock};

/// Post-pipeline finalize hook for the fastq stage: logs the record count and
/// the per-stage wall-time line.
pub(crate) struct FastqFinalizeHook {
    pub(crate) records_written: Arc<AtomicU64>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for FastqFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let FastqFinalizeHook { records_written, timer } = *self;
        let n = records_written.load(Ordering::Relaxed);
        info!("Fastq: completed ({n} FASTQ records written)");
        timer.log_completion(n);
        Ok(())
    }
}

/// Per-worker scratch for the FASTQ encoder: the reusable sequence and quality
/// decode buffers `write_fastq_record` fills per record.
pub(crate) struct FastqEncodeState {
    seq_buf: Vec<u8>,
    qual_buf: Vec<u8>,
}

impl crate::pipeline::core::item::HeapSize for FastqEncodeState {}

fn new_encode_state() -> FastqEncodeState {
    FastqEncodeState { seq_buf: Vec::with_capacity(512), qual_buf: Vec::with_capacity(512) }
}

/// Bump the shared record counter and log a progress line every million records.
fn bump_progress(counter: &AtomicU64, written: u64) {
    if written == 0 {
        return;
    }
    let prev = counter.fetch_add(written, Ordering::Relaxed);
    if (prev + written) / 1_000_000 > prev / 1_000_000 {
        info!("Wrote {} FASTQ records", prev + written);
    }
}

/// Interleaved FASTQ encode step: [`DecodedRecordBatch`] → one
/// [`DecompressedBlock`] of FASTQ text (all passing records, in input order).
pub(crate) fn build_fastq_encode_step(
    opts: Arc<FastqOptions>,
    per_step_byte_limit: u64,
    records_written: Arc<AtomicU64>,
) -> ProcessWithWorkerState<
    DecodedRecordBatch,
    DecompressedBlock,
    impl Fn(&mut FastqEncodeState, DecodedRecordBatch) -> io::Result<DecompressedBlock>
    + Send
    + Sync
    + 'static,
    FastqEncodeState,
    impl Fn() -> FastqEncodeState + Send + Sync + 'static,
> {
    process_with_worker_state::<DecodedRecordBatch, DecompressedBlock, _, FastqEncodeState, _>(
        "FastqEncode",
        per_step_byte_limit,
        new_encode_state,
        move |state: &mut FastqEncodeState,
              batch: DecodedRecordBatch|
              -> io::Result<DecompressedBlock> {
            let DecodedRecordBatch { batch_serial, records, .. } = batch;
            let mut out = Vec::with_capacity(records.len().saturating_mul(512));
            let mut written = 0u64;
            for record in &records {
                let flags = record.record().flags();
                if !opts.passes_filters(flags) {
                    continue;
                }
                write_fastq_record(
                    &mut out,
                    record.record(),
                    flags,
                    opts.no_suffix,
                    opts.umi_header.as_ref(),
                    &mut state.seq_buf,
                    &mut state.qual_buf,
                )
                .map_err(|e| io::Error::other(format!("fastq encode: {e:#}")))?;
                written += 1;
            }
            bump_progress(&records_written, written);
            Ok(DecompressedBlock { batch_serial, bytes: out })
        },
    )
}

/// Routes and encodes a single [`DecodedRecordBatch`] into R1/R2/other FASTQ
/// byte blocks.
///
/// Routing matches `samtools fastq -1/-2/-0`: R1 = `FIRST_SEGMENT` set /
/// `LAST_SEGMENT` clear, R2 = the reverse, "other" = single-end or ambiguous
/// (see [`classify_segment`]). Split mode always omits the `/1`/`/2` suffix.
///
/// **Emits a block on all three branches for every batch** — an empty
/// `Vec<u8>` where a branch has no records this batch — to keep each ordered
/// branch's `batch_serial` sequence gap-free (a skipped serial would wedge that
/// branch's reorder stage). Never returns `None` on a branch.
///
/// Extracted from [`build_fastq_paired_encode_step`]'s closure so the routing
/// logic is directly unit-testable without going through the pipeline
/// machinery (see `mod tests`).
fn encode_paired_batch(
    state: &mut FastqEncodeState,
    batch: DecodedRecordBatch,
    opts: &FastqOptions,
    records_written: &AtomicU64,
) -> io::Result<Process3Output<DecompressedBlock, DecompressedBlock, DecompressedBlock>> {
    let DecodedRecordBatch { batch_serial, records, .. } = batch;
    let mut r1: Vec<u8> = Vec::new();
    let mut r2: Vec<u8> = Vec::new();
    let mut other: Vec<u8> = Vec::new();
    let mut written = 0u64;
    for record in &records {
        let flags = record.record().flags();
        if !opts.passes_filters(flags) {
            continue;
        }
        // Split mode: the R1/R2 file already identifies the mate, so the
        // `/1`/`/2` suffix is always omitted (matches `samtools fastq -1/-2`).
        let out = match classify_segment(flags) {
            Segment::Read1 => &mut r1,
            Segment::Read2 => &mut r2,
            Segment::Other => &mut other,
        };
        write_fastq_record(
            out,
            record.record(),
            flags,
            /* no_suffix */ true,
            opts.umi_header.as_ref(),
            &mut state.seq_buf,
            &mut state.qual_buf,
        )
        .map_err(|e| io::Error::other(format!("fastq encode: {e:#}")))?;
        written += 1;
    }
    bump_progress(records_written, written);
    // Always emit all three branches (empty Vec where none) — no gaps.
    Ok(Process3Output::both3(
        DecompressedBlock { batch_serial, bytes: r1 },
        DecompressedBlock { batch_serial, bytes: r2 },
        DecompressedBlock { batch_serial, bytes: other },
    ))
}

/// Paired-split FASTQ encode step: [`DecodedRecordBatch`] → three
/// [`DecompressedBlock`]s of FASTQ text (R1, R2, other), one per output file.
///
/// See [`encode_paired_batch`] for the per-batch routing/encode logic; this
/// function only wires that logic into the `Process3WithWorkerState` step
/// shape.
#[allow(clippy::type_complexity)]
pub(crate) fn build_fastq_paired_encode_step(
    opts: Arc<FastqOptions>,
    per_step_byte_limit: u64,
    records_written: Arc<AtomicU64>,
) -> Process3WithWorkerState<
    DecodedRecordBatch,
    DecompressedBlock,
    DecompressedBlock,
    DecompressedBlock,
    FastqEncodeState,
    impl Fn() -> FastqEncodeState + Send + Sync + 'static,
    impl Fn(
        &mut FastqEncodeState,
        DecodedRecordBatch,
    ) -> io::Result<Process3Output<DecompressedBlock, DecompressedBlock, DecompressedBlock>>
    + Send
    + Sync
    + 'static,
> {
    process3_with_worker_state::<
        DecodedRecordBatch,
        DecompressedBlock,
        DecompressedBlock,
        DecompressedBlock,
        FastqEncodeState,
        _,
        _,
    >(
        "FastqEncodePaired",
        per_step_byte_limit,
        per_step_byte_limit,
        per_step_byte_limit,
        new_encode_state,
        move |state: &mut FastqEncodeState, batch: DecodedRecordBatch| {
            encode_paired_batch(state, batch, &opts, &records_written)
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paired_encode_step_profile_is_parallel_three_branches() {
        use crate::commands::fastq::FastqOptions;
        use crate::pipeline::core::step::{Step, StepKind};
        use std::sync::Arc;
        use std::sync::atomic::AtomicU64;

        let opts = Arc::new(FastqOptions {
            exclude_flags: 0,
            require_flags: 0,
            no_suffix: true,
            umi_header: None,
        });
        let step = build_fastq_paired_encode_step(opts, 1 << 20, Arc::new(AtomicU64::new(0)));
        let profile = step.profile();
        assert_eq!(profile.kind, StepKind::Parallel);
        assert_eq!(profile.output_queues.len(), 3);
        assert_eq!(profile.branch_ordering.len(), 3);
    }

    /// Build a `DecodedRecord` for a single-segment (unpaired-flags-wise
    /// simple) paired read: `PAIRED` plus exactly one of `FIRST_SEGMENT` /
    /// `LAST_SEGMENT`, matching how `classify_segment` routes R1 vs R2.
    fn decoded_record(name: &str, first_segment: bool) -> fgumi_bam_io::DecodedRecord {
        use fgumi_bam_io::{DecodedRecord, GroupKey};
        use fgumi_raw_bam::{SamBuilder as RawSamBuilder, flags as raw_flags};

        let mut flags = raw_flags::PAIRED;
        flags |= if first_segment { raw_flags::FIRST_SEGMENT } else { raw_flags::LAST_SEGMENT };
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags)
            .ref_id(0)
            .pos(99)
            .mapq(60);
        DecodedRecord::from_raw_bytes(b.build(), GroupKey::default())
    }

    fn no_filter_opts() -> FastqOptions {
        FastqOptions { exclude_flags: 0, require_flags: 0, no_suffix: true, umi_header: None }
    }

    /// Direct, deterministic regression test for the no-ordinal-gaps rule
    /// (see the doc comment on `encode_paired_batch`): a batch containing
    /// ONLY R1 records must still return `Some(empty block)` on the R2 and
    /// "other" branches, never `None`. Unlike the integration-level
    /// `test_fastq_paired_empty_r2_batch_no_hang` (which depends on the
    /// fixture fitting in a single BGZF block to avoid a real hang), this
    /// calls the routing/encode logic directly and asserts on `Process3Output`
    /// itself, so it fails immediately — no subprocess, no timing, no
    /// block-size heuristics — if a regression ever changes `b`/`c` to `None`.
    #[test]
    fn encode_paired_batch_emits_empty_block_not_none_when_branch_has_no_records() {
        use std::sync::atomic::AtomicU64;

        let opts = no_filter_opts();
        let mut state = new_encode_state();
        let records_written = AtomicU64::new(0);

        let batch = DecodedRecordBatch::new(0, vec![decoded_record("read1", true)]);
        let out = encode_paired_batch(&mut state, batch, &opts, &records_written)
            .expect("encode_paired_batch should succeed");

        let r1 = out.a.expect("R1 branch must be Some: this batch has an R1 record");
        assert!(!r1.bytes.is_empty(), "R1 branch must contain the encoded R1 record");

        let r2 = out.b.expect(
            "R2 branch must be Some (empty block) even though this batch has no R2 \
             records -- returning None here is exactly the regression that wedges \
             the R2 reorder stage",
        );
        assert!(r2.bytes.is_empty(), "R2 branch must be an empty block, not populated");

        let other = out.c.expect(
            "other branch must be Some (empty block) even though this batch has no \
             single-end/ambiguous records",
        );
        assert!(other.bytes.is_empty(), "other branch must be an empty block, not populated");
    }

    /// Sanity check for the routing side of `encode_paired_batch`: a batch
    /// with one record on each branch produces non-empty bytes on all three,
    /// so the "empty means absent" assertions above aren't trivially true for
    /// every branch regardless of content.
    #[test]
    fn encode_paired_batch_routes_mixed_batch_to_all_three_branches() {
        use std::sync::atomic::AtomicU64;

        let opts = no_filter_opts();
        let mut state = new_encode_state();
        let records_written = AtomicU64::new(0);

        // "Other" = neither/both segment flags set; use an unpaired single-end read.
        let mut other_builder = fgumi_raw_bam::SamBuilder::new();
        other_builder
            .read_name(b"single1")
            .sequence(b"GGGGCCCC")
            .qualities(&[30; 8])
            .flags(0)
            .ref_id(0)
            .pos(199)
            .mapq(60);
        let other_record = fgumi_bam_io::DecodedRecord::from_raw_bytes(
            other_builder.build(),
            fgumi_bam_io::GroupKey::default(),
        );

        let batch = DecodedRecordBatch::new(
            0,
            vec![decoded_record("read1", true), decoded_record("read1", false), other_record],
        );
        let out = encode_paired_batch(&mut state, batch, &opts, &records_written)
            .expect("encode_paired_batch should succeed");

        assert!(out.a.is_some_and(|b| !b.bytes.is_empty()), "R1 branch should be non-empty");
        assert!(out.b.is_some_and(|b| !b.bytes.is_empty()), "R2 branch should be non-empty");
        assert!(out.c.is_some_and(|b| !b.bytes.is_empty()), "other branch should be non-empty");
    }
}
