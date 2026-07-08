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

use crate::commands::fastq::{FastqOptions, write_fastq_record};
use crate::logging::OperationTimer;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::steps::process::{ProcessWithWorkerState, process_with_worker_state};
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
