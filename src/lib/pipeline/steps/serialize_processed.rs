//! `build_serialize_processed_groups_step`: shared builder for the
//! `BatchedProcessedPositionGroups → DecompressedBlock` serialize step.
//!
//! Standalone `fgumi group` and the runall `--stop-after group` fused
//! chain (`execute_group_pipeline` in `commands::runall`) both end on
//! the same shape: take an MI-id-assigned batch of position groups,
//! rewrite the MI tag bytes into each record's body, frame each body
//! with the BAM `block_size` prefix, and concatenate into a single
//! `DecompressedBlock` ready for `BgzfCompress`.
//!
//! Doing MI-tag-injection + framing in ONE closure pass (rather than
//! the three-step `BatchedProcessedPositionGroups → BatchedMiGroups →
//! RecordBatch → DecompressedBlock` consensus-tap detour) eliminates
//! two cross-step queue handoffs and two intermediate batch allocations
//! per BAM block. Measured ~11–15 % wall-time recovery on `--stop-after
//! group` benchmarks vs. the consensus-tap detour the runall path was
//! using before this builder existed.
//!
//! `emit_templates_raw_with_mi` lives in this module because both this
//! step AND the standalone single-threaded writer in `commands::group`
//! call it — keeping it pub(crate) at the step-library layer means
//! both callers can't drift in MI-injection semantics.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use log::info;

use crate::pipeline::steps::group::position::BatchedProcessedPositionGroups;
use crate::pipeline::steps::process::{ProcessWithWorkerState, process_with_worker_state};
use crate::pipeline::steps::types::DecompressedBlock;
use crate::template::Template;

/// Per-worker scratch state for the serialize step. `scratch` is the
/// record-aligned tag-update buffer; `mi_buf` is the formatted MI
/// string. Both are reused across batches for the same reason the
/// legacy `try_step_serialize` reused
/// `worker.core.{serialization_buffer, secondary_serialization_buffer}`.
///
/// Must be `pub` because [`build_serialize_processed_groups_step`]'s
/// return type names it via `ProcessWithWorkerState<…, SerializeState, …>` —
/// callers that bind the returned step to a `let` can't see a private
/// type in the public signature.
pub struct SerializeState {
    scratch: Vec<u8>,
    mi_buf: String,
}

impl crate::pipeline::core::item::HeapSize for SerializeState {}

/// Construct the `SerializeGroups` step: consumes
/// [`BatchedProcessedPositionGroups`], emits [`DecompressedBlock`].
///
/// The per-batch closure:
///   * walks every `ProcessedPositionGroup.templates`,
///   * calls [`emit_templates_raw_with_mi`] which writes the MI tag
///     bytes into the per-record `scratch` and appends framed bytes
///     into the batch's output `Vec<u8>`,
///   * updates `progress_counter` with the input record count and
///     logs a "Processed N records" line every million records.
///
/// The output `Vec<u8>` is sized at `total_records * 800` so a typical
/// BAM record (~256 B body + ~4 B framing + tag-update headroom) hits
/// zero `realloc`s mid-batch. Smaller starts (e.g. 64 KiB + grow)
/// landed in the medium → large mimalloc size class boundary and
/// regressed on CODEC 8M benches.
///
/// `assign_tag_bytes` is the 2-byte tag the MI string is written under
/// (`MI` by default; user-overridable via `--assign-tag`).
///
/// `progress_counter` is shared across all worker copies — bump
/// atomically with the batch's input record count.
pub fn build_serialize_processed_groups_step(
    per_step_byte_limit: u64,
    assign_tag_bytes: [u8; 2],
    progress_counter: Arc<AtomicU64>,
) -> ProcessWithWorkerState<
    BatchedProcessedPositionGroups,
    DecompressedBlock,
    impl Fn(&mut SerializeState, BatchedProcessedPositionGroups) -> io::Result<DecompressedBlock>
    + Send
    + Sync
    + 'static,
    SerializeState,
    impl Fn() -> SerializeState + Send + Sync + 'static,
> {
    process_with_worker_state::<
        BatchedProcessedPositionGroups,
        DecompressedBlock,
        _,
        SerializeState,
        _,
    >(
        "SerializeGroups",
        per_step_byte_limit,
        || SerializeState { scratch: Vec::with_capacity(512), mi_buf: String::with_capacity(16) },
        move |state: &mut SerializeState,
              item: BatchedProcessedPositionGroups|
              -> io::Result<DecompressedBlock> {
            let BatchedProcessedPositionGroups { batch_serial, groups } = item;
            // Size by record count, not template count: a paired-end template
            // holds ~2 records, so `templates.len()` under-sized the buffer ~2×
            // and forced reallocs. `input_record_count` is the per-group record
            // tally; summed once here for both the capacity hint and the
            // progress counter below.
            let total_input_records: u64 =
                groups.iter().fold(0u64, |acc, g| acc.saturating_add(g.input_record_count));
            let mut output = Vec::with_capacity(
                usize::try_from(total_input_records.saturating_mul(800)).unwrap_or(usize::MAX),
            );
            for processed in &groups {
                emit_templates_raw_with_mi(
                    &processed.templates,
                    0, // already mi_assign'd upstream — base offset is 0
                    assign_tag_bytes,
                    &mut state.scratch,
                    &mut state.mi_buf,
                    |bytes| {
                        output.extend_from_slice(bytes);
                        Ok(())
                    },
                )?;
            }

            let prev = progress_counter.fetch_add(total_input_records, Ordering::Relaxed);
            if (prev + total_input_records) / 1_000_000 > prev / 1_000_000 {
                info!("Processed {} records", prev + total_input_records);
            }

            Ok(DecompressedBlock { batch_serial, bytes: output })
        },
    )
}

/// Emit raw-byte primary reads for a batch of templates with MI tags injected.
///
/// For each template, if `has_mi` then the raw bytes are copied into `scratch`,
/// the MI tag is updated in place, and the result is emitted via `emit`; otherwise
/// the raw bytes are emitted unchanged. Each emitted record is prefixed with a
/// 4-byte little-endian `block_size`. `mi_buf` and `scratch` are reused across
/// templates to amortize allocations.
///
/// Shared by the threaded-pipeline `SerializeGroups` step (this module's
/// [`build_serialize_processed_groups_step`]) and the single-threaded writer
/// in [`crate::commands::group`] so they can't drift in MI-injection
/// semantics.
pub(crate) fn emit_templates_raw_with_mi(
    templates: &[Template],
    base_mi: u64,
    assign_tag_bytes: [u8; 2],
    scratch: &mut Vec<u8>,
    mi_buf: &mut String,
    mut emit: impl FnMut(&[u8]) -> io::Result<()>,
) -> io::Result<()> {
    for template in templates {
        let mi = template.mi;
        let has_mi = mi.is_assigned();
        if has_mi {
            // write_with_offset clears the buffer before writing.
            mi.write_with_offset(base_mi, mi_buf);
        }
        for raw in [template.r1(), template.r2()].into_iter().flatten() {
            if has_mi {
                scratch.clear();
                scratch.extend_from_slice(raw);
                fgumi_raw_bam::update_string_tag(scratch, assign_tag_bytes, mi_buf.as_bytes());
                let block_size = u32::try_from(scratch.len()).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("BAM record too large ({} bytes)", scratch.len()),
                    )
                })?;
                emit(&block_size.to_le_bytes())?;
                emit(scratch)?;
            } else {
                let block_size = u32::try_from(raw.len()).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("BAM record too large ({} bytes)", raw.len()),
                    )
                })?;
                emit(&block_size.to_le_bytes())?;
                emit(raw)?;
            }
        }
    }
    Ok(())
}
