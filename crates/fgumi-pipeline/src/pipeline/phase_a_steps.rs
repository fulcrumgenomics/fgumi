//! Phase A `try_step` functions for work-stealing pool stages.
//!
//! Each function attempts one unit of work for its stage. Returns a [`StepResult`] to
//! indicate success, output-full (item held in worker state), input-empty, or error.
//! The worker loop calls these in priority order determined by the scheduler.
//!
//! **Pool stages:**
//! - [`try_step_extract`]: Read one FASTQ template, produce an ordinal-tagged
//!   `(ordinal, Template, Vec<u8>)` triple. Multiple workers can compute in parallel
//!   after the FASTQ source mutex is released.
//! - [`try_step_correct`]: Pop an ordinal-tagged triple, correct the UMI, push the
//!   corrected triple to `q_corrected` (parallelizable, per-worker LRU caches).
//! - [`try_step_to_fastq`]: Pop an ordinal-tagged triple, insert into the reorder buffer,
//!   drain consecutive items in ordinal order to aligner stdin + `q_unmapped`.

use std::io::Write;
use std::sync::atomic::Ordering;

use anyhow::anyhow;
use fgumi_lib::correct;
use fgumi_lib::extract::extract_template_and_fastq;
use fgumi_lib::unified_pipeline::MemoryEstimate;
use noodles::sam;

use crate::pipeline::phase_a_graph::{OrdinalPair, PhaseAGraph, ToFastqState};
use crate::pipeline::phase_a_worker::PhaseAWorkerState;
use crate::pipeline::queue::WorkQueue;
use crate::pipeline::scheduler::StepResult;

// ============================================================================
// Extract
// ============================================================================

/// Attempt one Extract step: read one FASTQ template, produce an ordinal-tagged
/// `(ordinal, Template, Vec<u8>)` triple.
///
/// The FASTQ source mutex is held only during the read (serial for gzipped input).
/// After reading, the mutex is released and the compute-intensive
/// `extract_template_and_fastq` runs in parallel. The ordinal assigned at read time
/// tracks the original FASTQ order; the `ToFastq` reorder buffer restores this order
/// before writing to the aligner stdin.
///
/// On success, pushes the ordinal-tagged triple to `q_extracted`.
/// If the output queue is full, the triple is held in worker state for retry.
pub fn try_step_extract(graph: &PhaseAGraph, worker: &mut PhaseAWorkerState) -> StepResult {
    // 1. Try to push any held pair from a previous iteration.
    if let Some(pair) = worker.held_extract_pair.take() {
        let mem = worker.held_extract_pair_mem;
        if let Some(returned) = graph.q_extracted.push(pair, mem) {
            worker.held_extract_pair = Some(returned);
            return StepResult::OutputFull;
        }
        worker.held_extract_pair_mem = 0;
    }

    // 2. Check if extract is already done.
    if graph.extract_done.load(Ordering::Acquire) {
        return StepResult::InputEmpty;
    }

    // 3. Try to acquire the FASTQ source lock (non-blocking).
    let mut source_guard = match graph.extract_source.try_lock() {
        Ok(guard) => guard,
        Err(std::sync::TryLockError::WouldBlock) => return StepResult::InputEmpty,
        Err(std::sync::TryLockError::Poisoned(e)) => {
            return StepResult::Error(anyhow!("extract source lock poisoned: {e}"));
        }
    };

    // 4. Read one FASTQ read set (assigns a monotonic ordinal).
    let Some((ordinal, read_set)) = source_guard.next_read_set() else {
        // Source exhausted. Mark extract done (closes q_extracted).
        drop(source_guard);
        graph.mark_extract_done();
        return StepResult::InputEmpty;
    };
    // Release the lock before doing compute work — other workers can read the next
    // item while this worker runs extract_template_and_fastq in parallel.
    drop(source_guard);

    // 5. Extract template + FASTQ bytes (CPU-intensive, runs without any lock).
    let extracted = match extract_template_and_fastq(
        &read_set,
        graph.quality_encoding,
        &graph.extract_params,
    ) {
        Ok(e) => e,
        Err(e) => return StepResult::Error(e),
    };

    let template = extracted.template;
    let fastq_bytes = extracted.fastq_bytes;

    // 6. Push the ordinal-tagged triple to q_extracted.
    //    Items may arrive out of ordinal order due to parallel compute; the
    //    ToFastq reorder buffer will restore ordinal order downstream.
    let pair_mem = template.estimate_heap_size() + fastq_bytes.len();
    if let Some(returned) = graph.q_extracted.push((ordinal, template, fastq_bytes), pair_mem) {
        // Queue full — hold the triple for retry.
        worker.held_extract_pair = Some(returned);
        worker.held_extract_pair_mem = pair_mem;
        return StepResult::OutputFull;
    }

    StepResult::Success
}

// ============================================================================
// Correct
// ============================================================================

/// Attempt one Correct step: pop an ordinal-tagged `(ordinal, Template, Vec<u8>)` triple
/// from `q_extracted`, correct the UMI on the Template, and push the corrected triple
/// to `q_corrected`.
///
/// Correct is parallelizable — multiple workers can run it concurrently, each with
/// their own LRU cache. The ordinal and FASTQ bytes pass through unchanged.
///
/// When `q_extracted` is closed and empty, closes `q_corrected` to propagate the
/// close cascade to `ToFastq`.
pub fn try_step_correct(graph: &PhaseAGraph, worker: &mut PhaseAWorkerState) -> StepResult {
    let Some(ref correction_config) = graph.correction_config else {
        return StepResult::InputEmpty;
    };
    let Some(ref q_corrected) = graph.q_corrected else {
        return StepResult::InputEmpty;
    };

    // 1. Try to push any held pair from a previous iteration.
    if let Some(pair) = worker.held_correct_pair.take() {
        let mem = worker.held_correct_pair_mem;
        if let Some(returned) = q_corrected.push(pair, mem) {
            worker.held_correct_pair = Some(returned);
            return StepResult::OutputFull;
        }
        worker.held_correct_pair_mem = 0;
        graph.q_extracted.complete_in_flight();
        graph.q_extracted.release_memory(mem);
    }

    // 2. Pop from q_extracted.
    let Some((ordinal, mut template, fastq_bytes)) = graph.q_extracted.try_pop() else {
        if graph.q_extracted.is_closed_and_empty() {
            graph.mark_correct_done();
        }
        return StepResult::InputEmpty;
    };

    // 3. Extract UMI from the first record and correct it.
    let umi_tag = sam::alignment::record::data::field::Tag::new(
        correction_config.umi_tag[0],
        correction_config.umi_tag[1],
    );

    let umi_opt = template.records.first().and_then(|rec| {
        rec.data().get(&umi_tag).and_then(|v| {
            if let sam::alignment::record_buf::data::field::Value::String(s) = v {
                Some(String::from_utf8_lossy(s).to_string())
            } else {
                None
            }
        })
    });

    if let Some(umi) = umi_opt {
        let correction = correct::compute_template_correction(
            &umi,
            correction_config.umi_length,
            false, // revcomp not exposed in pipeline
            correction_config.max_mismatches,
            correction_config.min_distance_diff,
            &correction_config.encoded_umi_set,
            &mut worker.umi_cache,
        );

        if correction.matched {
            // Apply correction to all records in the template.
            if correction.needs_correction {
                for record in &mut template.records {
                    if correction.has_mismatches {
                        record.data_mut().insert(
                            sam::alignment::record::data::field::Tag::ORIGINAL_UMI_BARCODE_SEQUENCE,
                            sam::alignment::record_buf::data::field::Value::String(
                                correction.original_umi.clone().into(),
                            ),
                        );
                    }
                    if let Some(ref corrected) = correction.corrected_umi {
                        record.data_mut().insert(
                            umi_tag,
                            sam::alignment::record_buf::data::field::Value::String(
                                corrected.clone().into(),
                            ),
                        );
                    }
                }
            }
        } else {
            // Unmatched — drop the template (pipeline doesn't write rejects).
            graph.q_extracted.complete_in_flight();
            let mem = template.estimate_heap_size() + fastq_bytes.len();
            graph.q_extracted.release_memory(mem);
            return StepResult::Success;
        }
    }
    // If no UMI tag, pass template through unchanged.

    // 4. Push corrected triple (with ordinal) to q_corrected.
    let pair_mem = template.estimate_heap_size() + fastq_bytes.len();
    if let Some(returned) = q_corrected.push((ordinal, template, fastq_bytes), pair_mem) {
        worker.held_correct_pair = Some(returned);
        worker.held_correct_pair_mem = pair_mem;
        return StepResult::OutputFull;
    }

    // 5. Complete in-flight tracking on q_extracted.
    graph.q_extracted.complete_in_flight();
    graph.q_extracted.release_memory(pair_mem);

    StepResult::Success
}

// ============================================================================
// ToFastq
// ============================================================================

/// Attempt one `ToFastq` step: pop an ordinal-tagged triple from the input queue,
/// insert it into the reorder buffer, drain consecutive items in ordinal order,
/// and write them to the aligner stdin + push to `q_unmapped`.
///
/// Workers pop from the input queue *without* holding the `tofastq_state` lock. This
/// allows multiple workers to pop concurrently. The lock is then acquired to insert
/// into the reorder buffer, drain ready items, and perform the serial I/O (stdin write
/// and `q_unmapped` push) in ordinal order.
pub fn try_step_to_fastq(graph: &PhaseAGraph, _worker: &mut PhaseAWorkerState) -> StepResult {
    let input_q = graph.tofastq_input_queue();

    let Some((ordinal, template, fastq_bytes)) = input_q.try_pop() else {
        if input_q.is_closed_and_empty() {
            // A worker may have popped the last item and inserted it into the reorder
            // buffer but not yet drained it. Drain any stragglers before closing.
            if let Ok(mut state) = graph.tofastq_state.lock() {
                if let Err(e) = drain_reorder_buffer(graph, &mut state, input_q) {
                    return StepResult::Error(e);
                }
            }
            graph.close_aligner_stdin();
            graph.q_unmapped.close();
        }
        return StepResult::InputEmpty;
    };

    let mut state = match graph.tofastq_state.lock() {
        Ok(guard) => guard,
        Err(e) => {
            return StepResult::Error(anyhow!("tofastq_state lock poisoned: {e}"));
        }
    };

    state.reorder.insert(ordinal, (template, fastq_bytes));

    if let Err(e) = drain_reorder_buffer(graph, &mut state, input_q) {
        return StepResult::Error(e);
    }

    StepResult::Success
}

/// Drain consecutive items from the reorder buffer, writing FASTQ to stdin and
/// pushing Templates to `q_unmapped` in ordinal order.
///
/// Called while holding the `tofastq_state` lock.
fn drain_reorder_buffer(
    graph: &PhaseAGraph,
    state: &mut ToFastqState,
    input_q: &WorkQueue<OrdinalPair>,
) -> anyhow::Result<()> {
    while let Some((template, fastq_bytes)) = state.reorder.try_next() {
        let template_mem = template.estimate_heap_size();
        let pair_mem = template_mem + fastq_bytes.len();

        match write_fastq_to_stdin(state, &fastq_bytes) {
            Ok(()) | Err(WriteResult::StdinClosed) => {}
            Err(WriteResult::Error(e)) => {
                input_q.complete_in_flight();
                input_q.release_memory(pair_mem);
                return Err(e);
            }
        }

        if let Some(ref bs) = graph.batch_state {
            let num_bytes = fastq_bytes.len() as u64;
            #[expect(
                clippy::naive_bytecount,
                reason = "small buffer; adding bytecount dep not worthwhile"
            )]
            let newline_count = fastq_bytes.iter().filter(|&&b| b == b'\n').count() as u64;
            let num_records = newline_count / 4;
            // bwa -K counts bases, not bytes; paired-end FASTQ has ~2x bytes per base.
            bs.record_send(num_bytes / 2, num_records);
        }

        // Spin-push to q_unmapped — must complete before draining the next item
        // so that q_unmapped order matches stdin write order.
        let mut current_template = template;
        loop {
            match graph.q_unmapped.push(current_template, template_mem) {
                None => break,
                Some(returned) => {
                    if graph.cancel.load(Ordering::Acquire) {
                        input_q.complete_in_flight();
                        input_q.release_memory(pair_mem);
                        return Err(anyhow!("cancelled while pushing to q_unmapped"));
                    }
                    current_template = returned;
                    std::thread::yield_now();
                }
            }
        }

        input_q.complete_in_flight();
        input_q.release_memory(pair_mem);
    }

    Ok(())
}

/// Result of attempting to write FASTQ bytes to aligner stdin.
enum WriteResult {
    /// Aligner stdin has been closed (normal during shutdown).
    StdinClosed,
    /// I/O error writing to aligner stdin.
    Error(anyhow::Error),
}

/// Write FASTQ bytes to the aligner's stdin pipe via the `ToFastqState`.
fn write_fastq_to_stdin(state: &mut ToFastqState, data: &[u8]) -> Result<(), WriteResult> {
    let Some(stdin) = state.stdin.as_mut() else {
        return Err(WriteResult::StdinClosed);
    };

    stdin.write_all(data).map_err(|e| {
        if e.kind() == std::io::ErrorKind::BrokenPipe {
            WriteResult::StdinClosed
        } else {
            WriteResult::Error(anyhow!("failed to write FASTQ to aligner stdin: {e}"))
        }
    })
}

// ============================================================================
// Zipper (dedicated thread — not a pool stage)
// ============================================================================

/// Run the Zipper merge as a dedicated thread.
///
/// Sequentially pops one Template from `q_unmapped` and one from `q_mapped`, verifies
/// their querynames match (guaranteed by the ordering invariant), merges tags from the
/// unmapped into the mapped template, encodes the result to raw BAM bytes, and pushes
/// them to `q_zippered`.
///
/// Returns the total number of raw BAM records pushed.
///
/// # Errors
///
/// Returns an error if merge fails, encoding fails, names don't match (ordering
/// violation), or the mapped stream ends before the unmapped stream.
///
/// ## Ordering guarantee
///
/// `q_unmapped` is populated by the exclusive `ToFastq` stage in stdin write order.
/// `q_mapped` is populated by the stdout-reader in aligner output order. Since bwa
/// preserves input order within `-K` batches, these two queues arrive in the same
/// queryname order, enabling simple sequential matching without a `HashMap`.
pub fn zipper_thread_loop(graph: &PhaseAGraph) -> anyhow::Result<u64> {
    use fgumi_lib::progress::ProgressTracker;
    use fgumi_lib::vendored::encode_record_buf;

    let progress = ProgressTracker::new("Zipper merged records").with_interval(1_000_000);

    // Wait for the output header (set by stdout-reader after reading SAM header).
    let output_header = loop {
        if let Some(h) = graph.output_header.get() {
            break h;
        }
        if graph.cancel.load(Ordering::Acquire) {
            anyhow::bail!("zipper cancelled while waiting for output header");
        }
        // Also check if q_unmapped is already done (no templates to merge).
        if graph.q_unmapped.is_closed_and_empty() {
            log::info!(
                "zipper: q_unmapped closed before output header available, nothing to merge"
            );
            return Ok(0);
        }
        std::thread::sleep(std::time::Duration::from_millis(1));
    };

    let mut record_count: u64 = 0;

    loop {
        // Pop unmapped (blocking).
        let Some(unmapped) = graph.q_unmapped.pop() else {
            break; // q_unmapped closed+empty → done
        };
        graph.q_unmapped.complete_in_flight();

        // Pop mapped (blocking).
        let Some(mut mapped) = graph.q_mapped.pop() else {
            anyhow::bail!(
                "Zipper: mapped stream ended before unmapped. \
                 Unmapped template '{}' has no corresponding mapped template.",
                String::from_utf8_lossy(&unmapped.name)
            );
        };
        graph.q_mapped.complete_in_flight();

        // Verify names match — ordering guarantee from combined q_extracted queue.
        if unmapped.name != mapped.name {
            anyhow::bail!(
                "Zipper ordering violation: unmapped='{}' vs mapped='{}'. \
                 This indicates a bug in the pipeline ordering logic.",
                String::from_utf8_lossy(&unmapped.name),
                String::from_utf8_lossy(&mapped.name),
            );
        }

        // Merge unmapped metadata into the mapped template.
        (graph.merge_fn)(&unmapped, &mut mapped, &graph.tag_info, graph.skip_pa_tags)?;

        // Encode all records from the merged template to raw BAM bytes.
        let header = output_header.as_ref();
        for rec in mapped.all_reads() {
            let mut buf = Vec::new();
            encode_record_buf(&mut buf, header, rec)
                .map_err(|e| anyhow!("failed to encode merged record to BAM: {e}"))?;
            let mem = buf.len();
            graph.q_zippered.push_blocking(buf, mem)?;
            record_count += 1;
            progress.log_if_needed(1);
        }

        // Release memory for the unmapped template.
        let unmapped_mem = unmapped.estimate_heap_size();
        graph.q_unmapped.release_memory(unmapped_mem);
        let mapped_mem = mapped.estimate_heap_size();
        graph.q_mapped.release_memory(mapped_mem);
    }

    // Check for leftover mapped reads.
    if let Some(remaining) = graph.q_mapped.pop() {
        graph.q_mapped.complete_in_flight();
        anyhow::bail!(
            "Zipper: processed all unmapped reads but mapped reads remain. \
             Found template '{}'. Ensure unmapped and mapped inputs have matching read names.",
            String::from_utf8_lossy(&remaining.name)
        );
    }

    progress.log_final();
    Ok(record_count)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_result_variants() {
        // Ensure the enum compiles and can be matched.
        let _closed = WriteResult::StdinClosed;
        let _err = WriteResult::Error(anyhow!("test"));
    }
}
