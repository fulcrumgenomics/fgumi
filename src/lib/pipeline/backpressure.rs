//! Per-stage backpressure thresholds and the queue-memory-budget/high-water
//! mark relationship.
//!
//! Moved out of `unified_pipeline::base` (C5/R6a) — the definitions are pure
//! (a `u64` constant and a `fn(u64, u64) -> u64`) with no dependency on the
//! `unified_pipeline` scheduler/queue machinery, so they relocate cleanly to
//! the root crate's `pipeline` module. `unified_pipeline::base` keeps a
//! `pub use` shim so its own internal call sites and the still-present legacy
//! command paths keep resolving these names unchanged.

/// Default high-water mark for a reorder-buffered stage (512 MiB).
///
/// This is the point at which one stage stops taking on new work and the
/// scheduler starts prioritizing downstream steps (Process, Serialize,
/// Compress, Write) so the stage can drain.
///
/// It is reached only through `ReorderBufferState::effective_limit`, so it
/// applies to exactly three stages: Q2 and Q3 (whose counters span the queue
/// and its reorder buffer both) and the post-Group write reorder buffer. The
/// remaining queues — Q1, Q2b, Q4, Q6, Q7 — carry byte counters but no byte
/// mark of their own; they are bounded by their slot count and by the
/// pipeline-wide total gated at the Read step.
///
/// # Architecture — where the scheduler reads its drain signal
///
/// Only one of the three feeds the scheduler's `memory_high` / `memory_drained`
/// pair, and the two pipelines pick different ones:
/// - **BAM pipeline**: `q3_reorder_state`, the reorder buffer after Decode and
///   before Group (`BamPipelineState::is_memory_high`)
/// - **FASTQ pipeline**: `output.write_reorder_state`, the post-Compress write
///   reorder buffer (`FastqStepContext::get_backpressure`)
///
/// For the BAM pipeline that placement is deliberate: memory is tracked before
/// the exclusive Group step rather than after it, because tracking after Group
/// releases the pre-Group buffer's bytes before knowing whether the post-Group
/// queue can accept them, leaving data in an untracked intermediate buffer.
///
/// # Threshold behaviour
///
/// - When tracked memory >= the mark, `is_memory_high()` returns true: the
///   producing step declines new work and the scheduler enters "drain mode"
/// - When tracked memory < half the mark, `is_memory_drained()` returns true
///   and the scheduler exits drain mode (hysteresis prevents thrashing)
///
/// # Relationship to `--max-memory`
///
/// This is a *per-stage trigger*, not the user's capacity budget; see
/// [`stage_high_water_mark`] for why the two are separate and why a larger
/// budget does not raise this. The budget's capacity role is enforced on the
/// pipeline's total in-flight bytes at the Read step
/// (`BamPipelineState::read_admission_allowed`).
pub const BACKPRESSURE_THRESHOLD_BYTES: u64 = 512 * 1024 * 1024; // 512 MiB

/// Default high-water mark for Q5, the processed queue (256 MiB).
///
/// This is set lower than the Q3 mark (256 MiB vs 512 MiB) because items in Q5
/// are typically larger (e.g., `SimplexProcessedBatch` with `RecordBuf` vectors).
/// When Q5 memory reaches this mark, the Process step pauses to let
/// downstream steps (Serialize, Compress, Write) catch up.
///
/// Like [`BACKPRESSURE_THRESHOLD_BYTES`] it is applied through
/// [`stage_high_water_mark`], so a `--max-memory` below it tightens it and one
/// above it leaves it alone.
pub const Q5_BACKPRESSURE_THRESHOLD_BYTES: u64 = 256 * 1024 * 1024; // 256 MiB

/// The high-water mark one pipeline stage backs off at, given the declared
/// queue memory budget and that stage's default mark.
///
/// Returns `default_mark` when `queue_memory_limit` is 0 (meaning "unset"),
/// and `min(queue_memory_limit, default_mark)` otherwise.
///
/// # Why the budget cannot raise the mark (issue #765)
///
/// `--max-memory` is a *capacity* budget: how many bytes the pipeline may hold
/// in total. A stage's high-water mark is a *trigger*: the depth at which that
/// one stage stops pulling in work so the rest of the pipeline can catch up.
/// Historically the mark was the only place the budget reached, which fused the
/// two ideas into one number and made `--max-memory` above 512 MiB inert. The
/// capacity half now lives where it belongs — on the pipeline's total in-flight
/// bytes, gated at the Read step — and what is left here is only the trigger.
///
/// A trigger should not scale with the total budget, for three reasons:
///
/// - **It would buy no throughput.** The scheduler's response to `memory_high`
///   is a priority reordering, not a stop. Raising the mark only postpones the
///   moment output steps are favoured, so the stage parks more bytes and the
///   pipeline runs no faster.
/// - **Reorder capacity is bounded by serial skew, not by bytes.** The
///   pre-Group reorder buffers exist to absorb out-of-order completion between
///   workers. That skew is bounded by the number of batches in flight, which is
///   bounded by the queue slot count, so bytes past the skew window are dead
///   weight.
/// - **It would let one stage monopolize the budget.** These marks are the only
///   byte bound on Q2, Q3 and Q5 individually. Uncapped, a single stage could
///   park the entire budget and the total gate would fire before any other
///   stage buffered anything — converting a pipeline into a queue.
///
/// The downward direction is the useful one and is preserved: a budget smaller
/// than the default mark pulls the mark down with it, so no single stage can
/// hold more than the whole declared budget.
#[inline]
#[must_use]
pub fn stage_high_water_mark(queue_memory_limit: u64, default_mark: u64) -> u64 {
    if queue_memory_limit == 0 { default_mark } else { queue_memory_limit.min(default_mark) }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ========================================================================
    // Stage high-water mark tests (issue #765)
    // ========================================================================

    #[test]
    fn test_stage_high_water_mark_treats_zero_as_unset() {
        assert_eq!(
            stage_high_water_mark(0, BACKPRESSURE_THRESHOLD_BYTES),
            BACKPRESSURE_THRESHOLD_BYTES
        );
        assert_eq!(
            stage_high_water_mark(0, Q5_BACKPRESSURE_THRESHOLD_BYTES),
            Q5_BACKPRESSURE_THRESHOLD_BYTES
        );
    }

    #[test]
    fn test_stage_high_water_mark_tightens_under_a_small_budget() {
        // A budget below the default mark pulls the mark down with it, so a
        // single stage can never park more than the whole declared budget.
        assert_eq!(
            stage_high_water_mark(64 * 1024 * 1024, BACKPRESSURE_THRESHOLD_BYTES),
            64 * 1024 * 1024
        );
        assert_eq!(
            stage_high_water_mark(64 * 1024 * 1024, Q5_BACKPRESSURE_THRESHOLD_BYTES),
            64 * 1024 * 1024
        );
    }

    #[test]
    fn test_stage_high_water_mark_does_not_scale_above_the_default() {
        // Deliberate: the mark is a per-stage trigger, not the user's capacity
        // budget. See the `stage_high_water_mark` docs for why raising it buys
        // nothing. The budget's capacity role lives at the Read admission gate.
        let huge = 32 * 1024 * 1024 * 1024;
        assert_eq!(
            stage_high_water_mark(huge, BACKPRESSURE_THRESHOLD_BYTES),
            BACKPRESSURE_THRESHOLD_BYTES
        );
        assert_eq!(
            stage_high_water_mark(huge, Q5_BACKPRESSURE_THRESHOLD_BYTES),
            Q5_BACKPRESSURE_THRESHOLD_BYTES
        );
    }
}
