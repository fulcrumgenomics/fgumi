//! Phase A graph: queues and shared state for the alignment pipeline stages.
//!
//! The graph holds all inter-stage queues, shared configuration, and the mutex-protected
//! FASTQ source for Extract. Pool workers reference this graph immutably; mutable per-worker
//! state lives in [`PhaseAWorkerState`](super::phase_a_worker::PhaseAWorkerState).
//!
//! **Pool stages (work-stealing):** Extract, `ToFastq`.
//! **Dedicated threads (not in pool):** stdout-reader (SAM → `q_mapped`),
//! zipper (merge `q_unmapped` + `q_mapped` → `q_zippered`),
//! sort-accumulate (`q_zippered` → `SortedChunks`).

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, OnceLock};

use fgumi_lib::correct::EncodedUmiSet;
use fgumi_lib::extract::{ExtractParams, QualityEncoding};
use fgumi_lib::fastq::FastqSet;
use fgumi_lib::template::Template;
use fgumi_umi::TagInfo;
use noodles::sam::Header;

use crate::pipeline::batch_state::AlignerBatchState;
use crate::pipeline::queue::WorkQueue;
use crate::pipeline::reorder::ReorderBuffer;

/// Type alias for the zipper merge function.
///
/// The merge function transfers tags from the unmapped template to the mapped template.
/// This is a function pointer (not a closure) because the merge logic lives in the binary
/// crate (`commands::zipper::merge`) and is injected into the graph at construction time.
pub type MergeFn = fn(&Template, &mut Template, &TagInfo, bool) -> anyhow::Result<()>;

/// Type alias for an ordinal-tagged `(Template, FASTQ bytes)` pair.
///
/// The ordinal (`u64`) is assigned by [`ExtractSource`] at FASTQ read time and tracks the
/// original read order. This enables parallel Extract compute — workers can push to
/// `q_extracted` in any order, and the [`ToFastqState`] reorder buffer restores ordinal
/// order before writing to the aligner stdin.
pub type OrdinalPair = (u64, Template, Vec<u8>);

// ============================================================================
// ToFastqState
// ============================================================================

/// Shared state for the `ToFastq` stage: reorder buffer and aligner stdin pipe.
///
/// Protected by a single mutex in [`PhaseAGraph`]. Workers pop items from the input
/// queue *without* holding this lock, then acquire it to insert into the reorder buffer,
/// drain consecutive items, and write them to the aligner stdin in ordinal order.
///
/// This replaces the previous exclusive `tofastq_lock` approach: multiple workers can
/// pop from the input queue concurrently, and the reorder buffer restores the correct
/// FASTQ read order for the aligner.
pub struct ToFastqState {
    /// Reorder buffer that reassembles items in ordinal order.
    pub(super) reorder: ReorderBuffer<(Template, Vec<u8>)>,
    /// Aligner stdin writer. Set to `None` when the stdin pipe is closed (extract done).
    pub(super) stdin: Option<std::process::ChildStdin>,
}

// ============================================================================
// Queue capacity and memory defaults
// ============================================================================

/// Default slot capacity for Template queues (unmapped, corrected, mapped).
const TEMPLATE_QUEUE_CAPACITY: usize = 4096;

/// Default slot capacity for FASTQ bytes and zippered raw BAM bytes queues.
const BYTES_QUEUE_CAPACITY: usize = 4096;

/// Default memory limit per queue (256 MiB).
const DEFAULT_QUEUE_MEMORY_LIMIT: usize = 256 * 1024 * 1024;

// ============================================================================
// ExtractSource
// ============================================================================

/// Mutex-protected iterator over combined FASTQ read sets.
///
/// Extract is inherently serial for gzipped input, so only one worker runs it at a time
/// via `try_lock()`. Other workers that fail to acquire the lock move on to other stages.
pub struct ExtractSource {
    /// The FASTQ read-set iterator producing one [`FastqSet`] per template.
    /// `None` once the source is exhausted.
    inner: Option<Box<dyn Iterator<Item = FastqSet> + Send>>,
    /// Monotonic ordinal counter, incremented on each successful read.
    /// Ordinals track the original FASTQ read order so that downstream stages
    /// (via the reorder buffer) can restore this order after parallel compute.
    next_ordinal: u64,
}

impl ExtractSource {
    /// Create a new extract source from the given iterator.
    #[must_use]
    pub fn new(iter: Box<dyn Iterator<Item = FastqSet> + Send>) -> Self {
        Self { inner: Some(iter), next_ordinal: 0 }
    }

    /// Try to read the next FASTQ read set with its ordinal.
    ///
    /// Returns `Some((ordinal, read_set))` on success, `None` when exhausted.
    /// The ordinal is a monotonically increasing counter starting at 0, assigned
    /// in the order items are read from the FASTQ source.
    pub fn next_read_set(&mut self) -> Option<(u64, FastqSet)> {
        let iter = self.inner.as_mut()?;
        let item = iter.next();
        if let Some(read_set) = item {
            let ordinal = self.next_ordinal;
            self.next_ordinal += 1;
            Some((ordinal, read_set))
        } else {
            self.inner = None;
            None
        }
    }

    /// Returns `true` if the source has been exhausted.
    #[must_use]
    pub fn is_exhausted(&self) -> bool {
        self.inner.is_none()
    }
}

// ============================================================================
// CorrectionConfig
// ============================================================================

/// Configuration for in-pipeline UMI correction.
///
/// When present in the graph, the Correct stage is active and uses `q_corrected`
/// as an intermediate queue between Extract and `ToFastq`.
pub struct CorrectionConfig {
    /// Pre-encoded set of known UMI sequences for matching.
    pub encoded_umi_set: Arc<EncodedUmiSet>,
    /// UMI length (all segments must match this length).
    pub umi_length: usize,
    /// Maximum mismatches per UMI segment.
    pub max_mismatches: usize,
    /// Minimum distance difference between best and second-best match.
    pub min_distance_diff: usize,
    /// UMI tag name as 2-byte array (e.g. `*b"RX"`).
    pub umi_tag: [u8; 2],
    /// LRU cache capacity per worker thread (0 = no cache).
    pub cache_size: usize,
}

// ============================================================================
// PhaseAGraph
// ============================================================================

/// Phase A graph: queues and shared state for alignment pipeline stages.
///
/// Shared immutably across all pool workers. The graph connects two pool stages:
///
/// - **Extract**: reads FASTQ, produces `(Template, Vec<u8>)` pairs (template + FASTQ bytes)
/// - **`ToFastq`**: pops pairs from `q_extracted`, writes FASTQ bytes to aligner stdin,
///   pushes Templates to `q_unmapped` (maintaining insertion order = stdin order)
///
/// Dedicated threads handle the rest:
/// - **stdout-reader**: SAM from aligner → `q_mapped`
/// - **zipper**: joins `q_unmapped` + `q_mapped` → merge → encode → `q_zippered`
/// - **sort-accumulate**: `q_zippered` → `SortedChunks`
///
/// ## Ordering guarantee
///
/// Each FASTQ read set is assigned a monotonic ordinal by [`ExtractSource`] at read
/// time. After parallel Extract compute, items may arrive in `q_extracted` out of
/// ordinal order. The `ToFastq` stage uses a [`ReorderBuffer`] (inside
/// [`tofastq_state`](Self::tofastq_state)) to reassemble items in ordinal order
/// before writing to the aligner stdin and pushing to `q_unmapped`. Since bwa
/// preserves input order within `-K` batches, `q_unmapped` and `q_mapped` arrive
/// in the same queryname order, enabling simple sequential matching in the Zipper.
pub struct PhaseAGraph {
    // ---- Queues ----
    /// Ordinal-tagged `(ordinal, Template, FASTQ bytes)` triples from Extract.
    ///
    /// Items may arrive out of FASTQ-read order because Extract compute runs in parallel
    /// after the source mutex is released. The ordinal tracks the original read order so
    /// that the [`ToFastqState`] reorder buffer can restore it before writing to stdin.
    pub q_extracted: WorkQueue<OrdinalPair>,
    /// Corrected ordinal-tagged triples (feeds `ToFastq`).
    /// Only present when UMI correction is enabled.
    pub q_corrected: Option<WorkQueue<OrdinalPair>>,
    /// Unmapped Templates, populated by `ToFastq` after writing FASTQ to stdin.
    /// Order matches aligner stdin order (feeds the dedicated Zipper thread).
    pub q_unmapped: WorkQueue<Template>,
    /// Mapped templates from the aligner stdout reader thread (feeds Zipper).
    pub q_mapped: WorkQueue<Template>,
    /// Raw BAM bytes from Zipper merge (feeds the sort-accumulate thread).
    pub q_zippered: WorkQueue<Vec<u8>>,

    // ---- Extract state ----
    /// Shared FASTQ source, protected by mutex (Extract is serial for gzipped input).
    pub extract_source: Mutex<ExtractSource>,
    /// Extract configuration (tags, read group, read structures, etc.).
    pub extract_params: ExtractParams,
    /// Detected quality encoding for FASTQ input.
    pub quality_encoding: QualityEncoding,
    /// Set to `true` when the FASTQ source is exhausted.
    pub extract_done: AtomicBool,

    // ---- ToFastq reorder state ----
    /// Shared state for the `ToFastq` stage: reorder buffer + aligner stdin pipe.
    ///
    /// Workers pop ordinal-tagged items from the input queue without holding this lock,
    /// then acquire it to insert into the reorder buffer, drain consecutive items in
    /// ordinal order, write FASTQ bytes to stdin, and push Templates to `q_unmapped`.
    /// This replaces the previous exclusive `tofastq_lock` approach.
    pub tofastq_state: Mutex<ToFastqState>,

    // ---- Correction config ----
    /// UMI correction configuration. When `Some`, the Correct pool step is active.
    pub correction_config: Option<CorrectionConfig>,

    // ---- Zipper / output header ----
    /// Output SAM header used to encode merged records to BAM bytes.
    /// Set by the stdout-reader thread once the SAM header has been read from the aligner.
    pub output_header: OnceLock<Arc<Header>>,
    /// Tag manipulation info for the zipper merge (remove/reverse/revcomp).
    pub tag_info: Arc<TagInfo>,
    /// Whether to skip PA (primary alignment) tag generation in zipper.
    pub skip_pa_tags: bool,
    /// Zipper merge function: transfers tags from unmapped to mapped template.
    pub merge_fn: MergeFn,

    // ---- Aligner batch state (optional, when bwa -K is used) ----
    /// Batch coordination state for bwa `-K` flow control.
    pub batch_state: Option<Arc<AlignerBatchState>>,

    // ---- Cancellation ----
    /// Shared cancellation flag.
    pub cancel: Arc<AtomicBool>,
}

impl PhaseAGraph {
    /// Create a new Phase A graph with the given configuration.
    ///
    /// Queues are created with default capacities and memory limits.
    // The graph struct needs all its dependencies passed in at construction.
    #[expect(clippy::too_many_arguments, reason = "graph requires all dependencies")]
    #[must_use]
    pub fn new(
        extract_source: ExtractSource,
        extract_params: ExtractParams,
        quality_encoding: QualityEncoding,
        tag_info: Arc<TagInfo>,
        skip_pa_tags: bool,
        merge_fn: MergeFn,
        batch_state: Option<Arc<AlignerBatchState>>,
        aligner_stdin: std::process::ChildStdin,
        cancel: Arc<AtomicBool>,
        correction_config: Option<CorrectionConfig>,
    ) -> Self {
        let q_corrected = correction_config.as_ref().map(|_| {
            WorkQueue::new(TEMPLATE_QUEUE_CAPACITY, DEFAULT_QUEUE_MEMORY_LIMIT, "q_corrected")
        });

        Self {
            q_extracted: WorkQueue::new(
                TEMPLATE_QUEUE_CAPACITY,
                DEFAULT_QUEUE_MEMORY_LIMIT,
                "q_extracted",
            ),
            q_corrected,
            q_unmapped: WorkQueue::new(
                TEMPLATE_QUEUE_CAPACITY,
                DEFAULT_QUEUE_MEMORY_LIMIT,
                "q_unmapped",
            ),
            q_mapped: WorkQueue::new(
                TEMPLATE_QUEUE_CAPACITY,
                DEFAULT_QUEUE_MEMORY_LIMIT,
                "q_mapped",
            ),
            q_zippered: WorkQueue::new(
                BYTES_QUEUE_CAPACITY,
                DEFAULT_QUEUE_MEMORY_LIMIT,
                "q_zippered",
            ),
            extract_source: Mutex::new(extract_source),
            extract_params,
            quality_encoding,
            extract_done: AtomicBool::new(false),
            tofastq_state: Mutex::new(ToFastqState {
                reorder: ReorderBuffer::new(),
                stdin: Some(aligner_stdin),
            }),
            correction_config,
            output_header: OnceLock::new(),
            tag_info,
            skip_pa_tags,
            merge_fn,
            batch_state,
            cancel,
        }
    }

    /// Returns the input queue for the `ToFastq` stage.
    ///
    /// When correction is enabled, this is `q_corrected` (Correct pops from `q_extracted`
    /// and pushes to `q_corrected`). When correction is disabled, this is `q_extracted`
    /// directly (Extract → `ToFastq` with no intermediate queue).
    #[must_use]
    pub fn tofastq_input_queue(&self) -> &WorkQueue<OrdinalPair> {
        self.q_corrected.as_ref().unwrap_or(&self.q_extracted)
    }

    /// Sample backpressure signals from queue state.
    #[must_use]
    pub fn sample_backpressure(&self) -> super::phase_a::PhaseABackpressure {
        super::phase_a::PhaseABackpressure {
            extracted_high: self.q_extracted.fill_ratio() > super::BACKPRESSURE_HIGH_WATERMARK,
            tofastq_input_high: self.tofastq_input_queue().fill_ratio()
                > super::BACKPRESSURE_HIGH_WATERMARK,
            can_send_fastq: self.batch_state.as_ref().is_none_or(|bs| bs.can_send()),
        }
    }

    /// Returns `true` if all pool work is complete.
    ///
    /// The pool handles Extract, Correct (when enabled), and `ToFastq`. Work is complete
    /// when the `ToFastq` input queue is closed and empty (all pairs have been split
    /// by `ToFastq` and their Templates pushed to `q_unmapped`).
    #[must_use]
    pub fn is_pool_done(&self) -> bool {
        self.tofastq_input_queue().is_closed_and_empty()
    }

    /// Mark extract as done and close `q_extracted`.
    pub fn mark_extract_done(&self) {
        self.extract_done.store(true, Ordering::Release);
        self.q_extracted.close();
    }

    /// Close `q_corrected` after all extracted pairs have been corrected.
    ///
    /// Called by the Correct step when `q_extracted` is closed and empty.
    /// Only meaningful when correction is enabled.
    pub fn mark_correct_done(&self) {
        if let Some(ref q) = self.q_corrected {
            q.close();
        }
    }

    /// Close the aligner stdin pipe (signals EOF to the aligner).
    ///
    /// Called by `ToFastq` when the input queue is drained. Dropping the stdin handle
    /// closes the pipe and signals EOF to the aligner subprocess.
    pub fn close_aligner_stdin(&self) {
        if let Ok(mut guard) = self.tofastq_state.lock() {
            guard.stdin = None;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_source_exhaustion() {
        // We can't easily construct a real FastqSet, so test the wrapper logic
        // with a None-initialized source.
        let mut source = ExtractSource { inner: None, next_ordinal: 0 };
        assert!(source.is_exhausted());
        assert!(source.next_read_set().is_none());
    }

    #[test]
    fn test_backpressure_default() {
        let bp = super::super::phase_a::PhaseABackpressure::default();
        assert!(!bp.extracted_high);
        assert!(!bp.tofastq_input_high);
        assert!(bp.can_send_fastq);
    }
}
