//! `MetricsSink` — a terminal pipeline step that records per-coordinate-group
//! QC metrics (simplex/duplex family sizes, UMI counts, downsampling yield)
//! into a shared per-thread accumulator, emitting nothing downstream.
//!
//! This is the sink half of routing the standalone `simplex-metrics` /
//! `duplex-metrics` commands onto the typed-step pipeline. The parallel BAM
//! decode + `GroupByMi` front end (assembled by
//! [`ChainBuilder`](crate::pipeline::chains::builder::ChainBuilder)) feeds this
//! step `BatchedMiGroups`, and each worker records its batch into its own
//! `ConsensusMetricsSlot` via the SAME machinery the consensus commands'
//! metrics-on path uses (`push_mi_group_entries` → `split_into_runs` →
//! `classify_batch_runs` → `BoundaryReorder::submit`). The
//! [`ConsensusMetricsFinalizeHook`](crate::inline_metrics_collector::ConsensusMetricsFinalizeHook)
//! merges the per-thread slots and writes the byte-identical TSV files after
//! the pipeline drains.
//!
//! # Why `Outputs = ()`
//!
//! Metrics are a pure side effect recorded into the shared accumulator; there
//! is no BAM to write. A sink step (`type Outputs = ()`) terminates the chain
//! without a placeholder output, so the command's `SinkSpec::None` skips the
//! usual `BgzfCompress → WriteBgzfFile` tail entirely.
//!
//! # Parallelism and correctness
//!
//! `StepKind::Parallel`: every worker holds its own clone (an `Arc` of the
//! shared, mutex-sharded `ConsensusMetricsCaptures`). A coordinate group is the
//! atomic unit — `classify_batch_runs` records interior (fully-within-one-batch)
//! runs directly and defers batch-boundary runs to the shared `BoundaryReorder`,
//! which closes them in batch-serial order as the prefix becomes contiguous, so
//! a family split across a batch boundary is never double- or under-counted.
//! This is the exact contract the consensus T2 path relies on.

use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use crate::commands::shared_metrics::{consensus_guard_record, ensure_not_consensus_record};
use crate::inline_metrics_collector::{ConsensusMetricsCaptures, record_batch_metrics};
use crate::pipeline::core::item::HeapSize;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::group::mi::BatchedMiGroups;

/// Terminal metrics-recording step. Consumes [`BatchedMiGroups`], records each
/// batch into the shared per-thread accumulator, and produces no output.
///
/// Construct one per run and hand it to `PipelineBuilder::append_step`; the
/// framework clones it per worker via [`Step::new_worker_copy`] (every clone
/// shares the same `Arc`s, so all workers write into the one mutex-sharded
/// accumulator and the one boundary reorder).
pub(crate) struct MetricsSink {
    /// Shared per-thread accumulator + boundary reorder + intervals + output
    /// prefix. Cloned per worker (the `Arc` is shared; the shard mutexes inside
    /// keep per-worker writes contention-free in the common case).
    captures: Arc<ConsensusMetricsCaptures>,
    /// Input BAM header, used to resolve reference names for `TemplateInfo`.
    header: Arc<noodles::sam::Header>,
    /// Read-group → library index, for partitioning families by library.
    library_index: Arc<fgumi_bam_io::LibraryIndex>,
    /// Input path, for the consensus-BAM guard's error message.
    input: Arc<PathBuf>,
    /// Whether this worker has already run the consensus-BAM guard. Mirrors the
    /// serial `process_templates_from_bam`'s `consensus_checked` flag: the guard
    /// only needs to fire on the first qualifying record this worker sees (a
    /// consensus BAM carries the tags on every record), so it is a one-shot
    /// check, not a per-batch cost.
    consensus_checked: bool,
    name: &'static str,
}

impl MetricsSink {
    /// Build a `MetricsSink` from the shared captures + input header/library
    /// index. All three are `Arc`s so the per-worker clones are cheap and share
    /// the same backing state.
    #[must_use]
    pub(crate) fn new(
        captures: Arc<ConsensusMetricsCaptures>,
        header: Arc<noodles::sam::Header>,
        library_index: Arc<fgumi_bam_io::LibraryIndex>,
        input: Arc<PathBuf>,
    ) -> Self {
        Self {
            captures,
            header,
            library_index,
            input,
            consensus_checked: false,
            name: "MetricsSink",
        }
    }

    /// Record one batch of MI groups into this worker's accumulator slot.
    ///
    /// Mirrors the consensus metrics-on batch body
    /// (`run_simplex_consensus_batch_with_metrics`) exactly, minus the consensus
    /// calling: collect one `(TemplateInfo, ReadInfoKey)` entry per qualifying
    /// template across every family in the batch, split into maximal same-key
    /// runs, record interior runs directly, and submit boundary runs to the
    /// shared reorder (recording whatever groups that submission closes). The
    /// two `with_slot` acquisitions never nest around the reorder mutex (the
    /// reorder is taken/released inside `submit`, between them), so there is no
    /// lock-order cycle.
    fn record_batch(&mut self, item: BatchedMiGroups) -> io::Result<()> {
        let BatchedMiGroups { batch_serial, groups } = item;

        // Reject a consensus BAM the same way the serial
        // `process_templates_from_bam` does: on the first qualifying record this
        // worker sees, rather than by re-opening the input (which would break the
        // stdin/pipe inputs the metrics commands accept). The guard runs on the
        // grouped (post-`consensus_pregroup_keep_raw`) records, so it fires on a
        // mapped consensus BAM — the realistic misuse of feeding the aligned
        // consensus output back into a metrics tool.
        if !self.consensus_checked
            && let Some(guard) = groups.iter().find_map(|g| consensus_guard_record(&g.records))
        {
            ensure_not_consensus_record(guard, &self.input).map_err(io::Error::other)?;
            self.consensus_checked = true;
        }

        // The interior/boundary recording body is shared with the three consensus
        // metrics-on batch runners; see `record_batch_metrics`.
        record_batch_metrics(&self.captures, &self.header, &self.library_index, batch_serial, &groups)
    }
}

impl HeapSize for MetricsSink {}

impl Step for MetricsSink {
    type Input = BatchedMiGroups;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            // Parallel: every worker records its own batches into its own
            // accumulator slot; the group-atomicity contract (interior runs +
            // boundary reorder) makes the sharded accumulation order-independent.
            kind: StepKind::Parallel,
            sticky: false,
            // A sink has no output branches.
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(batch) => {
                self.record_batch(batch)?;
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }

    fn new_worker_copy(&self) -> Self {
        // Every per-worker copy shares the same `Arc`s: the mutex-sharded
        // accumulator distributes per-worker writes across its slots, and the
        // single `BoundaryReorder` serializes cross-batch group closing.
        Self {
            captures: Arc::clone(&self.captures),
            header: Arc::clone(&self.header),
            library_index: Arc::clone(&self.library_index),
            input: Arc::clone(&self.input),
            // Each worker re-runs the one-shot guard on its own first batch.
            consensus_checked: false,
            name: self.name,
        }
    }
}
