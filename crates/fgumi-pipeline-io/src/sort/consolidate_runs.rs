//! `ConsolidateRuns` — the merge of a consolidation job, the head of the
//! consolidation tail (`SpillReadPlanner` branch 2 → `ConsolidateRuns` →
//! `SpillBlockCompressConsolidation` → `WriteConsolidatedRun`).
//!
//! A `Detached` step on the `sort-io` driver (where `SpillWrite` waits for the
//! job). It merges one [`ConsolidationJob`]'s slots — whose reads the planner
//! plans and whose blocks `SortSpillDecompress` decompresses, exactly as for
//! the final merge — with a kind-dispatched merge driver
//! ([`fgumi_sort::merge_driver_for_kind`]), frames the merged records in the
//! spill framing, cuts them into raw blocks with the gather's rule and emits
//! each as a [`SpillBlockEvent::Block`] (`file_id` = the job id) for the pool to
//! compress. The output is therefore framed and compressed exactly as a run the
//! spill writer wrote.
//!
//! The merge is cooperative: a bounded number of records (and record bytes)
//! per `try_run`. On a stall it publishes the awaited slot on the consolidation
//! demand (whose deliveries unpark this driver), re-checks, decompresses the
//! slot's already-read blocks itself, and otherwise reports `Contention`. It
//! does not request a pool worker: supply reaches it through the planner's hot
//! set and the decompress step's directed notify.
//!
//! Ordinals are dense across jobs (one counter on the step), so the ordered
//! edge into the compressor never waits on an ordinal that will not come.

use std::collections::VecDeque;
use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{CounterSpec, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use fgumi_sort::{
    MergeCtxErased, MergeDemand, MergeDriverDyn, MergeStepFramed, QnClassTable, SpillKeyKind,
};

use super::merge::SelfServe;
use super::protocol::{ConsolidationJob, SpillBlockEvent};
use super::spill_block_compress::SpillBufPool;
use super::supply_ledger::SupplyLedger;

/// Counter slot: records merged.
const CONSOLIDATION_RECORDS: usize = 0;

/// Records merged per `try_run`. Large enough that per-call overhead is
/// negligible for short reads.
const RECORDS_PER_CALL: usize = 16 * 1024;

/// Record bytes merged per `try_run`, whichever of the two budgets is reached
/// first: long reads make a record budget alone unbounded in time, and a
/// `Detached` step must report progress promptly.
const BYTES_PER_CALL: usize = 8 * 1024 * 1024;

/// The job being merged.
struct Running {
    job_id: u32,
    key_kind: SpillKeyKind,
    driver: Box<dyn MergeDriverDyn + Send>,
    /// The job's slots, to name a failed input by its path.
    slots: Vec<Arc<fgumi_sort::SortMergeSlot>>,
    input_paths: Vec<PathBuf>,
    /// The block being framed.
    raw: Vec<u8>,
}

/// `Detached` consolidation merge (see the module docs).
pub struct ConsolidateRuns {
    running: Option<Running>,
    /// Dense across jobs.
    next_ordinal: u64,
    demand: Arc<MergeDemand>,
    class_table: Option<Arc<QnClassTable>>,
    block_size: usize,
    output_byte_limit: u64,
    outbox: VecDeque<SpillBlockEvent>,
    held: HeldSlot<Unpushed<SpillBlockEvent>>,
    records_per_call: usize,
    bytes_per_call: usize,
    /// The self-serve path's decompressor, built on first use.
    decomp: Option<fgumi_sort::SpillBlockDecompressor>,
    ledger: Option<Arc<SupplyLedger>>,
    window_budget: u64,
    records_merged: u64,
    /// Recycled block buffers, shared with the tail's compressor and writer.
    buf_pool: Option<Arc<SpillBufPool>>,
}

impl ConsolidateRuns {
    /// A step framing raw blocks of about `block_size` bytes (the gather's cut
    /// rule) into an output edge bounded at `output_byte_limit` bytes. Its
    /// demand is a fresh one until [`Self::with_merge_demand`].
    #[must_use]
    pub fn new(block_size: usize, output_byte_limit: u64) -> Self {
        Self {
            running: None,
            next_ordinal: 0,
            demand: Arc::new(MergeDemand::new()),
            class_table: None,
            block_size: block_size.max(1),
            output_byte_limit,
            outbox: VecDeque::new(),
            held: HeldSlot::new(),
            records_per_call: RECORDS_PER_CALL,
            bytes_per_call: BYTES_PER_CALL,
            decomp: None,
            ledger: None,
            window_budget: super::spill_decompress::reorder_window_budget(output_byte_limit),
            records_merged: 0,
            buf_pool: None,
        }
    }

    /// Frame blocks into buffers recycled through `pool`: the tail's compressor
    /// returns each raw block to it and its writer each compressed one, as the
    /// phase-1 spill path does, instead of allocating per block.
    #[must_use]
    pub fn with_buf_pool(mut self, pool: Arc<SpillBufPool>) -> Self {
        self.buf_pool = Some(pool);
        self
    }

    /// An empty block buffer with room for a full block and its overflow
    /// record: a recycled one when the pool has one.
    fn fresh_block(&self) -> Vec<u8> {
        let mut buf = self.buf_pool.as_ref().map(|p| p.take()).unwrap_or_default();
        buf.reserve(self.block_size + 1024);
        buf
    }

    /// The consolidation demand: shared with the read planner (its awaited and
    /// predicted slots are hot) and the decompress step (whose deliveries to
    /// the awaited slot unpark this step's driver).
    #[must_use]
    pub fn with_merge_demand(mut self, demand: Arc<MergeDemand>) -> Self {
        self.demand = demand;
        self
    }

    /// The sort's queryname class table, for queryname jobs' merge keys.
    #[must_use]
    pub fn with_class_table(mut self, table: Arc<QnClassTable>) -> Self {
        self.class_table = Some(table);
        self
    }

    /// The supply ledger the self-serve path releases stash charges to.
    #[must_use]
    pub fn with_supply_ledger(mut self, ledger: Arc<SupplyLedger>) -> Self {
        self.ledger = Some(ledger);
        self
    }

    /// Bound one `try_run` to `records` records and `bytes` record-body bytes.
    #[must_use]
    pub fn with_per_call(mut self, records: usize, bytes: usize) -> Self {
        self.records_per_call = records.max(1);
        self.bytes_per_call = bytes.max(1);
        self
    }

    /// The class table handed to queryname jobs (tests).
    #[cfg(any(test, feature = "test-utils"))]
    #[doc(hidden)]
    #[must_use]
    pub fn class_table_for_test(&self) -> Option<&Arc<QnClassTable>> {
        self.class_table.as_ref()
    }

    /// Records merged so far (tests).
    #[cfg(test)]
    pub(crate) fn records_merged_for_test(&self) -> u64 {
        self.records_merged
    }

    /// Queue a block of the running job.
    fn emit(&mut self, job_id: u32, key_kind: SpillKeyKind, bytes: Vec<u8>, last: bool) {
        self.outbox.push_back(SpillBlockEvent::Block {
            ordinal: self.next_ordinal,
            file_id: job_id,
            key_kind,
            is_last_in_file: last,
            records_ingested_so_far: 0,
            bytes,
        });
        self.next_ordinal += 1;
    }

    /// Start `job`'s merge.
    fn start(&mut self, job: ConsolidationJob) {
        let ConsolidationJob { job_id, slots, key_kind, input_paths } = job;
        let ctx = match (&self.class_table, key_kind) {
            (Some(t), SpillKeyKind::QuerynameLex | SpillKeyKind::QuerynameNatural) => {
                MergeCtxErased::Queryname(Arc::clone(t))
            }
            _ => MergeCtxErased::None,
        };
        let driver = fgumi_sort::merge_driver_for_kind(
            key_kind,
            slots.clone(),
            ctx,
            Arc::clone(&self.demand),
        );
        self.running =
            Some(Running { job_id, key_kind, driver, slots, input_paths, raw: self.fresh_block() });
    }

    /// Map a merge error to the step's error: a truncated input names its path.
    fn merge_error(run: &Running, e: &anyhow::Error) -> io::Error {
        let truncated =
            e.downcast_ref::<io::Error>().is_some_and(|e| e.kind() == io::ErrorKind::UnexpectedEof);
        if truncated
            && let Some(slot) = run.driver.failed_slot()
            && let Some(i) = run.slots.iter().position(|s| Arc::ptr_eq(s, slot))
        {
            return io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("consolidation input {} is truncated", run.input_paths[i].display()),
            );
        }
        io::Error::other(format!("ConsolidateRuns: merge step failed: {e:#}"))
    }

    /// Merge up to the per-call budgets of the running job. Returns the step
    /// outcome; the job's blocks are queued in `outbox`.
    fn advance(&mut self) -> io::Result<StepOutcome> {
        let mut run = self.running.take().expect("a running job");
        let key_bytes = run.key_kind.serialized_key_bytes();
        let (mut records, mut bytes) = (0usize, 0usize);
        while records < self.records_per_call && bytes < self.bytes_per_call {
            let before = run.raw.len();
            match run.driver.try_step_framed(&mut run.raw) {
                Err(e) => return Err(Self::merge_error(&run, &e)),
                Ok(MergeStepFramed::Produced) => {
                    records += 1;
                    bytes += run.raw.len() - before - 4 - key_bytes;
                    // The gather's cut rule (`frame_one_block`): a record that
                    // pushes a non-empty block past the size starts the next
                    // block; a block reaching the size closes, a lone oversized
                    // record forming its own.
                    if before != 0 && run.raw.len() > self.block_size {
                        // The overflow record moves to a full-capacity buffer
                        // (not a `split_off` sized to the record, which the
                        // next block would regrow from scratch).
                        let mut next = self.fresh_block();
                        next.extend_from_slice(&run.raw[before..]);
                        run.raw.truncate(before);
                        let block = std::mem::replace(&mut run.raw, next);
                        self.emit(run.job_id, run.key_kind, block, false);
                    }
                    if run.raw.len() >= self.block_size {
                        let block = std::mem::replace(&mut run.raw, self.fresh_block());
                        self.emit(run.job_id, run.key_kind, block, false);
                    }
                }
                Ok(MergeStepFramed::Stalled) => {
                    let Some(slot) = run.driver.stalled_slot().map(Arc::clone) else {
                        break;
                    };
                    // Declare the slot awaited and re-check (`MergeDemand`'s
                    // wake protocol): a block or EOF that already landed means
                    // keep merging.
                    if self.demand.await_slot(&slot) {
                        continue;
                    }
                    let mut serve = SelfServe {
                        decomp: &mut self.decomp,
                        ledger: self.ledger.as_deref(),
                        window_budget: self.window_budget,
                    };
                    let (served, ready) = serve.run(&slot)?;
                    self.demand.stats().record_self_serve(served, ready);
                    if ready {
                        self.demand.clear_awaited();
                        continue;
                    }
                    self.records_merged += records as u64;
                    self.running = Some(run);
                    return Ok(if records > 0 || !self.outbox.is_empty() {
                        StepOutcome::Progress
                    } else {
                        StepOutcome::Contention
                    });
                }
                Ok(MergeStepFramed::Done) => {
                    // The final block, empty when the job ended on a full block
                    // or had no records, so the writer always closes the file.
                    let block = std::mem::take(&mut run.raw);
                    self.emit(run.job_id, run.key_kind, block, true);
                    self.records_merged += records as u64;
                    return Ok(StepOutcome::Progress);
                }
            }
        }
        self.records_merged += records as u64;
        self.running = Some(run);
        Ok(StepOutcome::Progress)
    }
}

impl Step for ConsolidateRuns {
    type Input = ConsolidationJob;
    type Outputs = OrderedBytesSingle<SpillBlockEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ConsolidateRuns",
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_IO_GROUP)
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("consolidation_records", "records")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        let mut pushed = false;
        while let Some(ev) = self.outbox.pop_front() {
            pushed = true;
            if let Err(u) = ctx.outputs.push(ev) {
                self.held.put(u);
                return Ok(StepOutcome::Progress);
            }
        }
        if self.running.is_none() {
            // Check the input before popping (an empty poll is recorded).
            let Some(job) = ctx.input.pop() else {
                return Ok(if pushed {
                    StepOutcome::Progress
                } else if ctx.input.is_drained() {
                    StepOutcome::Finished
                } else {
                    StepOutcome::NoProgress
                });
            };
            self.start(job);
        }
        let before = self.records_merged;
        let outcome = self.advance()?;
        ctx.counters.add(CONSOLIDATION_RECORDS, self.records_merged - before);
        while let Some(ev) = self.outbox.pop_front() {
            if let Err(u) = ctx.outputs.push(ev) {
                self.held.put(u);
                return Ok(StepOutcome::Progress);
            }
        }
        Ok(if pushed { StepOutcome::Progress } else { outcome })
    }
}

/// The planner's consolidation branch when the chain has no consolidation tail
/// (a `--max-temp-files` limit that cannot trigger): `SpillWrite` then starts
/// no job, so nothing arrives here, and a job that did would be a wiring bug,
/// reported rather than dropped.
#[derive(Debug, Default)]
pub struct NoConsolidation;

impl Step for NoConsolidation {
    type Input = ConsolidationJob;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "NoConsolidation",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(job) = ctx.input.pop() {
            return Err(io::Error::other(format!(
                "consolidation job {} reached a sort chain built without a consolidation tail",
                job.job_id
            )));
        }
        Ok(if ctx.input.is_drained() { StepOutcome::Finished } else { StepOutcome::NoProgress })
    }
}

#[cfg(test)]
mod tests;
