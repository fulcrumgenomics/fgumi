//! `GroupByQueryname` step for the new typed-step pipeline framework.
//!
//! Mirrors the legacy pipeline's queryname grouping (`fgumi filter
//! --filter-by-template=true`, `fgumi sort --order queryname`):
//! consumes `DecodedRecordBatch` items, run-length groups consecutive
//! same-name records into `Template`s, and emits **batches** of
//! completed templates as one `BamTemplateBatch` per output ordinal.
//!
//! Design follows `GroupByMi` and `GroupByPosition`: `Serial` +
//! `ByItemOrdinal`, batched output to amortize downstream Serial
//! mutex cost, held-slot retry pattern for backpressure.

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BamTemplateBatch, DecodedRecordBatch};
use crate::template::Template;

/// Max input batches consumed per `try_run` invocation. Amortizes the
/// Serial mutex acquisition; mirrors `GroupByMi::MAX_BATCHES_PER_LOCK`.
const MAX_BATCHES_PER_LOCK: usize = 8;

/// Default target batch count. Mirrors legacy's
/// `TemplateGrouper::new(1000)` (the production batch size in
/// `Filter::execute_threads_mode_template`).
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 1000;

/// `Serial + ByItemOrdinal` queryname grouper. Records arriving in
/// queryname-grouped (or queryname-sorted) order are partitioned into
/// runs of consecutive same-name records; each run becomes one
/// [`Template`]. Completed templates accumulate until the target batch
/// size is reached, then emit as one [`BamTemplateBatch`].
///
/// State is held behind the framework's per-step mutex (the runtime
/// stores a `Serial` step inside `Arc<Mutex<...>>` and acquires
/// per-`try_run`). Held-slot retry pattern matches `GroupByMi`:
/// when an emit is rejected by downstream backpressure, the partial
/// state stays buffered until the held slot drains.
pub struct GroupByQueryname {
    /// Run-length state: queryname currently being accumulated.
    current_name: Option<Vec<u8>>,
    /// Run-length state: cached `name_hash` for fast pre-check.
    current_name_hash: Option<u64>,
    /// Run-length state: records accumulated for `current_name`.
    current_records: Vec<fgumi_raw_bam::RawRecord>,
    /// Self-managed monotonic ordinal — each emitted *batch* gets the
    /// next value.
    next_ordinal: u64,
    /// Accumulator: completed templates waiting to be packaged.
    accumulator: Vec<Template>,
    /// Held output slot when downstream rejected the most recent push.
    held: HeldSlot<Unpushed<BamTemplateBatch>>,
    target_batch_count: usize,
    output_byte_limit: u64,
    name: &'static str,
}

impl GroupByQueryname {
    /// Construct a `GroupByQueryname` step.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self::with_target_batch_count(output_byte_limit, DEFAULT_TARGET_BATCH_COUNT)
    }

    /// Construct with a custom target batch count.
    #[must_use]
    pub fn with_target_batch_count(output_byte_limit: u64, target_batch_count: usize) -> Self {
        Self {
            current_name: None,
            current_name_hash: None,
            current_records: Vec::new(),
            next_ordinal: 0,
            accumulator: Vec::with_capacity(target_batch_count),
            held: HeldSlot::new(),
            target_batch_count: target_batch_count.max(1),
            output_byte_limit,
            name: "GroupByQueryname",
        }
    }

    /// Flush the run-in-progress into the accumulator (if non-empty).
    fn flush_current_template(&mut self) -> io::Result<()> {
        if !self.current_records.is_empty() {
            let records = std::mem::take(&mut self.current_records);
            let template = Template::from_records(records)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            self.accumulator.push(template);
            self.current_name = None;
            self.current_name_hash = None;
        }
        Ok(())
    }

    /// Package the accumulator into a `BamTemplateBatch` and emit. On
    /// rejection, hold; the next `try_run` will retry via the held slot.
    fn emit_batch(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let serial = self.next_ordinal;
        self.next_ordinal += 1;
        let templates =
            std::mem::replace(&mut self.accumulator, Vec::with_capacity(self.target_batch_count));
        let out = BamTemplateBatch::new(serial, templates);
        if let Err(unpushed) = ctx.outputs.push(out) {
            self.held.put(unpushed);
        }
        StepOutcome::Progress
    }
}

impl Step for GroupByQueryname {
    type Input = DecodedRecordBatch;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held slot first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. If the accumulator is full enough, emit a batch first.
        if self.accumulator.len() >= self.target_batch_count {
            return Ok(self.emit_batch(ctx));
        }

        // 3. Process up to `MAX_BATCHES_PER_LOCK` input batches per call.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(batch) = ctx.input.pop() else { break };
            did_work = true;
            let DecodedRecordBatch { records, .. } = batch;
            for decoded in records {
                let name_hash = decoded.key.name_hash;
                let raw = decoded.into_raw_bytes();
                let read_name = fgumi_raw_bam::read_name(raw.as_ref());
                let same_template = match (self.current_name_hash, self.current_name.as_deref()) {
                    (Some(h), Some(name)) => h == name_hash && name == read_name,
                    _ => false,
                };
                if same_template {
                    self.current_records.push(raw);
                } else {
                    self.flush_current_template()?;
                    self.current_name = Some(read_name.to_vec());
                    self.current_name_hash = Some(name_hash);
                    self.current_records.push(raw);
                }
            }
            if self.accumulator.len() >= self.target_batch_count {
                break;
            }
        }

        if self.accumulator.len() >= self.target_batch_count {
            return Ok(self.emit_batch(ctx));
        }
        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // 4. No input this call. If upstream is drained, close the
        // run-in-progress, emit the final partial batch, and report
        // `Finished` once nothing remains. `held` is empty here (step 1
        // returned `Contention` otherwise); `emit_batch` parks a bounced
        // final push in `held` for step 1 to retry next pass.
        if ctx.input.is_drained() {
            self.flush_current_template()?;
            if !self.accumulator.is_empty() {
                return Ok(self.emit_batch(ctx));
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::core::item::Ordered;

    #[test]
    fn profile_advertises_serial_byordinal() {
        let s = GroupByQueryname::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "GroupByQueryname");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn batched_templates_carries_ordinal() {
        let templates = vec![];
        let wrapped = BamTemplateBatch::new(11, templates);
        assert_eq!(wrapped.ordinal(), 11);
    }
}
