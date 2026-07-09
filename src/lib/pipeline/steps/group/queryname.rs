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
    /// Run-length state: queryname currently being accumulated, kept in a
    /// reusable buffer so each template boundary refills (`clear` +
    /// `extend_from_slice`) rather than re-allocating a fresh `Vec`. Only
    /// meaningful while a run is active (`current_name_hash.is_some()`).
    current_name: Vec<u8>,
    /// Run-length state: cached `name_hash` for fast pre-check. `Some` iff a
    /// run is currently active; doubles as the run-active flag so the
    /// `current_name` allocation survives across template boundaries.
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
            current_name: Vec::new(),
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

    /// Flush the run-in-progress into the accumulator (if non-empty). The
    /// `current_name` buffer is retained (only the run-active flag is cleared)
    /// so its allocation is reused by the next template boundary.
    fn flush_current_template(&mut self) -> io::Result<()> {
        if !self.current_records.is_empty() {
            let records = std::mem::take(&mut self.current_records);
            let template = Template::from_records(records)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            self.accumulator.push(template);
            self.current_name_hash = None;
        }
        Ok(())
    }

    /// Accumulate one raw record into the current template run, opening a new
    /// run (and flushing the previous one) at a queryname boundary. The cached
    /// `name_hash` pre-check short-circuits most inequality before the byte
    /// compare; the `current_name` buffer is reused across boundaries.
    fn process_record(&mut self, raw: fgumi_raw_bam::RawRecord, name_hash: u64) -> io::Result<()> {
        let read_name = fgumi_raw_bam::read_name(raw.as_ref());
        // A run is active iff `current_name_hash` is `Some`.
        let same_template = match self.current_name_hash {
            Some(h) => h == name_hash && self.current_name == read_name,
            None => false,
        };
        if same_template {
            self.current_records.push(raw);
        } else {
            self.flush_current_template()?;
            // Reuse the buffer: clear + refill rather than allocate.
            self.current_name.clear();
            self.current_name.extend_from_slice(read_name);
            self.current_name_hash = Some(name_hash);
            self.current_records.push(raw);
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
            let records = batch.into_records();
            for decoded in records {
                let name_hash = decoded.key.name_hash;
                let raw = decoded.into_raw_bytes();
                self.process_record(raw, name_hash)?;
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
    use fgumi_bam_io::library::LibraryIndex;
    use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder};

    use fgumi_raw_bam::flags;

    /// Build a minimal mapped raw record with the given read name and flags.
    fn raw_named(name: &[u8], flag: u16) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name).sequence(b"ACGT").qualities(&[30; 4]).flags(flag);
        b.build()
    }

    /// Drive each `(name, flag)` record through [`GroupByQueryname::process_record`],
    /// close the final run, and return `(name, record_count)` per template in order.
    fn run_grouping(
        step: &mut GroupByQueryname,
        records: &[(&[u8], u16)],
    ) -> Vec<(Vec<u8>, usize)> {
        for (name, flag) in records {
            let raw = raw_named(name, *flag);
            let name_hash = LibraryIndex::hash_name(Some(name));
            step.process_record(raw, name_hash).expect("process_record");
        }
        step.flush_current_template().expect("flush");
        step.accumulator.iter().map(|t| (t.name.clone(), t.records().len())).collect()
    }

    const R1: u16 = flags::PAIRED | flags::FIRST_SEGMENT;
    const R2: u16 = flags::PAIRED | flags::LAST_SEGMENT;

    #[test]
    fn consecutive_same_name_forms_one_template() {
        let mut step = GroupByQueryname::new(1 << 20);
        let records: Vec<(&[u8], u16)> = vec![(b"read1", R1), (b"read1", R2)];
        let groups = run_grouping(&mut step, &records);
        assert_eq!(groups, vec![(b"read1".to_vec(), 2)]);
    }

    #[test]
    fn distinct_and_nonadjacent_names_split_runs() {
        let mut step = GroupByQueryname::new(1 << 20);
        // read1 appears again non-adjacently → a fresh template, not merged.
        let records: Vec<(&[u8], u16)> =
            vec![(b"read1", R1), (b"read1", R2), (b"read2", R1), (b"read1", R1)];
        let groups = run_grouping(&mut step, &records);
        assert_eq!(
            groups,
            vec![(b"read1".to_vec(), 2), (b"read2".to_vec(), 1), (b"read1".to_vec(), 1)]
        );
    }

    #[test]
    fn buffer_reuse_across_boundaries_preserves_grouping() {
        // Short-then-long-then-short name transitions exercise the reused
        // `current_name` buffer (clear + refill); grouping must stay exact.
        let mut step = GroupByQueryname::new(1 << 20);
        let records: Vec<(&[u8], u16)> = vec![
            (b"a", R1),
            (b"a", R2),
            (b"longer_name", R1),
            (b"b", R1),
            (b"b", R2),
            (b"another_longer_name", R1),
            (b"c", R1),
        ];
        let groups = run_grouping(&mut step, &records);
        assert_eq!(
            groups,
            vec![
                (b"a".to_vec(), 2),
                (b"longer_name".to_vec(), 1),
                (b"b".to_vec(), 2),
                (b"another_longer_name".to_vec(), 1),
                (b"c".to_vec(), 1),
            ]
        );
    }

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
