//! `GroupByMi` step for the new typed-step pipeline framework.
//!
//! Mirrors the legacy pipeline's MI-tag grouping (`fgumi simplex`,
//! `fgumi duplex`): consumes `DecodedRecordBatch` items in input order,
//! accumulates records by MI tag value via run-length grouping (assumes
//! input is already sorted by MI, which is the contract upstream of
//! consensus calling), and emits **batches** of completed `MiGroup`s as
//! a single `BatchedMiGroups` per output ordinal.
//!
//! Design follows `GroupByPosition`: `Serial` + `ByItemOrdinal`, batched
//! output to amortize downstream Serial mutex cost, held-slot retry
//! pattern for backpressure.

use std::io;

use fgumi_raw_bam::find_string_tag_in_record;

use crate::mi_group::MiGroup;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::DecodedRecordBatch;
use fgumi_bam_io::MemoryEstimate;

/// Max input batches consumed per `try_run` invocation. Amortizes the
/// Serial mutex acquisition; mirrors `GroupByPosition::MAX_BATCHES_PER_LOCK`.
const MAX_BATCHES_PER_LOCK: usize = 8;

/// Default target batch count. Mirrors legacy's per-batch group count
/// for simplex (50 groups per `MiGroupBatch`).
pub const DEFAULT_TARGET_BATCH_COUNT: usize = 50;

/// A batch of `MiGroup`s carrying its monotonic ordinal.
#[derive(Debug)]
pub struct BatchedMiGroups {
    pub batch_serial: u64,
    pub groups: Vec<MiGroup>,
}

impl BatchedMiGroups {
    #[must_use]
    pub fn new(batch_serial: u64, groups: Vec<MiGroup>) -> Self {
        Self { batch_serial, groups }
    }
}

impl HeapSize for BatchedMiGroups {
    fn heap_size(&self) -> usize {
        self.groups.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.groups.capacity() * std::mem::size_of::<MiGroup>()
    }
}

impl Ordered for BatchedMiGroups {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

/// Type alias for the optional MI-tag transform closure (mirrors
/// `MiGrouper`'s `mi_transform`). When set, the closure is invoked
/// with the raw MI tag bytes and produces a transformed key — used
/// e.g. by `duplex` to strip `/A` and `/B` suffixes so both strands
/// of a paired molecule group together.
type MiTransformFn = Box<dyn Fn(&[u8]) -> String + Send + Sync>;

/// Type alias for the optional record-filter closure (mirrors
/// `MiGrouper`'s `record_filter`). When set, records that fail the
/// filter are skipped without contributing to any MI group.
type RecordFilterFn = Box<dyn Fn(&[u8]) -> bool + Send + Sync>;

/// `Serial + ByItemOrdinal` MI-tag grouper. Records arriving in input
/// order are partitioned into runs of consecutive same-MI records;
/// each run becomes one [`MiGroup`]. Completed groups accumulate
/// until the target batch size is reached, then emit as one
/// [`BatchedMiGroups`].
///
/// State is held behind the framework's per-step mutex (the runtime
/// stores a `Serial` step inside `Arc<Mutex<...>>` and acquires
/// per-`try_run`). Held-slot retry pattern matches `GroupByPosition`:
/// when an emit is rejected by downstream backpressure, the partial
/// state stays buffered until the held slot drains.
pub struct GroupByMi {
    /// MI tag bytes to look up in each record's auxiliary data.
    tag: [u8; 2],
    /// Optional cell barcode tag. When set, the group key becomes
    /// `MI\tCELL`, so reads from different cells with the same MI
    /// tag are placed in separate groups (mirrors `MiGrouper`).
    cell_tag: Option<[u8; 2]>,
    /// Optional record filter (skip records that don't match).
    record_filter: Option<RecordFilterFn>,
    /// Optional MI-tag transformation. When set, the raw MI bytes are
    /// piped through this closure before being used as the group key
    /// (and stored on the emitted `MiGroup`). When unset, the raw MI
    /// bytes are used as a UTF-8 lossy string.
    mi_transform: Option<MiTransformFn>,
    /// Run-length state: the MI value currently being accumulated.
    current_mi: Option<String>,
    /// Run-length state: records accumulated for `current_mi`.
    current_records: Vec<fgumi_raw_bam::RawRecord>,
    /// Self-managed monotonic ordinal — each emitted *batch* gets the
    /// next value. Required for `BranchOrdering::ByItemOrdinal`'s
    /// `ReorderStage` to see contiguous serials.
    next_ordinal: u64,
    /// Accumulator: completed groups waiting to be packaged into a batch.
    accumulator: Vec<MiGroup>,
    /// Held output slot when downstream rejected the most recent push.
    held: HeldSlot<Unpushed<BatchedMiGroups>>,
    target_batch_count: usize,
    output_byte_limit: u64,
    name: &'static str,
}

impl GroupByMi {
    /// Construct a `GroupByMi` step grouping by the given 2-byte tag
    /// (typically `"MI"`).
    #[must_use]
    pub fn new(tag: [u8; 2], output_byte_limit: u64) -> Self {
        Self::with_target_batch_count(tag, output_byte_limit, DEFAULT_TARGET_BATCH_COUNT)
    }

    /// Construct with a custom target batch count.
    #[must_use]
    pub fn with_target_batch_count(
        tag: [u8; 2],
        output_byte_limit: u64,
        target_batch_count: usize,
    ) -> Self {
        Self {
            tag,
            cell_tag: None,
            record_filter: None,
            mi_transform: None,
            current_mi: None,
            current_records: Vec::new(),
            next_ordinal: 0,
            accumulator: Vec::with_capacity(target_batch_count),
            held: HeldSlot::new(),
            target_batch_count: target_batch_count.max(1),
            output_byte_limit,
            name: "GroupByMi",
        }
    }

    /// Set an optional cell-barcode tag for composite grouping.
    #[must_use]
    pub fn with_cell_tag(mut self, cell_tag: Option<[u8; 2]>) -> Self {
        self.cell_tag = cell_tag;
        self
    }

    /// Install an optional record filter. Records that fail are dropped
    /// before they contribute to any MI group. Mirrors
    /// `MiGrouper::with_filter_and_transform`'s filter closure.
    #[must_use]
    pub fn with_record_filter<F>(mut self, filter: F) -> Self
    where
        F: Fn(&[u8]) -> bool + Send + Sync + 'static,
    {
        self.record_filter = Some(Box::new(filter));
        self
    }

    /// Install an optional MI-tag transform closure. The raw MI tag
    /// bytes are piped through this closure to produce the per-record
    /// group key. Used by `duplex` to strip `/A`/`/B` suffixes so both
    /// strands of a paired molecule end up in the same MI group.
    #[must_use]
    pub fn with_mi_transform<F>(mut self, transform: F) -> Self
    where
        F: Fn(&[u8]) -> String + Send + Sync + 'static,
    {
        self.mi_transform = Some(Box::new(transform));
        self
    }

    /// Compute the composite group key for a record's bytes:
    /// `MI` alone, or `MI\tCELL` when `cell_tag` is set. When
    /// `mi_transform` is installed, the raw MI bytes are piped through
    /// it before forming the key.
    /// Returns `None` if the record lacks an MI tag.
    fn group_key(&self, bam: &[u8]) -> Option<String> {
        let mi_value = find_string_tag_in_record(bam, self.tag)?;
        let mut key = if let Some(ref transform) = self.mi_transform {
            transform(mi_value)
        } else {
            String::from_utf8_lossy(mi_value).into_owned()
        };
        if let Some(ct) = self.cell_tag {
            key.push('\t');
            if let Some(cell_value) = find_string_tag_in_record(bam, ct) {
                key.push_str(&String::from_utf8_lossy(cell_value));
            }
        }
        Some(key)
    }

    /// Run the optional filter against a record's bytes; returns `true`
    /// to keep, `false` to discard.
    #[inline]
    fn passes_filter(&self, bam: &[u8]) -> bool {
        match &self.record_filter {
            Some(filter) => filter(bam),
            None => true,
        }
    }

    /// Flush the run-in-progress into the accumulator (if non-empty).
    fn flush_current_group(&mut self) {
        if let Some(mi) = self.current_mi.take() {
            if !self.current_records.is_empty() {
                let records = std::mem::take(&mut self.current_records);
                self.accumulator.push(MiGroup::new(mi, records));
            }
        }
    }

    /// Package the accumulator into a `BatchedMiGroups` and emit. On
    /// rejection, hold; the next `try_run` will retry via the held slot.
    /// Always returns `Progress` because we either emitted or queued
    /// a batch for retry.
    fn emit_batch(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let serial = self.next_ordinal;
        self.next_ordinal += 1;
        let groups =
            std::mem::replace(&mut self.accumulator, Vec::with_capacity(self.target_batch_count));
        let out = BatchedMiGroups::new(serial, groups);
        if let Err(unpushed) = ctx.outputs.push(out) {
            self.held.put(unpushed);
        }
        StepOutcome::Progress
    }
}

impl Step for GroupByMi {
    type Input = DecodedRecordBatch;
    type Outputs = OrderedBytesSingle<BatchedMiGroups>;

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

        // 2. If the accumulator is full enough, emit a batch first;
        // hold off on more input until it lands.
        if self.accumulator.len() >= self.target_batch_count {
            return Ok(self.emit_batch(ctx));
        }

        // 3. Process up to `MAX_BATCHES_PER_LOCK` input batches per
        // call to amortize the Serial mutex acquisition.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(batch) = ctx.input.pop() else { break };
            did_work = true;
            let DecodedRecordBatch { records, .. } = batch;
            for decoded in records {
                let raw = decoded.into_raw_bytes();
                // Optional filter (mirrors `MiGrouper::should_keep`).
                if !self.passes_filter(raw.as_ref()) {
                    continue;
                }
                // Records without MI are skipped (matches `MiGrouper`).
                let Some(mi) = self.group_key(raw.as_ref()) else {
                    continue;
                };
                match &self.current_mi {
                    Some(current) if current == &mi => {
                        self.current_records.push(raw);
                    }
                    Some(_) => {
                        self.flush_current_group();
                        self.current_mi = Some(mi);
                        self.current_records.push(raw);
                    }
                    None => {
                        self.current_mi = Some(mi);
                        self.current_records.push(raw);
                    }
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
        // in-progress MI run, emit the final partial batch, and report
        // `Finished` once nothing remains. `held` is empty here (step 1
        // returned `Contention` otherwise); a bounced `emit_batch` push is
        // parked in `held` for step 1 to retry next pass.
        if ctx.input.is_drained() {
            self.flush_current_group();
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

    #[test]
    fn profile_advertises_serial_byordinal() {
        use crate::sam::SamTag;
        let s = GroupByMi::new(*SamTag::MI, 1024);
        let p = s.profile();
        assert_eq!(p.name, "GroupByMi");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn batched_mi_groups_carries_ordinal() {
        let groups = vec![MiGroup::new("0".to_string(), Vec::new())];
        let wrapped = BatchedMiGroups::new(7, groups);
        assert_eq!(wrapped.ordinal(), 7);
        assert_eq!(wrapped.groups.len(), 1);
    }
}
