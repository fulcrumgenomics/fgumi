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

use crate::mi_group::{MiGroup, MiKey, MiTransform};
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
    /// Run-length state: the grouping key currently being accumulated, stored as
    /// raw comparison bytes so run-boundary detection stays allocation-free on
    /// the common (no-transform) path. The owned display label is materialized
    /// only once per group (at the group boundary) via [`MiKey::label`].
    current_key: Option<MiKey>,
    /// Run-length state: records accumulated for `current_key`.
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
            current_key: None,
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

    /// Borrow the installed MI-transform closure (if any) as the
    /// `Option<&dyn Fn>` shape [`MiKey`] expects.
    #[inline]
    fn transform(&self) -> MiTransform<'_> {
        self.mi_transform.as_ref().map(|t| t.as_ref() as &dyn Fn(&[u8]) -> String)
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

    /// Accumulate one raw record into the current MI run, opening a new run (and
    /// flushing the previous one) at a group boundary. Run-boundary detection is
    /// an allocation-free byte comparison against the stored key on the common
    /// (no-transform) path; the owned key is built only when a new run begins.
    /// Records that fail the optional filter, or that carry no MI tag, are
    /// dropped (matching `MiGrouper`).
    fn process_record(&mut self, raw: fgumi_raw_bam::RawRecord) {
        // Optional filter (mirrors `MiGrouper::should_keep`).
        if !self.passes_filter(raw.as_ref()) {
            return;
        }
        let same_group = match &self.current_key {
            Some(key) => {
                match key.matches_record(raw.as_ref(), self.tag, self.cell_tag, self.transform()) {
                    Some(matches) => matches,
                    // No MI tag on this record: skip it.
                    None => return,
                }
            }
            None => false,
        };
        if same_group {
            self.current_records.push(raw);
        } else {
            // New run (or first record). Build the owned key once at the
            // boundary; `from_record` returns `None` (skip) without an MI tag.
            let Some(key) =
                MiKey::from_record(raw.as_ref(), self.tag, self.cell_tag, self.transform())
            else {
                return;
            };
            self.flush_current_group();
            self.current_key = Some(key);
            self.current_records.push(raw);
        }
    }

    /// Flush the run-in-progress into the accumulator (if non-empty). The owned
    /// display label is materialized here — once per group — from the stored key
    /// bytes, rather than per record.
    fn flush_current_group(&mut self) {
        if let Some(key) = self.current_key.take() {
            if !self.current_records.is_empty() {
                let records = std::mem::take(&mut self.current_records);
                self.accumulator.push(MiGroup::new(key.label(), records));
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
            let records = batch.into_records();
            for decoded in records {
                self.process_record(decoded.into_raw_bytes());
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
    use fgumi_raw_bam::RawRecord;

    /// Build a minimal unmapped raw BAM record carrying a single Z-type tag.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_with_tag(tag: &str, value: &str) -> RawRecord {
        raw_with_tags(&[(tag, value)])
    }

    /// Build a minimal unmapped raw BAM record carrying the given Z-type tags
    /// in the supplied aux-block order.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_with_tags(tags: &[(&str, &str)]) -> RawRecord {
        let name = b"read";
        let l_read_name: u8 = (name.len() + 1) as u8;
        let seq_len: u32 = 4;
        let seq_bytes = seq_len.div_ceil(2) as usize;

        let mut aux: Vec<u8> = Vec::new();
        for (tag, value) in tags {
            let t = tag.as_bytes();
            aux.extend_from_slice(&[t[0], t[1], b'Z']);
            aux.extend_from_slice(value.as_bytes());
            aux.push(0);
        }

        let total = 32 + l_read_name as usize + seq_bytes + seq_len as usize + aux.len();
        let mut buf = vec![0u8; total];
        buf[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        buf[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        buf[8] = l_read_name;
        buf[12..14].copy_from_slice(&0u16.to_le_bytes());
        buf[16..20].copy_from_slice(&seq_len.to_le_bytes());
        buf[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        buf[24..28].copy_from_slice(&(-1i32).to_le_bytes());
        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;
        let aux_start = 32 + l_read_name as usize + seq_bytes + seq_len as usize;
        buf[aux_start..aux_start + aux.len()].copy_from_slice(&aux);
        RawRecord::from(buf)
    }

    /// Build a minimal unmapped raw BAM record carrying no aux tags.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_without_tag() -> RawRecord {
        raw_with_tags(&[])
    }

    /// Drive a sequence of records through [`GroupByMi::process_record`], close
    /// the final run, and return `(label, record_count)` per completed group.
    fn run_grouping(step: &mut GroupByMi, records: Vec<RawRecord>) -> Vec<(String, usize)> {
        for raw in records {
            step.process_record(raw);
        }
        step.flush_current_group();
        step.accumulator.iter().map(|g| (g.mi.clone(), g.records.len())).collect()
    }

    #[test]
    fn single_mi_run_forms_one_group() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        let records: Vec<RawRecord> = (0..1000).map(|_| raw_with_tag("MI", "7")).collect();
        let groups = run_grouping(&mut step, records);
        assert_eq!(groups, vec![("7".to_string(), 1000)]);
    }

    #[test]
    fn distinct_mi_values_split_runs_and_label_once() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        let records = vec![
            raw_with_tag("MI", "1"),
            raw_with_tag("MI", "1"),
            raw_with_tag("MI", "2"),
            raw_with_tag("MI", "1"), // non-adjacent: a fresh run, not merged
        ];
        let groups = run_grouping(&mut step, records);
        assert_eq!(groups, vec![("1".to_string(), 2), ("2".to_string(), 1), ("1".to_string(), 1)]);
    }

    #[test]
    fn records_without_mi_are_skipped() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        let records = vec![
            raw_without_tag(), // skipped, no current run yet
            raw_with_tag("MI", "5"),
            raw_without_tag(), // skipped, current run preserved
            raw_with_tag("MI", "5"),
        ];
        let groups = run_grouping(&mut step, records);
        assert_eq!(groups, vec![("5".to_string(), 2)]);
    }

    #[test]
    fn cell_tag_composite_key_splits_and_labels() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20).with_cell_tag(Some(*SamTag::CB));
        // Same MI, different CB → distinct groups; aux order varies to exercise
        // the single-pass dual-tag lookup regardless of tag ordering.
        let records = vec![
            raw_with_tags(&[("MI", "1"), ("CB", "ACGT")]),
            raw_with_tags(&[("CB", "ACGT"), ("MI", "1")]),
            raw_with_tags(&[("MI", "1"), ("CB", "TGCA")]),
            raw_with_tags(&[("MI", "1"), ("CB", "TGCA")]),
        ];
        let groups = run_grouping(&mut step, records);
        assert_eq!(groups, vec![("1\tACGT".to_string(), 2), ("1\tTGCA".to_string(), 2)]);
    }

    #[test]
    fn missing_cell_tag_uses_empty_suffix() {
        use crate::sam::SamTag;
        // With a cell tag configured, a record that lacks `CB` keys on an empty
        // cell suffix (`MI\t`), distinct from a record that carries `CB`. Guards
        // against a future aux-lookup change silently merging or splitting MI
        // groups when the configured cell tag is absent.
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20).with_cell_tag(Some(*SamTag::CB));
        let records = vec![
            raw_with_tag("MI", "1"),                       // no CB → "1\t"
            raw_with_tag("MI", "1"),                       // no CB → "1\t"
            raw_with_tags(&[("MI", "1"), ("CB", "ACGT")]), // → "1\tACGT"
        ];
        let groups = run_grouping(&mut step, records);
        assert_eq!(groups, vec![("1\t".to_string(), 2), ("1\tACGT".to_string(), 1)]);
    }

    #[test]
    fn mi_transform_groups_strands_together() {
        use crate::sam::SamTag;
        // Strip a trailing "/A" or "/B" so both strands of a molecule group.
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20).with_mi_transform(|mi: &[u8]| {
            let s = String::from_utf8_lossy(mi);
            s.split('/').next().unwrap_or("").to_string()
        });
        let records =
            vec![raw_with_tag("MI", "1/A"), raw_with_tag("MI", "1/B"), raw_with_tag("MI", "2/A")];
        let groups = run_grouping(&mut step, records);
        // Label is the transformed key, materialized once per group.
        assert_eq!(groups, vec![("1".to_string(), 2), ("2".to_string(), 1)]);
    }

    #[test]
    fn record_filter_drops_records() {
        use crate::sam::SamTag;
        // Filter out any record whose MI tag value is "skip".
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20).with_record_filter(|bam: &[u8]| {
            fgumi_raw_bam::find_string_tag_in_record(bam, *SamTag::MI) != Some(b"skip")
        });
        let records = vec![
            raw_with_tag("MI", "1"),
            raw_with_tag("MI", "skip"), // dropped, run "1" preserved
            raw_with_tag("MI", "1"),
        ];
        let groups = run_grouping(&mut step, records);
        assert_eq!(groups, vec![("1".to_string(), 2)]);
    }

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
