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
//! output to amortize downstream `Serial` mutex cost, held-slot retry
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

/// Records retained for one MI key between long-run warnings.
///
/// `current_records` holds the whole run in progress and is not emitted until
/// the key changes or input drains, so its size is set by the largest MI family
/// in the input rather than by any configured limit: `MAX_BATCHES_PER_LOCK`
/// bounds input consumed *per call*, `target_batch_count` counts *completed*
/// groups, and the output queue's byte budget bounds what has been *pushed* —
/// which this has not been. Legacy `MiGrouper` (`crate::mi_group`) retains the
/// same way, so this is the ported behavior rather than a regression, but a
/// pathological family can still exhaust memory with no diagnostic at all.
///
/// This does not bound anything; it makes the growth *visible* before the run
/// dies, which is the cheap half of the fix and costs one comparison per record.
///
/// 100k is ~100x beyond any legitimate family — even deep amplicon data keeps a
/// single UMI at a single position in the thousands — so it should never fire on
/// real input, while still leaving roughly a 10x margin before the retained
/// records amount to serious memory (100k records is tens of MB; 1M is hundreds).
const LONG_RUN_WARN_INTERVAL: usize = 100_000;

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
    /// `current_records.len()` at which the next long-run warning fires, bumped
    /// by [`LONG_RUN_WARN_INTERVAL`] each time and reset at every run boundary.
    ///
    /// A running threshold rather than `len % INTERVAL == 0`: this is checked
    /// once per retained record, and a runtime modulo is an integer division on
    /// that path, where this is a single predictable compare.
    next_long_run_warn_at: usize,
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
    /// Records dropped because they carry no MI tag at all. Expected in small
    /// numbers (e.g. unmapped mates upstream filters missed); reported at
    /// end-of-stream so a silent whole-input drop cannot go unnoticed.
    skipped_no_mi: u64,
    /// Records dropped because MI is present but **not** a `Z` (string) tag.
    /// Almost always a malformed input rather than a legitimate skip: `fgumi
    /// group` and fgbio both write MI as a string (it can carry a `/A`,`/B`
    /// strand suffix), so an `MI:i:` integer silently matches nothing.
    skipped_non_string_mi: u64,
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
            next_long_run_warn_at: LONG_RUN_WARN_INTERVAL,
            next_ordinal: 0,
            accumulator: Vec::with_capacity(target_batch_count),
            held: HeldSlot::new(),
            target_batch_count: target_batch_count.max(1),
            output_byte_limit,
            skipped_no_mi: 0,
            skipped_non_string_mi: 0,
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
                // `None` means no usable MI tag on this record: count why, and skip.
                let Some(matches) =
                    key.matches_record(raw.as_ref(), self.tag, self.cell_tag, self.transform())
                else {
                    self.count_skipped_record(raw.as_ref());
                    return;
                };
                matches
            }
            None => false,
        };
        if same_group {
            self.current_records.push(raw);
            self.warn_if_run_is_long();
        } else {
            // New run (or first record). Build the owned key once at the
            // boundary; `from_record` returns `None` (skip) without an MI tag.
            let Some(key) =
                MiKey::from_record(raw.as_ref(), self.tag, self.cell_tag, self.transform())
            else {
                self.count_skipped_record(raw.as_ref());
                return;
            };
            self.flush_current_group();
            self.current_key = Some(key);
            self.current_records.push(raw);
            // A new run starts the ladder over, so a file with many large-but-
            // legal families warns once per family that crosses the threshold
            // rather than once and then never again.
            self.next_long_run_warn_at = LONG_RUN_WARN_INTERVAL;
        }
    }

    /// Warn every [`LONG_RUN_WARN_INTERVAL`] records retained for one MI key.
    ///
    /// Warning at each multiple rather than once keeps the signal proportional
    /// to the problem: a run that keeps growing keeps saying so, and because the
    /// gap between messages is a fixed number of records, the log cannot outpace
    /// the growth it is reporting.
    #[inline]
    fn warn_if_run_is_long(&mut self) {
        if self.current_records.len() < self.next_long_run_warn_at {
            return;
        }
        let label = self.current_key.as_ref().map_or_else(|| "<none>".to_string(), MiKey::label);
        log::warn!(
            "MI group '{label}' has reached {} records held in memory, which is far beyond a \
             typical family; this step holds a whole MI run before emitting it, so memory will \
             keep growing until the MI tag changes. Check that the input is grouped and sorted \
             by MI and that the MI tag is not constant across records.",
            self.current_records.len(),
        );
        self.next_long_run_warn_at =
            self.current_records.len().saturating_add(LONG_RUN_WARN_INTERVAL);
    }

    /// Classify a record that produced no usable MI key, so end-of-stream can tell
    /// "no MI tag" apart from "MI present but the wrong type".
    ///
    /// Only runs on the skip path, so it costs nothing for records that group
    /// normally.
    #[inline]
    fn count_skipped_record(&mut self, bam: &[u8]) {
        let aux = fgumi_raw_bam::aux_data_slice(bam);
        match fgumi_raw_bam::find_tag_type(aux, self.tag) {
            // Present but not `Z`: the caller almost certainly meant this record
            // to group, so it is tracked separately and warned about loudly.
            Some(kind) if kind != b'Z' => self.skipped_non_string_mi += 1,
            _ => self.skipped_no_mi += 1,
        }
    }

    /// Report records dropped for want of a usable MI tag, once, at end-of-stream.
    ///
    /// A wrong-typed MI is warned about unconditionally: it means the input was
    /// written with e.g. `MI:i:1` instead of `MI:Z:1`, which groups nothing and
    /// otherwise surfaces only as an inexplicably empty output.
    fn report_skipped_records(&self) {
        if self.skipped_non_string_mi > 0 {
            log::warn!(
                "{}: dropped {} record(s) whose {} tag is present but not a string (Z) tag. \
                 fgumi and fgbio write {} as a string because it may carry a /A or /B strand \
                 suffix; an integer tag such as `MI:i:1` matches no group. Re-run `fgumi group` \
                 on this input, or rewrite the tag as `MI:Z:`.",
                self.name,
                self.skipped_non_string_mi,
                String::from_utf8_lossy(&self.tag),
                String::from_utf8_lossy(&self.tag),
            );
        }
        if self.skipped_no_mi > 0 {
            log::warn!(
                "{}: dropped {} record(s) carrying no {} tag.",
                self.name,
                self.skipped_no_mi,
                String::from_utf8_lossy(&self.tag),
            );
        }
    }

    /// Flush the run-in-progress into the accumulator (if non-empty). The owned
    /// display label is materialized here — once per group — from the stored key
    /// bytes, rather than per record.
    fn flush_current_group(&mut self) {
        if let Some(key) = self.current_key.take()
            && !self.current_records.is_empty()
        {
            let records = std::mem::take(&mut self.current_records);
            self.accumulator.push(MiGroup::new(key.label(), records));
        }
    }

    /// Take at most `target_batch_count` groups off the accumulator, leaving
    /// any excess for the next batch.
    ///
    /// The last of the three grouping steps to get this cap; `GroupByPosition`
    /// and `GroupByQueryname` have the identical method for the identical
    /// reason. `process_record` appends per record and the `>=` checks run once
    /// per input batch, so one `DecodedRecordBatch` that closes many MI runs
    /// leaves the accumulator well past the target — handing it over wholesale
    /// would size the emitted batch by the input rather than by
    /// `target_batch_count`, defeating the threshold that exists to bound it.
    /// The `>=` re-entry in `try_run` and the drained path both loop back here
    /// until the accumulator falls below the target, so the excess is emitted
    /// rather than delayed.
    fn drain_one_batch(&mut self) -> Vec<MiGroup> {
        let take = self.accumulator.len().min(self.target_batch_count);
        self.accumulator.drain(..take).collect()
    }

    /// Package at most `target_batch_count` accumulated groups into a
    /// `BatchedMiGroups` and emit. On rejection, hold; the next `try_run` will
    /// retry via the held slot. Always returns `Progress` because we either
    /// emitted or queued a batch for retry. Leftovers above the cap stay in the
    /// accumulator, which is also what makes a rejected push safe.
    fn emit_batch(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let serial = self.next_ordinal;
        self.next_ordinal += 1;
        let groups = self.drain_one_batch();
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
            self.report_skipped_records();
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::RawRecord;
    use rstest::rstest;

    /// Build a minimal unmapped raw BAM record carrying a single Z-type tag.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_with_tag(tag: &str, value: &str) -> RawRecord {
        raw_with_tags(&[(tag, value)])
    }

    /// An emitted batch must never exceed `target_batch_count`.
    ///
    /// `process_record` appends per record while the `>=` checks run once per
    /// input batch, so a single `DecodedRecordBatch` that closes many MI runs
    /// leaves the accumulator well past the target. Taking it wholesale would
    /// size the emitted batch by the input rather than by the configured
    /// threshold. Pins that the excess is preserved rather than dropped, and
    /// that repeated draining converges — which is what makes the `>=` re-entry
    /// in `try_run` and the drained-path loop terminate.
    ///
    /// The third of three: `GroupByPosition` and `GroupByQueryname` carry the
    /// same test for the same guard.
    #[test]
    fn draining_caps_each_batch_and_preserves_the_excess() {
        use crate::sam::SamTag;
        const TARGET: usize = 64;
        let mut step = GroupByMi::with_target_batch_count(*SamTag::MI, 1 << 20, TARGET);

        // 150 completed groups, as one input batch closing many runs would leave.
        step.accumulator = (0..150).map(|i| MiGroup::new(format!("mi{i}"), Vec::new())).collect();

        let first = step.drain_one_batch();
        assert_eq!(first.len(), TARGET, "a batch must never exceed the target");
        assert_eq!(step.accumulator.len(), 86, "the excess stays queued, not dropped");

        let second = step.drain_one_batch();
        assert_eq!(second.len(), TARGET);
        assert_eq!(step.accumulator.len(), 22);

        // The tail comes out as a short final batch and the accumulator empties.
        let third = step.drain_one_batch();
        assert_eq!(third.len(), 22, "the remainder is emitted as a short batch");
        assert!(step.accumulator.is_empty());
        assert!(step.drain_one_batch().is_empty(), "draining an empty accumulator is a no-op");
    }

    /// The long-run warning re-arms at a fixed record interval and starts over
    /// at each run boundary.
    ///
    /// Drives the threshold from a small seeded value rather than pushing
    /// `LONG_RUN_WARN_INTERVAL` real records, so the arithmetic is pinned
    /// without a 100k-record test. What matters is that the next threshold is
    /// derived from the length *reached* (not the previous threshold), so the
    /// ladder cannot fall behind a run that grows by more than one record
    /// between checks, and that a new MI key resets it — otherwise a file of
    /// many large families would warn once and then stay silent.
    #[test]
    fn the_long_run_warning_re_arms_and_resets_at_a_run_boundary() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        assert_eq!(
            step.next_long_run_warn_at, LONG_RUN_WARN_INTERVAL,
            "a fresh step is armed at the interval",
        );

        // Open a run, then seed the threshold low enough to trip on the next
        // record rather than materializing 100k of them.
        step.process_record(raw_with_tag("MI", "fam"));
        step.next_long_run_warn_at = 2;
        step.process_record(raw_with_tag("MI", "fam"));
        assert_eq!(step.current_records.len(), 2, "both records are retained in one run");
        assert_eq!(
            step.next_long_run_warn_at,
            2 + LONG_RUN_WARN_INTERVAL,
            "the next warning is a full interval past the length just reached",
        );

        // A different MI value starts a new run, which re-arms from scratch.
        step.process_record(raw_with_tag("MI", "other"));
        assert_eq!(
            step.next_long_run_warn_at, LONG_RUN_WARN_INTERVAL,
            "a run boundary resets the ladder so each family warns on its own merits",
        );
    }

    /// Build a minimal unmapped raw BAM record carrying a single `i`-type
    /// (signed 32-bit integer) tag, e.g. the malformed `MI:i:1` that groups
    /// nothing because MI is defined as a string tag.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_with_int_tag(tag: &str, value: i32) -> RawRecord {
        let t = tag.as_bytes();
        let mut aux: Vec<u8> = vec![t[0], t[1], b'i'];
        aux.extend_from_slice(&value.to_le_bytes());
        raw_with_aux(&aux)
    }

    /// Build a minimal unmapped raw BAM record carrying the given Z-type tags
    /// in the supplied aux-block order.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_with_tags(tags: &[(&str, &str)]) -> RawRecord {
        let mut aux: Vec<u8> = Vec::new();
        for (tag, value) in tags {
            let t = tag.as_bytes();
            aux.extend_from_slice(&[t[0], t[1], b'Z']);
            aux.extend_from_slice(value.as_bytes());
            aux.push(0);
        }
        raw_with_aux(&aux)
    }

    /// Wrap a pre-encoded aux block in a minimal unmapped raw BAM record.
    #[allow(clippy::cast_possible_truncation)]
    fn raw_with_aux(aux: &[u8]) -> RawRecord {
        let name = b"read";
        let l_read_name: u8 = (name.len() + 1) as u8;
        let seq_len: u32 = 4;
        let seq_bytes = seq_len.div_ceil(2) as usize;

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
        buf[aux_start..aux_start + aux.len()].copy_from_slice(aux);
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

    /// Two MI values that are distinct raw bytes but render to the same
    /// replacement character must stay in separate groups.
    ///
    /// This is the regression test for the `MiKey` change: the previous shape
    /// ran `String::from_utf8_lossy` *before* comparing, so `\x80` and `\x81`
    /// both became U+FFFD and the three records below collapsed into one group
    /// of 3 — silently merging two distinct molecules. `MiKey` compares the raw
    /// tag bytes and only converts lossily to build the display label.
    ///
    /// The labels are asserted to be *identical* on purpose. `MiGroup.mi` is a
    /// diagnostic string, not an identity: nothing rekeys groups by it or writes
    /// it back into a BAM, and treating it as unique is exactly the mistake the
    /// old code made. The record counts are what carry the grouping claim.
    #[test]
    fn distinct_invalid_utf8_mi_values_stay_in_separate_groups() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        let records = vec![
            raw_with_aux(b"MIZ\x80\x00"),
            raw_with_aux(b"MIZ\x80\x00"),
            raw_with_aux(b"MIZ\x81\x00"),
        ];

        let groups = run_grouping(&mut step, records);

        assert_eq!(groups.len(), 2, "distinct raw MI bytes must not merge: {groups:?}");
        assert_eq!(groups[0].1, 2, "the two \\x80 records group together");
        assert_eq!(groups[1].1, 1, "the \\x81 record is its own group");
        assert_eq!(
            groups[0].0, groups[1].0,
            "both labels render to the same replacement character — the label is \
             diagnostic, and the separation above is what proves the keys differ",
        );
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

    /// A dropped record must be *counted*, and a wrong-typed MI counted apart from
    /// an absent one, so end-of-stream can warn instead of emitting a silently
    /// empty output.
    ///
    /// The `MI:i:` case is the one that bites: fgumi and fgbio write MI as a `Z`
    /// string (it may carry a `/A`,`/B` suffix), so an integer MI matches no group
    /// and the whole input vanishes with no diagnostic.
    #[rstest]
    #[case::integer_mi_is_wrong_type(vec![raw_with_int_tag("MI", 1), raw_with_int_tag("MI", 1)], 0, 2)]
    #[case::absent_mi(vec![raw_without_tag(), raw_without_tag()], 2, 0)]
    #[case::mixed(vec![raw_without_tag(), raw_with_int_tag("MI", 3)], 1, 1)]
    #[case::valid_string_mi_counts_nothing(vec![raw_with_tag("MI", "1"), raw_with_tag("MI", "1")], 0, 0)]
    fn skipped_records_are_counted_by_reason(
        #[case] records: Vec<RawRecord>,
        #[case] expected_no_mi: u64,
        #[case] expected_non_string_mi: u64,
    ) {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        run_grouping(&mut step, records);

        assert_eq!(step.skipped_no_mi, expected_no_mi, "records with no MI tag");
        assert_eq!(
            step.skipped_non_string_mi, expected_non_string_mi,
            "records whose MI is present but not a Z tag",
        );
    }

    /// An all-integer-MI input must produce zero groups — the regression that made
    /// four simplex parity tests compare one empty BAM against another.
    #[test]
    fn integer_mi_groups_nothing() {
        use crate::sam::SamTag;
        let mut step = GroupByMi::new(*SamTag::MI, 1 << 20);
        let records = (0..8).map(|_| raw_with_int_tag("MI", 1)).collect();
        assert!(run_grouping(&mut step, records).is_empty());
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
