//! One range of a range-parallel merge: read its run frames, decode them,
//! drop records outside the range, and k-way merge what is left.
//!
//! The merge is the serial merge restricted to `[lo, hi)`: the same leaves in
//! the same order (runs by `file_id`), the same loser
//! tree, the same tie rule (key, then leaf index). The bounds are points in
//! that order (key, leaf, entry; see [`RangeBound`]), so every record of a
//! source that sorts in `[lo, hi)` is kept in source order, and concatenating
//! every range's output in range order reproduces the serial merge byte for
//! byte.

use std::cell::RefCell;
use std::io;
use std::ops::Range;
use std::sync::Arc;

use crate::keys::RawSortKey;
use crate::loser_tree::LoserTree;
use crate::merge_driver::{slot_key_size, slot_parse_key_in};
use crate::merge_slots::SortMergeSlot;
use crate::range_merge::plan::{RangeBound, RangePlan, RangeSource};
use crate::run_index::RunIndex;
use crate::spill_block::RECORD_LEN_BYTES;
use crate::spill_block_reader::{ZSTD_LEN_PREFIX, ZstdFrameDecoder, checked_zstd_len};

thread_local! {
    /// Each pool thread's zstd frame decoder, reused across the range tasks it
    /// runs rather than built per task. zstd-only: it holds no BGZF state or
    /// scratch buffer, and decodes into each cursor's own exactly-sized buffer.
    static DECODER: RefCell<ZstdFrameDecoder> = RefCell::new(ZstdFrameDecoder::new());
}

/// Counters of one range merge.
#[derive(Clone, Copy, Debug, Default)]
pub struct RangeStats {
    /// Records the range emitted.
    pub records: u64,
    /// Decompressed bytes of the run frames the range decoded.
    pub decoded_bytes: u64,
    /// Decoded bytes (`[key][u32 len][body]`, as a spill frame holds them) of
    /// the records the range emitted, so `decoded_bytes / kept_run_bytes` is
    /// the decode amplification of the run data: 1.0 when no frame overlaps a
    /// range edge.
    pub kept_run_bytes: u64,
    /// Spill file bytes the range read.
    pub read_bytes: u64,
}

impl RangeStats {
    /// Add `other`'s counters to these.
    pub fn add(&mut self, other: &RangeStats) {
        self.records += other.records;
        self.decoded_bytes += other.decoded_bytes;
        self.kept_run_bytes += other.kept_run_bytes;
        self.read_bytes += other.read_bytes;
    }
}

/// A merged range.
#[derive(Debug, Default)]
pub struct RangeOut {
    /// The range's records as `[u32 LE len][body]`, in merge order.
    pub framed: Vec<u8>,
    /// `lens[i]` is the body length of record `i`.
    pub lens: Vec<u32>,
    /// The range's counters.
    pub stats: RangeStats,
}

/// The merge of one key range, free of the key type.
pub trait RangeTask: Send {
    /// Read every run source's frame span for the range from disk. Must be
    /// called once, before [`Self::merge`].
    ///
    /// # Errors
    ///
    /// An I/O error, a spill file shorter than its index says, or a second
    /// call.
    fn read(&mut self) -> io::Result<()>;

    /// Merge up to `max_records` more records. Returns `Some` with the whole
    /// range's output once the range is done, `None` while records remain.
    ///
    /// # Errors
    ///
    /// A corrupt frame, a record crossing a frame end, a call before
    /// [`Self::read`] or after the range is done.
    fn merge(&mut self, max_records: u64) -> io::Result<Option<RangeOut>>;

    /// Hand the task a finished range's output buffers to build its own
    /// output in (their contents are discarded; only their capacity is
    /// reused), so a large output buffer is not freed and page-faulted in
    /// again per range. Call before [`Self::read`]; later calls, and the
    /// default implementation, drop them.
    fn reuse_output_buffers(&mut self, framed: Vec<u8>, lens: Vec<u32>) {
        drop((framed, lens));
    }

    /// Merge up to `max_records` more records into `sink`, one body per
    /// record in merge order, instead of into the range's output buffer.
    /// Returns the range's counters once the range is done, `None` while
    /// records remain. A task merges either through [`Self::merge`] or
    /// through this, never both.
    ///
    /// # Errors
    ///
    /// As [`Self::merge`], a sink error, a task already merging through
    /// [`Self::merge`], or a task that cannot merge into a sink (the
    /// default implementation).
    fn merge_into(
        &mut self,
        sink: &mut dyn RecordSink,
        max_records: u64,
    ) -> io::Result<Option<RangeStats>> {
        let _ = (sink, max_records);
        Err(io::Error::new(io::ErrorKind::Unsupported, "this range task cannot merge into a sink"))
    }
}

/// Where [`RangeTask::merge_into`] puts the merged records.
pub trait RecordSink {
    /// Take the next record's body (no length prefix).
    ///
    /// # Errors
    ///
    /// The sink's own error; the merge stops with it.
    fn push_record(&mut self, body: &[u8]) -> io::Result<()>;
}

/// Clear `buf` and make room for exactly `need` elements without copying
/// its old contents: a buffer already large enough is kept, a smaller one is
/// replaced (growing it in place could copy its stale bytes).
fn reset_with_capacity<T>(buf: &mut Vec<T>, need: usize) {
    buf.clear();
    if buf.capacity() < need {
        *buf = Vec::with_capacity(need);
    }
}

/// A range plan seen through its ranges, free of the key type.
pub trait RangeMerge: Send + Sync {
    /// Number of ranges.
    fn range_count(&self) -> usize;
    /// Bytes range `j` holds while it runs ([`RangePlan::estimated_bytes`]).
    fn estimated_bytes(&self, j: usize) -> u64;
    /// Run bytes decoded over all ranges per run byte
    /// ([`RangePlan::planned_decode_amplification`]).
    fn planned_decode_amplification(&self) -> f64;
    /// Records across all ranges.
    fn total_records(&self) -> u64;
    /// The merge of range `j`.
    ///
    /// # Panics
    ///
    /// If `j >= self.range_count()`.
    fn task(self: Arc<Self>, j: usize) -> Box<dyn RangeTask>;
}

impl<K: RawSortKey + Default + Copy + Send + Sync + 'static> RangeMerge for RangePlan<K> {
    fn range_count(&self) -> usize {
        RangePlan::range_count(self)
    }

    fn estimated_bytes(&self, j: usize) -> u64 {
        RangePlan::estimated_bytes(self, j)
    }

    fn planned_decode_amplification(&self) -> f64 {
        RangePlan::planned_decode_amplification(self)
    }

    fn total_records(&self) -> u64 {
        RangePlan::total_records(self)
    }

    fn task(self: Arc<Self>, j: usize) -> Box<dyn RangeTask> {
        Box::new(PlanRangeTask::new(self, j))
    }
}

fn invalid_data(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Fill `buf` from `slot`'s file at `offset`
/// ([`fgumi_bam_io::pread::read_at_exact`], which retries short and
/// interrupted reads). An EOF before `buf` is full is a truncated run, named
/// by its slot.
fn pread_exact(slot: &SortMergeSlot, buf: &mut [u8], offset: u64) -> io::Result<()> {
    fgumi_bam_io::pread::read_at_exact(slot.source().as_ref(), buf, offset).map_err(|e| {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "spill slot {} truncated: {} bytes at offset {offset}: {e}",
                    slot.slot_id(),
                    buf.len()
                ),
            )
        } else {
            e
        }
    })
}

/// One run source's progress through its frame span.
struct RunCursor<K: RawSortKey> {
    /// Index of the source in the plan.
    source: usize,
    /// Index entries of the span not yet started.
    frames: Range<usize>,
    /// The span's file bytes (`[u32 len][zstd frame]` units), read by `read`.
    comp: Vec<u8>,
    /// Offset in `comp` of the next zstd frame to decode.
    comp_pos: usize,
    /// Offset in `comp` where the entry being walked ends (`comp_pos` once
    /// all its frames are decoded).
    entry_end: usize,
    /// The entry being walked.
    entry: usize,
    /// Decoded bytes of the entry being walked so far.
    entry_raw: u64,
    /// The decoded frame being walked.
    raw: Vec<u8>,
    /// Offset in `raw` of the next record.
    raw_pos: usize,
    /// The current record's body in `raw`.
    body: Range<usize>,
    /// Bytes of the serialized key ahead of each record's length prefix.
    key_size: usize,
    /// The key's per-source merge hint.
    hint: K::MergeHint,
}

/// Where a [`PlanRangeTask`] is.
enum Phase {
    /// `read` not called yet.
    Unread,
    /// Frames read; the merge not started.
    Read,
    /// Merging.
    Merging,
    /// The range was returned, or a call failed.
    Done,
}

/// [`RangeTask`] over range `j` of a [`RangePlan`].
struct PlanRangeTask<K: RawSortKey> {
    plan: Arc<RangePlan<K>>,
    j: usize,
    lo: Option<RangeBound<K>>,
    hi: Option<RangeBound<K>>,
    ctx: K::MergeCtx,
    phase: Phase,
    cursors: Vec<RunCursor<K>>,
    /// Loser-tree leaf -> cursor index (cursors with a first record, in
    /// source order).
    leaves: Vec<usize>,
    tree: Option<LoserTree<K>>,
    framed: Vec<u8>,
    lens: Vec<u32>,
    /// The read spans' raw bytes and records: what the output buffers
    /// reserve when the range merges into them.
    span: (usize, usize),
    /// Whether the merge runs into a sink (`merge_into`), once started.
    into_sink: Option<bool>,
    stats: RangeStats,
}

impl<K: RawSortKey + Default + Copy + Send + Sync + 'static> PlanRangeTask<K> {
    fn new(plan: Arc<RangePlan<K>>, j: usize) -> Self {
        let (lo, hi) = plan.bounds(j);
        let (lo, hi) = (lo.copied(), hi.copied());
        Self {
            plan,
            j,
            lo,
            hi,
            ctx: K::MergeCtx::default(),
            phase: Phase::Unread,
            cursors: Vec::new(),
            leaves: Vec::new(),
            tree: None,
            framed: Vec::new(),
            lens: Vec::new(),
            span: (0, 0),
            into_sink: None,
            stats: RangeStats::default(),
        }
    }

    fn read_spans(&mut self) -> io::Result<()> {
        let plan = Arc::clone(&self.plan);
        let key_size = slot_key_size::<K>().map_err(|e| invalid_data(e.to_string()))?;
        let (mut span_raw, mut span_records) = (0u64, 0u64);
        for (source, RangeSource { slot, index }) in plan.sources().iter().enumerate() {
            let frames = plan.frame_span(source, self.j);
            if frames.is_empty() {
                continue;
            }
            let (start, end) = (index.offset(frames.start), index.end(frames.end - 1));
            let len = usize::try_from(end - start).map_err(|_| {
                invalid_data(format!("spill slot {} span too large", slot.slot_id()))
            })?;
            let mut comp = vec![0u8; len];
            pread_exact(slot, &mut comp, start)?;
            self.stats.read_bytes += len as u64;
            for entry in frames.clone() {
                // A kept record's `[u32 len][body]` is never longer than its
                // spill framing, so the span's raw bytes and records bound the
                // output: reserving them up front means the output never
                // reallocates (the fit charges the span's raw bytes for it).
                span_raw += u64::from(index.raw_len(entry));
                span_records += u64::from(index.records(entry));
            }
            self.cursors.push(RunCursor {
                source,
                frames,
                comp,
                comp_pos: 0,
                entry_end: 0,
                entry: 0,
                entry_raw: 0,
                raw: Vec::new(),
                raw_pos: 0,
                body: 0..0,
                key_size,
                hint: K::MergeHint::default(),
            });
        }
        let too_large = || invalid_data(format!("range {} output too large", self.j));
        self.span = (
            usize::try_from(span_raw).map_err(|_| too_large())?,
            usize::try_from(span_records).map_err(|_| too_large())?,
        );
        Ok(())
    }

    /// Enter or continue the merge, into a sink or into the output buffers
    /// (`into_sink`): on the first call reserve the buffers (when merging into
    /// them) and prime the loser tree. Leaves the phase `Done` until the
    /// caller sets it back, so an error ends the task.
    fn start(&mut self, into_sink: bool) -> io::Result<()> {
        let phase = std::mem::replace(&mut self.phase, Phase::Done);
        match phase {
            Phase::Unread => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "range task merged before it was read",
            )),
            Phase::Done => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "range task merged after it finished",
            )),
            Phase::Read => {
                self.into_sink = Some(into_sink);
                if !into_sink {
                    reset_with_capacity(&mut self.framed, self.span.0);
                    reset_with_capacity(&mut self.lens, self.span.1);
                }
                self.prime()
            }
            Phase::Merging if self.into_sink == Some(into_sink) => Ok(()),
            Phase::Merging => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "range task switched between merging into a sink and into its buffers",
            )),
        }
    }

    /// Advance cursor `c` to its next record in range, returning its key, or
    /// `None` when the cursor has no more records in range.
    fn advance(&mut self, c: usize) -> io::Result<Option<K>> {
        // Borrow the fields apart, so the plan is read in place (no `Arc`
        // clone, an atomic shared by every range worker, per record).
        let Self { plan, cursors, ctx, lo, hi, stats, .. } = self;
        let run = &mut cursors[c];
        let RangeSource { slot, index } = &plan.sources()[run.source];
        next_run_record(run, slot, index, ctx, lo.as_ref(), hi.as_ref(), stats)
    }

    /// Append cursor `c`'s current record to `sink`, or to the output
    /// buffers when there is none.
    fn emit<S: RecordSink + ?Sized>(&mut self, c: usize, sink: Option<&mut S>) -> io::Result<()> {
        let run = &self.cursors[c];
        let body = &run.raw[run.body.clone()];
        // A body came from a `u32` length prefix.
        let len = u32::try_from(body.len()).expect("a record body fits u32");
        self.stats.kept_run_bytes += (run.key_size + RECORD_LEN_BYTES) as u64 + u64::from(len);
        if let Some(sink) = sink {
            sink.push_record(body)?;
        } else {
            self.framed.extend_from_slice(&len.to_le_bytes());
            self.framed.extend_from_slice(body);
            self.lens.push(len);
        }
        self.stats.records += 1;
        Ok(())
    }

    /// Build the loser tree from each cursor's first record in range.
    fn prime(&mut self) -> io::Result<()> {
        let mut keys = Vec::with_capacity(self.cursors.len());
        for c in 0..self.cursors.len() {
            if let Some(key) = self.advance(c)? {
                keys.push(key);
                self.leaves.push(c);
            }
        }
        if !keys.is_empty() {
            self.tree = Some(LoserTree::new(keys));
        }
        Ok(())
    }

    /// Merge up to `max_records` into `sink` (or the output buffers); `true`
    /// when the range is exhausted.
    fn merge_some<S: RecordSink + ?Sized>(
        &mut self,
        max_records: u64,
        mut sink: Option<&mut S>,
    ) -> io::Result<bool> {
        let Some(mut tree) = self.tree.take() else {
            return Ok(true);
        };
        let mut merged = 0u64;
        let done = loop {
            if !tree.winner_is_active() {
                break true;
            }
            if merged == max_records {
                break false;
            }
            let c = self.leaves[tree.winner()];
            self.emit(c, sink.as_deref_mut())?;
            merged += 1;
            match self.advance(c)? {
                Some(key) => tree.replace_winner(key),
                None => tree.remove_winner(),
            }
        };
        self.tree = Some(tree);
        Ok(done)
    }
}

/// The next record of `run` that sorts in `[lo, hi)`, decoding frames as the
/// walk reaches them. Records below `lo` are skipped; the first record at or
/// above `hi` ends the cursor (a run is sorted) and frees its buffers.
fn next_run_record<K: RawSortKey + Default + 'static>(
    run: &mut RunCursor<K>,
    slot: &SortMergeSlot,
    index: &RunIndex,
    ctx: &K::MergeCtx,
    lo: Option<&RangeBound<K>>,
    hi: Option<&RangeBound<K>>,
    stats: &mut RangeStats,
) -> io::Result<Option<K>> {
    let key_size = run.key_size;
    loop {
        if run.raw_pos >= run.raw.len() {
            if run.comp_pos >= run.entry_end {
                let Some(entry) = run.frames.next() else {
                    run.comp = Vec::new();
                    run.raw = Vec::new();
                    run.raw_pos = 0;
                    return Ok(None);
                };
                run.entry = entry;
                run.entry_end = run.comp_pos + comp_len_usize(slot, index, entry)?;
                run.entry_raw = 0;
            }
            decode_frame(run, slot, index)?;
            stats.decoded_bytes += run.raw.len() as u64;
            continue;
        }
        let p = run.raw_pos;
        let header_end = p + key_size + RECORD_LEN_BYTES;
        let not_aligned = || {
            invalid_data(format!(
                "spill frame is not record-aligned (slot {}, decoded offset {p})",
                slot.slot_id()
            ))
        };
        let len_bytes = run.raw.get(p + key_size..header_end).ok_or_else(not_aligned)?;
        let len = u32::from_le_bytes(len_bytes.try_into().expect("a 4-byte slice")) as usize;
        let body_end =
            header_end.checked_add(len).filter(|&e| e <= run.raw.len()).ok_or_else(not_aligned)?;
        let body = &run.raw[header_end..body_end];
        if K::EMBEDDED_IN_RECORD && body.len() < fgumi_raw_bam::MIN_BAM_RECORD_LEN {
            return Err(invalid_data(format!(
                "spill slot {} record of {} bytes is too short to hold its key",
                slot.slot_id(),
                body.len()
            )));
        }
        let key =
            slot_parse_key_in::<K>(ctx, &mut run.hint, &run.raw[p..p + key_size], body, false)
                .map_err(|e| invalid_data(format!("spill slot {}: {e}", slot.slot_id())))?;
        run.raw_pos = body_end;
        if lo.is_some_and(|lo| lo.has_below(&key, run.source, run.entry)) {
            continue;
        }
        if hi.is_some_and(|hi| !hi.has_below(&key, run.source, run.entry)) {
            run.frames.start = run.frames.end;
            run.entry_end = run.comp_pos;
            run.raw_pos = run.raw.len();
            continue;
        }
        run.body = header_end..body_end;
        return Ok(Some(key));
    }
}

/// Entry `entry`'s file span as a `usize`.
fn comp_len_usize(slot: &SortMergeSlot, index: &RunIndex, entry: usize) -> io::Result<usize> {
    usize::try_from(index.comp_len(entry))
        .map_err(|_| invalid_data(format!("spill slot {} entry too large", slot.slot_id())))
}

/// Decode the next zstd frame (`[u32 len][zstd]` at `run.comp_pos`) of the
/// entry being walked into `run.raw`. A sparse entry holds several frames; its
/// frames must tile it exactly, and decode to the raw length the index records.
fn decode_frame<K: RawSortKey>(
    run: &mut RunCursor<K>,
    slot: &SortMergeSlot,
    index: &RunIndex,
) -> io::Result<()> {
    let frame = run.entry;
    let corrupt = |what: &str| {
        invalid_data(format!("corrupt spill frame {frame} of slot {}: {what}", slot.slot_id()))
    };
    let remaining = run.entry_end - run.comp_pos;
    let prefix = run
        .comp
        .get(run.comp_pos..run.comp_pos + ZSTD_LEN_PREFIX.min(remaining))
        .filter(|p| p.len() == ZSTD_LEN_PREFIX)
        .ok_or_else(|| corrupt("no length prefix"))?;
    let frame_len = checked_zstd_len(prefix)?;
    if frame_len > remaining - ZSTD_LEN_PREFIX {
        return Err(corrupt(&format!(
            "length prefix {frame_len} does not match the index's {}",
            remaining - ZSTD_LEN_PREFIX
        )));
    }
    let unit_end = run.comp_pos + ZSTD_LEN_PREFIX + frame_len;
    let zstd = run
        .comp
        .get(run.comp_pos + ZSTD_LEN_PREFIX..unit_end)
        .ok_or_else(|| corrupt("span shorter than its index"))?;
    DECODER.with(|d| d.borrow_mut().decode_exact_into(zstd, &mut run.raw))?;
    run.comp_pos = unit_end;
    run.raw_pos = 0;
    run.entry_raw += run.raw.len() as u64;
    if run.comp_pos == run.entry_end && run.entry_raw != u64::from(index.raw_len(frame)) {
        return Err(corrupt(&format!(
            "decoded {} bytes, the index says {}",
            run.entry_raw,
            index.raw_len(frame)
        )));
    }
    Ok(())
}

impl<K: RawSortKey + Default + Copy + Send + Sync + 'static> RangeTask for PlanRangeTask<K> {
    fn read(&mut self) -> io::Result<()> {
        if !matches!(self.phase, Phase::Unread) {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "range task read twice"));
        }
        self.phase = Phase::Done;
        self.read_spans()?;
        self.phase = Phase::Read;
        Ok(())
    }

    fn merge(&mut self, max_records: u64) -> io::Result<Option<RangeOut>> {
        self.start(false)?;
        if !self.merge_some(max_records, None::<&mut dyn RecordSink>)? {
            self.phase = Phase::Merging;
            return Ok(None);
        }
        self.cursors = Vec::new();
        self.tree = None;
        Ok(Some(RangeOut {
            framed: std::mem::take(&mut self.framed),
            lens: std::mem::take(&mut self.lens),
            stats: self.stats,
        }))
    }

    fn reuse_output_buffers(&mut self, framed: Vec<u8>, lens: Vec<u32>) {
        if matches!(self.phase, Phase::Unread) {
            self.framed = framed;
            self.lens = lens;
        }
    }

    fn merge_into(
        &mut self,
        sink: &mut dyn RecordSink,
        max_records: u64,
    ) -> io::Result<Option<RangeStats>> {
        self.start(true)?;
        if !self.merge_some(max_records, Some(sink))? {
            self.phase = Phase::Merging;
            return Ok(None);
        }
        self.cursors = Vec::new();
        self.tree = None;
        Ok(Some(self.stats))
    }
}

#[cfg(test)]
mod tests;
