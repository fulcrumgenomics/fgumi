//! Deferred sort-key extraction — the ingest thread's largest serial cost,
//! moved onto the worker pool.
//!
//! Phase 1's ingest thread is a serial consumer, and the floor line names it as
//! the phase's binding limit: on a 16-thread whole-genome sort it is 137.5s of a
//! 145.4s read span, against a worker-capacity floor of 22.2s. Of that 137.5s,
//! **93.6s is `extract_template_key_inline`** — one aux-tag scan, one name hash
//! and one unclipped-5' computation per record, at 120 ns each over 780M
//! records. It is pure per-record CPU with no ordering requirement, so the only
//! reason it sat on the serial thread was that it ran where the record bytes
//! happened to be in hand.
//!
//! This module moves it. The ingest thread pushes a record's bytes into the
//! arena and records its extent, leaving the ref's key unset; batches of extents
//! are handed to the pool as [`KeyExtractionJob`]s, and the filled keys come
//! back over a channel to be spliced into the ref array. The pool has the
//! capacity to absorb it — 22.2s of worker busy across a 145s span — so the work
//! disappears rather than relocating.
//!
//! # Why batches are cut at arena segment boundaries
//!
//! A worker reads record bytes straight out of the arena while the ingest thread
//! is still appending to it. That is only sound against a segment the writer has
//! finished with, which is why [`SegmentedBuf`](crate::segmented_buf::SegmentedBuf)
//! splits its sealed segments from its live one: a batch names exactly one
//! sealed segment, holds it alive by `Arc`, and the writer never touches it
//! again. Records still in the live segment keep their keys unextracted until it
//! seals, or until the chunk barrier, whichever comes first.

use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender};

use fgumi_raw_bam::{RawRecordView, SamTag};

use crate::external::{
    DroppedLaneViolation, LibraryLookup, TemplateKeyVariant, extract_template_key_inline,
    verify_dropped_lanes,
};
use crate::inline::{TEMPLATE_HEADER_SIZE, TemplateKey, TemplateLaneKey, TemplateRecordBuffer};
use crate::worker_pool::SortWorkerPool;

/// Everything a deferred batch needs to reproduce the serial extraction exactly.
///
/// Built once per sort and shared by `Arc`, so a batch carries one pointer
/// rather than a copy of the library table and hasher.
pub(crate) struct KeyContext {
    /// Read-group to library-ordinal table, from the BAM header.
    pub(crate) lib_lookup: LibraryLookup,
    /// The cell-barcode tag to fold into the key, when the sort was asked for one.
    pub(crate) cell_tag: Option<SamTag>,
    /// Fixed-seed hasher for cell-barcode values. **Never reseed this** — it
    /// feeds the sort key, so a random seed would break byte-identity.
    pub(crate) cb_hasher: ahash::RandomState,
    /// The first record's full key: the baseline every later record's dropped
    /// lanes are verified against.
    pub(crate) first_key: TemplateKey,
    /// Which lanes the chosen narrowed key retains.
    pub(crate) variant: TemplateKeyVariant,
}

/// A dropped-lane violation, tagged with the record that carried it.
///
/// Carries the record's index within the whole chunk (not within the batch) so
/// the ingest thread can report the *first* offending record deterministically,
/// however the batches happened to be scheduled.
pub(crate) struct KeyViolation {
    /// Index of the offending record in the chunk's ref array.
    pub(crate) ref_index: usize,
    /// Which lane disagreed with the first record.
    pub(crate) violation: DroppedLaneViolation,
    /// The offending record's read name, for the error message.
    pub(crate) name: String,
}

/// One completed batch of keys, addressed by where they belong.
pub(crate) struct KeyBatchResult<K> {
    /// Index of this batch's first record in the chunk's ref array.
    pub(crate) first_ref: usize,
    /// Keys for `first_ref .. first_ref + keys.len()`, in ref order.
    pub(crate) keys: Box<[K]>,
    /// The lowest-indexed violation this batch found, if any.
    pub(crate) violation: Option<KeyViolation>,
}

/// A unit of deferred key extraction that any pool worker can run.
///
/// Type-erased on purpose: the worker pool is not generic over the sort key, and
/// making it so would thread `K` through every step, queue and stats array for
/// the benefit of one step. The batch carries its own sender instead, so the
/// pool only ever sees "run this".
pub(crate) trait KeyExtractionJob: Send {
    /// Extract every key in the batch and publish the result.
    ///
    /// Consumes the job: a batch runs exactly once.
    fn run(self: Box<Self>);
}

/// Extract one record's full key and check the lanes the chosen variant drops.
///
/// The single extraction implementation, shared by the batched path and the
/// immediate one, so the two cannot drift into producing different keys.
fn extract_one(ctx: &KeyContext, bam: &[u8]) -> (TemplateKey, Option<DroppedLaneViolation>) {
    let full = extract_template_key_inline(bam, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher);
    let violation = verify_dropped_lanes(&ctx.first_key, &full, ctx.variant);
    (full, violation)
}

/// A batch of template-coordinate records awaiting key extraction.
pub(crate) struct TemplateKeyBatch<K: TemplateLaneKey> {
    /// Shared extraction context.
    pub(crate) ctx: Arc<KeyContext>,
    /// The sealed arena segment holding every record body in this batch.
    pub(crate) segment: Arc<Vec<u8>>,
    /// Byte offset of `segment`'s first byte in the arena's global address space.
    pub(crate) segment_base: u64,
    /// Index of this batch's first record in the chunk's ref array.
    pub(crate) first_ref: usize,
    /// Per record, in ref order: the global offset of its BAM bytes (the inline
    /// header already skipped) and their length.
    pub(crate) extents: Box<[(u64, u32)]>,
    /// Where the filled keys are published.
    pub(crate) results: Sender<KeyBatchResult<K>>,
}

impl<K: TemplateLaneKey> KeyExtractionJob for TemplateKeyBatch<K>
where
    K: Send,
{
    fn run(self: Box<Self>) {
        let me = *self;
        let ctx = &me.ctx;
        let mut keys: Vec<K> = Vec::with_capacity(me.extents.len());
        let mut violation: Option<KeyViolation> = None;

        for (i, &(offset, len)) in me.extents.iter().enumerate() {
            let local = usize::try_from(offset - me.segment_base)
                .expect("record offset precedes its own segment");
            let bam = &me.segment[local..local + len as usize];

            let (full, lane) = extract_one(ctx, bam);

            // First violation wins: batches are scheduled in whatever order the
            // pool picks them up, so only the lowest index is reproducible.
            if violation.is_none()
                && let Some(v) = lane
            {
                let name =
                    String::from_utf8_lossy(RawRecordView::new(bam).read_name()).into_owned();
                violation = Some(KeyViolation { ref_index: me.first_ref + i, violation: v, name });
            }

            keys.push(K::from_full(&full));
        }

        // A closed channel means the ingest thread has already failed and gone
        // away; there is nobody left to report to and nothing to clean up.
        drop(me.results.send(KeyBatchResult {
            first_ref: me.first_ref,
            keys: keys.into_boxed_slice(),
            violation,
        }));
    }
}

// ============================================================================
// Ingest-side driver
// ============================================================================

/// What deferral actually achieved over a whole sort.
///
/// Reported because a deferral that silently stops engaging is invisible any
/// other way: the keys are correct whether a worker extracted them or the
/// barrier did, so no output check, equivalence test or record count can tell
/// the difference. Only this ratio can.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct KeyOverlapCensus {
    /// Records whose keys were extracted while the ingest thread kept reading.
    pub(crate) overlapped_records: u64,
    /// Records whose keys the ingest thread had to wait at a barrier for.
    pub(crate) barrier_records: u64,
    /// Exact seconds spent cutting and handing out batches.
    pub(crate) dispatch_secs: f64,
    /// Exact seconds spent waiting at barriers for keys.
    pub(crate) barrier_secs: f64,
}

impl KeyOverlapCensus {
    /// Share of records keyed off the serial thread's critical path, as a
    /// percentage. `None` when nothing was deferred at all.
    pub(crate) fn overlap_percent(self) -> Option<f64> {
        let total = self.total_records();
        if total == 0 {
            return None;
        }
        #[allow(clippy::cast_precision_loss)]
        Some(100.0 * self.overlapped_records as f64 / total as f64)
    }

    /// Total records accounted for.
    pub(crate) fn total_records(self) -> u64 {
        self.overlapped_records + self.barrier_records
    }
}

/// Records per deferred batch.
///
/// Big enough that the per-batch cost (one `Box`, one queue push, one channel
/// send, 48 KiB of extents) is noise against ~4k extractions — each batch is
/// still ~0.5 ms of work — and small enough that a batch's key array (4k x up to
/// 56 B = 224 KiB) stays inside L2, so splicing it back into the ref array hits
/// warm cache rather than a second trip to memory. Batches are also cut at arena
/// segment boundaries, so this is an upper bound rather than an exact size.
const KEY_BATCH_RECORDS: usize = 4 * 1024;

/// Drives deferred key extraction for one sort: cuts batches, hands them to the
/// pool, and splices the results back into the buffer's ref array.
///
/// One instance spans the whole sort and is [`reset`](Self::reset) per chunk, so
/// the channel and its context are allocated once rather than per spill.
pub(crate) struct DeferredKeys<K: TemplateLaneKey> {
    /// Shared extraction context, cloned into every batch by `Arc`.
    ctx: Arc<KeyContext>,
    /// Kept alive for the whole sort so the receiver never disconnects while
    /// batches are outstanding.
    tx: Sender<KeyBatchResult<K>>,
    /// Completed batches come back here.
    rx: Receiver<KeyBatchResult<K>>,
    /// Refs `0..dispatched` have been handed to some batch.
    dispatched: usize,
    /// Batches handed out whose result has not been absorbed yet.
    outstanding: usize,
    /// Refs whose key has actually been written back.
    filled: usize,
    /// Sealed arena segment count as of the last dispatch check.
    sealed_seen: usize,
    /// Lowest-indexed dropped-lane violation seen across all batches.
    violation: Option<KeyViolation>,
    /// Records dispatched *before* their chunk's barrier, summed over the whole
    /// sort. Not cleared by [`reset`](Self::reset).
    overlapped_records: u64,
    /// Records still unbatched when their chunk's barrier ran, so their keys
    /// were extracted with the ingest thread waiting. Not cleared by `reset`.
    barrier_records: u64,
    /// Seconds the ingest thread spent cutting and handing out batches.
    ///
    /// Timed exactly rather than sampled: dispatch fires once per sealed arena
    /// segment — a few thousand times across a whole sort, against hundreds of
    /// millions of records — so a 1-in-N per-record sample would be estimating a
    /// rare event with a method built for a uniform one. A few thousand clock
    /// pairs is unmeasurable overhead.
    dispatch_secs: f64,
    /// Seconds the ingest thread spent at a barrier waiting for keys.
    ///
    /// The number to watch: it is the part of extraction the pool did not manage
    /// to hide, and so the part that is still serial. Exact, per chunk.
    barrier_secs: f64,
    /// Whether to defer at all.
    ///
    /// False when the pool has fewer than two workers, because then there is no
    /// second thread for a batch to run on — deferral could only ever move the
    /// same extraction onto the same thread, later, having paid for an extent
    /// list, a key array and a splice on the way. That is not a judgement call
    /// about which is faster; with one worker the batched path is strictly more
    /// work. The single-threaded fast path (`--threads` absent) is exactly this
    /// case, so it keeps its original inline extraction.
    defer: bool,
}

impl<K: TemplateLaneKey> DeferredKeys<K> {
    /// Build a driver over `ctx`, deferring when Phase 1 has the workers to make
    /// deferral worth its bookkeeping.
    ///
    /// `phase1_threads` is the count of workers *active during Phase 1*, not the
    /// pool's size. The pool is sized to the wider of the two phases and capped
    /// per phase, so `--threads 1 --merge-threads 16` builds 16 workers of which
    /// exactly one is awake while ingest runs — gating on the pool size there
    /// would defer batches to a pool that is not allowed to run them.
    pub(crate) fn new(ctx: Arc<KeyContext>, phase1_threads: usize) -> Self {
        let defer = phase1_threads >= 2;
        let (tx, rx) = std::sync::mpsc::channel();
        Self {
            ctx,
            tx,
            rx,
            dispatched: 0,
            outstanding: 0,
            filled: 0,
            sealed_seen: 0,
            violation: None,
            overlapped_records: 0,
            barrier_records: 0,
            dispatch_secs: 0.0,
            barrier_secs: 0.0,
            defer,
        }
    }

    /// Get one record into the buffer with a sort key, now or later.
    ///
    /// The ingest loop's single entry point, so the deferred and immediate modes
    /// cannot diverge at the call site. Immediate mode extracts and pushes in one
    /// step exactly as the pre-deferral loop did; deferred mode pushes a
    /// placeholder and lets the pool fill it in.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be pushed into the arena, or (in
    /// immediate mode) if it carries a value in a lane the chosen key drops.
    pub(crate) fn push(
        &mut self,
        buffer: &mut TemplateRecordBuffer<K>,
        pool: &SortWorkerPool,
        bam: &[u8],
    ) -> anyhow::Result<()> {
        if !self.defer {
            let (full, lane) = extract_one(&self.ctx, bam);
            if let Some(v) = lane {
                let name =
                    String::from_utf8_lossy(RawRecordView::new(bam).read_name()).into_owned();
                return Err(crate::external::dropped_lane_error(&name, v));
            }
            buffer.push(bam, K::from_full(&full))?;
            // Counted as overlapped-with-nothing so the census still sums to the
            // record count and a single-worker run reads as 0% overlapped rather
            // than as no deferral information at all.
            self.barrier_records += 1;
            self.filled += 1;
            self.dispatched += 1;
            return Ok(());
        }
        buffer.push_deferred(bam)?;
        self.after_push(buffer, pool);
        Ok(())
    }

    /// Records whose keys were extracted while the ingest thread kept reading,
    /// and records whose keys it had to wait at the barrier for.
    ///
    /// This is how you tell the difference between deferral *working* and
    /// deferral silently not engaging: the second is indistinguishable from the
    /// first by output, by test, and by every equivalence check — the keys are
    /// right either way, they were just paid for serially. A run where
    /// `overlapped` is near zero has the cost of the old code plus the
    /// bookkeeping of the new.
    pub(crate) fn overlap_census(&self) -> KeyOverlapCensus {
        KeyOverlapCensus {
            overlapped_records: self.overlapped_records,
            barrier_records: self.barrier_records,
            dispatch_secs: self.dispatch_secs,
            barrier_secs: self.barrier_secs,
        }
    }

    /// Bytes currently held by in-flight batches, to be counted against the
    /// sort's memory limit.
    ///
    /// A record whose key is still outstanding is charged twice — once for the
    /// placeholder key already in its ref, once for the copy the batch is
    /// filling — plus its extent. Leaving this out would let the buffer overrun
    /// the limit by however much extraction happens to be lagging.
    pub(crate) fn in_flight_bytes(&self) -> usize {
        let per_record = std::mem::size_of::<K>() + std::mem::size_of::<(u64, u32)>();
        self.dispatched.saturating_sub(self.filled) * per_record
    }

    /// Hand out batches for every record that has become shareable, if the
    /// arena has sealed a segment since the last call.
    ///
    /// Cheap enough to call after every push: the common case is one integer
    /// comparison.
    pub(crate) fn after_push(
        &mut self,
        buffer: &mut TemplateRecordBuffer<K>,
        pool: &SortWorkerPool,
    ) {
        let sealed = buffer.sealed_segments();
        if !self.defer || sealed == self.sealed_seen {
            return;
        }
        self.sealed_seen = sealed;
        let started = std::time::Instant::now();
        // The record that triggered the seal landed in the new live segment, so
        // it is not shareable yet; everything before it is.
        let shareable = buffer.refs().len().saturating_sub(1);
        self.dispatch_upto(buffer, pool, shareable);
        self.absorb_ready(buffer);
        self.dispatch_secs += started.elapsed().as_secs_f64();
    }

    /// Cut and dispatch batches covering refs `self.dispatched..end`.
    ///
    /// Batches never span an arena segment: a batch names exactly one sealed
    /// segment and holds it by `Arc`, which is what makes reading it concurrent
    /// with the ingest thread's appends sound.
    fn dispatch_upto(
        &mut self,
        buffer: &TemplateRecordBuffer<K>,
        pool: &SortWorkerPool,
        end: usize,
    ) {
        let segment_size = buffer.segment_size() as u64;
        while self.dispatched < end {
            let refs = buffer.refs();
            let start = self.dispatched;
            let segment_index = refs[start].offset / segment_size;
            let limit = end.min(start + KEY_BATCH_RECORDS);
            let mut stop = start;
            while stop < limit && refs[stop].offset / segment_size == segment_index {
                stop += 1;
            }

            let segment_index = usize::try_from(segment_index).expect("segment index fits usize");
            let Some(segment) = buffer.sealed_segment(segment_index) else {
                // The segment is still live, so nothing from here on is
                // shareable yet. Leave it for the next seal or the barrier.
                return;
            };

            let extents: Box<[(u64, u32)]> = refs[start..stop]
                .iter()
                .map(|r| (r.offset + TEMPLATE_HEADER_SIZE as u64, r.len))
                .collect();
            let batch = TemplateKeyBatch::<K> {
                ctx: Arc::clone(&self.ctx),
                segment,
                segment_base: segment_index as u64 * segment_size,
                first_ref: start,
                extents,
                results: self.tx.clone(),
            };

            self.dispatched = stop;
            self.outstanding += 1;
            // A full queue means the pool is saturated, not that the batch can
            // be dropped: run it here rather than lose its keys.
            if let Err(batch) = pool.submit_key_job(Box::new(batch)) {
                batch.run();
            }
        }
    }

    /// Absorb every batch that has already finished, without waiting.
    fn absorb_ready(&mut self, buffer: &mut TemplateRecordBuffer<K>) {
        while let Ok(result) = self.rx.try_recv() {
            self.absorb(buffer, result);
        }
    }

    /// Write one batch's keys into the buffer and fold in its violation.
    fn absorb(&mut self, buffer: &mut TemplateRecordBuffer<K>, result: KeyBatchResult<K>) {
        buffer.fill_keys(result.first_ref, &result.keys);
        self.filled += result.keys.len();
        self.outstanding -= 1;
        if let Some(v) = result.violation
            && self.violation.as_ref().is_none_or(|seen| v.ref_index < seen.ref_index)
        {
            self.violation = Some(v);
        }
    }

    /// Wait until every record in the buffer has its real key, then report the
    /// first dropped-lane violation if there was one.
    ///
    /// Must be called before the buffer is sorted or spilled. Runs queued
    /// batches on the calling thread while it waits, so it cannot deadlock
    /// against a pool whose workers have already parked for the phase.
    ///
    /// # Errors
    ///
    /// Returns an error if a worker panicked while holding a batch. That batch's
    /// sender is dropped by the unwind without ever publishing, so its keys can
    /// never arrive — and because this driver holds its own sender, the channel
    /// never disconnects and a plain blocking wait would hang forever. The wait
    /// is therefore bounded and re-checks the pool's panic flag.
    ///
    /// # Panics
    ///
    /// Panics if the buffer still holds a record whose key was never filled —
    /// that would otherwise sort those records under a constant placeholder key
    /// and produce a wrong order that still looks like a valid BAM.
    pub(crate) fn finish(
        &mut self,
        buffer: &mut TemplateRecordBuffer<K>,
        pool: &SortWorkerPool,
    ) -> anyhow::Result<Option<KeyViolation>> {
        let started = std::time::Instant::now();
        let all = buffer.refs().len();
        if self.defer {
            // Everything left is in the live segment; seal it so it can be
            // shared. Skipped in immediate mode, where it would only churn a
            // fresh full-size segment allocation per chunk for no reader.
            buffer.seal_arena_segment();
            self.sealed_seen = buffer.sealed_segments();
            let overlapped = self.dispatched.min(all);
            self.overlapped_records += overlapped as u64;
            self.barrier_records += (all - overlapped) as u64;
            self.dispatch_upto(buffer, pool, all);
        }

        while self.outstanding > 0 {
            match self.rx.try_recv() {
                Ok(result) => {
                    self.absorb(buffer, result);
                    continue;
                }
                Err(std::sync::mpsc::TryRecvError::Empty) => {}
                // Unreachable: `self.tx` is alive for as long as `self` is.
                Err(std::sync::mpsc::TryRecvError::Disconnected) => break,
            }
            if pool.run_one_key_job() {
                continue;
            }
            // Nothing ready and nothing queued: every outstanding batch is
            // running on a worker right now, so waiting is the cheap move. The
            // timeout exists only so a panicked worker surfaces as an error
            // instead of an indefinite hang.
            match self.rx.recv_timeout(std::time::Duration::from_millis(50)) {
                Ok(result) => self.absorb(buffer, result),
                Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
                    anyhow::ensure!(
                        !pool.worker_panicked(),
                        "a sort worker panicked while extracting sort keys; {} key \
                         batches can never complete",
                        self.outstanding,
                    );
                }
                Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => break,
            }
        }

        self.barrier_secs += started.elapsed().as_secs_f64();
        assert_eq!(
            self.filled,
            all,
            "{} of {all} records were never given a sort key; sorting now would collate them \
             under the placeholder",
            all - self.filled,
        );
        Ok(self.violation.take())
    }

    /// Forget the previous chunk's accounting, after the buffer was cleared.
    pub(crate) fn reset(&mut self) {
        self.dispatched = 0;
        self.outstanding = 0;
        self.filled = 0;
        self.sealed_seen = 0;
        self.violation = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external::cb_hasher;
    use crate::inline::TemplateKey40;
    use noodles::sam::Header;

    /// Minimal mapped BAM record: paired, forward, mate at the same locus, with
    /// `aux` appended verbatim.
    #[allow(clippy::cast_possible_truncation)]
    fn mapped_bam(tid: i32, pos: i32, name: &[u8], aux: &[u8]) -> Vec<u8> {
        let mut bam = vec![0u8; 32];
        bam[0..4].copy_from_slice(&tid.to_le_bytes());
        bam[4..8].copy_from_slice(&pos.to_le_bytes());
        bam[8] = (name.len() + 1) as u8;
        bam[14..16].copy_from_slice(&3u16.to_le_bytes());
        bam[20..24].copy_from_slice(&tid.to_le_bytes());
        bam[24..28].copy_from_slice(&pos.to_le_bytes());
        bam.extend_from_slice(name);
        bam.push(0);
        bam.extend_from_slice(aux);
        bam
    }

    /// `CB:Z:<value>` aux bytes.
    fn cb_aux(value: &[u8]) -> Vec<u8> {
        let mut aux = b"CBZ".to_vec();
        aux.extend_from_slice(value);
        aux.push(0);
        aux
    }

    /// Pack `records` into one segment and describe them as a batch would see
    /// them, returning `(segment, extents)`.
    /// A packed fixture segment and the extents describing the records in it.
    type PackedSegment = (Arc<Vec<u8>>, Box<[(u64, u32)]>);

    fn pack(records: &[Vec<u8>]) -> PackedSegment {
        let mut segment = Vec::new();
        let mut extents = Vec::new();
        for r in records {
            let len = u32::try_from(r.len()).expect("fixture records are small");
            extents.push((segment.len() as u64, len));
            segment.extend_from_slice(r);
        }
        (Arc::new(segment), extents.into_boxed_slice())
    }

    /// A context whose baseline is `records[0]` and which retains every lane, so
    /// no record can trip the dropped-lane check.
    fn context(records: &[Vec<u8>], cell_tag: Option<SamTag>) -> Arc<KeyContext> {
        let lib_lookup = LibraryLookup::from_header(&Header::default());
        let hasher = cb_hasher();
        let first_key = extract_template_key_inline(&records[0], &lib_lookup, cell_tag, &hasher);
        Arc::new(KeyContext {
            lib_lookup,
            cell_tag,
            cb_hasher: hasher,
            first_key,
            variant: TemplateKeyVariant { cb: true, tertiary: true },
        })
    }

    fn run_batch<K: TemplateLaneKey + Send + 'static>(
        ctx: Arc<KeyContext>,
        segment: Arc<Vec<u8>>,
        extents: Box<[(u64, u32)]>,
        first_ref: usize,
    ) -> Option<KeyBatchResult<K>> {
        let (tx, rx) = std::sync::mpsc::channel();
        let job: Box<dyn KeyExtractionJob> = Box::new(TemplateKeyBatch::<K> {
            ctx,
            segment,
            segment_base: 0,
            first_ref,
            extents,
            results: tx,
        });
        job.run();
        rx.try_recv().ok()
    }

    #[test]
    fn test_a_batch_extracts_the_keys_serial_extraction_would_have() {
        // The whole change is only safe if deferring extraction cannot alter the
        // key, so this compares against the serial call the ingest loop used to
        // make, record for record.
        let records: Vec<Vec<u8>> =
            (0..8).map(|i| mapped_bam(1, 100 + i, format!("read{i}").as_bytes(), &[])).collect();
        let ctx = context(&records, None);
        let (segment, extents) = pack(&records);

        let result: KeyBatchResult<TemplateKey40> =
            run_batch(Arc::clone(&ctx), segment, extents, 0).expect("batch publishes a result");

        let expected: Vec<TemplateKey40> = records
            .iter()
            .map(|r| {
                let full =
                    extract_template_key_inline(r, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher);
                TemplateLaneKey::from_full(&full)
            })
            .collect();

        assert_eq!(result.first_ref, 0);
        assert_eq!(result.keys.len(), 8);
        assert_eq!(&result.keys[..], &expected[..]);
        assert!(result.violation.is_none());
    }

    #[test]
    fn test_a_batch_reports_the_lowest_indexed_violation_not_the_last() {
        // Batches run in whatever order the pool picks them up, so the reported
        // record has to be the first offending one or the error message is not
        // reproducible across runs. Two records violate here; the earlier must win.
        let mut records = vec![mapped_bam(1, 100, b"r0", &cb_aux(b"AAAA"))];
        records.push(mapped_bam(1, 101, b"r1", &cb_aux(b"AAAA")));
        records.push(mapped_bam(1, 102, b"r2", &cb_aux(b"CCCC")));
        records.push(mapped_bam(1, 103, b"r3", &cb_aux(b"AAAA")));
        records.push(mapped_bam(1, 104, b"r4", &cb_aux(b"GGGG")));

        // Force the CB lane to be dropped so the differing barcodes violate.
        let lib_lookup = LibraryLookup::from_header(&Header::default());
        let hasher = cb_hasher();
        let first_key =
            extract_template_key_inline(&records[0], &lib_lookup, Some(SamTag::CB), &hasher);
        let ctx = Arc::new(KeyContext {
            lib_lookup,
            cell_tag: Some(SamTag::CB),
            cb_hasher: hasher,
            first_key,
            variant: TemplateKeyVariant { cb: false, tertiary: false },
        });

        let (segment, extents) = pack(&records);
        // `first_ref` is deliberately non-zero: the reported index must be the
        // record's position in the chunk, not in the batch.
        let result: KeyBatchResult<crate::inline::TemplateKey24> =
            run_batch(ctx, segment, extents, 1000).expect("batch publishes a result");

        let v = result.violation.expect("a dropped CB lane must be reported");
        assert_eq!(v.ref_index, 1002, "the first violating record is r2, at chunk index 1002");
        assert_eq!(v.name, "r2");
    }

    #[test]
    fn test_a_multi_segment_ingest_keys_every_record_exactly_once() {
        // The end-to-end contract, and the only test that reaches the
        // dispatch-on-seal path: drive a real buffer and a real pool across
        // several arena segments, then check every key against what the serial
        // extraction this replaced would have produced. A batch that read the
        // wrong segment, double-counted a range, or silently skipped the live
        // segment's tail all show up here as a wrong or placeholder key.
        use crate::codec::SpillCodec;
        use crate::inline::TemplateRecordBuffer;
        use crate::worker_pool::SortWorkerPool;

        // Small enough that a few hundred records fill several segments.
        const SEGMENT: usize = 4096;

        let records: Vec<Vec<u8>> = (0..600i32)
            .map(|i| {
                let name = format!("read{i:05}");
                // Vary length so records straddle segment boundaries unevenly.
                let filler = vec![b'N'; usize::try_from(i % 37).expect("non-negative")];
                let mut aux = b"MCZ".to_vec();
                aux.extend_from_slice(&filler);
                aux.push(0);
                mapped_bam(1, 1000 + (i % 97), name.as_bytes(), &aux)
            })
            .collect();
        let ctx = context(&records, None);

        let pool = SortWorkerPool::new(4, 1, 6, SpillCodec::Bgzf);
        pool.set_phase(crate::worker_pool::phase::PHASE1);

        let mut buffer =
            TemplateRecordBuffer::<TemplateKey40>::with_segment_size(600, SEGMENT, SEGMENT);
        let mut deferred = DeferredKeys::<TemplateKey40>::new(Arc::clone(&ctx), 4);
        for record in &records {
            buffer.push_deferred(record).expect("push");
            deferred.after_push(&mut buffer, &pool);
        }
        assert!(
            buffer.sealed_segments() > 2,
            "the fixture must span several segments to exercise seal-driven dispatch, got {}",
            buffer.sealed_segments(),
        );

        // Most records must have been dispatched *before* the barrier. Without
        // this the test passes with mid-stream dispatch disabled entirely — the
        // barrier would extract everything and the keys would still be right,
        // which is exactly how this optimization could silently do nothing.
        assert_eq!(
            deferred.overlap_census().total_records(),
            0,
            "the census is only taken at the barrier"
        );

        let violation = deferred.finish(&mut buffer, &pool).expect("no worker panics");
        assert!(violation.is_none(), "no record drops a retained lane here");

        let census = deferred.overlap_census();
        let overlapped = census.overlapped_records;
        assert_eq!(census.total_records(), records.len() as u64);
        assert!(
            overlapped > (records.len() as u64) / 2,
            "most keys must be extracted while ingest is still reading; \
             only {overlapped} of {} were",
            records.len(),
        );

        let expected: Vec<TemplateKey40> = records
            .iter()
            .map(|r| {
                let full =
                    extract_template_key_inline(r, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher);
                TemplateLaneKey::from_full(&full)
            })
            .collect();
        let got: Vec<TemplateKey40> = buffer.refs().iter().map(|r| r.key).collect();
        assert_eq!(got, expected, "every record must carry the key serial extraction gives it");

        // And the bytes must still be retrievable at the offsets the refs claim.
        for (i, r) in buffer.refs().iter().enumerate() {
            assert_eq!(buffer.get_record(r), &records[i][..], "record {i} bytes");
        }
    }

    #[test]
    fn test_a_single_worker_pool_keys_inline_and_reports_no_overlap() {
        // With one worker there is no second thread for a batch to run on, so
        // deferral is skipped and extraction happens in the push, as it did
        // before this change. The keys must be identical either way, and the
        // census must say plainly that nothing overlapped rather than going
        // silent — a run reporting no deferral information at all would be
        // indistinguishable from deferral that broke.
        use crate::codec::SpillCodec;
        use crate::inline::TemplateRecordBuffer;
        use crate::worker_pool::SortWorkerPool;

        const SEGMENT: usize = 4096;

        let records: Vec<Vec<u8>> = (0..200i32)
            .map(|i| mapped_bam(1, 500 + i, format!("q{i:04}").as_bytes(), &[]))
            .collect();
        let ctx = context(&records, None);

        let pool = SortWorkerPool::new(1, 1, 6, SpillCodec::Bgzf);
        pool.set_phase(crate::worker_pool::phase::PHASE1);

        let mut buffer =
            TemplateRecordBuffer::<TemplateKey40>::with_segment_size(200, SEGMENT, SEGMENT);
        let mut deferred = DeferredKeys::<TemplateKey40>::new(Arc::clone(&ctx), 1);
        for record in &records {
            deferred.push(&mut buffer, &pool, record).expect("push");
        }
        let violation = deferred.finish(&mut buffer, &pool).expect("no worker panics");
        assert!(violation.is_none());

        let expected: Vec<TemplateKey40> = records
            .iter()
            .map(|r| {
                let full =
                    extract_template_key_inline(r, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher);
                TemplateLaneKey::from_full(&full)
            })
            .collect();
        let got: Vec<TemplateKey40> = buffer.refs().iter().map(|r| r.key).collect();
        assert_eq!(got, expected, "inline extraction must give the same keys");

        let census = deferred.overlap_census();
        assert_eq!(census.total_records(), records.len() as u64, "the census still sums");
        assert_eq!(census.overlap_percent(), Some(0.0), "nothing can overlap with one worker");
    }

    #[test]
    fn test_a_single_worker_pool_still_rejects_a_dropped_lane() {
        // The immediate path has its own error return, so it needs its own
        // coverage: a violation must surface as an error from the push rather
        // than silently keying the record under a lane the sort dropped.
        use crate::codec::SpillCodec;
        use crate::inline::{TemplateKey24, TemplateRecordBuffer};
        use crate::worker_pool::SortWorkerPool;

        let records = [
            mapped_bam(1, 100, b"okay", &cb_aux(b"AAAA")),
            mapped_bam(1, 101, b"bad", &cb_aux(b"TTTT")),
        ];
        let lib_lookup = LibraryLookup::from_header(&Header::default());
        let hasher = cb_hasher();
        let first_key =
            extract_template_key_inline(&records[0], &lib_lookup, Some(SamTag::CB), &hasher);
        let ctx = Arc::new(KeyContext {
            lib_lookup,
            cell_tag: Some(SamTag::CB),
            cb_hasher: hasher,
            first_key,
            variant: TemplateKeyVariant { cb: false, tertiary: false },
        });

        let pool = SortWorkerPool::new(1, 1, 6, SpillCodec::Bgzf);
        pool.set_phase(crate::worker_pool::phase::PHASE1);
        let mut buffer = TemplateRecordBuffer::<TemplateKey24>::with_capacity(2, 4096);
        let mut deferred = DeferredKeys::<TemplateKey24>::new(ctx, 1);

        deferred.push(&mut buffer, &pool, &records[0]).expect("the baseline record is fine");
        let err = deferred
            .push(&mut buffer, &pool, &records[1])
            .expect_err("the offending record must be rejected");
        assert!(err.to_string().contains("bad"), "the error names the record: {err}");
    }

    #[test]
    fn test_a_multi_segment_ingest_reports_a_violation_in_a_sealed_segment() {
        // A dropped-lane violation found by a *worker*, mid-stream, has to reach
        // the ingest thread — the serial path could return the error inline, and
        // this one cannot.
        use crate::codec::SpillCodec;
        use crate::inline::{TemplateKey24, TemplateRecordBuffer};
        use crate::worker_pool::SortWorkerPool;

        const SEGMENT: usize = 4096;

        let mut records: Vec<Vec<u8>> = (0..400)
            .map(|i| mapped_bam(1, 1000 + i, format!("r{i:04}").as_bytes(), &cb_aux(b"AAAA")))
            .collect();
        // One offender, early enough to land in a sealed segment.
        records[7] = mapped_bam(1, 1007, b"offender", &cb_aux(b"TTTT"));

        let lib_lookup = LibraryLookup::from_header(&Header::default());
        let hasher = cb_hasher();
        let first_key =
            extract_template_key_inline(&records[0], &lib_lookup, Some(SamTag::CB), &hasher);
        let ctx = Arc::new(KeyContext {
            lib_lookup,
            cell_tag: Some(SamTag::CB),
            cb_hasher: hasher,
            first_key,
            variant: TemplateKeyVariant { cb: false, tertiary: false },
        });

        let pool = SortWorkerPool::new(4, 1, 6, SpillCodec::Bgzf);
        pool.set_phase(crate::worker_pool::phase::PHASE1);

        let mut buffer =
            TemplateRecordBuffer::<TemplateKey24>::with_segment_size(400, SEGMENT, SEGMENT);
        let mut deferred = DeferredKeys::<TemplateKey24>::new(ctx, 4);
        for record in &records {
            buffer.push_deferred(record).expect("push");
            deferred.after_push(&mut buffer, &pool);
        }
        let violation = deferred
            .finish(&mut buffer, &pool)
            .expect("no worker panics")
            .expect("the offender must be reported");
        assert_eq!(violation.ref_index, 7);
        assert_eq!(violation.name, "offender");
    }

    #[test]
    fn test_a_batch_whose_receiver_is_gone_finishes_quietly() {
        // The ingest thread drops the receiver when it fails, and outstanding
        // batches must not panic on the way out or a real error is replaced by a
        // worker panic.
        let records = vec![mapped_bam(1, 100, b"r0", &[])];
        let ctx = context(&records, None);
        let (segment, extents) = pack(&records);

        let (tx, rx) = std::sync::mpsc::channel::<KeyBatchResult<TemplateKey40>>();
        drop(rx);
        let job: Box<dyn KeyExtractionJob> = Box::new(TemplateKeyBatch::<TemplateKey40> {
            ctx,
            segment,
            segment_base: 0,
            first_ref: 0,
            extents,
            results: tx,
        });
        job.run();
    }

    #[test]
    fn test_a_batch_reads_at_its_segments_base_not_the_arenas() {
        // Extents are global arena offsets while the segment handle is local, so
        // the batch has to subtract its base. Getting this wrong reads the wrong
        // record and silently produces a plausible-but-wrong key rather than
        // failing, which is why it is pinned separately.
        const BASE: u64 = 4096;

        let records: Vec<Vec<u8>> =
            (0..4).map(|i| mapped_bam(1, 200 + i, format!("s{i}").as_bytes(), &[])).collect();
        let ctx = context(&records, None);
        let (segment, local_extents) = pack(&records);
        let shifted: Box<[(u64, u32)]> =
            local_extents.iter().map(|&(o, l)| (o + BASE, l)).collect();

        let (tx, rx) = std::sync::mpsc::channel();
        let job: Box<dyn KeyExtractionJob> = Box::new(TemplateKeyBatch::<TemplateKey40> {
            ctx: Arc::clone(&ctx),
            segment,
            segment_base: BASE,
            first_ref: 0,
            extents: shifted,
            results: tx,
        });
        job.run();
        let result = rx.try_recv().expect("batch publishes a result");

        let expected: Vec<TemplateKey40> = records
            .iter()
            .map(|r| {
                let full =
                    extract_template_key_inline(r, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher);
                TemplateLaneKey::from_full(&full)
            })
            .collect();
        assert_eq!(&result.keys[..], &expected[..]);
    }
}
