//! `SortSpillDecompress` — Parallel typed step that reads spill chunk
//! files, decompresses them inline (codec-aware: BGZF blocks or zstd
//! frames per `slot.codec`, detected from the file magic at slot-open),
//! and pushes the decompressed bytes into per-slot bounded queues on
//! `SortMergeSlot`. Replaces the legacy +N `SortWorkerPool::Phase2`
//! OS threads with framework workers from the unified pipeline's
//! work-stealing pool.
//!
//! See `docs/design/sort-step-split-parity-fix.md` for the v4 design.
//!
//! ## Step body (`try_run`)
//!
//! 1. **Drain held output** if any (legacy `HeldSlot` pattern — pushes
//!    that failed last call due to a full output queue retry first).
//! 2. **Pop one input event** if available:
//!    * `SpillReady` → register the slot in the shared registry, AND
//!      forward verbatim to `SortMerge` for slot-table install.
//!    * `MemoryChunk` → forward verbatim to `SortMerge`.
//!    * `AllAnnounced` → forward verbatim.
//! 3. **Greedy slot fill** (`try_fill_some_slot`): walk a snapshot of
//!    the registry; for the first slot that has room in its
//!    `decompressed` queue (i.e., `decompressed.len() < CAP`) AND a
//!    free reader (`reader.try_lock()` succeeds), read up to
//!    `min(CAP - len, MAX_BATCH)` raw BGZF blocks from disk,
//!    decompress them inline on this worker, and push the batch into
//!    `decompressed`. If the read returns fewer than requested, set
//!    `queue_eof = true` (last batch).
//! 4. If no slot accepted reader/fill work: return `Contention` if
//!    ANY registered slot is still not `queue_eof` (consumer may
//!    drain space later); otherwise, once the input edge is drained
//!    too, report `Finished`.
//!
//! Producer is **non-blocking** — never waits on consumer activity.
//! The consumer (`MergeDriver::try_step` via `slot_try_load_block`) is
//! also non-blocking: when a slot's queue is empty and not yet
//! `queue_eof` it reports `WouldBlock` and the cooperative `SortMerge`
//! step yields, so the framework re-dispatches it after this producer
//! has pushed more blocks or set `queue_eof`. No condvar / notify.
//!
//! ## Atomic ordering for `queue_eof` / `decomp_error`
//!
//! Producers' BOTH success and error paths hold `slot.decompressed`
//! while storing the atomics. Consumer reads them under the same
//! mutex. The mutex release-acquire chain establishes happens-before
//! for both atomics simultaneously; consumer's load order does not
//! matter for correctness.

use std::io;
use std::sync::Arc;

use fgumi_sort::{PHASE2_DECOMP_CAP, SortMergeSlot, SpillBlockDecompressor};
use parking_lot::Mutex;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::Single;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::sort::protocol::{SortPhase1Event, SortPhase2Event};

/// Max raw blocks read+decompressed by a single `try_fill_some_slot`
/// call. Bounds the per-call work so the framework's round-robin
/// dispatch stays responsive (other steps don't starve while this
/// worker hogs a slot). At ~ms-per-block decompress and N=4, this is
/// ~4ms of work per call.
const MAX_BATCH_PER_CALL: usize = 4;

/// A slot announced via `SpillReady` and now eligible for greedy fill.
/// Stored in the step's shared registry; all Parallel worker copies
/// read it via the same `Arc`.
struct RegisteredSpill {
    slot: Arc<SortMergeSlot>,
}

/// Parallel step that reads + decompresses spill chunk files and
/// pushes results into per-slot queues. Forwards `SortPhase1Event`s
/// verbatim to `SortMerge`.
///
/// All Parallel worker copies share:
///
/// * `registry: Arc<Mutex<Vec<RegisteredSpill>>>` — files announced
///   via `SpillReady`. Append-only; per-worker `try_run` walks this
///   list when looking for fill work.
///
/// Per-worker:
///
/// * `block_dec: SpillBlockDecompressor` — fresh per worker copy; the
///   codec-aware (BGZF block / zstd frame) decompressor that reads +
///   inflates a batch of blocks and returns owned `Vec<u8>`s for the
///   per-slot queue (mimalloc size-class buffer reuse inside).
/// * `held: HeldSlot<Unpushed<SortPhase2Event>>` — backpressure retry
///   slot for the forwarded events (`SpillReady`/`MemoryChunk`/
///   `AllAnnounced`). The per-slot decompressed queue handles its
///   own backpressure separately (producer SKIPs when at cap).
pub struct SortSpillDecompress {
    registry: Arc<Mutex<Vec<RegisteredSpill>>>,
    /// Per-worker codec-aware decompressor (BGZF blocks or zstd frames).
    block_dec: SpillBlockDecompressor,
    held: HeldSlot<Unpushed<SortPhase2Event>>,
    /// Max in-flight forwarded events. `SortSpillDecompress` emits at
    /// most one event per `try_run` (forward from input pop). 64
    /// covers any realistic spill count.
    output_capacity: usize,
}

impl SortSpillDecompress {
    /// Construct a fresh step with an empty registry.
    #[must_use]
    pub fn new(output_capacity: usize) -> Self {
        Self {
            registry: Arc::new(Mutex::new(Vec::new())),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_capacity,
        }
    }

    /// Attempt to deliver any held output event. Returns `true` if
    /// the held slot is now empty.
    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        let Some(unpushed) = self.held.take() else {
            return true;
        };
        match ctx.outputs.retry(unpushed) {
            Ok(()) => true,
            Err(again) => {
                self.held.put(again);
                false
            }
        }
    }

    /// Push an output event, stashing it on backpressure. Returns
    /// `true` on success.
    fn push_or_hold(&mut self, ctx: &mut StepCtx<'_, Self>, event: SortPhase2Event) -> bool {
        match ctx.outputs.push(event) {
            Ok(()) => true,
            Err(unpushed) => {
                self.held.put(unpushed);
                false
            }
        }
    }

    /// Snapshot the registry's slots so the caller can walk them
    /// without holding the registry mutex across disk I/O or
    /// decompression.
    fn snapshot_registry(&self) -> Vec<Arc<SortMergeSlot>> {
        let registry = self.registry.lock();
        registry.iter().map(|e| Arc::clone(&e.slot)).collect()
    }

    /// Greedy slot-fill step. Walk registered slots; for the first
    /// one with `decompressed.len() < CAP` AND a free reader, read +
    /// decompress + push a batch of up to `MAX_BATCH_PER_CALL` blocks.
    /// Returns `Ok(true)` if any work happened.
    ///
    /// **Atomic ordering**: the success path sets `queue_eof = true`
    /// (when reader hits EOF) WHILE HOLDING the `decompressed` mutex,
    /// so the consumer's next lock-acquire synchronizes-with this
    /// release and sees both the pushed batch AND the eof flag.
    /// The error path same discipline (sets `decomp_error` +
    /// `queue_eof` while holding `decompressed`).
    fn try_fill_some_slot(&mut self) -> io::Result<bool> {
        for slot in self.snapshot_registry() {
            // Skip slots already done producing.
            if slot.queue_eof.load(std::sync::atomic::Ordering::Acquire) {
                continue;
            }

            // Try to acquire the slot's reader. If contended, move
            // to the next slot — another worker is feeding this one.
            let Ok(mut reader_guard) = slot.reader.try_lock() else {
                continue;
            };

            // Compute available queue space M = CAP - decompressed.len().
            // Brief decompressed-lock for the len read; skip slot if
            // queue is at cap (consumer must drain first).
            let m = {
                let dec = slot.decompressed.lock().expect("decompressed mutex poisoned");
                PHASE2_DECOMP_CAP.saturating_sub(dec.len())
            };
            if m == 0 {
                // Queue at cap. Drop reader and try another slot.
                continue;
            }
            let m = m.min(MAX_BATCH_PER_CALL);

            // Read + decompress up to m blocks from disk, codec-aware: BGZF
            // self-framed blocks or zstd `[len][frame]` records, per `slot.codec`
            // (detected from the file magic at slot-open). `read_blocks` returns
            // fewer than `m` (incl. 0) at EOF.
            //
            // **We hold `reader_guard` through the whole read+decompress + push**
            // for this batch — NOT just through the read. This is the
            // "I'm the producer for this slot" lock: while held, no sibling
            // worker can claim this slot, read another batch, and prematurely set
            // `queue_eof` while our batch is still in flight. Releasing the reader
            // before decompression would let a sibling worker see an empty read on
            // a slot we'd already exhausted and set `queue_eof` before our batch
            // lands — silently dropping our decompressed records.
            //
            // Cost: less per-slot parallelism (one worker decompresses per slot
            // at a time). Cross-slot parallelism preserved — different workers
            // handle different slots.
            let decompressed_batch =
                match self.block_dec.read_blocks(&mut reader_guard.inner, slot.codec, m) {
                    Ok(b) => b,
                    Err(e) => {
                        // Mark the slot as errored + EOF so the consumer exits
                        // cleanly with Err. Hold `reader_guard` through the flag
                        // block (lock order: reader → decompressed, same as the
                        // success paths below) and drop it AFTER. Dropping the
                        // reader first would let a sibling worker claim the slot,
                        // read a clean EOF, and set `queue_eof` WITHOUT
                        // `decomp_error` in the window — which the consumer (checks
                        // `decomp_error` then `queue_eof`) would see as a clean
                        // drain, silently dropping this error and the records after
                        // it.
                        {
                            let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
                            slot.decomp_error.store(true, std::sync::atomic::Ordering::Release);
                            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
                        }
                        drop(reader_guard);
                        return Err(e);
                    }
                };
            let got = decompressed_batch.len();
            let hit_eof = got < m;

            if got == 0 {
                // Empty read at EOF. Set queue_eof under the decompressed lock
                // (still holding reader so no sibling can race).
                {
                    let _g = slot.decompressed.lock().expect("decompressed mutex poisoned");
                    slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
                }
                drop(reader_guard);
                return Ok(true);
            }

            // Push all decompressed blocks under one decompressed-lock
            // acquisition. Set queue_eof under the same lock if this
            // batch is the last (hit_eof). Drop reader AFTER releasing
            // decompressed (lock order: reader → decompressed).
            {
                let mut dec = slot.decompressed.lock().expect("decompressed mutex poisoned");
                for b in decompressed_batch {
                    dec.push_back(b);
                }
                if hit_eof {
                    slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
                }
            }
            drop(reader_guard);
            return Ok(true);
        }
        Ok(false)
    }
}

impl Clone for SortSpillDecompress {
    fn clone(&self) -> Self {
        Self {
            // Per Parallel-step semantics: each worker copy shares
            // the registry so SpillReady popped by any worker is
            // immediately visible to all others.
            registry: Arc::clone(&self.registry),
            block_dec: SpillBlockDecompressor::new(),
            held: HeldSlot::new(),
            output_capacity: self.output_capacity,
        }
    }
}

impl Step for SortSpillDecompress {
    type Input = SortPhase1Event;
    type Outputs = Single<SortPhase2Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortSpillDecompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: self.output_capacity }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain held output first.
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // 2. Pop one input event, register if SpillReady, forward all.
        if let Some(event) = ctx.input.pop() {
            let forwarded = match event {
                SortPhase1Event::SpillReady { slot, path, records_ingested_so_far } => {
                    self.registry.lock().push(RegisteredSpill { slot: Arc::clone(&slot) });
                    SortPhase2Event::SpillReady { slot, path, records_ingested_so_far }
                }
                SortPhase1Event::MemoryChunk { chunk, records_ingested_so_far } => {
                    SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far }
                }
                SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, total_records } => {
                    SortPhase2Event::AllAnnounced { slot_count, memory_chunk_count, total_records }
                }
            };
            let _ = self.push_or_hold(ctx, forwarded);
            return Ok(StepOutcome::Progress);
        }

        // 3. Greedy slot-fill step.
        if self.try_fill_some_slot()? {
            return Ok(StepOutcome::Progress);
        }

        // 4. No fill work. Stay alive (return Contention) while ANY
        //    slot is still producing (i.e., !queue_eof). Otherwise, if the
        //    input edge is drained too, this step is done — report Finished.
        //
        //    See `docs/design/sort-step-split-parity-fix.md` Change
        //    2 — staying alive (Contention) while slots are producing is
        //    what prevents workers from Skipping prematurely during the
        //    cap-bound transient window that deadlocked v3.1.
        let any_alive = self
            .snapshot_registry()
            .iter()
            .any(|slot| !slot.queue_eof.load(std::sync::atomic::Ordering::Acquire));
        if any_alive {
            return Ok(StepOutcome::Contention);
        }

        // No fill work, no input this call, no slot still producing. If the
        // input edge is drained too, this step will never push again — report
        // Finished (only the last Parallel clone closes the shared output,
        // gated by the StepDrainCounter in the driver). Otherwise wait for
        // more SpillReady events.
        if ctx.input.is_drained() {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests;
