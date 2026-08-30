//! `SpillGather` — first step of the block-parallel spill-write split
//! (`SpillGather` → `SpillBlockCompress` → `SpillWrite`), replacing the monolithic
//! single-worker `CompressSpill`.
//!
//! `SpillGather` (`Serial`) consumes the [`SortChunkEvent`]s `SortBuffer`
//! emits and fans each `Spill` chunk into record-aligned **raw** (uncompressed)
//! [`SpillBlockEvent::Block`]s of ≤`BGZF_MAX_BLOCK_SIZE`, so the downstream
//! `Parallel` `SpillBlockCompress` can compress them across the framework pool. The
//! in-memory `Residual` and the terminal `AllAnnounced` pass straight through.
//!
//! # Run extension (spill coalescing)
//!
//! `SpillGather` is also the arena's **run-former**, the coalescing counterpart of
//! the owned engine's `RunFormer`. It sees spill chunks in seal order on a single
//! `Detached` thread, so it can compare each incoming chunk's minimum key against
//! the currently-open run's maximum key ([`RunBound::precedes_or_equal`]): when the
//! chunk is contiguous with the run, its blocks reuse the run's `file_id`, so a
//! whole run of contiguous chunks becomes **one physical spill file** (one merge
//! source) instead of one file per chunk. On an already-sorted input every chunk
//! is contiguous, so the N-way merge collapses to a trivial 1-way pass — matching
//! the owned engine's headline behaviour. The read side is unchanged: a slot is
//! still one file; there are simply fewer, longer files.
//!
//! Whether a run's final block is the *last block of its file* is only known once
//! the next chunk arrives (does it extend?) — so the run's final framed block is
//! **withheld** (`held_last_block`) and closed lazily: emitted with
//! `is_last_in_file = false` and the run's `file_id` if the next chunk extends, or
//! `true` if it does not (or on `Residual` / `AllAnnounced` / drain, which always
//! close the open run). The withheld block's ordinal is minted **at emission**, so
//! the ordinal stream stays dense for `ByItemOrdinal` (a minted-then-dropped
//! ordinal would wedge the single-cursor reorder). `AllAnnounced.slot_count` is
//! rewritten to the **run** count (`runs_written`), since that is how many
//! `SpillReady`s the merge will see and gate on.
//!
//! # Ordinal minting
//!
//! The step mints `ordinal` monotonically across **every** emitted item — every
//! block of every file plus the passthrough `Residual` / `AllAnnounced` — so the
//! output stream is dense and gap-free. That is what lets the framework's
//! single-cursor `ByItemOrdinal` reorder deliver the blocks to the `Serial`
//! `SpillWrite` in order without any per-file ordering primitive. Because
//! `SortBuffer` (Serial) emits `Spill` events one-at-a-time in `seq` order and
//! this step (Serial) drains them in order, each file's blocks are contiguous in
//! the ordinal stream.
//!
//! # Bounded memory (incremental framing)
//!
//! A spill chunk is large (≈ the per-thread sort budget). Framing it **all** into
//! `pending` at once would duplicate the whole chunk in memory and then drain
//! that copy slowly through the byte-bounded queue while `SortBuffer` races ahead
//! filling the next multi-GB buffer — a ~2× peak-RSS blow-up. Instead the chunk
//! is held in `active` and framed **incrementally**: each `try_run` frames at most
//! `MAX_EVENTS_PER_LOCK` blocks into `pending`, drains them, and only frees the
//! source chunk once its last record is framed. `pending` therefore holds ≤ a
//! handful of 64 KiB blocks (~½ MiB) regardless of chunk size, and the chunk
//! drains as fast as `SpillBlockCompress` consumes blocks.

use std::collections::VecDeque;
use std::io;

use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use fgumi_sort::{RunBound, frame_keyed_record_into};
use log::info;

use crate::sort::protocol::{MemoryChunkErased, SortChunkEvent, SpillBlockEvent};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Max staged events flushed to the output per `try_run` invocation. Matches
/// `SortBuffer::MAX_EVENTS_PER_LOCK` so a single fanned chunk drains in bounded
/// slices rather than holding the step lock for the whole spill.
const MAX_EVENTS_PER_LOCK: usize = 8;

/// Frame the `i`th record of a type-erased chunk into `out` in the spill layout
/// `[key?][u32 LE len][record]`, dispatching over the sort-key variant.
fn frame_record_at(chunk: &MemoryChunkErased, i: usize, out: &mut Vec<u8>) -> io::Result<()> {
    match chunk {
        MemoryChunkErased::Coordinate(c) => {
            frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i))
        }
        MemoryChunkErased::QuerynameLex(c) => {
            frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i))
        }
        MemoryChunkErased::QuerynameNatural(c) => {
            frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i))
        }
        MemoryChunkErased::TemplateCoordinate(c) => c.frame_record_into(i, out),
    }
}

/// Frame records `[start..)` of `chunk` into a single block `out` of at most
/// `block_size` bytes, returning the index of the next un-framed record (==
/// `chunk.len()` when the chunk is exhausted).
///
/// A record is never split across blocks: a record that would push a **non-empty**
/// block past `block_size` is left for the next block. The sole exception is a
/// single record larger than `block_size`, which forms its own oversized block
/// (BAM long reads can exceed 64 KiB; the codec re-blocks/​frames it safely,
/// matching the streaming `SyncSpillWriter`, so this is not an error).
fn frame_one_block(
    chunk: &MemoryChunkErased,
    start: usize,
    block_size: usize,
    out: &mut Vec<u8>,
) -> io::Result<usize> {
    let len = chunk.len();
    let mut i = start;
    while i < len {
        let before = out.len();
        frame_record_at(chunk, i, out)?;
        if before != 0 && out.len() > block_size {
            // This record overflowed a non-empty block — roll it back so it
            // starts the next block, and finish this one.
            out.truncate(before);
            break;
        }
        i += 1;
        if out.len() >= block_size {
            // Block full (the record fit exactly, or was a lone oversized record
            // framed into an empty block).
            break;
        }
    }
    Ok(i)
}

/// Frame an entire chunk into blocks (convenience for tests; production frames
/// incrementally via [`SpillGather::produce_blocks`]). Returns one `Vec<u8>`
/// per block; an empty chunk yields no blocks.
#[cfg(test)]
fn frame_chunk_into_blocks(
    chunk: &MemoryChunkErased,
    block_size: usize,
) -> io::Result<Vec<Vec<u8>>> {
    let mut blocks: Vec<Vec<u8>> = Vec::new();
    let len = chunk.len();
    let mut idx = 0;
    while idx < len {
        let mut block = Vec::with_capacity(block_size + 1024);
        idx = frame_one_block(chunk, idx, block_size, &mut block)?;
        if !block.is_empty() {
            blocks.push(block);
        }
    }
    Ok(blocks)
}

/// One spill chunk being framed incrementally into blocks across `try_run` calls.
struct ActiveSpill {
    /// The source sorted chunk, held until its last record is framed.
    chunk: MemoryChunkErased,
    /// Index of the next un-framed record.
    next_idx: usize,
    /// This chunk's run `file_id` — the open run's when the chunk extends it, or
    /// the chunk's own spill `seq` when it starts a new run.
    file_id: u32,
    records_ingested_so_far: u64,
}

/// A run's most-recently-framed block, withheld until we know whether it is the
/// last block *of its file*. Carries everything a [`SpillBlockEvent::Block`]
/// needs except `is_last_in_file` (decided at flush) and `ordinal` (minted at
/// flush, to keep the emitted stream dense).
struct HeldBlock {
    bytes: Vec<u8>,
    file_id: u32,
    records_ingested_so_far: u64,
}

/// `Serial` step that fans sorted spill chunks into raw blocks for `SpillBlockCompress`.
pub struct SpillGather {
    /// The spill chunk currently being framed incrementally (`None` between
    /// chunks). Holding it here — rather than materializing all its blocks into
    /// `pending` — is what keeps peak memory bounded.
    active: Option<ActiveSpill>,
    /// Staged block / passthrough events awaiting output. Bounded to ≤
    /// `MAX_EVENTS_PER_LOCK` blocks because framing only tops it up when it is
    /// already drained.
    pending: VecDeque<SpillBlockEvent>,
    /// Monotonic ordinal minted across every emitted item (dense, gap-free).
    next_ordinal: u64,
    held: HeldSlot<Unpushed<SpillBlockEvent>>,
    /// The currently-open run: its `file_id` and its maximum key so far. An
    /// incoming chunk whose minimum key is `>=` this maximum *extends* the run
    /// (reuses `file_id`); otherwise it starts a new run. `None` before the first
    /// spill and after a `Residual` / `AllAnnounced` closes the run.
    open_run: Option<(u32, RunBound)>,
    /// The open run's final framed block, withheld pending the extend/close
    /// decision of the next event (lazy run-close). At most one 64 KiB block.
    held_last_block: Option<HeldBlock>,
    /// Number of distinct runs (spill files) formed. Rewritten into
    /// `AllAnnounced.slot_count`, since that is how many `SpillReady`s the merge
    /// will receive and gate on.
    runs_written: u32,
    /// Number of non-empty spill chunks seen. With `runs_written`, gives the
    /// owned-style run-formation report (`extended = chunks_spilled - runs`).
    chunks_spilled: u32,
    /// Whether the run-formation summary has been logged (guards the single-shot
    /// log at the drained `Finished` transition).
    summary_logged: bool,
    /// Raw-block size threshold (records are cut into blocks of ≤ this many bytes).
    block_size: usize,
    output_byte_limit: u64,
}

impl SpillGather {
    /// Build a `SpillGather`. `output_byte_limit` byte-bounds the block-event
    /// output queue (its `Block` / `Residual` variants retain bytes/records).
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            active: None,
            pending: VecDeque::new(),
            next_ordinal: 0,
            held: HeldSlot::new(),
            open_run: None,
            held_last_block: None,
            runs_written: 0,
            chunks_spilled: 0,
            summary_logged: false,
            block_size: BGZF_MAX_BLOCK_SIZE,
            output_byte_limit,
        }
    }

    /// Take the next dense ordinal.
    fn next_ordinal(&mut self) -> u64 {
        let o = self.next_ordinal;
        self.next_ordinal += 1;
        o
    }

    /// The owned-engine-style run-formation summary line, or `None` when nothing
    /// spilled. Mirrors `fgumi_sort::external`'s `log_run_formation` so both sort
    /// engines report run formation identically: every spilled chunk either
    /// starts a run or extends one, so `extended = chunks_spilled - runs_written`.
    fn run_formation_summary(&self) -> Option<String> {
        if self.chunks_spilled == 0 {
            return None;
        }
        let extended = self.chunks_spilled.saturating_sub(self.runs_written);
        Some(format!(
            "Spill runs: {} from {} chunks ({extended} extended an existing run)",
            self.runs_written, self.chunks_spilled
        ))
    }

    /// Flush the withheld run-final block (if any) into `pending`, marking it
    /// `is_last_in_file = is_last` and minting its ordinal now (at emission, so
    /// the stream stays dense). No-op when nothing is held.
    fn flush_held_block(&mut self, is_last: bool) {
        let Some(hb) = self.held_last_block.take() else { return };
        let ordinal = self.next_ordinal();
        self.pending.push_back(SpillBlockEvent::Block {
            ordinal,
            file_id: hb.file_id,
            is_last_in_file: is_last,
            records_ingested_so_far: hb.records_ingested_so_far,
            bytes: hb.bytes,
        });
    }

    /// `StepCtx`-free core: stage one input event. A `Spill` chunk runs the
    /// run-former (extend-or-close-and-open) and is parked in `active` for
    /// incremental framing (no blocks produced yet). `Residual` / `AllAnnounced`
    /// close the open run (flushing its withheld final block) and pass straight
    /// through; `AllAnnounced.slot_count` is rewritten to the run count. Caller
    /// guarantees `active` is `None` (the previous chunk fully framed).
    fn stage_event(&mut self, event: SortChunkEvent) {
        match event {
            SortChunkEvent::Spill { seq, chunk, records_ingested_so_far } => {
                // `SortBuffer` only emits `Spill` for a non-empty buffer; skip an
                // (unexpected) empty chunk rather than open a zero-block file with
                // no `is_last_in_file` terminator. An empty chunk is a pure no-op:
                // it contributes no bounds and does NOT close the open run.
                if chunk.is_empty() {
                    return;
                }
                self.chunks_spilled += 1;
                let chunk_min = chunk.min_bound().expect("a non-empty chunk has a min bound");
                let chunk_max = chunk.max_bound().expect("a non-empty chunk has a max bound");
                // Extends the open run iff every record already in the run precedes
                // or equals this chunk's first record.
                let extends = match &self.open_run {
                    Some((_, open_max)) => open_max.precedes_or_equal(&chunk_min),
                    None => false,
                };
                // Close-or-continue the previous run's withheld final block: it is
                // the last block of its file iff this chunk does NOT extend it.
                self.flush_held_block(!extends);
                let file_id = if extends {
                    self.open_run.as_ref().expect("extends implies an open run").0
                } else {
                    self.runs_written += 1;
                    seq
                };
                self.open_run = Some((file_id, chunk_max));
                self.active =
                    Some(ActiveSpill { chunk, next_idx: 0, file_id, records_ingested_so_far });
            }
            SortChunkEvent::Residual { chunk, records_ingested_so_far } => {
                // A residual (in-memory, not spilled) ends any open spill run.
                self.flush_held_block(true);
                self.open_run = None;
                let ordinal = self.next_ordinal();
                self.pending.push_back(SpillBlockEvent::Residual {
                    ordinal,
                    chunk,
                    records_ingested_so_far,
                });
            }
            SortChunkEvent::AllAnnounced { memory_chunk_count, total_records, .. } => {
                // Terminal: close the final run, then announce the RUN count as the
                // slot count (the merge gates on `slots.len() == slot_count`, and
                // with coalescing #SpillReady = #runs, not #chunks).
                self.flush_held_block(true);
                self.open_run = None;
                let ordinal = self.next_ordinal();
                self.pending.push_back(SpillBlockEvent::AllAnnounced {
                    ordinal,
                    slot_count: self.runs_written,
                    memory_chunk_count,
                    total_records,
                });
            }
        }
    }

    /// Frame up to `MAX_EVENTS_PER_LOCK` non-final blocks of the `active` chunk
    /// into `pending`, minting dense ordinals. The chunk's **final** framed block
    /// is not pushed: it is withheld into `held_last_block` (its `is_last_in_file`
    /// is not yet known — it depends on whether the next chunk extends this run),
    /// and the chunk is freed (`active = None`). No-op when there is no active
    /// chunk.
    ///
    /// # Errors
    ///
    /// Propagates framing errors (e.g. a record too large for the `u32` length
    /// prefix).
    fn produce_blocks(&mut self) -> io::Result<()> {
        let block_size = self.block_size;
        while self.pending.len() < MAX_EVENTS_PER_LOCK {
            // Frame one block, scoping the borrow of `active` so `next_ordinal`
            // (a `&mut self` method) can run afterward.
            let framed = {
                let Some(active) = self.active.as_ref() else { return Ok(()) };
                let len = active.chunk.len();
                let mut bytes = Vec::with_capacity(block_size + 1024);
                let next = frame_one_block(&active.chunk, active.next_idx, block_size, &mut bytes)?;
                (bytes, next, next >= len, active.file_id, active.records_ingested_so_far)
            };
            let (bytes, next_idx, is_last, file_id, records_ingested_so_far) = framed;
            if is_last {
                // Withhold the run's final block: it is closed lazily once we know
                // whether the next chunk extends this run (see `stage_event` /
                // `flush_held_block`). The previous run's block was already flushed
                // before this chunk was staged, so nothing is being overwritten.
                debug_assert!(
                    self.held_last_block.is_none(),
                    "a prior held block must be flushed before withholding a new one"
                );
                self.held_last_block = Some(HeldBlock { bytes, file_id, records_ingested_so_far });
                self.active = None;
                return Ok(());
            }
            let ordinal = self.next_ordinal();
            self.pending.push_back(SpillBlockEvent::Block {
                ordinal,
                file_id,
                is_last_in_file: false,
                records_ingested_so_far,
                bytes,
            });
            // Advance the cursor for the next block.
            self.active.as_mut().expect("active present (not last)").next_idx = next_idx;
        }
        Ok(())
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    /// Push up to `MAX_EVENTS_PER_LOCK` staged events, parking the first that
    /// can't be pushed in `held`. Caller guarantees `held` is empty.
    fn emit_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let mut emitted = 0usize;
        while emitted < MAX_EVENTS_PER_LOCK {
            let Some(event) = self.pending.pop_front() else { break };
            if let Err(unpushed) = ctx.outputs.push(event) {
                self.held.put(unpushed);
                return StepOutcome::Progress;
            }
            emitted += 1;
        }
        if emitted > 0 { StepOutcome::Progress } else { StepOutcome::NoProgress }
    }
}

impl Step for SpillGather {
    type Input = SortChunkEvent;
    type Outputs = OrderedBytesSingle<SpillBlockEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillGather",
            // Off-pool on the coordination driver (N+2): the serial spill framing
            // + monotonic ordinal minting runs on the dedicated coordination
            // thread instead of a pool worker, so it never starves the parallel
            // spill compressors. Detached collapses `ByItemOrdinal` to `None`
            // exactly as `Serial` did (transport-identical).
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP)
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Drain staged events before producing/ingesting more — bounds peak memory.
        if !self.pending.is_empty() {
            return Ok(self.emit_pending(ctx));
        }

        // Continue framing the active chunk (incrementally), if any.
        if self.active.is_some() {
            self.produce_blocks()?;
            return Ok(self.emit_pending(ctx));
        }

        // No active chunk and nothing pending — ingest the next event.
        if let Some(event) = ctx.input.pop() {
            self.stage_event(event);
            // Frame the first blocks of a freshly-parked chunk so we make progress.
            if self.active.is_some() {
                self.produce_blocks()?;
            }
            if !self.pending.is_empty() {
                return Ok(self.emit_pending(ctx));
            }
            return Ok(StepOutcome::Progress);
        }

        // Drained: pending is empty and no chunk is mid-framing (checked above).
        if ctx.input.is_drained() {
            // `AllAnnounced` (always the terminal event) closes the open run and
            // flushes its withheld final block. A block still held at drain means
            // that close never happened — a lost run terminator would leave
            // `SpillWrite` with a dangling open file. Fail loudly rather than
            // silently drop a run.
            if self.held_last_block.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "SpillGather drained with an unclosed run (a withheld spill block was \
                     never flushed by AllAnnounced)",
                ));
            }
            // Report run formation once, matching the owned engine's line so
            // benchmark log parsing and cross-engine comparison see the same
            // "Spill runs: R from C chunks (E extended…)" wording.
            if !self.summary_logged {
                self.summary_logged = true;
                if let Some(line) = self.run_formation_summary() {
                    info!("{line}");
                }
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests;
