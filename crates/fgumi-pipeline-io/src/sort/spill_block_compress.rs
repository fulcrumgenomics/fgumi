//! `SpillBlockCompress` — middle step of the block-parallel spill-write split
//! (`SpillGather` → `SpillBlockCompress` → `SpillWrite`).
//!
//! `SpillBlockCompress` (`Parallel + ByItemOrdinal`) turns each planned
//! [`SpillBlockEvent::Range`] from `SpillGather` into a compressed
//! [`SpillBlockEvent::Block`]: it frames the range's records into a per-worker
//! scratch buffer and compresses that into a self-contained unit (framed BGZF
//! block(s) for bgzf, or a `[u32 len][zstd frame]` for zstd) via the shared
//! [`SpillBlockCompressor`] kernel. A raw `Block` (a run's final block, framed
//! by `SpillGather`) is compressed in place. `Residual` and `AllAnnounced` pass straight through. The
//! `ordinal` is preserved on every event, so the framework's `ByItemOrdinal`
//! reorder hands `SpillWrite` a dense, in-order stream regardless of which worker
//! compressed which block — exactly the output path's `BgzfCompress` idiom.
//!
//! Each `Parallel` worker holds its own [`SpillBlockCompressor`], built lazily on
//! the first block (the zstd compressor's construction is fallible, so it is
//! surfaced through `try_run`'s `io::Result` rather than the infallible
//! `new_worker_copy`).

use std::io;
use std::sync::Arc;

use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use fgumi_sort::{SpillBlockCompressor, SpillCodec};

use crate::sort::protocol::SpillBlockEvent;
use fgumi_pipeline_core::{
    HeldRetry, PhaseCap, Unpushed, admit_input,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Initial capacity of a worker's raw framing scratch: one full block plus
/// room for the record that ends it.
const RAW_SCRATCH_BYTES: usize = BGZF_MAX_BLOCK_SIZE + 1024;

/// Most compressed-block buffers [`SpillBufPool`] keeps for reuse.
const SPILL_BUF_POOL_CAP: usize = 256;

/// Free list of compressed spill-block buffers, shared by every
/// `SpillBlockCompress` worker (which take) and `SpillWrite` (which returns
/// each block's buffer once its bytes are in the file). Reuse removes the
/// per-block allocation and the cross-thread free on the writer. Bounded to
/// `SPILL_BUF_POOL_CAP` (256) buffers; past that a returned buffer is dropped.
#[derive(Debug, Default)]
pub struct SpillBufPool {
    free: parking_lot::Mutex<Vec<Vec<u8>>>,
}

impl SpillBufPool {
    /// A new, empty pool.
    #[must_use]
    pub fn new() -> Arc<Self> {
        Arc::new(Self::default())
    }

    /// A cleared buffer: a recycled one when available, else a new one.
    #[must_use]
    pub fn take(&self) -> Vec<u8> {
        self.free.lock().pop().unwrap_or_default()
    }

    /// Return `buf` for reuse (dropped if the pool is full).
    pub fn give(&self, mut buf: Vec<u8>) {
        if buf.capacity() == 0 {
            return;
        }
        buf.clear();
        let mut free = self.free.lock();
        if free.len() < SPILL_BUF_POOL_CAP {
            free.push(buf);
        }
    }
}

/// `Parallel + ByItemOrdinal` block compressor for the spill-write split: the
/// middle step of `SpillGather → SpillBlockCompress → SpillWrite`. It compresses
/// each block on a pool worker; the disk write is the separate `SpillWrite`
/// step.
pub struct SpillBlockCompress {
    codec: SpillCodec,
    compression: u32,
    /// Per-worker compressor, built lazily on the first block. `None` until then
    /// (and on fresh `new_worker_copy` clones).
    compressor: Option<SpillBlockCompressor>,
    /// Per-worker raw framing scratch, retained across blocks: a range is
    /// framed here and compressed from here, so no raw block is allocated per
    /// block. Grown by a lone oversized record and kept.
    raw: Vec<u8>,
    /// Per-worker scratch: where each frame of a multi-frame range ends in `raw`.
    frame_ends: Vec<usize>,
    held: HeldSlot<Unpushed<SpillBlockEvent>>,
    output_byte_limit: u64,
    /// Admission cap shared with the other phase-1 pool steps (set by the
    /// chain builder via `with_phase_cap`; `None` = uncapped).
    cap: Option<Arc<PhaseCap>>,
    /// Where compressed-block buffers come from (shared with `SpillWrite`).
    buf_pool: Option<Arc<SpillBufPool>>,
}

impl SpillBlockCompress {
    /// Build a `SpillBlockCompress` for `codec` at `compression`. `output_byte_limit`
    /// byte-bounds the compressed-block output queue.
    #[must_use]
    pub fn new(codec: SpillCodec, compression: u32, output_byte_limit: u64) -> Self {
        Self {
            codec,
            compression,
            compressor: None,
            raw: Vec::with_capacity(RAW_SCRATCH_BYTES),
            frame_ends: Vec::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            cap: None,
            buf_pool: None,
        }
    }

    /// The raw scratch's address and capacity (tests: pointer stability shows
    /// the buffer is reused, not reallocated).
    #[cfg(test)]
    pub(crate) fn raw_scratch(&self) -> (*const u8, usize) {
        (self.raw.as_ptr(), self.raw.capacity())
    }

    /// This worker's compressor, built on first use.
    fn compressor(&mut self) -> io::Result<&mut SpillBlockCompressor> {
        if self.compressor.is_none() {
            self.compressor = Some(SpillBlockCompressor::new(self.codec, self.compression)?);
        }
        Ok(self.compressor.as_mut().expect("compressor built above"))
    }

    /// Take compressed-block buffers from `pool` (shared with `SpillWrite`).
    #[must_use]
    pub fn with_buf_pool(mut self, pool: Arc<SpillBufPool>) -> Self {
        self.buf_pool = Some(pool);
        self
    }

    /// A buffer for one compressed block.
    fn out_buf(&self) -> Vec<u8> {
        self.buf_pool.as_ref().map(|p| p.take()).unwrap_or_default()
    }

    /// Share a phase admission cap with the other pool steps of this phase. At
    /// most `cap.max()` workers are inside any of those steps at once.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// One `try_run` with the phase cap moved out of `self` (see `try_run`).
    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Admission (see the pipeline-core `admission` module): a capped step
        // takes its phase permit BEFORE popping (an empty poll takes none and
        // is recorded as an empty pop), so a refused clone leaves the item for
        // a clone that may run and reports `Capped`. Uncapped, this is a no-op
        // and the pop below runs exactly as without a cap. The permit is held
        // to the end of this call.
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };

        if let Some(event) = ctx.input.pop() {
            let forwarded = self.compress_event(event)?;
            if let Err(unpushed) = ctx.outputs.push(forwarded) {
                self.held.put(unpushed);
            }
            return Ok(StepOutcome::Progress);
        }

        if ctx.input.is_drained() {
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        !matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld)
    }

    /// Frame and compress a `Range`, compress a raw `Block`, or pass a
    /// `Residual`/`AllAnnounced` through. The `ordinal` and routing fields are
    /// preserved. `StepCtx`-free for unit tests.
    ///
    /// # Errors
    ///
    /// Propagates framing, compressor-init or compression errors.
    pub(crate) fn compress_event(&mut self, event: SpillBlockEvent) -> io::Result<SpillBlockEvent> {
        match event {
            SpillBlockEvent::Range {
                ordinal,
                file_id,
                key_kind,
                is_last_in_file,
                records_ingested_so_far,
                chunk,
                start,
                end,
                framed_bytes,
                frame_bytes,
            } => {
                // `SpillGather` plans `start <= end <= chunk.len()`, and its cut
                // rule and the framing agree on the bytes. Both are checked in
                // every build, as the serial run-final path checks its block:
                // a disagreement fails the sort instead of panicking or
                // spilling a block the planner did not account for.
                let len = chunk.len() as u64;
                if start > end || end > len {
                    return Err(io::Error::other(format!(
                        "SpillBlockCompress: range [{start}, {end}) is not within the chunk's \
                         {len} records"
                    )));
                }
                let end_idx = usize::try_from(end).expect("a record index fits usize");
                let mut frame_start = usize::try_from(start).expect("a record index fits usize");
                // Frame every frame of the range into `raw`, noting where each
                // ends. A multi-frame range is re-cut with the planner's own rule,
                // so each frame is exactly the block it would have been as a
                // range of its own (and the spill bytes do not change).
                self.raw.clear();
                self.frame_ends.clear();
                while frame_start < end_idx {
                    let frame_end = if frame_bytes == 0 {
                        end_idx
                    } else {
                        crate::sort::spill_gather::plan_next_block(
                            &chunk,
                            frame_start,
                            frame_bytes as usize,
                        )
                        .0
                        .min(end_idx)
                    };
                    for i in frame_start..frame_end {
                        chunk.frame_record_into(i, &mut self.raw)?;
                    }
                    self.frame_ends.push(self.raw.len());
                    frame_start = frame_end;
                }
                if self.raw.len() != framed_bytes as usize {
                    return Err(io::Error::other(format!(
                        "SpillBlockCompress: range [{start}, {end}) framed {} bytes but the \
                         planner cut {framed_bytes}",
                        self.raw.len()
                    )));
                }
                // The framed bytes are all this needs from the chunk: release it
                // now, so the last range of a run returns the arena before the
                // compression, not after.
                drop(chunk);
                let raw = std::mem::take(&mut self.raw);
                let frame_ends = std::mem::take(&mut self.frame_ends);
                let mut out = self.out_buf();
                let compressed = self.compressor().and_then(|c| {
                    let mut at = 0;
                    for &frame_end in &frame_ends {
                        c.compress_block_into(&raw[at..frame_end], &mut out)?;
                        at = frame_end;
                    }
                    Ok(())
                });
                self.raw = raw;
                self.frame_ends = frame_ends;
                compressed?;
                Ok(SpillBlockEvent::Block {
                    ordinal,
                    file_id,
                    key_kind,
                    is_last_in_file,
                    records_ingested_so_far,
                    bytes: out,
                })
            }
            SpillBlockEvent::Block {
                ordinal,
                file_id,
                key_kind,
                is_last_in_file,
                records_ingested_so_far,
                bytes,
            } => {
                let mut compressed = self.out_buf();
                self.compressor()?.compress_block_into(&bytes, &mut compressed)?;
                Ok(SpillBlockEvent::Block {
                    ordinal,
                    file_id,
                    key_kind,
                    is_last_in_file,
                    records_ingested_so_far,
                    bytes: compressed,
                })
            }
            // Passthrough variants carry no compressible payload.
            other @ (SpillBlockEvent::Residual { .. } | SpillBlockEvent::AllAnnounced { .. }) => {
                Ok(other)
            }
        }
    }
}

impl Clone for SpillBlockCompress {
    fn clone(&self) -> Self {
        // Fresh per-worker compressor + held slot; shared config copied.
        Self {
            codec: self.codec,
            compression: self.compression,
            compressor: None,
            raw: Vec::with_capacity(RAW_SCRATCH_BYTES),
            frame_ends: Vec::new(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            cap: self.cap.clone(),
            buf_pool: self.buf_pool.clone(),
        }
    }
}

impl Step for SpillBlockCompress {
    type Input = SpillBlockEvent;
    type Outputs = OrderedBytesSingle<SpillBlockEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SpillBlockCompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Moved out for the call: the permit borrows the cap while `compress_event`
        // needs `&mut self`. A move, so no reference-count traffic.
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn phase_cap(&self) -> Option<&PhaseCap> {
        self.cap.as_deref()
    }
}

#[cfg(test)]
mod tests;
