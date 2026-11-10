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
//! With [`SpillBlockCompress::with_run_index`] (a zstd coordinate or
//! template-coordinate sort that can merge by range, on the main spill path
//! and the consolidation tail alike), each output `Block` of a run the phase-1
//! front still indexes (`index_stride > 0`) also carries a dense [`RunIndex`]
//! entry per frame — its compressed span, raw length, record count and first
//! key — which `SpillWrite` (or `WriteConsolidatedRun`) gathers into the run's
//! index.
//!
//! Each `Parallel` worker holds its own [`SpillBlockCompressor`], built lazily on
//! the first block (the zstd compressor's construction is fallible, so it is
//! surfaced through `try_run`'s `io::Result` rather than the infallible
//! `new_worker_copy`).

use std::io;
use std::sync::Arc;

use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use fgumi_sort::{
    RunIndex, SpillBlockCompressor, SpillCodec, SpillKeyKind, frame_meta, spill_key_width,
};

use crate::sort::protocol::SpillBlockEvent;
use fgumi_pipeline_core::{
    HeldRetry, PhaseCap, StepCap, StepCapRef, Unpushed, admit_input,
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
    /// The step's name (`SpillBlockCompress`, or the consolidation tail's
    /// `SpillBlockCompressConsolidation`): its stats row and phase.
    name: &'static str,
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
    cap: Option<StepCap>,
    /// Where compressed-block buffers come from (shared with `SpillWrite`).
    buf_pool: Option<Arc<SpillBufPool>>,
    /// Index each compressed block's frames ([`Self::with_run_index`]).
    run_index: bool,
    /// Per-worker scratch for a frame's first key while indexing.
    key_scratch: Vec<u8>,
}

impl SpillBlockCompress {
    /// Build a `SpillBlockCompress` for `codec` at `compression`. `output_byte_limit`
    /// byte-bounds the compressed-block output queue.
    #[must_use]
    pub fn new(codec: SpillCodec, compression: u32, output_byte_limit: u64) -> Self {
        Self {
            name: "SpillBlockCompress",
            codec,
            compression,
            compressor: None,
            raw: Vec::with_capacity(RAW_SCRATCH_BYTES),
            frame_ends: Vec::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            cap: None,
            buf_pool: None,
            run_index: false,
            key_scratch: Vec::new(),
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

    /// Name the step (its stats row and phase classification).
    #[must_use]
    pub fn named(mut self, name: &'static str) -> Self {
        self.name = name;
        self
    }

    /// Take compressed-block buffers from `pool` (shared with `SpillWrite`).
    #[must_use]
    pub fn with_buf_pool(mut self, pool: Arc<SpillBufPool>) -> Self {
        self.buf_pool = Some(pool);
        self
    }

    /// Index the frames of every block this step compresses (`frames` on each
    /// output `Block`) when `enabled`, the codec is zstd and the block's key
    /// kind is indexed ([`spill_key_width`]), for every block whose run the
    /// phase-1 front stamped a non-zero stride. Both the main spill path and
    /// the consolidation compressor enable it for a sort that can range.
    #[must_use]
    pub fn with_run_index(mut self, enabled: bool) -> Self {
        self.run_index = enabled;
        self
    }

    /// The key width of the index this step builds for a block keyed `kind`
    /// whose run the phase-1 front stamped `index_stride`, or `None` when
    /// that block is not indexed: indexing is on, the codec is zstd, the kind
    /// is indexed and the front still indexes the run (`index_stride > 0`).
    /// The index is always dense; `SpillWrite` applies the stride.
    fn index_width(&self, kind: SpillKeyKind, index_stride: u32) -> Option<u16> {
        if !self.run_index || self.codec != SpillCodec::Zstd || index_stride == 0 {
            return None;
        }
        spill_key_width(kind)
    }

    /// An empty index for a block keyed `kind` at `index_stride`, or `None`
    /// when this block's frames are not indexed ([`Self::index_width`]).
    fn new_index(&self, kind: SpillKeyKind, index_stride: u32) -> Option<RunIndex> {
        self.index_width(kind, index_stride).map(RunIndex::new)
    }

    /// Compress `raw`'s frames — each ending at the next of `frame_ends` — into
    /// one block buffer, indexing them when this block is indexed
    /// ([`Self::new_index`]).
    ///
    /// # Errors
    ///
    /// Propagates compressor-init, compression and frame-parse errors.
    fn compress_frames(
        &mut self,
        raw: &[u8],
        frame_ends: &[usize],
        key_kind: SpillKeyKind,
        index_stride: u32,
    ) -> io::Result<(Vec<u8>, Option<RunIndex>)> {
        let mut frames = self.new_index(key_kind, index_stride);
        let mut out = self.out_buf();
        let mut key_scratch = std::mem::take(&mut self.key_scratch);
        let compressed = self.compressor().and_then(|c| {
            let mut at = 0;
            for &frame_end in frame_ends {
                let frame = &raw[at..frame_end];
                compress_frame(c, frame, &mut out, key_kind, frames.as_mut(), &mut key_scratch)?;
                at = frame_end;
            }
            Ok(())
        });
        self.key_scratch = key_scratch;
        compressed.map(|()| (out, frames))
    }

    /// A buffer for one compressed block.
    fn out_buf(&self) -> Vec<u8> {
        self.buf_pool.as_ref().map(|p| p.take()).unwrap_or_default()
    }

    /// Admit through `cap` (`None` = uncapped): a fixed phase cap shared with
    /// the other pool steps of its phase, or a phased one that switches from
    /// phase 1 to phase 2 at `AllAnnounced`.
    #[must_use]
    pub fn with_step_cap(mut self, cap: impl Into<Option<StepCap>>) -> Self {
        self.cap = cap.into();
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
                index_stride,
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
                // Release the chunk (and with it, possibly, the run's arena)
                // before compressing.
                drop(chunk);
                let raw = std::mem::take(&mut self.raw);
                let frame_ends = std::mem::take(&mut self.frame_ends);
                let compressed = self.compress_frames(&raw, &frame_ends, key_kind, index_stride);
                self.raw = raw;
                self.frame_ends = frame_ends;
                let (out, frames) = compressed?;
                Ok(SpillBlockEvent::Block {
                    ordinal,
                    file_id,
                    key_kind,
                    is_last_in_file,
                    records_ingested_so_far,
                    bytes: out,
                    frames,
                    index_stride,
                })
            }
            SpillBlockEvent::Block {
                ordinal,
                file_id,
                key_kind,
                is_last_in_file,
                records_ingested_so_far,
                bytes,
                frames: _,
                index_stride,
            } => {
                let (compressed, frames) =
                    self.compress_frames(&bytes, &[bytes.len()], key_kind, index_stride)?;
                // The raw block is done with: recycle it for the producer that
                // draws its blocks from the same pool (`ConsolidateRuns`).
                if let Some(pool) = &self.buf_pool {
                    pool.give(bytes);
                }
                Ok(SpillBlockEvent::Block {
                    ordinal,
                    file_id,
                    key_kind,
                    is_last_in_file,
                    records_ingested_so_far,
                    bytes: compressed,
                    frames,
                    index_stride,
                })
            }
            // Passthrough variants carry no compressible payload.
            other @ (SpillBlockEvent::Residual { .. } | SpillBlockEvent::AllAnnounced { .. }) => {
                Ok(other)
            }
        }
    }
}

/// Compress one raw frame onto `out` and, when `index` is given, record it:
/// the bytes it added to `out` (its `u32` length prefix included), its raw
/// length, record count and first key. An empty frame writes nothing and is
/// not recorded.
///
/// # Errors
///
/// Propagates compression errors, and a frame whose records do not parse.
fn compress_frame(
    compressor: &mut SpillBlockCompressor,
    raw: &[u8],
    out: &mut Vec<u8>,
    key_kind: SpillKeyKind,
    index: Option<&mut RunIndex>,
    key_scratch: &mut Vec<u8>,
) -> io::Result<()> {
    let before = out.len();
    compressor.compress_block_into(raw, out)?;
    let Some(index) = index else { return Ok(()) };
    if raw.is_empty() {
        return Ok(());
    }
    let records = frame_meta(key_kind, raw, key_scratch)?;
    let comp_len = u32::try_from(out.len() - before)
        .map_err(|_| io::Error::other("SpillBlockCompress: a compressed frame exceeds 4 GiB"))?;
    let raw_len = u32::try_from(raw.len())
        .map_err(|_| io::Error::other("SpillBlockCompress: a raw frame exceeds 4 GiB"))?;
    index.push_frame(comp_len, raw_len, records, key_scratch);
    Ok(())
}

impl Clone for SpillBlockCompress {
    fn clone(&self) -> Self {
        // Fresh per-worker compressor + held slot; shared config copied.
        Self {
            name: self.name,
            codec: self.codec,
            compression: self.compression,
            compressor: None,
            raw: Vec::with_capacity(RAW_SCRATCH_BYTES),
            frame_ends: Vec::new(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            cap: self.cap.clone(),
            buf_pool: self.buf_pool.clone(),
            run_index: self.run_index,
            key_scratch: Vec::new(),
        }
    }
}

impl Step for SpillBlockCompress {
    type Input = SpillBlockEvent;
    type Outputs = OrderedBytesSingle<SpillBlockEvent>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
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
        let outcome = self.run_once(ctx, cap.as_ref().and_then(StepCap::active));
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn step_cap(&self) -> Option<StepCapRef<'_>> {
        self.cap.as_ref().map(StepCap::view)
    }
}

#[cfg(test)]
mod tests;
