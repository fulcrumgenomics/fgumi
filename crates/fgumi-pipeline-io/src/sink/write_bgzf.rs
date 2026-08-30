//! `WriteBgzfFile` sink step. `Serial + Affinity::Writer` by default, or
//! `StepKind::Detached` (its own dedicated thread, off the pool) when built
//! via [`WriteBgzfFile::with_detached`] — used only on the standalone-sort
//! terminal (lever 2, legacy "N + 2"). Receives pre-compressed `BgzfBlock`s
//! from `BgzfCompress` and writes them directly to disk.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use fgumi_bam_io::{BaiBuilder, write_bai_index};
use fgumi_bgzf::{BGZF_EOF, BGZF_MAX_BLOCK_SIZE, InlineBgzfCompressor};
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::types::BgzfBlock;
use fgumi_pipeline_core::{
    header::HeaderHandle,
    step::{Affinity, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// `Serial + sticky` BAM sink (or `Detached` — see [`Self::with_detached`])
/// that consumes pre-compressed `BgzfBlock`s.
pub struct WriteBgzfFile {
    state: Mutex<Option<WriterState>>,
    name: &'static str,
    /// When `Some`, advertise `StepKind::Detached` so the framework drives this
    /// sink on its own dedicated driver thread (off the work-stealing pool) in the
    /// given [`DetachedGroup`], instead of `Serial + Affinity::Writer`. The caller
    /// (e.g. the sort chain) supplies the group via [`Self::with_detached`], so
    /// this generic sink carries no chain-specific grouping. `None` (the default)
    /// keeps the pool-scheduled writer that every other chain uses.
    detached_group: Option<DetachedGroup>,
}

struct WriterState {
    out: BufWriter<File>,
    pending_header: Option<PendingHeader>,
    /// Cumulative compressed bytes written so far, header included. Used to
    /// attribute each joined block's constituent BGZF blocks to their
    /// physical (compressed) file offsets for inline BAI indexing.
    coffset: u64,
    /// Present only when [`WriteBgzfFile::with_bai_index`] was called. `None`
    /// for every writer that doesn't build an inline index.
    bai: Option<BaiSink>,
}

/// Per-writer inline BAI indexing state: the [`BaiBuilder`] accumulating
/// resolved chunks, the running physical block-number watermark (continuing
/// across batches so block numbers stay globally unique and monotonic), the
/// reference count to build the final index with, and the sidecar path to
/// write it to.
struct BaiSink {
    bai: BaiBuilder,
    next_block_no: u64,
    num_refs: usize,
    sidecar_path: PathBuf,
}

/// Transform applied to the aligner's runtime-resolved header before it is
/// written. Lets a post-`Align` stage (sort-order rewrite, consensus header,
/// clip) re-apply its header change on top of the resolved header — which
/// carries the aligner's runtime `@PG`/`@RG`/`@CO` — instead of discarding that
/// provenance by writing a build-time header.
pub type ResolvedHeaderTransform = Box<dyn Fn(Header) -> io::Result<Header> + Send + Sync>;

struct PendingHeader {
    handle: HeaderHandle,
    compression_level: u32,
    transform: Option<ResolvedHeaderTransform>,
}

impl WriteBgzfFile {
    /// Open `path`, BGZF-compress and write the BAM header bytes, return
    /// the sink ready to receive `BgzfBlock`s.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from path open or header write.
    pub fn new<P: AsRef<Path>>(
        path: P,
        header: &Header,
        compression_level: u32,
    ) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        let mut out = BufWriter::with_capacity(256 * 1024, file);

        let mut header_bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_bytes, header)
            .map_err(|e| io::Error::other(format!("write_bam_header: {e}")))?;

        let mut hc = InlineBgzfCompressor::new(compression_level);
        hc.write_all(&header_bytes)?;
        hc.flush()?;
        // Compress to a `Vec` first (rather than `write_blocks_to(&mut out)`
        // directly) so the header's compressed length is known up front — the
        // starting point for `coffset`, the cumulative compressed-byte counter
        // the inline BAI join uses to attribute blocks to physical offsets.
        let mut header_blocks = Vec::new();
        hc.write_blocks_to(&mut header_blocks)?;
        let coffset = header_blocks.len() as u64;
        out.write_all(&header_blocks)?;

        Ok(Self {
            state: Mutex::new(Some(WriterState { out, pending_header: None, coffset, bai: None })),
            name: "WriteBgzfFile",
            detached_group: None,
        })
    }

    /// Run this sink on a dedicated `StepKind::Detached` driver thread (in the
    /// given [`DetachedGroup`]) instead of as a pool-scheduled `Serial +
    /// Affinity::Writer` step. Used ONLY on the standalone-sort terminal (lever 2):
    /// it frees a pool worker for the compression-bound work, matching the legacy
    /// sort's dedicated writer thread. The caller chooses the group so this generic
    /// sink stays chain-agnostic. The `try_run` body and the bytes it writes are
    /// unchanged — the driver pops blocks in the same (reorder-stage-ordered)
    /// sequence — so the output BAM is byte-identical to the pool-scheduled writer.
    /// Affinity is ignored for `Detached`.
    #[must_use]
    pub fn with_detached(mut self, group: DetachedGroup) -> Self {
        self.detached_group = Some(group);
        self
    }

    /// Attach an inline BAI indexer: as each [`BgzfBlock`] carrying a manifest
    /// (see [`crate::types::BamIndexManifest`]) is written, its records are
    /// joined into a fresh [`BaiBuilder`] using this writer's own cumulative
    /// compressed offset, and the resulting `.bai` is written to
    /// `sidecar_path` once the BGZF stream is drained and finished.
    ///
    /// Requires the eager-header constructor ([`Self::new`]) — a
    /// deferred-header writer ([`Self::new_with_handle`]) has not yet written
    /// its header (and so has no fixed `coffset` starting point) when this is
    /// called, and joining a manifest against a `coffset` that later shifts
    /// underneath it would produce a `.bai` with wrong virtual offsets.
    /// Deferred-header inline indexing is unsupported for now.
    ///
    /// # Precondition: the record stream MUST already be coordinate-sorted
    ///
    /// The inline indexer position-bins each record's [`AlignmentContext`]
    /// (from its manifest entry) into the [`BaiBuilder`] purely by the order
    /// records arrive at this sink — it does **not** read or verify `@HD SO`,
    /// and it does not re-derive sort order from the coordinates themselves
    /// (unlike the retired `IndexBamFinalizeHook`, which built the `.bai` via
    /// `noodles::bam::fs::index` reading back a `SO:coordinate` BAM). Feeding
    /// a non-coordinate-sorted stream through this sink produces a `.bai` that
    /// parses and looks valid but is semantically wrong (queries against it
    /// will silently miss or misattribute records).
    ///
    /// Coordinate-only is enforced upstream of this sink today, not here: the
    /// command layer rejects non-coordinate `--write-index` requests, and the
    /// chain builder (`fgumi_lib::pipeline::chains::builder`) only constructs
    /// a `BamWithIndex` sink spec for coordinate sorts. Any future caller of
    /// `with_bai_index` must preserve (or re-verify) that guarantee before
    /// wiring a new record stream into this sink.
    ///
    /// [`AlignmentContext`]: fgumi_bam_io::AlignmentContext
    ///
    /// # Panics
    ///
    /// Panics if `state.pending_header` is `Some` (i.e. this was built via
    /// [`Self::new_with_handle`]).
    #[must_use]
    pub fn with_bai_index(self, sidecar_path: PathBuf, num_refs: usize) -> Self {
        {
            let mut guard = self.state.lock();
            let state = guard.as_mut().expect("state present");
            assert!(
                state.pending_header.is_none(),
                "with_bai_index requires the eager-header constructor (deferred-header inline indexing is unsupported)"
            );
            state.bai =
                Some(BaiSink { bai: BaiBuilder::new(), next_block_no: 0, num_refs, sidecar_path });
        }
        self
    }

    /// Open `path` and return the sink with the BAM header write
    /// deferred until an upstream step resolves `handle`.
    ///
    /// `transform`, when `Some`, is applied to the resolved header before it is
    /// written — see [`ResolvedHeaderTransform`] for why this exists; pass `None`
    /// when the resolved header should be written unchanged.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from path open. Header-write errors are
    /// surfaced from `try_run` once the handle resolves.
    pub fn new_with_handle<P: AsRef<Path>>(
        path: P,
        handle: HeaderHandle,
        compression_level: u32,
        transform: Option<ResolvedHeaderTransform>,
    ) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        let out = BufWriter::with_capacity(256 * 1024, file);
        Ok(Self {
            state: Mutex::new(Some(WriterState {
                out,
                pending_header: Some(PendingHeader { handle, compression_level, transform }),
                coffset: 0,
                bai: None,
            })),
            name: "WriteBgzfFile",
            detached_group: None,
        })
    }

    fn try_write_pending_header(state: &mut WriterState) -> io::Result<bool> {
        let Some(pending) = state.pending_header.as_ref() else {
            return Ok(true);
        };
        let header_clone = match pending.handle.try_get() {
            None => return Ok(false),
            Some(Err(e)) => return Err(e),
            Some(Ok(h)) => h.clone(),
        };
        let level = pending.compression_level;
        // Re-apply the post-align header change (sort order / consensus / clip)
        // on top of the aligner's runtime-resolved header so its @PG/@RG/@CO
        // survive into the written header.
        let header_clone = match &pending.transform {
            Some(transform) => transform(header_clone)?,
            None => header_clone,
        };

        let mut header_bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_bytes, &header_clone)
            .map_err(|e| io::Error::other(format!("write_bam_header: {e}")))?;
        let mut hc = InlineBgzfCompressor::new(level);
        hc.write_all(&header_bytes)?;
        hc.flush()?;
        // As in `new`: measure the header's compressed length via a `Vec`
        // first so `state.coffset` stays correct for a future fused
        // deferred-header + inline-index path (today `with_bai_index` asserts
        // against a pending header, so this is currently dead weight — but
        // cheap and correct to keep in sync).
        let mut header_blocks = Vec::new();
        hc.write_blocks_to(&mut header_blocks)?;
        state.coffset += header_blocks.len() as u64;
        state.out.write_all(&header_blocks)?;

        state.pending_header = None;
        Ok(true)
    }
}

impl Step for WriteBgzfFile {
    type Input = BgzfBlock;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            // Detached (own driver thread) when a group was set via
            // `with_detached`; otherwise the default pool-scheduled Serial + sticky
            // writer. `sticky` is irrelevant for Detached (it never enters a
            // worker's worklist).
            kind: if self.detached_group.is_some() { StepKind::Detached } else { StepKind::Serial },
            sticky: true,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        // The caller-supplied group (the sort chain passes `SORT_IO_GROUP`);
        // `PerStep` fallback is unreachable for a Serial (non-detached) writer
        // since `detached_group()` is only consulted for `Detached` steps.
        self.detached_group.unwrap_or(DetachedGroup::PerStep)
    }

    fn affinity(&self) -> Affinity {
        // Ignored for `Detached` (no pool worker drives it); kept for the
        // default Serial path where it pins the writer to the last worker.
        Affinity::Writer
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut guard = self.state.lock();
        let Some(state) = guard.as_mut() else {
            return Ok(StepOutcome::Finished);
        };

        let header_ready = Self::try_write_pending_header(state)?;
        if header_ready && let Some(block) = ctx.input.pop() {
            // Join this block's manifest into the inline BAI builder, using
            // this writer's own cumulative compressed offset, BEFORE writing
            // the block's bytes — `state.coffset` at this point is exactly the
            // physical file offset the block's first physical BGZF block will
            // land at.
            if let Some(sink) = state.bai.as_mut()
                && let Some(manifest) = &block.index
            {
                let base = sink.next_block_no;
                let mut off = state.coffset; // read by value; no borrow overlap with `sink`
                for (j, &clen) in manifest.phys_comp_len.iter().enumerate() {
                    sink.bai.note_block(base + j as u64, off);
                    off += u64::from(clen);
                }
                debug_assert_eq!(off, state.coffset + block.bytes.len() as u64);
                for e in &manifest.records {
                    let block_no = base + (e.uoffset as usize / BGZF_MAX_BLOCK_SIZE) as u64;
                    let within = e.uoffset as usize % BGZF_MAX_BLOCK_SIZE;
                    sink.bai.record_with_context(block_no, within, e.len as usize, e.ctx);
                }
                sink.next_block_no += manifest.phys_comp_len.len() as u64;
                sink.bai.resolve().map_err(io::Error::other)?;
                debug_assert_eq!(
                    sink.bai.pending(),
                    0,
                    "resolve drains fully per batch on the arena layout"
                );
                // With the cache fully drained, `resolve`'s own front-gated
                // prune did not run (nothing left in `entry_cache` to gate on)
                // — prune explicitly so `block_positions` stays bounded.
                sink.bai.prune_below(sink.next_block_no);
            }
            state.out.write_all(&block.bytes)?;
            state.coffset += block.bytes.len() as u64;
            return Ok(StepOutcome::Progress);
        }

        if ctx.input.is_drained() {
            if !header_ready {
                return Err(io::Error::other(
                    "WriteBgzfFile: input drained before HeaderHandle was resolved",
                ));
            }
            state.out.write_all(&BGZF_EOF)?;
            state.out.flush()?;
            if let Some(sink) = state.bai.take() {
                let index = sink.bai.build(sink.num_refs).map_err(io::Error::other)?;
                write_bai_index(&sink.sidecar_path, &index).map_err(io::Error::other)?;
                log::info!("Wrote BAM index: {}", sink.sidecar_path.display());
            }
            let _ = guard.take();
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

impl Drop for WriteBgzfFile {
    /// Cleanup-only drop. The BGZF EOF marker is written **exclusively** by the
    /// drained-finish path in `try_run` (which then takes the state so this drop
    /// is a no-op for a normally-finished sink). If state is still present here
    /// the sink was dropped before that path ran — i.e. an aborted/partial
    /// stream — so we deliberately do **not** append `BGZF_EOF`: stamping the
    /// EOF marker onto a truncated BAM would make it look like a complete stream
    /// and hide the truncation from downstream readers. We only flush whatever
    /// bytes were already buffered so the on-disk file reflects what was written
    /// (and stays detectably truncated). A still-pending header means nothing
    /// valid was written, so leave the file empty.
    fn drop(&mut self) {
        let mut guard = self.state.lock();
        if let Some(mut state) = guard.take() {
            if state.pending_header.is_some() {
                return;
            }
            let _ = state.out.flush();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_header() -> Header {
        Header::default()
    }

    #[test]
    fn profile_advertises_serial_writer_sink() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = empty_header();
        let step = WriteBgzfFile::new(&path, &header, 1).unwrap();
        let profile = step.profile();
        assert_eq!(profile.name, "WriteBgzfFile");
        assert_eq!(profile.kind, StepKind::Serial);
        assert!(profile.sticky);
        assert_eq!(step.affinity(), Affinity::Writer);
        assert_eq!(profile.output_queues.len(), 0);
        assert_eq!(profile.branch_ordering.len(), 0);
    }

    /// L2.6: `with_detached()` flips the profile kind to `Detached` (its own
    /// thread, off the pool) while leaving everything else — name, the
    /// (now-irrelevant) sticky flag, the empty output edges — unchanged. The
    /// default constructors stay `Serial`, so only the standalone-sort terminal
    /// that opts in is affected.
    #[test]
    fn with_detached_advertises_detached_kind() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = empty_header();
        let step = WriteBgzfFile::new(&path, &header, 1)
            .unwrap()
            .with_detached(DetachedGroup::Shared("test-io"));
        let profile = step.profile();
        assert_eq!(profile.name, "WriteBgzfFile");
        assert_eq!(profile.kind, StepKind::Detached);
        assert_eq!(step.detached_group(), DetachedGroup::Shared("test-io"));
        assert_eq!(profile.output_queues.len(), 0);
    }

    #[test]
    fn header_only_round_trip() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = empty_header();
        let step = WriteBgzfFile::new(&path, &header, 1).unwrap();
        let mut guard = step.state.lock();
        let mut state = guard.take().expect("state present");
        state.out.write_all(&BGZF_EOF).unwrap();
        state.out.flush().unwrap();
        drop(guard);

        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.len() >= 28, "BGZF EOF + header should be at least 28 bytes");
        assert_eq!(&bytes[0..2], &[0x1f, 0x8b], "BGZF/gzip magic at start");
        let tail = &bytes[bytes.len() - 28..];
        assert_eq!(tail, &BGZF_EOF, "file ends with BGZF EOF marker");
    }

    /// Positive coverage for the drained-finish branch of `try_run`, driven
    /// through a real 2-step pipeline (block source → sink): the sink must write
    /// exactly one trailing `BGZF_EOF` and retire its state. The negative branch
    /// (`drop_before_finish_does_not_append_eof_marker`) and the by-hand EOF
    /// (`header_only_round_trip`) left this positive path — the one that decides
    /// whether an output BAM is a valid, EOF-terminated stream — otherwise unpinned.
    #[test]
    fn try_run_drained_finish_writes_exactly_one_eof() {
        use fgumi_pipeline_core::{
            Unpushed,
            builder::{Pipeline, PipelineConfig},
            held::HeldSlot,
            outputs::OrderedBytesSingle,
            queues::QueueSpec,
            reorder::BranchOrdering,
        };

        /// Exclusive source draining a `Vec<BgzfBlock>`, one block per `try_run`.
        struct BlockSource {
            blocks: Vec<BgzfBlock>,
            held: HeldSlot<Unpushed<BgzfBlock>>,
        }
        impl Step for BlockSource {
            type Input = ();
            type Outputs = OrderedBytesSingle<BgzfBlock>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "BlockSource",
                    kind: StepKind::Exclusive,
                    sticky: true,
                    output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                    branch_ordering: vec![BranchOrdering::ByItemOrdinal],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                if let Some(unpushed) = self.held.take()
                    && let Err(again) = ctx.outputs.retry(unpushed)
                {
                    self.held.put(again);
                    return Ok(StepOutcome::Progress);
                }
                let Some(block) = self.blocks.pop() else {
                    return Ok(StepOutcome::Finished);
                };
                if let Err(unpushed) = ctx.outputs.push(block) {
                    self.held.put(unpushed);
                }
                Ok(StepOutcome::Progress)
            }
        }

        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let sink = WriteBgzfFile::new(&path, &empty_header(), 1).unwrap();

        // Two payload blocks with dense ordinals (0, 1). `pop()` drains the tail
        // first, so ordinal 0 is emitted before ordinal 1 for the `ByItemOrdinal`
        // reorder stage. The bytes are arbitrary non-EOF markers.
        let blocks = vec![
            BgzfBlock { batch_serial: 1, bytes: vec![0xAB; 16], uncompressed_size: 0, index: None },
            BgzfBlock { batch_serial: 0, bytes: vec![0xCD; 16], uncompressed_size: 0, index: None },
        ];
        let source = BlockSource { blocks, held: HeldSlot::new() };

        let builder = Pipeline::builder();
        builder.chain(source).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("pipeline builds");
        pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() })
            .expect("pipeline runs to completion");

        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.len() >= 28, "header + payload + EOF present");
        assert_eq!(
            &bytes[bytes.len() - 28..],
            &BGZF_EOF,
            "the drained-finish path terminates the stream with a BGZF EOF marker"
        );
        // Exactly one trailing EOF: stripping it must not reveal a second, which
        // is what a stray `Drop`-side append would produce.
        let without_eof = &bytes[..bytes.len() - 28];
        assert!(
            without_eof.len() < 28 || without_eof[without_eof.len() - 28..] != BGZF_EOF,
            "only the drained-finish path may append EOF; Drop must not add a second"
        );
        // Both payload blocks reached disk.
        assert!(bytes.windows(16).any(|w| w == [0xAB; 16]), "first block written");
        assert!(bytes.windows(16).any(|w| w == [0xCD; 16]), "second block written");
    }

    /// L5 positive coverage for the inline-BAI join: drives a real 2-step
    /// `BlockSource → WriteBgzfFile` pipeline (mirroring
    /// `try_run_drained_finish_writes_exactly_one_eof`) with `with_bai_index`
    /// attached and two payload blocks carrying manifests. Block B's record
    /// straddles a physical block boundary (`uoffset % B + len > B`) to
    /// exercise `compute_end_vpos`'s forward walk across noted blocks. The
    /// synthetic compressed bytes are not valid BGZF — validity of the BGZF
    /// stream is not under test here, only that the manifest offsets join
    /// correctly into a parseable `.bai`.
    #[test]
    #[allow(clippy::too_many_lines)] // inline test Step + exhaustive offset assertions
    fn inline_index_writes_valid_bai_with_straddling_record() {
        use fgumi_bam_io::AlignmentContext;
        use fgumi_pipeline_core::{
            Unpushed,
            builder::{Pipeline, PipelineConfig},
            held::HeldSlot,
            outputs::OrderedBytesSingle,
            queues::QueueSpec,
            reorder::BranchOrdering,
        };
        use noodles::core::Position;

        use crate::types::{BamIndexManifest, RecordIndexEntry};

        /// Exclusive source draining a `Vec<BgzfBlock>`, one block per `try_run`.
        struct BlockSource {
            blocks: Vec<BgzfBlock>,
            held: HeldSlot<Unpushed<BgzfBlock>>,
        }
        impl Step for BlockSource {
            type Input = ();
            type Outputs = OrderedBytesSingle<BgzfBlock>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "BlockSource",
                    kind: StepKind::Exclusive,
                    sticky: true,
                    output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                    branch_ordering: vec![BranchOrdering::ByItemOrdinal],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                if let Some(unpushed) = self.held.take()
                    && let Err(again) = ctx.outputs.retry(unpushed)
                {
                    self.held.put(again);
                    return Ok(StepOutcome::Progress);
                }
                let Some(block) = self.blocks.pop() else {
                    return Ok(StepOutcome::Finished);
                };
                if let Err(unpushed) = ctx.outputs.push(block) {
                    self.held.put(unpushed);
                }
                Ok(StepOutcome::Progress)
            }
        }

        const B: usize = BGZF_MAX_BLOCK_SIZE;

        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let bai_path = fgumi_bam_io::bai_sidecar_path(&path);
        let sink = WriteBgzfFile::new(&path, &empty_header(), 1)
            .unwrap()
            .with_bai_index(bai_path.clone(), 1);

        let ctx0 = Some(AlignmentContext {
            ref_id: 0,
            start: Position::try_from(1).unwrap(),
            end: Position::try_from(11).unwrap(),
            is_mapped: true,
        });
        let ctx1 = Some(AlignmentContext {
            ref_id: 0,
            start: Position::try_from(200).unwrap(),
            end: Position::try_from(210).unwrap(),
            is_mapped: true,
        });

        // Block A (ordinal 0): one physical block, one record fully inside it.
        let block_a = BgzfBlock {
            batch_serial: 0,
            bytes: vec![0xAA; 50],
            uncompressed_size: 30,
            index: Some(Box::new(BamIndexManifest {
                phys_comp_len: vec![50],
                records: vec![RecordIndexEntry { uoffset: 0, len: 20, ctx: ctx0 }],
            })),
        };
        // Block B (ordinal 1): two physical blocks (block numbers 1 and 2,
        // continuing from block A's single block 0); one record straddles the
        // boundary between them.
        let straddle_uoffset = u32::try_from(B - 10).unwrap();
        let block_b = BgzfBlock {
            batch_serial: 1,
            bytes: vec![0xBB; 90],
            uncompressed_size: u32::try_from(B).unwrap() + 20,
            index: Some(Box::new(BamIndexManifest {
                phys_comp_len: vec![45, 45],
                records: vec![RecordIndexEntry { uoffset: straddle_uoffset, len: 20, ctx: ctx1 }],
            })),
        };

        // `pop()` drains the tail first, so ordinal 0 (block_a) is emitted
        // before ordinal 1 (block_b) — see the comment in
        // `try_run_drained_finish_writes_exactly_one_eof`.
        let blocks = vec![block_b, block_a];
        let source = BlockSource { blocks, held: HeldSlot::new() };

        let builder = Pipeline::builder();
        builder.chain(source).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("pipeline builds");
        pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() })
            .expect("pipeline runs to completion");

        let index_bytes = std::fs::read(&bai_path).expect("bai sidecar written");
        let index = noodles::bam::bai::io::Reader::new(&index_bytes[..])
            .read_index()
            .expect("bai sidecar parses");
        let refseqs = index.reference_sequences();
        assert_eq!(refseqs.len(), 1, "num_refs propagated");

        // The join's virtual offsets are hand-computable from the fixed manifests.
        // The only runtime unknown is the compressed header length; recover it
        // independently from the file layout — header (H bytes) + block_a (50) +
        // block_b (90) + BGZF EOF — so nothing below is read from the code under test.
        let bam_bytes = std::fs::read(&path).expect("bam written");
        let header_len = bam_bytes.len() as u64 - 50 - 90 - BGZF_EOF.len() as u64;

        let vpos = |coffset: u64, uoffset: u16| {
            noodles::bgzf::VirtualPosition::try_from((coffset, uoffset)).expect("valid vpos")
        };
        // Both records fall in the same bin (positions 1 and 200 are both inside
        // the first 16 kbp bin), so their per-record chunks coalesce into one
        // chunk spanning the earliest start to the latest end. Its endpoints pin
        // the two offsets the join must get right:
        //   start = record 0's start — block A's single physical block at
        //           coffset `header_len`, uncompressed offset 0;
        //   end   = record 1's end — the straddling record ends 10 bytes into
        //           block B's *second* physical block (coffset header_len + 95),
        //           which is the forward walk `compute_end_vpos` performs. A
        //           dropped straddling record would end at record 0's end
        //           (vpos(header_len, 20)) instead.
        let expected_start = vpos(header_len, 0);
        let expected_end = vpos(header_len + 95, 10);

        let chunks: Vec<_> = refseqs[0]
            .bins()
            .iter()
            .flat_map(|(_, bin)| bin.chunks().iter().map(|c| (c.start(), c.end())))
            .collect();
        assert_eq!(
            chunks,
            vec![(expected_start, expected_end)],
            "the joined bin must hold exactly one chunk spanning record 0's start to the \
             straddling record 1's end"
        );
    }

    /// No payload blocks reach the sink (only the header). Drained-finish must
    /// still build and write a `.bai` whose references are all present but
    /// empty — matching a coordinate-sorted BAM with no records.
    #[test]
    fn inline_index_empty_input_writes_valid_empty_bai() {
        use fgumi_pipeline_core::{
            Unpushed,
            builder::{Pipeline, PipelineConfig},
            held::HeldSlot,
            outputs::OrderedBytesSingle,
            queues::QueueSpec,
            reorder::BranchOrdering,
        };
        /// Exclusive source that immediately finishes with no blocks.
        struct BlockSource {
            blocks: Vec<BgzfBlock>,
            held: HeldSlot<Unpushed<BgzfBlock>>,
        }
        impl Step for BlockSource {
            type Input = ();
            type Outputs = OrderedBytesSingle<BgzfBlock>;
            fn profile(&self) -> StepProfile {
                StepProfile {
                    name: "BlockSource",
                    kind: StepKind::Exclusive,
                    sticky: true,
                    output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                    branch_ordering: vec![BranchOrdering::ByItemOrdinal],
                }
            }
            fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
                if let Some(unpushed) = self.held.take()
                    && let Err(again) = ctx.outputs.retry(unpushed)
                {
                    self.held.put(again);
                    return Ok(StepOutcome::Progress);
                }
                let Some(block) = self.blocks.pop() else {
                    return Ok(StepOutcome::Finished);
                };
                if let Err(unpushed) = ctx.outputs.push(block) {
                    self.held.put(unpushed);
                }
                Ok(StepOutcome::Progress)
            }
        }

        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let bai_path = fgumi_bam_io::bai_sidecar_path(&path);
        let sink = WriteBgzfFile::new(&path, &empty_header(), 1)
            .unwrap()
            .with_bai_index(bai_path.clone(), 3);

        let source = BlockSource { blocks: vec![], held: HeldSlot::new() };

        let builder = Pipeline::builder();
        builder.chain(source).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("pipeline builds");
        pipeline
            .run(PipelineConfig { threads: 2, ..Default::default() })
            .expect("pipeline runs to completion");

        let index_bytes = std::fs::read(&bai_path).expect("bai sidecar written for empty input");
        let index = noodles::bam::bai::io::Reader::new(&index_bytes[..])
            .read_index()
            .expect("bai sidecar parses");
        assert_eq!(
            index.reference_sequences().len(),
            3,
            "empty index still carries an entry per reference"
        );
        for rs in index.reference_sequences() {
            assert!(rs.bins().is_empty(), "no records means no bins");
        }
    }

    #[test]
    fn new_with_handle_defers_header_until_resolved() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path, handle.clone(), 1, None).unwrap();

        let bytes_before = std::fs::read(&path).unwrap();
        assert_eq!(bytes_before.len(), 0, "no bytes written until header resolves");

        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            assert!(state.pending_header.is_some(), "handle still pending");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(!wrote, "unresolved handle should yield without writing");
            assert!(state.pending_header.is_some(), "still pending after no-op probe");
        }
        let bytes_mid = std::fs::read(&path).unwrap();
        assert_eq!(bytes_mid.len(), 0, "still nothing on disk after no-op probe");

        handle.set(empty_header()).expect("first set");
        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(wrote, "resolved handle should write");
            assert!(state.pending_header.is_none(), "pending slot cleared");
            state.out.flush().unwrap();
        }
        let bytes_after = std::fs::read(&path).unwrap();
        assert!(bytes_after.len() >= 2, "header BGZF block emitted");
        assert_eq!(&bytes_after[0..2], &[0x1f, 0x8b], "BGZF/gzip magic at start");
    }

    #[test]
    fn resolved_header_transform_runs_on_the_resolved_header() {
        use std::sync::{Arc, Mutex as StdMutex};

        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();

        // Capture the header the transform is handed, to prove it is the
        // runtime-resolved (aligner) header — not a build-time header.
        let seen: Arc<StdMutex<Option<Header>>> = Arc::new(StdMutex::new(None));
        let seen_for_closure = Arc::clone(&seen);
        let transform: ResolvedHeaderTransform = Box::new(move |resolved: Header| {
            *seen_for_closure.lock().unwrap() = Some(resolved.clone());
            Ok(resolved)
        });

        let step =
            WriteBgzfFile::new_with_handle(&path, handle.clone(), 1, Some(transform)).unwrap();

        // Nothing is written (and the transform must not run) until the handle
        // resolves.
        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            assert!(!WriteBgzfFile::try_write_pending_header(state).unwrap());
        }
        assert!(seen.lock().unwrap().is_none(), "transform must not run before resolution");

        // Resolve with a distinctive header standing in for the aligner's
        // runtime-resolved header.
        let resolved = Header::builder().add_comment("ALIGNER-PROVENANCE").build();
        handle.set(resolved.clone()).expect("set");
        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            assert!(WriteBgzfFile::try_write_pending_header(state).unwrap(), "resolved -> written");
            state.out.flush().unwrap();
        }

        // The transform ran, and it received the resolved header — so a
        // post-align stage's transform is applied on top of the aligner header.
        assert_eq!(
            seen.lock().unwrap().as_ref(),
            Some(&resolved),
            "transform must receive the resolved header, preserving aligner provenance",
        );
        assert!(std::fs::read(&path).unwrap().len() >= 2, "header block was written");
    }

    #[test]
    fn new_with_handle_propagates_poison() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path, handle.clone(), 1, None).unwrap();

        handle.poison(io::Error::new(io::ErrorKind::BrokenPipe, "aligner died")).unwrap();
        let mut guard = step.state.lock();
        let state = guard.as_mut().expect("state present");
        let err = WriteBgzfFile::try_write_pending_header(state).expect_err("poison");
        assert_eq!(err.kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(err.to_string(), "aligner died");
    }

    #[test]
    fn drop_before_finish_does_not_append_eof_marker() {
        // A sink dropped before the drained-finish path in `try_run` (e.g. a
        // pipeline abort) must NOT append the BGZF EOF marker. Appending it
        // would stamp a "complete stream" signature onto a truncated BAM,
        // hiding the truncation from downstream readers. See the `Drop` doc.
        let path = tempfile::NamedTempFile::new().unwrap();
        let path_buf = path.path().to_path_buf();
        let header = empty_header();
        // `new` eagerly writes the header (pending_header is None), so the
        // only thing standing between this state and a valid EOF marker is
        // the `try_run` drained-finish path, which we never reach.
        let step = WriteBgzfFile::new(&path_buf, &header, 1).unwrap();
        drop(step);

        let bytes = std::fs::read(&path_buf).unwrap();
        assert!(bytes.len() >= 28, "header bytes should be on disk");
        let tail = &bytes[bytes.len() - 28..];
        assert_ne!(tail, &BGZF_EOF, "aborted output must not end with a valid BGZF EOF marker");
    }

    #[test]
    fn drop_with_unresolved_handle_leaves_empty_file() {
        let path = tempfile::NamedTempFile::new().unwrap();
        let path_buf = path.path().to_path_buf();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path_buf, handle, 1, None).unwrap();
        drop(step);

        let bytes = std::fs::read(&path_buf).unwrap();
        assert_eq!(bytes.len(), 0, "Drop with unresolved handle must skip EOF — see Drop doc");
    }
}
