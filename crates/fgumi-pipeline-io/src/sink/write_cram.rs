//! CRAM output sink (feature `cram`): three steps in place of
//! `BgzfCompress → WriteBgzfFile`.
//!
//! ```text
//! DecompressedBlock / OutputBlock (framed BAM records, in order)
//!   ↓ CramPack    Serial    records → htslib containers (no compression)
//!   ↓ CramEncode  Parallel  one container per job, on the pipeline's workers
//!   ↓ CramWrite   Serial    containers in order (+ CRAI), then EOF and close
//! ```
//!
//! This is htslib's own three-stage CRAM writer with fgumi's pool in place of
//! htslib's thread pool (see the `fgumi-htslib` crate): container encoding is
//! the expensive part and the only parallel one, exactly as for BGZF blocks.

use std::collections::VecDeque;
use std::io;
use std::marker::PhantomData;
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_htslib::{CramFinisher, CramWriteOptions, CramWriter, EncodeJob, EncodedContainer};
use fgumi_pipeline_core::{
    BranchOrdering, CounterSpec, HeapSize, HeldRetry, HeldSlot, Ordered, OrderedBytesSingle,
    PhaseCap, QueueSpec, StepCapRef, Unpushed, admit_input,
    header::HeaderHandle,
    step::{Affinity, DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::output_compress::{CompressInput, OutputBlock, OutputDepth};
use crate::sink::write_bgzf::{ResolvedHeaderTransform, resolve_deferred_header};

/// A container awaiting compression.
pub struct CramJob {
    job: EncodeJob,
}

impl HeapSize for CramJob {
    fn heap_size(&self) -> usize {
        usize::try_from(self.job.approx_bytes()).unwrap_or(usize::MAX)
    }
}
impl Ordered for CramJob {
    fn ordinal(&self) -> u64 {
        self.job.serial()
    }
}

/// A compressed container awaiting its in-order write.
pub struct CramEncoded {
    container: EncodedContainer,
    approx_bytes: u64,
}

impl HeapSize for CramEncoded {
    fn heap_size(&self) -> usize {
        // A heuristic, not a measurement: compressed CRAM is typically several
        // times smaller than the raw BAM records, so charge a quarter of the
        // job's approximate raw size. htslib's encoded container size is not
        // exposed by `fgumi_htslib::EncodedContainer`; this only sizes the
        // queue's backpressure, so a wrong ratio costs memory headroom, not
        // correctness.
        usize::try_from(self.approx_bytes / 4).unwrap_or(usize::MAX)
    }
}
impl Ordered for CramEncoded {
    fn ordinal(&self) -> u64 {
        self.container.serial()
    }
}

/// Where the output header comes from.
pub enum CramHeader {
    /// Known when the chain is built.
    Eager(Box<Header>),
    /// Resolved at run time (an aligner contributes `@PG`/`@RG`/`@CO`), then
    /// optionally transformed — the same contract as `WriteBgzfFile`.
    Deferred { handle: HeaderHandle, transform: Option<ResolvedHeaderTransform> },
}

/// Everything needed to create the file, consumed when it is opened.
pub struct CramSetup {
    /// Where the file is created: the destination, or a staged temporary
    /// renamed into place.
    pub path: PathBuf,
    /// The user's destination, named in errors (equal to `path` unless the
    /// output is staged).
    pub display: PathBuf,
    pub header: CramHeader,
    /// The reference FASTA, if one is configured.
    pub reference: Option<PathBuf>,
    /// Without a reference: store bases verbatim rather than let htslib look
    /// one up (see `fgumi_htslib::CramWriteOptions::no_ref`).
    pub no_ref: bool,
    /// The CRAM file-definition ID written in the file header (the
    /// destination's name when `path` is a staged temporary).
    pub file_id: Option<String>,
    /// CRAI path, when an index was requested (coordinate-sorted output only).
    pub crai: Option<PathBuf>,
}

/// The finisher handed from [`CramPack`] (which seals the writer) to
/// [`CramWrite`] (which closes it after the last container).
pub type FinisherSlot = Arc<Mutex<Option<CramFinisher>>>;

fn header_text(header: &Header) -> io::Result<String> {
    String::from_utf8(fgumi_bam_io::sam_header_text(header)?)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "SAM header is not UTF-8"))
}

// ---------------------------------------------------------------------------

/// Where the output file is in its life.
enum PackState {
    /// Not yet created: waiting for the header (a deferred header resolves at
    /// run time).
    Pending(CramSetup),
    /// Open: records are being packed.
    Open(CramWriter),
    /// Sealed: the last containers are queued and the finisher handed off.
    Sealed,
}

impl PackState {
    /// Open the file once the header is known, and return its writer.
    /// `Ok(None)`: the header is still pending (or the writer is sealed).
    fn ensure_open(&mut self) -> io::Result<Option<&mut CramWriter>> {
        if let PackState::Pending(setup) = self {
            let header = match &setup.header {
                CramHeader::Eager(h) => Some(header_text(h.as_ref())?),
                CramHeader::Deferred { handle, transform } => {
                    resolve_deferred_header(handle, transform.as_ref())?
                        .map(|h| header_text(&h))
                        .transpose()?
                }
            };
            let Some(text) = header else {
                return Ok(None);
            };
            *self = PackState::Open(create_writer(setup, &text)?);
        }
        Ok(match self {
            PackState::Open(writer) => Some(writer),
            PackState::Pending(_) | PackState::Sealed => None,
        })
    }
}

/// Create the CRAM file described by `setup`, with SAM header `text`.
fn create_writer(setup: &CramSetup, text: &str) -> io::Result<CramWriter> {
    let opts = CramWriteOptions {
        reference: setup.reference.as_deref(),
        no_ref: setup.no_ref,
        crai: setup.crai.as_deref(),
        file_id: setup.file_id.as_deref(),
        ..CramWriteOptions::default()
    };
    fgumi_pipeline_core::external_wait(|| CramWriter::create(&setup.path, text, &opts)).map_err(
        |e| {
            let what = if setup.path == setup.display {
                setup.display.display().to_string()
            } else {
                format!("{} (staged at {})", setup.display.display(), setup.path.display())
            };
            io::Error::new(e.kind(), format!("creating CRAM output {what}: {e}"))
        },
    )
}

/// `Serial` step: packs framed BAM records into CRAM containers.
pub struct CramPack<I: CompressInput = crate::types::DecompressedBlock> {
    state: PackState,
    pending: VecDeque<CramJob>,
    held: HeldSlot<Unpushed<CramJob>>,
    finisher: FinisherSlot,
    depth: Option<Arc<OutputDepth>>,
    output_byte_limit: u64,
    _input: PhantomData<fn(I)>,
}

/// Counter slots for [`CramPack`].
const RECORDS: usize = 0;

/// Input blocks [`CramPack`] packs per `try_run` before yielding, so the
/// serial step amortizes its dispatch over several blocks.
const MAX_BLOCKS_PER_RUN: usize = 8;

impl<I: CompressInput> CramPack<I> {
    #[must_use]
    pub fn new(setup: CramSetup, finisher: FinisherSlot, output_byte_limit: u64) -> Self {
        Self {
            state: PackState::Pending(setup),
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            finisher,
            depth: None,
            output_byte_limit,
            _input: PhantomData,
        }
    }

    /// Count pops in the producer's queue-depth gate (the standalone sort's
    /// range merge gates its pool claims on it, as for `BgzfCompress`).
    #[must_use]
    pub fn with_output_depth(mut self, depth: Option<Arc<OutputDepth>>) -> Self {
        self.depth = depth;
        self
    }

    /// Push the held job, then queued jobs, in order. `false`: the output
    /// refused one (held).
    fn flush(
        pending: &mut VecDeque<CramJob>,
        held: &mut HeldSlot<Unpushed<CramJob>>,
        ctx: &mut StepCtx<'_, Self>,
    ) -> bool {
        if matches!(ctx.outputs.retry_held(held), HeldRetry::StillHeld) {
            return false;
        }
        while let Some(job) = pending.pop_front() {
            if let Err(unpushed) = ctx.outputs.push(job) {
                held.put(unpushed);
                return false;
            }
        }
        true
    }

    /// Seal the open writer: queue its last containers and hand the finisher
    /// to [`CramWrite`].
    fn seal(&mut self) -> io::Result<()> {
        if let PackState::Open(writer) = std::mem::replace(&mut self.state, PackState::Sealed) {
            let (jobs, finisher) = writer.seal()?;
            self.pending.extend(jobs.into_iter().map(|job| CramJob { job }));
            *self.finisher.lock() = Some(finisher);
        }
        Ok(())
    }
}

impl<I: CompressInput> Step for CramPack<I> {
    type Input = I;
    type Outputs = OrderedBytesSingle<CramJob>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CramPack",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("records", "records")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !Self::flush(&mut self.pending, &mut self.held, ctx) {
            return Ok(StepOutcome::Contention);
        }
        if matches!(self.state, PackState::Sealed) {
            return Ok(StepOutcome::Finished);
        }
        if self.state.ensure_open()?.is_none() {
            if !ctx.input.is_drained() {
                return Ok(StepOutcome::NoProgress);
            }
            // Observing drain synchronizes with the producer's final header
            // publication, which may have followed our first header probe.
            if self.state.ensure_open()?.is_none() {
                return Err(io::Error::other(
                    "CramPack: input drained before the output header was resolved",
                ));
            }
        }
        let PackState::Open(writer) = &mut self.state else {
            unreachable!("ensure_open resolved the writer");
        };

        let mut packed = 0;
        let mut drained = false;
        while packed < MAX_BLOCKS_PER_RUN {
            let Some(item) = ctx.input.pop() else {
                drained = ctx.input.is_drained();
                break;
            };
            packed += 1;
            if let Some(depth) = &self.depth {
                depth.note_popped();
            }
            let block = match item.into_output_block() {
                OutputBlock::Raw(block) => block,
                OutputBlock::Compressed(_) => {
                    return Err(io::Error::other(
                        "CramPack: received a BGZF-compressed block; the producer must not \
                         compress output bound for CRAM",
                    ));
                }
            };
            let n = crate::types::walk_framed_records(&block.bytes, |body| writer.write_raw(body))?;
            ctx.counters.add(RECORDS, n);
            self.pending.extend(writer.drain_jobs().map(|job| CramJob { job }));
            if !Self::flush(&mut self.pending, &mut self.held, ctx) {
                // Backpressure: yield with the rest of the input unread.
                return Ok(StepOutcome::Progress);
            }
        }
        if drained {
            self.seal()?;
            Self::flush(&mut self.pending, &mut self.held, ctx);
            return Ok(StepOutcome::Progress);
        }
        Ok(if packed > 0 { StepOutcome::Progress } else { StepOutcome::NoProgress })
    }
}

// ---------------------------------------------------------------------------

/// `Parallel` step: compresses one CRAM container per call on a pool worker.
pub struct CramEncode {
    held: HeldSlot<Unpushed<CramEncoded>>,
    output_byte_limit: u64,
    cap: Option<Arc<PhaseCap>>,
}

impl Clone for CramEncode {
    fn clone(&self) -> Self {
        Self {
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            cap: self.cap.clone(),
        }
    }
}

impl CramEncode {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit, cap: None }
    }

    /// Share a phase admission cap (the standalone sort's phase-2 cap), as
    /// `BgzfCompress` does.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };
        let Some(CramJob { job }) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let approx_bytes = job.approx_bytes();
        let container = job.encode()?;
        if let Err(unpushed) = ctx.outputs.push(CramEncoded { container, approx_bytes }) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

impl Step for CramEncode {
    type Input = CramJob;
    type Outputs = OrderedBytesSingle<CramEncoded>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CramEncode",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn step_cap(&self) -> Option<StepCapRef<'_>> {
        self.cap.as_deref().map(StepCapRef::Fixed)
    }
}

// ---------------------------------------------------------------------------

/// `Serial` writer step (or `Detached`, see [`Self::with_detached`]): writes
/// compressed containers in order, then closes.
pub struct CramWrite {
    finisher: FinisherSlot,
    done: bool,
    /// When `Some`, the step is `Detached` in this group (its own driver thread,
    /// off the pool), as `WriteBgzfFile::with_detached` does for the BAM sink.
    detached_group: Option<DetachedGroup>,
}

const CONTAINERS_WRITTEN: usize = 0;

impl CramWrite {
    #[must_use]
    pub fn new(finisher: FinisherSlot) -> Self {
        Self { finisher, done: false, detached_group: None }
    }

    /// Run the writer on its own driver thread in `group` instead of as a
    /// pool-scheduled `Serial` step. The standalone sort uses this (as it does
    /// for `WriteBgzfFile`) so the writer's I/O does not take a pool turn. The
    /// containers are written in the same order either way.
    #[must_use]
    pub fn with_detached(mut self, group: DetachedGroup) -> Self {
        self.detached_group = Some(group);
        self
    }
}

impl Step for CramWrite {
    type Input = CramEncoded;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CramWrite",
            kind: if self.detached_group.is_some() { StepKind::Detached } else { StepKind::Serial },
            sticky: true,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        // Only consulted for a `Detached` step, i.e. when a group was set.
        self.detached_group.unwrap_or(DetachedGroup::PerStep)
    }

    fn affinity(&self) -> Affinity {
        // Ignored when `Detached`.
        Affinity::Writer
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("containers_written", "containers")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if self.done {
            return Ok(StepOutcome::Finished);
        }
        if let Some(encoded) = ctx.input.pop() {
            fgumi_pipeline_core::external_wait(|| encoded.container.write())?;
            ctx.counters.add(CONTAINERS_WRITTEN, 1);
            return Ok(StepOutcome::Progress);
        }
        if ctx.input.is_drained() {
            let Some(finisher) = self.finisher.lock().take() else {
                return Err(io::Error::other(
                    "CramWrite: input drained but the writer was never sealed",
                ));
            };
            fgumi_pipeline_core::external_wait(|| finisher.finish())?;
            self.done = true;
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Resolve the header exactly when the consumer observes drain, reproducing
    /// an empty producer finishing between the first header probe and drain.
    struct DrainHeader {
        header: Option<HeaderHandle>,
    }

    impl fgumi_pipeline_core::step::InputHandle<crate::types::DecompressedBlock> for DrainHeader {
        fn pop(&self) -> Option<crate::types::DecompressedBlock> {
            None
        }

        fn is_drained(&self) -> bool {
            if let Some(header) = &self.header
                && !header.is_set()
            {
                header.set(Header::default()).unwrap();
            }
            true
        }

        fn is_empty(&self) -> bool {
            true
        }
    }

    struct DrainProbe {
        pack: CramPack,
        input: DrainHeader,
    }

    impl Step for DrainProbe {
        type Input = ();
        type Outputs = OrderedBytesSingle<CramJob>;

        fn profile(&self) -> StepProfile {
            self.pack.profile()
        }

        fn counters(&self) -> &'static [CounterSpec] {
            self.pack.counters()
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.pack.try_run(&mut StepCtx {
                input: &self.input,
                outputs: ctx.outputs,
                counters: ctx.counters,
                pool: ctx.pool,
                signal: ctx.signal,
            })
        }
    }

    #[rstest::rstest]
    #[case::resolved_at_drain(true)]
    #[case::never_resolved(false)]
    fn deferred_empty_stream_rechecks_header_after_drain(#[case] resolve: bool) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.cram");
        let header = HeaderHandle::new();
        let finisher = FinisherSlot::default();
        let setup = CramSetup {
            path: path.clone(),
            display: path.clone(),
            header: CramHeader::Deferred { handle: header.clone(), transform: None },
            reference: None,
            no_ref: true,
            file_id: None,
            crai: None,
        };
        let mut step = DrainProbe {
            pack: CramPack::new(setup, finisher.clone(), 1 << 20),
            input: DrainHeader { header: resolve.then_some(header) },
        };
        let probe = fgumi_pipeline_core::testing::StepProbe::new(&step);
        let outcome = probe.try_run(&mut step);
        if resolve {
            assert_eq!(outcome.unwrap(), StepOutcome::Progress);
            finisher.lock().take().expect("empty writer sealed").finish().unwrap();
            let mut reader = fgumi_htslib::CramReader::open(&path, None).unwrap();
            assert!(reader.next_slice().unwrap().is_none());
        } else {
            let error = outcome.expect_err("an orphaned header remains an error");
            assert!(error.to_string().contains("before the output header was resolved"));
            assert!(!path.exists());
        }
    }

    /// `CramWrite` is a pool-scheduled `Serial` writer by default; the
    /// standalone sort's `with_detached` moves it to its own driver thread in
    /// the given group, as it does `WriteBgzfFile`.
    #[rstest::rstest]
    #[case::pool(None, StepKind::Serial, DetachedGroup::PerStep)]
    #[case::detached(
        Some(DetachedGroup::Shared("test-io")),
        StepKind::Detached,
        DetachedGroup::Shared("test-io")
    )]
    fn cram_write_is_serial_unless_detached(
        #[case] group: Option<DetachedGroup>,
        #[case] expected_kind: StepKind,
        #[case] expected_group: DetachedGroup,
    ) {
        let writer = CramWrite::new(FinisherSlot::default());
        let writer = group
            .map_or(writer, |group| CramWrite::new(FinisherSlot::default()).with_detached(group));
        let profile = writer.profile();
        assert_eq!(profile.name, "CramWrite");
        assert_eq!(profile.kind, expected_kind);
        assert_eq!(writer.detached_group(), expected_group);
    }
}
