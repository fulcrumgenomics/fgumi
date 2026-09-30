//! `SubprocessAlignStep` — the subprocess align backend.
//!
//! A `Serial` typed `Step` that wraps an aligner subprocess and pairs aligner
//! output with the original unmapped tags, emitting a `ZipperBatch` stream for
//! the shared `MergeAlignedStep` to merge.
//!
//! This is the former `AlignAndMergeStep` with its inline `merge_zipper_batch`
//! call removed: instead of merging on the dispatching worker, `try_run` now
//! pushes the raw `ZipperBatch` to an
//! `OrderedBytesSingle<ZipperBatch>` output. The reader still releases the
//! in-flight gate at token consumption, so the memory profile is unchanged. The
//! writer/reader threads, the byte-budget gate, header handling, and
//! finalize/Drop teardown are unchanged from the pre-split step.
//!
//! ## Shape
//!
//! ```text
//!  upstream (BamTemplateBatch)
//!     │
//!     ▼
//!  ┌──────────────── SubprocessAlignStep ────────────────┐
//!  │  Owns:                                               │
//!  │    • AlignerProcess (subprocess + stderr ring)       │
//!  │    • stdin-writer thread (FASTQ → aligner stdin)     │
//!  │    • stdout-reader thread (SAM/BAM parse; emits      │
//!  │      ZipperBatch)                                    │
//!  │    • Three channels (in_chan, token_chan, out_chan)  │
//!  │      and a shared error slot                         │
//!  │                                                      │
//!  │  Step::try_run (Serial; any worker can dispatch):    │
//!  │    upstream → in_chan → (writer) → aligner stdin     │
//!  │    aligner stdout → (reader) → out_chan → push       │
//!  │                              ZipperBatch downstream   │
//!  └──────────────────────────────────────────────────────┘
//! ```
//!
//! Why `Serial` instead of `Exclusive`: `Exclusive` pinned dispatch to
//! one worker, which meant any other step's blocking work on that
//! worker (e.g. a held-slot flush spin-retrying on a full downstream)
//! could starve this step and deadlock the pipeline.
//! `Serial` lets any free worker dispatch it; the framework mutex
//! preserves the "one dispatcher at a time" invariant that
//! `in_tx`/`out_rx`/`held_*` require.
//!
//! ## `BatchToken` / `ZipperBatch` protocol (index-based pairing)
//!
//! Pairing of unmapped templates with their alignments is **structural
//! by index, not queryname**. The writer thread keeps a per-batch
//! template count and ships a `BatchToken { unmapped, n_templates,
//! serial }` through `token_chan` to the reader thread. The reader
//! pops one token at a time, reads exactly `n_templates` templates
//! from aligner stdout, packages them into a `ZipperBatch { serial,
//! mapped, unmapped }`, and sends that on `out_chan`.
//! `try_run` pops `ZipperBatch`es and pushes them to the shared merge step.
//!
//! This relies on the aligner preserving input record order (bwa-mem
//! and bwa-mem3 do this when `-K` chunk size is set, which our presets
//! require). A per-template queryname-equality check — always on, not a
//! debug assertion — rejects out-of-order aligner output with a loud
//! error rather than silently mismerging tags onto the wrong reads.

use std::io::{self, BufReader, BufWriter, Read, Write};
use std::process::ChildStdout;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, Sender, SyncSender, TrySendError, channel, sync_channel};
use std::thread::{self, JoinHandle};

use anyhow::anyhow;
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::aligner::AlignerProcess;
use crate::commands::fastq::{FastqRecordBuffers, write_fastq_record};
use crate::pipeline::core::builder::PipelineBuilder;
use crate::pipeline::core::header::HeaderHandle;
use crate::pipeline::core::item::HeapSize;
use crate::pipeline::core::topology::{BranchIdx, StepIdx};
use crate::pipeline::core::{
    BranchOrdering, HeldSlot, OrderedBytesSingle, QueueSpec, Step, StepCtx, StepKind, StepOutcome,
    StepProfile, Unpushed,
};
use crate::pipeline::steps::align::merge::MergeAlignedStep;
use crate::pipeline::steps::align::{
    AlignBackend, AlignWired, AlignWiringCtx, ConsumerGoneGuard, InFlightGate, ZipperBatch,
    in_flight_budget_for_chunk_size, is_primary_for_alignment, merge_aligner_header,
    no_primary_records_message, split_pair_into_singles, validate_sq_consistency,
};
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::template::Template;

/// Number of aligner stderr lines retained for failure diagnostics.
const ALIGNER_STDERR_RING_SIZE: usize = 50;

/// Bound for `in_chan` (dispatcher → writer thread). Two batches in
/// flight match bwa's `kt_pipeline` `p_nt=2` double-buffer: one batch
/// being aligned, one queued.
///
/// Note: the writer→reader **token** channel is intentionally
/// unbounded (`std::sync::mpsc::channel()`), NOT bounded by this
/// constant. See the comment at the `token_chan` construction site
/// in [`SubprocessAlignStep::new`] for the deadlock rationale: bwa's
/// `-K` flag lets it buffer many input batches before emitting any
/// output, so the writer must be able to keep pushing tokens past
/// the reader without blocking. Backpressure on the writer-side flow is
/// enforced by the OS pipe buffer plus bwa's `-K`-bounded internal buffer.
const IN_CHAN_DEPTH: usize = 2;

/// Bound for `out_chan` (reader thread → dispatcher). Two batches lets
/// the reader stay one ahead while the dispatcher pushes the previous
/// one downstream; backpressure cascades back when downstream stalls.
const OUT_CHAN_DEPTH: usize = 2;

/// FASTQ-writer [`BufWriter`] capacity. 1 MiB: bigger than the OS pipe buffer
/// (~64 KiB) so we amortize `write` syscalls, but not so large that we delay the
/// aligner's first chunk perceptibly.
const FASTQ_WRITER_BUF_BYTES: usize = 1 << 20;

/// Stdout-side [`BufReader`] capacity for the aligner pipe. 64 KiB
/// pairs with one OS pipe buffer's worth of pending bytes.
const STDOUT_READER_BUF_BYTES: usize = 64 * 1024;

/// Stable error-message fragment for the "aligner exited before
/// emitting any bytes" branch in [`reader_loop_inner`]. Exposed so
/// tests can assert on a token (not the full prose) that's
/// guaranteed stable across rewordings of the user-facing message.
pub(crate) const ERR_ALIGNER_EXITED_BEFORE_OUTPUT: &str = "exited before emitting any output";

// ──────────────────────────────────────────────────────────────────────────
// Configuration
// ──────────────────────────────────────────────────────────────────────────

/// Immutable, Arc-friendly configuration for [`SubprocessAlignStep`].
///
/// One `Arc`-cheap clone of the inner is made per spawned thread.
#[derive(Clone)]
pub(crate) struct SubprocessConfig {
    /// Partial output header built at construction time
    /// (dict-derived `@SQ` + unmapped-derived `@HD`/`@CO`/`@RG`/`@PG` +
    /// fgumi's own `@PG`). The reader thread merges the aligner's emitted
    /// `@PG`/`@CO`/`@RG` lines into this and resolves `header_handle`.
    pub(crate) partial_output_header: Arc<Header>,

    /// One-shot handle the reader thread resolves with the merged
    /// header. The downstream `WriteBgzfFile` (constructed via
    /// `new_with_handle`) polls this.
    pub(crate) header_handle: HeaderHandle,

    /// Byte limit for the **downstream** output queue
    /// (`OrderedBytesSingle<ZipperBatch>` `ByteBounded`) this step pushes onto.
    pub(crate) output_byte_limit: u64,

    /// Byte budget for the **internal** in-flight unmapped reads (the
    /// writer→reader token backlog). The writer blocks once in-flight
    /// unmapped bytes reach this, bounding the otherwise-unbounded backlog
    /// so a fast-draining aligner can't accumulate the whole input's
    /// unmapped reads in RAM (issue #382). Derive it from the aligner's
    /// `-K` chunk size via [`in_flight_budget_for_chunk_size`].
    pub(crate) in_flight_unmapped_budget: u64,

    /// Whether a same-queryname output group with two unpaired primaries is
    /// accepted as bwa's mid-pair split (see [`AlignedGroup::MidPairSplit`]).
    /// Set only for the subprocess presets (`mem -p -K`); otherwise the group
    /// fails `Template::from_records`'s "multiple primaries" check, as it did
    /// before the split was recognized.
    pub(crate) accept_mid_pair_split: bool,
}

// ──────────────────────────────────────────────────────────────────────────
// Internal types
// ──────────────────────────────────────────────────────────────────────────

/// Metadata carried alongside an in-flight unmapped batch from the
/// stdin-writer to the stdout-reader.
///
/// The writer pushes one `BatchToken` after writing all FASTQ records
/// for an input batch. The reader pops one token, reads exactly
/// `n_templates` templates from the aligner's stdout, and emits one
/// `ZipperBatch { mapped, unmapped, serial }`.
struct BatchToken {
    unmapped: BamTemplateBatch,
    n_templates: usize,
    /// Serial number assigned by the writer in push order. Carried
    /// through into the `ZipperBatch` and ultimately the emitted
    /// `BamTemplateBatch::batch_serial` so downstream's
    /// `ByItemOrdinal` consumer sees a monotonic sequence.
    serial: u64,
}

/// Shared state visible to both threads and `try_run`. Used to
/// surface a thread-local IO error to the dispatcher loop.
///
/// The error slot is set-once: the first thread to fail wins; later
/// thread errors are dropped because they're likely cascading
/// consequences of the first (e.g., closed pipe → writer error
/// downstream of a crashed reader).
struct SharedState {
    error_slot: Mutex<Option<io::Error>>,
}

impl SharedState {
    fn new() -> Self {
        Self { error_slot: Mutex::new(None) }
    }

    /// Record an error if none has been recorded yet. First writer wins.
    fn record_error(&self, err: io::Error) {
        let mut guard = self.error_slot.lock();
        if guard.is_none() {
            *guard = Some(err);
        }
    }

    /// Take any recorded error.
    fn take_error(&self) -> Option<io::Error> {
        self.error_slot.lock().take()
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Step
// ──────────────────────────────────────────────────────────────────────────

/// `Serial` Step that wraps an aligner subprocess + two internal
/// threads and emits [`ZipperBatch`]es for the shared merge step.
///
/// See module doc for design rationale (including why `Serial`, not
/// `Exclusive`).
pub(crate) struct SubprocessAlignStep {
    cfg: SubprocessConfig,

    /// Owned aligner subprocess. Taken out in `finalize` (or `Drop` on the
    /// error path) for `wait()`.
    aligner: Option<AlignerProcess>,

    /// Sender side of `in_chan`. Dropped in `try_run`'s close-stdin phase to
    /// signal EOF to the writer thread (which closes aligner stdin
    /// in turn).
    in_tx: Option<SyncSender<BamTemplateBatch>>,

    /// Receiver side of `out_chan`. The reader thread sends raw
    /// `ZipperBatch`es here; `try_run` pushes them downstream.
    out_rx: Option<Receiver<ZipperBatch>>,

    writer_thread: Option<JoinHandle<io::Result<()>>>,
    reader_thread: Option<JoinHandle<io::Result<()>>>,

    shared: Arc<SharedState>,

    /// Held slot for a batch that was popped from `ctx.input` but
    /// couldn't be `try_send`'d to `in_chan` (writer thread is slow).
    held_in: HeldSlot<BamTemplateBatch>,

    /// Held slot for an `Unpushed<ZipperBatch>` we received from
    /// the reader but couldn't `push` to `ctx.outputs` (downstream
    /// backpressure).
    held_out: HeldSlot<Unpushed<ZipperBatch>>,

    name: &'static str,
}

impl SubprocessAlignStep {
    /// Spawn the aligner subprocess + the two I/O threads and return
    /// a ready-to-dispatch Step.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - the subprocess spawn fails,
    /// - the aligner doesn't provide stdin/stdout pipes (should be
    ///   impossible given `Stdio::piped()` is used by `AlignerProcess`),
    /// - either I/O thread spawn fails.
    pub(crate) fn new(cfg: SubprocessConfig, aligner_command: &str) -> io::Result<Self> {
        let mut aligner = AlignerProcess::spawn(aligner_command, ALIGNER_STDERR_RING_SIZE)
            .map_err(|e| io::Error::other(format!("SubprocessAlignStep::new: spawn: {e:#}")))?;

        let aligner_stdin = aligner
            .take_stdin()
            .ok_or_else(|| io::Error::other("aligner stdin pipe not available"))?;
        let aligner_stdout = aligner
            .take_stdout()
            .ok_or_else(|| io::Error::other("aligner stdout pipe not available"))?;

        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(IN_CHAN_DEPTH);
        // The `token_chan` itself stays unbounded: a bounded *channel*
        // would deadlock the writer once full (writer blocks on
        // `token_tx.send(...)`, never observes `in_rx` disconnection,
        // never drops its BufWriter, the aligner never sees stdin EOF,
        // never flushes). Instead the writer→reader unmapped backlog is
        // bounded by *bytes* via `InFlightGate`: the writer reserves
        // before feeding a batch and the reader releases when it consumes
        // the token. The budget is sized ≥ the aligner's `-K` buffering,
        // so a streaming aligner always emits before the writer blocks
        // (transient block, not a deadlock); bwa's `-K`-bounded internal
        // buffer + the OS pipe still provide the first line of throttling.
        // This bound is what keeps a fast-draining (non-`-K`-throttling)
        // aligner from accumulating the whole input's unmapped reads in
        // RAM — issue #382.
        let (token_tx, token_rx) = channel::<BatchToken>();
        let (out_tx, out_rx) = sync_channel::<ZipperBatch>(OUT_CHAN_DEPTH);

        let shared = Arc::new(SharedState::new());
        let gate = Arc::new(InFlightGate::new(cfg.in_flight_unmapped_budget));

        let writer_thread = {
            let shared = Arc::clone(&shared);
            let gate = Arc::clone(&gate);
            thread::Builder::new()
                .name("aam-fastq-writer".into())
                .spawn(move || writer_loop(in_rx, token_tx, aligner_stdin, shared, &gate))
                .map_err(|e| {
                    io::Error::other(format!(
                        "SubprocessAlignStep::new: failed to spawn writer thread: {e}"
                    ))
                })?
        };

        let reader_thread = {
            let cfg = cfg.clone();
            let shared = Arc::clone(&shared);
            let gate = Arc::clone(&gate);
            thread::Builder::new()
                .name("aam-sam-reader".into())
                .spawn(move || {
                    // Latch consumer-gone on EVERY reader exit — normal return,
                    // error, or a panic in `reader_loop` — via the RAII guard, so
                    // a writer blocked in `gate.acquire` always wakes and bails
                    // instead of hanging (and `Drop`'s `writer_thread.join()`
                    // cannot deadlock on a panicked reader).
                    let _consumer_gone = ConsumerGoneGuard(&gate);
                    reader_loop(token_rx, out_tx, aligner_stdout, cfg, shared, &gate)
                })
                .map_err(|e| {
                    io::Error::other(format!(
                        "SubprocessAlignStep::new: failed to spawn reader thread: {e}"
                    ))
                })?
        };

        Ok(Self {
            cfg,
            aligner: Some(aligner),
            in_tx: Some(in_tx),
            out_rx: Some(out_rx),
            writer_thread: Some(writer_thread),
            reader_thread: Some(reader_thread),
            shared,
            held_in: HeldSlot::new(),
            held_out: HeldSlot::new(),
            name: "SubprocessAlign",
        })
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Writer thread
// ──────────────────────────────────────────────────────────────────────────

/// Drain `in_rx` to the aligner's stdin as interleaved FASTQ, pushing
/// one `BatchToken` per batch to `token_tx` so the reader can pair
/// alignments back by count.
///
/// Lifecycle:
/// - When `in_rx` closes (sender dropped), flush + drop the
///   [`BufWriter`] so the aligner sees stdin EOF. `token_tx` drops
///   at the same time, signalling EOF to the reader after the last
///   token is consumed.
/// - On any IO error, record it in `shared.error_slot` and exit.
///   Channels' drop closes them; reader / `try_run` notice via
///   `Disconnected`.
///
/// The function intentionally owns its channel ends + `ChildStdin` so
/// they all get dropped at return (closing the aligner's stdin pipe).
#[allow(clippy::needless_pass_by_value)]
fn writer_loop(
    in_rx: Receiver<BamTemplateBatch>,
    token_tx: Sender<BatchToken>,
    aligner_stdin: std::process::ChildStdin,
    shared: Arc<SharedState>,
    gate: &InFlightGate,
) -> io::Result<()> {
    let mut writer = BufWriter::with_capacity(FASTQ_WRITER_BUF_BYTES, aligner_stdin);
    let mut fastq_buffers = FastqRecordBuffers::with_capacity(512);
    let mut next_serial: u64 = 0;

    while let Ok(batch) = in_rx.recv() {
        // Reserve this batch's unmapped bytes against the in-flight budget
        // before feeding it to the aligner. Blocks (transiently) if the
        // backlog is at budget; returns `false` only if the reader has
        // exited, in which case we stop (the channel send below would also
        // fail). Bounds the otherwise-unbounded token backlog (issue #382).
        if !gate.acquire(batch.heap_size() as u64) {
            return Ok(());
        }
        // Count only templates that contribute at least one record
        // to the FASTQ stream. A template whose every record is
        // filtered by `FASTQ_WRITER_EXCLUDE_FLAGS` produces zero
        // aligner output; if we counted it, the reader would expect
        // one more mapped template than the aligner emits and
        // surface a misleading "aligner emitted fewer alignments"
        // error. Expected input is an unmapped BAM (from
        // `fgumi extract`), which has no secondaries — a
        // fully-filtered template here means the user passed a
        // re-aligned BAM by mistake. Error out loudly with the
        // queryname so the misconfiguration is obvious.
        let mut n_templates = 0usize;
        for template in batch.templates() {
            let mut wrote_any = false;
            for record in &template.records {
                let flags = record.flags();
                // Align primary reads only (skip SECONDARY/SUPPLEMENTARY).
                if !is_primary_for_alignment(flags) {
                    continue;
                }
                if let Err(e) = write_fastq_record(
                    &mut writer,
                    record,
                    flags,
                    /* no_suffix */ true,
                    &mut fastq_buffers,
                    /* umi_header */ None,
                ) {
                    // Record + return the same full error message so
                    // `finalize`'s join path and the shared
                    // error_slot fallback path both see the same
                    // text — no information loss either way. Same
                    // pattern as the no-primary path below.
                    let msg = format!(
                        "align-and-merge writer: write_fastq_record for template '{name}': {e:#}",
                        name = String::from_utf8_lossy(template.name()),
                    );
                    shared.record_error(io::Error::other(msg.clone()));
                    return Err(io::Error::other(msg));
                }
                wrote_any = true;
            }
            if !wrote_any {
                let msg = format!(
                    "align-and-merge writer: {}",
                    no_primary_records_message(template.name())
                );
                shared.record_error(io::Error::other(msg.clone()));
                return Err(io::Error::other(msg));
            }
            // KNOWN DIVERGENCE — resolve in the align-and-merge WIRING PR.
            // This guard only rejects a template where EVERY record is filtered
            // (no primary at all). A *partially* filtered template — e.g. a
            // re-aligned BAM by mistake where one mate has a primary and the
            // other is only SECONDARY/SUPPLEMENTARY — writes fewer FASTQ records
            // than the unmapped template physically carries, and nothing
            // downstream re-checks the per-template record count, so the reader
            // pairs a short mapped template against the full unmapped one and
            // merge_raw can mismerge. The queryname hard-check in the reader
            // (see the pairing guard) already converts the common reorder case
            // into a loud error, so this is a narrow residual edge — and it is
            // unreachable on real input (`fgumi extract` output carries no
            // secondaries). The correct fix compares *primaries written* against
            // the unmapped template's expected mate count
            // (`r1.is_some() + r2.is_some()`), NOT `records.len()` (which counts
            // secondaries and would false-positive on the exact re-aligned-BAM
            // shape the SECONDARY/SUPPLEMENTARY filter is designed to tolerate).
            // Deferred so the wiring PR's byte-parity tests gate it.
            n_templates += 1;
        }

        let token = BatchToken { unmapped: batch, n_templates, serial: next_serial };
        next_serial = next_serial.wrapping_add(1);

        // Send is bounded; block until the reader (or backpressure)
        // makes room. If the reader has died, the channel disconnect
        // surfaces here and we exit.
        if token_tx.send(token).is_err() {
            // Reader thread is gone — we cannot make progress.
            // Don't record an error: the reader's own error (if any)
            // is the root cause; we just shut down quietly.
            return Ok(());
        }
    }

    // in_rx closed → upstream is done. Flush + drop the writer to
    // close the aligner's stdin → aligner flushes its output.
    if let Err(e) = writer.flush() {
        // Record and return the SAME detailed message so `finalize`'s join path
        // and the error_slot fallback both see the actionable text (matching the
        // write_fastq_record / no-primary paths above).
        let msg = format!("align-and-merge writer: flush: {e}");
        shared.record_error(io::Error::other(msg.clone()));
        return Err(io::Error::other(msg));
    }
    drop(writer);
    Ok(())
}

// ──────────────────────────────────────────────────────────────────────────
// Reader thread
// ──────────────────────────────────────────────────────────────────────────

/// Parse the aligner's stdout, pair each emitted template with its
/// corresponding unmapped via the token channel, and emit
/// `ZipperBatch` to `out_tx`. Merging happens on the shared
/// `MergeAlignedStep`; the reader stays I/O-bound so bwa's stdout never
/// stalls behind CPU-heavy merge work.
///
/// First action is to read the aligner's SAM/BAM header, merge with
/// the partial header, and resolve `cfg.header_handle`. Once that
/// resolves, the downstream writer can start.
///
/// Takes its channel ends + `ChildStdout` + config by value so they
/// drop at return; the wrapping helper translates the `Result` into
/// an error-slot + handle-poison action.
#[allow(clippy::needless_pass_by_value)]
fn reader_loop(
    token_rx: Receiver<BatchToken>,
    out_tx: SyncSender<ZipperBatch>,
    aligner_stdout: ChildStdout,
    cfg: SubprocessConfig,
    shared: Arc<SharedState>,
    gate: &InFlightGate,
) -> io::Result<()> {
    let result = reader_loop_inner(token_rx, &out_tx, aligner_stdout, &cfg, gate);
    match result {
        Ok(()) => Ok(()),
        Err(e) => {
            // Poison the HeaderHandle so the downstream writer's
            // `try_get` resolves to Err instead of looping forever.
            // It's harmless if the handle was already set (e.g.,
            // error happened mid-record after header was emitted) —
            // `poison` returns AlreadySetError, ignored.
            let _ = cfg
                .header_handle
                .poison(io::Error::new(e.kind(), format!("align-and-merge reader: {e}")));
            shared.record_error(io::Error::new(e.kind(), format!("align-and-merge reader: {e}")));
            Err(e)
        }
    }
}

// Long but linear: the reader thread's body is a five-step
// procedure (peek format → set up BAM reader → validate @SQ →
// resolve HeaderHandle → process tokens) that fits naturally in
// one function. Splitting it would require threading the BAM
// reader through helper signatures whose generics are already
// noisy enough (see `AlignerBamReader`).
#[allow(clippy::needless_pass_by_value, clippy::too_many_lines)]
fn reader_loop_inner(
    token_rx: Receiver<BatchToken>,
    out_tx: &SyncSender<ZipperBatch>,
    aligner_stdout: ChildStdout,
    cfg: &SubprocessConfig,
    gate: &InFlightGate,
) -> io::Result<()> {
    // 1. Peek 4 bytes to detect BGZF (BAM) vs SAM text. Returns the
    //    peeked bytes prepended to a chained Read so we don't lose
    //    the prefix.
    let (is_bgzf, peek_filled, stdout) = peek_aligner_format(aligner_stdout)?;

    if peek_filled == 0 {
        // Aligner exited before emitting any output — usually an
        // aligner crash or misconfiguration (wrong binary path,
        // missing index files, segfault before reading FASTQ).
        // Distinguished from "non-BGZF content" so the message
        // doesn't misleadingly point at SAM/BAM format conversion.
        let marker = ERR_ALIGNER_EXITED_BEFORE_OUTPUT;
        return Err(io::Error::other(format!(
            "SubprocessAlignStep: aligner {marker}. Likely causes: subprocess \
             startup error (missing binary, missing index files), aligner \
             argument error, or the aligner segfaulted before reading FASTQ. \
             Check the aligner's stderr above for the specific failure."
        )));
    }

    // 2. Set up the appropriate record source. BAM = BGZF →
    //    noodles BAM reader → unwrap back into raw bytes reader.
    //    SAM = text-mode noodles SAM reader. Both paths parse the
    //    aligner-emitted header first; we feed it through the same
    //    `@SQ` validation + merge + `HeaderHandle` resolution before
    //    constructing the `TemplateStream` variant.
    //
    //    The noodles BAM reader does *not* read past the header
    //    into the body, so `bam_reader.into_inner()` returns a BGZF
    //    reader positioned at the first record. Same property
    //    holds for the SAM reader.
    let (mut template_stream, aligner_header) = if is_bgzf {
        let buffered = BufReader::with_capacity(STDOUT_READER_BUF_BYTES, stdout);
        let bgzf = noodles::bgzf::io::Reader::new(buffered);
        let mut bam_reader = noodles::bam::io::Reader::from(bgzf);
        let aligner_header = bam_reader
            .read_header()
            .map_err(|e| io::Error::other(format!("aligner BAM header: {e}")))?;
        let bgzf = bam_reader.into_inner();
        let stream = TemplateStream::Bam(BamTemplateStream::new(
            fgumi_raw_bam::RawBamReader::new(bgzf),
            cfg.accept_mid_pair_split,
        ));
        (stream, aligner_header)
    } else {
        let buffered = BufReader::with_capacity(STDOUT_READER_BUF_BYTES, stdout);
        let mut sam_reader = noodles::sam::io::Reader::new(buffered);
        let aligner_header = sam_reader
            .read_header()
            .map_err(|e| io::Error::other(format!("aligner SAM header: {e}")))?;
        let header_arc = Arc::new(aligner_header.clone());
        let stream = TemplateStream::Sam(SamTemplateStream::new(
            sam_reader,
            header_arc,
            cfg.accept_mid_pair_split,
        ));
        (stream, aligner_header)
    };

    // 3. Verify the aligner's `@SQ` table matches the partial
    //    header's (which came from the reference dict). If they
    //    don't, the aligner was built against a different FASTA
    //    than the dict and the resulting BAM's `tid` integers
    //    would silently point at the wrong contig. Catch it here
    //    rather than emit a corrupt BAM.
    validate_sq_consistency(&cfg.partial_output_header, &aligner_header)?;

    // 4. Merge aligner header into the partial output header, then
    //    resolve the HeaderHandle so the downstream writer can start.
    //    `merge_aligner_header` is intentionally permissive: aligner
    //    `@PG` / `@RG` / `@CO` lines are added; `@SQ` from the
    //    aligner is discarded because the partial header already has
    //    the dict-derived canonical reference list.
    let merged = merge_aligner_header(&cfg.partial_output_header, &aligner_header);
    // The handle should be unresolved at this point in production —
    // this step is the sole producer. A pre-set handle indicates either a
    // test (where `from_header` was used) or a wiring bug. We log
    // the latter case at warn level so it surfaces in CI / prod
    // pipelines while letting tests proceed unchanged. The merged
    // header is dropped in either case; the existing value wins.
    if cfg.header_handle.set(merged).is_err() {
        log::warn!(
            "SubprocessAlignStep: HeaderHandle was already resolved before the aligner \
             emitted its header — aligner @PG/@RG/@CO contributions will not appear \
             in the output. This is expected in tests using HeaderHandle::from_header \
             but indicates a wiring bug in production."
        );
    }

    // 5. Process tokens. One token = one upstream batch. We read
    //    exactly token.n_templates templates from the aligner output
    //    for each token, via the `TemplateStream` iterator (which
    //    encapsulates the scratch + peeked state).
    while let Ok(token) = token_rx.recv() {
        // The token's unmapped bytes have left the (bounded) token backlog;
        // release them against the in-flight budget so a writer blocked in
        // `gate.acquire` can resume. The bytes still live briefly in the
        // `ZipperBatch` flowing through the depth-bounded `out_chan` +
        // `held_out`, which is separately bounded — so releasing
        // here is what keeps the *unbounded* token backlog near budget.
        gate.release(token.unmapped.heap_size() as u64);

        if token.n_templates == 0 {
            // Empty batch from upstream — write nothing FASTQ-side,
            // so the aligner emits nothing for this token. Emit an
            // empty `ZipperBatch` to keep ordinals contiguous so
            // downstream `ByItemOrdinal` consumers don't stall.
            // `GroupByQueryname` never emits empty batches in
            // practice, but the framework's batch contract doesn't
            // forbid them and emitting empty here is robust to
            // future upstream changes.
            let empty_unmapped = BamTemplateBatch::new(token.serial, Vec::new());
            let zb =
                ZipperBatch { serial: token.serial, mapped: Vec::new(), unmapped: empty_unmapped };
            if out_tx.send(zb).is_err() {
                return Ok(()); // downstream gone
            }
            continue;
        }

        let mut mapped_templates: Vec<Template> = Vec::with_capacity(token.n_templates);
        // Unmapped-template indices bwa split mid-pair; their unmapped halves are
        // split to match after the loop (rare: needs mixed single/paired input).
        let mut split_indices: Vec<usize> = Vec::new();

        for i in 0..token.n_templates {
            let mapped = template_stream.next_template()?;

            let mapped = mapped.ok_or_else(|| {
                io::Error::other(format!(
                    "align-and-merge reader: alignerstdout EOF after {emitted} of {expected} \
                     templates in batch {serial} — aligner emitted fewer alignments \
                     than input reads",
                    emitted = i,
                    expected = token.n_templates,
                    serial = token.serial,
                ))
            })?;

            // Verify the i-th aligner-emitted template really is the read we
            // paired it with by position. Order preservation is an *external*
            // contract (bwa-mem/bwa-mem3 with `-K` honor it), but
            // `aligner_command` is a free-form user command, so a reordering
            // aligner could keep the per-batch count while shuffling templates —
            // silently merging one read's tags/UMI onto another's alignment.
            // The compare is one byte-slice comparison per template (both names
            // are already in hand), negligible against the per-record merge_raw
            // it guards, so it stays on in release rather than as a debug_assert.
            if token.unmapped.templates()[i].name() != mapped.name() {
                return Err(io::Error::other(format!(
                    "align-and-merge reader: queryname mismatch — unmapped[{i}]='{u}' but \
                     mapped='{m}'; the aligner emitted templates out of input order (does the \
                     aligner command preserve input order, e.g. bwa `-K`?)",
                    u = String::from_utf8_lossy(token.unmapped.templates()[i].name()),
                    m = String::from_utf8_lossy(mapped.name()),
                )));
            }

            match mapped {
                AlignedGroup::Template(t) => mapped_templates.push(t),
                AlignedGroup::MidPairSplit(first, second) => {
                    mapped_templates.push(first);
                    mapped_templates.push(second);
                    split_indices.push(i);
                }
            }
        }

        let unmapped = if split_indices.is_empty() {
            token.unmapped
        } else {
            split_unmapped_pairs(token.unmapped, &split_indices)?
        };
        let zb = ZipperBatch { serial: token.serial, mapped: mapped_templates, unmapped };
        if out_tx.send(zb).is_err() {
            // Downstream gone — bail. The dispatcher will observe
            // out_rx disconnected on its next try_run.
            return Ok(());
        }
    }

    // token_rx closed → writer is done. If TemplateStream has a
    // peeked record still in flight, the aligner emitted at least
    // one more template than the sum of n_templates across all
    // tokens — fatal mismatch (input reads != output reads).
    if template_stream.has_peeked() {
        return Err(io::Error::other(
            "align-and-merge reader: aligneremitted at least one more template than the sum of \
             input batch sizes (residual record carried after final token)",
        ));
    }

    // Probe for one more record past the final token. If the
    // aligner emitted multiple extra templates this only detects
    // the first; the count isn't useful for the error message, so
    // we report "at least one" rather than a plural that would be
    // wrong for exactly-one cases.
    if template_stream.probe_trailing().map_err(|e| {
        io::Error::other(format!("align-and-merge reader: probing trailing aligner output: {e}"))
    })? {
        return Err(io::Error::other(
            "align-and-merge reader: aligneremitted at least one more record than the sum of \
             input batch sizes (extra record at EOF)",
        ));
    }

    Ok(())
}

/// Rebuild `unmapped` with each template at `split_indices` (ascending) split
/// into its two single-read halves, matching the aligner's mid-pair split so the
/// mapped and unmapped halves still zip one-to-one.
fn split_unmapped_pairs(
    unmapped: BamTemplateBatch,
    split_indices: &[usize],
) -> io::Result<BamTemplateBatch> {
    let (serial, templates) = unmapped.into_parts();
    let mut out = Vec::with_capacity(templates.len() + split_indices.len());
    let mut splits = split_indices.iter().peekable();
    for (i, template) in templates.into_iter().enumerate() {
        if splits.next_if_eq(&&i).is_some() {
            let (first, second) = split_pair_into_singles(template)?;
            out.push(first);
            out.push(second);
        } else {
            out.push(template);
        }
    }
    Ok(BamTemplateBatch::new(serial, out))
}

/// Peek the leading bytes of `stdout` to classify the aligner's output
/// format. Returns `(is_bgzf, peek_filled, chained_reader)` where:
/// - `is_bgzf` is true iff the peeked prefix is a valid BGZF block
///   header per the shared [`fgumi_bgzf::is_bgzf_header`] classifier,
///   which validates the full header (magic + FEXTRA + the `BC`
///   subfield), not just the 4 magic bytes — so a plain-gzip SAM stream
///   is correctly NOT treated as BAM.
/// - `peek_filled` is the number of bytes actually peeked
///   (`0..=BGZF_HEADER_SIZE`). `0` means the aligner exited before
///   producing any output, which the caller distinguishes from a
///   "wrong format" condition for a clearer error message.
/// - `chained_reader` re-emits the peeked bytes followed by the rest
///   of the original stream — callers must drain the returned reader,
///   not the original.
///
/// The return type is `Box<dyn Read + Send>` so the same boxed reader
/// can feed both the BGZF (BAM) path and the SAM-text path without
/// changing every downstream type signature.
fn peek_aligner_format(mut stdout: ChildStdout) -> io::Result<(bool, usize, Box<dyn Read + Send>)> {
    // Peek a full BGZF header's worth so the shared validator has enough bytes
    // for a real verdict (a genuine BGZF header is exactly `BGZF_HEADER_SIZE`
    // bytes). A shorter (or empty) prefix simply fails validation and falls to
    // the SAM-text path; `peek_filled == 0` is handled specially upstream as
    // "aligner exited before emitting any output".
    let mut peek_buf = [0u8; fgumi_bgzf::BGZF_HEADER_SIZE];
    let mut peek_filled = 0usize;
    while peek_filled < peek_buf.len() {
        match stdout.read(&mut peek_buf[peek_filled..]) {
            Ok(0) => break, // aligner exited before filling the peek buffer
            Ok(n) => peek_filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    let is_bgzf = fgumi_bgzf::is_bgzf_header(&peek_buf[..peek_filled]);
    let prefix = std::io::Cursor::new(peek_buf[..peek_filled].to_vec());
    let chained: Box<dyn Read + Send> = Box::new(prefix.chain(stdout));
    Ok((is_bgzf, peek_filled, chained))
}

/// Concrete type of the BAM reader we use over the aligner's piped
/// stdout. Defined to keep [`TemplateStream`]'s type signature
/// manageable and to make the (single-threaded, pipe-friendly)
/// BGZF stack explicit at this call site.
type AlignerBamReader =
    fgumi_raw_bam::RawBamReader<noodles::bgzf::io::Reader<BufReader<Box<dyn Read + Send>>>>;

/// Concrete type of the SAM reader we use over the aligner's piped
/// stdout when the peeked magic bytes show plain text (not BGZF).
/// Same dyn-erasure boundary as `AlignerBamReader`.
type AlignerSamReader = noodles::sam::io::Reader<BufReader<Box<dyn Read + Send>>>;

/// Reads templates one at a time from the aligner's emitted output
/// (BAM-on-pipe or SAM-text-on-pipe), groups records by queryname,
/// and yields fully-assembled [`Template`]s.
///
/// Two variants for the two stdout formats the aligner can produce.
/// They share the same `scratch` + `peeked` + `name_buf`
/// state-machine shape — only the per-record read primitive differs:
/// * BAM: `RawBamReader::read_record(&mut RawRecord)` — zero-copy
///   into a reusable [`fgumi_raw_bam::RawRecord`] buffer.
/// * SAM: `sam::io::Reader::read_record_buf(&Header, &mut RecordBuf)`
///   followed by [`fgumi_raw_bam::encode_record_buf_to_raw`] to bring
///   the record into the same `RawRecord` representation as BAM.
enum TemplateStream {
    Bam(BamTemplateStream),
    Sam(SamTemplateStream),
}

impl TemplateStream {
    fn next_template(&mut self) -> io::Result<Option<AlignedGroup>> {
        match self {
            Self::Bam(s) => s.next_template(),
            Self::Sam(s) => s.next_template(),
        }
    }

    /// `true` if the stream has a peeked record from a prior
    /// `next_template` call — i.e., the aligner emitted at least
    /// one record past the last consumed template boundary.
    fn has_peeked(&self) -> bool {
        match self {
            Self::Bam(s) => s.peeked.is_some(),
            Self::Sam(s) => s.peeked.is_some(),
        }
    }

    /// Probe one more record from the underlying stream. Returns
    /// `true` if a record was read (which counts as the aligner
    /// having emitted extra records past the last token's
    /// boundary). Caller is expected to have already checked
    /// [`Self::has_peeked`]; this method does NOT take a peeked
    /// record into account, only the underlying reader.
    fn probe_trailing(&mut self) -> io::Result<bool> {
        match self {
            Self::Bam(s) => s.probe_trailing(),
            Self::Sam(s) => s.probe_trailing(),
        }
    }
}

/// One group of consecutive same-queryname records from the aligner.
enum AlignedGroup {
    /// The template's alignments (the usual case).
    Template(Template),
    /// bwa's mid-pair split: with mixed single/paired input a `-K` chunk
    /// boundary fell between a pair's two reads, so bwa aligned them as two
    /// unpaired reads and emitted them back to back under the pair's name —
    /// the first read's records, then the second's. The reader splits the
    /// pair's unmapped template the same way (see [`split_pair_into_singles`]).
    /// Produced only when [`SubprocessConfig::accept_mid_pair_split`] is set
    /// (the subprocess presets).
    MidPairSplit(Template, Template),
}

impl AlignedGroup {
    /// The group's queryname.
    fn name(&self) -> &[u8] {
        match self {
            Self::Template(t) | Self::MidPairSplit(t, _) => t.name(),
        }
    }
}

/// Where a same-queryname record group splits if it is bwa's mid-pair split:
/// no record is flagged paired and there are exactly two primaries, in which
/// case the second read's records start at the second primary. `None` for any
/// other group.
fn mid_pair_split_point(records: &[fgumi_raw_bam::RawRecord]) -> Option<usize> {
    use fgumi_raw_bam::flags::{PAIRED, SECONDARY, SUPPLEMENTARY};
    if records.iter().any(|r| r.flags() & PAIRED != 0) {
        return None;
    }
    let mut primaries =
        records.iter().enumerate().filter(|(_, r)| r.flags() & (SECONDARY | SUPPLEMENTARY) == 0);
    let (first, second) = (primaries.next(), primaries.next());
    match (first, second, primaries.next()) {
        (Some(_), Some((split_at, _)), None) => Some(split_at),
        _ => None,
    }
}

/// Group consecutive same-queryname records from `read` into one [`AlignedGroup`].
///
/// Shared by [`BamTemplateStream`] and [`SamTemplateStream`]: only the
/// per-record read primitive differs (zero-copy BAM read vs SAM decode+encode),
/// so it is injected as the `read` closure. `peeked` carries the first record of
/// the *next* template (read one past this template's boundary) across calls,
/// and `name_buf` is a reusable queryname buffer. A group is returned as
/// [`AlignedGroup::MidPairSplit`] only when `accept_mid_pair_split` is set.
/// Returns `Ok(None)` at end of stream.
fn assemble_next_template(
    peeked: &mut Option<fgumi_raw_bam::RawRecord>,
    name_buf: &mut Vec<u8>,
    accept_mid_pair_split: bool,
    mut read: impl FnMut() -> io::Result<Option<fgumi_raw_bam::RawRecord>>,
) -> io::Result<Option<AlignedGroup>> {
    let mut records: Vec<fgumi_raw_bam::RawRecord> = Vec::with_capacity(2);

    if let Some(first) = peeked.take() {
        records.push(first);
    } else {
        match read()? {
            Some(r) => records.push(r),
            None => return Ok(None),
        }
    }

    name_buf.clear();
    name_buf.extend_from_slice(records[0].read_name());

    while let Some(r) = read()? {
        if r.read_name() == name_buf.as_slice() {
            records.push(r);
        } else {
            *peeked = Some(r);
            break;
        }
    }

    // Build the owned queryname only on the error path — the success path (every
    // template but a malformed-record edge case) is the hot loop.
    let build = |records| {
        Template::from_records(records).map_err(|e| {
            io::Error::other(format!(
                "Template::from_records for aligner-emitted queryname '{name}': {e}",
                name = String::from_utf8_lossy(name_buf),
            ))
        })
    };
    // Without `accept_mid_pair_split`, a two-unpaired-primary group falls
    // through to `build`, whose "multiple primaries" error is the loud failure
    // an aligner that never pairs (a free-form `--aligner::command`) needs.
    if let Some(split_at) = mid_pair_split_point(&records).filter(|_| accept_mid_pair_split) {
        let second = records.split_off(split_at);
        return Ok(Some(AlignedGroup::MidPairSplit(build(records)?, build(second)?)));
    }
    Ok(Some(AlignedGroup::Template(build(records)?)))
}

struct BamTemplateStream {
    reader: AlignerBamReader,
    scratch: fgumi_raw_bam::RawRecord,
    peeked: Option<fgumi_raw_bam::RawRecord>,
    /// Reusable per-template name buffer. Re-cleared on each
    /// `next_template` call; avoids ~one `Vec<u8>` allocation
    /// per template (~40 bytes typical) on hot loops.
    name_buf: Vec<u8>,
    /// See [`SubprocessConfig::accept_mid_pair_split`].
    accept_mid_pair_split: bool,
}

impl BamTemplateStream {
    fn new(reader: AlignerBamReader, accept_mid_pair_split: bool) -> Self {
        Self {
            reader,
            scratch: fgumi_raw_bam::RawRecord::new(),
            peeked: None,
            name_buf: Vec::with_capacity(64),
            accept_mid_pair_split,
        }
    }

    fn next_template(&mut self) -> io::Result<Option<AlignedGroup>> {
        // Disjoint field capture (Rust 2021+): the closure borrows
        // `reader`/`scratch` while `peeked`/`name_buf` pass as sibling args, so
        // there's no whole-`self` aliasing conflict. `scratch` is reused across
        // reads, so we clone it out.
        assemble_next_template(
            &mut self.peeked,
            &mut self.name_buf,
            self.accept_mid_pair_split,
            || {
                let n = self.reader.read_record(&mut self.scratch)?;
                Ok(if n == 0 { None } else { Some(self.scratch.clone()) })
            },
        )
    }

    fn probe_trailing(&mut self) -> io::Result<bool> {
        let n = self.reader.read_record(&mut self.scratch)?;
        Ok(n > 0)
    }
}

struct SamTemplateStream {
    reader: AlignerSamReader,
    /// Aligner-emitted header; needed by
    /// [`fgumi_raw_bam::encode_record_buf_to_raw`] per record.
    /// Owned (not a reference) so the stream is self-contained and
    /// the noodles SAM reader's `read_record_buf(&header, ...)` call
    /// can borrow it freely on each iteration.
    header: Arc<Header>,
    /// Scratch `RecordBuf` reused across `read_record_buf` calls.
    scratch: noodles::sam::alignment::RecordBuf,
    /// First record of the next template, encoded to `RawRecord`.
    /// Stashed when `next_template` reads past the current
    /// template's last record.
    peeked: Option<fgumi_raw_bam::RawRecord>,
    name_buf: Vec<u8>,
    /// See [`SubprocessConfig::accept_mid_pair_split`].
    accept_mid_pair_split: bool,
}

impl SamTemplateStream {
    fn new(reader: AlignerSamReader, header: Arc<Header>, accept_mid_pair_split: bool) -> Self {
        Self {
            reader,
            header,
            scratch: noodles::sam::alignment::RecordBuf::default(),
            peeked: None,
            name_buf: Vec::with_capacity(64),
            accept_mid_pair_split,
        }
    }

    /// Read one record-buf from the SAM stream and encode it into a
    /// fresh `RawRecord`. Returns `Ok(None)` on EOF.
    fn read_next_raw(&mut self) -> io::Result<Option<fgumi_raw_bam::RawRecord>> {
        let n = self.reader.read_record_buf(&self.header, &mut self.scratch)?;
        if n == 0 {
            return Ok(None);
        }
        // RecordBuf → RawRecord. `encode_record_buf_to_raw` allocates
        // a fresh output Vec per record (vs zero-copy on the BAM
        // path). For SAM-emitting aligners this is the cost-of-doing
        // business; aligners that can emit BAM avoid it entirely.
        let raw = fgumi_raw_bam::encode_record_buf_to_raw(&self.scratch, &self.header)
            .map_err(|e| io::Error::other(format!("encode_record_buf_to_raw: {e}")))?;
        Ok(Some(raw))
    }

    fn next_template(&mut self) -> io::Result<Option<AlignedGroup>> {
        // Disjoint field capture (see BAM sibling): the closure borrows
        // `reader`/`header`/`scratch` while `peeked`/`name_buf` pass as sibling
        // args. The read primitive is inlined here rather than calling
        // `read_next_raw` (which would borrow all of `self` and conflict with
        // the `&mut self.peeked`/`&mut self.name_buf` args); `read_next_raw`
        // remains as the primitive for `probe_trailing`.
        assemble_next_template(
            &mut self.peeked,
            &mut self.name_buf,
            self.accept_mid_pair_split,
            || {
                let n = self.reader.read_record_buf(&self.header, &mut self.scratch)?;
                if n == 0 {
                    Ok(None)
                } else {
                    fgumi_raw_bam::encode_record_buf_to_raw(&self.scratch, &self.header)
                        .map(Some)
                        .map_err(|e| io::Error::other(format!("encode_record_buf_to_raw: {e}")))
                }
            },
        )
    }

    fn probe_trailing(&mut self) -> io::Result<bool> {
        Ok(self.read_next_raw()?.is_some())
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Step impl
// ──────────────────────────────────────────────────────────────────────────

impl Step for SubprocessAlignStep {
    type Input = BamTemplateBatch;
    type Outputs = OrderedBytesSingle<ZipperBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            // Serial — single-dispatcher-at-a-time is the real
            // requirement (in_tx/out_rx/held_* are not Sync); the
            // framework's mutex enforces it. `Exclusive` would pin
            // dispatch to one worker, which lets any blocking call
            // on that worker (e.g. another step's held-slot flush
            // spin) starve this step.
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.cfg.output_byte_limit }],
            // FIFO: the Serial step emits `ZipperBatch`es in dense serial order,
            // and the downstream (Parallel) `MergeAlignedStep` restores input
            // order from each batch's serial via its own `ByItemOrdinal` output,
            // so no reorder is needed on this edge.
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Fast-fail on any async error from the threads.
        if let Some(e) = self.shared.take_error() {
            return Err(e);
        }

        // Already finalized (the threads were joined and `out_rx` taken on a
        // prior pass). The shared `finished` latch normally stops other
        // workers from re-dispatching a finished Serial step before this is
        // reached, but guard defensively so a re-entry is an idempotent
        // `Finished` rather than an `expect` panic.
        if self.out_rx.is_none() {
            return Ok(StepOutcome::Finished);
        }

        let mut did_work = false;
        let outcome = |did_work: bool| {
            if did_work { StepOutcome::Progress } else { StepOutcome::NoProgress }
        };

        // 2. Drain held output slot first. Backpressure cases below
        //    return `NoProgress` (not `Contention`): `Contention` is
        //    reserved by the framework for Serial-step mutex contention,
        //    not for full output queues.
        if let Some(unpushed) = self.held_out.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => did_work = true,
                Err(again) => {
                    self.held_out.put(again);
                    return Ok(StepOutcome::NoProgress);
                }
            }
        }

        // 3. Pull `ZipperBatch`es from the reader and push them downstream. The
        //    loop exits (without returning) on `Empty`/`Disconnected`. This is
        //    the steady-state pump; it does NOT by itself guarantee `out_rx` is
        //    empty at completion (the reader could buffer a final batch after
        //    this loop's last `try_recv` saw `Empty`). The completion gate
        //    below re-drains `out_rx` after observing `reader.is_finished()`
        //    to close that window — see step 5.
        {
            let out_rx = self.out_rx.as_ref().expect("out_rx Some (checked above)");
            // Loop exits on `Empty` (nothing right now) or `Disconnected`
            // (reader done — if it errored, the slot is set and the next pass
            // surfaces it). Either way the receiver is drained for this call.
            while let Ok(zb) = out_rx.try_recv() {
                match ctx.outputs.push(zb) {
                    Ok(()) => did_work = true,
                    Err(unpushed) => {
                        self.held_out.put(unpushed);
                        return Ok(outcome(did_work));
                    }
                }
            }
        }

        // 4. Ingest phase: while `in_tx` is open, feed the aligner. Held_out
        //    has priority over held_in deliberately: downstream backpressure
        //    must propagate upstream before we accept more input, otherwise the
        //    internal in_chan grows unbounded. The `in_tx` borrow is scoped to
        //    the feed so it ends before 4c may `take()` it.
        let mut close_stdin = false;
        if let Some(in_tx) = self.in_tx.as_ref() {
            let input_drained = ctx.input.is_drained();
            // 4a. Drain held input slot.
            if let Some(batch) = self.held_in.take() {
                match in_tx.try_send(batch) {
                    Ok(()) => did_work = true,
                    Err(TrySendError::Full(b)) => {
                        self.held_in.put(b);
                        return Ok(outcome(did_work));
                    }
                    Err(TrySendError::Disconnected(_)) => {
                        // Prefer the writer's own recorded error (the root cause)
                        // over a generic disconnect message; the writer records
                        // to `error_slot` before dropping its channel end, so a
                        // disconnect here usually means it already failed for a
                        // specific, actionable reason.
                        return Err(self.shared.take_error().unwrap_or_else(|| {
                            io::Error::other(
                                "align-and-merge: writer thread is gone (channel disconnected)",
                            )
                        }));
                    }
                }
            }
            // 4b. Pump from ctx.input → in_chan until either side stalls.
            loop {
                let Some(batch) = ctx.input.pop() else {
                    break;
                };
                match in_tx.try_send(batch) {
                    Ok(()) => did_work = true,
                    Err(TrySendError::Full(b)) => {
                        self.held_in.put(b);
                        return Ok(outcome(did_work));
                    }
                    Err(TrySendError::Disconnected(_)) => {
                        // Prefer the writer's own recorded error (the root cause)
                        // over a generic disconnect message; the writer records
                        // to `error_slot` before dropping its channel end, so a
                        // disconnect here usually means it already failed for a
                        // specific, actionable reason.
                        return Err(self.shared.take_error().unwrap_or_else(|| {
                            io::Error::other(
                                "align-and-merge: writer thread is gone (channel disconnected)",
                            )
                        }));
                    }
                }
            }
            // 4c-guard: ready to close once input is drained and `held_in` is
            // flushed (closing earlier would lose the last input batch).
            close_stdin = input_drained && !self.held_in.is_held();
        }
        if self.in_tx.is_some() {
            // 4c. Close `in_tx` so the writer sees disconnect → aligner stdin
            //     EOF → aligner flushes its remaining output. Re-dispatch then
            //     enters the drain phase below to pump the aligner's tail.
            if close_stdin {
                drop(self.in_tx.take());
                did_work = true;
            }
            return Ok(outcome(did_work));
        }

        // 5. Drain phase: `in_tx` is closed; the aligner is flushing its tail
        //    through the reader. Reaching here means step 3 drained `out_rx`
        //    empty this call. Complete only once the reader has also exited AND
        //    nothing is held — then join the threads (which cannot block, since
        //    the reader is finished and the receiver is drained), wait on the
        //    subprocess, and surface errors.
        let reader_finished = self.reader_thread.as_ref().is_none_or(JoinHandle::is_finished);
        if !self.held_out.is_held() && reader_finished {
            // The reader has exited, so it will send no more batches — but it
            // may have buffered a final batch in `out_rx` in the window between
            // step 3's last `try_recv` (which saw `Empty`) and `is_finished()`
            // flipping true. A finished reader has dropped its `out_tx`, so
            // drain `out_rx` to completion NOW (it yields any remaining batches
            // then `Disconnected`) before joining — otherwise that final batch
            // would be lost when `finalize` drops the receiver. If a batch
            // can't be pushed, park it in `held_out` and re-dispatch.
            {
                let out_rx = self.out_rx.as_ref().expect("out_rx Some (checked above)");
                while let Ok(zb) = out_rx.try_recv() {
                    if let Err(unpushed) = ctx.outputs.push(zb) {
                        self.held_out.put(unpushed);
                        return Ok(StepOutcome::Progress);
                    }
                }
            }
            // out_rx fully drained, reader finished, held_out empty → done.
            return self.finalize();
        }
        Ok(outcome(did_work))
    }
}

impl SubprocessAlignStep {
    /// Completion barrier, reached from `try_run` once the reader has exited,
    /// `out_rx` is empty, and `held_out` is empty (so no aligner output can be
    /// lost). Joins the writer + reader daemons, waits for the aligner
    /// subprocess, and surfaces any error in the documented order
    /// (aligner-exit → reader → writer → `error_slot` fallback). Returns
    /// `StepOutcome::Finished`. The joins cannot deadlock here: the gate
    /// guarantees the reader returned (so the bounded `out_chan` is closed) and
    /// the writer has already seen `in_tx` disconnect.
    fn finalize(&mut self) -> io::Result<StepOutcome> {
        let writer_handle =
            self.writer_thread.take().expect("writer_thread is Some until finalize");
        let reader_handle =
            self.reader_thread.take().expect("reader_thread is Some until finalize");
        // Drop the receiver; the reader has already exited.
        let _ = self.out_rx.take();

        // Join threads. Both have exited by now — the writer when in_rx
        // disconnected (we dropped in_tx), the reader when aligner stdout
        // reached EOF (observed via `is_finished()` in the gate).
        let writer_result = writer_handle.join();
        let reader_result = reader_handle.join();

        // Wait for the subprocess to exit. `AlignerProcess::wait` consumes it.
        let aligner_result = self.aligner.take().expect("aligner is Some until finalize").wait();

        // Surface errors. Order: aligner exit first because its stderr ring is
        // the most actionable diagnostic in the common failure mode
        // (aligner-crash). Reader / writer errors typically cascade from an
        // aligner failure (broken pipes, mismatched template counts).
        aligner_result
            .map_err(|e| io::Error::other(format!("aligner subprocess failed: {e:#}")))?;

        // Then reader: it observes aligner stdout, so its errors (header parse,
        // mismatched template count, etc.) describe aligner output shape.
        match reader_result {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_panic) => return Err(io::Error::other("align-and-merge reader thread panicked")),
        }

        // Writer last: its only failure mode (modulo pipe-broken from an
        // aligner crash, already covered above) is a write_fastq_record bug.
        match writer_result {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_panic) => return Err(io::Error::other("align-and-merge writer thread panicked")),
        }

        // Fallback: surface any error recorded in `shared.error_slot` that
        // didn't come through a thread's `Result`. Both threads return Err
        // when they record to the slot, so a non-empty slot here indicates a
        // future bug — surface it loudly rather than swallow.
        if let Some(e) = self.shared.take_error() {
            return Err(e);
        }

        Ok(StepOutcome::Finished)
    }
}

impl Drop for SubprocessAlignStep {
    /// Best-effort resource cleanup for error paths that bypass the clean
    /// `try_run` completion (`finalize`). The framework's `PipelineSignal`
    /// records the first step error and other workers stop on their next loop
    /// iteration; the step never reaches `finalize` when this happens, so the
    /// threads/aligner are still live and must be torn down here.
    ///
    /// **This is also the cancel-response path.** When
    /// `PipelineSignal::cancel` fires, workers exit at their next
    /// loop iteration, `Pipeline::run` returns, the step's
    /// `Arc` refcount hits zero, and `Drop` runs. `aligner.kill()`
    /// here SIGKILLs bwa; its pipes close immediately so the reader
    /// thread's stdout `read` unblocks within microseconds and the
    /// daemon cascade completes.
    ///
    /// Order is **disconnect-channels-then-kill-then-join** — the
    /// opposite of the clean `finalize` path. Two distinct blocking points
    /// have to be unblocked before we can join:
    ///
    /// * Reader blocked on `aligner_stdout.read(...)` → unblocked by
    ///   killing the subprocess (closing its stdout fd).
    /// * Reader blocked on `out_tx.send(...)` → unblocked by
    ///   dropping `out_rx` (the receiver) so `send` returns
    ///   `SendError`. This case arises when downstream stopped
    ///   accepting batches before drain fired.
    /// * Writer blocked on `in_rx.recv()` → unblocked by dropping
    ///   `in_tx`.
    /// * Writer blocked on `token_tx.send(...)` → unblocked by the
    ///   reader exiting (drops its `token_rx`). The reader exiting
    ///   is gated on the two cases above.
    ///
    /// Steps:
    /// 1. Drop `in_tx` (unblocks writer's `recv`).
    /// 2. Drop `out_rx` (unblocks reader's `send`).
    /// 3. Kill + wait subprocess (unblocks reader's stdout `read`).
    /// 4. Join both threads.
    ///
    /// All errors are silently ignored — we're already on a teardown
    /// path with nowhere to surface them.
    fn drop(&mut self) {
        // Poison the header handle first. If the reader thread *panicked*
        // before resolving the header (the panic unwinds past `reader_loop`'s
        // `Err`-path poison, and the reader is a manual std::thread not covered
        // by the framework's resume_unwind), the downstream WriteBgzfFile would
        // otherwise block forever waiting on the handle. Poisoning here makes it
        // fail with an error instead of hanging. Idempotent: if the header was
        // already resolved on the normal path, `poison` returns `Err` (ignored).
        let _ = self.cfg.header_handle.poison(io::Error::other(
            "SubprocessAlignStep dropped before resolving the output header",
        ));

        drop(self.in_tx.take());
        drop(self.out_rx.take());

        if let Some(mut aligner) = self.aligner.take() {
            // Bounded teardown: only issue the unbounded `wait()` once the child
            // has actually been reaped. If `kill()` timed out (stuck child, e.g.
            // uninterruptible I/O), calling `wait()` would block teardown
            // forever — defeating the 1s kill deadline. In that case we skip
            // `wait()` and let `AlignerProcess::Drop` (itself bounded) clean up.
            if aligner.kill() {
                let _ = aligner.wait();
            }
        }

        if let Some(h) = self.writer_thread.take() {
            let _ = h.join();
        }
        if let Some(h) = self.reader_thread.take() {
            let _ = h.join();
        }
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Backend
// ──────────────────────────────────────────────────────────────────────────

/// The subprocess align backend: wires a [`SubprocessAlignStep`] after the
/// queryname-grouped input, then the shared `Parallel` [`MergeAlignedStep`]
/// that zips its `ZipperBatch` stream, and returns the merged tail.
pub(crate) struct SubprocessBackend {
    /// The resolved aligner shell command (e.g. `bwa-mem3 mem -p -K … <ref> …`).
    pub(crate) command: String,
    /// The aligner's `-K` chunk size (bases), used to derive the in-flight
    /// byte budget.
    pub(crate) chunk_size: u64,
    /// See [`SubprocessConfig::accept_mid_pair_split`].
    pub(crate) accept_mid_pair_split: bool,
}

impl SubprocessBackend {
    /// Minimum pool workers this backend needs for steady-state progress. The
    /// subprocess step spawns two daemon threads (FASTQ writer + SAM/BAM
    /// reader) plus the aligner subprocess, so 4 framework workers are the
    /// floor: source preamble, dispatch, downstream, plus a spare.
    pub(crate) const MIN_WORKERS: usize = 4;
    /// Whether this backend prefers the chain builder's drain-first scheduler.
    /// `false`: the subprocess align path keeps its pre-split scheduling
    /// behavior.
    pub(crate) const PREFERS_DRAIN_FIRST: bool = false;
}

impl AlignBackend for SubprocessBackend {
    fn describe(&self) -> String {
        format!("subprocess aligner: {}", self.command)
    }

    fn wire(
        self: Box<Self>,
        pipeline: &PipelineBuilder,
        input: (StepIdx, BranchIdx),
        ctx: &AlignWiringCtx,
    ) -> anyhow::Result<AlignWired> {
        // fgumi's side of the subprocess route churns the same per-batch
        // buffers (FASTQ out, BAM in, merge); the aligner child gets the same
        // setting through its environment (see `AlignerProcess::spawn`).
        super::retain_freed_memory_unless_user_set();
        let step = SubprocessAlignStep::new(
            SubprocessConfig {
                partial_output_header: Arc::clone(&ctx.partial_output_header),
                header_handle: ctx.header_handle.clone(),
                output_byte_limit: ctx.per_step_byte_limit,
                in_flight_unmapped_budget: in_flight_budget_for_chunk_size(self.chunk_size),
                accept_mid_pair_split: self.accept_mid_pair_split,
            },
            &self.command,
        )
        .map_err(|e| anyhow!("SubprocessAlignStep::new: {e:#}"))?;

        let zipper_tail = pipeline.append_step(step, input);
        // The shared merge step zips the aligner's records with the unmapped
        // halves.
        let tail = pipeline
            .append_step(MergeAlignedStep::from_shared(Arc::clone(&ctx.merge)), zipper_tail);
        Ok(AlignWired {
            tail,
            min_workers: Self::MIN_WORKERS,
            prefers_drain_first: Self::PREFERS_DRAIN_FIRST,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED, REVERSE, SUPPLEMENTARY};

    fn make_test_cfg() -> SubprocessConfig {
        SubprocessConfig {
            partial_output_header: Arc::new(Header::default()),
            header_handle: HeaderHandle::new(),
            output_byte_limit: 1024 * 1024,
            in_flight_unmapped_budget: in_flight_budget_for_chunk_size(150_000_000),
            accept_mid_pair_split: true,
        }
    }

    /// `bash -c 'true'` is a stand-in for a subprocess that exits
    /// cleanly. Lets us verify the spawn path doesn't panic and that
    /// Drop tears down without leaking.
    #[test]
    fn spawn_with_trivial_command_does_not_panic() {
        let cfg = make_test_cfg();
        // The aligner exits before any FASTQ is written. The reader
        // thread will fail to find BGZF magic at stdout EOF, which
        // we expect — we just verify spawn + Drop.
        let step = SubprocessAlignStep::new(cfg, "true").expect("spawn true");
        let _ = step;
    }

    #[test]
    fn profile_advertises_serial_nonsticky_fifo() {
        // Serial + !sticky so any free worker can dispatch when an owner is
        // stuck elsewhere; the `ZipperBatch` output is FIFO (`None`) because the
        // downstream Parallel `MergeAlignedStep` restores order via its own
        // `ByItemOrdinal` output.
        let cfg = make_test_cfg();
        let step = SubprocessAlignStep::new(cfg, "cat").expect("spawn cat");
        let p = step.profile();
        assert_eq!(p.name, "SubprocessAlign");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky, "non-sticky so any free worker can dispatch");
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
        assert_eq!(p.output_queues.len(), 1);
    }

    #[test]
    fn drop_kills_long_running_subprocess() {
        // Spawn a "sleep" subprocess; Drop should kill it within
        // AlignerProcess::Drop's 1-second deadline.
        let cfg = make_test_cfg();
        let start = std::time::Instant::now();
        {
            let _step = SubprocessAlignStep::new(cfg, "sleep 60").expect("spawn sleep");
        } // Drop runs here.
        let elapsed = start.elapsed();
        // AlignerProcess::kill waits up to 1s; thread joins should
        // be near-instant once channels close. Generous bound.
        assert!(
            elapsed < std::time::Duration::from_secs(5),
            "Drop took {elapsed:?}; expected <5s — long-running subprocess not killed?"
        );
    }

    #[test]
    fn empty_stdout_poisons_header_handle() {
        // When the subprocess emits nothing, the reader can't find BGZF
        // magic, returns Err, and *poisons* the header handle rather than
        // resolving it — so a downstream writer blocked on the handle fails
        // instead of hanging.
        let cfg = make_test_cfg();
        let handle = cfg.header_handle.clone();

        // `true` exits immediately — no stdout. Reader sees empty
        // stdout, peek finds 0 bytes, is_bgzf=false → Err →
        // handle.poison.
        let step = SubprocessAlignStep::new(cfg, "true").expect("spawn true");

        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
        while !handle.is_set() && std::time::Instant::now() < deadline {
            thread::yield_now();
            thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(handle.is_set(), "header handle was not resolved (set or poisoned) within 2s");
        let result = handle.try_get().expect("set or poisoned");
        // We expect Err (poisoned by reader's "requires BAM" error
        // path when stdout is empty/non-BGZF).
        assert!(result.is_err(), "expected handle to be poisoned on empty stdout");

        let _ = step;
    }

    /// Build a `RawRecord` for a single read. `flags` controls whether
    /// the record is primary / secondary / supplementary. Used to
    /// construct test fixtures for the writer-loop filter behavior.
    fn make_record(qname: &[u8], flags: u16) -> fgumi_raw_bam::RawRecord {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(qname).flags(flags).sequence(b"ACGT").qualities(b"IIII");
        b.build()
    }

    /// Build a Template whose every record is filtered by
    /// `FASTQ_WRITER_EXCLUDE_FLAGS`. The writer should reject this
    /// input loudly so the user sees the "expected unmapped input"
    /// guidance rather than a misleading "aligner emitted fewer
    /// alignments" downstream error.
    fn make_all_secondary_template(qname: &[u8]) -> Template {
        let rec1 = make_record(qname, fgumi_raw_bam::flags::SECONDARY);
        let rec2 = make_record(qname, fgumi_raw_bam::flags::SUPPLEMENTARY);
        Template::from_records(vec![rec1, rec2]).expect("template with secondary records")
    }

    /// Build a paired-end primary template (one R1 + one R2, both
    /// primary alignments) — the shape `writer_loop` expects from a
    /// realistic unmapped BAM (`fgumi extract` output).
    fn make_paired_primary_template(qname: &[u8]) -> Template {
        let r1 = make_record(qname, PAIRED | FIRST_SEGMENT);
        let r2 = make_record(qname, PAIRED | LAST_SEGMENT);
        Template::from_records(vec![r1, r2]).expect("paired primary template")
    }

    #[test]
    fn writer_loop_errors_on_template_with_no_primary_records() {
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;

        // `cat > /dev/null` consumes our FASTQ writes (if any) so we
        // don't SIGPIPE. We expect zero FASTQ to actually be
        // written before the writer errors out.
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn cat");
        let cat_stdin = cat.stdin.take().expect("cat stdin");

        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(1);
        // Unbounded token channel mirrors production (see comment at
        // the `token_chan` construction site in `SubprocessAlignStep::new`).
        let (token_tx, _token_rx) = channel::<BatchToken>();
        let shared = Arc::new(SharedState::new());

        let bad_template = make_all_secondary_template(b"secondary_only");
        let batch = BamTemplateBatch::new(0, vec![bad_template]);
        in_tx.send(batch).expect("send batch");
        drop(in_tx);

        // Non-engaging budget: this test has no reader to release the gate,
        // so size it so the single batch never blocks (the gate's blocking
        // behavior is covered by `in_flight_gate_*` unit tests).
        let gate = InFlightGate::new(u64::MAX);
        let result = writer_loop(in_rx, token_tx, cat_stdin, Arc::clone(&shared), &gate);

        assert!(result.is_err(), "writer_loop should reject all-filtered template");
        let msg = result.unwrap_err().to_string();
        assert!(
            msg.contains("no primary records") && msg.contains("secondary_only"),
            "error should name the offending queryname: {msg}"
        );

        // Also surfaced into the shared error slot for try_run to
        // observe on its next iteration.
        let slot_err = shared.take_error().expect("error_slot populated");
        assert!(slot_err.to_string().contains("no primary records"));

        let _ = cat.kill();
        let _ = cat.wait();
    }

    /// Regression for the AAM deadlock fixed in commit `069c265`.
    ///
    /// **Invariant**: `token_chan` is unbounded; the writer must
    /// drain `in_rx` regardless of whether anyone is recv'ing tokens.
    #[test]
    fn writer_loop_does_not_deadlock_when_reader_does_not_drain_tokens() {
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;
        use std::time::{Duration, Instant};

        // `cat > /dev/null` consumes our FASTQ writes silently —
        // standing in for an aligner that buffers stdin without
        // emitting anything on stdout.
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn cat");
        let cat_stdin = cat.stdin.take().expect("cat stdin");

        // Match production: in_chan small + bounded (mirrors
        // `IN_CHAN_DEPTH = 2`); token channel unbounded.
        let n_batches: usize = 16;
        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(2);
        let (token_tx, token_rx) = channel::<BatchToken>();
        let shared = Arc::new(SharedState::new());

        let shared_for_worker = Arc::clone(&shared);
        // Non-engaging budget: this test deliberately never drains tokens, so
        // size the in-flight gate so it can't block — the assertion under test
        // is that `token_tx.send` (the unbounded channel) never blocks the
        // writer, which the gate's byte bound must not change for in-flight
        // below budget. The gate's blocking is covered by `in_flight_gate_*`.
        let gate = Arc::new(InFlightGate::new(u64::MAX));
        let gate_for_worker = Arc::clone(&gate);
        let writer_handle = std::thread::Builder::new()
            .name("writer-loop-deadlock-test".into())
            .spawn(move || {
                writer_loop(in_rx, token_tx, cat_stdin, shared_for_worker, &gate_for_worker)
            })
            .expect("spawn writer worker");

        let producer_handle = std::thread::Builder::new()
            .name("producer-deadlock-test".into())
            .spawn(move || {
                for i in 0..n_batches {
                    let qname = format!("read_{i}").into_bytes();
                    let template = make_paired_primary_template(&qname);
                    let batch = BamTemplateBatch::new(i as u64, vec![template]);
                    if in_tx.send(batch).is_err() {
                        // writer dropped early — surface via join.
                        return Err::<(), String>("in_rx closed before producer finished".into());
                    }
                }
                drop(in_tx);
                Ok(())
            })
            .expect("spawn producer worker");

        // Generous timeout: the test should complete in <1 s on a
        // healthy system; 30 s leaves headroom for slow CI runners
        // while still failing fast on a deadlock.
        let deadline = Instant::now() + Duration::from_secs(30);
        loop {
            if writer_handle.is_finished() && producer_handle.is_finished() {
                break;
            }
            assert!(
                Instant::now() < deadline,
                "writer_loop did not exit within 30s — likely a deadlock regression \
                 (token_chan is bounded and reader isn't draining)"
            );
            std::thread::sleep(Duration::from_millis(50));
        }
        producer_handle
            .join()
            .expect("producer thread did not panic")
            .expect("producer did not finish");
        let result = writer_handle.join().expect("writer thread did not panic");
        result.expect("writer_loop returned Err — unexpected for valid input");

        // All N tokens should have accumulated in the (unbounded)
        // token channel — the writer pushed one token per batch
        // without ever blocking.
        let mut n_tokens = 0;
        while token_rx.try_recv().is_ok() {
            n_tokens += 1;
        }
        assert_eq!(
            n_tokens, n_batches,
            "expected {n_batches} tokens pushed into the unbounded channel"
        );

        // No error should have been recorded.
        assert!(shared.take_error().is_none(), "writer_loop must not record an error");

        let _ = cat.kill();
        let _ = cat.wait();
    }

    /// Build a minimal valid header-only BAM file at `path`. Used as
    /// the output of a fake-aligner shell command so the reader can
    /// exercise its real-BAM code path (header parse + `@SQ`
    /// validation + handle resolve + EOF probe) without depending on
    /// a real aligner.
    fn write_minimal_header_only_bam(path: &std::path::Path) {
        use fgumi_bgzf::{BGZF_EOF, InlineBgzfCompressor};
        use std::fs::File;

        let mut header_bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_bytes, &Header::default())
            .expect("write_bam_header");
        let mut hc = InlineBgzfCompressor::new(1);
        hc.write_all(&header_bytes).expect("compress header");
        hc.flush().expect("flush header");

        let mut f = File::create(path).expect("create fixture");
        hc.write_blocks_to(&mut f).expect("write compressed header blocks");
        f.write_all(&BGZF_EOF).expect("write BGZF EOF");
        f.flush().expect("flush fixture");
    }

    #[test]
    fn header_handle_resolves_to_ok_when_aligner_emits_valid_bam() {
        let tmp = tempfile::TempDir::new().expect("tempdir");
        let fixture = tmp.path().join("empty.bam");
        write_minimal_header_only_bam(&fixture);

        let cfg = make_test_cfg();
        let handle = cfg.header_handle.clone();

        // Fake-aligner: emit the pre-staged header-only BAM bytes
        // and exit. The test never pushes FASTQ input, so the
        // aligner doesn't need to consume stdin — `cat <file>`
        // emits the file's bytes and exits when it's done.
        let cmd = format!("cat {}", fixture.display());
        let step = SubprocessAlignStep::new(cfg, &cmd).expect("spawn fake aligner");

        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while !handle.is_set() && std::time::Instant::now() < deadline {
            thread::sleep(std::time::Duration::from_millis(20));
        }
        assert!(handle.is_set(), "handle not resolved within 5s");
        let result = handle.try_get().expect("set");
        assert!(result.is_ok(), "header should resolve to Ok, not poison: {:?}", result.err());

        let _ = step;
    }

    #[test]
    fn empty_stdout_surfaces_aligner_exited_before_emitting() {
        // Confirms the `peek_filled == 0` branch in
        // `peek_aligner_format` produces the stable "exited before
        // emitting" marker rather than going down the BAM or SAM
        // reader paths (both of which would fail later with less
        // actionable errors).
        let cfg = make_test_cfg();
        let handle = cfg.header_handle.clone();
        let _step = SubprocessAlignStep::new(cfg, "true").expect("spawn true");

        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
        while !handle.is_set() && std::time::Instant::now() < deadline {
            thread::sleep(std::time::Duration::from_millis(10));
        }
        assert!(handle.is_set());
        let err = handle.try_get().unwrap().expect_err("poisoned");
        let msg = err.to_string();
        assert!(
            msg.contains(ERR_ALIGNER_EXITED_BEFORE_OUTPUT),
            "expected zero-bytes error path, got: {msg}"
        );
    }

    /// Build a minimal valid SAM-text fixture (header only) at `path`.
    fn write_minimal_header_only_sam(path: &std::path::Path) {
        use std::fs::File;
        let mut f = File::create(path).expect("create sam fixture");
        // Bare-minimum SAM header: one `@HD` line.
        f.write_all(b"@HD\tVN:1.6\n").expect("write SAM header");
        f.flush().expect("flush SAM fixture");
    }

    #[test]
    fn header_handle_resolves_to_ok_when_aligner_emits_sam_text() {
        // Companion to `..._when_aligner_emits_valid_bam` — same
        // assertion but exercising the SAM-text branch through
        // `peek_aligner_format` + `SamTemplateStream`.
        let tmp = tempfile::TempDir::new().expect("tempdir");
        let fixture = tmp.path().join("empty.sam");
        write_minimal_header_only_sam(&fixture);

        let cfg = make_test_cfg();
        let handle = cfg.header_handle.clone();

        let cmd = format!("cat {}", fixture.display());
        let _step = SubprocessAlignStep::new(cfg, &cmd).expect("spawn fake SAM aligner");

        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while !handle.is_set() && std::time::Instant::now() < deadline {
            thread::sleep(std::time::Duration::from_millis(20));
        }
        assert!(handle.is_set(), "handle not resolved within 5s on SAM-text input");
        let result = handle.try_get().expect("set");
        assert!(result.is_ok(), "SAM-text header should resolve handle to Ok: {:?}", result.err(),);
    }

    /// How one group of same-queryname aligner records should assemble: one
    /// template of `n` records, or bwa's mid-pair split into two halves.
    #[derive(Debug, PartialEq, Eq)]
    enum Assembled {
        One(usize),
        Split(usize, usize),
    }

    /// A group of same-name records assembles as one template unless it is bwa's
    /// mid-pair split: two primaries and no record flagged paired, which happens
    /// when a `-K` chunk boundary falls between a pair's reads. A split is cut at
    /// the second primary, so each half keeps its own supplementaries.
    #[rstest::rstest]
    #[case::proper_pair(&[PAIRED | FIRST_SEGMENT, PAIRED | LAST_SEGMENT], Assembled::One(2))]
    #[case::single_with_supplementary(&[0, SUPPLEMENTARY], Assembled::One(2))]
    #[case::mid_pair_split(&[0, REVERSE], Assembled::Split(1, 1))]
    #[case::mid_pair_split_with_supplementaries(
        &[0, SUPPLEMENTARY, REVERSE, SUPPLEMENTARY],
        Assembled::Split(2, 2)
    )]
    #[case::mid_pair_split_first_half_supplementary(&[0, SUPPLEMENTARY, REVERSE], Assembled::Split(2, 1))]
    fn assemble_next_template_recognizes_a_mid_pair_split(
        #[case] flags: &[u16],
        #[case] expected: Assembled,
    ) {
        let recs: Vec<_> = flags.iter().map(|&f| make_record(b"pe10", f)).collect();
        let mut it = recs.into_iter();
        let mut peeked: Option<fgumi_raw_bam::RawRecord> = None;
        let mut name_buf: Vec<u8> = Vec::new();
        let group = assemble_next_template(&mut peeked, &mut name_buf, true, || Ok(it.next()))
            .expect("group ok")
            .expect("group present");
        assert_eq!(group.name(), b"pe10");
        let assembled = match group {
            AlignedGroup::Template(t) => Assembled::One(t.records().len()),
            AlignedGroup::MidPairSplit(a, b) => {
                Assembled::Split(a.records().len(), b.records().len())
            }
        };
        assert_eq!(assembled, expected);
    }

    /// `mid_pair_split_point` fires only for bwa's mid-pair split shape — no
    /// record flagged paired and exactly two primaries — and splits at the
    /// second primary.
    #[rstest::rstest]
    #[case::proper_pair(&[PAIRED | FIRST_SEGMENT, PAIRED | LAST_SEGMENT], None)]
    #[case::one_mate_paired(&[PAIRED | FIRST_SEGMENT, 0], None)]
    #[case::two_unpaired_primaries(&[0, REVERSE], Some(1))]
    #[case::supplementary_between(&[0, SUPPLEMENTARY, REVERSE], Some(2))]
    #[case::three_unpaired_primaries(&[0, 0, REVERSE], None)]
    #[case::single_with_supplementary(&[0, SUPPLEMENTARY], None)]
    fn mid_pair_split_point_matches_only_two_unpaired_primaries(
        #[case] flags: &[u16],
        #[case] expected: Option<usize>,
    ) {
        let recs: Vec<_> = flags.iter().map(|&f| make_record(b"pe10", f)).collect();
        assert_eq!(mid_pair_split_point(&recs), expected);
    }

    /// Without `accept_mid_pair_split` (a free-form `--aligner::command`), two
    /// unpaired primaries under one name are not reinterpreted as a split: the
    /// group fails `Template::from_records` loudly, naming the queryname,
    /// instead of silently splitting a pair an aligner never paired.
    #[test]
    fn two_unpaired_primaries_error_when_the_split_is_not_accepted() {
        let recs = vec![make_record(b"pe10", 0), make_record(b"pe10", REVERSE)];
        let mut it = recs.into_iter();
        let mut peeked: Option<fgumi_raw_bam::RawRecord> = None;
        let mut name_buf: Vec<u8> = Vec::new();
        let Err(err) = assemble_next_template(&mut peeked, &mut name_buf, false, || Ok(it.next()))
        else {
            panic!("two unpaired primaries must error when the split is not accepted");
        };
        let msg = err.to_string();
        assert!(
            msg.contains("Template::from_records for aligner-emitted queryname 'pe10'"),
            "error names the queryname: {msg}"
        );
    }

    /// `split_unmapped_pairs` splits exactly the templates at `split_indices`,
    /// in place and in order, so the unmapped halves line up one-to-one with the
    /// aligner's split output; untouched templates pass through unchanged.
    #[test]
    fn split_unmapped_pairs_splits_only_the_listed_templates_in_order() {
        use crate::pipeline::core::item::Ordered;
        let single = |name: &[u8]| {
            Template::from_records(vec![make_record(name, fgumi_raw_bam::flags::UNMAPPED)])
                .expect("single")
        };
        let pair = |name: &[u8]| {
            Template::from_records(vec![
                make_record(name, PAIRED | FIRST_SEGMENT | fgumi_raw_bam::flags::UNMAPPED),
                make_record(name, PAIRED | LAST_SEGMENT | fgumi_raw_bam::flags::UNMAPPED),
            ])
            .expect("pair")
        };
        let batch = BamTemplateBatch::new(7, vec![single(b"se1"), pair(b"pe2"), pair(b"pe3")]);
        let out = split_unmapped_pairs(batch, &[1]).expect("split ok");
        assert_eq!(out.ordinal(), 7, "the batch serial is preserved");
        let shape: Vec<(Vec<u8>, usize)> =
            out.templates().iter().map(|t| (t.name().to_vec(), t.records().len())).collect();
        assert_eq!(
            shape,
            vec![
                (b"se1".to_vec(), 1),
                (b"pe2".to_vec(), 1),
                (b"pe2".to_vec(), 1),
                (b"pe3".to_vec(), 2),
            ],
            "only pe2 is split, into two single-record halves, in place"
        );
        for half in &out.templates()[1..3] {
            assert_eq!(
                half.records()[0].flags() & (PAIRED | FIRST_SEGMENT | LAST_SEGMENT),
                0,
                "each split half is an unpaired read"
            );
        }
    }

    /// The shared queryname-grouping state machine (delegated to by both
    /// `BamTemplateStream` and `SamTemplateStream`) must group consecutive
    /// same-queryname records into one `Template` and stash the first record of
    /// the *next* template across calls. A mock `read` closure feeds canned
    /// records: readA (R1+R2), then readB (single), then EOF.
    #[test]
    fn assemble_next_template_groups_consecutive_querynames() {
        let recs = vec![
            make_record(b"readA", PAIRED | FIRST_SEGMENT),
            make_record(b"readA", PAIRED | LAST_SEGMENT),
            make_record(b"readB", 0),
        ];
        let mut it = recs.into_iter();
        let mut peeked: Option<fgumi_raw_bam::RawRecord> = None;
        let mut name_buf: Vec<u8> = Vec::new();

        let t1 = assemble_next_template(&mut peeked, &mut name_buf, true, || Ok(it.next()))
            .expect("t1 ok")
            .expect("t1 present");
        let AlignedGroup::Template(t1) = t1 else { panic!("t1 is not a mid-pair split") };
        assert_eq!(t1.name(), b"readA");
        assert_eq!(t1.read_count(), 2, "readA's R1+R2 group into one template");

        let t2 = assemble_next_template(&mut peeked, &mut name_buf, true, || Ok(it.next()))
            .expect("t2 ok")
            .expect("t2 present");
        let AlignedGroup::Template(t2) = t2 else { panic!("t2 is not a mid-pair split") };
        assert_eq!(t2.name(), b"readB");
        assert_eq!(t2.read_count(), 1, "readB is a solo template");

        assert!(
            assemble_next_template(&mut peeked, &mut name_buf, true, || Ok(it.next()))
                .expect("t3 ok")
                .is_none(),
            "stream is exhausted after the last template"
        );
    }

    /// Write a BAM (empty header) holding one record per `(qname, flags)` to
    /// `path`: the staged output of a fake aligner.
    fn write_bam_fixture(path: &std::path::Path, records: &[(&[u8], u16)]) {
        use fgumi_bgzf::{BGZF_EOF, InlineBgzfCompressor};
        use std::fs::File;

        let mut bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut bytes, &Header::default()).expect("write_bam_header");
        for &(qname, flags) in records {
            fgumi_raw_bam::write_raw_record(&mut bytes, &make_record(qname, flags))
                .expect("write record");
        }
        let mut c = InlineBgzfCompressor::new(1);
        c.write_all(&bytes).expect("compress");
        c.flush().expect("flush");
        let mut f = File::create(path).expect("create fixture");
        c.write_blocks_to(&mut f).expect("write compressed blocks");
        f.write_all(&BGZF_EOF).expect("write BGZF EOF");
    }

    /// Write a SAM (`@HD` only) holding one unmapped `ACGT` record per
    /// `(qname, flags)` to `path`: the staged output of a fake aligner.
    fn write_sam_fixture(path: &std::path::Path, records: &[(&[u8], u16)]) {
        let mut text = b"@HD\tVN:1.6\n".to_vec();
        for &(qname, flags) in records {
            let qname = String::from_utf8_lossy(qname);
            writeln!(text, "{qname}\t{flags}\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").expect("format");
        }
        std::fs::write(path, text).expect("write sam fixture");
    }

    /// End to end through `reader_loop_inner` for both aligner output formats: a
    /// pair the aligner split mid-pair (two unpaired primaries under one name)
    /// zips as two single-read mapped templates against its unmapped template
    /// split to match, a single-end template passes through, and an empty batch
    /// still yields an (empty) `ZipperBatch` to keep serials contiguous.
    #[rstest::rstest]
    fn reader_loop_zips_a_mid_pair_split_from_aligner_output(#[values(false, true)] sam: bool) {
        use fgumi_raw_bam::flags::UNMAPPED;
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;

        let tmp = tempfile::TempDir::new().expect("tempdir");
        let aligned: [(&[u8], u16); 3] =
            [(b"pe1", UNMAPPED), (b"pe1", UNMAPPED | REVERSE), (b"se1", UNMAPPED)];
        let fixture = tmp.path().join(if sam { "out.sam" } else { "out.bam" });
        if sam {
            write_sam_fixture(&fixture, &aligned);
        } else {
            write_bam_fixture(&fixture, &aligned);
        }
        let mut cat = Command::new("cat")
            .arg(&fixture)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .spawn()
            .expect("spawn cat");
        let stdout = cat.stdout.take().expect("cat stdout");

        let (token_tx, token_rx) = channel::<BatchToken>();
        let (out_tx, out_rx) = sync_channel::<ZipperBatch>(4);
        let empty = BamTemplateBatch::new(0, Vec::new());
        token_tx.send(BatchToken { unmapped: empty, n_templates: 0, serial: 0 }).unwrap();
        let single = Template::from_records(vec![make_record(b"se1", UNMAPPED)]).unwrap();
        let unmapped = BamTemplateBatch::new(1, vec![make_paired_primary_template(b"pe1"), single]);
        token_tx.send(BatchToken { unmapped, n_templates: 2, serial: 1 }).unwrap();
        drop(token_tx);

        let cfg = make_test_cfg();
        reader_loop_inner(token_rx, &out_tx, stdout, &cfg, &InFlightGate::new(u64::MAX))
            .expect("reader_loop_inner");
        drop(out_tx);
        let _ = cat.wait();

        let batches: Vec<ZipperBatch> = out_rx.iter().collect();
        assert_eq!(batches.iter().map(|b| b.serial).collect::<Vec<_>>(), vec![0, 1]);
        assert!(batches[0].mapped.is_empty() && batches[0].unmapped.templates().is_empty());
        let shape = |ts: &[Template]| {
            ts.iter().map(|t| (t.name().to_vec(), t.read_count())).collect::<Vec<_>>()
        };
        let expected = vec![(b"pe1".to_vec(), 1), (b"pe1".to_vec(), 1), (b"se1".to_vec(), 1)];
        assert_eq!(shape(&batches[1].mapped), expected, "mapped side");
        assert_eq!(shape(batches[1].unmapped.templates()), expected, "unmapped side");
        for half in &batches[1].unmapped.templates()[..2] {
            assert_eq!(half.records()[0].flags() & PAIRED, 0, "split half is unpaired");
        }
    }
}
