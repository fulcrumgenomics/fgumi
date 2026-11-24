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
use std::time::{Duration, Instant};

use anyhow::anyhow;
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::aligner::{AlignerProcess, AlignerScratchDir, SubprocessRuntime};
use crate::commands::fastq::{
    FastqRecordBuffers, UNALIGNED_SAM_HEADER, UnalignedReadFormat, write_fastq_record,
    write_unaligned_sam_record,
};
use crate::pipeline::core::builder::PipelineBuilder;
use crate::pipeline::core::header::HeaderHandle;
use crate::pipeline::core::item::HeapSize;
use crate::pipeline::core::topology::{BranchIdx, StepIdx};
use crate::pipeline::core::{
    BranchOrdering, ExternalWaitClaim, HeldSlot, OrderedBytesSingle, QueueSpec, Step, StepCtx,
    StepKind, StepOutcome, StepProfile, Unpushed,
};
use crate::pipeline::steps::align::merge::MergeAlignedStep;
use crate::pipeline::steps::align::{
    AlignBackend, AlignWired, AlignWiringCtx, ConsumerGoneGuard, IN_FLIGHT_BUDGET_ENV,
    IN_FLIGHT_CEILING_TIMEOUT, IN_FLIGHT_CEILING_TIMEOUT_ENV, IN_FLIGHT_STALL_GRACE, InFlightGate,
    ZipperBatch, ceiling_timeout_override, check_meth_aligner_version,
    in_flight_budget_for_chunk_size, in_flight_budget_override, is_primary_for_alignment,
    merge_aligner_header, no_primary_records_message, read_ahead_ceiling, split_pair_into_singles,
    validate_sq_consistency,
};
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::template::Template;

/// Number of aligner stderr lines retained for failure diagnostics.
const ALIGNER_STDERR_RING_SIZE: usize = 50;

/// Grace period for the reader to publish a root cause after a writer cascade.
/// A failed stdin pipe does not guarantee stdout EOF, so only this failure path
/// is bounded; healthy alignment can take as long as it needs.
const CASCADE_DRAIN_TIMEOUT: Duration = Duration::from_secs(1);

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
    /// `-K` chunk size via [`in_flight_budget_for_chunk_size`]. It is where
    /// the gate starts: a stalled aligner grows it up to
    /// `in_flight_read_ahead_ceiling` (see [`InFlightGate`]).
    pub(crate) in_flight_unmapped_budget: u64,

    /// Hard ceiling the in-flight budget may grow to when the aligner needs
    /// more input than the budget before it emits; past it the run fails.
    /// From [`read_ahead_ceiling`] (`--aligner::max-read-ahead`).
    pub(crate) in_flight_read_ahead_ceiling: u64,

    /// How long a batch at the ceiling may see no aligner output consumed
    /// before the run fails: [`IN_FLIGHT_CEILING_TIMEOUT`], unless the hidden
    /// [`IN_FLIGHT_CEILING_TIMEOUT_ENV`] overrides it.
    pub(crate) in_flight_ceiling_timeout: std::time::Duration,

    /// Whether a same-queryname output group with two unpaired primaries is
    /// accepted as bwa's mid-pair split (see [`AlignedGroup::MidPairSplit`]).
    /// Set only for the subprocess presets (`mem -p -K`); otherwise the group
    /// fails `Template::from_records`'s "multiple primaries" check, as it did
    /// before the split was recognized.
    pub(crate) accept_mid_pair_split: bool,

    /// Whether the aligner command passes bwa-mem3's `--meth`
    /// ([`crate::aligner::command_passes_meth`]): the reader thread then
    /// refuses a bwa-mem3 older than the `--meth` `NM` fix
    /// ([`check_meth_aligner_version`]).
    pub(crate) passes_meth: bool,
    /// Format the writer thread writes the reads in (`--aligner::input-format`):
    /// interleaved FASTQ, or unaligned SAM (an `@HD` line, then one record per
    /// read). Either way R1 precedes R2 and secondary/supplementary records are
    /// skipped, so the reader's pairing is unchanged.
    pub(crate) input_format: UnalignedReadFormat,

    /// Require complete primary pairs before writing a template (STAR SAM PE).
    pub(crate) require_complete_pairs: bool,

    /// Maximum time from spawning a slow-starting aligner to its first output
    /// (see [`SubprocessRuntime::startup_timeout`]). Until output arrives its
    /// silence is an external wait for the watchdog, not a deadlock. Expiry
    /// fails the run; `None` uses the normal watchdog from the outset.
    pub(crate) startup_timeout: Option<std::time::Duration>,

    /// A scratch directory the aligner writes into, held until the step (and
    /// its reader thread) are gone, after the aligner has been reaped. Never
    /// read: holding it is what keeps the directory alive.
    #[allow(dead_code)]
    pub(crate) work_dir: Option<AlignerScratchDir>,
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
/// consequences of the first. The one exception is a *cascade* error (see
/// [`Self::record_cascade_error`]): the writer stopping because the other end
/// went away, either a broken pipe to the aligner's stdin or the reader no
/// longer taking batch tokens. That is a symptom, never a cause, and the reader
/// that stopped the aligner (by closing its stdout on bad output) records its
/// error only after its stdout and token receiver have dropped, so a cascade
/// error is held aside and a root-cause error recorded during the bounded drain
/// still wins.
struct SharedState {
    error_slot: Mutex<ErrorSlot>,
    /// The aligner's start-up wait, for a [`SubprocessConfig::startup_timeout`]
    /// aligner.
    startup: Mutex<StartupWait>,
}

/// The external wait a [`SubprocessConfig::startup_timeout`] aligner's start-up is
/// counted as. Only a pipeline worker can claim the watchdog's gauge, so
/// `try_run` claims it and the reader thread, which sees the first output byte
/// arrive, releases it.
#[derive(Default)]
struct StartupWait {
    /// The aligner has written its first byte, or the reader has exited.
    over: bool,
    /// Time of the first successful stdout read, independent of dispatch timing.
    first_output_at: Option<std::time::Instant>,
    claim: Option<ExternalWaitClaim>,
}

/// See [`SharedState`].
#[derive(Default)]
struct ErrorSlot {
    /// The first root-cause error.
    error: Option<io::Error>,
    /// When the first cascade was recorded and its error. The deadline starts
    /// here, even if dispatch has already closed input or is backpressured.
    cascade: Option<(Instant, io::Error)>,
    /// The first returned dispatch error, retained before the Serial guard is
    /// released so another worker cannot enter a successful or blocking drain.
    reported_error: Option<(io::ErrorKind, String)>,
}

impl SharedState {
    fn new() -> Self {
        Self {
            error_slot: Mutex::new(ErrorSlot::default()),
            startup: Mutex::new(StartupWait::default()),
        }
    }

    /// From a pipeline worker: count the aligner's start-up as an external
    /// wait, unless it is over or already counted.
    fn claim_startup_wait(
        &self,
        started_at: std::time::Instant,
        timeout: std::time::Duration,
    ) -> io::Result<()> {
        let mut startup = self.startup.lock();
        let observed_elapsed = startup.first_output_at.map(|at| at.duration_since(started_at));
        let elapsed = observed_elapsed.or_else(|| (!startup.over).then(|| started_at.elapsed()));
        if elapsed.is_some_and(|elapsed| elapsed >= timeout) {
            startup.over = true;
            startup.claim = None;
            return Err(io::Error::new(
                io::ErrorKind::TimedOut,
                format!(
                    "aligner startup timed out after {timeout:?} without output; increase \
                     --aligner::startup-timeout (seconds) for a slower genome load"
                ),
            ));
        }
        if !startup.over && startup.claim.is_none() {
            startup.claim = ExternalWaitClaim::begin();
        }
        Ok(())
    }

    /// Record first output as soon as a read returns bytes, not after the
    /// format probe has filled or the dispatcher runs again.
    fn observe_output(&self) {
        let observed_at = std::time::Instant::now();
        let mut startup = self.startup.lock();
        startup.first_output_at.get_or_insert(observed_at);
        startup.over = true;
        startup.claim = None;
    }

    /// End the aligner's start-up (first byte written, or reader gone), releasing any
    /// external wait it was counted as.
    fn end_startup_wait(&self) {
        let mut startup = self.startup.lock();
        startup.over = true;
        startup.claim = None;
    }

    /// Reconstruct the failure already returned by a previous dispatch.
    fn reported_dispatch_error(&self) -> Option<io::Error> {
        self.error_slot
            .lock()
            .reported_error
            .as_ref()
            .map(|(kind, message)| io::Error::new(*kind, message.clone()))
    }

    /// Retain the selected failure before the dispatcher releases its Serial
    /// guard, preserving its kind and message if another worker publishes first.
    fn record_dispatch_error(&self, error: &io::Error) {
        self.error_slot
            .lock()
            .reported_error
            .get_or_insert_with(|| (error.kind(), error.to_string()));
    }

    /// Record a root-cause error if none has been recorded yet. First writer
    /// wins.
    fn record_error(&self, err: io::Error) {
        let mut guard = self.error_slot.lock();
        if guard.error.is_none() {
            guard.error = Some(err);
        }
    }

    /// Record a cascade error: the writer stopped because the aligner stopped
    /// reading its stdin or the reader stopped taking tokens. It is reported
    /// only if no root-cause error arrives before the drain deadline.
    fn record_cascade_error(&self, err: io::Error) {
        let mut guard = self.error_slot.lock();
        if guard.cascade.is_none() {
            guard.cascade = Some((Instant::now(), err));
        }
    }

    /// Take the recorded root-cause error, if any, leaving a cascade error.
    fn take_root_cause(&self) -> Option<io::Error> {
        self.error_slot.lock().error.take()
    }

    /// Whether a cascade error has been recorded.
    fn has_cascade(&self) -> bool {
        self.error_slot.lock().cascade.is_some()
    }

    /// Take a root cause immediately, or a cascade once its reader grace period
    /// expires. One lock makes a root cause published by the deadline win.
    fn take_dispatch_error(&self, now: Instant) -> Option<io::Error> {
        let mut guard = self.error_slot.lock();
        if let Some(error) = guard.error.take() {
            return Some(error);
        }
        if guard.cascade.as_ref().is_some_and(|(recorded_at, _)| {
            now.saturating_duration_since(*recorded_at) >= CASCADE_DRAIN_TIMEOUT
        }) {
            return guard.cascade.take().map(|(_, error)| error);
        }
        None
    }

    /// Take the recorded root-cause error, else the cascade error.
    fn take_error(&self) -> Option<io::Error> {
        let mut guard = self.error_slot.lock();
        guard.error.take().or_else(|| guard.cascade.take().map(|(_, error)| error))
    }
}

/// Whether a writer error is a broken pipe to the aligner's stdin: the aligner
/// stopped reading (it exited, or was killed because its output was rejected),
/// so the writer's failure is a consequence, not a cause.
fn is_broken_pipe(err: &anyhow::Error) -> bool {
    err.chain().any(|cause| {
        cause.downcast_ref::<io::Error>().is_some_and(|e| e.kind() == io::ErrorKind::BrokenPipe)
    })
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
    /// Startup deadline is measured from spawn, not from the first dispatch.
    started_at: std::time::Instant,

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
        let started_at = std::time::Instant::now();
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
        let gate = Arc::new(InFlightGate::new(
            cfg.in_flight_unmapped_budget,
            cfg.in_flight_read_ahead_ceiling,
            IN_FLIGHT_STALL_GRACE,
            cfg.in_flight_ceiling_timeout,
        ));

        let writer_thread = {
            let input_format = cfg.input_format;
            let require_complete_pairs = cfg.require_complete_pairs;
            let shared = Arc::clone(&shared);
            let gate = Arc::clone(&gate);
            thread::Builder::new()
                .name("aam-fastq-writer".into())
                .spawn(move || {
                    writer_loop(
                        in_rx,
                        token_tx,
                        aligner_stdin,
                        shared,
                        &gate,
                        input_format,
                        require_complete_pairs,
                    )
                })
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
            started_at,
            held_in: HeldSlot::new(),
            held_out: HeldSlot::new(),
            name: "SubprocessAlign",
        })
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Writer thread
// ──────────────────────────────────────────────────────────────────────────

/// Drain `in_rx` to the aligner's stdin as interleaved FASTQ or, per
/// `input_format`, unaligned SAM (its `@HD` line first), pushing one
/// `BatchToken` per batch to `token_tx` so the reader can pair alignments back
/// by count.
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
    input_format: UnalignedReadFormat,
    require_complete_pairs: bool,
) -> io::Result<()> {
    let mut writer = BufWriter::with_capacity(FASTQ_WRITER_BUF_BYTES, aligner_stdin);
    let mut fastq_buffers = FastqRecordBuffers::with_capacity(512);
    let mut next_serial: u64 = 0;

    if input_format == UnalignedReadFormat::Sam
        && let Err(e) = writer.write_all(UNALIGNED_SAM_HEADER)
    {
        let msg = format!("align-and-merge writer: SAM header: {e}");
        shared.record_error(io::Error::other(msg.clone()));
        return Err(io::Error::other(msg));
    }

    while let Ok(batch) = in_rx.recv() {
        // Reserve this batch's unmapped bytes against the in-flight budget
        // before feeding it to the aligner. Blocks (transiently) if the
        // backlog is at budget; returns `false` only if the reader has
        // exited, in which case we stop (the channel send below would also
        // fail). Bounds the otherwise-unbounded token backlog (issue #382).
        // Errs when the aligner reads past the read-ahead ceiling without
        // emitting; nothing of this batch has been written then.
        let batch_bytes = batch.heap_size() as u64;
        // Before acquiring, hand the aligner everything already written: if the
        // gate makes the writer wait, the tail of the previous batch must not
        // sit in the FASTQ buffer, since an aligner that needs it to finish a
        // chunk would emit nothing, so the gate would wait out its grace and
        // grow for no reason (or fail at the ceiling). Flushing every batch
        // costs at most one extra write per batch.
        if let Err(e) = writer.flush() {
            let msg = format!("align-and-merge writer: flush: {e}");
            return Err(record_writer_error(&shared, e.kind() == io::ErrorKind::BrokenPipe, msg));
        }
        match gate.acquire(batch_bytes) {
            Ok(true) => {}
            Ok(false) => return Ok(()),
            Err(e) => {
                let msg = format!("align-and-merge writer: {e}");
                shared.record_error(io::Error::other(msg.clone()));
                return Err(io::Error::other(msg));
            }
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
            if require_complete_pairs {
                use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED};
                let mask = PAIRED | FIRST_SEGMENT | LAST_SEGMENT;
                let complete =
                    template.r1().is_some_and(|r| r.flags() & mask == PAIRED | FIRST_SEGMENT)
                        && template.r2().is_some_and(|r| r.flags() & mask == PAIRED | LAST_SEGMENT);
                if !complete {
                    let msg = format!(
                        "align-and-merge writer: STAR requires a complete primary pair for template '{}'",
                        String::from_utf8_lossy(template.name()),
                    );
                    shared.record_error(io::Error::other(msg.clone()));
                    return Err(io::Error::other(msg));
                }
            }
            let mut wrote_any = false;
            for record in &template.records {
                let flags = record.flags();
                // Shared with the in-process prepare step so both backends select
                // byte-identical primary reads (skip SECONDARY/SUPPLEMENTARY).
                if !is_primary_for_alignment(flags) {
                    continue;
                }
                let result = match input_format {
                    UnalignedReadFormat::Fastq => write_fastq_record(
                        &mut writer,
                        record,
                        flags,
                        /* no_suffix */ true,
                        &mut fastq_buffers,
                        /* umi_header */ None,
                    ),
                    UnalignedReadFormat::Sam => {
                        write_unaligned_sam_record(&mut writer, record, flags, &mut fastq_buffers)
                    }
                };
                if let Err(e) = result {
                    // Record + return the same full error message so
                    // `finalize`'s join path and the shared
                    // error_slot fallback path both see the same
                    // text — no information loss either way. Same
                    // pattern as the no-primary path below.
                    let format = match input_format {
                        UnalignedReadFormat::Fastq => "FASTQ",
                        UnalignedReadFormat::Sam => "SAM",
                    };
                    let msg = format!(
                        "align-and-merge writer: writing {format} for template '{name}': {e:#}",
                        name = String::from_utf8_lossy(template.name()),
                    );
                    return Err(record_writer_error(&shared, is_broken_pipe(&e), msg));
                }
                wrote_any = true;
            }
            if !wrote_any {
                // Shared with the in-process prepare step via
                // `no_primary_records_message` so the two backends' hard-error
                // wording cannot drift; the "writer:" prefix stays subprocess-local.
                let msg = format!(
                    "align-and-merge writer: {}",
                    no_primary_records_message(template.name())
                );
                shared.record_error(io::Error::other(msg.clone()));
                return Err(io::Error::other(msg));
            }
            // KNOWN DIVERGENCE — resolve in the align-and-merge WIRING PR.
            // Outside STAR (which requires complete pairs above), this guard
            // only rejects a template where EVERY record is filtered
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
            // Reader thread is gone — we cannot make progress. Its own error
            // (if any) is the root cause, but it records that only after
            // dropping the token receiver, so record a cascade error: the step
            // then drains and waits for the reader rather than failing at once
            // with a generic "writer is gone". If the reader records nothing,
            // the drain deadline reports that input was left unaligned.
            shared.record_cascade_error(io::Error::other(
                "align-and-merge writer: the reader stopped while input remained",
            ));
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
        return Err(record_writer_error(&shared, e.kind() == io::ErrorKind::BrokenPipe, msg));
    }
    drop(writer);
    Ok(())
}

/// Record a writer I/O failure in `shared` and return it. A broken pipe to the
/// aligner's stdin is a cascade error (see [`SharedState::record_cascade_error`])
/// and keeps its [`io::ErrorKind::BrokenPipe`] kind, so `finalize` can rank it
/// below the error that stopped the aligner; anything else is a root cause.
fn record_writer_error(shared: &SharedState, broken_pipe: bool, msg: String) -> io::Error {
    if broken_pipe {
        shared.record_cascade_error(io::Error::new(io::ErrorKind::BrokenPipe, msg.clone()));
        io::Error::new(io::ErrorKind::BrokenPipe, msg)
    } else {
        shared.record_error(io::Error::other(msg.clone()));
        io::Error::other(msg)
    }
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
    /// Ends the aligner's start-up wait on every reader exit, a panic
    /// included, so the watchdog never keeps treating a dead reader's
    /// start-up as an external wait.
    struct EndStartup<'a>(&'a SharedState);
    impl Drop for EndStartup<'_> {
        fn drop(&mut self) {
            self.0.end_startup_wait();
        }
    }
    let _end_startup = EndStartup(&shared);
    let result = reader_loop_inner(token_rx, &out_tx, aligner_stdout, &cfg, &shared, gate);
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
    shared: &SharedState,
    gate: &InFlightGate,
) -> io::Result<()> {
    // 1. Peek a BGZF header to detect BAM vs SAM text. Returns the
    //    peeked bytes prepended to a chained Read so we don't lose
    //    the prefix.
    let (is_bgzf, peek_filled, stdout) = peek_aligner_format(aligner_stdout, shared)?;
    // The aligner has written (or exited without writing): its start-up is
    // over, and from here on the watchdog judges the align stage by its
    // progress.
    shared.end_startup_wait();

    if peek_filled == 0 {
        // Aligner exited before emitting any output — usually an
        // aligner crash or misconfiguration (wrong binary path,
        // missing index files, segfault before reading its input).
        // Distinguished from "non-BGZF content" so the message
        // doesn't misleadingly point at SAM/BAM format conversion.
        let marker = ERR_ALIGNER_EXITED_BEFORE_OUTPUT;
        return Err(io::Error::other(format!(
            "SubprocessAlignStep: aligner {marker}. Likely causes: subprocess \
             startup error (missing binary, missing index files), aligner \
             argument error, the aligner reads a different input format than \
             --aligner::input-format, or the aligner segfaulted before reading \
             its input. Check the aligner's stderr above for the specific failure."
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
    if cfg.passes_meth {
        check_meth_aligner_version(&aligner_header)?;
    }

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
    // header is dropped in either case; the existing value wins. A handle
    // poisoned by teardown (the step failed, e.g. at the read-ahead ceiling,
    // while this reader still waited for the header) is not a wiring bug, so
    // it is not warned about.
    if cfg.header_handle.set(merged).is_err()
        && !matches!(cfg.header_handle.try_get(), Some(Err(_)))
    {
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
            // Blocked downstream is not an aligner stall (see `downstream_wait`).
            let sent = {
                let _downstream = gate.downstream_wait();
                out_tx.send(zb)
            };
            if sent.is_err() {
                return Ok(()); // downstream gone
            }
            continue;
        }

        let mut mapped_templates: Vec<Template> = Vec::with_capacity(token.n_templates);
        // Unmapped-template indices bwa split mid-pair; their unmapped halves are
        // split to match after the loop (rare: needs mixed single/paired input).
        let mut split_indices: Vec<usize> = Vec::new();

        for i in 0..token.n_templates {
            let mapped = template_stream.next_template(gate)?;

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
        // Blocked downstream is not an aligner stall (see `downstream_wait`).
        let sent = {
            let _downstream = gate.downstream_wait();
            out_tx.send(zb)
        };
        if sent.is_err() {
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
    if template_stream.probe_trailing(gate).map_err(|e| {
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
fn peek_aligner_format(
    mut stdout: ChildStdout,
    shared: &SharedState,
) -> io::Result<(bool, usize, Box<dyn Read + Send>)> {
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
            Ok(n) => {
                if peek_filled == 0 {
                    shared.observe_output();
                }
                peek_filled += n;
            }
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
    fn next_template(&mut self, gate: &InFlightGate) -> io::Result<Option<AlignedGroup>> {
        match self {
            Self::Bam(s) => s.next_template(gate),
            Self::Sam(s) => s.next_template(gate),
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
    fn probe_trailing(&mut self, gate: &InFlightGate) -> io::Result<bool> {
        match self {
            Self::Bam(s) => s.probe_trailing(gate),
            Self::Sam(s) => s.probe_trailing(gate),
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
    /// the first read's records, then the second's. The in-process backend
    /// splits such a pair the same way (see [`split_pair_into_singles`]).
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

    fn next_template(&mut self, gate: &InFlightGate) -> io::Result<Option<AlignedGroup>> {
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
                Ok(if n == 0 {
                    None
                } else {
                    gate.record_progress();
                    Some(self.scratch.clone())
                })
            },
        )
    }

    fn probe_trailing(&mut self, gate: &InFlightGate) -> io::Result<bool> {
        let n = self.reader.read_record(&mut self.scratch)?;
        if n > 0 {
            gate.record_progress();
        }
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
    fn read_next_raw(
        &mut self,
        gate: &InFlightGate,
    ) -> io::Result<Option<fgumi_raw_bam::RawRecord>> {
        let n = self.reader.read_record_buf(&self.header, &mut self.scratch)?;
        if n == 0 {
            return Ok(None);
        }
        gate.record_progress();
        // RecordBuf → RawRecord. `encode_record_buf_to_raw` allocates
        // a fresh output Vec per record (vs zero-copy on the BAM
        // path). For SAM-emitting aligners this is the cost-of-doing
        // business; aligners that can emit BAM avoid it entirely.
        let raw = fgumi_raw_bam::encode_record_buf_to_raw(&self.scratch, &self.header)
            .map_err(|e| io::Error::other(format!("encode_record_buf_to_raw: {e}")))?;
        Ok(Some(raw))
    }

    fn next_template(&mut self, gate: &InFlightGate) -> io::Result<Option<AlignedGroup>> {
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
                    gate.record_progress();
                    fgumi_raw_bam::encode_record_buf_to_raw(&self.scratch, &self.header)
                        .map(Some)
                        .map_err(|e| io::Error::other(format!("encode_record_buf_to_raw: {e}")))
                }
            },
        )
    }

    fn probe_trailing(&mut self, gate: &InFlightGate) -> io::Result<bool> {
        Ok(self.read_next_raw(gate)?.is_some())
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
        if let Some(error) = self.shared.reported_dispatch_error() {
            return Err(error);
        }
        let result = self.try_run_inner(ctx);
        if let Err(error) = &result {
            self.shared.record_dispatch_error(error);
        }
        result
    }
}

impl SubprocessAlignStep {
    /// Run one dispatch, with the outer entry point retaining any returned error.
    fn try_run_inner(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Fast-fail on a root-cause error from the threads. This is how a
        //    reader error is normally reported: the reader records it as it
        //    exits, before the step reaches `finalize`. A writer cascade error
        //    alone gets a bounded grace period for the reader's error to arrive.
        //    Check before all backpressure returns, even after input was closed.
        if let Some(err) = self.shared.take_root_cause() {
            return Err(err);
        }
        // An already-finished aligner can supply the cause of a writer cascade.
        // Collect it before returning the cascade, but never wait on a live child,
        // an I/O thread, or its stderr relay from this failure path.
        if self.shared.has_cascade()
            && self.reader_thread.as_ref().is_some_and(JoinHandle::is_finished)
            && self.writer_thread.as_ref().is_some_and(JoinHandle::is_finished)
            && let Some(aligner) = self.aligner.as_mut()
            && aligner.is_finished()?
        {
            return self.finalize();
        }
        if let Some(err) = self.shared.take_dispatch_error(Instant::now()) {
            return Err(err);
        }
        if let Some(timeout) = self.cfg.startup_timeout {
            self.shared.claim_startup_wait(self.started_at, timeout)?;
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
                    Err(TrySendError::Disconnected(_)) => return self.writer_gone(),
                }
            }
            // 4b. Pump from ctx.input → in_chan until either side stalls.
            while let Some(batch) = ctx.input.pop() {
                match in_tx.try_send(batch) {
                    Ok(()) => did_work = true,
                    Err(TrySendError::Full(b)) => {
                        self.held_in.put(b);
                        return Ok(outcome(did_work));
                    }
                    Err(TrySendError::Disconnected(_)) => return self.writer_gone(),
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
        let writer_finished = self.writer_thread.as_ref().is_none_or(JoinHandle::is_finished);
        // A finished reader does not imply a finished writer or child. In the
        // cascade case, keep dispatching until step 1 reports the reader's root
        // cause or the expired cascade; never enter finalize's blocking waits.
        // Check the slot AFTER writer completion so its final flush cannot race
        // this gate and publish a cascade after we have committed to finalize.
        if !self.held_out.is_held()
            && reader_finished
            && writer_finished
            && !self.shared.has_cascade()
        {
            // Both threads are now finished, so no later root cause can race
            // this final check or leave us waiting on the child after a failure.
            if let Some(err) = self.shared.take_root_cause() {
                return Err(err);
            }
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
    /// The writer thread exited while input remained. A root-cause error it (or
    /// the reader) recorded is returned at once. A cascade error alone (a
    /// broken pipe to the aligner's stdin, or the reader no longer taking
    /// tokens) usually follows the reader rejecting the aligner's output and
    /// closing its stdout, an error the reader may still be recording; so the
    /// step stops feeding the aligner and enters the drain phase. The reader's
    /// error then surfaces at `try_run`'s step 1, and the cascade error only
    /// if no root cause arrives before the grace period expires.
    ///
    /// Unlike a successful drain, this waits at most [`CASCADE_DRAIN_TIMEOUT`]
    /// from the first recorded cascade. If stdout stays open, `try_run` reports
    /// the cascade and Drop kills the aligner rather than waiting for EOF.
    fn writer_gone(&mut self) -> io::Result<StepOutcome> {
        if let Some(err) = self.shared.take_root_cause() {
            return Err(err);
        }
        if self.shared.has_cascade() {
            drop(self.in_tx.take());
            drop(self.held_in.take());
            return Ok(StepOutcome::Progress);
        }
        Err(io::Error::other("align-and-merge: writer thread is gone (channel disconnected)"))
    }

    /// Completion barrier after a successful output drain, or after a cascade
    /// once both I/O threads, child, and stderr relay have finished. A successful
    /// drain requires `out_rx` and `held_out` to be empty; the cascade path always
    /// fails, so any buffered output is discarded.
    /// Joins the writer + reader daemons, waits for the aligner
    /// subprocess, and surfaces the error [`rank_align_errors`] puts first (the
    /// reader's, with the aligner's exit after it; else the aligner's exit;
    /// else the writer's), with the `error_slot` as fallback. A reader error
    /// usually never gets here: `try_run` returns it as soon as it is recorded.
    /// Returns
    /// `StepOutcome::Finished`. Both callers guarantee that the I/O threads have
    /// finished, so these joins cannot wait on live pipes or channels.
    fn finalize(&mut self) -> io::Result<StepOutcome> {
        let writer_handle =
            self.writer_thread.take().expect("writer_thread is Some until finalize");
        let reader_handle =
            self.reader_thread.take().expect("reader_thread is Some until finalize");
        // Drop the receiver; the reader has already exited.
        let _ = self.out_rx.take();

        // Both I/O threads have exited: after input/output EOF on a successful
        // drain, or after a cascade has stopped the writer and reader.
        let writer_result = writer_handle.join();
        let reader_result = reader_handle.join();

        // Wait for the subprocess to exit. `AlignerProcess::wait` consumes it.
        //
        // If the reader failed (an error, or a panic), nothing reads the
        // aligner's output any more and the aligner may not exit on its own:
        // tear its process group down (bounded) instead of waiting without
        // bound. An aligner that had already exited keeps its own exit status,
        // reported first as below; one this kill ended reports only a signal,
        // which would mask the reader's diagnostic, so its status is dropped.
        let mut aligner = self.aligner.take().expect("aligner is Some until finalize");
        let aligner_result = if matches!(reader_result, Ok(Ok(()))) {
            aligner.wait()
        } else {
            let exited_on_its_own = aligner.has_exited();
            match (aligner.kill(), exited_on_its_own) {
                (true, true) => aligner.wait(),
                (true, false) => {
                    let _ = aligner.wait();
                    Ok(())
                }
                // Not confirmed gone: skip the unbounded `wait()`, as `Drop`
                // does; `AlignerProcess::drop` skips the stderr join too.
                (false, _) => Ok(()),
            }
        };

        let reader_result = match reader_result {
            Ok(result) => result,
            Err(_panic) => Err(io::Error::other("align-and-merge reader thread panicked")),
        };
        let writer_result = match writer_result {
            Ok(result) => result,
            Err(_panic) => Err(io::Error::other("align-and-merge writer thread panicked")),
        };
        if let Some(err) = rank_align_errors(aligner_result, reader_result, writer_result) {
            return Err(err);
        }

        // Fallback: surface any error recorded in `shared.error_slot` that
        // didn't come through a thread's `Result`. The writer returns `Ok`
        // after recording a cascade error when the reader stopped taking
        // tokens, so a cascade here means input was left unaligned; anything
        // else would be a future bug. Surface it rather than swallow it.
        if let Some(e) = self.shared.take_error() {
            return Err(e);
        }

        Ok(StepOutcome::Finished)
    }
}

/// The error `finalize` reports for the aligner's exit, the reader and the
/// writer, or `None` if all three succeeded.
///
/// A reader error is reported first, as `try_run` reports it when it arrives
/// before `finalize`: when the reader rejected the aligner's output (out of
/// input order, a missing template, a malformed record) it is the cause, and
/// when the aligner crashed it already says the aligner's output ended early.
/// The aligner's exit status and stderr tail follow it when the aligner also
/// failed. How the aligner died cannot tell these apart: an aligner that
/// ignores SIGPIPE exits with its own status when its stdout is closed. Then
/// the aligner's failure on its own, and last the writer's error (a broken
/// pipe included).
fn rank_align_errors(
    aligner: anyhow::Result<()>,
    reader: io::Result<()>,
    writer: io::Result<()>,
) -> Option<io::Error> {
    match (aligner, reader) {
        (Err(aligner), Err(reader)) => Some(io::Error::new(
            reader.kind(),
            format!("{reader}; the aligner exited: {aligner:#}"),
        )),
        (Ok(()), Err(reader)) => Some(reader),
        (Err(aligner), Ok(())) => {
            Some(io::Error::other(format!("aligner subprocess failed: {aligner:#}")))
        }
        (Ok(()), Ok(())) => writer.err(),
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
    /// here kills the aligner's whole process group (SIGTERM, then
    /// SIGKILL after a short grace), not just its `/bin/bash -c`
    /// wrapper, so every process holding the pipes exits, the reader
    /// thread's stdout `read` unblocks, and the daemon cascade completes.
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
            // forever — defeating kill()'s bound (the SIGTERM grace plus the
            // reap deadline, about 1.5 s). In that case we skip
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
    /// See [`SubprocessConfig::input_format`].
    pub(crate) input_format: UnalignedReadFormat,
    /// `--aligner::max-read-ahead`, or `None` for the default ceiling (see
    /// [`read_ahead_ceiling`]).
    pub(crate) max_read_ahead: Option<u64>,
    /// What the preset knows about running its aligner (budget floor, slow
    /// start-up, scratch directory).
    pub(crate) runtime: SubprocessRuntime,
}

/// The in-flight budget: the [`IN_FLIGHT_BUDGET_ENV`] override (`env`) when it
/// is set, else the larger of the `--aligner::chunk-size`-derived budget and
/// the preset's read-ahead floor (`floor`, `star`'s estimate of what it reads
/// before its first output). The override is a diagnostic knob and bypasses
/// the floor.
fn in_flight_budget(env: Option<&str>, chunk_size: u64, floor: u64) -> u64 {
    in_flight_budget_override(env)
        .unwrap_or_else(|| in_flight_budget_for_chunk_size(chunk_size).max(floor))
}

/// The warning for an explicit `--aligner::max-read-ahead` below the in-flight
/// budget: the gate starts at the ceiling then, so the budget itself is
/// lowered, not just the room to grow. The budget comes from
/// `--aligner::chunk-size` or a preset's read-ahead estimate (for `star`,
/// `--aligner::threads` and `--limitIObufferSize`), so the message names both.
fn max_read_ahead_lowers_budget_warning(budget: u64, requested: Option<u64>) -> Option<String> {
    requested.filter(|&ceiling| ceiling < budget).map(|ceiling| {
        format!(
            "--aligner::max-read-ahead {} is below the in-flight budget {} (set by \
             --aligner::chunk-size or the preset's read-ahead estimate), so the budget is \
             lowered to it and the gate cannot grow",
            bytesize::ByteSize(ceiling),
            bytesize::ByteSize(budget)
        )
    })
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
        crate::aligner::retain_freed_memory_unless_user_set("align");
        let budget = in_flight_budget(
            std::env::var(IN_FLIGHT_BUDGET_ENV).ok().as_deref(),
            self.chunk_size,
            self.runtime.min_in_flight_budget,
        );
        if let Some(dir) = &self.runtime.work_dir {
            log::info!("align: aligner scratch directory {}", dir.path().display());
        }
        let ceiling = read_ahead_ceiling(
            budget,
            self.max_read_ahead,
            crate::system::detect_total_memory() as u64,
        );
        if let Some(warning) = max_read_ahead_lowers_budget_warning(budget, self.max_read_ahead) {
            log::warn!("{warning}");
        }
        log::debug!(
            "align: in-flight budget {}, read-ahead ceiling {}",
            bytesize::ByteSize(budget),
            bytesize::ByteSize(ceiling)
        );
        let step = SubprocessAlignStep::new(
            SubprocessConfig {
                partial_output_header: Arc::clone(&ctx.partial_output_header),
                header_handle: ctx.header_handle.clone(),
                output_byte_limit: ctx.per_step_byte_limit,
                in_flight_unmapped_budget: budget,
                in_flight_read_ahead_ceiling: ceiling,
                in_flight_ceiling_timeout: ceiling_timeout_override(
                    std::env::var(IN_FLIGHT_CEILING_TIMEOUT_ENV).ok().as_deref(),
                )
                .unwrap_or(IN_FLIGHT_CEILING_TIMEOUT),
                accept_mid_pair_split: self.accept_mid_pair_split,
                passes_meth: crate::aligner::command_passes_meth(&self.command),
                input_format: self.input_format,
                startup_timeout: self.runtime.startup_timeout,
                require_complete_pairs: self.runtime.require_complete_pairs,
                work_dir: self.runtime.work_dir.clone(),
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
            refill: None,
        })
    }
}

#[cfg(test)]
mod tests {
    /// An explicit `--aligner::max-read-ahead` below the derived budget is
    /// warned about (it lowers the budget); at or above it, or unset, it is not.
    #[rstest::rstest]
    #[case::below(100, Some(50), true)]
    #[case::equal(100, Some(100), false)]
    #[case::above(100, Some(200), false)]
    #[case::unset(100, None, false)]
    fn max_read_ahead_below_the_budget_is_warned(
        #[case] budget: u64,
        #[case] requested: Option<u64>,
        #[case] warned: bool,
    ) {
        assert_eq!(max_read_ahead_lowers_budget_warning(budget, requested).is_some(), warned);
    }

    /// The budget is the larger of the chunk-size-derived budget and the
    /// preset's floor, unless the environment override is set, which wins
    /// outright (floor included).
    #[rstest::rstest]
    #[case::floor_above(None, 1, u64::MAX / 4, u64::MAX / 4)]
    #[case::floor_below(None, 1_000_000, 1, in_flight_budget_for_chunk_size(1_000_000))]
    #[case::floor_equal(
        None,
        1_000_000,
        in_flight_budget_for_chunk_size(1_000_000),
        in_flight_budget_for_chunk_size(1_000_000)
    )]
    #[case::no_floor(None, 1_000_000, 0, in_flight_budget_for_chunk_size(1_000_000))]
    #[case::override_bypasses_floor(Some("4096"), 1_000_000, u64::MAX / 4, 4096)]
    fn in_flight_budget_takes_the_floor_unless_overridden(
        #[case] env: Option<&str>,
        #[case] chunk_size: u64,
        #[case] floor: u64,
        #[case] expected: u64,
    ) {
        assert_eq!(in_flight_budget(env, chunk_size, floor), expected);
    }

    use super::*;
    use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED, REVERSE, SUPPLEMENTARY};

    #[rstest::rstest]
    #[case::before_deadline(9, None, false)]
    #[case::at_deadline(10, None, true)]
    #[case::after_deadline(11, None, true)]
    #[case::timely_output_late_dispatch(11, Some(9), false)]
    #[case::output_at_deadline(11, Some(10), true)]
    #[case::late_output_late_dispatch(12, Some(11), true)]
    fn startup_wait_deadline(
        #[case] elapsed_secs: u64,
        #[case] output_secs: Option<u64>,
        #[case] expires: bool,
    ) {
        let shared = SharedState::new();
        // Use seconds for a generous margin without sleeping: all observation
        // and spawn instants are synthetic, on the same monotonic clock.
        let started_at = std::time::Instant::now()
            .checked_sub(std::time::Duration::from_secs(elapsed_secs))
            .unwrap();
        if let Some(output_secs) = output_secs {
            shared.startup.lock().first_output_at =
                Some(started_at + std::time::Duration::from_secs(output_secs));
            shared.end_startup_wait();
        }
        let result = shared.claim_startup_wait(started_at, std::time::Duration::from_secs(10));
        assert_eq!(result.is_err(), expires);
        if let Err(error) = result {
            assert_eq!(error.kind(), io::ErrorKind::TimedOut);
            let startup = shared.startup.lock();
            assert!(startup.over);
            assert!(startup.claim.is_none(), "expiry releases the watchdog claim");
        }
    }

    /// `finalize` reports the reader's error first (with the aligner's exit
    /// after it when the aligner also failed), else the aligner's failure, and
    /// the writer's error only when nothing else failed.
    #[rstest::rstest]
    #[case::reader_and_aligner(
        Some("exited with status 1"),
        Some("bad order"),
        Some(Pipe::Broken),
        Some("bad order; the aligner exited: exited with status 1")
    )]
    #[case::reader_only(None, Some("bad order"), Some(Pipe::Broken), Some("bad order"))]
    #[case::aligner_only(
        Some("exited with status 1"),
        None,
        Some(Pipe::Broken),
        Some("aligner subprocess failed: exited with status 1")
    )]
    #[case::writer_broken_pipe_only(None, None, Some(Pipe::Broken), Some("broken pipe"))]
    #[case::writer_only(None, None, Some(Pipe::Other), Some("writer failed"))]
    #[case::all_ok(None, None, None, None)]
    fn rank_align_errors_puts_the_cause_first(
        #[case] aligner: Option<&str>,
        #[case] reader: Option<&str>,
        #[case] writer: Option<Pipe>,
        #[case] expected: Option<&str>,
    ) {
        let aligner = aligner.map_or(Ok(()), |msg| Err(anyhow!("{msg}")));
        let reader =
            reader.map_or(Ok(()), |msg| Err(io::Error::new(io::ErrorKind::InvalidData, msg)));
        let writer = writer.map_or(Ok(()), |pipe| match pipe {
            Pipe::Broken => Err(io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe")),
            Pipe::Other => Err(io::Error::other("writer failed")),
        });
        let reader_failed = reader.is_err();
        let ranked = rank_align_errors(aligner, reader, writer);
        if reader_failed {
            let kind = ranked.as_ref().map(io::Error::kind);
            assert_eq!(kind, Some(io::ErrorKind::InvalidData), "the reader's kind is kept");
        }
        assert_eq!(ranked.map(|e| e.to_string()).as_deref(), expected);
    }

    /// How the stdin writer failed.
    #[derive(Debug, Clone, Copy)]
    enum Pipe {
        Broken,
        Other,
    }

    /// A cascade error recorded before the reader's error does not mask it: the
    /// root cause wins whichever thread records first, and the cascade is
    /// reported only when nothing else is.
    #[rstest::rstest]
    #[case::cascade_first(true)]
    #[case::cause_first(false)]
    fn shared_state_prefers_the_root_cause(#[case] cascade_first: bool) {
        let shared = SharedState::new();
        let cascade = || io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe");
        if cascade_first {
            shared.record_cascade_error(cascade());
            assert!(shared.has_cascade(), "the cascade is recorded");
            shared.record_error(io::Error::other("queryname mismatch"));
        } else {
            shared.record_error(io::Error::other("queryname mismatch"));
            shared.record_cascade_error(cascade());
        }
        assert_eq!(shared.take_error().expect("an error").to_string(), "queryname mismatch");

        let only_cascade = SharedState::new();
        only_cascade.record_cascade_error(cascade());
        assert!(only_cascade.take_root_cause().is_none(), "a cascade is not a root cause");
        assert_eq!(only_cascade.take_error().expect("the cascade").to_string(), "broken pipe");
    }

    /// Only cascade errors expire; a root cause wins even at the deadline.
    #[rstest::rstest]
    #[case::healthy(false, None, 10_000, None)]
    #[case::root_only(false, Some("reader failed"), 0, Some(("reader failed", io::ErrorKind::InvalidData)))]
    #[case::before_deadline(true, None, 999, None)]
    #[case::at_deadline(true, None, 1000, Some(("broken pipe", io::ErrorKind::BrokenPipe)))]
    #[case::after_deadline(true, None, 1001, Some(("broken pipe", io::ErrorKind::BrokenPipe)))]
    #[case::reader_before_deadline(true, Some("reader failed"), 999, Some(("reader failed", io::ErrorKind::InvalidData)))]
    #[case::reader_at_deadline(true, Some("reader failed"), 1000, Some(("reader failed", io::ErrorKind::InvalidData)))]
    #[case::reader_after_deadline(true, Some("reader failed"), 1001, Some(("reader failed", io::ErrorKind::InvalidData)))]
    fn shared_state_bounds_only_cascade_errors(
        #[case] cascade: bool,
        #[case] root: Option<&str>,
        #[case] elapsed_ms: u64,
        #[case] expected: Option<(&str, io::ErrorKind)>,
    ) {
        let shared = SharedState::new();
        if cascade {
            shared.record_cascade_error(io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe"));
        }
        let recorded_at = shared
            .error_slot
            .lock()
            .cascade
            .as_ref()
            .map_or_else(Instant::now, |(recorded_at, _)| *recorded_at);
        if let Some(msg) = root {
            shared.record_error(io::Error::new(io::ErrorKind::InvalidData, msg));
        }
        let error = shared.take_dispatch_error(recorded_at + Duration::from_millis(elapsed_ms));
        assert_eq!(
            error.map(|error| (error.to_string(), error.kind())),
            expected.map(|(msg, kind)| (msg.to_string(), kind))
        );
    }

    /// A writer broken pipe is recorded as a cascade error and keeps its kind;
    /// any other writer failure is a root cause.
    #[rstest::rstest]
    #[case::broken_pipe(true, io::ErrorKind::BrokenPipe)]
    #[case::other(false, io::ErrorKind::Other)]
    fn record_writer_error_routes_a_broken_pipe_to_the_cascade(
        #[case] broken_pipe: bool,
        #[case] expected_kind: io::ErrorKind,
    ) {
        let shared = SharedState::new();
        let returned = record_writer_error(&shared, broken_pipe, "writer failed".to_string());
        assert_eq!(returned.kind(), expected_kind);
        assert_eq!(shared.has_cascade(), broken_pipe);
        assert_eq!(shared.take_root_cause().is_some(), !broken_pipe);
    }

    /// The writer's broken-pipe check sees through `anyhow` context.
    #[test]
    fn is_broken_pipe_finds_a_wrapped_io_error() {
        let wrapped = anyhow::Error::from(io::Error::new(io::ErrorKind::BrokenPipe, "pipe"))
            .context("writing a record");
        assert!(is_broken_pipe(&wrapped));
        assert!(!is_broken_pipe(&anyhow::anyhow!("no primary records")));
    }

    /// What the slot holds when the writer thread is found gone.
    #[derive(Debug, Clone, Copy)]
    enum Recorded {
        RootCause,
        Cascade,
        Nothing,
    }

    /// `writer_gone` returns a recorded root cause at once, drains on a cascade
    /// alone (closing the aligner's stdin and dropping held input), and fails
    /// with a generic error when nothing was recorded.
    #[rstest::rstest]
    #[case::root_cause(Recorded::RootCause, Err("queryname mismatch"))]
    #[case::cascade_only(Recorded::Cascade, Ok(()))]
    #[case::nothing(Recorded::Nothing, Err("writer thread is gone"))]
    fn writer_gone_drains_only_on_a_cascade(
        #[case] recorded: Recorded,
        #[case] expected: Result<(), &str>,
    ) {
        let mut step = SubprocessAlignStep::new(make_test_cfg(), "cat").expect("spawn cat");
        step.held_in.put(BamTemplateBatch::new(0, vec![make_paired_primary_template(b"q")]));
        match recorded {
            Recorded::RootCause => step.shared.record_error(io::Error::other("queryname mismatch")),
            Recorded::Cascade => step
                .shared
                .record_cascade_error(io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe")),
            Recorded::Nothing => {}
        }
        match (step.writer_gone(), expected) {
            (Ok(outcome), Ok(())) => {
                assert!(matches!(outcome, StepOutcome::Progress), "drains: {outcome:?}");
                assert!(step.in_tx.is_none(), "the aligner's stdin is closed");
                assert!(!step.held_in.is_held(), "held input is dropped");
                assert!(step.shared.has_cascade(), "the cascade stays for the bounded drain");
            }
            (Err(err), Err(needle)) => {
                assert!(err.to_string().contains(needle), "expected {needle:?}, got: {err}");
            }
            (got, want) => panic!("expected {want:?}, got {got:?}"),
        }
    }

    /// A published root cause or expired cascade wins over a startup timeout.
    #[rstest::rstest]
    #[case::root_cause(false, io::ErrorKind::InvalidData, "reader rejected output")]
    #[case::expired_cascade(true, io::ErrorKind::BrokenPipe, "broken pipe")]
    fn dispatch_error_precedes_startup_timeout(
        #[case] cascade: bool,
        #[case] expected_kind: io::ErrorKind,
        #[case] expected_message: &str,
    ) {
        use crate::pipeline::core::testing::StepProbe;

        let mut cfg = make_test_cfg();
        cfg.startup_timeout = Some(Duration::ZERO);
        let mut step = SubprocessAlignStep::new(cfg, "exec sleep 8").expect("spawn sleep");
        if cascade {
            step.shared
                .record_cascade_error(io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe"));
            step.shared.error_slot.lock().cascade.as_mut().expect("cascade").0 =
                Instant::now().checked_sub(CASCADE_DRAIN_TIMEOUT).expect("backdate cascade");
        } else {
            step.shared
                .record_error(io::Error::new(io::ErrorKind::InvalidData, "reader rejected output"));
        }
        let probe = StepProbe::new(&step);
        let error = probe.try_run(&mut step).expect_err("the recorded failure must win");
        assert_eq!(error.kind(), expected_kind);
        assert_eq!(error.to_string(), expected_message);
    }

    /// An expired cascade must beat a held-output retry that downstream refuses.
    #[test]
    fn cascade_deadline_precedes_output_backpressure() {
        use crate::pipeline::core::testing::StepProbe;

        let mut cfg = make_test_cfg();
        cfg.output_byte_limit = 1;
        let mut step = SubprocessAlignStep::new(cfg, "exec sleep 8").expect("spawn sleep");
        let (out_tx, out_rx) = sync_channel(OUT_CHAN_DEPTH);
        step.out_rx = Some(out_rx);
        for serial in 0..2 {
            out_tx
                .send(ZipperBatch {
                    serial,
                    mapped: Vec::new(),
                    unmapped: BamTemplateBatch::new(
                        serial,
                        vec![make_paired_primary_template(b"q")],
                    ),
                })
                .expect("queue output");
        }
        drop(out_tx);
        let probe = StepProbe::new(&step);
        assert_eq!(probe.try_run(&mut step).expect("fill downstream"), StepOutcome::Progress);
        assert!(step.held_out.is_held(), "the second batch must be held");
        assert_eq!(probe.try_run(&mut step).expect("refused retry"), StepOutcome::NoProgress);

        step.shared.record_cascade_error(io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe"));
        step.shared.error_slot.lock().cascade.as_mut().expect("cascade").0 =
            Instant::now().checked_sub(CASCADE_DRAIN_TIMEOUT).expect("backdate cascade");
        let error = probe.try_run(&mut step).expect_err("deadline must beat backpressure");
        assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(error.to_string(), "broken pipe");
    }

    /// An exited aligner's status and stderr must beat a cascade, even after
    /// its deadline, when the writer failed before sending any batch token.
    #[rstest::rstest]
    #[case::before_deadline(0)]
    #[case::after_deadline(1100)]
    fn cascade_preserves_an_exited_aligners_failure(
        #[case] elapsed_ms: u64,
        #[values(UnalignedReadFormat::Fastq, UnalignedReadFormat::Sam)]
        input_format: UnalignedReadFormat,
    ) {
        use crate::pipeline::core::testing::StepProbe;

        let mut cfg = make_test_cfg();
        cfg.input_format = input_format;
        cfg.require_complete_pairs = true;
        let mut step = SubprocessAlignStep::new(
            cfg,
            "exec 0<&-; printf '@HD\\tVN:1.6\\n'; printf 'startup-failed\\n' >&2; exit 42",
        )
        .expect("spawn failing aligner");
        let deadline = Instant::now() + Duration::from_secs(2);
        while !step.shared.startup.lock().over {
            assert!(Instant::now() < deadline, "the aligner must close stdin and write its header");
            thread::sleep(Duration::from_millis(1));
        }
        while !step.aligner.as_mut().expect("aligner").is_finished().expect("completion probe") {
            assert!(Instant::now() < deadline, "the failing aligner and stderr relay must finish");
            thread::sleep(Duration::from_millis(1));
        }
        // Force a write before token publication rather than fitting into the
        // buffer and letting a token turn the reader's clean EOF into an error.
        let sequence = vec![b'A'; FASTQ_WRITER_BUF_BYTES];
        let qualities = vec![b'I'; FASTQ_WRITER_BUF_BYTES];
        let records = [FIRST_SEGMENT, LAST_SEGMENT]
            .into_iter()
            .map(|mate| {
                let mut builder = fgumi_raw_bam::SamBuilder::new();
                builder
                    .read_name(b"early_exit")
                    .flags(PAIRED | mate)
                    .sequence(&sequence)
                    .qualities(&qualities);
                builder.build()
            })
            .collect();
        let template = Template::from_records(records).expect("complete primary pair");
        step.in_tx
            .take()
            .expect("writer input")
            .send(BamTemplateBatch::new(0, vec![template]))
            .expect("queue input and close it");
        while !step.reader_thread.as_ref().expect("reader").is_finished()
            || !step.writer_thread.as_ref().expect("writer").is_finished()
        {
            assert!(Instant::now() < deadline, "I/O threads must finish after the failed write");
            thread::sleep(Duration::from_millis(1));
        }
        assert!(step.shared.take_root_cause().is_none(), "no token reached the reader");
        step.shared.error_slot.lock().cascade.as_mut().expect("writer cascade").0 = Instant::now()
            .checked_sub(Duration::from_millis(elapsed_ms))
            .expect("backdate cascade");
        let probe = StepProbe::new(&step);
        probe.close_input();
        let error = probe.try_run(&mut step).expect_err("dispatch must report the exited child");
        let message = error.to_string();
        assert!(
            message.contains("aligner subprocess failed: aligner process exited with status 42."),
            "{message}"
        );
        assert!(message.contains("startup-failed"), "the stderr tail must survive: {message}");
    }

    /// A second worker may enter after the Serial guard is released but before
    /// the first worker publishes its error to the pipeline signal.
    #[rstest::rstest]
    #[case::expired_cascade(false, io::ErrorKind::BrokenPipe, "Broken pipe")]
    #[case::reader_root(true, io::ErrorKind::InvalidData, "reader rejected output")]
    fn dispatch_failure_survives_redispatch(
        #[case] reader_root: bool,
        #[case] expected_kind: io::ErrorKind,
        #[case] expected_message: &str,
        #[values(UnalignedReadFormat::Fastq, UnalignedReadFormat::Sam)]
        input_format: UnalignedReadFormat,
    ) {
        use crate::pipeline::core::testing::StepProbe;

        let mut cfg = make_test_cfg();
        cfg.input_format = input_format;
        cfg.require_complete_pairs = true;
        let mut step = SubprocessAlignStep::new(
            cfg,
            "exec 0<&-; printf '@HD\\tVN:1.6\\n'; exec 1>&-; exec sleep 4",
        )
        .expect("spawn live closed-pipe aligner");
        let deadline = Instant::now() + Duration::from_secs(2);
        while !step.shared.startup.lock().over {
            assert!(Instant::now() < deadline, "the aligner must close stdin and write its header");
            thread::sleep(Duration::from_millis(1));
        }
        let sequence = vec![b'A'; FASTQ_WRITER_BUF_BYTES];
        let qualities = vec![b'I'; FASTQ_WRITER_BUF_BYTES];
        let records = [FIRST_SEGMENT, LAST_SEGMENT]
            .into_iter()
            .map(|mate| {
                let mut builder = fgumi_raw_bam::SamBuilder::new();
                builder
                    .read_name(b"closed_pipes")
                    .flags(PAIRED | mate)
                    .sequence(&sequence)
                    .qualities(&qualities);
                builder.build()
            })
            .collect();
        step.in_tx
            .take()
            .expect("writer input")
            .send(BamTemplateBatch::new(
                0,
                vec![Template::from_records(records).expect("complete primary pair")],
            ))
            .expect("queue input and close it");
        while !step.reader_thread.as_ref().expect("reader").is_finished()
            || !step.writer_thread.as_ref().expect("writer").is_finished()
        {
            assert!(Instant::now() < deadline, "I/O threads must finish while the child is live");
            thread::sleep(Duration::from_millis(1));
        }
        assert!(step.shared.take_root_cause().is_none(), "the real write failed before a token");
        step.shared.error_slot.lock().cascade.as_mut().expect("writer cascade").0 =
            Instant::now().checked_sub(CASCADE_DRAIN_TIMEOUT).expect("backdate cascade");
        if reader_root {
            step.shared
                .record_error(io::Error::new(io::ErrorKind::InvalidData, "reader rejected output"));
        }
        let probe = StepProbe::new(&step);
        probe.close_input();
        let first = probe.try_run(&mut step).expect_err("first worker observes the failure");
        assert_eq!(first.kind(), expected_kind);
        let first_message = first.to_string();
        assert!(first_message.contains(expected_message), "{first_message}");
        for _ in 0..3 {
            let start = Instant::now();
            let error = probe.try_run(&mut step).expect_err("redispatch must stay failed");
            assert!(
                start.elapsed() < Duration::from_secs(1),
                "redispatch must not wait on the child"
            );
            assert_eq!(error.kind(), expected_kind);
            assert_eq!(
                error.to_string(),
                first_message,
                "the original failure must win publication"
            );
            assert!(step.aligner.is_some(), "a failed step must not enter finalize");
        }
    }

    #[rstest::rstest]
    #[case::healthy(0, None)]
    #[case::aligner_failure(
        42,
        Some("aligner subprocess failed: aligner process exited with status 42.")
    )]
    fn completed_drain_reentry_preserves_its_outcome(
        #[case] exit_code: u8,
        #[case] expected_error: Option<&str>,
    ) {
        use crate::pipeline::core::testing::StepProbe;

        let mut step = SubprocessAlignStep::new(
            make_test_cfg(),
            &format!("printf '@HD\\tVN:1.6\\n'; exit {exit_code}"),
        )
        .expect("spawn header-only aligner");
        drop(step.in_tx.take());
        let deadline = Instant::now() + Duration::from_secs(2);
        while !step.reader_thread.as_ref().expect("reader").is_finished()
            || !step.writer_thread.as_ref().expect("writer").is_finished()
            || !step.aligner.as_mut().expect("aligner").is_finished().expect("completion probe")
        {
            assert!(Instant::now() < deadline, "the aligner and all I/O threads must finish");
            thread::sleep(Duration::from_millis(1));
        }
        let probe = StepProbe::new(&step);
        probe.close_input();
        for _ in 0..3 {
            let result = probe.try_run(&mut step);
            if let Some(expected_error) = expected_error {
                let error = result.expect_err("a failed drain must stay failed");
                assert_eq!(error.kind(), io::ErrorKind::Other);
                assert!(error.to_string().starts_with(expected_error), "{error}");
            } else {
                assert_eq!(result.expect("successful drain"), StepOutcome::Finished);
            }
        }
    }

    /// A successful reader can finish while the child stays alive. A fresh
    /// cascade must not enter finalize's wait; an expired one must fail dispatch.
    #[rstest::rstest]
    #[case::before_deadline(false)]
    #[case::at_deadline(true)]
    fn cascade_drain_does_not_finalize_a_live_child(#[case] expired: bool) {
        use crate::pipeline::core::testing::StepProbe;

        let mut step = SubprocessAlignStep::new(
            make_test_cfg(),
            "printf '@HD\\tVN:1.6\\n'; exec 1>&-; exec sleep 8",
        )
        .expect("spawn header-only aligner");
        drop(step.in_tx.take());
        let deadline = Instant::now() + Duration::from_secs(2);
        while !step.reader_thread.as_ref().expect("reader").is_finished()
            || !step.writer_thread.as_ref().expect("writer").is_finished()
        {
            assert!(Instant::now() < deadline, "I/O threads must finish before the child exits");
            thread::sleep(Duration::from_millis(1));
        }
        assert!(step.shared.take_error().is_none(), "both I/O threads must succeed");
        step.shared.record_cascade_error(io::Error::new(io::ErrorKind::BrokenPipe, "broken pipe"));
        if expired {
            step.shared.error_slot.lock().cascade.as_mut().expect("cascade").0 =
                Instant::now().checked_sub(CASCADE_DRAIN_TIMEOUT).expect("backdate cascade");
        }
        let probe = StepProbe::new(&step);
        probe.close_input();
        let start = Instant::now();
        let result = probe.try_run(&mut step);
        assert!(start.elapsed() < Duration::from_secs(1), "dispatch must not wait on the child");
        if expired {
            let error = result.expect_err("expired cascade must fail");
            assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
            assert_eq!(error.to_string(), "broken pipe");
        } else {
            assert_eq!(result.expect("bounded drain"), StepOutcome::NoProgress);
        }
    }

    /// A flush before gate acquisition keeps a broken pipe as a cascade, so a
    /// reader error published afterward still wins. FASTQ fails flushing batch
    /// one before batch two; SAM fails flushing its buffered header before batch one.
    #[rstest::rstest]
    #[case::fastq(UnalignedReadFormat::Fastq, 1)]
    #[case::sam(UnalignedReadFormat::Sam, 0)]
    fn writer_pre_batch_flush_keeps_a_broken_pipe_as_a_cascade(
        #[case] input_format: UnalignedReadFormat,
        #[case] expected_tokens: usize,
        #[values(false, true)] require_complete_pairs: bool,
    ) {
        use std::io::BufRead;

        let mut aligner = AlignerProcess::spawn(
            "exec 0<&-; printf 'ready\\n'; exec sleep 8",
            ALIGNER_STDERR_RING_SIZE,
        )
        .expect("spawn closed-stdin aligner");
        let mut stdout = BufReader::new(aligner.take_stdout().expect("stdout"));
        let mut ready = String::new();
        stdout.read_line(&mut ready).expect("wait for stdin closure");
        assert_eq!(ready, "ready\n");

        let (in_tx, in_rx) = sync_channel(2);
        for serial in 0..2 {
            in_tx
                .send(BamTemplateBatch::new(serial, vec![make_paired_primary_template(b"q")]))
                .expect("queue batch");
        }
        drop(in_tx);
        let (token_tx, token_rx) = channel();
        let shared = Arc::new(SharedState::new());
        let gate =
            InFlightGate::new(u64::MAX, u64::MAX, IN_FLIGHT_STALL_GRACE, IN_FLIGHT_CEILING_TIMEOUT);
        let result = writer_loop(
            in_rx,
            token_tx,
            aligner.take_stdin().expect("stdin"),
            Arc::clone(&shared),
            &gate,
            input_format,
            require_complete_pairs,
        );
        let cascade_recorded = shared.has_cascade();
        shared.record_error(io::Error::new(io::ErrorKind::InvalidData, "reader rejected output"));
        let reported = shared.take_error().expect("recorded failure");
        drop(stdout);
        drop(aligner);

        let writer_error = result.expect_err("the pre-batch flush must fail");
        assert_eq!(writer_error.kind(), io::ErrorKind::BrokenPipe);
        assert!(writer_error.to_string().starts_with("align-and-merge writer: flush:"));
        assert_eq!(token_rx.try_iter().count(), expected_tokens, "the intended flush must fail");
        assert!(cascade_recorded, "the broken pipe must not occupy the root slot");
        assert_eq!(reported.kind(), io::ErrorKind::InvalidData);
        assert_eq!(reported.to_string(), "reader rejected output");
    }

    /// When the reader has stopped taking tokens, the writer exits `Ok` but
    /// records a cascade error, so the step waits for the reader's own error
    /// and the drain deadline still reports unaligned input if there is none.
    #[test]
    fn writer_loop_records_a_cascade_when_the_reader_is_gone() {
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;

        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn cat");
        let cat_stdin = cat.stdin.take().expect("cat stdin");

        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(1);
        let (token_tx, token_rx) = channel::<BatchToken>();
        drop(token_rx);
        let shared = Arc::new(SharedState::new());
        in_tx
            .send(BamTemplateBatch::new(0, vec![make_paired_primary_template(b"q")]))
            .expect("send batch");
        drop(in_tx);

        let gate =
            InFlightGate::new(u64::MAX, u64::MAX, IN_FLIGHT_STALL_GRACE, IN_FLIGHT_CEILING_TIMEOUT);
        let result = writer_loop(
            in_rx,
            token_tx,
            cat_stdin,
            Arc::clone(&shared),
            &gate,
            UnalignedReadFormat::Fastq,
            false,
        );

        assert!(result.is_ok(), "the reader's error, not the writer's, is the cause: {result:?}");
        assert!(shared.take_root_cause().is_none(), "no root cause is recorded");
        let cascade = shared.take_error().expect("a cascade is recorded");
        assert!(cascade.to_string().contains("the reader stopped"), "got: {cascade}");

        let _ = cat.kill();
        let _ = cat.wait();
    }

    fn make_test_cfg() -> SubprocessConfig {
        SubprocessConfig {
            partial_output_header: Arc::new(Header::default()),
            header_handle: HeaderHandle::new(),
            output_byte_limit: 1024 * 1024,
            in_flight_unmapped_budget: in_flight_budget_for_chunk_size(150_000_000),
            in_flight_read_ahead_ceiling: u64::MAX,
            in_flight_ceiling_timeout: IN_FLIGHT_CEILING_TIMEOUT,
            accept_mid_pair_split: true,
            passes_meth: false,
            input_format: UnalignedReadFormat::Fastq,
            startup_timeout: None,
            require_complete_pairs: false,
            work_dir: None,
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
        // AlignerProcess::kill waits up to the SIGTERM grace plus the reap
        // deadline (about 1.5 s, or about 2.5 s if a group member outlives the
        // SIGKILL); thread joins should
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

    #[rstest::rstest]
    #[case::star_unpaired(true, &[0], false)]
    #[case::star_missing_r1(true, &[PAIRED | LAST_SEGMENT], false)]
    #[case::star_missing_r2(true, &[PAIRED | FIRST_SEGMENT], false)]
    #[case::star_secondary_mate(true, &[PAIRED | FIRST_SEGMENT, PAIRED | LAST_SEGMENT | SUPPLEMENTARY], false)]
    #[case::star_unpaired_r1(true, &[0, PAIRED | LAST_SEGMENT], false)]
    #[case::star_pair(true, &[PAIRED | FIRST_SEGMENT, PAIRED | LAST_SEGMENT], true)]
    #[case::command_single(false, &[0], true)]
    #[case::command_orphan(false, &[PAIRED | FIRST_SEGMENT], true)]
    fn writer_requires_complete_pairs_only_when_configured(
        #[case] require_complete_pairs: bool,
        #[case] flags: &[u16],
        #[case] accepted: bool,
    ) {
        use std::process::{Command, Stdio};
        let dir = tempfile::tempdir().unwrap();
        let captured = dir.path().join("stdin.sam");
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(std::fs::File::create(&captured).unwrap())
            .spawn()
            .unwrap();
        let records = flags.iter().map(|&f| make_record(b"pair_check", f)).collect();
        let template = Template::from_records(records).unwrap();
        let (in_tx, in_rx) = sync_channel(1);
        let (token_tx, token_rx) = channel();
        in_tx.send(BamTemplateBatch::new(0, vec![template])).unwrap();
        drop(in_tx);
        let shared = Arc::new(SharedState::new());
        let gate =
            InFlightGate::new(u64::MAX, u64::MAX, IN_FLIGHT_STALL_GRACE, IN_FLIGHT_CEILING_TIMEOUT);
        let result = writer_loop(
            in_rx,
            token_tx,
            cat.stdin.take().unwrap(),
            Arc::clone(&shared),
            &gate,
            UnalignedReadFormat::Sam,
            require_complete_pairs,
        );
        assert!(cat.wait().unwrap().success());
        assert_eq!(result.is_ok(), accepted);
        let captured = std::fs::read_to_string(captured).unwrap();
        if accepted {
            assert_eq!(token_rx.recv().unwrap().n_templates, 1);
            assert!(captured.contains("pair_check"));
        } else {
            let error = result.unwrap_err().to_string();
            assert!(
                error.contains("complete primary pair") && error.contains("pair_check"),
                "{error}"
            );
            assert_eq!(shared.take_error().unwrap().to_string(), error);
            assert_eq!(
                captured,
                std::str::from_utf8(UNALIGNED_SAM_HEADER).unwrap(),
                "no record of the incomplete template may reach STAR"
            );
            assert!(token_rx.recv().is_err());
        }
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
        let gate =
            InFlightGate::new(u64::MAX, u64::MAX, IN_FLIGHT_STALL_GRACE, IN_FLIGHT_CEILING_TIMEOUT);
        let result = writer_loop(
            in_rx,
            token_tx,
            cat_stdin,
            Arc::clone(&shared),
            &gate,
            UnalignedReadFormat::Fastq,
            false,
        );

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

    /// The bytes `writer_loop` writes to the aligner's stdin for one template
    /// built from records given out of order (R2, a supplementary R1, then R1),
    /// in each `--aligner::input-format`: R1 before R2, the supplementary record
    /// skipped, no mate suffix, and for SAM the `@HD` line first.
    #[rstest::rstest]
    #[case::fastq(UnalignedReadFormat::Fastq, "@q1\nACGT\n+\nIIII\n@q1\nTTGC\n+\n5555\n")]
    #[case::sam(
        UnalignedReadFormat::Sam,
        "@HD\tVN:1.6\tSO:unsorted\n\
         q1\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
         q1\t141\t*\t0\t0\t*\t*\t0\t0\tTTGC\t5555\n"
    )]
    fn writer_loop_writes_primaries_in_input_format(
        #[case] input_format: UnalignedReadFormat,
        #[case] expected: &str,
    ) {
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;

        let dir = tempfile::tempdir().expect("tempdir");
        let out_path = dir.path().join("stdin.txt");
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(std::fs::File::create(&out_path).expect("create capture file"))
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn cat");
        let cat_stdin = cat.stdin.take().expect("cat stdin");

        let r2 = {
            let mut b = fgumi_raw_bam::SamBuilder::new();
            b.read_name(b"q1")
                .flags(PAIRED | LAST_SEGMENT)
                .sequence(b"TTGC")
                .qualities(&[20, 20, 20, 20]);
            b.build()
        };
        let supplementary = make_record(b"q1", PAIRED | FIRST_SEGMENT | SUPPLEMENTARY);
        let r1 = {
            let mut b = fgumi_raw_bam::SamBuilder::new();
            b.read_name(b"q1")
                .flags(PAIRED | FIRST_SEGMENT)
                .sequence(b"ACGT")
                .qualities(&[40, 40, 40, 40]);
            b.build()
        };
        let template = Template::from_records(vec![r2, supplementary, r1]).expect("template");

        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(1);
        let (token_tx, token_rx) = channel::<BatchToken>();
        in_tx.send(BamTemplateBatch::new(0, vec![template])).expect("send batch");
        drop(in_tx);
        let gate =
            InFlightGate::new(u64::MAX, u64::MAX, IN_FLIGHT_STALL_GRACE, IN_FLIGHT_CEILING_TIMEOUT);
        writer_loop(
            in_rx,
            token_tx,
            cat_stdin,
            Arc::new(SharedState::new()),
            &gate,
            input_format,
            false,
        )
        .expect("writer_loop");
        cat.wait().expect("cat exits");

        let token = token_rx.recv().expect("one token");
        assert_eq!(token.n_templates, 1);
        assert_eq!(std::fs::read_to_string(&out_path).expect("read capture"), expected);
    }

    /// Before waiting on the gate the writer flushes what it has buffered, also
    /// once the gate is stalled and the next batch would pass the ceiling: the
    /// aligner may need that tail to finish a chunk, and holding it back would
    /// fail the run at the ceiling. Batch A fills an empty gate; B waits out
    /// the grace and grows the gate to exactly A + B (stalled); C would pass
    /// that ceiling, so the writer waits again, and fails once no output has
    /// been consumed for the ceiling timeout. During that wait B must already
    /// have reached the aligner, and the failure is recorded for `try_run`.
    #[test]
    fn writer_loop_flushes_before_waiting_at_the_ceiling() {
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;
        use std::time::{Duration, Instant};

        let dir = tempfile::tempdir().expect("tempdir");
        let captured = dir.path().join("stdin.fq");
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(std::fs::File::create(&captured).expect("create capture file"))
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn cat");
        let cat_stdin = cat.stdin.take().expect("cat stdin");

        let batches: Vec<BamTemplateBatch> = ["read_A", "read_B", "read_C"]
            .iter()
            .enumerate()
            .map(|(i, name)| {
                BamTemplateBatch::new(i as u64, vec![make_paired_primary_template(name.as_bytes())])
            })
            .collect();
        let ceiling = batches[0].heap_size() as u64 + batches[1].heap_size() as u64;
        let grace = Duration::from_secs(2);
        let gate = Arc::new(InFlightGate::new(1, ceiling, grace, grace + Duration::from_secs(1)));

        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(4);
        let (token_tx, _token_rx) = channel::<BatchToken>();
        for batch in batches {
            in_tx.send(batch).expect("queue batch");
        }
        drop(in_tx);
        let writer_gate = Arc::clone(&gate);
        let shared = Arc::new(SharedState::new());
        let writer_shared = Arc::clone(&shared);
        let writer = std::thread::spawn(move || {
            writer_loop(
                in_rx,
                token_tx,
                cat_stdin,
                writer_shared,
                &writer_gate,
                UnalignedReadFormat::Fastq,
                false,
            )
        });

        // B's stall grows the gate to the ceiling; C then waits on it.
        let deadline = Instant::now() + Duration::from_secs(30);
        while gate.budget() != ceiling {
            assert!(Instant::now() < deadline, "the gate never grew for batch B");
            std::thread::sleep(Duration::from_millis(10));
        }
        std::thread::sleep(grace / 4);
        let seen = std::fs::read_to_string(&captured).expect("read capture");
        assert!(
            seen.contains("@read_B"),
            "batch B must reach the aligner before the writer waits at the ceiling; saw:\n{seen}"
        );

        let err = writer.join().expect("writer thread").expect_err("C passes the ceiling");
        assert!(err.to_string().contains("read ahead more than"), "unexpected: {err}");
        let recorded = shared.take_error().expect("the ceiling error is recorded for try_run");
        assert!(recorded.to_string().contains("read ahead more than"), "recorded: {recorded}");
        cat.wait().expect("cat exits");
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
        let gate = Arc::new(InFlightGate::new(
            u64::MAX,
            u64::MAX,
            IN_FLIGHT_STALL_GRACE,
            IN_FLIGHT_CEILING_TIMEOUT,
        ));
        let gate_for_worker = Arc::clone(&gate);
        let writer_handle = std::thread::Builder::new()
            .name("writer-loop-deadlock-test".into())
            .spawn(move || {
                writer_loop(
                    in_rx,
                    token_tx,
                    cat_stdin,
                    shared_for_worker,
                    &gate_for_worker,
                    UnalignedReadFormat::Fastq,
                    false,
                )
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
        assert!(result.is_ok(), "SAM-text header should resolve handle to Ok: {:?}", result.err());
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

    /// Encode a header and one same-queryname record per independently readable
    /// chunk, without EOF, so a pipe producer can pace the reader's progress.
    fn same_queryname_output_chunks(sam: bool, record_count: usize) -> Vec<Vec<u8>> {
        use fgumi_bgzf::InlineBgzfCompressor;
        use fgumi_raw_bam::flags::UNMAPPED;

        // Separate BGZF blocks let the real BAM reader consume each record
        // immediately, just as newline-delimited SAM does.
        let encode = |bytes: Vec<u8>| {
            if sam {
                bytes
            } else {
                let mut c = InlineBgzfCompressor::new(1);
                c.write_all(&bytes).unwrap();
                c.flush().unwrap();
                let mut compressed = Vec::new();
                c.write_blocks_to(&mut compressed).unwrap();
                compressed
            }
        };
        let mut header = Vec::new();
        if sam {
            header.extend_from_slice(b"@HD\tVN:1.6\n");
        } else {
            fgumi_bam_io::write_bam_header(&mut header, &Header::default()).unwrap();
        }
        let mut chunks = vec![encode(header)];
        for i in 0..record_count {
            let flags = UNMAPPED | if i == 0 { 0 } else { SUPPLEMENTARY };
            let mut bytes = Vec::new();
            if sam {
                writeln!(bytes, "large\t{flags}\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").unwrap();
            } else {
                fgumi_raw_bam::write_raw_record(&mut bytes, &make_record(b"large", flags)).unwrap();
            }
            chunks.push(encode(bytes));
        }
        chunks
    }

    /// A large same-queryname group must keep the ceiling waiter alive before
    /// even one template completes. Removing either format's record-progress
    /// signal must fail the active-output case; EOF is deliberately held back
    /// in the stopped-output case so only the quiet timeout can end the wait.
    #[rstest::rstest]
    #[case::bam_active_output_outlasts_the_timeout(false, 30, false)]
    #[case::sam_active_output_outlasts_the_timeout(true, 30, false)]
    #[case::bam_output_stops_before_the_template_finishes(false, 30, true)]
    #[case::sam_output_stops_before_the_template_finishes(true, 30, true)]
    fn reader_loop_reports_record_progress_at_the_ceiling(
        #[case] sam: bool,
        #[case] record_count: usize,
        #[case] stops: bool,
    ) {
        use fgumi_bgzf::BGZF_EOF;
        use fgumi_raw_bam::flags::UNMAPPED;
        use std::process::{Command, Stdio};
        use std::sync::mpsc::{TryRecvError, sync_channel};
        use std::time::{Duration, Instant};

        let chunks = same_queryname_output_chunks(sam, record_count);
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .expect("spawn streaming output");
        let mut stdin = cat.stdin.take().unwrap();
        let stdout = cat.stdout.take().unwrap();
        let fail_after = Duration::from_millis(500);
        let (active_tx, active_rx) = channel();
        let producer = thread::spawn(move || {
            for chunk in chunks {
                stdin.write_all(&chunk).unwrap();
                stdin.flush().unwrap();
                thread::sleep(Duration::from_millis(50));
            }
            active_tx.send(()).unwrap();
            if stops {
                thread::sleep(fail_after * 3);
            }
            if !sam {
                stdin.write_all(&BGZF_EOF).unwrap();
            }
        });

        let unmapped = BamTemplateBatch::new(
            0,
            vec![Template::from_records(vec![make_record(b"large", UNMAPPED)]).unwrap()],
        );
        let ceiling = unmapped.heap_size() as u64 + 100;
        let gate =
            Arc::new(InFlightGate::new(ceiling, ceiling, Duration::from_millis(50), fail_after));
        assert!(gate.acquire(ceiling).unwrap());
        let (token_tx, token_rx) = channel();
        token_tx.send(BatchToken { unmapped, n_templates: 1, serial: 0 }).unwrap();
        drop(token_tx);
        let (out_tx, out_rx) = sync_channel(1);
        let cfg = make_test_cfg();
        let header_handle = cfg.header_handle.clone();
        let reader_gate = Arc::clone(&gate);
        let reader = thread::spawn(move || {
            reader_loop_inner(token_rx, &out_tx, stdout, &cfg, &SharedState::new(), &reader_gate)
        });
        let deadline = Instant::now() + Duration::from_secs(5);
        while !header_handle.is_set() || gate.in_flight() != 100 {
            assert!(Instant::now() < deadline, "reader did not start its token");
            thread::sleep(Duration::from_millis(5));
        }
        let waiter_gate = Arc::clone(&gate);
        let (wait_tx, wait_rx) = channel();
        let waiter = thread::spawn(move || {
            wait_tx.send(waiter_gate.acquire(ceiling)).unwrap();
        });

        active_rx.recv_timeout(Duration::from_secs(10)).expect("output keeps flowing");
        let while_active = wait_rx.try_recv();
        producer.join().unwrap();
        reader.join().unwrap().expect("reader completes the template");
        cat.wait().unwrap();
        let while_reserved = wait_rx.try_recv();
        let reserved = gate.in_flight();
        gate.release(100);
        waiter.join().unwrap();

        assert!(
            matches!(while_active, Err(TryRecvError::Empty)),
            "output must stay active longer than the timeout without failing: {while_active:?}"
        );
        assert_eq!(reserved, 100, "record progress must not release unmapped bytes");
        if stops {
            let err = while_reserved.expect("stopped output must time out").unwrap_err();
            assert!(err.to_string().contains("read ahead more than"), "unexpected: {err}");
        } else {
            assert!(
                matches!(while_reserved, Err(TryRecvError::Empty)),
                "active output must neither time out nor admit bytes: {while_reserved:?}"
            );
            assert!(wait_rx.recv().unwrap().unwrap(), "release admits the waiting batch");
        }
        let batch = out_rx.recv().unwrap();
        assert_eq!(batch.mapped[0].records().len(), record_count);
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
        reader_loop_inner(
            token_rx,
            &out_tx,
            stdout,
            &cfg,
            &SharedState::new(),
            &InFlightGate::new(
                u64::MAX,
                u64::MAX,
                IN_FLIGHT_STALL_GRACE,
                IN_FLIGHT_CEILING_TIMEOUT,
            ),
        )
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
