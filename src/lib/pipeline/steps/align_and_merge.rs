//! `AlignAndMergeStep` — typed `Step` that wraps an aligner
//! subprocess and pairs aligner output with original unmapped tags.
//!
//! Replaces the legacy hand-rolled orchestrator in
//! `src/lib/align_and_merge.rs`: today's orchestrator drains aligner
//! stdout to a tempfile and then re-opens the input BAM via
//! `Zipper::execute`. The new Step streams template-by-template with
//! no tempfile bridge.
//!
//! ## Shape
//!
//! ```text
//!  upstream (BamTemplateBatch)
//!     │
//!     ▼
//!  ┌──────────────── AlignAndMergeStep ─────────────────┐
//!  │  Owns:                                              │
//!  │    • AlignerProcess (subprocess + stderr ring)      │
//!  │    • stdin-writer thread (FASTQ → aligner stdin)    │
//!  │    • stdout-reader thread (SAM/BAM parse only —     │
//!  │      no merge; emits ZipperBatch for try_run)       │
//!  │    • Three bounded channels (in_chan, token_chan,   │
//!  │      out_chan) and a shared error slot              │
//!  │                                                     │
//!  │  Step::try_run (Serial; any worker can dispatch):   │
//!  │    upstream → in_chan → (writer) → aligner stdin    │
//!  │    aligner stdout → (reader) → out_chan ──┐         │
//!  │                                            ▼         │
//!  │                              merge_zipper_batch      │
//!  │                              (merge_raw + bisulfite) │
//!  │                                            │         │
//!  │                                            ▼         │
//!  │                                       downstream     │
//!  └──────────────────────────────────────────────────────┘
//! ```
//!
//! Why `Serial` instead of `Exclusive`: `Exclusive` pinned dispatch to
//! one worker, which meant any other step's blocking work on that
//! worker (e.g. a held-slot flush spin-retrying on a full downstream)
//! could starve AAM and deadlock the pipeline.
//! `Serial` lets any free worker dispatch AAM; the framework mutex
//! preserves the "one dispatcher at a time" invariant that
//! `in_tx`/`out_rx`/`held_*` require. See
//! `docs/design/aam-bridge-refactor.md` for the full diagnosis.
//!
//! Why merge happens in `try_run`, not on the reader: the reader is
//! the sole consumer of bwa's stdout, and `merge_raw` per template
//! is CPU-heavy (especially with bisulfite restore). Doing the merge
//! on the reader makes bwa stall on stdout flush when merge is the
//! bottleneck. Emitting raw `ZipperBatch` (paired mapped + unmapped)
//! and merging in `try_run` keeps the reader I/O-bound and lets
//! whichever worker dispatches AAM amortize merge across the same
//! `OUT_CHAN_DEPTH` of in-flight batches.
//!
//! ## [`BatchToken`] / [`ZipperBatch`] protocol (index-based pairing)
//!
//! Pairing of unmapped templates with their alignments is **structural
//! by index, not queryname**. The writer thread keeps a per-batch
//! template count and ships a `BatchToken { unmapped, n_templates,
//! serial }` through `token_chan` to the reader thread. The reader
//! pops one token at a time, reads exactly `n_templates` templates
//! from aligner stdout, packages them into a `ZipperBatch { serial,
//! mapped, unmapped }`, and sends that on `out_chan`. `try_run`
//! pops `ZipperBatch`es and calls `merge_zipper_batch` (which runs
//! `merge_raw` + optional bisulfite restore per template-pair) to
//! produce the merged `BamTemplateBatch` pushed downstream.
//!
//! This relies on the aligner preserving input record order (bwa-mem
//! and bwa-mem3 do this when `-K` chunk size is set, which our presets
//! require). A debug-only queryname-equality assertion catches order
//! violations during development; production runs trust the contract.
//!
//! ## Header propagation
//!
//! The reader thread's first action (before any record) is to parse
//! the aligner's emitted SAM header, merge it with the
//! construction-time partial header (dict-derived `@SQ` + unmapped
//! `@HD`/`@CO`/`@RG`/`@PG` + fgumi's own `@PG`), and resolve the
//! supplied [`HeaderHandle`]. The downstream `WriteBgzfFile`
//! (constructed via `new_with_handle`) polls the handle and proceeds
//! once it's resolved.
//!
//! ## Lifecycle
//!
//! - `new` spawns the subprocess + both daemon threads + sets up
//!   channels.
//! - `try_run` is dispatcher-driven, non-blocking, with two `HeldSlot`s
//!   (`held_in` for upstream→`in_chan` backpressure, `held_out` for
//!   merged-batch→downstream backpressure). It runs three phases across
//!   re-dispatches: ingest (pump `ctx.input` → `in_tx` while merging the
//!   aligner's output), close-stdin (once `ctx.input` is drained and
//!   `held_in` is flushed, drop `in_tx` → aligner stdin EOF), and drain
//!   (pump the aligner's tail through `out_rx` until the reader has exited
//!   and the receiver is empty, then `finalize`). Returning
//!   `StepOutcome::Finished` only at the end of the drain phase means the
//!   final flush yields and is re-dispatched rather than blocking the worker.
//! - `finalize` joins both daemons (which cannot block — the reader has
//!   exited and the receiver is drained), waits the subprocess, and surfaces
//!   the stderr ring on non-zero exit in order (aligner → reader → writer).
//! - `Drop` is the resource-cleanup backstop on error paths — the clean
//!   `finalize` is never reached when any step returns `Err` (the
//!   framework signals teardown via `PipelineSignal::is_done`).
//!
//! ## Output formats
//!
//! The reader auto-detects BGZF magic on the aligner's first 4
//! bytes:
//! * BGZF → noodles BAM reader, zero-copy per-record via
//!   [`fgumi_raw_bam::RawBamReader`].
//! * Anything else → noodles SAM text reader, one
//!   `RecordBuf → RawRecord` encode per record via
//!   [`fgumi_raw_bam::encode_record_buf_to_raw`].
//!
//! BAM is the higher-throughput path because there's no encode
//! cost; aligners that can emit BAM (e.g. via a `samtools view -bu -`
//! pipe in the user's command) avoid the SAM encode entirely.

use std::io::{self, BufReader, BufWriter, Read, Write};
use std::process::ChildStdout;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, Sender, SyncSender, TrySendError, channel, sync_channel};
use std::thread::{self, JoinHandle};

use noodles::sam::Header;
use parking_lot::Mutex;

use crate::aligner::AlignerProcess;
use crate::commands::fastq::write_fastq_record;
use crate::commands::zipper::merge_one_template;
use crate::pipeline::core::header::HeaderHandle;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::{
    BranchOrdering, HeldSlot, OrderedBytesSingle, QueueSpec, Step, StepCtx, StepKind, StepOutcome,
    StepProfile, Unpushed,
};
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::reference::ReferenceReader;
use crate::template::Template;
use crate::umi::TagInfo;

/// Number of aligner stderr lines retained for failure diagnostics.
/// Same value the legacy orchestrator uses (`align_and_merge.rs:74`).
const ALIGNER_STDERR_RING_SIZE: usize = 50;

/// First four bytes of every BGZF/BAM stream. Used to auto-detect
/// whether the aligner emits BAM (BGZF-compressed) or SAM (text). See
/// `align_and_merge.rs:79`.
const BGZF_MAGIC: [u8; 4] = [0x1f, 0x8b, 0x08, 0x04];

/// Flag bits filtered from records written as FASTQ to the aligner.
/// `0x900` = `SECONDARY | SUPPLEMENTARY`. Same rationale as the legacy
/// orchestrator (`align_and_merge.rs:91`): if a previously-aligned BAM
/// is passed by mistake, dropping these prevents duplicate FASTQ
/// emissions for the same template.
const FASTQ_WRITER_EXCLUDE_FLAGS: u16 =
    fgumi_raw_bam::flags::SECONDARY | fgumi_raw_bam::flags::SUPPLEMENTARY;

/// Bound for `in_chan` (dispatcher → writer thread). Two batches in
/// flight match bwa's `kt_pipeline` `p_nt=2` double-buffer: one batch
/// being aligned, one queued. See the design doc
/// `align-and-merge-as-step.md` for the derivation from bwa-mem2's
/// `fastmap.cpp`.
///
/// Note: the writer→reader **token** channel is intentionally
/// unbounded (`std::sync::mpsc::channel()`), NOT bounded by this
/// constant. See the comment at the `token_chan` construction site
/// in [`AlignAndMergeStep::new`] for the deadlock rationale: bwa's
/// `-K` flag lets it buffer many input batches before emitting any
/// output, so the writer must be able to keep pushing tokens past
/// the reader without blocking — otherwise the writer is stuck on
/// `token_tx.send` and never closes bwa's stdin, and bwa never
/// flushes. Backpressure on the writer-side flow is enforced
/// naturally by the OS pipe buffer (~64 KiB on macOS) plus bwa's
/// `-K`-bounded internal buffer, both of which throttle the
/// FASTQ-write path well before the token channel ever fills up.
const IN_CHAN_DEPTH: usize = 2;

/// Bound for `out_chan` (reader thread → dispatcher). Two batches lets
/// the reader stay one ahead while the dispatcher pumps the previous
/// out to downstream; backpressure cascades back when downstream
/// stalls.
const OUT_CHAN_DEPTH: usize = 2;

/// Estimated in-flight bytes per aligner-input base, used to derive the
/// in-flight-unmapped byte budget from the aligner's `-K` chunk size.
/// One base costs ~1 byte of sequence + ~1 byte of quality held in the
/// unmapped `RawRecord`, plus read-name / tag overhead — 3 is a
/// deliberately generous estimate so the budget never under-shoots the
/// aligner's real buffering (under-shooting risks a writer stall).
const IN_FLIGHT_BYTES_PER_BASE: u64 = 3;

/// Multiplier over a single `-K` chunk for the in-flight budget. bwa's
/// `kt_pipeline` double-buffers (`p_nt=2`: one chunk aligning, one
/// queued), so the aligner can hold ~2 chunks before emitting; 4× adds
/// slack for the FASTQ-write buffer and pipeline jitter so the writer
/// never blocks before the aligner has a full chunk to emit (which
/// would deadlock — see `InFlightGate`).
const IN_FLIGHT_CHUNK_MULTIPLIER: u64 = 4;

/// Floor for the in-flight-unmapped budget. Keeps the budget workable
/// for small `-K` values (custom aligner commands, tests) where the
/// `-K`-derived value would be tiny, and covers aligners whose internal
/// buffering exceeds their nominal `-K`.
const IN_FLIGHT_MIN_BUDGET: u64 = 512 * 1024 * 1024;

/// Derive the in-flight-unmapped byte budget for AAM from the aligner's
/// `-K` chunk size (bases per batch). The budget bounds the otherwise
/// unbounded writer→reader token backlog: the writer blocks once the
/// in-flight unmapped bytes reach this, so a fast-draining aligner can't
/// accumulate the whole input's unmapped reads in RAM (issue #382).
///
/// Sized `≥` the aligner's real buffering so a **streaming** aligner
/// (one that emits output after ≤ `-K` bases — every real aligner)
/// always has a full chunk to emit before the writer blocks, so the
/// block is transient, not a deadlock. A non-streaming aligner (one that
/// reads all stdin before emitting) is out of contract: it would stall
/// the writer here rather than OOM. For the default `-K` of 150M bases
/// this is ~1.8 GiB.
#[must_use]
pub fn in_flight_budget_for_chunk_size(chunk_size_bases: u64) -> u64 {
    chunk_size_bases
        .saturating_mul(IN_FLIGHT_BYTES_PER_BASE)
        .saturating_mul(IN_FLIGHT_CHUNK_MULTIPLIER)
        .max(IN_FLIGHT_MIN_BUDGET)
}

/// Byte-budget gate bounding AAM's in-flight unmapped reads (the
/// writer→reader token backlog). The writer [`acquire`](InFlightGate::acquire)s
/// before feeding a batch to the aligner and the reader
/// [`release`](InFlightGate::release)s when it consumes the matching
/// token, so resident in-flight unmapped bytes stay near `budget`.
///
/// Deadlock-safety: the gate blocks the writer only while in-flight is at
/// budget AND non-empty — a single oversized batch always passes when the
/// gate is empty. Because the budget is sized `≥` a streaming aligner's
/// `-K` buffering (see [`in_flight_budget_for_chunk_size`]), the aligner
/// always has a full chunk to emit before the writer blocks, so the
/// reader drains and the block lifts. If the reader exits (EOF or error),
/// it latches `consumer_gone` so a blocked writer wakes and bails rather
/// than hanging.
struct InFlightGate {
    budget: u64,
    inner: std::sync::Mutex<GateInner>,
    cond: std::sync::Condvar,
}

struct GateInner {
    in_flight: u64,
    consumer_gone: bool,
}

impl InFlightGate {
    fn new(budget: u64) -> Self {
        Self {
            budget,
            inner: std::sync::Mutex::new(GateInner { in_flight: 0, consumer_gone: false }),
            cond: std::sync::Condvar::new(),
        }
    }

    /// Reserve `n` in-flight bytes, blocking while the gate is full and
    /// non-empty. Returns `false` if the consumer (reader) has gone — the
    /// caller (writer) should then stop.
    fn acquire(&self, n: u64) -> bool {
        let mut g = self.inner.lock().expect("InFlightGate mutex poisoned");
        // Block while non-empty AND this reservation would exceed budget.
        // The `in_flight != 0` guard lets a single oversized batch through
        // when the gate is empty (it can't be split, so holding it is
        // unavoidable) — without it the writer would deadlock on a batch
        // larger than the whole budget.
        while !g.consumer_gone && g.in_flight != 0 && g.in_flight.saturating_add(n) > self.budget {
            g = self.cond.wait(g).expect("InFlightGate mutex poisoned");
        }
        if g.consumer_gone {
            return false;
        }
        g.in_flight = g.in_flight.saturating_add(n);
        true
    }

    /// Release `n` in-flight bytes (reader consumed a token) and wake the
    /// writer if it is blocked.
    fn release(&self, n: u64) {
        let mut g = self.inner.lock().expect("InFlightGate mutex poisoned");
        g.in_flight = g.in_flight.saturating_sub(n);
        drop(g);
        self.cond.notify_all();
    }

    /// Latch that the consumer (reader) has exited, so a blocked writer
    /// wakes and bails instead of hanging forever.
    fn mark_consumer_gone(&self) {
        let mut g = self.inner.lock().expect("InFlightGate mutex poisoned");
        g.consumer_gone = true;
        drop(g);
        self.cond.notify_all();
    }
}

/// RAII guard that latches consumer-gone on drop — on normal return, on an
/// error return, AND on a panic unwinding out of `reader_loop`. Without the
/// guard the panic path would skip `mark_consumer_gone`, leaving a writer
/// parked in `gate.acquire` blocked forever and `Drop`'s `writer_thread.join()`
/// hung with it.
struct ConsumerGoneGuard<'a>(&'a InFlightGate);

impl Drop for ConsumerGoneGuard<'_> {
    fn drop(&mut self) {
        self.0.mark_consumer_gone();
    }
}

/// FASTQ-writer [`BufWriter`] capacity. 1 MiB matches the legacy
/// orchestrator's choice (`align_and_merge.rs:259`): bigger than the
/// OS pipe buffer (~64 KiB) so we amortize `write` syscalls, but not
/// so large that we delay the aligner's first chunk perceptibly.
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

/// Immutable, Arc-friendly configuration shared by [`AlignAndMergeStep`].
///
/// One `Arc` clone of the inner is made per spawned thread; the cost is
/// negligible compared to the data flowing through.
#[derive(Clone)]
pub struct AlignAndMergeConfig {
    /// Tag-merge rules (remove / reverse / revcomp) plumbed to
    /// `merge_raw`.
    pub tag_info: Arc<TagInfo>,

    /// Whether to skip TC (template-coordinate) tag handling. Mirrors
    /// `ZipperMergeConfig::skip_tc_tags`.
    pub skip_tc_tags: bool,

    /// Optional reference reader used by
    /// `restore_unconverted_bases_in_raw_template` (bisulfite path).
    /// `None` for normal alignment.
    pub reference: Option<Arc<ReferenceReader>>,

    /// Partial output header built at construction time
    /// (dict-derived `@SQ` + unmapped-derived
    /// `@HD`/`@CO`/`@RG`/`@PG` + fgumi's own `@PG`). The reader
    /// thread merges the aligner's emitted `@PG`/`@CO`/`@RG` lines
    /// into this and resolves `header_handle`.
    pub partial_output_header: Arc<Header>,

    /// One-shot handle the reader thread resolves with the merged
    /// header. The downstream `WriteBgzfFile` (constructed via
    /// `new_with_handle`) polls this.
    pub header_handle: HeaderHandle,

    /// Counter for records emitted to downstream. Exposed back to the
    /// caller after `Pipeline::run` returns so summary logging can
    /// report a real throughput number.
    pub records_emitted: Arc<AtomicU64>,

    /// Byte limit for the **downstream** output queue
    /// (`OrderedBytesSingle` `ByteBounded`) that AAM pushes merged
    /// `BamTemplateBatch`es onto.
    ///
    /// AAM also has an *internal* `out_chan` (`SyncSender<ZipperBatch>`,
    /// depth `OUT_CHAN_DEPTH = 2`) carrying the reader-to-`try_run`
    /// hand-off. Because `ZipperBatch` carries both `mapped` and
    /// `unmapped` halves (paired template Vecs), its per-item heap
    /// is roughly **2×** the per-item heap of the old design's
    /// `BamTemplateBatch` (which only carried merged templates).
    /// Combined with `held_out` the worst-case in-flight footprint
    /// is ≈ 3 × 2 × `batch_bytes`.
    ///
    /// For default 1000-template batches this is ~600 KB, negligible.
    /// For production CODEC pipelines with multi-MB batches it can
    /// approach ~1.5 GB additional in-flight; if memory is tight,
    /// drop `OUT_CHAN_DEPTH` to 1 (one rebuild) — the reader
    /// throughput cost is paid only when bwa stalls on stdout for
    /// one batch worth of time, which is rare on production hardware.
    /// This field gates the downstream queue, not the internal one.
    pub output_byte_limit: u64,

    /// Byte budget for AAM's **internal** in-flight unmapped reads (the
    /// writer→reader token backlog). The writer blocks once in-flight
    /// unmapped bytes reach this, bounding the otherwise-unbounded backlog
    /// so a fast-draining aligner can't accumulate the whole input's
    /// unmapped reads in RAM (issue #382). Derive it from the aligner's
    /// `-K` chunk size via [`in_flight_budget_for_chunk_size`].
    pub in_flight_unmapped_budget: u64,
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
/// `ZipperBatch { mapped, unmapped, serial }` for `try_run` to merge.
struct BatchToken {
    unmapped: BamTemplateBatch,
    n_templates: usize,
    /// Serial number assigned by the writer in push order. Carried
    /// through into the `ZipperBatch` and ultimately the emitted
    /// `BamTemplateBatch::batch_serial` so downstream's
    /// `ByItemOrdinal` consumer sees a monotonic sequence.
    serial: u64,
}

/// Reader-thread → `try_run` intermediate. Aligner-parsed but
/// pre-merge templates paired with their original unmapped halves.
/// Merging (`merge_raw` + optional bisulfite restore) happens later
/// on whatever worker dispatches AAM; the reader stays I/O-bound.
#[derive(Debug)]
struct ZipperBatch {
    serial: u64,
    mapped: Vec<Template>,
    unmapped: BamTemplateBatch,
}

impl HeapSize for ZipperBatch {
    fn heap_size(&self) -> usize {
        // Matches `BamTemplateBatch::total_bytes` semantics: sum of
        // `Template::heap_size` only; no `Vec`-allocation overhead.
        // `Template` carries only an inherent `heap_size` method (not
        // a `HeapSize` trait impl) to keep `template.rs` framework-agnostic,
        // so the sum is hand-rolled here.
        let mapped_heap: usize = self.mapped.iter().map(Template::heap_size).sum();
        mapped_heap + self.unmapped.heap_size()
    }
}

impl Ordered for ZipperBatch {
    fn ordinal(&self) -> u64 {
        self.serial
    }
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

/// `Exclusive` Step that wraps an aligner subprocess + two internal
/// threads + the inline `merge_raw` of tags from unmapped templates
/// onto the aligned mapped templates.
///
/// See module doc for design rationale.
pub struct AlignAndMergeStep {
    cfg: AlignAndMergeConfig,

    /// Owned aligner subprocess. Taken out in `finalize` (or `Drop` on the
    /// error path) for `wait()`.
    aligner: Option<AlignerProcess>,

    /// Sender side of `in_chan`. Dropped in `try_run`'s close-stdin phase to
    /// signal EOF to the writer thread (which closes aligner stdin
    /// in turn).
    in_tx: Option<SyncSender<BamTemplateBatch>>,

    /// Receiver side of `out_chan`. The reader thread sends raw
    /// (pre-merge) `ZipperBatch`es here; `try_run` does the per-
    /// template `merge_raw` + bisulfite restore, then pushes a
    /// merged `BamTemplateBatch` to `ctx.outputs`.
    out_rx: Option<Receiver<ZipperBatch>>,

    writer_thread: Option<JoinHandle<io::Result<()>>>,
    reader_thread: Option<JoinHandle<io::Result<()>>>,

    shared: Arc<SharedState>,

    /// Held slot for a batch that was popped from `ctx.input` but
    /// couldn't be `try_send`'d to `in_chan` (writer thread is slow).
    held_in: HeldSlot<BamTemplateBatch>,

    /// Held slot for an `Unpushed<BamTemplateBatch>` we received from
    /// the reader but couldn't `push` to `ctx.outputs` (downstream
    /// backpressure).
    held_out: HeldSlot<Unpushed<BamTemplateBatch>>,

    name: &'static str,
}

impl AlignAndMergeStep {
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
    pub fn new(cfg: AlignAndMergeConfig, aligner_command: &str) -> io::Result<Self> {
        let mut aligner = AlignerProcess::spawn(aligner_command, ALIGNER_STDERR_RING_SIZE)
            .map_err(|e| io::Error::other(format!("AlignAndMergeStep::new: spawn: {e:#}")))?;

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
                        "AlignAndMergeStep::new: failed to spawn writer thread: {e}"
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
                        "AlignAndMergeStep::new: failed to spawn reader thread: {e}"
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
            name: "AlignAndMerge",
        })
    }
}

// ──────────────────────────────────────────────────────────────────────────
// Writer thread
// ──────────────────────────────────────────────────────────────────────────

/// Drain `in_rx` to the aligner's stdin as interleaved FASTQ, pushing
/// one [`BatchToken`] per batch to `token_tx` so the reader can pair
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
    let mut seq_buf: Vec<u8> = Vec::with_capacity(512);
    let mut qual_buf: Vec<u8> = Vec::with_capacity(512);
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
        // error. AAM's expected input is an unmapped BAM (from
        // `fgumi extract`), which has no secondaries — a
        // fully-filtered template here means the user passed a
        // re-aligned BAM by mistake. Error out loudly with the
        // queryname so the misconfiguration is obvious.
        let mut n_templates = 0usize;
        for template in batch.templates() {
            let mut wrote_any = false;
            for record in &template.records {
                let flags = record.flags();
                if (flags & FASTQ_WRITER_EXCLUDE_FLAGS) != 0 {
                    continue;
                }
                if let Err(e) = write_fastq_record(
                    &mut writer,
                    record,
                    flags,
                    /* no_suffix */ true,
                    /* umi_header */ None,
                    &mut seq_buf,
                    &mut qual_buf,
                ) {
                    // Record + return the same full error message so
                    // `finalize`'s join path and the shared
                    // error_slot fallback path both see the same
                    // text — no information loss either way. Same
                    // pattern as the no-primary path below.
                    let msg = format!(
                        "AAM writer: write_fastq_record for template '{name}': {e:#}",
                        name = String::from_utf8_lossy(&template.name),
                    );
                    shared.record_error(io::Error::other(msg.clone()));
                    return Err(io::Error::other(msg));
                }
                wrote_any = true;
            }
            if !wrote_any {
                let msg = format!(
                    "AAM writer: template '{name}' has no primary records (all flagged \
                     SECONDARY/SUPPLEMENTARY). AAM expects unmapped BAM input — did you pass \
                     a re-aligned BAM by mistake? Run `fgumi extract` first or pre-filter \
                     the input.",
                    name = String::from_utf8_lossy(&template.name),
                );
                shared.record_error(io::Error::other(msg.clone()));
                return Err(io::Error::other(msg));
            }
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
        let io_err = io::Error::other(format!("AAM writer: flush: {e}"));
        shared.record_error(io_err);
        return Err(io::Error::other("AAM writer: flush failed"));
    }
    drop(writer);
    Ok(())
}

// ──────────────────────────────────────────────────────────────────────────
// Reader thread
// ──────────────────────────────────────────────────────────────────────────

/// Parse the aligner's stdout, pair each emitted template with its
/// corresponding unmapped via the token channel, and emit
/// `ZipperBatch` to `out_tx`. Merging (`merge_raw` and any bisulfite
/// restore) happens in `AlignAndMergeStep::try_run`, run by whatever
/// worker the framework dispatches; the reader stays I/O-bound so
/// bwa's stdout never stalls behind CPU-heavy merge work.
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
    cfg: AlignAndMergeConfig,
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
            let _ = cfg.header_handle.poison(io::Error::new(e.kind(), format!("AAM reader: {e}")));
            shared.record_error(io::Error::new(e.kind(), format!("AAM reader: {e}")));
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
    cfg: &AlignAndMergeConfig,
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
            "AlignAndMergeStep: aligner {marker}. Likely causes: subprocess \
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
        let stream =
            TemplateStream::Bam(BamTemplateStream::new(fgumi_raw_bam::RawBamReader::new(bgzf)));
        (stream, aligner_header)
    } else {
        let buffered = BufReader::with_capacity(STDOUT_READER_BUF_BYTES, stdout);
        let mut sam_reader = noodles::sam::io::Reader::new(buffered);
        let aligner_header = sam_reader
            .read_header()
            .map_err(|e| io::Error::other(format!("aligner SAM header: {e}")))?;
        let header_arc = Arc::new(aligner_header.clone());
        let stream = TemplateStream::Sam(SamTemplateStream::new(sam_reader, header_arc));
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
    // AAM is the sole producer. A pre-set handle indicates either a
    // test (where `from_header` was used) or a wiring bug. We log
    // the latter case at warn level so it surfaces in CI / prod
    // pipelines while letting tests proceed unchanged. The merged
    // header is dropped in either case; the existing value wins.
    if cfg.header_handle.set(merged).is_err() {
        log::warn!(
            "AlignAndMergeStep: HeaderHandle was already resolved before the aligner \
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
        // `held_out` + merge, which is separately bounded — so releasing
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

        for i in 0..token.n_templates {
            let mapped = template_stream.next_template()?;

            let mapped = mapped.ok_or_else(|| {
                io::Error::other(format!(
                    "AAM reader: aligner stdout EOF after {emitted} of {expected} \
                     templates in batch {serial} — aligner emitted fewer alignments \
                     than input reads",
                    emitted = i,
                    expected = token.n_templates,
                    serial = token.serial,
                ))
            })?;

            // Cheap dev-time integrity check; production trusts
            // bwa's order preservation. Done here (in the reader,
            // pre-merge) because the queryname comparison only needs
            // both templates' `name` fields, which both halves
            // already carry.
            debug_assert_eq!(
                token.unmapped.templates()[i].name,
                mapped.name,
                "AAM reader: queryname mismatch — unmapped[{i}]='{u}' but mapped='{m}' \
                 (aligner reordered output?)",
                i = i,
                u = String::from_utf8_lossy(&token.unmapped.templates()[i].name),
                m = String::from_utf8_lossy(&mapped.name),
            );

            mapped_templates.push(mapped);
        }

        let zb = ZipperBatch {
            serial: token.serial,
            mapped: mapped_templates,
            unmapped: token.unmapped,
        };
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
            "AAM reader: aligner emitted at least one more template than the sum of \
             input batch sizes (residual record carried after final token)",
        ));
    }

    // Probe for one more record past the final token. If the
    // aligner emitted multiple extra templates this only detects
    // the first; the count isn't useful for the error message, so
    // we report "at least one" rather than a plural that would be
    // wrong for exactly-one cases.
    if template_stream.probe_trailing().map_err(|e| {
        io::Error::other(format!("AAM reader: probing trailing aligner output: {e}"))
    })? {
        return Err(io::Error::other(
            "AAM reader: aligner emitted at least one more record than the sum of \
             input batch sizes (extra record at EOF)",
        ));
    }

    Ok(())
}

/// Merge one `ZipperBatch` into a `BamTemplateBatch`. Per-template
/// `merge_raw` plus optional bisulfite restore; folds record-count
/// and heap-size accounting into the same single pass so the
/// resulting `BamTemplateBatch` doesn't re-walk the templates to
/// compute `total_bytes`.
///
/// **Single-threaded by construction.** This runs synchronously on
/// whichever worker dispatches AAM (Serial); merge throughput is
/// effectively single-core. If profiling reveals merge as a
/// bottleneck (most likely on bisulfite workloads where
/// `restore_unconverted_bases_in_raw_template` does meaningful
/// per-record CPU), the right move is to extract a `Parallel`
/// `MergeAligned` step that consumes `ZipperBatch` from a multi-
/// consumer queue and emits `BamTemplateBatch`. The `ZipperBatch`
/// intermediate is shaped for that promotion — no fixture/test
/// refactor required at the call site. See
/// `docs/design/aam-bridge-refactor.md` §8 (out-of-scope follow-ups).
fn merge_zipper_batch(zb: ZipperBatch, cfg: &AlignAndMergeConfig) -> io::Result<BamTemplateBatch> {
    let ZipperBatch { serial, mapped, unmapped } = zb;
    debug_assert_eq!(
        mapped.len(),
        unmapped.templates().len(),
        "ZipperBatch invariant: mapped and unmapped must have equal length",
    );

    let mut merged: Vec<Template> = Vec::with_capacity(mapped.len());
    let mut total_records: u64 = 0;
    let mut total_bytes: usize = 0;
    for (mut mapped_template, unmapped_template) in mapped.into_iter().zip(unmapped.templates()) {
        merge_one_template(
            unmapped_template,
            &mut mapped_template,
            &cfg.tag_info,
            cfg.skip_tc_tags,
            cfg.reference.as_deref(),
            &cfg.partial_output_header,
        )
        .map_err(|e| io::Error::other(format!("AAM: {e:#}")))?;

        total_records += mapped_template.records.len() as u64;
        total_bytes += mapped_template.heap_size();
        merged.push(mapped_template);
    }

    cfg.records_emitted.fetch_add(total_records, Ordering::Relaxed);
    Ok(BamTemplateBatch::from_parts(serial, merged, total_bytes))
}

/// Peek the first 4 bytes of `stdout` to detect BGZF magic. Returns
/// `(is_bgzf, peek_filled, chained_reader)` where:
/// - `is_bgzf` is true iff exactly 4 bytes were peeked and they match
///   `BGZF_MAGIC`.
/// - `peek_filled` is the number of bytes actually peeked (0..=4).
///   `0` means the aligner exited before producing any output, which
///   the caller can distinguish from a "wrong format" condition for
///   a clearer error message.
/// - `chained_reader` re-emits the peeked bytes followed by the rest
///   of the original stream — callers must drain the returned reader,
///   not the original.
///
/// The return type is `Box<dyn Read + Send>` so the same boxed reader
/// can feed both the BGZF path (BAM) and any future SAM-text path
/// without changing every downstream type signature.
///
/// Lifted from `align_and_merge.rs:165-187`.
fn peek_aligner_format(mut stdout: ChildStdout) -> io::Result<(bool, usize, Box<dyn Read + Send>)> {
    let mut peek_buf = [0u8; 4];
    let mut peek_filled = 0usize;
    while peek_filled < peek_buf.len() {
        match stdout.read(&mut peek_buf[peek_filled..]) {
            Ok(0) => break, // aligner exited with < 4 bytes of output
            Ok(n) => peek_filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    let is_bgzf = peek_filled == 4 && peek_buf == BGZF_MAGIC;
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
///   The SAM path costs one `RecordBuf` → `RawRecord` encode per
///   record (vs zero-copy on BAM); aligners emitting BAM are
///   preferred for high-throughput workloads.
///
/// The prior caller-threaded `(scratch, peeked)` API is gone —
/// centralising the state inside this enum makes the
/// "forgot-to-thread-peeked-and-lost-a-record" footgun
/// syntactically impossible.
enum TemplateStream {
    Bam(BamTemplateStream),
    Sam(SamTemplateStream),
}

impl TemplateStream {
    fn next_template(&mut self) -> io::Result<Option<Template>> {
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

struct BamTemplateStream {
    reader: AlignerBamReader,
    scratch: fgumi_raw_bam::RawRecord,
    peeked: Option<fgumi_raw_bam::RawRecord>,
    /// Reusable per-template name buffer. Re-cleared on each
    /// `next_template` call; avoids ~one `Vec<u8>` allocation
    /// per template (~40 bytes typical) on hot loops.
    name_buf: Vec<u8>,
}

impl BamTemplateStream {
    fn new(reader: AlignerBamReader) -> Self {
        Self {
            reader,
            scratch: fgumi_raw_bam::RawRecord::new(),
            peeked: None,
            name_buf: Vec::with_capacity(64),
        }
    }

    fn next_template(&mut self) -> io::Result<Option<Template>> {
        let mut records: Vec<fgumi_raw_bam::RawRecord> = Vec::with_capacity(2);

        if let Some(first) = self.peeked.take() {
            records.push(first);
        } else {
            let n = self.reader.read_record(&mut self.scratch)?;
            if n == 0 {
                return Ok(None);
            }
            // `scratch.clone()` here (and below) is intentional: the
            // same `scratch` buffer is reused for the *next*
            // `read_record` call (this call's loop and the following
            // `next_template` call's first read), so we can't
            // `mem::take` it without surrendering the capacity reuse
            // that the scratch buffer exists to provide — each taken
            // record would start a fresh allocation. The accepted
            // follow-up (deferred, no profile points here today) is a
            // small `RawRecord` free-list so `next_template` can swap a
            // recycled buffer in rather than clone; at 10k–100k
            // records/sec the clones are cheap. See S5a1-007.
            records.push(self.scratch.clone());
        }

        self.name_buf.clear();
        self.name_buf.extend_from_slice(records[0].read_name());

        loop {
            let n = self.reader.read_record(&mut self.scratch)?;
            if n == 0 {
                break;
            }
            if self.scratch.read_name() == self.name_buf.as_slice() {
                records.push(self.scratch.clone());
            } else {
                self.peeked = Some(self.scratch.clone());
                break;
            }
        }

        let queryname = String::from_utf8_lossy(&self.name_buf).to_string();
        let template = Template::from_records(records).map_err(|e| {
            io::Error::other(format!(
                "Template::from_records for aligner-emitted queryname '{queryname}': {e}"
            ))
        })?;
        Ok(Some(template))
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
    /// Scratch [`RecordBuf`] reused across `read_record_buf` calls.
    scratch: noodles::sam::alignment::RecordBuf,
    /// First record of the next template, encoded to `RawRecord`.
    /// Stashed when [`next_template`] reads past the current
    /// template's last record.
    peeked: Option<fgumi_raw_bam::RawRecord>,
    name_buf: Vec<u8>,
}

impl SamTemplateStream {
    fn new(reader: AlignerSamReader, header: Arc<Header>) -> Self {
        Self {
            reader,
            header,
            scratch: noodles::sam::alignment::RecordBuf::default(),
            peeked: None,
            name_buf: Vec::with_capacity(64),
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

    fn next_template(&mut self) -> io::Result<Option<Template>> {
        let mut records: Vec<fgumi_raw_bam::RawRecord> = Vec::with_capacity(2);

        if let Some(first) = self.peeked.take() {
            records.push(first);
        } else {
            match self.read_next_raw()? {
                Some(r) => records.push(r),
                None => return Ok(None),
            }
        }

        self.name_buf.clear();
        self.name_buf.extend_from_slice(records[0].read_name());

        while let Some(r) = self.read_next_raw()? {
            if r.read_name() == self.name_buf.as_slice() {
                records.push(r);
            } else {
                self.peeked = Some(r);
                break;
            }
        }

        let queryname = String::from_utf8_lossy(&self.name_buf).to_string();
        let template = Template::from_records(records).map_err(|e| {
            io::Error::other(format!(
                "Template::from_records for aligner-emitted queryname '{queryname}': {e}"
            ))
        })?;
        Ok(Some(template))
    }

    fn probe_trailing(&mut self) -> io::Result<bool> {
        Ok(self.read_next_raw()?.is_some())
    }
}

/// Compare the aligner's emitted `@SQ` table to the partial output
/// header's. The partial header's `@SQ` came from the reference dict
/// at construction time; the aligner's `@SQ` comes from whatever
/// FASTA the aligner was indexed against. If they don't match, the
/// aligner's per-record `tid` integers index into a different
/// `@SQ` ordering and the merged BAM would silently have corrupt
/// reference IDs.
///
/// Comparison is name + length only (M5/UR/AS/SP are dict-specific
/// fields the aligner doesn't propagate, so we don't require them).
fn validate_sq_consistency(partial: &Header, aligner: &Header) -> io::Result<()> {
    let partial_refs = partial.reference_sequences();
    let aligner_refs = aligner.reference_sequences();

    if partial_refs.len() != aligner_refs.len() {
        return Err(io::Error::other(format!(
            "AAM: aligner @SQ count ({}) does not match reference dict @SQ count ({}). \
             The aligner was indexed against a different FASTA than the supplied --ref. \
             Re-index or fix the --ref path.",
            aligner_refs.len(),
            partial_refs.len(),
        )));
    }

    for (idx, ((p_name, p_map), (a_name, a_map))) in
        partial_refs.iter().zip(aligner_refs.iter()).enumerate()
    {
        if p_name != a_name {
            return Err(io::Error::other(format!(
                "AAM: aligner @SQ name mismatch at position {idx}: dict='{dict_name}' \
                 aligner='{aln_name}'. The aligner was indexed against a different \
                 FASTA than the supplied --ref.",
                dict_name = String::from_utf8_lossy(p_name),
                aln_name = String::from_utf8_lossy(a_name),
            )));
        }
        if p_map.length() != a_map.length() {
            return Err(io::Error::other(format!(
                "AAM: aligner @SQ length mismatch for '{name}': dict={dict_len} \
                 aligner={aln_len}. The aligner was indexed against a different \
                 FASTA than the supplied --ref.",
                name = String::from_utf8_lossy(p_name),
                dict_len = p_map.length(),
                aln_len = a_map.length(),
            )));
        }
    }

    Ok(())
}

/// Merge aligner-emitted header lines into the partial output header.
///
/// The aligner contributes:
/// - `@PG` lines (with bwa version, command line, etc.) — appended.
///   On duplicate ID, the partial's existing PG wins (the aligner's
///   PG with that ID is dropped). **Limitation:** this does not
///   maintain the BAM spec's `@PG.PP` chain when the aligner adds
///   its own PG to a chain that already has one with the same ID
///   (rare in practice — bwa's PG ID is "bwa", which won't be in
///   `fgumi extract`'s unmapped BAM output). Maintaining a proper
///   PP chain requires constructing fresh IDs and threading the
///   PP pointer; deferred to a follow-up commit if real-world data
///   surfaces a collision.
/// - `@RG` lines (if `-R` was passed to the aligner) — appended;
///   partial's existing RG wins on duplicate ID.
/// - `@CO` comment lines — concatenated (partial first, then
///   aligner).
///
/// The aligner's `@SQ` lines are deliberately discarded:
/// [`validate_sq_consistency`] runs first to ensure the two `@SQ`
/// tables agree, after which the partial's (dict-derived) version
/// is authoritative.
fn merge_aligner_header(partial: &Header, aligner: &Header) -> Header {
    use bstr::BString;

    let mut builder = Header::builder();

    if let Some(hd) = partial.header() {
        builder = builder.set_header(hd.clone());
    }

    for (name, map) in partial.reference_sequences() {
        builder = builder.add_reference_sequence(name.clone(), map.clone());
    }

    let mut rg_seen: std::collections::HashSet<BString> = std::collections::HashSet::new();
    for (id, rg) in partial.read_groups() {
        builder = builder.add_read_group(id.clone(), rg.clone());
        rg_seen.insert(id.clone());
    }
    for (id, rg) in aligner.read_groups() {
        if !rg_seen.contains(id) {
            builder = builder.add_read_group(id.clone(), rg.clone());
        }
    }

    let mut pg_seen: std::collections::HashSet<BString> = std::collections::HashSet::new();
    for (id, pg) in partial.programs().as_ref() {
        builder = builder.add_program(id.clone(), pg.clone());
        pg_seen.insert(id.clone());
    }
    for (id, pg) in aligner.programs().as_ref() {
        if !pg_seen.contains(id) {
            builder = builder.add_program(id.clone(), pg.clone());
        }
    }

    for c in partial.comments() {
        builder = builder.add_comment(c.clone());
    }
    for c in aligner.comments() {
        builder = builder.add_comment(c.clone());
    }

    builder.build()
}

// ──────────────────────────────────────────────────────────────────────────
// Step impl
// ──────────────────────────────────────────────────────────────────────────

impl Step for AlignAndMergeStep {
    type Input = BamTemplateBatch;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            // Serial — single-dispatcher-at-a-time is the real
            // requirement (in_tx/out_rx/held_* are not Sync); the
            // framework's mutex enforces it. `Exclusive` would pin
            // dispatch to one worker, which lets any blocking call
            // on that worker (e.g. another step's held-slot flush
            // spin) starve AAM.
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.cfg.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
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
        //    reserved by the framework (see `core/step.rs:10-12`) for
        //    Serial-step mutex contention, not for full output
        //    queues. Using the wrong variant skews dispatcher stats
        //    and can change scheduling.
        if let Some(unpushed) = self.held_out.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => did_work = true,
                Err(again) => {
                    self.held_out.put(again);
                    return Ok(StepOutcome::NoProgress);
                }
            }
        }

        // 3. Pull `ZipperBatch`es from the reader, merge per-template,
        //    push merged `BamTemplateBatch` to `ctx.outputs`. The loop exits
        //    (without returning) on `Empty`/`Disconnected`. This is the
        //    steady-state pump; it does NOT by itself guarantee `out_rx` is
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
                let merged = merge_zipper_batch(zb, &self.cfg)?;
                match ctx.outputs.push(merged) {
                    Ok(()) => did_work = true,
                    Err(unpushed) => {
                        self.held_out.put(unpushed);
                        return Ok(outcome(did_work));
                    }
                }
            }
        }

        // 4. Ingest phase: while `in_tx` is open, feed the aligner. Held_out
        //    has priority over held_in deliberately (see ZipperMergeStep
        //    precedent): downstream backpressure must propagate upstream before
        //    we accept more input, otherwise the internal in_chan grows
        //    unbounded. The `in_tx` borrow is scoped to the feed so it ends
        //    before 4c may `take()` it.
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
                        return Err(io::Error::other(
                            "AAM: writer thread is gone (channel disconnected)",
                        ));
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
                        return Err(io::Error::other(
                            "AAM: writer thread is gone (channel disconnected)",
                        ));
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
            // would be lost when `finalize` drops the receiver. If a merged
            // batch can't be pushed, park it in `held_out` and re-dispatch.
            {
                let out_rx = self.out_rx.as_ref().expect("out_rx Some (checked above)");
                while let Ok(zb) = out_rx.try_recv() {
                    let merged = merge_zipper_batch(zb, &self.cfg)?;
                    if let Err(unpushed) = ctx.outputs.push(merged) {
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

impl AlignAndMergeStep {
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
            Err(_panic) => return Err(io::Error::other("AAM reader thread panicked")),
        }

        // Writer last: its only failure mode (modulo pipe-broken from an
        // aligner crash, already covered above) is a write_fastq_record bug.
        match writer_result {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_panic) => return Err(io::Error::other("AAM writer thread panicked")),
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

impl Drop for AlignAndMergeStep {
    /// Best-effort resource cleanup for error paths that bypass the clean
    /// `try_run` completion (`finalize`). The framework's `PipelineSignal`
    /// records the first step error and other workers stop on their next loop
    /// iteration; the step never reaches `finalize` when this happens, so the
    /// threads/aligner are still live and must be torn down here.
    ///
    /// **This is also the cancel-response path.** When
    /// `PipelineSignal::cancel` fires, workers exit at their next
    /// loop iteration, `Pipeline::run` returns, the `AlignAndMergeStep`
    /// `Arc` refcount hits zero, and `Drop` runs. `aligner.kill()`
    /// here SIGKILLs bwa; its pipes close immediately so the reader
    /// thread's stdout `read` unblocks within microseconds and the
    /// daemon cascade completes. Total cancel latency is bounded by:
    /// (one worker-loop iteration to observe `is_done`) + (longest
    /// in-flight `try_run` wall — single-digit ms for AAM's merge
    /// step) + (Drop+join, ≤ µs). Routing `is_done` through `try_run`
    /// to kill earlier would shave ~ms but requires either a
    /// framework change (`StepCtx` exposing `signal`) or threading
    /// the signal through construction; not worth it for the
    /// observable cancel response we ship today.
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
            "AlignAndMergeStep dropped before resolving the output header",
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
// Tests
// ──────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_cfg() -> AlignAndMergeConfig {
        AlignAndMergeConfig {
            tag_info: Arc::new(TagInfo::new(vec![], vec![], vec![])),
            skip_tc_tags: true,
            reference: None,
            partial_output_header: Arc::new(Header::default()),
            header_handle: HeaderHandle::new(),
            records_emitted: Arc::new(AtomicU64::new(0)),
            output_byte_limit: 1024 * 1024,
            in_flight_unmapped_budget: in_flight_budget_for_chunk_size(150_000_000),
        }
    }

    #[test]
    fn in_flight_budget_derivation() {
        // Default -K (150M bases) → 4 × 150M × 3 = 1.8 GiB.
        assert_eq!(in_flight_budget_for_chunk_size(150_000_000), 4 * 150_000_000 * 3);
        // Small / zero -K floors at IN_FLIGHT_MIN_BUDGET.
        assert_eq!(in_flight_budget_for_chunk_size(0), IN_FLIGHT_MIN_BUDGET);
        assert_eq!(in_flight_budget_for_chunk_size(1_000_000), IN_FLIGHT_MIN_BUDGET);
    }

    #[test]
    fn consumer_gone_latched_even_when_reader_scope_panics() {
        // A panic unwinding out of the reader scope must still latch
        // consumer-gone via `ConsumerGoneGuard`, so a writer parked in
        // `acquire` bails (returns false) instead of blocking forever. Without
        // the guard the panic would skip `mark_consumer_gone`, deadlocking the
        // writer and the `writer_thread.join()` in `Drop`.
        let gate = InFlightGate::new(1024);
        assert!(gate.acquire(10), "acquire succeeds before the consumer is gone");

        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _consumer_gone = ConsumerGoneGuard(&gate);
            panic!("simulated reader_loop panic");
        }));
        assert!(panicked.is_err(), "the panic must propagate out of the guarded scope");

        // The guard's `Drop` latched consumer-gone during unwind, so a
        // subsequent writer reservation bails instead of blocking.
        assert!(!gate.acquire(10), "a writer must bail once the reader scope has panicked");
    }

    #[test]
    fn in_flight_gate_oversized_single_batch_passes_when_empty() {
        // A batch larger than the whole budget must still pass when the gate
        // is empty — it can't be split, so holding it is unavoidable, and
        // blocking would deadlock.
        let gate = InFlightGate::new(100);
        assert!(gate.acquire(1000), "oversized batch must pass when gate empty");
    }

    #[test]
    fn in_flight_gate_blocks_until_release() {
        use std::time::Duration;
        let gate = Arc::new(InFlightGate::new(100));
        assert!(gate.acquire(100), "first acquire fills the budget");

        let g2 = Arc::clone(&gate);
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            // in_flight is 100 (full, non-empty) → this blocks until release.
            assert!(g2.acquire(50));
            tx.send(()).unwrap();
        });
        // Still blocked.
        assert!(
            rx.recv_timeout(Duration::from_millis(150)).is_err(),
            "second acquire must block while the gate is full"
        );
        gate.release(100);
        rx.recv_timeout(Duration::from_secs(5)).expect("acquire must unblock after release");
        h.join().unwrap();
    }

    #[test]
    fn in_flight_gate_consumer_gone_unblocks_and_bails() {
        use std::time::Duration;
        let gate = Arc::new(InFlightGate::new(100));
        assert!(gate.acquire(100), "fill the budget");

        let g2 = Arc::clone(&gate);
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            let ok = g2.acquire(50); // blocks (full) until consumer-gone
            tx.send(ok).unwrap();
        });
        assert!(rx.recv_timeout(Duration::from_millis(150)).is_err(), "must be blocked");
        gate.mark_consumer_gone();
        let ok = rx.recv_timeout(Duration::from_secs(5)).expect("must wake on consumer-gone");
        assert!(!ok, "acquire must return false once the consumer is gone");
        h.join().unwrap();
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
        let step = AlignAndMergeStep::new(cfg, "true");
        // `true` may or may not have stdin/stdout pipes open by the
        // time we get here — both outcomes are valid; we just check
        // spawn doesn't panic.
        let _ = step;
    }

    #[test]
    fn profile_advertises_serial_nonsticky_byordinal() {
        // v2.2 refactor: AAM was `Exclusive`+`sticky` (which caused
        // the worker-affinity hang); is now `Serial`+`!sticky` so any
        // free worker can dispatch when an owner is stuck elsewhere.
        // See docs/design/aam-bridge-refactor.md §0.
        let cfg = make_test_cfg();
        let step = AlignAndMergeStep::new(cfg, "cat").expect("spawn cat");
        let p = step.profile();
        assert_eq!(p.name, "AlignAndMerge");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky, "post-v2.2: non-sticky so any free worker can dispatch");
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        assert_eq!(p.output_queues.len(), 1);
        // Don't strictly assert the exact QueueSpec variant to avoid
        // brittle coupling to the framework's internal layout.
    }

    #[test]
    fn drop_kills_long_running_subprocess() {
        // Spawn a "sleep" subprocess; Drop should kill it within
        // AlignerProcess::Drop's 1-second deadline.
        let cfg = make_test_cfg();
        let start = std::time::Instant::now();
        {
            let _step = AlignAndMergeStep::new(cfg, "sleep 60").expect("spawn sleep");
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
    fn header_handle_is_set_after_clean_drain_on_empty_input() {
        // Spawn a command that emits a minimal BAM header (no
        // records). The reader's first action is to read the
        // header and resolve the handle.
        //
        // `samtools view -bH /dev/null` would do this, but we can't
        // rely on samtools in unit tests. Instead, use a hand-crafted
        // minimal BAM: BGZF-compressed `BAM\1` + 0 header text + 0
        // references + BGZF EOF. We emit this via `printf | gzip`
        // — wait, BAM uses BGZF specifically, not plain gzip.
        //
        // Easier path: just check the handle is NOT set when the
        // subprocess emits nothing (the reader fails to find BGZF
        // magic, returns Err, poisons the handle). This still
        // verifies the poison path.
        let cfg = make_test_cfg();
        let handle = cfg.header_handle.clone();

        // `true` exits immediately — no stdout. Reader sees empty
        // stdout, peek finds 0 bytes, is_bgzf=false → Err →
        // handle.poison.
        let step = AlignAndMergeStep::new(cfg, "true").expect("spawn true");

        // Drive the step a few times to let the reader thread run
        // and propagate. The exact threading interleaving means we
        // can't be deterministic — but we can give it a generous
        // window.
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

        // Subprocess + threads tear down at end-of-scope when
        // `step` is dropped.
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
    /// `FASTQ_WRITER_EXCLUDE_FLAGS`. AAM's writer should reject this
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
        use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED};
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
        // the `token_chan` construction site in `AlignAndMergeStep::new`).
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

        // `cat` self-exits once `writer_loop` drops its BufWriter
        // (and the wrapped `ChildStdin`) at function return —
        // closing the pipe makes `cat` see stdin EOF. The explicit
        // kill+wait below is belt-and-braces in case `writer_loop`
        // bailed before dropping (e.g. on an unexpected panic).
        let _ = cat.kill();
        let _ = cat.wait();
    }

    /// Regression for the AAM deadlock fixed in commit `069c265`.
    ///
    /// **Setup**: feed `writer_loop` more batches than the old
    /// `IN_CHAN_DEPTH=2` token-channel bound would have allowed,
    /// while keeping the reader side parked — `token_rx` is held
    /// but never receives. With the old bounded `token_chan`, the
    /// writer would block on its third `token_tx.send(...)` inside
    /// the `while let Ok(batch) = in_rx.recv()` body and never
    /// exit, causing `writer_loop` to hang here.
    ///
    /// **Invariant**: `token_chan` is unbounded; the writer must
    /// drain `in_rx` regardless of whether anyone is recv'ing tokens.
    /// The test runs `writer_loop` on a worker thread with a 30 s
    /// join timeout — a regression that re-introduces the bound
    /// shows up as a `join` timeout (test failure), not as silent
    /// flakiness.
    ///
    /// `cat > /dev/null` plays the role of bwa: it absorbs the
    /// FASTQ writes without emitting anything. With the unbounded
    /// `token_chan` the writer pushes 16 tokens and then exits
    /// cleanly when `in_tx` is dropped.
    #[test]
    fn writer_loop_does_not_deadlock_when_reader_does_not_drain_tokens() {
        use std::process::{Command, Stdio};
        use std::sync::mpsc::sync_channel;
        use std::time::{Duration, Instant};

        // `cat > /dev/null` consumes our FASTQ writes silently —
        // standing in for an aligner that buffers stdin without
        // emitting anything on stdout (the exact condition that
        // triggered the production deadlock on bwa with
        // `-K 150000000` and a small input).
        let mut cat = Command::new("cat")
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .expect("spawn cat");
        let cat_stdin = cat.stdin.take().expect("cat stdin");

        // Match production: in_chan small + bounded (mirrors
        // `IN_CHAN_DEPTH = 2`); token channel unbounded.
        // Pre-fix this would have deadlocked at batch 3; use 16 so a
        // partial regression (e.g. someone bumping the bound to 8 or
        // 16 without making it unbounded) is also caught.
        let n_batches: usize = 16;
        let (in_tx, in_rx) = sync_channel::<BamTemplateBatch>(2);
        let (token_tx, token_rx) = channel::<BatchToken>();
        let shared = Arc::new(SharedState::new());

        // Drive the writer on a worker thread so the test can apply a
        // join timeout — a deadlock surfaces as a timeout rather than
        // hanging the test runner.
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

        // Pump N batches of one paired template each from a separate
        // producer thread. If the writer deadlocks on
        // `token_tx.send` (the failure mode this test guards), the
        // bounded `in_chan` (depth 2) will also fill up and block
        // this producer — but the main test thread stays responsive
        // and surfaces the timeout below. Pushing the sends in-line
        // on the main thread would let the same `in_chan`-full
        // condition wedge the test runner itself.
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
                // Closing `in_tx` flips `in_rx.recv()` to `Err` after
                // the last buffered batch is consumed. The writer must
                // observe that and exit; with the old bounded
                // `token_chan`, the writer would be stuck on
                // `token_tx.send` for batch 3+ and never get back to
                // the `recv` site.
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

    fn make_sq_header(refs: &[(&str, usize)]) -> Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        let mut b = Header::builder();
        for (name, length) in refs {
            let map: Map<ReferenceSequence> = Map::<ReferenceSequence>::new(
                std::num::NonZeroUsize::new(*length).expect("nonzero ref length"),
            );
            b = b.add_reference_sequence(bstr::BString::from(*name), map);
        }
        b.build()
    }

    #[test]
    fn validate_sq_consistency_accepts_matching_refs() {
        let partial = make_sq_header(&[("chr1", 1000), ("chr2", 2000)]);
        let aligner = make_sq_header(&[("chr1", 1000), ("chr2", 2000)]);
        validate_sq_consistency(&partial, &aligner).expect("matching refs are ok");
    }

    #[test]
    fn validate_sq_consistency_rejects_count_mismatch() {
        let partial = make_sq_header(&[("chr1", 1000), ("chr2", 2000)]);
        let aligner = make_sq_header(&[("chr1", 1000)]);
        let err =
            validate_sq_consistency(&partial, &aligner).expect_err("count mismatch must reject");
        let msg = err.to_string();
        assert!(msg.contains("@SQ count") && msg.contains("does not match"));
    }

    #[test]
    fn validate_sq_consistency_rejects_length_mismatch() {
        let partial = make_sq_header(&[("chr1", 1000)]);
        let aligner = make_sq_header(&[("chr1", 999)]);
        let err =
            validate_sq_consistency(&partial, &aligner).expect_err("length mismatch must reject");
        let msg = err.to_string();
        assert!(msg.contains("length mismatch") && msg.contains("chr1"));
    }

    #[test]
    fn validate_sq_consistency_rejects_name_mismatch() {
        // Two refs so the position-index reporting can be exercised
        // (round-1 had a format-args bug that put the dict name in
        // the position slot; a single-ref test wouldn't have caught
        // it).
        let partial = make_sq_header(&[("chr1", 1000), ("chr2", 1000)]);
        let aligner = make_sq_header(&[("chr1", 1000), ("chrZ", 1000)]);
        let err =
            validate_sq_consistency(&partial, &aligner).expect_err("name mismatch must reject");
        let msg = err.to_string();
        assert!(msg.contains("name mismatch"));
        assert!(msg.contains("position 1"), "position index in message: {msg}");
        assert!(msg.contains("dict='chr2'"), "dict name in message: {msg}");
        assert!(msg.contains("aligner='chrZ'"), "aligner name in message: {msg}");
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
        // emits the file's bytes and exits when it's done. The
        // writer thread sits in `in_rx.recv()` until the step is
        // dropped at the end of the test.
        let cmd = format!("cat {}", fixture.display());
        let step = AlignAndMergeStep::new(cfg, &cmd).expect("spawn fake aligner");

        // Wait up to 5s for the reader thread to parse the header
        // and resolve the handle. The shell command takes a moment
        // to start.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while !handle.is_set() && std::time::Instant::now() < deadline {
            thread::sleep(std::time::Duration::from_millis(20));
        }
        assert!(handle.is_set(), "handle not resolved within 5s");
        let result = handle.try_get().expect("set");
        assert!(result.is_ok(), "header should resolve to Ok, not poison: {:?}", result.err());

        // Subprocess + threads tear down at end-of-scope when
        // `step` is dropped.
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
        let _step = AlignAndMergeStep::new(cfg, "true").expect("spawn true");

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
    /// Used as the output of a fake-aligner shell command so the
    /// reader can exercise the SAM-text code path through
    /// `peek_aligner_format` → `SamTemplateStream`.
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
        let _step = AlignAndMergeStep::new(cfg, &cmd).expect("spawn fake SAM aligner");

        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);
        while !handle.is_set() && std::time::Instant::now() < deadline {
            thread::sleep(std::time::Duration::from_millis(20));
        }
        assert!(handle.is_set(), "handle not resolved within 5s on SAM-text input");
        let result = handle.try_get().expect("set");
        assert!(result.is_ok(), "SAM-text header should resolve handle to Ok: {:?}", result.err(),);
    }

    /// `merge_aligner_header` must keep the *partial* header's `@PG` on a
    /// duplicate ID (the aligner's PG with that ID is dropped, not merged)
    /// while still appending aligner PGs whose IDs are unique to the
    /// aligner. This pins the dedup behavior documented on the function.
    #[test]
    fn merge_aligner_header_keeps_partial_pg_on_duplicate_id() {
        use bstr::BString;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::Program;
        use noodles::sam::header::record::value::map::program::tag as pg_tag;

        // Helper: build a @PG map carrying a distinguishing PN field so we
        // can tell the partial's PG apart from the aligner's.
        let program_with_name = |name: &str| -> Map<Program> {
            Map::<Program>::builder().insert(pg_tag::NAME, name).build().expect("valid @PG")
        };

        // Partial header: a "dup" PG (PN=partial) plus a partial-only PG.
        let partial = Header::builder()
            .add_program(BString::from("dup"), program_with_name("partial"))
            .add_program(BString::from("onlypartial"), program_with_name("partial"))
            .build();

        // Aligner header: a "dup" PG (PN=aligner, must be dropped) plus a
        // unique "bwa" PG (must be appended).
        let aligner = Header::builder()
            .add_program(BString::from("dup"), program_with_name("aligner"))
            .add_program(BString::from("bwa"), program_with_name("aligner"))
            .build();

        let merged = merge_aligner_header(&partial, &aligner);
        let programs = merged.programs();
        let programs = programs.as_ref();

        // Exactly three PGs survive: dup (partial's), onlypartial, bwa.
        assert_eq!(programs.len(), 3, "expected dup + onlypartial + bwa");

        // Read back the PN field for a given @PG ID.
        let pn_of = |id: &str| -> Option<String> {
            programs
                .get(&BString::from(id))
                .and_then(|pg| pg.other_fields().get(&pg_tag::NAME).map(ToString::to_string))
        };

        // The duplicate-ID PG is the PARTIAL's (PN=partial), not the
        // aligner's — the aligner's same-ID PG was dropped.
        assert_eq!(
            pn_of("dup").as_deref(),
            Some("partial"),
            "duplicate @PG ID must retain the partial's PG, not the aligner's"
        );
        // The partial-only PG is preserved.
        assert_eq!(pn_of("onlypartial").as_deref(), Some("partial"));
        // The aligner-unique PG is appended.
        assert_eq!(
            pn_of("bwa").as_deref(),
            Some("aligner"),
            "aligner @PG with a non-duplicate ID must be appended"
        );
    }
}
