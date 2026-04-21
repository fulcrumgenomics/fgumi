//! Unified progress tracker for fgumi.
//!
//! Every fgumi command emits progress events through this module. A dedicated
//! tracker thread receives events via a lock-free channel and renders them to
//! stderr based on the environment:
//!
//! - TTY + interactive → dashboard (indicatif bars + status line)
//! - Non-TTY or batch   → heartbeat (periodic logfmt lines) + end summary
//! - `Mode::None`       → tracker thread is not spawned; every producer call
//!   is a single branch-predicted atomic load, zero allocations
//!
//! Tasks 1-3 scope: module skeleton, Event enum, producer API, and the three
//! renderers (Dashboard, Heartbeat, Summary). Wiring into runall and the
//! standalone commands lands in Tasks 4-5.
//!
//! See `docs/superpowers/specs/2026-04-18-progress-tracker-design.md`.

// Clippy allowances for the rendering module:
//
// - `cast_precision_loss` / `cast_possible_truncation` / `cast_sign_loss`:
//   rate and ETA math intentionally uses f64 and then truncates back to u64
//   for display. The magnitudes involved (records/sec, seconds remaining)
//   never approach u64::MAX, so the lossy conversions are benign.
// - `format_push_string`: the summary builder assembles a single buffer with
//   multi-line format!() appends for readability; the alternative (write!
//   against a Formatter) is harder to follow.
// - `inline_always`: `with_sender` is the producer hot path; inlining it
//   into every caller removes the channel existence check on `Mode::None`.
// - `too_many_lines`: `tracker_loop` and friends orchestrate multiple event
//   types in a single state machine; splitting would hurt readability.
// - `map_unwrap_or`: idiomatic here for `Option<Instant>` time math.
// - `match_same_arms`: the cases genuinely represent the same behavior and
//   the arm bodies document the intent.
// - `doc_markdown`: documentation references `fgumi` / `tracker`; they are
//   not identifiers that need backticks.
// - `missing_panics_doc`: `init` only panics in a test helper path.
// - `needless_pass_by_value` on `tracker_loop`: the loop owns the receiver
//   for its lifetime; taking by value communicates ownership transfer.
// - `used_underscore_items`: `_send_error` / `_send_warning` are the macro
//   implementation details; the macros export a clean public surface.
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::format_push_string,
    clippy::inline_always,
    clippy::too_many_lines,
    clippy::map_unwrap_or,
    clippy::match_same_arms,
    clippy::doc_markdown,
    clippy::missing_panics_doc,
    clippy::needless_pass_by_value,
    clippy::used_underscore_items
)]

use std::collections::BTreeMap;
use std::io::IsTerminal;
use std::sync::{Arc, OnceLock};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use crossbeam_channel::{Sender, bounded};
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};

// ─── Public API: init / guard / mode ────────────────────────────────────

/// Rendering mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// TTY → Dashboard; non-TTY → Heartbeat. Default.
    Auto,
    /// Force dashboard (debug aid).
    Dashboard,
    /// Force heartbeat (for tailing logs even on a TTY).
    Heartbeat,
    /// No tracker. Every producer call is a no-op.
    None,
}

/// RAII handle. Dropping shuts down the tracker thread.
/// Idempotent: a second drop is a no-op.
pub struct TrackerGuard {
    handle: Option<JoinHandle<()>>,
    shutdown_tx: Option<Sender<Event>>,
}

impl Drop for TrackerGuard {
    fn drop(&mut self) {
        if let Some(tx) = self.shutdown_tx.take() {
            let _ = tx.send(Event::Shutdown);
        }
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

/// Spawn the tracker thread (or return a no-op guard).
///
/// Second invocation in the same process is a no-op: the global channel
/// is set exactly once. Calling `init(Mode::None)` leaves the channel
/// empty forever; every producer call in the whole binary becomes a no-op.
#[must_use = "drop the guard to shut down the tracker cleanly"]
pub fn init(mode: Mode) -> TrackerGuard {
    if mode == Mode::None {
        return TrackerGuard { handle: None, shutdown_tx: None };
    }

    let (tx, rx) = bounded::<Event>(10_000);

    if CHANNEL.set(tx.clone()).is_err() {
        return TrackerGuard { handle: None, shutdown_tx: None };
    }

    let handle = std::thread::Builder::new()
        .name("fgumi-progress".to_string())
        .spawn(move || tracker_loop(rx, mode))
        .expect("spawn fgumi-progress thread");

    TrackerGuard { handle: Some(handle), shutdown_tx: Some(tx) }
}

// ─── Internal state ──────────────────────────────────────────────────────

static CHANNEL: OnceLock<Sender<Event>> = OnceLock::new();

#[allow(dead_code)] // Used by the Task 2 producer API.
#[inline(always)]
pub(crate) fn with_sender(f: impl FnOnce(&Sender<Event>)) {
    if let Some(tx) = CHANNEL.get() {
        f(tx);
    }
}

// ─── Event enum ──────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
#[allow(dead_code)] // Variants are emitted starting from Task 2.
pub(crate) enum Event {
    PipelineStarted { name: Arc<str>, total: Option<u64> },
    PipelineFinished { ok: bool },
    StageStarted { stage: Arc<str>, total: Option<u64> },
    StageFinished { stage: Arc<str> },

    // Pipeline-level (drive header bar + ETA)
    RecordsRead { n: u64 },
    RecordsWritten { n: u64 },
    ReadsComplete { total: u64 },

    // Per-stage (drive per-stage bars + bottleneck)
    RecordsIn { stage: Arc<str>, n: u64 },
    RecordsOut { stage: Arc<str>, n: u64 },
    QueueDepth { stage: Arc<str>, bytes: u64 },

    Warning { stage: Arc<str>, msg: String },
    Error { stage: Arc<str>, msg: String },

    Shutdown,
}

// ─── Pipeline-level counters (drive header bar + ETA) ────────────────────

/// Record a batch of records consumed from the source (FASTQ or BAM reader).
/// Accumulates into the tracker's `total_read` counter.
#[inline]
pub fn records_read(n: u64) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::RecordsRead { n });
    });
}

/// Record a batch of records written to the final output. Drives the
/// pipeline's header progress bar.
#[inline]
pub fn records_written(n: u64) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::RecordsWritten { n });
    });
}

/// Called exactly once when the source reader reaches EOF on all input.
/// Unlocks the ETA computation: before this, ETA is rendered as `???`;
/// after, ETA = (total - records_written) / bottleneck_rate.
pub fn reads_complete(total: u64) {
    tracing::info!(total, "reads complete");
    with_sender(|tx| {
        let _ = tx.try_send(Event::ReadsComplete { total });
    });
}

// ─── Per-stage counters (drive per-stage bars + bottleneck) ──────────────

/// Record a batch of records entering `stage`.
#[inline]
pub fn records_in(stage: &str, n: u64) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::RecordsIn { stage: stage.into(), n });
    });
}

/// Record a batch of records leaving `stage`. Used to advance that stage's
/// progress bar in the dashboard.
#[inline]
pub fn records_out(stage: &str, n: u64) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::RecordsOut { stage: stage.into(), n });
    });
}

/// Report queue-held bytes for `stage` (snapshot). Used to detect the
/// current bottleneck.
#[inline]
pub fn queue_depth(stage: &str, bytes: u64) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::QueueDepth { stage: stage.into(), bytes });
    });
}

// ─── Lifecycle (tracker + tracing::info) ─────────────────────────────────

/// Emit once per fgumi invocation, before orchestration begins.
pub fn pipeline_started(name: &str, total: Option<u64>) {
    tracing::info!(pipeline = name, total = total, "pipeline started");
    with_sender(|tx| {
        let _ = tx.try_send(Event::PipelineStarted { name: name.into(), total });
    });
}

/// Emit once per fgumi invocation, after `run()` returns (ok or err).
pub fn pipeline_finished(ok: bool) {
    tracing::info!(ok, "pipeline finished");
    with_sender(|tx| {
        let _ = tx.try_send(Event::PipelineFinished { ok });
    });
}

/// Emit before `stage` begins processing. `total` is Some(n) when the
/// total record count is known ahead of time (usually None).
pub fn stage_started(stage: &str, total: Option<u64>) {
    tracing::info!(stage, total, "stage started");
    with_sender(|tx| {
        let _ = tx.try_send(Event::StageStarted { stage: stage.into(), total });
    });
}

/// Emit after `stage` has drained all input and closed its output.
pub fn stage_finished(stage: &str) {
    tracing::info!(stage, "stage finished");
    with_sender(|tx| {
        let _ = tx.try_send(Event::StageFinished { stage: stage.into() });
    });
}

// ─── Message helpers (invoked by the macros below) ───────────────────────

#[doc(hidden)]
pub fn _send_error(stage: &str, msg: String) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::Error { stage: stage.into(), msg });
    });
}

#[doc(hidden)]
pub fn _send_warning(stage: &str, msg: String) {
    with_sender(|tx| {
        let _ = tx.try_send(Event::Warning { stage: stage.into(), msg });
    });
}

// ─── Message macros ──────────────────────────────────────────────────────
//
// error!/warn!: flow to both tracing AND the tracker (surfaced in dashboard
// status line; counted in heartbeat and summary)
// info!/debug!/trace!: flow to tracing only

/// Emit an error. Flows to `tracing::error!` AND the tracker.
#[macro_export]
macro_rules! progress_error {
    (stage = $stage:expr, $($arg:tt)+) => {{
        let msg = format!($($arg)+);
        ::tracing::error!(stage = $stage, "{}", msg);
        $crate::progress::_send_error($stage, msg);
    }};
    ($($arg:tt)+) => {{
        ::tracing::error!($($arg)+);
    }};
}

/// Emit a warning. Flows to `tracing::warn!` AND the tracker when `stage` is set.
#[macro_export]
macro_rules! progress_warn {
    (stage = $stage:expr, $($arg:tt)+) => {{
        let msg = format!($($arg)+);
        ::tracing::warn!(stage = $stage, "{}", msg);
        $crate::progress::_send_warning($stage, msg);
    }};
    ($($arg:tt)+) => {{
        ::tracing::warn!($($arg)+);
    }};
}

/// Emit info-level message. Tracing only.
#[macro_export]
macro_rules! progress_info {
    (stage = $stage:expr, $($arg:tt)+) => {
        ::tracing::info!(stage = $stage, $($arg)+)
    };
    ($($arg:tt)+) => {
        ::tracing::info!($($arg)+)
    };
}

/// Emit debug-level message. Tracing only.
#[macro_export]
macro_rules! progress_debug {
    (stage = $stage:expr, $($arg:tt)+) => {
        ::tracing::debug!(stage = $stage, $($arg)+)
    };
    ($($arg:tt)+) => {
        ::tracing::debug!($($arg)+)
    };
}

/// Emit trace-level message. Tracing only.
#[macro_export]
macro_rules! progress_trace {
    (stage = $stage:expr, $($arg:tt)+) => {
        ::tracing::trace!(stage = $stage, $($arg)+)
    };
    ($($arg:tt)+) => {
        ::tracing::trace!($($arg)+)
    };
}

// Re-export macros at the module path so callers can write `progress::warn!`.
pub use progress_debug as debug;
pub use progress_error as error;
pub use progress_info as info;
pub use progress_trace as trace;
pub use progress_warn as warn;

// ─── Per-stage state ─────────────────────────────────────────────────────

/// Per-stage running state.
#[derive(Debug)]
struct StageState {
    name: Arc<str>,
    #[allow(dead_code)] // Reserved for future use (e.g. absolute progress).
    total: Option<u64>,
    records_in: u64,
    records_out: u64,
    queue_bytes: u64,
    /// Exponentially-weighted moving average of records_out per second.
    rate_ema: f64,
    last_rate_tick: Instant,
    last_rate_records: u64,
    bar: Option<ProgressBar>,
}

impl StageState {
    fn new(name: Arc<str>, total: Option<u64>) -> Self {
        Self {
            name,
            total,
            records_in: 0,
            records_out: 0,
            queue_bytes: 0,
            rate_ema: 0.0,
            last_rate_tick: Instant::now(),
            last_rate_records: 0,
            bar: None,
        }
    }

    /// Update the EMA rate if at least 0.5s has elapsed since last tick.
    fn update_rate(&mut self, now: Instant) {
        let dt = now.duration_since(self.last_rate_tick).as_secs_f64();
        if dt >= 0.5 {
            let delta = self.records_out.saturating_sub(self.last_rate_records) as f64;
            let current = delta / dt;
            // EMA with alpha=0.3 — responsive but smoothed.
            self.rate_ema = 0.7 * self.rate_ema + 0.3 * current;
            self.last_rate_tick = now;
            self.last_rate_records = self.records_out;
        }
    }
}

// ─── Tracker state ───────────────────────────────────────────────────────

/// Rendering mode after resolving `Mode::Auto` against the actual environment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EffectiveMode {
    Dashboard,
    Heartbeat,
}

/// Tracker-thread-local state.
struct TrackerState {
    #[allow(dead_code)] // Requested mode; the effective mode drives rendering.
    mode: Mode,
    effective: EffectiveMode,
    pipeline_name: Option<Arc<str>>,
    pipeline_started_at: Option<Instant>,
    // Pipeline-level totals (drive header bar + ETA).
    records_read: u64,
    records_written: u64,
    /// Set by `ReadsComplete`; enables ETA computation.
    total_read: Option<u64>,
    header_bar: Option<ProgressBar>,
    // Per-stage state.
    stages: BTreeMap<Arc<str>, StageState>,
    /// Insertion order for deterministic dashboard rendering.
    stage_order: Vec<Arc<str>>,
    warnings: u64,
    errors: u64,
    last_warning: Option<(Arc<str>, String, Instant)>,
    last_error: Option<(Arc<str>, String, Instant)>,
    multi: Option<MultiProgress>,
    last_heartbeat: Instant,
    heartbeat_interval: Duration,
}

impl TrackerState {
    /// Compute ETA. Returns `None` before `reads_complete` fires, or before a
    /// bottleneck stage has accumulated enough rate signal.
    fn eta(&self) -> Option<Duration> {
        let total = self.total_read?;
        let bn = self.stages.values().max_by_key(|s| s.queue_bytes)?;
        if bn.rate_ema < 1.0 {
            return None;
        }
        let remaining = total.saturating_sub(self.records_written) as f64;
        Some(Duration::from_secs_f64(remaining / bn.rate_ema))
    }
}

/// Resolve the requested mode into the concrete rendering mode.
fn resolve_mode(requested: Mode) -> EffectiveMode {
    match requested {
        Mode::Dashboard => EffectiveMode::Dashboard,
        Mode::Heartbeat => EffectiveMode::Heartbeat,
        Mode::Auto => {
            if std::io::stderr().is_terminal() {
                EffectiveMode::Dashboard
            } else {
                EffectiveMode::Heartbeat
            }
        }
        Mode::None => unreachable!("resolve_mode should not be called for Mode::None"),
    }
}

// ─── Tracker thread loop ──────────────────────────────────────────────────

fn tracker_loop(rx: crossbeam_channel::Receiver<Event>, mode: Mode) {
    let effective = resolve_mode(mode);
    let multi = if effective == EffectiveMode::Dashboard {
        let m = MultiProgress::new();
        m.set_draw_target(ProgressDrawTarget::stderr_with_hz(10));
        Some(m)
    } else {
        None
    };

    let mut state = TrackerState {
        mode,
        effective,
        pipeline_name: None,
        pipeline_started_at: None,
        records_read: 0,
        records_written: 0,
        total_read: None,
        header_bar: None,
        stages: BTreeMap::new(),
        stage_order: Vec::new(),
        warnings: 0,
        errors: 0,
        last_warning: None,
        last_error: None,
        multi,
        last_heartbeat: Instant::now(),
        heartbeat_interval: Duration::from_secs(30),
    };

    // Emit an initial heartbeat 5s after start (non-TTY only).
    let first_heartbeat = if effective == EffectiveMode::Heartbeat {
        crossbeam_channel::after(Duration::from_secs(5))
    } else {
        crossbeam_channel::never()
    };
    let mut first_heartbeat_done = false;
    let ticker = crossbeam_channel::tick(Duration::from_millis(500));

    loop {
        crossbeam_channel::select! {
            recv(rx) -> ev => {
                match ev {
                    Ok(Event::Shutdown) => break,
                    Ok(ev) => handle_event(&mut state, ev),
                    Err(_) => break,
                }
            }
            recv(ticker) -> _ => {
                tick(&mut state);
            }
            recv(first_heartbeat) -> _ => {
                if !first_heartbeat_done && state.effective == EffectiveMode::Heartbeat {
                    emit_heartbeat(&state);
                    first_heartbeat_done = true;
                }
            }
        }
    }

    // Drain any remaining events (short deadline) before the final summary.
    let deadline = Instant::now() + Duration::from_millis(100);
    while Instant::now() < deadline {
        match rx.try_recv() {
            Ok(Event::Shutdown) => break,
            Ok(ev) => handle_event(&mut state, ev),
            Err(_) => break,
        }
    }

    emit_summary(&state);
}

// ─── Event handling ──────────────────────────────────────────────────────

fn handle_event(state: &mut TrackerState, ev: Event) {
    let now = Instant::now();
    match ev {
        Event::PipelineStarted { name, total: _ } => {
            state.pipeline_name = Some(name.clone());
            state.pipeline_started_at = Some(now);
            // Install a header spinner at the top. It becomes a real progress
            // bar once `ReadsComplete` fires.
            if let Some(multi) = &state.multi {
                let pb = multi.add(ProgressBar::new_spinner());
                pb.set_style(
                    ProgressStyle::with_template("{prefix:14} {spinner} {msg}")
                        .expect("static template compiles"),
                );
                pb.set_prefix((*name).to_string());
                pb.enable_steady_tick(Duration::from_millis(120));
                state.header_bar = Some(pb);
            }
        }
        Event::PipelineFinished { ok: _ } => {
            // Finalize any stage bars still in flight.
            let names: Vec<Arc<str>> = state.stage_order.clone();
            for name in names {
                if let Some(s) = state.stages.get(&name)
                    && let Some(bar) = &s.bar
                {
                    bar.finish();
                }
            }
            if let Some(bar) = &state.header_bar {
                bar.finish();
            }
        }
        Event::StageStarted { stage, total } => {
            if !state.stages.contains_key(&stage) {
                state.stage_order.push(stage.clone());
            }
            let mut ss = StageState::new(stage.clone(), total);
            if state.effective == EffectiveMode::Dashboard
                && let Some(multi) = &state.multi
            {
                let style =
                    ProgressStyle::with_template("  {prefix:14} [{bar:20}] {pos}/{len} {msg}")
                        .expect("static template compiles")
                        .progress_chars("██▓");
                let pb = if let Some(t) = total {
                    multi.add(ProgressBar::new(t))
                } else {
                    multi.add(ProgressBar::new_spinner())
                };
                pb.set_style(style);
                pb.set_prefix((*stage).to_string());
                ss.bar = Some(pb);
            }
            state.stages.insert(stage, ss);
        }
        Event::StageFinished { stage } => {
            if let Some(s) = state.stages.get(&stage)
                && let Some(bar) = &s.bar
            {
                bar.finish();
            }
        }
        Event::RecordsRead { n } => {
            state.records_read = state.records_read.saturating_add(n);
            // Before reads_complete: header is a spinner; show read count in msg.
            if state.total_read.is_none()
                && let Some(bar) = &state.header_bar
            {
                bar.set_message(format!("{} read", format_count(state.records_read)));
            }
        }
        Event::RecordsWritten { n } => {
            state.records_written = state.records_written.saturating_add(n);
            if let Some(bar) = &state.header_bar {
                // After ReadsComplete, this advances the real bar.
                // Before ReadsComplete, this is a no-op on the spinner.
                bar.set_position(state.records_written);
            }
        }
        Event::ReadsComplete { total } => {
            state.total_read = Some(total);
            // Rebuild the header as a real progress bar now that the total is known.
            if let (Some(multi), Some(old)) = (&state.multi, state.header_bar.take()) {
                old.finish_and_clear();
                let style = ProgressStyle::with_template(
                    "{prefix:14} [{bar:30}] {pos}/{len}  ETA {eta_precise}",
                )
                .expect("static template compiles")
                .progress_chars("██▓");
                let pb = multi.insert(0, ProgressBar::new(total));
                pb.set_style(style);
                pb.set_prefix("overall");
                pb.set_position(state.records_written);
                state.header_bar = Some(pb);
            }
            if let Some(multi) = &state.multi {
                multi.suspend(|| {
                    tracing::info!(total, "reads complete; ETA now computable");
                });
            } else {
                tracing::info!(total, "reads complete; ETA now computable");
            }
        }
        Event::RecordsIn { stage, n } => {
            if let Some(s) = state.stages.get_mut(&stage) {
                s.records_in = s.records_in.saturating_add(n);
            }
        }
        Event::RecordsOut { stage, n } => {
            if let Some(s) = state.stages.get_mut(&stage) {
                s.records_out = s.records_out.saturating_add(n);
                if let Some(bar) = &s.bar {
                    bar.set_position(s.records_out);
                }
            }
        }
        Event::QueueDepth { stage, bytes } => {
            if let Some(s) = state.stages.get_mut(&stage) {
                s.queue_bytes = bytes;
            }
        }
        Event::Warning { stage, msg } => {
            state.warnings = state.warnings.saturating_add(1);
            state.last_warning = Some((stage, msg, now));
        }
        Event::Error { stage, msg } => {
            state.errors = state.errors.saturating_add(1);
            state.last_error = Some((stage, msg, now));
        }
        Event::Shutdown => unreachable!("handled in tracker_loop"),
    }
}

fn tick(state: &mut TrackerState) {
    let now = Instant::now();
    for s in state.stages.values_mut() {
        s.update_rate(now);
    }
    if state.effective == EffectiveMode::Heartbeat
        && now.duration_since(state.last_heartbeat) >= state.heartbeat_interval
    {
        emit_heartbeat(state);
        state.last_heartbeat = now;
    }
    if state.effective == EffectiveMode::Dashboard {
        update_dashboard_status_line(state);
    }
}

fn bottleneck_stage(state: &TrackerState) -> Option<&Arc<str>> {
    state.stages.values().max_by_key(|s| s.queue_bytes).map(|s| &s.name)
}

/// Refresh the per-stage bar messages with rate + bottleneck marker.
fn update_dashboard_status_line(state: &TrackerState) {
    let bn = bottleneck_stage(state).cloned();
    for s in state.stages.values() {
        if let Some(bar) = &s.bar {
            let rate_k = (s.rate_ema as u64) / 1_000;
            let is_bn = bn.as_ref().is_some_and(|b| b == &s.name);
            let msg = if is_bn {
                format!("{rate_k}k/s  🚦 bottleneck")
            } else {
                format!("{rate_k}k/s")
            };
            bar.set_message(msg);
        }
    }
}

/// Emit a logfmt heartbeat (one `runall.heartbeat` line + one per stage).
fn emit_heartbeat(state: &TrackerState) {
    let elapsed = state
        .pipeline_started_at
        .map(|t| Instant::now().duration_since(t).as_secs_f64())
        .unwrap_or(0.0);
    let mem_peak: u64 = state.stages.values().map(|s| s.queue_bytes).sum();
    let eta_s: Option<u64> = state.eta().map(|d| d.as_secs());

    let write = || {
        tracing::info!(
            pipeline = state.pipeline_name.as_deref().unwrap_or(""),
            elapsed_s = elapsed,
            read = state.records_read,
            written = state.records_written,
            total_read = state.total_read,
            eta_s = eta_s,
            mem_bytes = mem_peak,
            warnings = state.warnings,
            errors = state.errors,
            "runall.heartbeat"
        );

        let bn = bottleneck_stage(state).cloned();
        for s in state.stages.values() {
            let is_bn = bn.as_ref().is_some_and(|b| b == &s.name);
            tracing::info!(
                stage = &*s.name,
                records_in = s.records_in,
                records_out = s.records_out,
                rate = s.rate_ema as u64,
                queued_bytes = s.queue_bytes,
                bottleneck = is_bn,
                "runall.stage"
            );
        }
    };

    if let Some(multi) = &state.multi {
        multi.suspend(write);
    } else {
        write();
    }
}

/// Emit a human-readable summary table to stderr.
fn emit_summary(state: &TrackerState) {
    use std::io::Write;

    let elapsed = state
        .pipeline_started_at
        .map(|t| Instant::now().duration_since(t))
        .unwrap_or(Duration::ZERO);
    let pipeline = state.pipeline_name.as_deref().unwrap_or("fgumi");

    // Build the summary in a single buffer then write atomically. If the
    // dashboard is still active, suspend it so bars don't interleave.
    let mut buf = String::new();
    buf.push('\n');
    buf.push_str(&format!("{} completed in {}\n\n", pipeline, format_duration(elapsed)));

    let ratio = if state.records_written > 0 {
        state.records_read as f64 / state.records_written as f64
    } else {
        0.0
    };
    buf.push_str(&format!(
        "  {} records read     {} records written  (ratio {:.2}x)\n\n",
        format_count(state.records_read),
        format_count(state.records_written),
        ratio,
    ));
    buf.push_str("  Stage         records in      records out    rate\n");

    let bn = bottleneck_stage(state).cloned();
    for stage in &state.stage_order {
        if let Some(s) = state.stages.get(stage) {
            let is_bn = bn.as_ref().is_some_and(|b| b == stage);
            let bn_mark = if is_bn { "  <- bottleneck" } else { "" };
            buf.push_str(&format!(
                "  {:12}  {:>12}  {:>12}    {}/s{}\n",
                stage,
                s.records_in,
                s.records_out,
                format_rate(s.rate_ema as u64),
                bn_mark
            ));
        }
    }
    buf.push_str(&format!("\n  {} warnings  {} errors\n\n", state.warnings, state.errors));

    let write_buf = || {
        let mut out = std::io::stderr().lock();
        let _ = out.write_all(buf.as_bytes());
        let _ = out.flush();
    };

    if let Some(multi) = &state.multi {
        multi.suspend(write_buf);
    } else {
        write_buf();
    }
}

// ─── Formatting helpers ──────────────────────────────────────────────────

fn format_duration(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    if h > 0 {
        format!("{h}:{m:02}:{s:02}")
    } else if m > 0 {
        format!("{m}:{s:02}")
    } else {
        format!("{s}s")
    }
}

fn format_rate(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{}k", n / 1_000)
    } else {
        format!("{n}")
    }
}

fn format_count(n: u64) -> String {
    if n >= 1_000_000_000 {
        format!("{:.2}B", n as f64 / 1e9)
    } else if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1e6)
    } else if n >= 1_000 {
        format!("{:.0}k", n as f64 / 1e3)
    } else {
        format!("{n}")
    }
}

// ============================================================================
// Legacy `ProgressTracker` — simple atomic-counter + EMA rate estimator.
//
// Retained for standalone commands (`fgumi filter`, `fgumi dedup`, etc.) that
// log progress at interval boundaries rather than emitting structured events
// for the runall tracker. The event-based API above (`records_read`,
// `stage_started`, etc.) is the preferred surface for runall-orchestrated
// work; this tracker is for single-command callers that want simple periodic
// "processed N records" logging with optional ETA.
// ============================================================================

use std::sync::Mutex as LegacyMutex;
use std::sync::atomic::{AtomicU64 as LegacyAtomicU64, Ordering as LegacyOrdering};
use std::time::Instant as LegacyInstant;

/// Smoothing constant for the EMA rate estimator.
/// 0.3 balances responsiveness to rate changes with stability (same default as tqdm).
const EMA_ALPHA: f64 = 0.3;

struct EmaState {
    smoothed_rate: f64,
    calls: u32,
    last_count: u64,
    last_time: LegacyInstant,
}

impl EmaState {
    fn new() -> Self {
        Self { smoothed_rate: 0.0, calls: 0, last_count: 0, last_time: LegacyInstant::now() }
    }

    fn update(&mut self, current_count: u64) -> f64 {
        if current_count <= self.last_count {
            return self.corrected_rate();
        }
        let now = LegacyInstant::now();
        let dt = now.duration_since(self.last_time).as_secs_f64();
        if dt > 0.0 {
            #[allow(clippy::cast_precision_loss)]
            let dn = (current_count - self.last_count) as f64;
            let instantaneous_rate = dn / dt;
            self.smoothed_rate =
                EMA_ALPHA * instantaneous_rate + (1.0 - EMA_ALPHA) * self.smoothed_rate;
            self.calls += 1;
            self.last_count = current_count;
            self.last_time = now;
        }
        self.corrected_rate()
    }

    fn corrected_rate(&self) -> f64 {
        if self.calls == 0 {
            return 0.0;
        }
        let beta = 1.0 - EMA_ALPHA;
        let correction = 1.0 - beta.powi(self.calls.cast_signed());
        if correction <= 0.0 { 0.0 } else { self.smoothed_rate / correction }
    }
}

fn fmt_duration(secs: f64) -> String {
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    format_duration(std::time::Duration::from_secs(secs.round() as u64))
}

/// Thread-safe progress tracker that logs at interval boundaries.
///
/// Used by standalone commands to periodically emit "Processed N records"
/// messages. Optional total enables percentage + ETA via an EMA rate.
pub struct ProgressTracker {
    interval: u64,
    message: String,
    count: LegacyAtomicU64,
    total: Option<u64>,
    start_time: LegacyInstant,
    ema: LegacyMutex<EmaState>,
}

impl ProgressTracker {
    /// Create a tracker with default interval of 10,000.
    #[must_use]
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            interval: 10_000,
            message: message.into(),
            count: LegacyAtomicU64::new(0),
            total: None,
            start_time: LegacyInstant::now(),
            ema: LegacyMutex::new(EmaState::new()),
        }
    }

    /// Set the logging interval.
    #[must_use]
    pub fn with_interval(mut self, interval: u64) -> Self {
        self.interval = interval;
        self
    }

    /// Set the total expected count (enables percentage + ETA).
    #[must_use]
    pub fn with_total(mut self, total: u64) -> Self {
        self.total = (total > 0).then_some(total);
        self
    }

    /// Add to the count and log at crossed interval boundaries.
    /// Returns true if the final count is exactly on an interval.
    #[allow(clippy::cast_precision_loss)]
    pub fn log_if_needed(&self, additional: u64) -> bool {
        if additional == 0 {
            let count = self.count.load(LegacyOrdering::Relaxed);
            return count > 0 && count.is_multiple_of(self.interval);
        }
        let prev = self.count.fetch_add(additional, LegacyOrdering::Relaxed);
        let new_count = prev + additional;
        let prev_intervals = prev / self.interval;
        let new_intervals = new_count / self.interval;
        if new_intervals > prev_intervals {
            let rate = if self.total.is_some() {
                if let Ok(mut ema) = self.ema.lock() { ema.update(new_count) } else { 0.0 }
            } else {
                0.0
            };
            for i in (prev_intervals + 1)..=new_intervals {
                let milestone = i * self.interval;
                if let Some(total) = self.total {
                    let pct = (milestone as f64 / total as f64) * 100.0;
                    let eta_suffix = if rate > 0.0 {
                        let remaining = total.saturating_sub(milestone) as f64;
                        format!(", ETA ~{}", fmt_duration(remaining / rate))
                    } else {
                        String::new()
                    };
                    tracing::info!(
                        "{} {} / {} ({:.1}%{})",
                        self.message,
                        milestone,
                        total,
                        pct,
                        eta_suffix,
                    );
                } else {
                    tracing::info!("{} {}", self.message, milestone);
                }
            }
        }
        new_count.is_multiple_of(self.interval)
    }

    /// Log the final count. Emits even if off-interval when a total is set.
    pub fn log_final(&self) {
        let count = self.count.load(LegacyOrdering::Relaxed);
        if count == 0 && self.total.is_none() {
            return;
        }
        if self.total.is_some() {
            let elapsed = self.start_time.elapsed().as_secs_f64();
            tracing::info!("{} {} (complete, {})", self.message, count, fmt_duration(elapsed));
        } else if !self.log_if_needed(0) {
            tracing::info!("{} {} (complete)", self.message, count);
        }
    }

    /// Current count.
    #[must_use]
    pub fn count(&self) -> u64 {
        self.count.load(LegacyOrdering::Relaxed)
    }
}
// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
use std::sync::Mutex;

// Serialize tests that touch the global CHANNEL — OnceLock can't reset.
#[cfg(test)]
static TEST_LOCK: Mutex<()> = Mutex::new(());

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_none_is_no_op() {
        let _lock = TEST_LOCK.lock().unwrap();
        let _g = init(Mode::None);
    }

    #[test]
    fn init_auto_spawns_and_shuts_down_cleanly() {
        let _lock = TEST_LOCK.lock().unwrap();
        let g = init(Mode::Auto);
        drop(g);
    }

    #[test]
    fn tracker_state_eta_none_before_reads_complete() {
        // Can't easily test against a real tracker thread without spawning it,
        // so test the TrackerState::eta() method against a fabricated state.
        // (TrackerState is private, so this test must live inside the module.)
        let state = TrackerState {
            mode: Mode::Heartbeat,
            effective: EffectiveMode::Heartbeat,
            pipeline_name: None,
            pipeline_started_at: None,
            records_read: 1000,
            records_written: 500,
            total_read: None, // not complete yet
            header_bar: None,
            stages: BTreeMap::new(),
            stage_order: vec![],
            warnings: 0,
            errors: 0,
            last_warning: None,
            last_error: None,
            multi: None,
            last_heartbeat: std::time::Instant::now(),
            heartbeat_interval: std::time::Duration::from_secs(30),
        };
        assert!(state.eta().is_none(), "ETA must be None before reads_complete");
    }

    #[test]
    fn tracker_state_eta_some_after_reads_complete_with_rate() {
        let mut state = TrackerState {
            mode: Mode::Heartbeat,
            effective: EffectiveMode::Heartbeat,
            pipeline_name: None,
            pipeline_started_at: None,
            records_read: 10_000,
            records_written: 5_000,
            total_read: Some(10_000),
            header_bar: None,
            stages: BTreeMap::new(),
            stage_order: vec![],
            warnings: 0,
            errors: 0,
            last_warning: None,
            last_error: None,
            multi: None,
            last_heartbeat: std::time::Instant::now(),
            heartbeat_interval: std::time::Duration::from_secs(30),
        };
        // Insert a bottleneck stage with a known rate.
        let name: Arc<str> = Arc::from("sort");
        let mut ss = StageState::new(name.clone(), None);
        ss.queue_bytes = 1024; // pick as bottleneck
        ss.rate_ema = 1000.0; // 1000 records/sec
        state.stages.insert(name, ss);
        let eta = state.eta().expect("eta should compute");
        // 5000 remaining records / 1000 per sec = 5 seconds
        assert_eq!(eta.as_secs(), 5);
    }

    #[test]
    fn tracker_state_eta_none_with_zero_rate() {
        let mut state = TrackerState {
            mode: Mode::Heartbeat,
            effective: EffectiveMode::Heartbeat,
            pipeline_name: None,
            pipeline_started_at: None,
            records_read: 10_000,
            records_written: 5_000,
            total_read: Some(10_000),
            header_bar: None,
            stages: BTreeMap::new(),
            stage_order: vec![],
            warnings: 0,
            errors: 0,
            last_warning: None,
            last_error: None,
            multi: None,
            last_heartbeat: std::time::Instant::now(),
            heartbeat_interval: std::time::Duration::from_secs(30),
        };
        let name: Arc<str> = Arc::from("sort");
        let mut ss = StageState::new(name.clone(), None);
        ss.queue_bytes = 1024;
        ss.rate_ema = 0.5; // below 1.0 threshold
        state.stages.insert(name, ss);
        assert!(state.eta().is_none(), "ETA should be None when rate < 1.0");
    }

    #[test]
    fn producer_calls_are_no_op_without_init() {
        let _lock = TEST_LOCK.lock().unwrap();
        // Before init(Mode::Auto) runs in some other test, CHANNEL may be
        // unset. If it's been set by a prior test, these calls become
        // real sends (the tracker is a zombie by then since its guard was
        // dropped). Either way, these must not panic.
        records_read(100);
        records_written(50);
        records_in("sort", 1);
        records_out("sort", 2);
        queue_depth("sort", 1024);
        pipeline_started("test", Some(100));
        stage_started("sort", Some(100));
        stage_finished("sort");
        reads_complete(100);
        pipeline_finished(true);
        _send_warning("sort", "probe".to_string());
        _send_error("sort", "probe".to_string());
    }
}
