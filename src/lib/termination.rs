//! Process-wide handling of the termination signals SIGINT, SIGTERM and SIGHUP.
//!
//! Without a handler, those signals kill fgumi with no destructors run: sort
//! temp directories stay on disk, staged outputs are left beside their
//! destination, and a subprocess aligner (which runs in its own process group,
//! so a terminal's Ctrl-C does not reach it) keeps running until it next
//! touches a closed pipe.
//!
//! [`install`], called once from `main`, routes those signals into the
//! cancellation fgumi already has for failures:
//!
//! 1. Code that holds state needing cleanup opens a `Scope` *before* creating
//!    that state, and attaches what a termination request should invoke with
//!    `Scope::on_cancel` (for a pipeline, its `CancelHandle::cancel`).
//! 2. On the first signal, if any scope is open, every attached target is
//!    invoked and the request is remembered; a target attached later is invoked
//!    as it is attached, so a signal that lands while a chain is being built is
//!    not lost. With no scope open the process dies by the signal at once, as it
//!    did before this module existed.
//! 3. The command unwinds through its normal error path, running every
//!    destructor (the aligner's group is stopped, temp directories and
//!    unpublished staged outputs are removed). `main` hands the result to
//!    [`conclude`], which restores the signal's default action and re-raises it,
//!    so the parent sees a real death by that signal, not an exit code.
//!
//! Two escape hatches keep fgumi killable if cleanup stalls: a second signal
//! arriving at least `SECOND_SIGNAL_DEBOUNCE` after the first terminates at
//! once (from signal-hook's own handler, so no fgumi thread needs to run), and a
//! watchdog terminates with the first signal `DEFAULT_GRACE` after it was
//! received (`FGUMI_SIGNAL_GRACE_SECS` overrides the bound; `0` disables all of
//! this handling). Signals ignored when fgumi started (`nohup`) stay ignored.
//!
//! Invariant: once a request is accepted, the `fgumi-signal` thread only
//! sleeps, reads signal-hook's pipe, takes the coordinator's leaf mutex and
//! stores atomics. It never logs, joins or calls a cancel target (those run on
//! a one-shot thread), so a blocked stderr or a wedged worker cannot disable the
//! watchdog.
//!
//! Library callers never call [`install`]; their scopes are inert registry
//! entries and signals keep their default action.

use std::fmt;
use std::sync::atomic::{AtomicBool, AtomicI32, Ordering};
use std::sync::{Arc, LazyLock};
use std::time::{Duration, Instant};

use parking_lot::Mutex;

/// How long cleanup may take after the first signal before the watchdog
/// terminates the process. Above the aligner teardown's ~4.5 s worst case and
/// below the 10 s a `docker stop` waits before sending SIGKILL.
pub(crate) const DEFAULT_GRACE: Duration = Duration::from_secs(8);

/// A second signal arriving sooner than this after the first is treated as a
/// duplicate of it (GNU `timeout`, a closing terminal and signal-forwarding
/// wrappers all deliver one request twice within microseconds), not as a
/// request to stop immediately.
pub(crate) const SECOND_SIGNAL_DEBOUNCE: Duration = Duration::from_millis(500);

/// Environment variable overriding [`DEFAULT_GRACE`], in seconds; `0` turns
/// signal handling off.
pub(crate) const GRACE_ENV: &str = "FGUMI_SIGNAL_GRACE_SECS";

/// How often the signal thread checks the watchdog while a request is pending.
const TICK: Duration = Duration::from_millis(10);

/// The signals fgumi turns into cancellation.
#[cfg(unix)]
const HANDLED: [i32; 3] =
    [signal_hook::consts::SIGINT, signal_hook::consts::SIGTERM, signal_hook::consts::SIGHUP];

/// Marker error: the work stopped because a termination signal cancelled it.
///
/// Attach it as the source of a cancelled command's error so [`conclude`] can
/// tell a cancellation (reported once, when the signal arrived) from a real
/// failure (printed before the signal death).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Cancelled;

impl fmt::Display for Cancelled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("cancelled")
    }
}

impl std::error::Error for Cancelled {}

/// Whether `err`, or anything in its chain, is a [`Cancelled`].
#[must_use]
pub fn is_cancellation(err: &anyhow::Error) -> bool {
    err.chain().any(<dyn std::error::Error>::is::<Cancelled>)
}

/// What a termination request invokes: must be idempotent and must not block.
type CancelTarget = Arc<dyn Fn() + Send + Sync>;

// ---------------------------------------------------------------------------
// The coordinator: a pure state machine (no threads, signals or I/O).
// ---------------------------------------------------------------------------

/// Where the process is in handling a termination request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    /// No request accepted yet.
    Idle,
    /// A request is being honoured: cleanup is running.
    Cancelling {
        /// The first accepted signal (the one the process dies by).
        signal: i32,
        /// When it was accepted.
        since: Instant,
        /// When the watchdog terminates the process (`None`: a grace too
        /// large to represent, so no watchdog).
        deadline: Option<Instant>,
        /// Whether the handler-level second-signal hatch is armed.
        second_armed: bool,
    },
    /// Default termination by `signal` is under way.
    Forced { signal: i32 },
}

/// What the signal thread must do after a signal.
enum SignalAction {
    /// Die by this signal now, without cleanup.
    Terminate(i32),
    /// Invoke these targets (the first accepted signal).
    Cancel(Vec<CancelTarget>),
    /// Nothing (a duplicate inside the debounce).
    Ignore,
}

/// What the signal thread must do on a watchdog tick.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TickAction {
    Wait,
    /// The debounce has elapsed: let a further signal terminate from its handler.
    ArmSecondSignal,
    /// The grace period has elapsed: die by this signal.
    Terminate(i32),
}

struct Inner {
    phase: Phase,
    grace: Duration,
    /// Open scopes, by id, with their target once attached.
    scopes: Vec<(u64, Option<CancelTarget>)>,
    next_id: u64,
}

/// The registry of open scopes plus the termination state machine.
///
/// The mutex is a leaf: targets are cloned out and called after it is
/// released, so no foreign code ever runs under it.
struct Coordinator {
    inner: Mutex<Inner>,
    /// The accepted signal (`0`: none), mirrored from `inner.phase` on every
    /// transition so [`Coordinator::requested`] (polled per record) is one
    /// atomic load.
    pending: AtomicI32,
}

impl Coordinator {
    fn new(grace: Duration) -> Self {
        Self {
            inner: Mutex::new(Inner { phase: Phase::Idle, grace, scopes: Vec::new(), next_id: 0 }),
            pending: AtomicI32::new(0),
        }
    }

    fn set_grace(&self, grace: Duration) {
        self.inner.lock().grace = grace;
    }

    /// The accepted signal, if a request is pending.
    fn requested(&self) -> Option<i32> {
        match self.pending.load(Ordering::SeqCst) {
            0 => None,
            signal => Some(signal),
        }
    }

    /// Move to `phase`, keeping `pending` in step. Called with the lock held.
    fn set_phase(&self, inner: &mut Inner, phase: Phase) {
        inner.phase = phase;
        let signal = match phase {
            Phase::Idle => 0,
            Phase::Cancelling { signal, .. } | Phase::Forced { signal } => signal,
        };
        self.pending.store(signal, Ordering::SeqCst);
    }

    /// Register a new scope with no target, returning its id.
    fn enter(&self) -> u64 {
        let mut inner = self.inner.lock();
        let id = inner.next_id;
        inner.next_id += 1;
        inner.scopes.push((id, None));
        id
    }

    /// Attach (or replace) scope `id`'s target. Returns whether a request is
    /// already pending, in which case the caller must invoke `target` itself:
    /// the snapshot the signal took at acceptance did not include it.
    fn set_target(&self, id: u64, target: CancelTarget) -> bool {
        let mut inner = self.inner.lock();
        if let Some(entry) = inner.scopes.iter_mut().find(|(scope_id, _)| *scope_id == id) {
            entry.1 = Some(target);
        }
        !matches!(inner.phase, Phase::Idle)
    }

    /// Deregister scope `id`. The phase is untouched: once accepted, a request
    /// ends only in [`conclude`] or by the watchdog.
    fn exit(&self, id: u64) {
        self.inner.lock().scopes.retain(|(scope_id, _)| *scope_id != id);
    }

    fn on_signal(&self, signal: i32, now: Instant) -> SignalAction {
        let mut inner = self.inner.lock();
        match inner.phase {
            Phase::Idle if inner.scopes.is_empty() || inner.grace.is_zero() => {
                self.set_phase(&mut inner, Phase::Forced { signal });
                SignalAction::Terminate(signal)
            }
            Phase::Idle => {
                let deadline = now.checked_add(inner.grace);
                self.set_phase(
                    &mut inner,
                    Phase::Cancelling { signal, since: now, deadline, second_armed: false },
                );
                SignalAction::Cancel(
                    inner.scopes.iter().filter_map(|(_, target)| target.clone()).collect(),
                )
            }
            Phase::Cancelling { since, .. }
                if now.saturating_duration_since(since) < SECOND_SIGNAL_DEBOUNCE =>
            {
                SignalAction::Ignore
            }
            Phase::Cancelling { .. } | Phase::Forced { .. } => {
                self.set_phase(&mut inner, Phase::Forced { signal });
                SignalAction::Terminate(signal)
            }
        }
    }

    fn on_tick(&self, now: Instant) -> TickAction {
        let mut inner = self.inner.lock();
        match &mut inner.phase {
            Phase::Forced { signal } => TickAction::Terminate(*signal),
            Phase::Cancelling { signal, deadline: Some(deadline), .. } if now >= *deadline => {
                let signal = *signal;
                self.set_phase(&mut inner, Phase::Forced { signal });
                TickAction::Terminate(signal)
            }
            Phase::Cancelling { since, second_armed, .. }
                if !*second_armed
                    && now.saturating_duration_since(*since) >= SECOND_SIGNAL_DEBOUNCE =>
            {
                *second_armed = true;
                TickAction::ArmSecondSignal
            }
            Phase::Idle | Phase::Cancelling { .. } => TickAction::Wait,
        }
    }
}

/// The process's coordinator. Its grace is set by [`install`].
static COORDINATOR: LazyLock<Coordinator> = LazyLock::new(|| Coordinator::new(DEFAULT_GRACE));

// ---------------------------------------------------------------------------
// Scopes and polling: usable from any code, inert unless `install` ran.
// ---------------------------------------------------------------------------

/// Declares that the current work holds state needing cleanup and will return
/// promptly once cancelled. While any scope is open, a termination signal
/// cancels instead of killing the process outright.
///
/// Open it before creating what needs cleanup and declare it so that it drops
/// after that state (e.g. as the last field of the owning struct).
#[must_use = "a scope protects the work only while it is alive"]
pub(crate) struct Scope {
    id: u64,
}

/// Open a [`Scope`].
pub(crate) fn scope() -> Scope {
    Scope { id: COORDINATOR.enter() }
}

impl Scope {
    /// Attach (or replace) what a termination request invokes for this scope.
    /// Runs `cancel` at once if a request is already pending.
    ///
    /// `cancel` must be idempotent and must not block: it is called from a
    /// helper thread while the rest of the process keeps running.
    pub(crate) fn on_cancel<F: Fn() + Send + Sync + 'static>(&self, cancel: F) {
        let target: CancelTarget = Arc::new(cancel);
        if COORDINATOR.set_target(self.id, Arc::clone(&target)) {
            target();
        }
    }
}

impl Drop for Scope {
    fn drop(&mut self) {
        COORDINATOR.exit(self.id);
    }
}

impl fmt::Debug for Scope {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Scope").field("id", &self.id).finish()
    }
}

/// The accepted termination signal, if a request is pending.
#[must_use]
pub(crate) fn requested() -> Option<i32> {
    COORDINATOR.requested()
}

/// Poll point for loops outside the pipeline framework.
///
/// # Errors
///
/// [`Cancelled`] once a termination request is pending.
pub(crate) fn check() -> Result<(), Cancelled> {
    #[cfg(test)]
    if TEST_REQUEST.get() {
        return Err(Cancelled);
    }
    match requested() {
        Some(_) => Err(Cancelled),
        None => Ok(()),
    }
}

#[cfg(test)]
thread_local! {
    /// A termination request visible to [`check`] on this thread only, so a
    /// unit test can drive a polling loop's cancellation without touching the
    /// process-wide coordinator that concurrently running tests share.
    static TEST_REQUEST: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Make [`check`] report a pending request on the calling thread (test only).
#[cfg(test)]
pub(crate) fn set_test_request_on_this_thread(pending: bool) {
    TEST_REQUEST.set(pending);
}

/// Hand the command's result back to the termination handling.
///
/// Returns `result` unchanged when no request is pending. Otherwise every
/// destructor inside the command has already run, so this dies by the accepted
/// signal and does not return: a real failure is printed first (it is not what
/// the user asked for); a cancellation is reported only by the one-time
/// `Received <signal>` notice. As PID 1 (where a self-sent signal cannot kill)
/// it exits with status `128 + signal` instead.
///
/// # Errors
///
/// Returns `result`'s error when no termination request is pending.
pub fn conclude(result: anyhow::Result<()>) -> anyhow::Result<()> {
    let Some(signal) = requested() else {
        return result;
    };
    match &result {
        // The notice may not be out yet if the work noticed the request before
        // the helper thread logged it.
        Err(e) if is_cancellation(e) => log_notice_once(signal),
        Err(e) => {
            // Not `eprintln!`, which panics if stderr is gone (a Ctrl-C also
            // kills a `| tee` reading it); the signal death must still happen.
            use std::io::Write;
            let _ = writeln!(std::io::stderr(), "Error: {e:?}");
        }
        Ok(()) => log::debug!("the command completed before the termination request took effect"),
    }
    terminate(signal)
}

// ---------------------------------------------------------------------------
// Process glue: handler installation and the signal thread (Unix only).
// ---------------------------------------------------------------------------

/// Install the termination-signal handling for this process. Call once, from
/// `main`, after the logger is initialised.
///
/// Does nothing when `FGUMI_SIGNAL_GRACE_SECS` is `0`, for signals fgumi was
/// started with ignored, and on non-Unix platforms.
///
/// # Errors
///
/// A handler could not be registered or the signal thread could not start.
/// Any signal already registered then terminates by its default action (as
/// PID 1, it exits with status `128 + signal`).
#[cfg(unix)]
pub fn install() -> std::io::Result<()> {
    let grace = match std::env::var(GRACE_ENV) {
        Ok(value) => parse_grace(&value).unwrap_or_else(|| {
            log::warn!(
                "{GRACE_ENV}={value:?} is not a non-negative number of seconds; using {}s",
                DEFAULT_GRACE.as_secs()
            );
            DEFAULT_GRACE
        }),
        Err(_) => DEFAULT_GRACE,
    };
    if grace.is_zero() {
        return Ok(());
    }
    let ignored = inherited_ignored();
    let handled: Vec<i32> = HANDLED.iter().copied().filter(|sig| !ignored.contains(sig)).collect();
    if handled.is_empty() {
        return Ok(());
    }
    COORDINATOR.set_grace(grace);

    // Registered before the iterator so it runs first: once armed, a signal
    // terminates from the handler itself, with no fgumi thread involved. It
    // starts armed, so a signal arriving before the signal thread is running
    // (or after a failed installation) still terminates by its default action
    // rather than reaching a handler with nothing behind it.
    let force_next = Arc::new(AtomicBool::new(true));
    // As PID 1 a self-raised default signal is ignored, so exit with the
    // shell's code for the signal instead, as `terminate` does.
    let is_init = rustix::process::getpid().is_init();
    for &signal in &handled {
        if is_init {
            signal_hook::flag::register_conditional_shutdown(
                signal,
                128 + signal,
                Arc::clone(&force_next),
            )?;
        } else {
            signal_hook::flag::register_conditional_default(signal, Arc::clone(&force_next))?;
        }
    }
    let signals = signal_hook::iterator::Signals::new(&handled)?;
    let thread_force_next = Arc::clone(&force_next);
    std::thread::Builder::new()
        .name("fgumi-signal".to_string())
        .spawn(move || signal_loop(signals, &thread_force_next))?;
    force_next.store(false, Ordering::SeqCst);
    Ok(())
}

/// Non-Unix platforms keep the default signal behaviour.
///
/// # Errors
///
/// Never.
#[cfg(not(unix))]
pub fn install() -> std::io::Result<()> {
    Ok(())
}

/// Arms default termination when dropped: if the signal thread ever exits
/// (including by a panic), signals terminate by their default action instead
/// of reaching handlers with nothing behind them.
#[cfg(unix)]
struct ArmOnExit<'a>(&'a AtomicBool);

#[cfg(unix)]
impl Drop for ArmOnExit<'_> {
    fn drop(&mut self) {
        self.0.store(true, Ordering::SeqCst);
    }
}

/// The `fgumi-signal` thread. See the module docs for what it may do once a
/// request is pending.
#[cfg(unix)]
#[allow(clippy::needless_pass_by_value)] // the thread owns its `Signals`
fn signal_loop(mut signals: signal_hook::iterator::Signals, force_next: &AtomicBool) {
    let _arm_on_exit = ArmOnExit(force_next);
    // Phase 1: parked in the kernel until the first signal, zero cost for a
    // normal run.
    while COORDINATOR.requested().is_none() {
        for signal in signals.wait() {
            act(COORDINATOR.on_signal(signal, Instant::now()), signal);
        }
    }
    // Phase 2: a request is pending; run the watchdog.
    loop {
        for signal in signals.pending() {
            act(COORDINATOR.on_signal(signal, Instant::now()), signal);
        }
        match COORDINATOR.on_tick(Instant::now()) {
            TickAction::Wait => {}
            TickAction::ArmSecondSignal => {
                force_next.store(true, Ordering::SeqCst);
            }
            TickAction::Terminate(signal) => terminate(signal),
        }
        std::thread::sleep(TICK);
    }
}

#[cfg(unix)]
fn act(action: SignalAction, signal: i32) {
    match action {
        SignalAction::Ignore => {}
        SignalAction::Terminate(signal) => terminate(signal),
        SignalAction::Cancel(targets) => {
            // Targets and the log line run off the signal thread so neither a
            // slow target nor a blocked stderr can stall the watchdog.
            let spawned =
                std::thread::Builder::new().name("fgumi-signal-cancel".to_string()).spawn({
                    let targets = targets.clone();
                    move || {
                        for target in &targets {
                            target();
                        }
                        log_notice_once(signal);
                    }
                });
            // With no thread to run the targets on, terminate as an
            // unhandled signal would: running them here could stall the
            // signal thread, and with it both escape hatches.
            if spawned.is_err() {
                terminate(signal);
            }
        }
    }
}

/// Whether the first-signal notice has been logged.
static NOTICE_LOGGED: AtomicBool = AtomicBool::new(false);

/// Log the first-signal notice, once, from whichever thread gets there first.
#[cfg_attr(not(unix), allow(dead_code))]
fn log_notice_once(signal: i32) {
    if NOTICE_LOGGED.swap(true, Ordering::SeqCst) {
        return;
    }
    #[cfg(unix)]
    let name = signal_hook::low_level::signal_name(signal).unwrap_or("signal");
    #[cfg(not(unix))]
    let name = signal.to_string();
    // A separate statement, so the coordinator's guard is released before the
    // log call, which can block on stderr.
    let grace = COORDINATOR.inner.lock().grace.as_secs_f64();
    log::warn!(
        "Received {name}: stopping and cleaning up (up to {grace:.1}s). Send it again to stop \
         immediately."
    );
}

/// Die by `signal`'s default action.
#[cfg(unix)]
fn terminate(signal: i32) -> ! {
    // A pid-namespace init (fgumi as PID 1 in a container without an init)
    // cannot be killed by a signal it sends itself, so `raise` would return
    // and `abort` would be ignored too. Exit with the shell's code for that
    // signal instead.
    if rustix::process::getpid().is_init() {
        std::process::exit(128 + signal);
    }
    // Restores SIG_DFL, unblocks and raises; for a terminating signal it does
    // not return (it falls back to `abort` itself). `abort`, not
    // `process::exit`, keeps `atexit` handlers (which can block) out of it.
    let _ = signal_hook::low_level::emulate_default_handler(signal);
    std::process::abort()
}

#[cfg(not(unix))]
fn terminate(_signal: i32) -> ! {
    // Unreachable: no request is ever accepted without `install`'s handlers.
    std::process::abort()
}

/// Parse a grace period in seconds: a finite, non-negative decimal.
fn parse_grace(value: &str) -> Option<Duration> {
    let seconds: f64 = value.trim().parse().ok()?;
    Duration::try_from_secs_f64(seconds).ok()
}

/// The handled signals whose disposition was inherited as ignored
/// (`nohup`, a background job in a non-interactive shell). Installing a
/// handler would override that, so those signals are left alone.
///
/// No safe API reads a signal disposition, so this probes: on Linux the
/// `SigIgn:` mask of `/proc/self/status`; elsewhere (or if `/proc` is
/// unreadable) a `/bin/sh` child that signals itself, which survives exactly
/// when the signal is ignored (an ignored disposition survives `exec`, and a
/// shell cannot reset a signal ignored on entry). A failed probe treats the
/// signal as not ignored. Must run before any handler is installed: `exec`
/// resets handled signals to their default.
#[cfg(unix)]
fn inherited_ignored() -> Vec<i32> {
    #[cfg(target_os = "linux")]
    if let Ok(status) = std::fs::read_to_string("/proc/self/status")
        && let Some(mask) = parse_sig_ign(&status)
    {
        return HANDLED.iter().copied().filter(|&sig| mask_contains(mask, sig)).collect();
    }
    probe_ignored_with_shell()
}

#[cfg(unix)]
fn probe_ignored_with_shell() -> Vec<i32> {
    use std::os::unix::process::ExitStatusExt;
    use std::process::{Command, Stdio};

    let children: Vec<(i32, std::io::Result<std::process::Child>)> = HANDLED
        .iter()
        .map(|&signal| {
            let name = signal_hook::low_level::signal_name(signal)
                .and_then(|name| name.strip_prefix("SIG"))
                .unwrap_or("0");
            let child = Command::new("/bin/sh")
                .arg("-c")
                .arg(format!("kill -{name} $$"))
                .stdin(Stdio::null())
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn();
            (signal, child)
        })
        .collect();
    children
        .into_iter()
        .filter_map(|(signal, child)| match child.and_then(|mut child| child.wait()) {
            Ok(status) if status.success() => Some(signal),
            Ok(status) if status.signal() == Some(signal) => None,
            Ok(status) => {
                log::debug!("signal probe for {signal} was inconclusive ({status})");
                None
            }
            Err(e) => {
                log::debug!("could not probe whether signal {signal} is ignored: {e}");
                None
            }
        })
        .collect()
}

/// The `SigIgn:` mask from a `/proc/<pid>/status` file.
#[cfg_attr(not(target_os = "linux"), allow(dead_code))]
fn parse_sig_ign(status: &str) -> Option<u64> {
    status
        .lines()
        .find_map(|line| line.strip_prefix("SigIgn:"))
        .and_then(|hex| u64::from_str_radix(hex.trim(), 16).ok())
}

/// Whether `signal` is set in a kernel signal mask (bit `signal - 1`).
#[cfg_attr(not(target_os = "linux"), allow(dead_code))]
fn mask_contains(mask: u64, signal: i32) -> bool {
    u32::try_from(signal - 1).ok().and_then(|bit| mask.checked_shr(bit)).is_some_and(|m| m & 1 == 1)
}

#[cfg(test)]
mod tests {
    use std::sync::atomic::{AtomicUsize, Ordering};

    use rstest::rstest;

    use super::*;

    const INT: i32 = 2;
    const TERM: i32 = 15;
    const GRACE: Duration = Duration::from_secs(8);

    fn ms(n: u64) -> Duration {
        Duration::from_millis(n)
    }

    /// A coordinator with `scopes` open scopes, each with a counting target.
    fn coordinator(scopes: usize) -> (Coordinator, Arc<AtomicUsize>) {
        let coordinator = Coordinator::new(GRACE);
        let calls = Arc::new(AtomicUsize::new(0));
        for _ in 0..scopes {
            let id = coordinator.enter();
            let calls = Arc::clone(&calls);
            assert!(!coordinator.set_target(
                id,
                Arc::new(move || {
                    calls.fetch_add(1, Ordering::SeqCst);
                })
            ));
        }
        (coordinator, calls)
    }

    /// A signal action reduced to something comparable.
    #[derive(Debug, PartialEq, Eq)]
    enum Act {
        Terminate(i32),
        Cancel(usize),
        Ignore,
    }

    fn act_of(action: SignalAction) -> Act {
        match action {
            SignalAction::Terminate(signal) => Act::Terminate(signal),
            SignalAction::Cancel(targets) => Act::Cancel(targets.len()),
            SignalAction::Ignore => Act::Ignore,
        }
    }

    #[rstest]
    #[case::idle_no_scope_terminates(0, &[], INT, Act::Terminate(INT), Some(INT))]
    #[case::idle_with_scope_cancels(2, &[], INT, Act::Cancel(2), Some(INT))]
    #[case::duplicate_inside_debounce_ignored(1, &[(INT, 0)], TERM, Act::Ignore, Some(INT))]
    #[case::second_after_debounce_terminates(1, &[(INT, 0)], TERM, Act::Terminate(TERM), Some(TERM))]
    fn on_signal_transitions(
        #[case] scopes: usize,
        #[case] earlier: &[(i32, u64)],
        #[case] signal: i32,
        #[case] expected: Act,
        #[case] requested_after: Option<i32>,
    ) {
        let (coordinator, _) = coordinator(scopes);
        let t0 = Instant::now();
        for &(earlier_signal, at) in earlier {
            let _ = coordinator.on_signal(earlier_signal, t0 + ms(at));
        }
        let at = if expected == Act::Ignore { ms(100) } else { SECOND_SIGNAL_DEBOUNCE };
        assert_eq!(act_of(coordinator.on_signal(signal, t0 + at)), expected);
        assert_eq!(coordinator.requested(), requested_after);
    }

    #[test]
    fn signal_after_forced_terminates_again() {
        let (coordinator, _) = coordinator(0);
        let now = Instant::now();
        assert_eq!(act_of(coordinator.on_signal(INT, now)), Act::Terminate(INT));
        assert_eq!(act_of(coordinator.on_signal(TERM, now)), Act::Terminate(TERM));
    }

    #[test]
    fn zero_grace_terminates_even_with_a_scope() {
        let (coordinator, _) = coordinator(1);
        coordinator.set_grace(Duration::ZERO);
        assert_eq!(act_of(coordinator.on_signal(INT, Instant::now())), Act::Terminate(INT));
    }

    #[rstest]
    #[case::before_debounce(ms(100), TickAction::Wait)]
    #[case::at_debounce(SECOND_SIGNAL_DEBOUNCE, TickAction::ArmSecondSignal)]
    #[case::at_deadline(GRACE, TickAction::Terminate(INT))]
    #[case::past_deadline(GRACE + ms(1), TickAction::Terminate(INT))]
    fn on_tick_after_a_signal(#[case] elapsed: Duration, #[case] expected: TickAction) {
        let (coordinator, _) = coordinator(1);
        let t0 = Instant::now();
        let _ = coordinator.on_signal(INT, t0);
        assert_eq!(coordinator.on_tick(t0 + elapsed), expected);
    }

    #[test]
    fn on_tick_arms_once_then_waits_then_terminates() {
        let (coordinator, _) = coordinator(1);
        let t0 = Instant::now();
        let _ = coordinator.on_signal(TERM, t0);
        assert_eq!(coordinator.on_tick(t0 + ms(600)), TickAction::ArmSecondSignal);
        assert_eq!(coordinator.on_tick(t0 + ms(700)), TickAction::Wait);
        assert_eq!(coordinator.on_tick(t0 + GRACE), TickAction::Terminate(TERM));
        // Forced is terminal.
        assert_eq!(coordinator.on_tick(t0 + GRACE + ms(1)), TickAction::Terminate(TERM));
    }

    /// A grace too large to add to an `Instant` means no watchdog, not a
    /// panic on the signal thread (which would leave signals swallowed).
    #[test]
    fn unrepresentable_grace_disables_the_watchdog() {
        let (coordinator, _) = coordinator(1);
        coordinator.set_grace(Duration::MAX);
        let t0 = Instant::now();
        assert_eq!(act_of(coordinator.on_signal(INT, t0)), Act::Cancel(1));
        assert_eq!(coordinator.on_tick(t0 + SECOND_SIGNAL_DEBOUNCE), TickAction::ArmSecondSignal);
        assert_eq!(coordinator.on_tick(t0 + Duration::from_secs(1_000_000)), TickAction::Wait);
    }

    #[test]
    fn on_tick_idle_waits() {
        let (coordinator, _) = coordinator(1);
        assert_eq!(coordinator.on_tick(Instant::now()), TickAction::Wait);
    }

    #[test]
    fn target_attached_while_idle_is_in_the_snapshot() {
        let (coordinator, calls) = coordinator(1);
        let SignalAction::Cancel(targets) = coordinator.on_signal(INT, Instant::now()) else {
            panic!("expected Cancel");
        };
        for target in &targets {
            target();
        }
        assert_eq!(calls.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn target_attached_while_cancelling_is_reported_not_snapshotted() {
        let coordinator = Coordinator::new(GRACE);
        let id = coordinator.enter();
        // A scope with no target yet: the request is accepted and sticky.
        let SignalAction::Cancel(targets) = coordinator.on_signal(INT, Instant::now()) else {
            panic!("expected Cancel");
        };
        assert!(targets.is_empty());
        assert!(coordinator.set_target(id, Arc::new(|| {})), "caller must invoke it");
    }

    #[test]
    fn exiting_the_last_scope_while_cancelling_keeps_the_request() {
        let coordinator = Coordinator::new(GRACE);
        let id = coordinator.enter();
        let _ = coordinator.on_signal(INT, Instant::now());
        coordinator.exit(id);
        assert_eq!(coordinator.requested(), Some(INT));
    }

    #[test]
    fn signal_after_every_scope_exited_terminates() {
        let coordinator = Coordinator::new(GRACE);
        let id = coordinator.enter();
        coordinator.exit(id);
        assert_eq!(act_of(coordinator.on_signal(TERM, Instant::now())), Act::Terminate(TERM));
    }

    #[test]
    fn global_scope_is_inert_without_install() {
        let scope = scope();
        let ran = Arc::new(AtomicUsize::new(0));
        let counter = Arc::clone(&ran);
        scope.on_cancel(move || {
            counter.fetch_add(1, Ordering::SeqCst);
        });
        assert_eq!(ran.load(Ordering::SeqCst), 0);
        assert_eq!(requested(), None);
        assert!(check().is_ok());
        assert!(conclude(Err(anyhow::anyhow!("boom"))).is_err());
    }

    #[rstest]
    #[case::integer("8", Some(Duration::from_secs(8)))]
    #[case::zero("0", Some(Duration::ZERO))]
    #[case::fraction("1.5", Some(ms(1500)))]
    #[case::padded(" 2 ", Some(Duration::from_secs(2)))]
    #[case::empty("", None)]
    #[case::negative("-1", None)]
    #[case::text("abc", None)]
    #[case::infinite("inf", None)]
    #[case::nan("NaN", None)]
    fn parse_grace_cases(#[case] value: &str, #[case] expected: Option<Duration>) {
        assert_eq!(parse_grace(value), expected);
    }

    #[rstest]
    #[case::none("Name:\tfgumi\nSigIgn:\t0000000000000000\n", Some(0))]
    #[case::hup("SigIgn:\t0000000000000001\n", Some(0x1))]
    #[case::int("SigIgn:\t0000000000000002\n", Some(0x2))]
    #[case::hup_term("SigBlk:\t0\nSigIgn:\t0000000000004001\nSigCgt:\t0\n", Some(0x4001))]
    #[case::missing("Name:\tfgumi\n", None)]
    fn parse_sig_ign_cases(#[case] status: &str, #[case] expected: Option<u64>) {
        assert_eq!(parse_sig_ign(status), expected);
    }

    #[rstest]
    #[case::hup_set(0x1, 1, true)]
    #[case::int_set(0x2, 2, true)]
    #[case::term_set(0x4000, 15, true)]
    #[case::term_clear(0x4001 & !0x4000, 15, false)]
    #[case::zero_signal(u64::MAX, 0, false)]
    #[case::out_of_range(u64::MAX, 65, false)]
    fn mask_contains_cases(#[case] mask: u64, #[case] signal: i32, #[case] expected: bool) {
        assert_eq!(mask_contains(mask, signal), expected);
    }

    #[rstest]
    #[case::bare(anyhow::Error::new(Cancelled), true)]
    #[case::under_context(anyhow::Error::new(Cancelled).context("a").context("b"), true)]
    #[case::unrelated(anyhow::anyhow!("boom"), false)]
    fn is_cancellation_cases(#[case] err: anyhow::Error, #[case] expected: bool) {
        assert_eq!(is_cancellation(&err), expected);
    }

    #[test]
    fn documented_bounds() {
        assert_eq!(DEFAULT_GRACE, Duration::from_secs(8));
        assert_eq!(SECOND_SIGNAL_DEBOUNCE, ms(500));
    }
}
