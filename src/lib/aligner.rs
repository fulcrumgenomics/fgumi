//! Aligner subprocess management for the runall `AlignAndMerge` stage.
//!
//! Two halves:
//!
//! * [`AlignerProcess`] spawns an aligner as a child process via
//!   `/bin/bash -c <command>` and exposes managed stdin / stdout / stderr
//!   handles. A background thread continuously drains the child's stderr
//!   to fgumi's own stderr (so aligner progress / warnings are visible
//!   live) AND retains the last `ring_size` lines in a ring buffer that
//!   gets surfaced in the error message if the subprocess exits non-zero.
//!
//! * [`AlignerPreset`] builds well-formed aligner shell commands for the
//!   two presets fgumi knows about (`bwa-mem3`, `bwa`) and validates
//!   that the reference's index files exist before any input is read.
//!   For `--aligner::command "..."` users (free-form mode), the
//!   [`substitute_template`] helper fills the `{ref}` / `{threads}`
//!   placeholders.
//!
//! Used by the subprocess align backend's `SubprocessAlignStep` in
//! `src/lib/pipeline/steps/align/subprocess.rs`. This module is the
//! framework-agnostic subprocess primitive; the typed `Step` impl that owns the
//! I/O threads lives in that step module.

use std::collections::VecDeque;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::Duration;

use anyhow::{Context, Result, bail};
use clap::Args;
use fgumi_cli_macros::multi_options;
use fgumi_consensus::MethylationMode;

use crate::commands::common::{EM_SEQ_MODE_NAME, TAPS_MODE_NAME};
use crate::commands::fastq::UnalignedReadFormat;

mod scoring;
mod star;
pub(crate) use scoring::{BwaScoring, command_passes_meth};
use scoring::{mem_option_words, mem_word_end};

/// A scratch directory an aligner subprocess writes into (`star`'s logs and
/// per-chunk files), kept until the aligner has exited and then removed when
/// the last holder drops it.
pub(crate) type AlignerScratchDir = std::sync::Arc<tempfile::TempDir>;

// ============================================================================
// AlignerProcess
// ============================================================================

/// How a child that exited without a code died: `signal N` on Unix, where
/// that is the only other way to exit, else a bare `signal`.
fn killed_by(status: std::process::ExitStatus) -> String {
    #[cfg(unix)]
    if let Some(signal) = std::os::unix::process::ExitStatusExt::signal(&status) {
        return format!("signal {signal}");
    }
    let _ = status;
    "signal".to_string()
}

/// A running aligner subprocess with managed stdin, stdout, and stderr.
///
/// The aligner is spawned via `/bin/bash -c`, allowing the command string
/// to contain pipes, redirects, and other shell constructs. Stderr from
/// the child process is relayed to fgumi's own stderr in real time, and
/// the last `ring_size` lines are retained for inclusion in error
/// messages if the process exits non-zero.
///
/// # Drain contract (avoid deadlock)
///
/// The child is spawned with all three standard pipes piped. The stderr
/// pipe is drained automatically by a background thread. Stdin and
/// stdout, however, are NOT drained by the framework — the caller MUST
/// arrange to:
///
///   * write FASTQ to `take_stdin()`'s handle, then drop it to signal EOF,
///   * concurrently read SAM/BAM from `take_stdout()`'s handle,
///
/// before calling [`Self::wait`]. If the caller writes stdin without
/// concurrently draining stdout, the kernel's ~64 KB stdout pipe buffer
/// fills, the aligner blocks on `write`, the caller blocks on `write` to
/// stdin, and the pipeline deadlocks with no error. The intended caller
/// is Step 2 of the runall `AlignAndMerge` chain, which spawns paired
/// reader/writer threads precisely for this reason.
///
/// [`Self::wait`] additionally drops any not-yet-taken stdin handle
/// before waiting so a misuse (caller forgets to take stdin) surfaces
/// as a fast subprocess exit rather than a hang.
pub struct AlignerProcess {
    child: Child,
    /// Wrapper status observed before teardown can have caused its exit.
    exit_before_teardown: Option<std::process::ExitStatus>,
    /// Once signals start, their exit status must not mask a reader failure.
    teardown_started: bool,
    /// The child's stdin pipe, available until [`Self::take_stdin`] is called.
    stdin: Option<ChildStdin>,
    /// The child's stdout pipe, available until [`Self::take_stdout`] is called.
    stdout: Option<ChildStdout>,
    /// Background thread that relays child stderr into [`Self::stderr_capture`].
    stderr_thread: Option<JoinHandle<()>>,
    /// What the relay thread has captured of the aligner's stderr. Shared rather
    /// than returned by the thread, so it can be read even when the relay is
    /// still blocked on a stderr pipe that a process the aligner started holds
    /// open (see [`STDERR_DRAIN_GRACE`]).
    stderr_capture: Arc<parking_lot::Mutex<StderrCapture>>,
    /// The aligner's process group (its id is the wrapper shell's PID), which
    /// [`Self::kill`] signals so the aligner's own processes die with the
    /// wrapper. `None` once the group has been torn down or the wrapper exited
    /// through [`Self::wait`]: after that the id may name unrelated processes,
    /// so it is never signalled again.
    #[cfg(unix)]
    process_group: Option<rustix::process::Pid>,
    /// What the group teardown in [`Self::kill`] concluded: `Some(false)` when
    /// the wrapper or another group member could not be confirmed gone, so
    /// something may still hold the stderr pipe. `Drop` reads it to skip the
    /// stderr-relay join instead of re-deriving exit from the (reaped)
    /// wrapper alone.
    #[cfg(unix)]
    group_teardown_confirmed: Option<bool>,
}

/// How long [`AlignerProcess::kill`] waits after SIGTERM to the aligner's
/// process group before escalating to SIGKILL.
#[cfg(unix)]
const PROCESS_GROUP_TERM_GRACE: Duration = Duration::from_millis(500);

/// How long [`AlignerProcess::kill`] waits for the wrapper to be reaped after
/// SIGKILL.
const KILL_REAP_DEADLINE: Duration = Duration::from_secs(1);

/// How long [`AlignerProcess::wait`] and `Drop` let the stderr relay drain after
/// the aligner has exited. The relay ends at EOF on the stderr pipe, which a
/// process the aligner started (and left running) can hold open indefinitely;
/// past this grace the relay is left to finish on its own instead of hanging
/// teardown.
const STDERR_DRAIN_GRACE: Duration = Duration::from_secs(2);

/// Interval at which the bounded stderr-relay join polls for the relay's end.
/// Short, since the relay usually ends within moments of the aligner's exit.
const STDERR_JOIN_POLL: Duration = Duration::from_millis(5);

/// Interval at which [`AlignerProcess::kill`] polls for exit.
const KILL_POLL_INTERVAL: Duration = Duration::from_millis(50);

/// The mimalloc option environment variables (`MIMALLOC_PURGE_DELAY` and its
/// legacy name) a user can set to choose the purge delay themselves.
const MIMALLOC_PURGE_ENV: [&str; 2] = ["MIMALLOC_PURGE_DELAY", "MIMALLOC_RESET_DELAY"];

/// Whether the user set mimalloc's purge delay in the environment.
pub(crate) fn user_set_mimalloc_purge() -> bool {
    names_set_mimalloc_purge(std::env::vars_os().map(|(name, _)| name))
}

/// Delay mimalloc returning freed pages to the OS (by
/// [`fgumi_sort::RETAINED_PURGE_DELAY_MS`]), unless the user set the purge
/// delay. The align stage frees and reallocates large per-batch buffers on
/// every pool thread, and mimalloc's default purge (decommit after 1 s) turns
/// that reuse into millions of page faults.
///
/// The setting is process-wide and is never restored, so it covers every stage
/// of the `runall` process, not just align. The delay is bounded rather than
/// never-purge on purpose: with `-1`, memory a fused sort frees at each spill
/// is never returned, and peak RSS grows by roughly the sort's memory budget
/// (+5.5 GB in-process on 30M pairs, enough to exceed a 30 GB limit). Measured
/// on 3M pairs at 32 threads (extract through simplex consensus, subprocess and
/// in-process aligners, with and without a spilling sort), the 60 s delay
/// matches never-purging on wall time, CPU and page faults (~30k vs ~750k at
/// mimalloc's default, which is ~3% slower); on 30M pairs it gives back most of
/// never-purging's extra memory (fgumi 13.4 GB vs 17.6 GB, 11.6 GB at the
/// default). Set `MIMALLOC_PURGE_DELAY` to choose otherwise.
///
/// Both align backends call it at wiring; `stage` labels the debug line that
/// reports the resulting delay. A sort does not: v0.7.0 never set the delay,
/// and on a production-sized sort the 60 s delay cost 3–5 GiB of peak RSS with
/// no speed benefit.
pub(crate) fn retain_freed_memory_unless_user_set(stage: &str) {
    retain_unless_user_set(
        std::env::vars_os().map(|(name, _)| name),
        fgumi_sort::retain_freed_memory,
    );
    log::debug!("{stage}: mimalloc purge delay {} ms", fgumi_sort::mi_purge_delay_ms());
}

/// Call `retain` unless `names` (environment variable names) set the purge
/// delay; returns whether it was called.
fn retain_unless_user_set<I, S>(names: I, retain: impl FnOnce()) -> bool
where
    I: IntoIterator<Item = S>,
    S: AsRef<std::ffi::OsStr>,
{
    let retained = !names_set_mimalloc_purge(names);
    if retained {
        retain();
    }
    retained
}

/// Whether `names` (environment variable names) include either
/// [`MIMALLOC_PURGE_ENV`] name. The comparison ignores ASCII case because
/// mimalloc's own lookup does, so `mimalloc_purge_delay=250` is the user's
/// choice too. Any value counts, including one mimalloc cannot parse: that one
/// is mimalloc's to reject, not ours to replace.
fn names_set_mimalloc_purge<I, S>(names: I) -> bool
where
    I: IntoIterator<Item = S>,
    S: AsRef<std::ffi::OsStr>,
{
    names.into_iter().any(|name| {
        let name = name.as_ref();
        MIMALLOC_PURGE_ENV.iter().any(|option| name.eq_ignore_ascii_case(option))
    })
}

impl AlignerProcess {
    /// Spawn an aligner subprocess with the given shell command.
    ///
    /// The command is run via `/bin/bash -c`, so it supports pipes and
    /// redirects. A background thread is started immediately to relay
    /// the child's stderr to fgumi's stderr in real time; the last
    /// `ring_size` lines are captured for failure diagnostics.
    ///
    /// # Arguments
    ///
    /// * `command`   - Shell command to run (passed to `/bin/bash -c`).
    /// * `ring_size` - Number of stderr lines to retain for error reporting.
    ///
    /// # Errors
    ///
    /// Returns an error if the subprocess cannot be spawned.
    ///
    /// # Panics
    ///
    /// Panics if the OS does not provide the stderr pipe after
    /// configuring `Stdio::piped()`, which should never happen in
    /// practice.
    pub fn spawn(command: &str, ring_size: usize) -> Result<Self> {
        let mut cmd = Command::new("/bin/bash");
        cmd.args(["-c", command])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        // bwa-mem3 links mimalloc, whose default purge decommits freed pages
        // after 1 s and costs the aligner page faults on every batch, so the
        // child gets the same bounded delay as fgumi itself
        // (`fgumi_sort::RETAINED_PURGE_DELAY_MS`). A user-set value (either
        // name, any case) is inherited unchanged; aligners without mimalloc
        // ignore the variable. The variable reaches every process in the shell
        // command, so with `--aligner::command` any mimalloc-linked tool in the
        // user's pipeline also delays returning freed pages; set
        // `MIMALLOC_PURGE_DELAY` to choose otherwise.
        if !user_set_mimalloc_purge() {
            cmd.env("MIMALLOC_PURGE_DELAY", fgumi_sort::RETAINED_PURGE_DELAY_MS.to_string());
        }
        // Run the aligner in its own process group, so teardown can signal
        // every process the shell started (`kill()`). Killing only the
        // `/bin/bash -c` wrapper leaves the real aligner, a grandchild, holding
        // the pipes open, and fgumi never exits after a failure.
        //
        // The cost is that no signal sent to fgumi's process group reaches the
        // aligner any more: not the terminal's SIGINT (Ctrl-C) or SIGTSTP
        // (Ctrl-Z), and not a group-directed kill from GNU `timeout`, a shell's
        // `kill %job`, or a supervisor that signals the whole group. fgumi
        // installs no handler for these, so such a signal kills fgumi by the
        // default action without running any teardown; the aligner then sees
        // EOF on stdin, and EPIPE/SIGPIPE when it next writes to stdout, which
        // can be a whole batch later. fgumi's own failure and cancellation
        // paths kill the group through `kill()`.
        //
        // As a background group in fgumi's session, the aligner is also stopped
        // (SIGTTIN/SIGTTOU) if it reads or configures the controlling terminal
        // itself, by opening `/dev/tty` as an ssh passphrase or sudo password
        // prompt does; pipes cover only its standard streams. An
        // `--aligner::command` must therefore not prompt on the terminal (use
        // e.g. `ssh -o BatchMode=yes`, `sudo -n`).
        #[cfg(unix)]
        std::os::unix::process::CommandExt::process_group(&mut cmd, 0);
        let mut child =
            cmd.spawn().with_context(|| format!("failed to spawn aligner command: {command}"))?;
        // `process_group(0)` makes the group id the child's PID.
        #[cfg(unix)]
        let process_group = Some(rustix::process::Pid::from_child(&child));

        let stdin = child.stdin.take();
        let stdout = child.stdout.take();
        let child_stderr = child.stderr.take().expect("stderr was configured as piped");

        let stderr_capture = Arc::new(parking_lot::Mutex::new(StderrCapture::new(ring_size)));
        let relay_capture = Arc::clone(&stderr_capture);
        let stderr_thread = thread::spawn(move || relay_stderr(child_stderr, &relay_capture));

        Ok(Self {
            child,
            exit_before_teardown: None,
            teardown_started: false,
            stdin,
            stdout,
            stderr_thread: Some(stderr_thread),
            stderr_capture,
            #[cfg(unix)]
            process_group,
            #[cfg(unix)]
            group_teardown_confirmed: None,
        })
    }

    /// Take the stdin pipe for writing FASTQ data.
    ///
    /// Returns `None` if stdin has already been taken. **The caller
    /// must also drain stdout concurrently (via [`Self::take_stdout`]
    /// on a separate thread) — otherwise the kernel pipe buffer fills
    /// and the pipeline deadlocks.** See the struct-level "Drain
    /// contract" section.
    pub fn take_stdin(&mut self) -> Option<ChildStdin> {
        self.stdin.take()
    }

    /// Take the stdout pipe for reading SAM data.
    ///
    /// Returns `None` if stdout has already been taken. **The caller
    /// must also feed stdin concurrently (via [`Self::take_stdin`] on
    /// a separate thread) — otherwise the aligner produces no output
    /// and the reader blocks indefinitely.** See the struct-level
    /// "Drain contract" section.
    pub fn take_stdout(&mut self) -> Option<ChildStdout> {
        self.stdout.take()
    }

    /// Check without blocking whether both the child and its stderr relay
    /// have finished, so [`Self::wait`] can collect their results immediately.
    ///
    /// # Errors
    ///
    /// Returns an error if querying the child's exit status fails.
    pub(crate) fn is_finished(&mut self) -> std::io::Result<bool> {
        Ok(self.child.try_wait()?.is_some()
            && self.stderr_thread.as_ref().is_none_or(JoinHandle::is_finished))
    }

    /// Wait for the aligner process to finish.
    ///
    /// Drops any not-yet-taken stdin/stdout handle before waiting — most
    /// aligners (bwa-mem3, bwa) read stdin until EOF, so an undropped
    /// stdin would hang the child forever. A caller who legitimately
    /// wants to feed stdin from a thread should call `take_stdin()`
    /// first, write+drop in that thread, then call `wait()`.
    ///
    /// Joins the stderr relay thread and includes the last captured
    /// stderr lines in the error message if the process exits with a
    /// non-zero status. The join waits at most `STDERR_DRAIN_GRACE` (2 s) after
    /// the wrapper exits: a process the wrapper started that still holds the
    /// stderr pipe after the wrapper exits (with any status) is left running,
    /// since `wait()` never signals the group once the wrapper is reaped, and
    /// its stderr is still relayed. The lines reported are then those captured
    /// by the deadline (including a last line not yet ended by a newline), and
    /// the message says later lines are not shown. If the detached relay still
    /// holds the capture lock, the error reports the status without a tail
    /// rather than waiting for that lock. Successful exits never copy the tail.
    ///
    /// # Errors
    ///
    /// Returns an error if the process exits with a non-zero status
    /// code or if waiting on the process fails.
    pub fn wait(mut self) -> Result<()> {
        // Close any retained stdin/stdout BEFORE waiting. If a caller
        // never called `take_stdin()`, the child would block on
        // `read(0, ...)` forever and `child.wait()` would deadlock.
        // Dropping signals EOF to the child. Same logic applies to
        // stdout (less common) — an undrained stdout pipe blocks the
        // child on write. The intended caller drains both via spawned
        // threads before calling `wait()`; this is the safety net for
        // direct / partial use.
        drop(self.stdin.take());
        drop(self.stdout.take());

        let status = match self.child.wait() {
            Ok(status) => status,
            Err(err) => {
                // `ECHILD`: the wrapper was reaped elsewhere, so its PID (the
                // group id) may already be reused; never signal the group
                // again, including from `Drop`. Any other error leaves the
                // wrapper possibly alive, so `Drop` still tears the group down.
                #[cfg(unix)]
                if is_already_reaped(&err) {
                    self.process_group = None;
                }
                return Err(err).context("failed to wait for aligner process");
            }
        };
        // The wrapper is reaped, so its PID (the group id) can be reused once
        // the group empties: never signal the group from here on, including
        // from `Drop`.
        #[cfg(unix)]
        {
            self.process_group = None;
        }

        // Let the stderr relay drain (bounded: a process the aligner started may
        // still hold the pipe open), then read the lines it captured.
        self.exit_result(status)
    }

    /// Kill after a reader failure, retaining diagnostics from a natural exit.
    /// Unlike `wait`, this never waits without a deadline for a live wrapper or
    /// requires the process group to disappear before reporting a cached status.
    pub(crate) fn finish_teardown(&mut self) -> Result<()> {
        self.kill();
        match self.exit_before_teardown {
            Some(status) => self.exit_result(status),
            None => Ok(()),
        }
    }

    /// Report an available wrapper status with a bounded stderr snapshot.
    fn exit_result(&mut self, status: std::process::ExitStatus) -> Result<()> {
        let drained = self.join_stderr_relay();
        if status.success() {
            return Ok(());
        }
        // A detached relay can still be decoding a large line while holding
        // the capture mutex. Never wait for that lock after the drain deadline.
        let snapshot = self.stderr_capture.try_lock().map(|capture| capture.tail());
        let incomplete = !drained || snapshot.is_none();
        let mut last_lines = snapshot.unwrap_or_default();
        if incomplete {
            last_lines.push(
                "(the aligner's stderr capture was incomplete when fgumi stopped waiting; \
                 later lines are not shown)"
                    .to_string(),
            );
        }

        let code = status.code().map_or_else(|| killed_by(status), |c| c.to_string());
        let tail = if last_lines.is_empty() {
            "(no stderr captured)".to_string()
        } else {
            last_lines.join("\n")
        };
        bail!("aligner process exited with status {code}. Last stderr lines:\n{tail}");
    }

    /// Join the stderr relay once it reaches EOF, waiting at most
    /// [`STDERR_DRAIN_GRACE`], and return whether it did. Called only after the
    /// aligner has exited, when the relay ends as soon as every holder of the
    /// stderr pipe has closed it and the relay has printed what it read. A
    /// process the aligner started and left running (e.g. a backgrounded
    /// helper) can hold the pipe open indefinitely, and a slow consumer of
    /// fgumi's own stderr can hold the relay back; then the relay is left
    /// detached (it keeps relaying and ends at EOF), rather than hanging.
    /// Neither timeout nor panic logs on the caller: that sink may be blocked
    /// too. Both mark the returned error's stderr snapshot as incomplete.
    fn join_stderr_relay(&mut self) -> bool {
        let Some(handle) = self.stderr_thread.take() else {
            return true;
        };
        let deadline = std::time::Instant::now() + STDERR_DRAIN_GRACE;
        while !handle.is_finished() {
            if std::time::Instant::now() >= deadline {
                return false;
            }
            thread::sleep(STDERR_JOIN_POLL);
        }
        handle.join().is_ok()
    }

    /// Kill the aligner process.
    ///
    /// On Unix the whole process group is killed, not just the `/bin/bash -c`
    /// wrapper: the real aligner is usually the wrapper's child, and while it
    /// lives it holds the stdin/stdout/stderr pipes open, so the I/O threads
    /// and the stderr relay never finish. The group gets SIGTERM first, so an
    /// aligner (or a wrapper script) that cleans up on SIGTERM, e.g. removes
    /// its scratch files or flushes its log, gets the chance; whatever is still
    /// alive after a 500 ms grace (`PROCESS_GROUP_TERM_GRACE`) gets SIGKILL. Only the first
    /// call signals the group, and a call after [`Self::wait`] reaped the
    /// wrapper, or after the wrapper was reaped elsewhere (`ECHILD`), never does
    /// (the group id may have been reused). Later calls return the recorded
    /// group result; a reaped wrapper alone cannot confirm group exit. On
    /// other platforms only the wrapper is killed (SIGKILL).
    ///
    /// Then waits up to 1 s (`KILL_REAP_DEADLINE`) for the wrapper to be reaped
    /// and, on Unix, for the rest of the group to be gone. Returns `true` if
    /// both are confirmed (or were already the case), and `false` otherwise:
    /// the wrapper was still running at the deadline (e.g. stuck in
    /// uninterruptible I/O), or another group member was, which may still hold
    /// the pipes. A group whose orphaned members are never reaped (fgumi as PID
    /// 1 in a container without an init) also reports `false`, since a zombie
    /// still counts as a member, and always waits out the grace. Callers MUST
    /// honor a `false` return by NOT issuing a subsequent unbounded
    /// [`wait`](Self::wait) on the child, which would re-introduce the very
    /// teardown hang this bounded kill exists to prevent.
    pub fn kill(&mut self) -> bool {
        #[cfg(unix)]
        if let Some(confirmed) = self.group_teardown_confirmed {
            return confirmed;
        }
        #[cfg(unix)]
        if let Some(group) = self.process_group.take() {
            let confirmed = self.kill_process_group(group).unwrap_or(false);
            self.group_teardown_confirmed = Some(confirmed);
            return confirmed;
        }
        self.kill_wrapper(std::time::Instant::now() + KILL_REAP_DEADLINE)
    }

    /// SIGTERM `group`, give it [`PROCESS_GROUP_TERM_GRACE`] to exit, then
    /// SIGKILL whatever is left, reap the wrapper, and confirm the group is
    /// gone. An external wrapper reap stops signalling without confirming the
    /// group is gone, since the group identity can no longer be trusted.
    ///
    /// The group's id is the wrapper's PID. While the wrapper is unreaped that
    /// PID cannot be reused, and once it is reaped (here) the id still cannot
    /// be handed out as a new PID while the group has members, so after the
    /// reap the group is signalled only after a probe found it non-empty. The
    /// signal follows that probe immediately, without an intervening poll sleep.
    /// A wrapper reaped elsewhere (`ECHILD`) gives no such guarantee, so its
    /// group is never signalled.
    #[cfg(unix)]
    fn kill_process_group(&mut self, group: rustix::process::Pid) -> Result<bool, WrapperState> {
        use rustix::process::Signal;

        if self.wait_for_group(group, std::time::Instant::now())? {
            return Ok(true);
        }
        self.teardown_started |= signal_process_group(group, Signal::TERM);
        if self.wait_for_group(group, std::time::Instant::now() + PROCESS_GROUP_TERM_GRACE)? {
            return Ok(true);
        }
        if self.wait_for_group(group, std::time::Instant::now())? {
            return Ok(true);
        }
        let deadline = std::time::Instant::now() + KILL_REAP_DEADLINE;
        self.teardown_started |= signal_process_group(group, Signal::KILL);
        self.reap_killed_group(group, deadline)
    }

    /// Reap a killed wrapper and confirm group exit within the escalation deadline.
    #[cfg(unix)]
    fn reap_killed_group(
        &mut self,
        group: rustix::process::Pid,
        deadline: std::time::Instant,
    ) -> Result<bool, WrapperState> {
        if !self.kill_wrapper(deadline) {
            return Ok(false);
        }
        // The wrapper is reaped; confirm no other member (the real aligner)
        // survived the SIGKILL, within the same reap deadline.
        self.wait_for_group(group, deadline)
    }

    /// Poll until `deadline` for the wrapper to exit (reaping it) and for
    /// `group` to become empty. An external reap returns an error: signalling must
    /// stop, but the group's exit is not confirmed by the wrapper alone.
    #[cfg(unix)]
    fn wait_for_group(
        &mut self,
        group: rustix::process::Pid,
        deadline: std::time::Instant,
    ) -> Result<bool, WrapperState> {
        loop {
            // Probe the group only once the wrapper is reaped: an unreaped
            // (zombie) wrapper is still a member, so the probe would succeed.
            match self.wrapper_state() {
                WrapperState::ReapedElsewhere => return Err(WrapperState::ReapedElsewhere),
                WrapperState::Exited if process_group_is_empty(group) => return Ok(true),
                WrapperState::Exited | WrapperState::Running | WrapperState::Unknown => {}
            }
            if std::time::Instant::now() >= deadline {
                return Ok(false);
            }
            thread::sleep(KILL_POLL_INTERVAL);
        }
    }

    /// The wrapper's state, reaping it if it has exited.
    #[cfg(unix)]
    fn wrapper_state(&mut self) -> WrapperState {
        match self.child.try_wait() {
            Ok(Some(status)) => {
                if !self.teardown_started {
                    self.exit_before_teardown = Some(status);
                }
                WrapperState::Exited
            }
            Ok(None) => WrapperState::Running,
            Err(err) if is_already_reaped(&err) => WrapperState::ReapedElsewhere,
            Err(_) => WrapperState::Unknown,
        }
    }

    /// SIGKILL the wrapper only and wait until `deadline` for it
    /// to be reaped. Returns as [`Self::kill`] does.
    fn kill_wrapper(&mut self, deadline: std::time::Instant) -> bool {
        if !self.teardown_started
            && let Ok(Some(status)) = self.child.try_wait()
        {
            self.exit_before_teardown = Some(status);
            return true;
        }
        self.teardown_started |= self.child.kill().is_ok();
        loop {
            match self.child.try_wait() {
                Ok(Some(status)) => {
                    if !self.teardown_started {
                        self.exit_before_teardown = Some(status);
                    }
                    return true;
                }
                Ok(None) => {
                    if std::time::Instant::now() >= deadline {
                        return false;
                    }
                    thread::sleep(KILL_POLL_INTERVAL);
                }
                // ECHILD confirms a prior reap; other errors leave exit unknown.
                Err(error) => return is_already_reaped(&error),
            }
        }
    }

    /// Whether the aligner's process group has not yet been torn down by
    /// [`Self::kill`] or released by a clean [`Self::wait`].
    fn process_group_pending(&self) -> bool {
        #[cfg(unix)]
        {
            self.process_group.is_some()
        }
        #[cfg(not(unix))]
        {
            false
        }
    }

    /// For `Drop`: whether the wrapper's exit is confirmed, killing it (bounded)
    /// if it is still running or its state cannot be read.
    fn kill_unless_exited(&mut self) -> bool {
        match self.child.try_wait() {
            // Already exited: stderr fd is closed, relay will finish.
            Ok(Some(_)) => true,
            // Still running: kill() reports `true` only on a confirmed reap
            // (clean exit or ECHILD); `false` means it may still be alive.
            Ok(None) => self.kill(),
            // `try_wait` errored. ECHILD (already reaped) confirms exit. For any
            // other errno the child may still be alive, so fall back to the
            // bounded `kill()` — returning `is_already_reaped` alone here would
            // drop the handle without ever attempting SIGKILL and leak the
            // subprocess. `kill()` re-classifies ECHILD as a confirmed reap.
            Err(err) => is_already_reaped(&err) || self.kill(),
        }
    }

    /// Whether the aligner's `/bin/bash -c` wrapper has exited, reaping it if
    /// so. A wrapper reaped elsewhere (`ECHILD`) counts as exited; any other
    /// `try_wait` error does not.
    pub fn has_exited(&mut self) -> bool {
        match self.child.try_wait() {
            Ok(status) => status.is_some(),
            Err(err) => is_already_reaped(&err),
        }
    }

    /// Return the process ID of the aligner child process.
    #[must_use]
    pub fn pid(&self) -> u32 {
        self.child.id()
    }
}

/// `ECHILD` errno. POSIX assigns it the value 10 on every Unix fgumi runs an
/// external aligner on (Linux, macOS, the BSDs); `std::io` exposes no
/// `ErrorKind` for it, and `nix`/`libc` are not available on all targets here,
/// so we match the raw errno directly. On non-Unix targets nothing returns this
/// value, so [`is_already_reaped`] simply never matches — the conservative
/// "not confirmed reaped" default.
const ECHILD: i32 = 10;

/// Whether a `try_wait` error means the child was *already reaped* (and so its
/// exit is confirmed). Only `ECHILD` qualifies — it is what the OS returns once
/// a SIGCHLD handler (or a prior wait) has already collected the child. Any
/// other errno is unexpected and must NOT be treated as a confirmed exit.
fn is_already_reaped(err: &std::io::Error) -> bool {
    err.raw_os_error() == Some(ECHILD)
}

/// The `/bin/bash -c` wrapper's state as [`AlignerProcess::kill`]'s group
/// teardown sees it.
#[cfg(unix)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum WrapperState {
    /// Still running.
    Running,
    /// Exited and reaped by this process: its PID stays reserved while the
    /// group has members, so the group may be probed and signalled.
    Exited,
    /// Reaped elsewhere (`ECHILD`): its PID may already be reused, so the
    /// group must never be signalled.
    ReapedElsewhere,
    /// `try_wait` failed for another reason; treated as possibly running.
    Unknown,
}

impl Drop for AlignerProcess {
    fn drop(&mut self) {
        // Determine whether the child's exit is *confirmed*. Only then is the
        // stderr relay joined: that thread blocks reading the child's stderr
        // fd, which stays open while the child lives, so a join on a still-alive
        // child could only end at the drain grace, spending it on a relay that
        // cannot finish.
        //
        // Without a clean `wait()`, the aligner's process group is killed even
        // if the wrapper itself has exited: an orphaned grandchild holding the
        // stderr pipe would otherwise keep the relay running, so the join below
        // would wait out the drain grace and leave the relay detached.
        //
        // A group teardown that already ran in `kill()` is authoritative: the
        // wrapper is reaped then, so re-checking it would report exit even
        // when another group member still holds the pipe.
        #[cfg(unix)]
        let recorded = self.group_teardown_confirmed;
        #[cfg(not(unix))]
        let recorded: Option<bool> = None;
        let exit_confirmed = match recorded {
            Some(false) => {
                // A late wrapper exit still needs reaping, but cannot confirm
                // group exit or justify another wait on its pipe holders.
                let deadline = std::time::Instant::now() + KILL_REAP_DEADLINE;
                while matches!(self.child.try_wait(), Ok(None)) {
                    if std::time::Instant::now() >= deadline {
                        break;
                    }
                    thread::sleep(KILL_POLL_INTERVAL);
                }
                false
            }
            Some(true) => true,
            None if self.process_group_pending() => self.kill(),
            None => self.kill_unless_exited(),
        };
        // Join the stderr relay only when exit is confirmed, and then only for
        // the drain grace. If the child could not be confirmed dead (kill timed
        // out / unclassified errno), skip the join: the relay cannot finish
        // while the child holds the fd, so it is left detached at once.
        if exit_confirmed {
            let _ = self.join_stderr_relay();
        }
    }
}

/// Send `signal` to every process in `group`, returning whether the kernel
/// accepted it. Do not log on bounded teardown paths: the sink can block too.
#[cfg(unix)]
fn signal_process_group(group: rustix::process::Pid, signal: rustix::process::Signal) -> bool {
    rustix::process::kill_process_group(group, signal).is_ok()
}

/// Whether `group` has no processes left (`kill(-group, 0)` reports `ESRCH`).
/// A zombie member still counts until its parent reaps it.
#[cfg(unix)]
fn process_group_is_empty(group: rustix::process::Pid) -> bool {
    rustix::process::test_kill_process_group(group) == Err(rustix::io::Errno::SRCH)
}

/// What the stderr relay has captured, shared with [`AlignerProcess`] so it can
/// be read while the relay is still running.
struct StderrCapture {
    /// How many complete lines to keep; 0 disables capture (the relay still
    /// drains and prints the pipe).
    ring_size: usize,
    /// The last `ring_size` complete lines.
    lines: VecDeque<String>,
    /// The bytes of the line being read, before its newline arrives. Part of
    /// the tail too, so a last line without a newline is reported even while a
    /// process the aligner started keeps the pipe open.
    partial: Vec<u8>,
}

impl StderrCapture {
    fn new(ring_size: usize) -> Self {
        Self { ring_size, lines: VecDeque::with_capacity(ring_size), partial: Vec::new() }
    }

    /// End the line in `partial` (stripping a `\r` before a newline), keep it
    /// in the ring, and return it for printing.
    fn finish_line(&mut self, newline: bool) -> String {
        let mut bytes = std::mem::take(&mut self.partial);
        if newline && bytes.last() == Some(&b'\r') {
            bytes.pop();
        }
        let line = String::from_utf8_lossy(&bytes).into_owned();
        if self.ring_size > 0 {
            if self.lines.len() == self.ring_size {
                self.lines.pop_front();
            }
            self.lines.push_back(line.clone());
        }
        line
    }

    /// The captured lines, then the unfinished line, if any.
    fn tail(&self) -> Vec<String> {
        let mut tail: Vec<String> = self.lines.iter().cloned().collect();
        if self.ring_size > 0 && !self.partial.is_empty() {
            tail.push(String::from_utf8_lossy(&self.partial).into_owned());
        }
        tail
    }
}

/// Read lines from `stderr`, print each to fgumi's stderr, and retain the last
/// lines in `capture`, until the child closes its stderr fd. The capture is
/// shared, so what was read so far can be read before this returns: each line
/// is kept before it is printed, and a line still waiting for its newline is
/// kept as it arrives.
fn relay_stderr(stderr: impl std::io::Read, capture: &parking_lot::Mutex<StderrCapture>) {
    let mut reader = BufReader::new(stderr);

    // Byte-oriented drain: `BufRead::lines()` yields `Err` on a non-UTF-8 line,
    // and breaking on it would STOP draining the aligner's stderr — the aligner
    // then blocks once its stderr pipe fills (~64 KiB) and the whole align
    // pipeline deadlocks. Read raw bytes and decode lossily so invalid UTF-8
    // never halts the relay.
    loop {
        let chunk = match reader.fill_buf() {
            Ok(chunk) => chunk,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            // An unrecoverable pipe read error: stop draining.
            Err(_) => break,
        };
        // EOF: the child (and everything that inherited the pipe) closed it.
        if chunk.is_empty() {
            break;
        }
        let (used, line) = {
            let mut capture = capture.lock();
            if let Some(end) = chunk.iter().position(|&b| b == b'\n') {
                capture.partial.extend_from_slice(&chunk[..end]);
                (end + 1, Some(capture.finish_line(true)))
            } else {
                capture.partial.extend_from_slice(chunk);
                (chunk.len(), None)
            }
        };
        reader.consume(used);
        if let Some(line) = line {
            eprintln!("{line}");
        }
    }
    // A last line without a newline.
    let mut capture = capture.lock();
    if !capture.partial.is_empty() {
        let line = capture.finish_line(false);
        drop(capture);
        eprintln!("{line}");
    }
}

// ============================================================================
// AlignerPreset
// ============================================================================

/// Preset aligner configurations with command-building and index validation.
///
/// Four presets — `bwa-mem3` (the project's primary subprocess aligner),
/// `bwa` (legacy reference), `bwa-mem3-inproc` (bwa-mem3 linked directly into
/// fgumi, gated by the `aligner-bwa-mem3` build feature; see
/// `AlignerOptions::resolve`), and `star` (the splice-aware STAR aligner for
/// RNA-seq, which takes its own index directory; see the `star` module). Under
/// runall's `--methylation-mode`, both bwa-mem3 presets align conversion-aware
/// with `bwa-mem3 mem --meth` (EM-seq) or `--meth=taps` (TAPS; see
/// `bwa_mem3_meth_flags` and `AlignerPreset::validate_with_methylation`);
/// `bwa` and `star` have no such mode, so they are rejected under EM-seq and
/// align TAPS reads plain;
/// other bisulfite aligners (e.g. bwameth) go through `--aligner::command`.
///
/// The CLI value for each subprocess variant matches the on-disk binary
/// name (`bwa`, `bwa-mem3`). Variant identifiers (`Bwa`, `BwaMem3`) are
/// chosen so `rename_all = "kebab-case"` round-trips through the
/// identifier ↔ binary-name mapping without divergence: `Bwa` ↔
/// `"bwa"`, `BwaMem3` ↔ `"bwa-mem3"`. `Display` emits the same
/// strings so error messages match what the user typed.
///
/// `BwaMem3InProc` exists in every build regardless of the
/// `aligner-bwa-mem3` feature — so `--help` lists it and clap accepts
/// the value everywhere — and only `resolve` gates whether it can
/// actually be used, so the feature-off error names the precise fix
/// (`--features aligner-bwa-mem3`) instead of clap rejecting an
/// unrecognized value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum AlignerPreset {
    /// Classic BWA aligner (`bwa mem`). CLI value: `bwa`.
    Bwa,
    /// BWA-MEM3 aligner (`bwa-mem3 mem`). CLI value: `bwa-mem3`.
    BwaMem3,
    /// BWA-MEM3 linked in-process (no subprocess, no shell command). CLI
    /// value: `bwa-mem3-inproc`. Requires fgumi built with the
    /// `aligner-bwa-mem3` feature; see `AlignerOptions::resolve`.
    ///
    /// `#[value(name = ...)]` overrides the derived `rename_all =
    /// "kebab-case"` value: `heck`'s kebab-case splits `Mem3In` at the
    /// digit→uppercase boundary too, so the derived value would be
    /// `bwa-mem3-in-proc`, not `bwa-mem3-inproc`.
    #[value(name = "bwa-mem3-inproc")]
    BwaMem3InProc,
    /// STAR, the splice-aware RNA-seq aligner (`STAR --runMode alignReads`).
    /// CLI value: `star`. Needs `--aligner::index <DIR>` (the STAR genome
    /// directory); `--ref` stays the FASTA whose `.dict` names the output's
    /// sequences. Reads are fed as unaligned SAM and returned in input order;
    /// the input must be paired-end (STAR reads it with `--readFilesType SAM
    /// PE`). The wait for STAR's genome load is reported, never failed, by the
    /// deadlock watchdog. `runall` only: `fgumi realign` and `runall --realign`
    /// reject it. STAR's MAPQ is kept: 255 for a unique read (`--outSAMmapqUnique` in
    /// `--aligner::extra-args` changes it), 3 for two loci, 1 for three or
    /// four, 0 for more.
    Star,
}

impl std::fmt::Display for AlignerPreset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.cli_name())
    }
}

impl AlignerPreset {
    /// Default binary name (without path) for this preset: the kebab-case CLI
    /// value (`bwa`, `bwa-mem3`) except for `star`, whose binary is `STAR`.
    /// Internal helper — callers outside this module should use
    /// [`Display`](std::fmt::Display) for the preset's name.
    #[must_use]
    pub(crate) fn binary_name(self) -> &'static str {
        match self {
            Self::Star => "STAR",
            other => other.cli_name(),
        }
    }

    /// The kebab-case CLI value (`bwa`, `bwa-mem3`, `bwa-mem3-inproc`, `star`),
    /// as the user typed it.
    #[must_use]
    fn cli_name(self) -> &'static str {
        match self {
            Self::Bwa => "bwa",
            Self::BwaMem3 => "bwa-mem3",
            Self::BwaMem3InProc => "bwa-mem3-inproc",
            Self::Star => "star",
        }
    }

    /// Whether this preset needs external-binary discovery (`which` on
    /// `PATH`, or a `--aligner-bin` override). `false` only for
    /// [`Self::BwaMem3InProc`], which links bwa-mem3 directly into the
    /// fgumi process and has no separate executable to find; `validate`
    /// skips the binary-discovery step for it but still requires the
    /// index files.
    #[must_use]
    pub(crate) fn requires_binary(self) -> bool {
        !matches!(self, Self::BwaMem3InProc)
    }

    /// Binary name to suggest in an index-missing fix-it hint. Always the
    /// real `bwa-mem3` executable — [`Self::BwaMem3InProc`] links bwa-mem3
    /// in-process for alignment, but the index itself is still built with
    /// the standalone `bwa-mem3 index` command, not a nonexistent
    /// `bwa-mem3-inproc` binary.
    #[must_use]
    fn index_binary_hint(self) -> &'static str {
        match self {
            Self::BwaMem3InProc => "bwa-mem3",
            other => other.binary_name(),
        }
    }

    /// Index file extensions this preset requires alongside the reference
    /// FASTA. All must be present for `validate` to succeed.
    #[must_use]
    pub(crate) fn index_extensions(self) -> &'static [&'static str] {
        match self {
            // BWA classic: amb / ann / bwt / pac / sa
            Self::Bwa => &[".amb", ".ann", ".bwt", ".pac", ".sa"],
            // BWA-MEM3 (subprocess and in-process share the same on-disk index):
            // amb / ann / bwt.2bit.64 / pac. `bwa-mem3 0.4.0+` pac-fetches the
            // reference from `.pac` on demand and no longer writes `.0123` by
            // default (fg-labs/bwa-mem3#177); `.0123` is now opt-in via `index
            // --emit-unpacked-ref` and `mem` ignores any present, so requiring it
            // here rejects a valid 0.4.0 index. `.bwt.2bit.64` is the canonical
            // index sentinel. Indexes built by older bwa-mem3 still carry these
            // four files, so dropping `.0123` is backward-compatible.
            Self::BwaMem3 | Self::BwaMem3InProc => &[".amb", ".ann", ".bwt.2bit.64", ".pac"],
            // STAR's index is a directory of its own (`--aligner::index`),
            // checked by `star::check_index`, not files beside the FASTA.
            Self::Star => &[],
        }
    }

    /// The conversion-aware `mem` flags this preset aligns `methylation` with:
    /// [`bwa_mem3_meth_flags`] for the bwa-mem3 presets, and `None` for `bwa`,
    /// which has no conversion-aware mode (it aligns TAPS plain, and
    /// [`Self::validate_with_methylation`] rejects it for EM-seq).
    #[must_use]
    pub(crate) fn meth_flags(self, methylation: MethylationMode) -> Option<&'static str> {
        match self {
            Self::BwaMem3 | Self::BwaMem3InProc => bwa_mem3_meth_flags(methylation),
            Self::Bwa | Self::Star => None,
        }
    }

    /// Whether this preset consumes `methylation` (`--methylation-mode`): it
    /// aligns conversion-aware ([`Self::meth_flags`]), or is rejected by
    /// [`Self::validate_methylation`] (`bwa` or `star` under EM-seq). Under
    /// TAPS they align plain, so do not consume it.
    #[must_use]
    pub(crate) fn consumes_methylation(self, methylation: MethylationMode) -> bool {
        self.meth_flags(methylation).is_some()
            || (matches!(self, Self::Bwa | Self::Star) && methylation == MethylationMode::EmSeq)
    }

    /// Build the aligner shell command for this preset.
    ///
    /// The command reads FASTQ from `/dev/stdin` and writes its alignments to
    /// stdout, allowing it to be connected to fgumi's pipeline via
    /// stdin/stdout pipes. `bwa-mem3` writes uncompressed BAM (`--bam=0`):
    /// fgumi's single reader thread then takes records zero-copy instead of
    /// parsing SAM text, which at 32 threads kept bwa-mem3 blocked on its
    /// output for most of the run (a fused extract → correct → align ran 36%
    /// faster). `bwa` has no BAM output and writes SAM.
    ///
    /// # Arguments
    ///
    /// * `reference` - Path to the reference FASTA (index must already exist).
    /// * `threads` - Number of threads to pass to the aligner (`-t`).
    /// * `chunk_size` - Chunk size in bases to pass to the aligner (`-K`).
    /// * `binary_override` - Optional explicit path to the aligner binary.
    ///   When `Some`, replaces the preset's default binary name; when
    ///   `None`, the bare binary name is used and the shell resolves it
    ///   via `PATH`.
    /// * `methylation` - The chemistry from `--methylation-mode`. The bwa-mem3
    ///   presets add bwa-mem3's `--meth` flags for EM-seq or TAPS (see
    ///   [`bwa_mem3_meth_flags`]). `bwa` has no conversion-aware mode: it aligns
    ///   TAPS plain, and [`Self::validate_with_methylation`] rejects it for
    ///   EM-seq first.
    // `pub(crate)`, not `pub`: this interpolates `reference` / `binary_override`
    // into a string later run via `/bin/bash -c`, and it does not itself validate
    // those paths — it trusts that `validate_with_methylation` (which runs
    // `check_shell_safe_path`; `validate` and `AlignerOptions::resolve` go through
    // it) was called first. Restricting the method to in-crate callers keeps a direct
    // external caller from bypassing that validation and injecting shell syntax.
    #[must_use]
    pub(crate) fn build_command(
        self,
        reference: &Path,
        threads: usize,
        chunk_size: u64,
        binary_override: Option<&Path>,
        methylation: MethylationMode,
    ) -> String {
        debug_assert!(
            self != Self::Star,
            "star builds its command with star::build_star_command, not build_command"
        );
        let binary = match binary_override {
            Some(path) => path.display().to_string(),
            None => self.binary_name().to_string(),
        };
        let ref_str = reference.display();
        let output = match self {
            Self::BwaMem3 | Self::BwaMem3InProc => " --bam=0",
            // `star` builds its command with `star::build_star_command`.
            Self::Bwa | Self::Star => "",
        };
        let meth = self.meth_flags(methylation).map(|flag| format!(" {flag}")).unwrap_or_default();
        format!("{binary} mem{output}{meth} -p -K {chunk_size} -t {threads} {ref_str} /dev/stdin")
    }

    /// Validate that the aligner binary and required index files are present
    /// for plain (non-bisulfite) alignment.
    ///
    /// Equivalent to the in-crate `validate_with_methylation` with
    /// [`MethylationMode::Disabled`]. Checks, in order:
    /// 1. The reference FASTA itself exists (so a typo in `--ref`
    ///    surfaces with a clean error rather than "index file not
    ///    found").
    /// 2. The reference path does not contain shell-unsafe characters
    ///    (preset mode hands the path through `/bin/bash -c`, so a
    ///    path like `/data/my refs/genome.fa` would be word-split into
    ///    two arguments). Presets where `Self::requires_binary` is
    ///    `false` ([`Self::BwaMem3InProc`]) never build a shell command,
    ///    so they reject only control characters, which would corrupt
    ///    the `@PG CL:` header line the path is written into.
    /// 3. The aligner binary is reachable: `--aligner-bin` override
    ///    path (must be an existing file), else `which::which("<binary_name>")`
    ///    on `PATH`. `--aligner-bin` paths are also checked for
    ///    shell-unsafe characters. Skipped entirely for presets where
    ///    `Self::requires_binary` is `false` ([`Self::BwaMem3InProc`]
    ///    links bwa-mem3 in-process, so there is no binary to find).
    /// 4. All expected index file(s) exist alongside the reference.
    ///
    /// `star` is rejected: its index is a separate directory, which only
    /// `AlignerOptions::resolve` (given `--aligner::index`) can check.
    ///
    /// # Errors
    ///
    /// Returns an error with a fix-it hint if any of the above fail, or for
    /// `star`.
    pub fn validate(self, reference: &Path, binary_override: Option<&Path>) -> Result<()> {
        if self == Self::Star {
            bail!(
                "AlignerPreset::validate cannot check a star index (a separate directory); \
                 resolve the options with AlignerOptions::resolve and --aligner::index"
            );
        }
        self.validate_with_methylation(reference, binary_override, MethylationMode::Disabled)
    }

    /// Validate that the aligner binary and the index files needed to align
    /// for `methylation` (`--methylation-mode`) are present.
    ///
    /// Runs [`Self::validate`]'s checks, except that the index files required
    /// beside the reference depend on the chemistry: when the preset aligns
    /// conversion-aware ([`Self::meth_flags`]: the bwa-mem3 presets under
    /// EM-seq or TAPS) it needs the `bwa-mem3 index --meth` files
    /// [`Self::validate_methylation`] checks; otherwise (disabled, or `bwa`
    /// under TAPS) the plain index ([`Self::index_extensions`]). Under EM-seq,
    /// `bwa` is rejected before anything else, since no index or binary fixes
    /// it.
    ///
    /// # Errors
    ///
    /// Returns an error with a fix-it hint if any check fails.
    pub(crate) fn validate_with_methylation(
        self,
        reference: &Path,
        binary_override: Option<&Path>,
        methylation: MethylationMode,
    ) -> Result<()> {
        self.reject_without_bisulfite_mode(methylation)?;
        self.validate_reference_and_binary(reference, binary_override)?;
        if self.meth_flags(methylation).is_some() {
            self.require_meth_index_files(reference, methylation)
        } else {
            Self::require_index_files(reference, self.index_extensions(), || {
                format!("run `{} index {}`", self.index_binary_hint(), reference.display())
            })
        }
    }

    /// Steps 1-3 of [`Self::validate`]: the reference FASTA exists, its path
    /// (and any `--aligner-bin` path) is safe to interpolate, and the aligner
    /// binary is reachable when the preset needs one.
    fn validate_reference_and_binary(
        self,
        reference: &Path,
        binary_override: Option<&Path>,
    ) -> Result<()> {
        let binary_name = self.binary_name();

        // 1. Reference FASTA presence (separated from the index-loop
        //    error so users see "ref doesn't exist" before "index
        //    missing" when the path is wrong).
        if !reference.is_file() {
            bail!(
                "reference FASTA not found or not a regular file: {} \
                 (preset `{self}` requires the FASTA + its index files)",
                reference.display()
            );
        }

        // 2. Reference path shell-safety (preset mode argv flows through
        //    /bin/bash -c, so paths with spaces, quotes, or shell metas
        //    silently break). The in-process backend loads the index
        //    through the library, so only the `@PG CL:` line constrains it.
        match self {
            Self::Bwa | Self::BwaMem3 => check_shell_safe_path(reference, "--ref")?,
            Self::BwaMem3InProc => check_header_safe_path(reference, "--ref")?,
            // STAR aligns against `--aligner::index`; `--ref` only names the
            // dictionary and never reaches a command line or `@PG`.
            Self::Star => {}
        }

        // 3. Binary discovery — skipped for presets that don't need one.
        if self.requires_binary() {
            match binary_override {
                Some(path) => {
                    check_shell_safe_path(path, "--aligner-bin")?;
                    if !path.is_file() {
                        bail!(
                            "--aligner-bin path is not a regular file: {} \
                             (preset `{self}` requires the binary)",
                            path.display()
                        );
                    }
                    // Verify executability up front so the failure stays
                    // actionable here, instead of being deferred until
                    // `/bin/bash -c` tries to exec the path.
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::PermissionsExt;
                        let mode = std::fs::metadata(path)
                            .with_context(|| {
                                format!(
                                    "reading metadata for --aligner-bin path: {}",
                                    path.display()
                                )
                            })?
                            .permissions()
                            .mode();
                        if mode & 0o111 == 0 {
                            bail!(
                                "--aligner-bin path is not executable: {} \
                                 (preset `{self}` requires an executable binary; \
                                 `chmod +x` it or pass a different path)",
                                path.display()
                            );
                        }
                    }
                }
                None => {
                    which::which(binary_name).with_context(|| {
                        format!(
                            "aligner binary `{binary_name}` not found on PATH \
                             (preset `{self}`; pass --aligner-bin <path> to override)"
                        )
                    })?;
                }
            }
        }

        Ok(())
    }

    /// Require `<prefix><ext>` to be a regular file for every `ext` in
    /// `extensions`, naming the first missing file and appending `fix()` (a
    /// fix-it hint) to the error.
    fn require_index_files<F: Fn() -> String>(
        prefix: &Path,
        extensions: &[&str],
        fix: F,
    ) -> Result<()> {
        for ext in extensions {
            let index_path = append_extension(prefix, ext);
            if !index_path.is_file() {
                bail!("required index file not found: {} ({})", index_path.display(), fix());
            }
        }
        Ok(())
    }

    /// Reject `bwa` and `star` under EM-seq: neither has a bisulfite mode. (They
    /// align TAPS plain, since TAPS converts only the methylated cytosines.)
    fn reject_without_bisulfite_mode(self, methylation: MethylationMode) -> Result<()> {
        if matches!(self, Self::Bwa | Self::Star) && methylation == MethylationMode::EmSeq {
            bail!(
                "--methylation-mode {EM_SEQ_MODE_NAME} needs a bisulfite-aware aligner, and `{self}` has \
                 none; use --aligner::preset bwa-mem3 or bwa-mem3-inproc (which align with \
                 bwa-mem3's --meth), or --aligner::command"
            );
        }
        Ok(())
    }

    /// Index file extensions the original (un-converted) `<ref>` prefix needs
    /// when aligning with `bwa-mem3 mem --meth` (EM-seq) or `--meth=taps`
    /// (TAPS).
    ///
    /// `bwa-mem3 mem --meth` seeds against the converted `<ref>.meth` index and
    /// reads only the original reference's contigs and sequence (`.amb`,
    /// `.ann`, `.pac`), never the plain `<ref>.bwt.2bit.64`. A
    /// `bwa-mem3 index --meth --meth-only` index omits that file, so the
    /// subprocess preset does not require it: its list is
    /// [`Self::index_extensions`] minus the plain FM-index `.bwt.2bit.64`.
    ///
    /// The in-process backend still requires the full plain index. This
    /// mirrors bwa-mem3-rs's index check in `BwaIndex::load_meth`, which
    /// requires the plain `<ref>.bwt.2bit.64` on the original prefix before
    /// loading a dual index, so it cannot load a `--meth-only` index. That
    /// list is a private constant in bwa-mem3-rs, so it cannot be imported
    /// here.
    #[must_use]
    fn meth_original_index_extensions(self) -> &'static [&'static str] {
        match self {
            Self::BwaMem3 => &[".amb", ".ann", ".pac"],
            // `Bwa` and `Star` are unreachable here: they have no `--meth` flags
            // (`meth_flags`), so no dual index is ever checked for them.
            Self::BwaMem3InProc | Self::Bwa | Self::Star => self.index_extensions(),
        }
    }

    /// The fix-it hint for a missing file of the `bwa-mem3 index --meth` dual
    /// index. The subprocess preset also accepts a `--meth-only` index, but the
    /// hint suggests the full dual index: a `--meth-only` index cannot align
    /// plain (duplex re-alignment, runs without `--methylation-mode`), and the
    /// in-process backend needs the full dual index
    /// ([`Self::meth_original_index_extensions`]). From 0.18.0, bwa-mem3's
    /// `index --meth` also builds the plain index and refuses to replace
    /// existing index files without `--force`, which older releases do not
    /// accept, so the hint passes `--force` and says when to drop it.
    ///
    /// TAPS converts only the methylated cytosines, so its reads also align
    /// plain; the TAPS hint names that fallback.
    fn meth_index_hint(self, reference: &Path, methylation: MethylationMode) -> String {
        let binary = self.index_binary_hint();
        let mode = methylation_mode_name(methylation);
        let flag = bwa_mem3_meth_flag_name(methylation);
        let reference = reference.display();
        let run = format!(
            "run `{binary} index --meth --force {reference}`, which also rebuilds the plain \
             index; drop `--force` for bwa-mem3 before 0.18.0, which does not accept it"
        );
        let plain_fallback = if methylation == MethylationMode::Taps {
            "; TAPS reads also align plain without it, with --aligner::preset bwa (which \
             needs a bwa index) or an --aligner::command running plain `bwa-mem3 mem`"
        } else {
            ""
        };
        match self {
            Self::BwaMem3InProc => format!(
                "--methylation-mode {mode} aligns with bwa-mem3 {flag}, and the in-process \
                 backend needs the full dual index, including the plain \
                 `{reference}.bwt.2bit.64`; {run}{plain_fallback}"
            ),
            // `Bwa` and `Star` are unreachable here: they have no `--meth` flags
            // (`meth_flags`), so no dual index is ever checked for them.
            Self::BwaMem3 | Self::Bwa | Self::Star => format!(
                "--methylation-mode {mode} aligns with bwa-mem3 {flag}, which needs the dual \
                 index; {run}{plain_fallback}"
            ),
        }
    }

    /// Validate that this preset can align for `methylation` (`--methylation-mode`).
    ///
    /// EM-seq reads must be aligned conversion-aware: every unmethylated C reads
    /// as T, so a plain aligner scores most cytosines as mismatches, heavily
    /// converted reads map worse, and the methylation calls are biased. TAPS
    /// converts only the methylated cytosines, so its reads align plain too, but
    /// bwa-mem3's `--meth=taps` places them better. The bwa-mem3 presets align
    /// both with bwa-mem3's `--meth` flags (see `bwa_mem3_meth_flags`) against a
    /// dual index built by `bwa-mem3 index --meth`: the converted seed index
    /// `<ref>.meth.{amb,ann,bwt.2bit.64,pac}` plus the original reference's
    /// `<ref>.{amb,ann,pac}`. This checks for those files. The plain
    /// `<ref>.bwt.2bit.64` is not required for the `bwa-mem3` preset, so a
    /// `bwa-mem3 index --meth --meth-only` index (which omits it) is accepted;
    /// `bwa-mem3-inproc` still requires it (see `meth_original_index_extensions`).
    ///
    /// `bwa` has no conversion-aware mode: it is rejected for EM-seq, and aligns
    /// TAPS plain (no `.meth` index needed). A no-op for
    /// [`MethylationMode::Disabled`].
    ///
    /// # Errors
    ///
    /// Returns an error if the preset is `bwa` under EM-seq, or a bwa-mem3 preset
    /// under EM-seq or TAPS is missing a dual index file.
    ///
    /// # Pairing with `validate`
    ///
    /// This checks the dual index files only, and is meant to be paired with
    /// [`Self::validate`], which requires the plain index (including
    /// `<ref>.bwt.2bit.64`). The pair therefore rejects a `--meth-only` index.
    /// In-crate resolution uses `validate_with_methylation` instead, which
    /// picks the required index files by chemistry and accepts a `--meth-only`
    /// index for the `bwa-mem3` preset under EM-seq or TAPS.
    pub fn validate_methylation(
        self,
        reference: &Path,
        methylation: MethylationMode,
    ) -> Result<()> {
        self.reject_without_bisulfite_mode(methylation)?;
        if self.meth_flags(methylation).is_none() {
            return Ok(());
        }
        self.require_meth_index_files(reference, methylation)
    }

    /// Require the `bwa-mem3 index --meth` files this preset reads under
    /// `methylation` (EM-seq or TAPS): the original reference's files
    /// ([`Self::meth_original_index_extensions`]) and the converted
    /// `<ref>.meth` seed index.
    fn require_meth_index_files(
        self,
        reference: &Path,
        methylation: MethylationMode,
    ) -> Result<()> {
        let fix = || self.meth_index_hint(reference, methylation);
        Self::require_index_files(reference, self.meth_original_index_extensions(), fix)?;
        let seed_prefix = append_extension(reference, BWA_MEM3_METH_INDEX_SUFFIX);
        Self::require_index_files(&seed_prefix, self.index_extensions(), fix)
    }
}

/// The suffix bwa-mem3 appends to the reference for the converted seed index
/// of a `bwa-mem3 index --meth` dual index (`<ref>.meth.bwt.2bit.64`, ...).
pub(crate) const BWA_MEM3_METH_INDEX_SUFFIX: &str = ".meth";

/// bwa-mem3's `mem` flags for `methylation`, as one space-separated string to
/// splice into a command line: `--meth --meth-scoring genomic` (conversion-aware
/// alignment, EM-seq chemistry, genomic scoring) for EM-seq,
/// `--meth=taps --meth-scoring neutral` for TAPS, `None` when disabled. This is
/// the one place that decides whether a chemistry aligns conversion-aware with
/// the bwa-mem3 presets; the index check, the in-process options and the `@PG`
/// record all follow it.
///
/// The scoring is bwa-mem3's default for each chemistry from 0.17.0 (EM-seq
/// was collapsed, `-B 2`, before), and the in-process backend sets it
/// explicitly too. Spelling it out pins the scoring mode, though not every
/// other `--meth` default, for any subprocess `bwa-mem3` from 0.4.0 (the first
/// release with `--meth-scoring`), so it scores as the in-process backend does
/// and as realign rescores.
///
/// TAPS converts only the methylated cytosines, so its reads also align plain,
/// but `--meth=taps` (same `--meth` index, seeding and `XR`/`XG` tags, neutral
/// conversion scoring by default) places them better. It must be `--meth=taps`:
/// a bare `--meth` inverts every `XM` call on TAPS reads, and bwa-mem3 binds the
/// chemistry only in the `=` form (in `--meth taps`, `taps` is a positional).
#[must_use]
pub(crate) fn bwa_mem3_meth_flags(methylation: MethylationMode) -> Option<&'static str> {
    match methylation {
        MethylationMode::EmSeq => Some("--meth --meth-scoring genomic"),
        MethylationMode::Taps => Some("--meth=taps --meth-scoring neutral"),
        MethylationMode::Disabled => None,
    }
}

/// The bare `--meth` flag (without `--meth-scoring`) `methylation` aligns
/// with, for messages: `--meth` for EM-seq, `--meth=taps` for TAPS.
///
/// # Panics
///
/// On [`MethylationMode::Disabled`], which aligns with no `--meth` flag; the
/// only caller, `AlignerPreset::meth_index_hint`, runs only for a chemistry
/// [`bwa_mem3_meth_flags`] aligns conversion-aware.
fn bwa_mem3_meth_flag_name(methylation: MethylationMode) -> &'static str {
    match methylation {
        MethylationMode::EmSeq => "--meth",
        MethylationMode::Taps => "--meth=taps",
        MethylationMode::Disabled => {
            unreachable!("no dual index is checked for --methylation-mode disabled")
        }
    }
}

/// The `--methylation-mode` value that selects `methylation`, for messages.
///
/// # Panics
///
/// On [`MethylationMode::Disabled`], which no `--methylation-mode` value
/// selects; the only caller, `AlignerPreset::meth_index_hint`, runs only for a
/// chemistry [`bwa_mem3_meth_flags`] aligns conversion-aware.
fn methylation_mode_name(methylation: MethylationMode) -> &'static str {
    match methylation {
        MethylationMode::EmSeq => EM_SEQ_MODE_NAME,
        MethylationMode::Taps => TAPS_MODE_NAME,
        MethylationMode::Disabled => unreachable!(
            "no dual index is checked for --methylation-mode disabled, which aligns with \
             no conversion-aware flag"
        ),
    }
}

// ============================================================================
// Command-mode template substitution
// ============================================================================

/// Substitute `{ref}` and `{threads}` placeholders in a free-form
/// aligner command template.
///
/// Used by `--aligner::command "..."` (command mode) to fill in the
/// reference path and thread count from runall's `--ref` / `--threads`
/// flags before the command is handed to `/bin/bash -c`.
///
/// Substitution is literal text replacement — there is no shell
/// escaping. A user who wants to put a literal `{ref}` in their
/// command must work around the collision themselves.
///
/// # Errors
///
/// Returns an error if `template` does not contain `{ref}` — without
/// a reference path the aligner has nothing to align against and
/// would fail at run time with a less-clear error.
pub fn substitute_template(template: &str, reference: &Path, threads: usize) -> Result<String> {
    if !template.contains("{ref}") {
        bail!(
            "--aligner::command template does not contain `{{ref}}`; \
             the aligner needs a reference path. Example: \
             \"bwa-mem3 mem --bam=0 -p -K 150000000 -t {{threads}} {{ref}} /dev/stdin\""
        );
    }
    // Single left-to-right pass so injected text is never re-scanned. Two
    // sequential `str::replace` calls would corrupt a reference path that
    // literally contains `{threads}`: the `{ref}` pass injects it, then the
    // `{threads}` pass rewrites the path's own `{threads}` into the thread
    // count. Consuming each placeholder exactly once and never re-examining
    // emitted text closes that collision in both directions.
    let ref_str = reference.display().to_string();
    let threads_str = threads.to_string();
    let mut out = String::with_capacity(template.len());
    let mut rest = template;
    while !rest.is_empty() {
        if let Some(after) = rest.strip_prefix("{ref}") {
            out.push_str(&ref_str);
            rest = after;
        } else if let Some(after) = rest.strip_prefix("{threads}") {
            out.push_str(&threads_str);
            rest = after;
        } else {
            let mut chars = rest.chars();
            // `rest` is non-empty (loop guard), so `next()` is always `Some`;
            // the `None` arm is unreachable but keeps this panic-free.
            match chars.next() {
                Some(c) => out.push(c),
                None => break,
            }
            rest = chars.as_str();
        }
    }
    Ok(out)
}

// ============================================================================
// AlignerOptions — runall CLI surface
// ============================================================================

/// Default chunk size in bases per aligner batch (`-K` flag). Matches the
/// `bwa mem -K` examples in the bwa-mem3 documentation; large enough that
/// the aligner sees full batches but small enough that ~2 batches in
/// flight fit in a few GB of RAM.
pub const DEFAULT_ALIGNER_CHUNK_SIZE: u64 = 150_000_000;

/// `--aligner::dedup-reads`: whether the in-process bwa-mem3 backend skips
/// seeding exact duplicate read pairs, copying each one's alignment from one
/// copy (its representative) in its `-K` cohort (bwa-mem3's `--dedup-reads`
/// memo). Output is byte-identical either way; `on` saves the seeding work on
/// PCR-duplicate-rich input (UMI and amplicon libraries).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum DedupReads {
    /// Seed each distinct read pair once per cohort (the default).
    On,
    /// Seed every read pair.
    Off,
}

/// Default sub-batch size, in templates, for the in-process bwa-mem3
/// backend's internal batching. Only meaningful with
/// `--aligner::preset bwa-mem3-inproc`; see
/// `AlignerOptions::sub_batch_templates` and
/// [`ResolvedBackend::InProcessBwaMem3`].
///
/// One sub-batch of pairs fills exactly one of bwa-mem3's `BATCH_SIZE` read
/// batches ([`bwa_mem3_rs::kernel_batch_size`]: 1024 reads on aarch64, 512
/// elsewhere), so every seed, extension and mate-rescue kernel call in the
/// shim runs on a full batch, as the `bwa-mem3` CLI's own workers do. A
/// half-filled batch costs ~1.2% CPU on Graviton4. Larger sub-batches gain
/// nothing and leave fewer work items to balance across the pool at high
/// thread counts.
#[cfg(feature = "aligner-bwa-mem3")]
pub(crate) fn default_sub_batch_templates() -> usize {
    bwa_mem3_rs::kernel_batch_size() / 2
}

/// Per-stage aligner tuning knobs.
///
/// Annotated with `#[multi_options("aligner", "Aligner Options")]` so
/// `runall` exposes each field as `--aligner::<flag>` via the generated
/// `MultiAlignerOptions` companion struct. Mutual-exclusion between
/// `preset` and `command` is enforced inside `Self::resolve` rather
/// than at clap-parse time — the macro doesn't rewrite clap's
/// `conflicts_with` field-name references when the prefixed identifiers
/// shift, so we validate logically after `MultiAlignerOptions::validate`.
///
/// Field shapes:
/// - `preset` — `Option<AlignerPreset>` so clap parses it as
///   `--aligner::preset {bwa-mem3|bwa}` (no default; if unset and
///   `command` is also unset, `Self::resolve` errors).
/// - `command` — `Option<String>` for the free-form mode.
/// - `threads` — `Option<usize>` so `Self::resolve` can default it
///   from `std::thread::available_parallelism()` at runtime (the macro
///   can't represent "all cores" as a const default).
/// - `chunk_size` — `u64` with a const default; drives the `-K` flag
///   in preset mode and the Step-1 batch size in both modes.
/// - `sub_batch_templates` — `Option<usize>`, hidden from `--help`.
///   Only meaningful with `--aligner::preset bwa-mem3-inproc`;
///   `Self::resolve` rejects it with any other preset or with command
///   mode, and defaults it from `default_sub_batch_templates` when unset.
/// - `dedup_reads` — `Option<DedupReads>`, likewise in-process only; `None`
///   means `on`.
/// - `max_read_ahead` — `Option<u64>` bytes; `None` means the default ceiling,
///   which depends on the host's memory (resolved when the align stage is wired).
///   Subprocess aligners only (preset or command).
/// - `index` / `extra_args` — `--aligner::preset star` only (its genome
///   directory and extra STAR parameters); `Self::resolve` rejects them
///   otherwise.
#[multi_options("aligner", "Aligner Options")]
#[derive(Args, Debug, Clone)]
pub struct AlignerOptions {
    /// Named aligner preset. Mutually exclusive with `--aligner::command`.
    /// Every preset but `bwa-mem3-inproc` runs the aligner as a subprocess in
    /// its own process group, which a Ctrl-C or `timeout` sent to fgumi does
    /// not reach. After fgumi stops, the aligner may continue running until it
    /// exits or accesses the closed stdin or stdout pipe.
    #[arg(long = "preset", value_enum)]
    pub preset: Option<AlignerPreset>,

    /// Free-form aligner command template. Supports `{ref}` and
    /// `{threads}` placeholders. Mutually exclusive with
    /// `--aligner::preset` and every other `--aligner::*` flag except
    /// `--aligner::chunk-size`, `--aligner::input-format`, and
    /// `--aligner::max-read-ahead`. The command runs in its own process group,
    /// so it must not prompt on the terminal (e.g. use `ssh -o BatchMode=yes`,
    /// `sudo -n`): a terminal read from it stops the aligner and stalls the run.
    /// A Ctrl-C or `timeout` sent to fgumi does not reach it either. After
    /// fgumi stops, the command may continue running indefinitely; EOF or a
    /// broken pipe affects it only when it accesses the corresponding stdin
    /// or stdout pipe.
    #[arg(long = "command")]
    pub command: Option<String>,

    /// Thread count passed to the aligner via `-t` (`--runThreadN` for `star`;
    /// preset mode only).
    /// Defaults to runall's top-level `--threads` value when unset.
    /// The aligner uses these threads concurrently with the rest of
    /// the AAM pipeline; oversubscription is the OS scheduler's
    /// problem.
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    /// Bases per aligner batch (the bwa presets' `-K` flag; also sizes the
    /// in-flight budget in every subprocess mode, which `star` raises to cover
    /// STAR's read-ahead).
    #[arg(long = "chunk-size", default_value_t = DEFAULT_ALIGNER_CHUNK_SIZE)]
    pub chunk_size: u64,

    /// Sub-batch size, in templates, for the in-process bwa-mem3
    /// backend's internal batching. Hidden from `--help` — only valid
    /// with `--aligner::preset bwa-mem3-inproc`; `Self::resolve` rejects
    /// it for every other preset or for command mode, and caps it at
    /// `MAX_SUB_BATCH_TEMPLATES`. `None` uses `default_sub_batch_templates()`.
    #[arg(long = "sub-batch-templates", hide = true)]
    pub sub_batch_templates: Option<usize>,

    /// Skip seeding exact duplicate read pairs (same bases in both mates) and
    /// copy their alignment from one copy in the cohort, as bwa-mem3's
    /// `--dedup-reads` does. Output is byte-identical either way; `on` is faster
    /// on PCR-duplicate-rich input. Only valid with `--aligner::preset
    /// bwa-mem3-inproc`. [default: on]
    #[arg(long = "dedup-reads", value_enum)]
    pub dedup_reads: Option<DedupReads>,

    /// Format of the unaligned reads written to the aligner's stdin: `fastq`
    /// (interleaved FASTQ) or `sam` (unaligned SAM: an `@HD` line, then one
    /// record per read with FLAG 77 for R1, 141 for R2 and 4 for an unpaired
    /// read, no tags), for aligners that read SAM input. Only valid with
    /// `--aligner::command`; the bwa presets read FASTQ and reject
    /// `--aligner::input-format sam`, and `star` reads SAM. [default: fastq]
    #[arg(long = "input-format", value_enum)]
    pub input_format: Option<UnalignedReadFormat>,

    /// The aligner's index directory, for an aligner that does not index
    /// beside the `--ref` FASTA. Required with `--aligner::preset star`, and
    /// only valid with it: the STAR genome directory, built by `STAR --runMode
    /// genomeGenerate` from the same FASTA as `--ref`. Its `chrName.txt` and
    /// `chrLength.txt`, and the combined `chrNameLength.txt`, must match the
    /// `--ref` `.dict`; mandatory index files are checked before STAR starts.
    /// Index and scratch paths must be UTF-8 (spaces are supported).
    /// Not supported by `fgumi realign`.
    #[arg(long = "index", value_name = "DIR")]
    pub index: Option<PathBuf>,

    /// Extra parameters appended to the STAR command line, as one
    /// whitespace-separated string, e.g. `--outSAMmapqUnique 60
    /// --outSAMattributes NH HI AS nM NM MD`. Only valid with `--aligner::preset
    /// star`. Rejected: the parameters the preset sets itself, `--twopassMode`
    /// (two-pass re-reads the input), `--outFilterType BySJout`,
    /// `--outSAMprimaryFlag AllBestScore`, `--outSAMfilter`, `--readMapNumber`,
    /// `--parametersFiles`, `--outSAMreadID` other than `Standard`,
    /// `--outSAMmode` other than `Full`, `--genomeLoad` other than `NoSharedMemory`,
    /// `--outSAMmultNmax` other than -1 or a
    /// positive limit, flag masks that change pairing, mapping, strand or
    /// alignment-classification bits, and the parameters that write side
    /// files (`--quantMode`, `--outReadsUnmapped`, chimeric output, `--soloType`,
    /// `--sjdbInsertSave All`). Both `--name value` and `--name=value` are
    /// accepted and normalized. Words are split on whitespace: quotes are
    /// rejected and a value cannot contain a space. A
    /// `--limitIObufferSize` here also sizes the in-flight budget. Not supported
    /// by `fgumi realign`.
    #[arg(long = "extra-args", value_name = "ARGS", allow_hyphen_values = true)]
    pub extra_args: Option<String>,

    /// Maximum seconds from spawning STAR to its first output (1 or greater).
    /// Only valid with `--aligner::preset star`. Expiry fails the run and kills
    /// the subprocess. Increase this for a slow genome load. [default: 900]
    #[arg(long = "startup-timeout", value_name = "SECONDS", value_parser = clap::value_parser!(u64).range(1..))]
    pub startup_timeout: Option<u64>,
    /// Most memory the unmapped reads sent to the aligner but not yet returned
    /// by it may take, for an aligner that reads further ahead before its first
    /// output than `--aligner::chunk-size` provides for. fgumi then admits more
    /// input after a short pause instead of waiting on the aligner, up to this
    /// limit; past it, fgumi waits for the aligner and fails only if no output
    /// comes for five minutes. A value below the in-flight budget
    /// `--aligner::chunk-size` sets lowers that budget too (with a warning).
    /// Takes a size such as `8GiB` (a plain number is MiB). Not valid with
    /// `--aligner::preset bwa-mem3-inproc`. [default: a quarter of the
    /// machine's memory, and at least the in-flight budget]
    #[arg(long = "max-read-ahead", value_name = "SIZE", value_parser = parse_max_read_ahead)]
    pub max_read_ahead: Option<u64>,

    /// The chemistry to align for. `#[arg(skip)]`: it is runall's top-level
    /// `--methylation-mode`, which runall copies here so the align stage and
    /// the consensus callers it feeds agree. A bwa-mem3 preset aligns with
    /// bwa-mem3's `--meth` under EM-seq and `--meth=taps` under TAPS.
    #[arg(skip)]
    pub methylation_mode: MethylationMode,
}

/// Hand-rolled `Default` impl. **Must** match each field's clap
/// `default_value_t` exactly: the `multi_options` macro emits
/// `default_value_t = AlignerOptions::default().<field>` on the
/// generated `MultiAlignerOptions`, so a `#[derive(Default)]` here
/// would set `chunk_size = 0` and propagate that as `--aligner::chunk-size
/// [default: 0]` — and from there as `bwa-mem3 -K 0` at runtime. The
/// macro's smoke tests document this trap explicitly; this impl is
/// the antidote.
impl Default for AlignerOptions {
    fn default() -> Self {
        Self {
            preset: None,
            command: None,
            threads: None,
            chunk_size: DEFAULT_ALIGNER_CHUNK_SIZE,
            sub_batch_templates: None,
            dedup_reads: None,
            input_format: None,
            index: None,
            extra_args: None,
            startup_timeout: None,
            max_read_ahead: None,
            methylation_mode: MethylationMode::Disabled,
        }
    }
}

/// Parse `--aligner::max-read-ahead`: a memory size as `--max-memory` takes it
/// (a plain number is MiB; zero is rejected).
fn parse_max_read_ahead(s: &str) -> std::result::Result<u64, String> {
    crate::validation::parse_memory_size(s).map_err(|e| e.to_string())
}

/// The alignment backend a [`ResolvedAligner`] will run: an external
/// aligner subprocess reached via a shell command, or — behind the
/// `aligner-bwa-mem3` build feature — bwa-mem3 linked directly into the
/// fgumi process.
///
/// Exists in every build regardless of the feature (so [`AlignerOptions::resolve`]
/// can emit the feature-off error uniformly before constructing anything, and
/// the align stage's `backend_for` has one match with no `#[cfg]` on the enum
/// itself).
#[derive(Debug, Clone)]
pub(crate) enum ResolvedBackend {
    /// `--aligner::preset` (subprocess presets) or `--aligner::command`:
    /// the shell command to spawn via [`AlignerProcess::spawn`].
    Subprocess {
        command: String,
        /// Whether the aligner's output may carry bwa's mid-pair split (a
        /// pair aligned as two unpaired reads because a `-K` chunk cut fell
        /// between them). `true` only for the subprocess presets, which run
        /// `mem -p -K` and so produce exactly that shape; a free-form
        /// `--aligner::command` gets the loud "multiple primaries" error
        /// instead, because the same shape from an aligner that never pairs
        /// would otherwise split every pair silently.
        accept_mid_pair_split: bool,
        /// Format of the reads written to the aligner's stdin
        /// (`--aligner::input-format`). FASTQ for the bwa presets, SAM for
        /// `star`.
        input_format: UnalignedReadFormat,
        /// `--aligner::max-read-ahead`: the ceiling the in-flight budget may
        /// grow to, or `None` for the default. An explicit value below the
        /// budget `--aligner::chunk-size` derives also lowers the budget.
        max_read_ahead: Option<u64>,
        /// What a preset knows about running its aligner beyond the command.
        runtime: SubprocessRuntime,
    },
    /// `--aligner::preset bwa-mem3-inproc`. Constructed only when
    /// `resolve` runs under the `aligner-bwa-mem3` feature (it bails
    /// before reaching this arm otherwise), but the variant itself is
    /// unconditional — see the enum-level doc.
    #[allow(dead_code)] // constructed only with the feature; matched everywhere
    InProcessBwaMem3 {
        reference: PathBuf,
        sub_batch_templates: usize,
        /// `--aligner::dedup-reads` resolved (`on` unless set to `off`).
        dedup_reads: bool,
        /// The bisulfite chemistry to align for (`--methylation-mode`), or
        /// [`MethylationMode::Disabled`].
        methylation: MethylationMode,
        /// bwa's `-h <primary>,<alt>` (`MemOpts::set_xa_max_hits`), or `None`
        /// for bwa-mem3's default. Set only by the realign stage, via
        /// [`ResolvedAligner::with_xa_max_hits`].
        xa_max_hits: Option<(u32, u32)>,
    },
}

/// What a subprocess preset knows about running its aligner beyond the command
/// line. The default (every field off) is what the bwa presets and
/// `--aligner::command` run with.
#[derive(Debug, Clone, Default)]
pub(crate) struct SubprocessRuntime {
    /// STAR's paired SAM input requires a primary R1 and R2 on every template.
    /// Off for bwa and free-form commands, which may accept single reads.
    pub require_complete_pairs: bool,
    /// A scratch directory the command writes into, kept until the aligner has
    /// exited and then removed (`star`'s logs and per-chunk files).
    pub work_dir: Option<AlignerScratchDir>,
    /// A floor for the in-flight unmapped budget, in fgumi heap bytes: the
    /// aligner's known read-ahead before its first output (`star`, sized by
    /// `star::star_in_flight_budget`). The budget is the larger of this and the
    /// `--aligner::chunk-size`-derived one.
    pub min_in_flight_budget: u64,
    /// Maximum startup wait for a known slow-starting aligner (`star`). Until
    /// its first output, the wait is external to the deadlock watchdog, but
    /// exceeding this deadline fails the run. `None` uses the normal watchdog.
    pub startup_timeout: Option<std::time::Duration>,
}

/// Result of [`AlignerOptions::resolve`] — a ready-to-construct aligner
/// backend plus the parameters the downstream chain needs.
///
/// `pub(crate)` because only the runall AAM dispatch consumes it;
/// promote to `pub` if a cross-crate caller materializes.
///
/// Every field is read when the chain builder wires the AAM stage
/// (`chains/builder.rs`): `backend` seeds the align step via
/// the align stage's `backend_for`, `chunk_size` derives the step's
/// `in_flight_unmapped_budget` (via `in_flight_budget_for_chunk_size`), and
/// `mode`/`threads` are info-logged.
#[derive(Debug, Clone)]
pub(crate) struct ResolvedAligner {
    /// The backend to construct (by the align stage's `backend_for`).
    pub backend: ResolvedBackend,
    /// The aligner's `-K` chunk size, in bases. The chain builder converts
    /// it into the align backend's in-flight byte budget so the
    /// resident unmapped backlog stays proportional to one aligner chunk.
    pub chunk_size: u64,
    /// Resolved thread count (default-substituted for preset mode;
    /// `None` if command mode let the user hardcode their own, or if
    /// the in-process backend shares runall's single thread budget).
    pub threads: Option<usize>,
    /// Which mode produced this — preset (with which preset) or
    /// command. Used in info-logging.
    pub mode: ResolvedAlignerMode,
}

impl ResolvedAligner {
    /// Set bwa's `-h <primary>,<alt>`: how many alternative hits are reported
    /// in `XA` (alternative hits on ALT contigs use `alt`). Used only by the
    /// realign stage, and only when the user sets it (`option` names the flag
    /// they used, for messages), so the raw-read align stage and an unset value
    /// leave the aligner's own default and the command untouched.
    ///
    /// The subprocess backend gets `-h <primary>,<alt>` inserted right after the
    /// command's first `mem` word (so it precedes the positional arguments
    /// whatever `getopt` the aligner uses), or appended, with a warning, when the
    /// command has none. A command that already passes `-h` is left as is, with a
    /// warning that its own `-h` wins. Preset paths are shell-safe (no
    /// whitespace), so that word is always the subcommand. The in-process backend
    /// sets it on its `MemOpts`.
    #[must_use]
    pub(crate) fn with_xa_max_hits(mut self, primary: u32, alt: u32, option: &str) -> Self {
        match &mut self.backend {
            ResolvedBackend::Subprocess { command, .. } => {
                let flag = format!("-h {primary},{alt}");
                if command_passes_h(command) {
                    log::warn!(
                        "realign: --aligner::command already passes bwa's -h, which wins: \
                         {option} {primary},{alt} is not applied"
                    );
                } else {
                    *command = insert_after_mem(command, &flag).unwrap_or_else(|| {
                        log::warn!(
                            "realign: --aligner::command has no `mem` word to put bwa's `{flag}` \
                             ({option}) after, so it is appended to the end of the command; \
                             make sure the aligner still receives it"
                        );
                        format!("{command} {flag}")
                    });
                }
            }
            ResolvedBackend::InProcessBwaMem3 { xa_max_hits, .. } => {
                *xa_max_hits = Some((primary, alt));
            }
        }
        self
    }
}

/// Whether `command` passes bwa's `-h` among the words after its `mem` word (the
/// whole command when it has none), up to the first `|`, `;`, `&&` or `||`:
/// as `-h` or with its value attached (`-h3`), and also at the end of a cluster
/// of bwa-mem3's argument-less short flags (`-Yh`, `-ph5`), as `getopt` reads it.
/// Words are split on whitespace without interpreting shell quoting (the tree has
/// no shell-word lexer), so a `-h` inside a quoted argument such as `-R '... -h
/// ...'` counts; `fgumi realign --help` documents this.
fn command_passes_h(command: &str) -> bool {
    /// bwa-mem3 `mem`'s short flags that take no argument, so `-h` may follow
    /// them in one cluster.
    const NO_ARG_FLAGS: &str = "aCMYpPSVj5qu";
    let passes_h = |word: &str| {
        let Some(cluster) = word.strip_prefix('-').filter(|c| !c.starts_with('-')) else {
            return false;
        };
        let flags_end = cluster.find(|c: char| !NO_ARG_FLAGS.contains(c)).unwrap_or(cluster.len());
        cluster[flags_end..]
            .strip_prefix('h')
            .is_some_and(|v| v.is_empty() || v.starts_with(|c: char| c.is_ascii_digit()))
    };
    mem_option_words(command).any(passes_h)
}

/// Insert `flag` after the first `mem` word of `command` ([`mem_word_end`]), or
/// `None` when there is none.
fn insert_after_mem(command: &str, flag: &str) -> Option<String> {
    let end = mem_word_end(command)?;
    Some(format!("{} {flag}{}", &command[..end], &command[end..]))
}

/// How a [`ResolvedAligner`] was produced. The AAM chain builder logs the mode
/// (`chains/builder.rs`), so the `Preset` payload is consumed only through the
/// derived `Debug` — the `allow(dead_code)` suppresses the never-read lint,
/// which fires for a field read solely via a derived trait.
#[allow(dead_code)] // `Preset`'s payload is diagnostic-only (Debug-logged).
#[derive(Debug, Clone, Copy)]
pub(crate) enum ResolvedAlignerMode {
    /// `--aligner::preset <preset>`. Carries the chosen preset so the
    /// caller can log it.
    Preset(AlignerPreset),
    /// `--aligner::command "..."`.
    Command,
}

/// The largest `--aligner::sub-batch-templates` the in-process backend accepts.
/// Each sub-batch pre-sizes a `Vec` of this many templates (about 13 MB at the
/// cap) and counts its pairs and single-end reads in `u32`s. The cap is 128 to
/// 256 times the default (`default_sub_batch_templates`, 512 or 256), far past
/// any size that helps throughput.
const MAX_SUB_BATCH_TEMPLATES: usize = 1 << 16;

/// The largest `--aligner::chunk-size` a preset accepts: `i32::MAX`, since bwa
/// and bwa-mem3 parse `-K` with `atoi` into an `int`.
const MAX_PRESET_CHUNK_SIZE: u64 = 2_147_483_647;

impl AlignerOptions {
    /// Reject the flags that only one kind of backend reads, before the
    /// preset/command dispatch in [`Self::resolve`], so each is rejected
    /// uniformly regardless of mode and of whether this binary was built with
    /// `aligner-bwa-mem3`.
    fn check_backend_specific_flags(&self) -> Result<()> {
        // `--aligner::sub-batch-templates` is a hidden, in-process-only knob.
        if self.sub_batch_templates.is_some()
            && !matches!(self.preset, Some(AlignerPreset::BwaMem3InProc))
        {
            bail!(
                "--aligner::sub-batch-templates is only valid with --aligner::preset \
                 bwa-mem3-inproc; it configures the in-process backend's internal \
                 batching and has no effect on the subprocess aligner or \
                 --aligner::command"
            );
        }
        // Likewise `--aligner::dedup-reads`: it configures the in-process
        // backend's read-pair memo. A subprocess bwa-mem3 runs its own
        // `--dedup-reads` default, which `--aligner::command` can set.
        if self.dedup_reads.is_some() && !matches!(self.preset, Some(AlignerPreset::BwaMem3InProc))
        {
            bail!(
                "--aligner::dedup-reads is only valid with --aligner::preset \
                 bwa-mem3-inproc; it configures the in-process backend's read-pair \
                 memo (a subprocess bwa-mem3 applies its own --dedup-reads, which \
                 --aligner::command can set)"
            );
        }
        // `--aligner::max-read-ahead` bounds the subprocess backend's in-flight
        // gate; the in-process backend never blocks on aligner read-ahead.
        if self.max_read_ahead.is_some()
            && matches!(self.preset, Some(AlignerPreset::BwaMem3InProc))
        {
            bail!(
                "--aligner::max-read-ahead is not valid with --aligner::preset bwa-mem3-inproc; \
                 it bounds the input fgumi holds for a subprocess aligner"
            );
        }
        // `--aligner::input-format` changes what is written to an aligner
        // subprocess's stdin. The bwa presets align FASTQ (`mem -p`, or in
        // process) and `star` reads SAM, so only command mode may choose.
        let star = self.preset == Some(AlignerPreset::Star);
        match (self.input_format, self.preset) {
            (Some(UnalignedReadFormat::Sam), Some(preset)) if !star => bail!(
                "--aligner::input-format sam is only valid with --aligner::command; \
                 --aligner::preset {preset} reads FASTQ"
            ),
            (Some(UnalignedReadFormat::Fastq), Some(preset)) if star => bail!(
                "--aligner::input-format fastq is not valid with --aligner::preset {preset}, \
                 which reads unaligned SAM"
            ),
            _ => {}
        }
        // `--aligner::index` and `--aligner::extra-args` configure the STAR
        // command; the bwa presets index beside `--ref` and a command carries
        // its own arguments.
        if !star {
            for (set, flag) in [
                (self.index.is_some(), "--aligner::index"),
                (self.extra_args.is_some(), "--aligner::extra-args"),
                (self.startup_timeout.is_some(), "--aligner::startup-timeout"),
            ] {
                if set {
                    bail!(
                        "{flag} is only valid with --aligner::preset star (bwa presets read \
                         the index beside --ref; --aligner::command carries its own arguments)"
                    );
                }
            }
        }
        Ok(())
    }

    /// Reject an aligner the realign stage cannot use: `star`. Realign chooses
    /// each read's primary among bwa's `XA` alternatives and rescores a
    /// promoted alignment with bwa scoring, and STAR emits neither.
    ///
    /// # Errors
    ///
    /// Returns an error for `--aligner::preset star`.
    pub(crate) fn check_realign_support(&self) -> Result<()> {
        if self.preset == Some(AlignerPreset::Star) {
            bail!(
                "--aligner::preset star cannot realign: realign picks alternatives from bwa's \
                 XA tag and rescores them with bwa scoring, which STAR does not emit; use a \
                 bwa-mem3 preset"
            );
        }
        Ok(())
    }

    /// Validate the option combination and produce a [`ResolvedAligner`]
    /// ready for the align stage's `backend_for` to construct.
    ///
    /// # Arguments
    ///
    /// * `reference` - From `runall --ref`. Required in both modes.
    /// * `top_threads` - From `runall --threads`. Used as the preset
    ///   default if `--aligner::threads` is unset.
    /// * `aligner_bin` - From `runall --aligner-bin`. Preset-mode-only
    ///   binary override; must be `None` in command mode (rejected
    ///   here with a clear error).
    ///
    /// # Errors
    ///
    /// - `--aligner::chunk-size` is zero, or exceeds `i32::MAX` with a preset.
    /// - Neither `--aligner::preset` nor `--aligner::command` was set.
    /// - Both were set (mutual exclusion violation).
    /// - `--aligner::sub-batch-templates` set with any preset other than
    ///   `bwa-mem3-inproc`, or with command mode; or set to zero or more than
    ///   `MAX_SUB_BATCH_TEMPLATES`.
    /// - `--aligner::dedup-reads` set with any preset other than
    ///   `bwa-mem3-inproc`, or with command mode.
    /// - `--aligner::input-format sam` with a bwa preset, or `fastq` with `star`.
    /// - `--aligner::index` or `--aligner::extra-args` with anything but `star`.
    /// - `star` without `--aligner::index`, with an index whose
    ///   chromosome metadata disagrees with the `--ref` `.dict`, or with an
    ///   extra STAR parameter the streaming preset cannot run with (see
    ///   `star::StarExtraArgs::parse`).
    /// - `--aligner::max-read-ahead` set with `bwa-mem3-inproc`.
    /// - Command mode + a preset-only flag (`--aligner-bin` or
    ///   `--aligner::threads`).
    /// - Command mode template missing `{ref}`.
    /// - `--aligner::preset bwa-mem3-inproc` + `--aligner::threads` or
    ///   `--aligner-bin` (the in-process backend has no subprocess and
    ///   shares runall's thread budget).
    /// - `--aligner::preset bwa-mem3-inproc` built without the
    ///   `aligner-bwa-mem3` feature.
    /// - Preset-mode index files / binary missing (delegated to
    ///   [`AlignerPreset::validate_with_methylation`]).
    /// - `methylation_mode` EM-seq with the `bwa` preset, or EM-seq or TAPS
    ///   with a bwa-mem3 preset whose `bwa-mem3 index --meth` files are missing
    ///   (also delegated to [`AlignerPreset::validate_with_methylation`]).
    ///   Command mode is not checked: the command picks its own aligner and
    ///   flags.
    pub(crate) fn resolve(
        self,
        reference: &Path,
        top_threads: usize,
        aligner_bin: Option<&Path>,
    ) -> Result<ResolvedAligner> {
        // A zero chunk size reaches the aligner as `-K 0` in preset mode and
        // zeroes the pipeline's in-flight unmapped budget in both modes;
        // attribute the failure to the flag rather than to a cryptic aligner
        // error downstream.
        if self.chunk_size == 0 {
            bail!(
                "--aligner::chunk-size must be greater than 0 (got 0); it sets the \
                 aligner's -K batch size and the pipeline's in-flight unmapped budget"
            );
        }
        // Every preset hands the chunk size to bwa's `-K`, which bwa and
        // bwa-mem3 parse with `atoi` into an `int`: a larger value overflows
        // there, so the subprocess preset would fall back to thread-scaled
        // batching while the in-process cutter used the exact value — the two
        // backends would cut cohorts differently. Command mode sets its own
        // `-K` (this flag only sizes its in-flight budget), so it is exempt.
        if self.preset.is_some_and(|p| p != AlignerPreset::Star)
            && self.chunk_size > MAX_PRESET_CHUNK_SIZE
        {
            bail!(
                "--aligner::chunk-size must be at most {MAX_PRESET_CHUNK_SIZE} (got {}) with \
                 --aligner::preset: bwa parses -K as a 32-bit int",
                self.chunk_size
            );
        }
        self.check_backend_specific_flags()?;
        match (self.preset, self.command.clone()) {
            (None, None) => bail!(
                "aligning requires one of `--aligner::preset` \
                 (e.g. bwa-mem3, bwa) or `--aligner::command \"...\"`"
            ),
            (Some(_), Some(_)) => bail!(
                "--aligner::preset and --aligner::command are mutually exclusive; \
                 pass only one"
            ),
            (Some(AlignerPreset::Star), None) => {
                star::resolve_star(&self, reference, top_threads, aligner_bin)
            }
            (Some(AlignerPreset::BwaMem3InProc), None) => resolve_inproc(
                reference,
                self.threads,
                aligner_bin,
                self.sub_batch_templates,
                self.dedup_reads,
                self.chunk_size,
                self.methylation_mode,
            ),
            (Some(preset), None) => {
                // Preset mode: validate indexes + binary, build argv.
                preset.validate_with_methylation(reference, aligner_bin, self.methylation_mode)?;
                // An explicit `--aligner::threads 0` would reach the aligner as
                // `-t 0`. `None` means "default to runall's --threads" (not
                // zero), so only reject an explicit zero.
                if self.threads == Some(0) {
                    bail!("--aligner::threads must be greater than 0 (got 0)");
                }
                let threads = self.threads.unwrap_or(top_threads);
                let command = preset.build_command(
                    reference,
                    threads,
                    self.chunk_size,
                    aligner_bin,
                    self.methylation_mode,
                );
                Ok(ResolvedAligner {
                    // Both subprocess presets run `mem -p -K`, whose smart
                    // pairing splits a pair across a chunk cut (bwa and
                    // bwa-mem3 share `bseq_read`'s even-count cut and
                    // `bseq_classify`).
                    backend: ResolvedBackend::Subprocess {
                        command,
                        accept_mid_pair_split: true,
                        input_format: UnalignedReadFormat::Fastq,
                        max_read_ahead: self.max_read_ahead,
                        runtime: SubprocessRuntime::default(),
                    },
                    chunk_size: self.chunk_size,
                    threads: Some(threads),
                    mode: ResolvedAlignerMode::Preset(preset),
                })
            }
            (None, Some(template)) => {
                // Command mode: only `--aligner-bin` and
                // `--aligner::threads` are preset-only knobs; rejecting
                // them surfaces the misuse loud.
                //
                // `--aligner::chunk-size` is NOT preset-only — the chain
                // builder derives the AAM step's
                // `in_flight_unmapped_budget` from it regardless of mode
                // (see `in_flight_budget_for_chunk_size`), and
                // command-mode users may legitimately want to align it
                // with their `-K` choice in the template.
                if aligner_bin.is_some() {
                    bail!(
                        "--aligner-bin is only valid with --aligner::preset; \
                         pass the binary path inside --aligner::command \"...\" instead"
                    );
                }
                if self.threads.is_some() {
                    bail!(
                        "--aligner::threads is only valid with --aligner::preset; \
                         hardcode the thread count in --aligner::command \"...\" \
                         or use the `{{threads}}` placeholder"
                    );
                }
                // Substitute {ref} / {threads}. {ref} is required;
                // {threads} is optional in command mode (the user may
                // hardcode any thread count).
                let command = substitute_template(&template, reference, top_threads)?;
                Ok(ResolvedAligner {
                    backend: ResolvedBackend::Subprocess {
                        command,
                        accept_mid_pair_split: false,
                        input_format: self.input_format.unwrap_or_default(),
                        max_read_ahead: self.max_read_ahead,
                        runtime: SubprocessRuntime::default(),
                    },
                    chunk_size: self.chunk_size,
                    threads: None,
                    mode: ResolvedAlignerMode::Command,
                })
            }
        }
    }
}

/// Resolve `--aligner::preset bwa-mem3-inproc` into a [`ResolvedAligner`].
/// Split out of [`AlignerOptions::resolve`] (which dispatches to this
/// function for that one preset) to keep `resolve` under clippy's
/// line-count ceiling; see `resolve`'s doc comment for the full error
/// contract this preset participates in.
///
/// The flag-misuse checks (`threads`, `aligner_bin`, a zero or oversized
/// `sub_batch_templates`) run BEFORE the feature-off bail so a caller
/// passing both a bad flag and running a feature-off binary sees the
/// flag-specific error rather than the generic "rebuild fgumi" one — and
/// so these checks (and their error text) are exercised identically
/// whether or not this binary was built with `aligner-bwa-mem3`.
///
/// `reference` and `chunk_size` are read only inside the
/// `#[cfg(feature = "aligner-bwa-mem3")]` arm below, so a feature-off
/// build never references them; the `cfg_attr` below covers that build
/// instead of underscore-prefixing names that ARE used once the feature
/// is on.
#[cfg_attr(not(feature = "aligner-bwa-mem3"), allow(unused_variables))]
fn resolve_inproc(
    reference: &Path,
    threads: Option<usize>,
    aligner_bin: Option<&Path>,
    sub_batch_templates: Option<usize>,
    dedup_reads: Option<DedupReads>,
    chunk_size: u64,
    methylation: MethylationMode,
) -> Result<ResolvedAligner> {
    if threads.is_some() {
        bail!(
            "--aligner::threads is not valid with --aligner::preset \
             bwa-mem3-inproc: the in-process aligner has a single \
             budget with the rest of runall (set --threads instead)"
        );
    }
    if aligner_bin.is_some() {
        bail!(
            "--aligner-bin is not valid with --aligner::preset \
             bwa-mem3-inproc: the in-process backend links bwa-mem3 \
             directly, so there is no binary to override"
        );
    }
    if sub_batch_templates == Some(0) {
        bail!("--aligner::sub-batch-templates must be greater than 0 (got 0)");
    }
    if let Some(n) = sub_batch_templates
        && n > MAX_SUB_BATCH_TEMPLATES
    {
        bail!("--aligner::sub-batch-templates must be at most {MAX_SUB_BATCH_TEMPLATES} (got {n})");
    }

    #[cfg(not(feature = "aligner-bwa-mem3"))]
    {
        bail!(
            "--aligner::preset bwa-mem3-inproc requires fgumi built with \
             `--features aligner-bwa-mem3`; this binary was built without \
             it (use --aligner::preset bwa-mem3 for the subprocess aligner)"
        )
    }

    #[cfg(feature = "aligner-bwa-mem3")]
    {
        // Binary discovery is skipped for this preset (`requires_binary()`
        // is false); index files are still required.
        AlignerPreset::BwaMem3InProc.validate_with_methylation(reference, None, methylation)?;
        let sub_batch_templates = sub_batch_templates.unwrap_or_else(default_sub_batch_templates);
        Ok(ResolvedAligner {
            backend: ResolvedBackend::InProcessBwaMem3 {
                reference: reference.to_path_buf(),
                sub_batch_templates,
                dedup_reads: dedup_reads != Some(DedupReads::Off),
                methylation,
                xa_max_hits: None,
            },
            chunk_size,
            threads: None,
            mode: ResolvedAlignerMode::Preset(AlignerPreset::BwaMem3InProc),
        })
    }
}

/// Append a suffix to a path without replacing the existing extension.
///
/// For example, `append_extension("/ref/genome.fa", ".bwt")` →
/// `/ref/genome.fa.bwt`.
fn append_extension(path: &Path, suffix: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// Reject paths containing characters that bash interprets specially
/// inside an unquoted command argument. Used by preset-mode `validate`
/// — the argv it builds is interpolated into `/bin/bash -c "..."` so
/// any whitespace, quote, or metacharacter would break word-splitting.
///
/// The check is conservative: we allow only `[A-Za-z0-9_./:-]` plus
/// non-ASCII (UTF-8 bytes that aren't whitespace or metas). `:` is
/// permitted because it is not a shell metacharacter in argument position
/// (it only matters to bash inside parameter expansions / as the `PATH`
/// separator, neither of which applies to an interpolated argv word).
/// Users with paths outside this set can switch to `--aligner::command
/// "..."` (free-form mode), where they own the quoting.
///
/// `flag_label` is interpolated into the error message so the user
/// sees which flag they need to clean up.
fn check_shell_safe_path(path: &Path, flag_label: &str) -> Result<()> {
    let s = path.to_string_lossy();
    // A leading '-' makes the interpolated word look like a flag to the aligner
    // (e.g. `-genome.fa` parsed as an option). The per-character allowlist below
    // permits '-' (it is legal mid-path), so guard the leading position here.
    if s.starts_with('-') {
        bail!(
            "{flag_label} path {:?} starts with '-', which the aligner would parse \
             as a flag; prefix it with `./` or use `--aligner::command \"...\"` \
             (free-form mode) instead.",
            path.display().to_string()
        );
    }
    for ch in s.chars() {
        let safe = ch.is_alphanumeric()
            || matches!(ch, '_' | '.' | '/' | '-' | ':')
            || (!ch.is_ascii() && !ch.is_whitespace());
        if !safe {
            bail!(
                "{flag_label} path {:?} contains a shell-unsafe character ({ch:?}); \
                 preset-mode commands flow through `/bin/bash -c` and cannot safely \
                 handle this path. Use `--aligner::command \"...\"` (free-form mode, \
                 with your own quoting) for paths with spaces or shell metacharacters.",
                path.display().to_string()
            );
        }
    }
    Ok(())
}

/// Reject a path containing a control character (tab, newline, ...): the
/// in-process backend writes the reference path into the `@PG CL:` header
/// line, where a tab would start a new field and a newline a new record.
fn check_header_safe_path(path: &Path, flag_label: &str) -> Result<()> {
    if let Some(ch) = path.to_string_lossy().chars().find(|c| c.is_control()) {
        bail!(
            "{flag_label} path {:?} contains a control character ({ch:?}), which \
             cannot be written into the output's @PG CL: header line",
            path.display().to_string()
        );
    }
    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

/// The seven files `bwa-mem3 index --meth --meth-only` writes: the original
/// reference's `.amb/.ann/.pac` and the full `.meth` seed index, with no
/// plain `<ref>.bwt.2bit.64`.
#[cfg(test)]
pub(crate) const METH_ONLY_INDEX_FILES: [&str; 7] =
    [".amb", ".ann", ".pac", ".meth.amb", ".meth.ann", ".meth.bwt.2bit.64", ".meth.pac"];

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    use rstest::rstest;

    use super::*;

    /// A user-set purge delay (either name, any case) wins: the delay is not
    /// set.
    #[rstest]
    #[case::unset(&["PATH", "MIMALLOC_VERBOSE"], true)]
    #[case::purge_delay(&["MIMALLOC_PURGE_DELAY"], false)]
    #[case::legacy_name(&["MIMALLOC_RESET_DELAY"], false)]
    #[case::any_case(&["mimalloc_purge_delay"], false)]
    fn the_purge_delay_is_retained_unless_the_user_set_it(
        #[case] names: &[&str],
        #[case] expected: bool,
    ) {
        let mut called = false;
        assert_eq!(retain_unless_user_set(names, || called = true), expected);
        assert_eq!(called, expected);
    }

    #[rstest]
    #[case::current_name(&["MIMALLOC_PURGE_DELAY"], true)]
    #[case::legacy_name(&["MIMALLOC_RESET_DELAY"], true)]
    #[case::lowercase(&["mimalloc_purge_delay"], true)]
    #[case::mixed_case_legacy(&["Mimalloc_Reset_Delay"], true)]
    #[case::among_others(&["PATH", "HOME", "mimalloc_purge_delay"], true)]
    #[case::unset(&["PATH", "HOME"], false)]
    #[case::other_mimalloc_option(&["MIMALLOC_VERBOSE", "MIMALLOC_ARENA_EAGER_COMMIT"], false)]
    #[case::prefix_only(&["MIMALLOC_PURGE_DELAY_MS", "X_MIMALLOC_PURGE_DELAY"], false)]
    fn names_set_mimalloc_purge_matches_like_mimalloc(
        #[case] names: &[&str],
        #[case] expected: bool,
    ) {
        assert_eq!(names_set_mimalloc_purge(names), expected);
    }

    /// A sub-batch of pairs must fill exactly one bwa-mem3 kernel batch: two
    /// reads per template, `BATCH_SIZE` reads per kernel call.
    #[cfg(feature = "aligner-bwa-mem3")]
    #[test]
    fn default_sub_batch_fills_one_kernel_batch() {
        assert_eq!(2 * default_sub_batch_templates(), bwa_mem3_rs::kernel_batch_size());
        let expected = if cfg!(target_arch = "aarch64") { 512 } else { 256 };
        assert_eq!(default_sub_batch_templates(), expected);
    }

    /// A non-UTF-8 stderr line must NOT halt the relay: draining has to
    /// continue to EOF, otherwise the aligner blocks once its stderr pipe fills
    /// (~64 KiB) and the whole align pipeline deadlocks. Regression test for the
    /// `lines()`-breaks-on-invalid-UTF-8 hang.
    #[test]
    fn relay_stderr_keeps_draining_past_non_utf8_line() {
        let mut data: Vec<u8> = Vec::new();
        data.extend_from_slice(b"first line\n");
        data.extend_from_slice(&[0xff, 0xfe, b'\n']); // invalid UTF-8 line
        data.extend_from_slice(b"third line\n");

        let capture = parking_lot::Mutex::new(StderrCapture::new(8));
        relay_stderr(std::io::Cursor::new(data), &capture);
        let ring = capture.into_inner().lines;

        // All three lines were read (the invalid one lossily decoded), proving
        // the drain continued past the bad line to EOF.
        assert_eq!(ring.len(), 3, "drain must not stop at the non-UTF-8 line");
        assert_eq!(ring[0], "first line");
        // Pin the "lossily decoded" claim itself: 0xff and 0xfe are each individually
        // ill-formed in UTF-8 (never a valid lead or continuation byte), so
        // `String::from_utf8_lossy` emits one U+FFFD per byte — the bad bytes are
        // replaced, not dropped.
        assert_eq!(
            ring[1], "\u{FFFD}\u{FFFD}",
            "invalid bytes must be lossily decoded, not dropped"
        );
        assert_eq!(ring[2], "third line");
    }

    /// Spawn a simple `echo` command and verify that stdout can be read.
    #[test]
    fn test_spawn_echo() {
        let mut proc = AlignerProcess::spawn("echo hello", 10).expect("spawn should succeed");
        let mut stdout = proc.take_stdout().expect("stdout should be available");

        let mut output = String::new();
        stdout.read_to_string(&mut output).expect("should read stdout");
        assert_eq!(output.trim(), "hello");

        proc.wait().expect("process should exit successfully");
    }

    /// The aligner child inherits fgumi's bounded mimalloc purge delay, not
    /// mimalloc's 1 s default and not a never-purge `-1`. Skipped when the
    /// environment already sets the delay, which the child then inherits as is.
    #[test]
    fn test_spawn_passes_bounded_purge_delay_to_child() {
        if user_set_mimalloc_purge() {
            return;
        }
        let mut proc = AlignerProcess::spawn("printf %s \"${MIMALLOC_PURGE_DELAY-unset}\"", 10)
            .expect("spawn should succeed");
        let mut stdout = proc.take_stdout().expect("stdout should be available");

        let mut output = String::new();
        stdout.read_to_string(&mut output).expect("should read stdout");
        assert_eq!(output, fgumi_sort::RETAINED_PURGE_DELAY_MS.to_string());

        proc.wait().expect("process should exit successfully");
    }

    /// Spawn a command that writes to both stdout and stderr; verify the stderr
    /// relay does not corrupt or interleave into stdout. The stderr ring itself
    /// is pinned by `test_nonzero_exit_surfaces_stderr`, since `wait()` discards
    /// the ring on a successful exit.
    #[test]
    fn test_stderr_does_not_leak_into_stdout() {
        let mut proc = AlignerProcess::spawn("bash -c 'echo err >&2; echo out'", 10)
            .expect("spawn should succeed");

        let mut stdout = proc.take_stdout().expect("stdout should be available");
        let mut stdout_buf = String::new();
        stdout.read_to_string(&mut stdout_buf).expect("should read stdout");
        assert_eq!(stdout_buf.trim(), "out");

        proc.wait().expect("process should exit successfully");
    }

    /// Spawn a command that exits with a non-zero status and verify
    /// that `wait` returns an error containing the captured stderr.
    #[test]
    fn test_nonzero_exit_surfaces_stderr() {
        let proc = AlignerProcess::spawn("bash -c 'echo oh-no >&2; exit 1'", 10).expect("spawn ok");
        let result = proc.wait();
        assert!(result.is_err(), "wait() should return Err for non-zero exit");
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("status"), "error message should mention status: {msg}");
        assert!(msg.contains("oh-no"), "error should include stderr ring: {msg}");
    }

    /// Child exit is insufficient while the stderr relay still has work.
    #[rstest]
    #[case::relay_pending(false)]
    #[case::relay_finished(true)]
    fn completion_probe_does_not_wait_for_stderr(#[case] release_relay: bool) {
        let mut process =
            AlignerProcess::spawn("printf 'failure-detail\\n' >&2; exit 42", 10).expect("spawn");
        process.child.wait().expect("child exit");
        process.stderr_thread.take().expect("relay").join().expect("stderr capture");
        let (release_tx, release_rx) = std::sync::mpsc::channel();
        process.stderr_thread = Some(thread::spawn(move || {
            let _ = release_rx.recv_timeout(std::time::Duration::from_secs(2));
        }));
        if release_relay {
            release_tx.send(()).expect("release relay");
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(2);
            while !process.stderr_thread.as_ref().expect("relay").is_finished() {
                assert!(std::time::Instant::now() < deadline, "relay must finish");
                thread::yield_now();
            }
        }
        let start = std::time::Instant::now();
        let finished = process.is_finished().expect("nonblocking completion probe");
        assert!(start.elapsed() < std::time::Duration::from_secs(1), "probe must not join stderr");
        assert_eq!(finished, release_relay);
        drop(release_tx);
        let error = process.wait().expect_err("collect the failed child's result").to_string();
        assert!(error.contains("exited with status 42."), "{error}");
        assert!(error.contains("failure-detail"), "{error}");
    }

    /// An aligner killed by a signal is reported with the signal's number, so
    /// e.g. an OOM kill (9) is not mistaken for any other death.
    #[cfg(unix)]
    #[test]
    fn test_killed_aligner_names_the_signal() {
        let proc = AlignerProcess::spawn("kill -KILL $$", 10).expect("spawn ok");
        let msg = proc.wait().expect_err("a killed aligner fails").to_string();
        assert!(msg.contains("exited with status signal 9."), "error should name SIGKILL: {msg}");
    }

    /// Verify that stdin can be written to and the subprocess reads it.
    #[test]
    fn test_stdin_write() {
        let mut proc = AlignerProcess::spawn("bash -c 'cat'", 10).expect("spawn should succeed");
        let mut stdin = proc.take_stdin().expect("stdin should be available");
        let mut stdout = proc.take_stdout().expect("stdout should be available");

        stdin.write_all(b"hello pipe\n").expect("should write to stdin");
        drop(stdin); // close stdin so the child process can exit

        let mut output = String::new();
        stdout.read_to_string(&mut output).expect("should read stdout");
        assert_eq!(output.trim(), "hello pipe");

        proc.wait().expect("process should exit successfully");
    }

    /// `kill()` reports `true` when the child is reaped within the deadline.
    /// A plain `sleep` responds to SIGKILL immediately, so the bounded kill
    /// reaps it well inside the reap deadline — and a `true` return is what lets
    /// `SubprocessAlignStep::drop` safely issue its follow-up `wait()` without
    /// risking a teardown hang.
    #[test]
    fn test_kill_reports_reaped_for_killable_child() {
        let mut proc = AlignerProcess::spawn("sleep 30", 10).expect("spawn should succeed");
        assert!(proc.kill(), "a killable child must be reaped within the deadline");
        // A second kill on an already-gone child still reports reaped, never
        // a spurious timeout.
        assert!(proc.kill(), "killing an already-reaped child must still report reaped");
    }

    /// `pid` as a rustix [`Pid`](rustix::process::Pid). Panics on 0 or a value
    /// past `i32::MAX`, neither of which a test's child process can have.
    #[cfg(unix)]
    fn as_pid(pid: u32) -> rustix::process::Pid {
        i32::try_from(pid).ok().and_then(rustix::process::Pid::from_raw).expect("a valid pid")
    }

    /// Whether `pid` is running, excluding zombies, which hold no pipes.
    /// Keep this separate from process-group emptiness: zombies remain members.
    #[cfg(unix)]
    fn process_exists(pid: u32) -> bool {
        let probe = rustix::process::test_kill_process(as_pid(pid));
        if probe == Err(rustix::io::Errno::SRCH) {
            return false;
        }
        let state = std::process::Command::new("/bin/ps")
            .args(["-o", "stat=", "-p", &pid.to_string()])
            .output()
            .expect("probe process state");
        if !state.status.success() {
            assert_eq!(
                rustix::process::test_kill_process(as_pid(pid)),
                Err(rustix::io::Errno::SRCH),
                "ps failed without confirming PID absence: {state:?}"
            );
            return false;
        }
        let state = String::from_utf8_lossy(&state.stdout);
        assert!(!state.trim().is_empty(), "ps returned no state");
        !state.trim_start().starts_with('Z')
    }

    /// An exited but unreaped child must not count as a live pipe holder.
    #[cfg(unix)]
    #[rstest]
    #[case::running(false, true)]
    #[case::zombie(true, false)]
    fn test_process_liveness_excludes_zombies(#[case] exited: bool, #[case] expected: bool) {
        let mut child = std::process::Command::new("sleep").arg("30").spawn().expect("spawn sleep");
        if exited {
            child.kill().expect("kill child without reaping");
            let deadline = std::time::Instant::now() + Duration::from_secs(5);
            loop {
                let state = std::process::Command::new("ps")
                    .args(["-o", "stat=", "-p", &child.id().to_string()])
                    .output()
                    .expect("ps state");
                if String::from_utf8_lossy(&state.stdout).trim_start().starts_with('Z') {
                    break;
                }
                assert!(std::time::Instant::now() < deadline, "child did not become a zombie");
                thread::sleep(Duration::from_millis(10));
            }
        }
        let alive = process_exists(child.id());
        let _ = child.kill();
        child.wait().expect("reap child");
        assert_eq!(alive, expected);
    }

    /// Poll until `pid` is gone, for at most `timeout`; returns whether it went.
    #[cfg(unix)]
    fn wait_for_exit(pid: u32, timeout: Duration) -> bool {
        let deadline = std::time::Instant::now() + timeout;
        while process_exists(pid) {
            if std::time::Instant::now() >= deadline {
                return false;
            }
            thread::sleep(Duration::from_millis(20));
        }
        true
    }

    /// Retrying a failed group teardown cannot turn a cached wrapper exit
    /// into confirmation that every group member is gone.
    #[cfg(unix)]
    #[test]
    fn test_kill_retry_preserves_unconfirmed_group_result() {
        let mut proc = AlignerProcess::spawn("exit 0", 10).unwrap();
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        while !proc.has_exited() {
            assert!(std::time::Instant::now() < deadline);
            thread::sleep(Duration::from_millis(10));
        }
        proc.process_group = None;
        proc.group_teardown_confirmed = Some(false);
        assert!(!proc.kill(), "a reaped wrapper does not confirm group teardown");
    }

    /// An unconfirmed group result must not prevent reaping a wrapper that
    /// exits during Drop's bounded cleanup.
    #[cfg(unix)]
    #[test]
    fn test_drop_reaps_wrapper_after_unconfirmed_group_teardown() {
        let mut proc = AlignerProcess::spawn("sleep 0.1", 10).unwrap();
        let pid = as_pid(proc.pid());
        proc.process_group = None;
        proc.group_teardown_confirmed = Some(false);
        drop(proc);
        let probe = rustix::process::waitpid(Some(pid), rustix::process::WaitOptions::NOHANG);
        // Reap a regressed implementation's zombie, without hiding its result.
        if !matches!(probe, Err(rustix::io::Errno::CHILD)) {
            let _ = rustix::process::waitpid(Some(pid), rustix::process::WaitOptions::empty());
        }
        assert!(
            matches!(probe, Err(rustix::io::Errno::CHILD)),
            "Drop left the wrapper unreaped: {probe:?}"
        );
    }

    /// Read the PID a test aligner command wrote to `path`, polling until the
    /// shell has written it.
    #[cfg(unix)]
    fn read_pid_file(path: &Path) -> u32 {
        let deadline = std::time::Instant::now() + Duration::from_secs(10);
        loop {
            if let Ok(text) = std::fs::read_to_string(path)
                && let Ok(pid) = text.trim().parse()
            {
                return pid;
            }
            assert!(std::time::Instant::now() < deadline, "no pid written to {}", path.display());
            thread::sleep(Duration::from_millis(20));
        }
    }

    /// SIGKILLs a test's stray process on drop, so a failing assertion does not
    /// leave a long `sleep` behind.
    #[cfg(unix)]
    struct KillOnDrop(Option<u32>);

    #[cfg(unix)]
    impl KillOnDrop {
        fn new(pid: u32) -> Self {
            Self(Some(pid))
        }

        /// Stop signalling once exit is confirmed, before the PID can be reused.
        fn disarm(&mut self) {
            self.0 = None;
        }
    }

    #[cfg(unix)]
    impl Drop for KillOnDrop {
        fn drop(&mut self) {
            if let Some(pid) = self.0 {
                let _ = rustix::process::kill_process(as_pid(pid), rustix::process::Signal::KILL);
            }
        }
    }

    /// Teardown must kill the aligner's whole process group, not only the
    /// `/bin/bash -c` wrapper. The command backgrounds a long-lived grandchild
    /// that inherits the stdout and stderr pipes; with only the wrapper killed,
    /// the grandchild keeps both pipes open, so the stdout reader never sees EOF
    /// (and the stderr relay is left running past its drain grace). The teardown here
    /// is the one `SubprocessAlignStep::drop` runs (`kill()` then `wait()`),
    /// bounded by a timeout so a regression fails instead of hanging. The
    /// `ignores_sigterm` case needs the SIGKILL that follows the SIGTERM grace.
    #[cfg(unix)]
    #[rstest]
    #[case::exits_on_sigterm("-", false)]
    #[case::ignores_sigterm("''", true)]
    fn test_kill_terminates_grandchild_holding_pipes(
        #[case] term_trap: &str,
        #[case] requires_escalation: bool,
    ) {
        let dir = tempfile::tempdir().expect("tempdir");
        let pid_file = dir.path().join("grandchild.pid");
        let ready_file = dir.path().join("grandchild.ready");
        // Publish readiness from the helper only after installing the trap.
        // The wrapper publishes $!, which also works with macOS's Bash 3.2.
        let command = format!(
            "(trap {term_trap} TERM; echo ready > '{}'; exec sleep 600) & echo $! > '{}'; wait",
            ready_file.display(),
            pid_file.display()
        );
        let mut proc = AlignerProcess::spawn(&command, 10).expect("spawn should succeed");
        let mut stdout = proc.take_stdout().expect("stdout should be available");
        let grandchild_pid = read_pid_file(&pid_file);
        let mut cleanup = KillOnDrop::new(grandchild_pid);
        let ready_by = std::time::Instant::now() + Duration::from_secs(10);
        while !ready_file.exists() {
            assert!(std::time::Instant::now() < ready_by, "helper never installed its TERM trap");
            thread::sleep(Duration::from_millis(10));
        }

        // As in the pipeline, one thread blocks reading stdout while another
        // tears the aligner down.
        let reader = thread::spawn(move || {
            let mut sink = Vec::new();
            stdout.read_to_end(&mut sink).is_ok()
        });
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        thread::spawn(move || {
            let start = std::time::Instant::now();
            let group_gone = proc.kill();
            let elapsed = start.elapsed();
            let wrapper_reaped = proc.child.try_wait().expect("poll wrapper").is_some();
            if group_gone {
                let _ = proc.wait();
            }
            let stdout_closed = reader.join().unwrap_or(false);
            let _ = done_tx.send((wrapper_reaped, stdout_closed, elapsed));
        });
        let (reaped, stdout_closed, elapsed) = done_rx
            .recv_timeout(Duration::from_secs(30))
            .expect("teardown hung: the grandchild kept the aligner's pipes open");
        assert!(reaped, "kill() must reap the wrapper even if zombies remain in the group");
        if requires_escalation {
            assert!(
                elapsed >= PROCESS_GROUP_TERM_GRACE,
                "TERM-ignoring helper must require escalation: {elapsed:?}"
            );
        }
        assert!(stdout_closed, "stdout must reach EOF once the group is killed");
        let exited = wait_for_exit(grandchild_pid, Duration::from_secs(10));
        if exited {
            cleanup.disarm();
        }
        assert!(exited, "grandchild (pid {grandchild_pid}) survived teardown");
    }

    /// Dropping an `AlignerProcess` whose wrapper has exited but was never
    /// `wait()`ed still kills its process group: a grandchild left holding the
    /// stdout and stderr pipes would otherwise outlive fgumi, and keep `Drop`'s
    /// stderr-relay join waiting out its drain grace.
    #[cfg(unix)]
    #[test]
    fn test_drop_kills_group_after_wrapper_exits_without_wait() {
        let dir = tempfile::tempdir().expect("tempdir");
        let pid_file = dir.path().join("grandchild.pid");
        let command = format!("sleep 600 & echo $! > '{}'", pid_file.display());
        let mut proc = AlignerProcess::spawn(&command, 10).expect("spawn should succeed");
        let grandchild_pid = read_pid_file(&pid_file);
        let mut cleanup = KillOnDrop::new(grandchild_pid);

        let deadline = std::time::Instant::now() + Duration::from_secs(10);
        while proc.child.try_wait().expect("try_wait").is_none() {
            assert!(std::time::Instant::now() < deadline, "the wrapper never exited");
            thread::sleep(Duration::from_millis(20));
        }
        let (done_tx, done_rx) = std::sync::mpsc::channel();
        thread::spawn(move || {
            drop(proc);
            let _ = done_tx.send(());
        });
        done_rx
            .recv_timeout(Duration::from_secs(30))
            .expect("Drop hung: the grandchild kept the stderr pipe open");
        let exited = wait_for_exit(grandchild_pid, Duration::from_secs(10));
        if exited {
            cleanup.disarm();
        }
        assert!(exited, "grandchild (pid {grandchild_pid}) survived Drop");
    }

    /// A clean exit never signals the process group: by then the group's id
    /// may belong to unrelated processes. A grandchild that detached from the
    /// pipes and outlives a successful `wait()` is left running.
    #[cfg(unix)]
    #[test]
    fn test_clean_wait_does_not_signal_process_group() {
        let dir = tempfile::tempdir().expect("tempdir");
        let pid_file = dir.path().join("grandchild.pid");
        let command = format!(
            "(exec sleep 600 </dev/null >/dev/null 2>&1) & echo $! > '{}'",
            pid_file.display()
        );
        let proc = AlignerProcess::spawn(&command, 10).expect("spawn should succeed");
        let grandchild_pid = read_pid_file(&pid_file);
        let _cleanup = KillOnDrop::new(grandchild_pid);

        proc.wait().expect("the wrapper should exit successfully");
        thread::sleep(Duration::from_millis(200));
        assert!(process_exists(grandchild_pid), "a clean exit must not kill the process group");
    }

    /// A clean exit whose stderr pipe is still held open by a process the
    /// aligner started (here a backgrounded `sleep` that inherited stderr)
    /// must not hang `wait()` joining the stderr relay: it returns within the
    /// drain grace, with the lines captured so far, and leaves that process
    /// running. A failing exit reports the stderr tail captured by then
    /// (including a last line without a newline) and says it may be incomplete.
    #[cfg(unix)]
    #[rstest]
    #[case::clean_exit("echo started >&2", None)]
    #[case::failing_exit("echo started >&2; exit 3", Some("started"))]
    #[case::failing_exit_unterminated_last_line(
        "echo started >&2; printf 'fatal: bad index' >&2; exit 3",
        Some("started\nfatal: bad index\n")
    )]
    fn test_wait_is_bounded_when_a_grandchild_holds_stderr(
        #[case] wrapper_tail: &str,
        #[case] expected_tail: Option<&str>,
    ) {
        let dir = tempfile::tempdir().expect("tempdir");
        let pid_file = dir.path().join("grandchild.pid");
        let command = format!(
            "(exec sleep 600 </dev/null >/dev/null) & echo $! > '{}'; {wrapper_tail}",
            pid_file.display()
        );
        let proc = AlignerProcess::spawn(&command, 10).expect("spawn should succeed");
        let grandchild_pid = read_pid_file(&pid_file);
        let _cleanup = KillOnDrop::new(grandchild_pid);

        let (done_tx, done_rx) = std::sync::mpsc::channel();
        thread::spawn(move || {
            let _ = done_tx.send(proc.wait().map_err(|e| e.to_string()));
        });
        let result = done_rx
            .recv_timeout(STDERR_DRAIN_GRACE + Duration::from_secs(20))
            .expect("wait() hung on a stderr pipe held open by the aligner's child");
        match expected_tail {
            None => result.expect("a clean exit must succeed"),
            Some(tail) => {
                let msg = result.expect_err("a failing exit must be reported");
                assert!(msg.contains("status 3"), "got: {msg}");
                assert!(msg.contains(tail), "expected the tail {tail:?}, got: {msg}");
                assert!(msg.contains("later lines are not shown"), "got: {msg}");
            }
        }
        assert!(process_exists(grandchild_pid), "wait() must not kill the aligner's child");
    }

    /// A capture lock held past relay timeout must not delay the exit result.
    #[rstest]
    #[case::success_complete(true, true, false, None)]
    #[case::success_pending_ready(true, false, false, None)]
    #[case::success_pending_locked(true, false, true, None)]
    #[case::failure_complete(false, true, false, Some((true, false)))]
    #[case::failure_pending_ready(false, false, false, Some((true, true)))]
    #[case::failure_pending_locked(false, false, true, Some((false, true)))]
    fn test_exit_result_does_not_wait_for_capture_lock(
        #[case] success: bool,
        #[case] relay_complete: bool,
        #[case] capture_locked: bool,
        #[case] expected_error: Option<(bool, bool)>,
    ) {
        let mut proc = AlignerProcess::spawn(if success { "true" } else { "exit 7" }, 10)
            .expect("spawn wrapper");
        let status = proc.child.wait().expect("reap wrapper");
        #[cfg(unix)]
        {
            proc.process_group = None;
        }
        proc.stderr_thread.take().unwrap().join().expect("original relay");
        proc.stderr_capture.lock().lines.push_back("fatal-index".to_string());
        let capture = Arc::clone(&proc.stderr_capture);
        let (ready_tx, ready_rx) = std::sync::mpsc::channel();
        let (release_tx, release_rx) = std::sync::mpsc::channel();
        let (relay_done_tx, relay_done_rx) = std::sync::mpsc::channel();
        let relay = thread::spawn(move || {
            let _guard = capture_locked.then(|| capture.lock());
            ready_tx.send(()).unwrap();
            if !relay_complete {
                let _ = release_rx.recv();
            }
            let _ = relay_done_tx.send(());
        });
        ready_rx.recv_timeout(Duration::from_secs(5)).expect("capture fixture ready");
        if relay_complete {
            relay.join().expect("completed relay");
        } else {
            proc.stderr_thread = Some(relay);
        }
        let (tx, rx) = std::sync::mpsc::channel();
        let worker = thread::spawn(move || {
            let _ = tx.send(proc.exit_result(status));
        });
        let result = rx.recv_timeout(STDERR_DRAIN_GRACE + Duration::from_secs(2));
        // Release even a regressed blocking snapshot before asserting.
        let _ = release_tx.send(());
        let relay_done = relay_done_rx.recv_timeout(Duration::from_secs(5));
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        while !worker.is_finished() && std::time::Instant::now() < deadline {
            thread::sleep(Duration::from_millis(5));
        }
        let joined = worker.is_finished().then(|| worker.join());
        relay_done.expect("relay cleanup");
        assert!(joined.is_some_and(|result| result.is_ok()), "exit result worker did not finish");
        let result = result.expect("exit result waited for the capture lock");
        if let Some((expected_tail, expected_incomplete)) = expected_error {
            let error = result.expect_err("failed wrapper").to_string();
            assert!(error.contains("status 7"), "{error}");
            assert_eq!(error.contains("fatal-index"), expected_tail, "{error}");
            assert_eq!(error.contains("later lines are not shown"), expected_incomplete, "{error}");
        } else {
            result.expect("successful wrapper");
        }
    }

    /// The post-SIGKILL group probe must not restart a reap budget already
    /// consumed by wrapper reaping or scheduling delays.
    #[cfg(unix)]
    #[rstest]
    #[case::empty_group_after_budget_elapsed(false, Duration::ZERO, true)]
    #[case::unreaped_member_after_budget_elapsed(true, Duration::ZERO, false)]
    #[case::unreaped_member_with_short_budget(true, Duration::from_millis(100), false)]
    fn test_post_kill_reaping_shares_one_deadline(
        #[case] retain_member: bool,
        #[case] remaining: Duration,
        #[case] expected_empty: bool,
    ) {
        use std::os::unix::process::CommandExt;

        let mut proc = AlignerProcess::spawn("sleep 30", 0).expect("spawn wrapper");
        let group = proc.process_group.expect("wrapper process group");
        let mut member = retain_member.then(|| {
            let mut child = std::process::Command::new("sleep")
                .arg("30")
                .process_group(i32::try_from(proc.pid()).unwrap())
                .spawn()
                .expect("spawn group member");
            child.kill().expect("kill member without reaping it");
            child
        });
        proc.child.kill().expect("kill wrapper");
        proc.child.wait().expect("reap wrapper");

        let start = std::time::Instant::now();
        let result = proc.reap_killed_group(group, start + remaining);
        let elapsed = start.elapsed();
        if let Some(member) = &mut member {
            member.wait().expect("reap group member");
        }
        proc.process_group = None;
        assert_eq!(result, Ok(expected_empty));
        assert!(elapsed < Duration::from_millis(500), "reap budget restarted: {elapsed:?}");
        if retain_member {
            assert!(elapsed >= remaining, "unreaped group confirmed too soon: {elapsed:?}");
        }
    }

    /// A blocked logging sink must not extend bounded failure cleanup.
    #[cfg(unix)]
    #[rstest]
    #[case::stderr_drain(false)]
    #[case::unconfirmed_group(true)]
    fn test_stderr_timeout_does_not_write_to_logging_sink(#[case] zombie_group: bool) {
        use std::os::unix::process::CommandExt;

        struct BlockingLogger {
            released: std::sync::Mutex<bool>,
            ready: std::sync::Condvar,
        }
        impl log::Log for BlockingLogger {
            fn enabled(&self, _: &log::Metadata<'_>) -> bool {
                true
            }
            fn log(&self, _: &log::Record<'_>) {
                let mut released = self.released.lock().unwrap();
                while !*released {
                    released = self.ready.wait(released).unwrap();
                }
            }
            fn flush(&self) {}
        }
        static LOGGER: BlockingLogger = BlockingLogger {
            released: std::sync::Mutex::new(false),
            ready: std::sync::Condvar::new(),
        };
        // libtest runs cases in one process; isolate the process-global logger
        // so neither another case nor an unrelated test can install or release it.
        let current_thread = thread::current();
        let test_name = current_thread.name().expect("named test thread");
        if std::env::var("FGUMI_BLOCKING_LOGGER_CASE").as_deref() != Ok(test_name) {
            let output = std::process::Command::new(std::env::current_exe().expect("test binary"))
                .args(["--exact", test_name, "--nocapture"])
                .env("FGUMI_BLOCKING_LOGGER_CASE", test_name)
                .output()
                .expect("run logger case in its own process");
            assert!(output.status.success(), "logger case failed: {output:?}");
            return;
        }

        log::set_logger(&LOGGER).expect("install logger in isolated test process");
        log::set_max_level(log::LevelFilter::Warn);
        let dir = tempfile::tempdir().expect("tempdir");
        let pid_file = dir.path().join("holder.pid");
        let release = dir.path().join("release");
        let command = if zombie_group {
            format!("while [ ! -f '{}' ]; do sleep 0.01; done; exit 3", release.display())
        } else {
            format!("sleep 30 </dev/null >/dev/null & echo $! > '{}'; exit 3", pid_file.display())
        };
        let mut proc = AlignerProcess::spawn(&command, 10).expect("spawn aligner");
        let mut member = zombie_group.then(|| {
            let mut child = std::process::Command::new("sleep")
                .arg("30")
                .process_group(i32::try_from(proc.pid()).unwrap())
                .spawn()
                .unwrap();
            child.kill().unwrap();
            std::fs::write(&release, b"ready").unwrap();
            let deadline = std::time::Instant::now() + Duration::from_secs(5);
            while !proc.has_exited() {
                assert!(std::time::Instant::now() < deadline);
                thread::sleep(Duration::from_millis(10));
            }
            child
        });
        let _cleanup = (!zombie_group).then(|| KillOnDrop::new(read_pid_file(&pid_file)));
        let (tx, rx) = std::sync::mpsc::channel();
        let worker = thread::spawn(move || {
            let result = if zombie_group { proc.finish_teardown() } else { proc.wait() };
            let _ = tx.send(result);
        });
        let result = rx.recv_timeout(STDERR_DRAIN_GRACE + Duration::from_secs(2));
        *LOGGER.released.lock().unwrap() = true;
        LOGGER.ready.notify_all();
        worker.join().expect("wait thread");
        if let Some(member) = &mut member {
            member.wait().unwrap();
        }
        let error = result
            .expect("stderr timeout blocked on the logging sink")
            .expect_err("wrapper failure")
            .to_string();
        assert!(error.contains("status 3"), "{error}");
        if !zombie_group {
            assert!(error.contains("later lines are not shown"), "{error}");
        }
    }

    /// How an aligner is torn down, for
    /// [`test_teardown_is_bounded_when_a_process_outside_the_group_holds_stderr`].
    #[derive(Debug, Clone, Copy)]
    enum Teardown {
        /// The wrapper exits on its own and is dropped without `wait()`.
        DropAfterExit,
        /// The running wrapper is killed, then waited for.
        KillThenWait,
    }

    /// `Drop` and `kill()` + `wait()` join the stderr relay for at most the drain
    /// grace too. A process the aligner started in its own process group (here a
    /// `set -m` background job) survives the group teardown and keeps the
    /// inherited stderr open, so the relay never reaches EOF.
    #[cfg(unix)]
    #[rstest]
    #[case::drop_after_exit(Teardown::DropAfterExit)]
    #[case::kill_then_wait(Teardown::KillThenWait)]
    fn test_teardown_is_bounded_when_a_process_outside_the_group_holds_stderr(
        #[case] teardown: Teardown,
    ) {
        let dir = tempfile::tempdir().expect("tempdir");
        let pid_file = dir.path().join("holder.pid");
        let wrapper_tail = match teardown {
            Teardown::DropAfterExit => "",
            Teardown::KillThenWait => "; set +m; exec sleep 600",
        };
        let command = format!(
            "set -m; (exec sleep 600 </dev/null >/dev/null) & echo $! > '{}'{wrapper_tail}",
            pid_file.display()
        );
        let mut proc = AlignerProcess::spawn(&command, 10).expect("spawn should succeed");
        let holder_pid = read_pid_file(&pid_file);
        let _cleanup = KillOnDrop::new(holder_pid);
        if matches!(teardown, Teardown::DropAfterExit) {
            let deadline = std::time::Instant::now() + Duration::from_secs(10);
            while !proc.has_exited() {
                assert!(std::time::Instant::now() < deadline, "the wrapper did not exit");
                thread::sleep(Duration::from_millis(10));
            }
        }

        let (done_tx, done_rx) = std::sync::mpsc::channel();
        thread::spawn(move || {
            match teardown {
                Teardown::DropAfterExit => drop(proc),
                Teardown::KillThenWait => {
                    proc.kill();
                    let _ = proc.wait();
                }
            }
            let _ = done_tx.send(());
        });
        done_rx
            .recv_timeout(STDERR_DRAIN_GRACE + Duration::from_secs(20))
            .expect("teardown hung on a stderr pipe held open outside the aligner's group");
        assert!(process_exists(holder_pid), "a process outside the group is not signalled");
    }

    /// `is_already_reaped` confirms exit only for `ECHILD`. An unrelated errno
    /// (e.g. `EINVAL`) must NOT be treated as a confirmed reap, so `kill()`
    /// reports `false` and the caller avoids a follow-up unbounded `wait()`.
    #[rstest]
    #[case::already_reaped(Some(ECHILD), true)]
    #[case::unrelated_errno(Some(22), false)]
    #[case::non_os_error(None, false)]
    fn test_is_already_reaped_only_for_echild(#[case] errno: Option<i32>, #[case] expected: bool) {
        let error = errno.map_or_else(
            || std::io::Error::other("synthetic non-os error"),
            std::io::Error::from_raw_os_error,
        );
        assert_eq!(is_already_reaped(&error), expected);
    }

    /// An external reap removes the safe process-group identity. Cleanup must
    /// tolerate ECHILD without signalling a potentially reused group.
    #[cfg(unix)]
    #[rstest]
    #[case::wait("wait")]
    #[case::kill_wrapper("kill_wrapper")]
    #[case::drop("drop")]
    fn test_cleanup_after_external_wrapper_reap(#[case] operation: &str) {
        let mut proc = AlignerProcess::spawn("exit 0", 10).unwrap();
        rustix::process::waitpid(Some(as_pid(proc.pid())), rustix::process::WaitOptions::empty())
            .expect("external reap");
        match operation {
            "wait" => assert!(proc.wait().unwrap_err().to_string().contains("failed to wait")),
            "kill_wrapper" => {
                proc.process_group = None;
                assert!(proc.kill());
            }
            "drop" => {
                proc.process_group = None;
                drop(proc);
            }
            _ => unreachable!(),
        }
    }

    /// ECHILD prevents safe group signalling, but cannot confirm the group's
    /// surviving members have exited.
    #[cfg(unix)]
    #[test]
    fn test_external_wrapper_reap_does_not_confirm_group_exit() {
        use std::os::unix::process::CommandExt;
        let dir = tempfile::tempdir().unwrap();
        let release = dir.path().join("release");
        let command = format!("while [ ! -f '{}' ]; do sleep 0.01; done", release.display());
        let mut proc = AlignerProcess::spawn(&command, 10).unwrap();
        let mut member = std::process::Command::new("sleep")
            .arg("30")
            .process_group(i32::try_from(proc.pid()).unwrap())
            .spawn()
            .unwrap();
        std::fs::write(&release, b"ready").unwrap();
        rustix::process::waitpid(Some(as_pid(proc.pid())), rustix::process::WaitOptions::empty())
            .unwrap();
        let confirmed = proc.kill();
        let survived = process_exists(member.id());
        let _ = member.kill();
        member.wait().unwrap();
        assert!(!confirmed, "an external wrapper reap cannot confirm group exit");
        assert!(survived, "an externally reaped wrapper's group must not be signalled");
    }

    /// A relay panic must report an incomplete snapshot without logging on the
    /// caller thread or implying that the stderr pipe necessarily stayed open.
    #[test]
    fn test_relay_panic_marks_stderr_snapshot_incomplete() {
        let mut proc = AlignerProcess::spawn("exit 3", 10).unwrap();
        proc.stderr_thread = Some(thread::spawn(|| panic!("intentional relay panic")));
        let error = proc.wait().unwrap_err().to_string();
        assert!(error.contains("status 3") && error.contains("capture was incomplete"), "{error}");
    }

    /// A wrapper exit observed by kill's own initial probe must retain its
    /// diagnostic even if the caller's earlier running probe returned false.
    #[cfg(unix)]
    #[rstest]
    #[case::natural_failure("echo fatal-probe >&2; exit 7", Some("status 7"))]
    #[case::natural_success("exit 0", None)]
    fn test_teardown_observes_exit_after_prior_running_probe(
        #[case] tail: &str,
        #[case] expected: Option<&str>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let release = dir.path().join("release");
        let command =
            format!("while [ ! -f '{}' ]; do sleep 0.01; done; {tail}", release.display());
        let mut proc = AlignerProcess::spawn(&command, 10).unwrap();
        assert!(!proc.has_exited());
        std::fs::write(&release, b"ready").unwrap();
        // Observe exit without reaping or using has_exited again: the next
        // wrapper probe belongs to kill itself, as in the reviewed race.
        let deadline = std::time::Instant::now() + Duration::from_secs(5);
        loop {
            let state = std::process::Command::new("/bin/ps")
                .args(["-o", "stat=", "-p", &proc.pid().to_string()])
                .output()
                .unwrap();
            if String::from_utf8_lossy(&state.stdout).trim_start().starts_with('Z') {
                break;
            }
            assert!(std::time::Instant::now() < deadline);
            thread::sleep(Duration::from_millis(10));
        }
        let result = proc.finish_teardown();
        match expected {
            Some(status) => {
                let error = result.unwrap_err().to_string();
                assert!(error.contains(status) && error.contains("fatal-probe"), "{error}");
            }
            None => result.unwrap(),
        }
    }

    /// Preset command strings with the default-PATH binary (bare binary name,
    /// resolved by the shell via `PATH`).
    #[rstest]
    #[case::bwa(
        AlignerPreset::Bwa,
        8,
        10_000_000,
        "bwa mem -p -K 10000000 -t 8 /ref/genome.fa /dev/stdin"
    )]
    #[case::bwa_mem3(
        AlignerPreset::BwaMem3,
        4,
        5_000_000,
        "bwa-mem3 mem --bam=0 -p -K 5000000 -t 4 /ref/genome.fa /dev/stdin"
    )]
    fn preset_build_command_default_binary(
        #[case] preset: AlignerPreset,
        #[case] threads: usize,
        #[case] chunk_size: u64,
        #[case] expected: &str,
    ) {
        let reference = Path::new("/ref/genome.fa");
        let cmd =
            preset.build_command(reference, threads, chunk_size, None, MethylationMode::Disabled);
        assert_eq!(cmd, expected);
    }

    /// `--aligner-bin` override replaces the bare binary name.
    #[test]
    fn test_build_command_with_binary_override() {
        let reference = Path::new("/data/hg38.fa");
        let override_path = PathBuf::from("/opt/bwa-mem3-2.2.1/bwa-mem3");
        let cmd = AlignerPreset::BwaMem3.build_command(
            reference,
            16,
            100_000,
            Some(&override_path),
            MethylationMode::Disabled,
        );
        assert!(cmd.contains("/opt/bwa-mem3-2.2.1/bwa-mem3"), "override binary not in cmd: {cmd}");
        assert!(!cmd.starts_with("bwa-mem3 "), "bare name should not appear: {cmd}");
    }

    /// `index_extensions` shape per preset.
    #[test]
    fn test_index_extensions() {
        assert_eq!(AlignerPreset::Bwa.index_extensions(), &[".amb", ".ann", ".bwt", ".pac", ".sa"]);
        // `bwa-mem3 0.4.0+` no longer writes `.0123` (fg-labs/bwa-mem3#177), so it
        // is not a required index file.
        assert_eq!(
            AlignerPreset::BwaMem3.index_extensions(),
            &[".amb", ".ann", ".bwt.2bit.64", ".pac"]
        );
    }

    /// `AlignerPreset` CLI parsing — both kebab-case variants
    /// (`bwa`, `bwa-mem3`) round-trip cleanly through `clap::ValueEnum`.
    #[test]
    fn test_aligner_preset_parses() {
        use clap::ValueEnum;
        assert_eq!(AlignerPreset::from_str("bwa", true).unwrap(), AlignerPreset::Bwa);
        assert_eq!(AlignerPreset::from_str("BWA", true).unwrap(), AlignerPreset::Bwa);
        assert_eq!(AlignerPreset::from_str("bwa-mem3", true).unwrap(), AlignerPreset::BwaMem3);
        assert_eq!(AlignerPreset::from_str("BWA-MEM3", true).unwrap(), AlignerPreset::BwaMem3);
        // `bwa-mem` would be the kebab of the old `BwaMem` variant —
        // ensure the variant rename ($BwaMem → Bwa) cleared this
        // stale CLI value.
        assert!(AlignerPreset::from_str("bwa-mem", true).is_err());
    }

    /// `Display` round-trips through the CLI value — what the user
    /// types matches what error messages print.
    #[test]
    fn test_aligner_preset_display_roundtrips() {
        use clap::ValueEnum;
        for preset in [AlignerPreset::Bwa, AlignerPreset::BwaMem3] {
            let printed = preset.to_string();
            let parsed = AlignerPreset::from_str(&printed, false).unwrap();
            assert_eq!(parsed, preset, "Display value {printed:?} should re-parse to {preset:?}");
        }
    }

    /// `AlignerOptions::default()` must produce field values that
    /// match each `default_value_t` on the clap derive. The macro
    /// emits `default_value_t = AlignerOptions::default().<field>`
    /// for the `Multi`-side, so any divergence here would surface as
    /// `--aligner::<flag> [default: 0]` in `--help` (and worse,
    /// `-K 0` reaching the aligner). Pin the contract.
    #[test]
    fn test_aligner_options_default_pins_clap_defaults() {
        let d = AlignerOptions::default();
        assert!(d.preset.is_none());
        assert!(d.command.is_none());
        assert!(d.threads.is_none());
        assert_eq!(
            d.chunk_size, DEFAULT_ALIGNER_CHUNK_SIZE,
            "chunk_size must default to DEFAULT_ALIGNER_CHUNK_SIZE, NOT to u64::default() (0). \
             A `#[derive(Default)]` here would zero this and propagate -K 0 to the aligner."
        );
    }

    /// `resolve` returns the right `ResolvedAligner` for command mode.
    #[test]
    fn test_resolve_command_mode_basic() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            preset: None,
            command: Some("bwa-mem3 mem -p -K 1000 -t {threads} {ref} /dev/stdin".to_string()),
            threads: None,
            chunk_size: DEFAULT_ALIGNER_CHUNK_SIZE,
            sub_batch_templates: None,
            dedup_reads: None,
            input_format: None,
            index: None,
            extra_args: None,
            startup_timeout: None,
            max_read_ahead: None,
            methylation_mode: MethylationMode::Disabled,
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        let ResolvedBackend::Subprocess { command, .. } = resolved.backend else {
            panic!("command mode must resolve to ResolvedBackend::Subprocess");
        };
        assert!(command.contains(&ref_path.display().to_string()));
        assert!(command.contains("-t 4"));
        assert!(matches!(resolved.mode, ResolvedAlignerMode::Command));
    }

    /// `--aligner::input-format` reaches the subprocess backend in command mode,
    /// defaulting to FASTQ.
    #[rstest]
    #[case::unset(None, UnalignedReadFormat::Fastq)]
    #[case::fastq(Some(UnalignedReadFormat::Fastq), UnalignedReadFormat::Fastq)]
    #[case::sam(Some(UnalignedReadFormat::Sam), UnalignedReadFormat::Sam)]
    fn test_resolve_command_mode_carries_input_format(
        #[case] input_format: Option<UnalignedReadFormat>,
        #[case] expected: UnalignedReadFormat,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            command: Some("aligner {ref}".to_string()),
            input_format,
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        let ResolvedBackend::Subprocess { input_format, .. } = resolved.backend else {
            panic!("command mode must resolve to ResolvedBackend::Subprocess");
        };
        assert_eq!(input_format, expected);
    }

    /// `--aligner::input-format sam` is rejected with every preset (they read
    /// FASTQ), before the preset's index or binary is looked for.
    #[rstest]
    #[case::bwa(AlignerPreset::Bwa)]
    #[case::bwa_mem3(AlignerPreset::BwaMem3)]
    #[case::bwa_mem3_inproc(AlignerPreset::BwaMem3InProc)]
    fn test_resolve_rejects_sam_input_format_with_preset(#[case] preset: AlignerPreset) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            preset: Some(preset),
            input_format: Some(UnalignedReadFormat::Sam),
            ..AlignerOptions::default()
        };
        let msg = opts.resolve(&ref_path, 4, None).unwrap_err().to_string();
        assert!(
            msg.contains("--aligner::input-format sam is only valid with --aligner::command")
                && msg.contains(&format!("--aligner::preset {preset}")),
            "got: {msg}"
        );
    }

    /// An explicit `--aligner::input-format fastq` is accepted with every bwa
    /// subprocess preset: it names the format those presets already read.
    #[rstest]
    #[case::bwa(AlignerPreset::Bwa)]
    #[case::bwa_mem3(AlignerPreset::BwaMem3)]
    fn test_resolve_accepts_fastq_input_format_with_preset(#[case] preset: AlignerPreset) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in preset.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let bin = make_existing_binary(tmp.path());
        let opts = AlignerOptions {
            preset: Some(preset),
            input_format: Some(UnalignedReadFormat::Fastq),
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, Some(&bin)).unwrap();
        let ResolvedBackend::Subprocess { input_format, .. } = resolved.backend else {
            panic!("expected the subprocess backend");
        };
        assert_eq!(input_format, UnalignedReadFormat::Fastq);
    }

    /// `resolve` rejects `--aligner-bin` with command mode.
    #[test]
    fn test_resolve_command_mode_rejects_aligner_bin() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let bin = tmp.path().join("mock-aligner");
        std::fs::write(&bin, b"#!/bin/sh\n").unwrap();
        let opts = AlignerOptions {
            preset: None,
            command: Some("bwa-mem3 mem {ref} /dev/stdin".to_string()),
            threads: None,
            chunk_size: DEFAULT_ALIGNER_CHUNK_SIZE,
            sub_batch_templates: None,
            dedup_reads: None,
            input_format: None,
            index: None,
            extra_args: None,
            startup_timeout: None,
            max_read_ahead: None,
            methylation_mode: MethylationMode::Disabled,
        };
        let err = opts.resolve(&ref_path, 4, Some(&bin)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("--aligner-bin is only valid with --aligner::preset"), "got: {msg}");
    }

    /// `resolve` rejects neither preset nor command set.
    #[test]
    fn test_resolve_rejects_no_mode_selected() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions::default();
        let err = opts.resolve(&ref_path, 4, None).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("requires one of"), "got: {msg}");
    }

    /// Test-only helper for the `bwa-mem3-inproc` resolve-error table below.
    /// Parses a preset name plus a flat list of `--flag value` pairs into an
    /// `AlignerOptions`, then resolves it against a freshly created (empty,
    /// unindexed) reference FASTA. This is deliberately NOT a general CLI
    /// parser — it only understands the handful of flags the in-process
    /// resolve-error cases exercise (`--aligner::threads`, `--aligner-bin`,
    /// `--aligner::sub-batch-templates`, `--aligner::chunk-size`); anything
    /// else panics loudly so a typo in a case's `extra` list fails fast
    /// instead of silently no-op'ing.
    fn resolve_from_args(preset: &str, extra: &[&str]) -> Result<ResolvedAligner> {
        use clap::ValueEnum;

        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();

        let mut opts = AlignerOptions {
            preset: Some(
                AlignerPreset::from_str(preset, false)
                    .unwrap_or_else(|e| panic!("resolve_from_args: bad preset {preset:?}: {e}")),
            ),
            ..AlignerOptions::default()
        };
        let mut aligner_bin: Option<PathBuf> = None;

        let mut iter = extra.iter().copied();
        while let Some(flag) = iter.next() {
            let value = iter
                .next()
                .unwrap_or_else(|| panic!("resolve_from_args: flag {flag} is missing a value"));
            match flag {
                "--aligner::threads" => {
                    opts.threads = Some(
                        value
                            .parse()
                            .unwrap_or_else(|e| panic!("resolve_from_args: bad threads: {e}")),
                    );
                }
                "--aligner-bin" => aligner_bin = Some(PathBuf::from(value)),
                "--aligner::chunk-size" => {
                    opts.chunk_size = value
                        .parse()
                        .unwrap_or_else(|e| panic!("resolve_from_args: bad chunk-size: {e}"));
                }
                "--aligner::sub-batch-templates" => {
                    opts.sub_batch_templates = Some(value.parse().unwrap_or_else(|e| {
                        panic!("resolve_from_args: bad sub-batch-templates: {e}")
                    }));
                }
                "--aligner::dedup-reads" => {
                    use clap::ValueEnum;
                    opts.dedup_reads = Some(
                        DedupReads::from_str(value, false)
                            .unwrap_or_else(|e| panic!("resolve_from_args: bad dedup-reads: {e}")),
                    );
                }
                "--aligner::max-read-ahead" => {
                    opts.max_read_ahead =
                        Some(parse_max_read_ahead(value).unwrap_or_else(|e| {
                            panic!("resolve_from_args: bad max-read-ahead: {e}")
                        }));
                }
                other => panic!("resolve_from_args: unrecognized flag {other}"),
            }
        }

        opts.resolve(&ref_path, 4, aligner_bin.as_deref())
    }

    /// `resolve` for `--aligner::preset bwa-mem3-inproc`: the
    /// preset-specific flag rejections (`--aligner::threads`,
    /// `--aligner-bin`, a zero or oversized `--aligner::sub-batch-templates`) fire with
    /// their own needle regardless of whether this binary was built with
    /// `aligner-bwa-mem3` — see the ordering rationale in `resolve` itself.
    ///
    /// The feature-OFF-only "names the feature" case lives in its own
    /// `#[cfg(not(feature = "aligner-bwa-mem3"))]` test
    /// (`resolve_inproc_feature_off_names_feature` below), NOT in this
    /// table: under the feature, `resolve_inproc` skips the feature bail
    /// entirely and falls through to `validate`'s index-missing error
    /// instead, so asserting the feature-off needle here would fail under
    /// `cargo nextest run --features aligner-bwa-mem3` (confirmed by
    /// running it).
    #[rstest]
    #[case::inproc_rejects_threads(
        "bwa-mem3-inproc",
        &["--aligner::threads", "8"],
        "single budget"
    )]
    #[case::inproc_rejects_bin(
        "bwa-mem3-inproc",
        &["--aligner-bin", "/x"],
        "no binary to override"
    )]
    #[case::inproc_rejects_zero_sub_batch_templates(
        "bwa-mem3-inproc",
        &["--aligner::sub-batch-templates", "0"],
        "--aligner::sub-batch-templates must be greater than 0"
    )]
    #[case::inproc_rejects_sub_batch_templates_past_cap(
        "bwa-mem3-inproc",
        &["--aligner::sub-batch-templates", "65537"],
        "--aligner::sub-batch-templates must be at most 65536 (got 65537)"
    )]
    #[case::dedup_reads_rejected_with_subprocess_preset(
        "bwa-mem3",
        &["--aligner::dedup-reads", "off"],
        "--aligner::dedup-reads is only valid with --aligner::preset bwa-mem3-inproc"
    )]
    #[case::sub_batch_templates_rejected_with_subprocess_preset(
        "bwa-mem3",
        &["--aligner::sub-batch-templates", "10"],
        "only valid with --aligner::preset bwa-mem3-inproc"
    )]
    #[case::inproc_rejects_chunk_size_past_i32(
        "bwa-mem3-inproc",
        &["--aligner::chunk-size", "2147483648"],
        "--aligner::chunk-size must be at most 2147483647 (got 2147483648)"
    )]
    #[case::bwa_mem3_rejects_chunk_size_past_i32(
        "bwa-mem3",
        &["--aligner::chunk-size", "2147483648"],
        "--aligner::chunk-size must be at most 2147483647 (got 2147483648)"
    )]
    #[case::bwa_rejects_chunk_size_past_i32(
        "bwa",
        &["--aligner::chunk-size", "4294967296"],
        "--aligner::chunk-size must be at most 2147483647 (got 4294967296)"
    )]
    #[case::inproc_rejects_max_read_ahead(
        "bwa-mem3-inproc",
        &["--aligner::max-read-ahead", "8GiB"],
        "--aligner::max-read-ahead is not valid with --aligner::preset bwa-mem3-inproc"
    )]
    fn resolve_inproc_errors(#[case] preset: &str, #[case] extra: &[&str], #[case] needle: &str) {
        let err = resolve_from_args(preset, extra).unwrap_err().to_string();
        assert!(err.contains(needle), "got: {err}");
    }

    /// `resolve` for `--aligner::preset bwa-mem3-inproc` under a
    /// feature-OFF build: names the exact rebuild fix. Gated
    /// `#[cfg(not(feature = "aligner-bwa-mem3"))]` because under the
    /// feature, `resolve_inproc` never reaches this bail — it falls
    /// through to `validate` instead (see
    /// `resolve_inproc_skips_binary_discovery` below for the feature-on
    /// counterpart). Split out of `resolve_inproc_errors` per code review:
    /// an earlier, ungated version of this case passed only by accident of
    /// which feature set `cargo ci-test` happens to build with, and would
    /// fail under the `aligner-ffi` CI job (which builds with the feature).
    #[test]
    #[cfg(not(feature = "aligner-bwa-mem3"))]
    fn resolve_inproc_feature_off_names_feature() {
        let err = resolve_from_args("bwa-mem3-inproc", &[]).unwrap_err().to_string();
        assert!(
            err.contains("requires fgumi built with `--features aligner-bwa-mem3`"),
            "got: {err}"
        );
    }

    /// Under `aligner-bwa-mem3`, `bwa-mem3-inproc` skips binary discovery
    /// entirely (`AlignerPreset::requires_binary()` is `false` for it): with
    /// the reference's index files present, `resolve` must succeed even
    /// though there is no `bwa-mem3-inproc` executable anywhere on `PATH`
    /// to find — proving `validate` never calls `which::which` for this
    /// preset. If `requires_binary()` regressed to `true`, this would fail
    /// with "aligner binary `bwa-mem3-inproc` not found on PATH", not with
    /// a compile error, so this test is the only thing that would catch
    /// that regression.
    #[test]
    #[cfg(feature = "aligner-bwa-mem3")]
    fn resolve_inproc_skips_binary_discovery() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in AlignerPreset::BwaMem3InProc.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let opts =
            AlignerOptions { preset: Some(AlignerPreset::BwaMem3InProc), ..Default::default() };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        assert!(
            matches!(resolved.backend, ResolvedBackend::InProcessBwaMem3 { .. }),
            "got: {:?}",
            resolved.backend
        );
    }

    /// `--aligner::max-read-ahead` reaches the resolved subprocess aligner, as
    /// the in-flight gate's ceiling; unset leaves the default to the backend.
    #[rstest]
    #[case::unset(None)]
    #[case::set(Some(8 << 30))]
    fn resolve_command_mode_carries_max_read_ahead(#[case] max_read_ahead: Option<u64>) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            command: Some("bwa-mem3 mem -t {threads} {ref} /dev/stdin".to_string()),
            max_read_ahead,
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        let ResolvedBackend::Subprocess { max_read_ahead: resolved_max, .. } = resolved.backend
        else {
            panic!("command mode must resolve to ResolvedBackend::Subprocess");
        };
        assert_eq!(resolved_max, max_read_ahead);
    }

    /// `--aligner::max-read-ahead` takes a memory size as `--max-memory` does: a
    /// plain number is MiB, and zero is rejected.
    #[rstest]
    #[case::binary_units("8GiB", Some(8 << 30))]
    #[case::plain_mebibytes("768", Some(768 << 20))]
    #[case::zero("0", None)]
    #[case::not_a_size("lots", None)]
    fn parse_max_read_ahead_sizes(#[case] value: &str, #[case] expected: Option<u64>) {
        assert_eq!(parse_max_read_ahead(value).ok(), expected);
    }

    /// `--aligner::dedup-reads` with command mode is rejected like the other
    /// in-process-only knobs.
    #[test]
    fn resolve_command_mode_rejects_dedup_reads() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            command: Some("bwa-mem3 mem -t {threads} {ref} /dev/stdin".to_string()),
            dedup_reads: Some(DedupReads::Off),
            ..AlignerOptions::default()
        };
        let err = opts.resolve(&ref_path, 4, None).unwrap_err().to_string();
        assert!(
            err.contains(
                "--aligner::dedup-reads is only valid with --aligner::preset bwa-mem3-inproc"
            ),
            "got: {err}"
        );
    }

    /// Under `aligner-bwa-mem3`, `--aligner::dedup-reads` resolves to on unless
    /// set to `off`.
    #[rstest]
    #[case::default_is_on(&[], true)]
    #[case::on(&["--aligner::dedup-reads", "on"], true)]
    #[case::off(&["--aligner::dedup-reads", "off"], false)]
    #[cfg(feature = "aligner-bwa-mem3")]
    fn resolve_inproc_dedup_reads(#[case] extra: &[&str], #[case] expected: bool) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in AlignerPreset::BwaMem3InProc.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let mut opts =
            AlignerOptions { preset: Some(AlignerPreset::BwaMem3InProc), ..Default::default() };
        if let [_, value] = extra {
            use clap::ValueEnum;
            opts.dedup_reads = Some(DedupReads::from_str(value, false).unwrap());
        }
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        let ResolvedBackend::InProcessBwaMem3 { dedup_reads, .. } = resolved.backend else {
            panic!("bwa-mem3-inproc must resolve to ResolvedBackend::InProcessBwaMem3");
        };
        assert_eq!(dedup_reads, expected);
    }

    /// Helper: create a known-existing binary path for tests that
    /// want validate's binary check to pass. Uses a file inside the
    /// tempdir so the test does not depend on `/bin/bash` (or any
    /// system path).
    fn make_existing_binary(dir: &std::path::Path) -> PathBuf {
        let p = dir.join("mock-aligner");
        std::fs::write(&p, b"#!/bin/sh\n").unwrap();
        // `validate` now verifies executability as part of the binary check,
        // so mark the mock +x. Tests using this helper want the binary check
        // to pass and to exercise the *next* check (e.g. index files).
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = std::fs::metadata(&p).unwrap().permissions();
            perms.set_mode(0o755);
            std::fs::set_permissions(&p, perms).unwrap();
        }
        p
    }

    /// `validate` errors out on missing index files (binary + FASTA
    /// checks are arranged to pass so we exercise only the index check).
    #[test]
    fn test_validate_missing_indexes() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let bin = make_existing_binary(tmp.path());
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(&bin)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("index file not found"), "expected index error, got: {msg}");
        assert!(msg.contains("bwa-mem3 index"), "expected fix-it hint, got: {msg}");
    }

    /// A preset that passes validation resolves to the subprocess backend with
    /// bwa's mid-pair split accepted: both presets run `mem -p -K`.
    #[rstest]
    #[case::bwa_mem3(AlignerPreset::BwaMem3)]
    #[case::bwa(AlignerPreset::Bwa)]
    fn resolve_preset_accepts_mid_pair_split(
        #[case] preset: AlignerPreset,
        #[values(None, Some(8 << 30))] max_read_ahead: Option<u64>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in preset.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let bin = make_existing_binary(tmp.path());
        let opts =
            AlignerOptions { preset: Some(preset), max_read_ahead, ..AlignerOptions::default() };
        let resolved = opts.resolve(&ref_path, 4, Some(&bin)).unwrap();
        let ResolvedBackend::Subprocess {
            command,
            accept_mid_pair_split,
            input_format,
            max_read_ahead: resolved_max,
            ..
        } = resolved.backend
        else {
            panic!("a subprocess preset must resolve to ResolvedBackend::Subprocess");
        };
        assert!(accept_mid_pair_split, "preset mode must accept bwa's mid-pair split");
        assert_eq!(input_format, UnalignedReadFormat::Fastq, "presets read FASTQ");
        assert_eq!(resolved_max, max_read_ahead, "--aligner::max-read-ahead reaches the preset");
        assert!(command.starts_with(&bin.display().to_string()), "got: {command}");
    }

    /// `with_xa_max_hits` puts `-h` right after the `mem` subcommand, ahead of the
    /// positional arguments, or appends it when the command has no `mem` word.
    #[rstest]
    #[case::preset(
        "bwa-mem3 mem --bam=0 -p -K 10 -t 4 ref.fa /dev/stdin",
        "bwa-mem3 mem -h 5,200 --bam=0 -p -K 10 -t 4 ref.fa /dev/stdin"
    )]
    #[case::mem_last("bwa mem", "bwa mem -h 5,200")]
    #[case::memory_path_not_mem(
        "/opt/memx/bwa mem ref.fa /dev/stdin",
        "/opt/memx/bwa mem -h 5,200 ref.fa /dev/stdin"
    )]
    #[case::mem_prefix_word(
        "tool memo x; bwa mem ref.fa -",
        "tool memo x; bwa mem -h 5,200 ref.fa -"
    )]
    #[case::no_mem("minimap2 -ax sr ref.fa -", "minimap2 -ax sr ref.fa - -h 5,200")]
    #[case::tab_before_mem("bwa\tmem ref.fa -", "bwa\tmem -h 5,200 ref.fa -")]
    // The command's own `-h` wins: the command is left unchanged.
    #[case::user_h("bwa mem -h 10,10 ref.fa -", "bwa mem -h 10,10 ref.fa -")]
    #[case::user_h_attached("bwa mem -h3 ref.fa -", "bwa mem -h3 ref.fa -")]
    #[case::user_h_in_cluster("bwa mem -Yh 10 ref.fa -", "bwa mem -Yh 10 ref.fa -")]
    #[case::user_h_in_cluster_attached("bwa mem -ph5 ref.fa -", "bwa mem -ph5 ref.fa -")]
    // `-th` would be `-t h` (a value), not `-h`; `--help` is not `-h`.
    #[case::h_as_value_of_t("bwa mem -th ref.fa -", "bwa mem -h 5,200 -th ref.fa -")]
    #[case::long_option("bwa mem --hello ref.fa -", "bwa mem -h 5,200 --hello ref.fa -")]
    #[case::h_after_pipe_ignored(
        "bwa mem ref.fa - | tool -h 2",
        "bwa mem -h 5,200 ref.fa - | tool -h 2"
    )]
    fn with_xa_max_hits_inserts_after_mem(#[case] command: &str, #[case] expected: &str) {
        let resolved = ResolvedAligner {
            backend: ResolvedBackend::Subprocess {
                command: command.to_string(),
                accept_mid_pair_split: true,
                input_format: UnalignedReadFormat::Fastq,
                max_read_ahead: None,
                runtime: SubprocessRuntime::default(),
            },
            chunk_size: 10,
            threads: Some(4),
            mode: ResolvedAlignerMode::Command,
        }
        .with_xa_max_hits(5, 200, "--xa-max-hits");
        let ResolvedBackend::Subprocess { command, .. } = resolved.backend else {
            panic!("backend kind must be unchanged");
        };
        assert_eq!(command, expected);
    }

    /// `insert_after_mem` finds the ` mem` subcommand word, or reports there is
    /// none (the caller then appends `-h` and warns).
    #[rstest]
    #[case::preset("bwa-mem3 mem -p ref.fa -", Some("bwa-mem3 mem -h 1,2 -p ref.fa -"))]
    #[case::mem_prefix_word("tool memo x", None)]
    #[case::no_mem("minimap2 -ax sr ref.fa -", None)]
    fn insert_after_mem_finds_the_subcommand(
        #[case] command: &str,
        #[case] expected: Option<&str>,
    ) {
        assert_eq!(insert_after_mem(command, "-h 1,2").as_deref(), expected);
    }

    /// The in-process backend records `-h` for its `MemOpts` instead.
    #[test]
    fn with_xa_max_hits_sets_inproc_option() {
        let resolved = ResolvedAligner {
            backend: ResolvedBackend::InProcessBwaMem3 {
                reference: PathBuf::from("ref.fa"),
                sub_batch_templates: 8,
                dedup_reads: true,
                methylation: MethylationMode::Disabled,
                xa_max_hits: None,
            },
            chunk_size: 10,
            threads: None,
            mode: ResolvedAlignerMode::Preset(AlignerPreset::BwaMem3InProc),
        }
        .with_xa_max_hits(20, 300, "--xa-max-hits");
        let ResolvedBackend::InProcessBwaMem3 { xa_max_hits, .. } = resolved.backend else {
            panic!("backend kind must be unchanged");
        };
        assert_eq!(xa_max_hits, Some((20, 300)));
    }

    /// The in-process preset links bwa-mem3, so `validate` never looks for a
    /// binary, but its index is still built by the standalone tool: the
    /// missing-index hint names `bwa-mem3 index`, not `bwa-mem3-inproc index`.
    #[test]
    fn test_validate_inproc_missing_index_hints_standalone_binary() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let msg = AlignerPreset::BwaMem3InProc.validate(&ref_path, None).unwrap_err().to_string();
        assert!(msg.contains("(run `bwa-mem3 index "), "got: {msg}");
    }

    /// The in-process preset never builds a shell command, so a reference path
    /// with a space validates for it while the subprocess presets reject it; a
    /// control character, which would corrupt the `@PG CL:` line, is rejected
    /// by every preset.
    #[rstest]
    #[case::inproc_space(AlignerPreset::BwaMem3InProc, "my ref.fa", None)]
    #[case::bwa_mem3_space(AlignerPreset::BwaMem3, "my ref.fa", Some("shell-unsafe character"))]
    #[case::bwa_space(AlignerPreset::Bwa, "my ref.fa", Some("shell-unsafe character"))]
    #[case::inproc_tab(AlignerPreset::BwaMem3InProc, "my\tref.fa", Some("control character"))]
    #[case::bwa_mem3_tab(AlignerPreset::BwaMem3, "my\tref.fa", Some("shell-unsafe character"))]
    fn test_validate_reference_path_safety(
        #[case] preset: AlignerPreset,
        #[case] file_name: &str,
        #[case] needle: Option<&str>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join(file_name);
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in preset.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let bin = make_existing_binary(tmp.path());
        let bin = preset.requires_binary().then_some(bin.as_path());
        let result = preset.validate(&ref_path, bin);
        match needle {
            None => result.unwrap(),
            Some(needle) => {
                let msg = result.unwrap_err().to_string();
                assert!(msg.contains(needle), "got: {msg}");
            }
        }
    }

    /// Without `--aligner-bin`, `validate` looks the preset's binary up on
    /// `PATH`: absent, it fails naming the binary and the override flag;
    /// present, it gets as far as the (missing) index files.
    #[test]
    fn test_validate_looks_up_binary_on_path() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let msg = AlignerPreset::BwaMem3.validate(&ref_path, None).unwrap_err().to_string();
        if which::which("bwa-mem3").is_ok() {
            assert!(msg.contains("index file not found"), "got: {msg}");
        } else {
            assert!(
                msg.contains("aligner binary `bwa-mem3` not found on PATH")
                    && msg.contains("--aligner-bin"),
                "got: {msg}"
            );
        }
    }

    /// `validate` errors out on a nonexistent `--aligner-bin` override path.
    #[test]
    fn test_validate_override_path_missing() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let missing = tmp.path().join("not-a-real-binary");
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(&missing)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("--aligner-bin path"), "got: {msg}");
        assert!(msg.contains("not a regular file"), "got: {msg}");
    }

    /// `validate` errors out when `--aligner-bin` points at a directory.
    #[test]
    fn test_validate_override_path_is_dir() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(tmp.path())).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("not a regular file"), "got: {msg}");
    }

    /// `validate` reports the FASTA-missing error BEFORE the
    /// index-missing error so users see the right fix-it hint.
    #[test]
    fn test_validate_missing_reference() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("does-not-exist.fa");
        let bin = make_existing_binary(tmp.path());
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(&bin)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("reference FASTA not found"), "got: {msg}");
        // Must NOT mention indexes — that error would be confusing.
        assert!(!msg.contains("index file not found"), "got: {msg}");
    }

    /// `validate` rejects a reference path containing a space (preset
    /// mode flows through `/bin/bash -c` so unquoted spaces break
    /// word-splitting).
    #[test]
    fn test_validate_rejects_shell_unsafe_reference() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path().join("dir with space");
        std::fs::create_dir(&ref_dir).unwrap();
        let ref_path = ref_dir.join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let bin = make_existing_binary(tmp.path());
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(&bin)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("shell-unsafe"), "got: {msg}");
        assert!(msg.contains("--aligner::command"), "should point to escape hatch: {msg}");
    }

    /// `relay_stderr` with `ring_size = 0` still drains the pipe (so the
    /// child doesn't block on stderr-write) but captures nothing, so `wait()`
    /// reports "(no stderr captured)".
    #[test]
    fn test_stderr_ring_size_zero_disables_capture() {
        let proc =
            AlignerProcess::spawn("bash -c 'echo lost-line >&2; exit 1'", 0).expect("spawn ok");
        let err = proc.wait().unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("(no stderr captured)"), "ring=0 should disable capture: {msg}");
    }

    /// `substitute_template` fills both placeholders.
    #[test]
    fn test_substitute_template_both_placeholders() {
        let reference = Path::new("/ref/genome.fa");
        let template = "bwa-mem3 mem -p -K 150000000 -t {threads} {ref} /dev/stdin";
        let result = substitute_template(template, reference, 16).unwrap();
        assert_eq!(result, "bwa-mem3 mem -p -K 150000000 -t 16 /ref/genome.fa /dev/stdin");
    }

    /// `substitute_template` fills `{ref}` even if `{threads}` is missing.
    #[test]
    fn test_substitute_template_threads_optional() {
        let reference = Path::new("/ref/genome.fa");
        let template = "bwa mem -t 8 {ref} /dev/stdin";
        let result = substitute_template(template, reference, 16).unwrap();
        assert_eq!(result, "bwa mem -t 8 /ref/genome.fa /dev/stdin");
    }

    /// `substitute_template` errors if `{ref}` is missing.
    #[test]
    fn test_substitute_template_requires_ref() {
        let reference = Path::new("/ref/genome.fa");
        let template = "bwa mem -t {threads} /dev/stdin";
        let err = substitute_template(template, reference, 16).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("{ref}"), "error should mention {{ref}}: {msg}");
    }

    /// A reference path that literally contains `{threads}` must be emitted
    /// verbatim, not rewritten into the thread count: the single-pass
    /// substitution never re-scans injected text. Regression for the
    /// replace-`{ref}`-then-replace-`{threads}` collision.
    #[test]
    fn test_substitute_template_ref_path_containing_threads_token_is_not_mangled() {
        let reference = Path::new("/data/{threads}/genome.fa");
        let template = "bwa-mem3 mem -t {threads} {ref} /dev/stdin";
        let result = substitute_template(template, reference, 8).unwrap();
        assert_eq!(result, "bwa-mem3 mem -t 8 /data/{threads}/genome.fa /dev/stdin");
    }

    /// Command mode sets its own `-K`, so a chunk size past `i32::MAX` (which
    /// only sizes its in-flight budget) is accepted.
    #[test]
    fn resolve_accepts_command_mode_chunk_size_past_i32() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            command: Some("bwa-mem3 mem -t {threads} {ref} /dev/stdin".to_string()),
            chunk_size: 2_147_483_648,
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        assert_eq!(resolved.chunk_size, 2_147_483_648);
    }

    /// `resolve` rejects `--aligner::chunk-size 0` (would reach the aligner as
    /// `-K 0` and zero the in-flight budget) with a flag-attributed error.
    #[test]
    fn test_resolve_rejects_zero_chunk_size() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let opts = AlignerOptions {
            preset: None,
            command: Some("bwa-mem3 mem {ref} /dev/stdin".to_string()),
            threads: None,
            chunk_size: 0,
            sub_batch_templates: None,
            dedup_reads: None,
            input_format: None,
            index: None,
            extra_args: None,
            startup_timeout: None,
            max_read_ahead: None,
            methylation_mode: MethylationMode::Disabled,
        };
        let err = opts.resolve(&ref_path, 4, None).unwrap_err();
        assert!(err.to_string().contains("--aligner::chunk-size must be greater than 0"));
    }

    /// `resolve` rejects an explicit `--aligner::threads 0` in preset mode
    /// (`None` means "default to runall --threads", so only a literal 0 is
    /// rejected).
    #[test]
    fn test_resolve_rejects_zero_preset_threads() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let bin = make_existing_binary(tmp.path());
        // Index files present so we reach the threads check, not the index bail.
        for ext in AlignerPreset::BwaMem3.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let opts = AlignerOptions {
            preset: Some(AlignerPreset::BwaMem3),
            command: None,
            threads: Some(0),
            chunk_size: DEFAULT_ALIGNER_CHUNK_SIZE,
            sub_batch_templates: None,
            dedup_reads: None,
            input_format: None,
            index: None,
            extra_args: None,
            startup_timeout: None,
            max_read_ahead: None,
            methylation_mode: MethylationMode::Disabled,
        };
        let err = opts.resolve(&ref_path, 4, Some(&bin)).unwrap_err();
        assert!(err.to_string().contains("--aligner::threads must be greater than 0"));
    }

    /// `validate` errors when the `--aligner-bin` override exists but is not
    /// executable (mode & 0o111 == 0).
    #[test]
    #[cfg(unix)]
    fn test_validate_rejects_non_executable_binary() {
        use std::os::unix::fs::PermissionsExt;
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let bin = tmp.path().join("mock-aligner");
        std::fs::write(&bin, b"#!/bin/sh\n").unwrap();
        let mut perms = std::fs::metadata(&bin).unwrap().permissions();
        perms.set_mode(0o644); // readable but not executable
        std::fs::set_permissions(&bin, perms).unwrap();
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(&bin)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("not executable"), "got: {msg}");
    }

    /// `check_shell_safe_path` flags the `--aligner-bin` path too (not just
    /// `--ref`), interpolating the flag label into the message.
    #[test]
    fn test_validate_rejects_shell_unsafe_aligner_bin() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        // A binary path with a space is shell-unsafe.
        let bin_dir = tmp.path().join("dir with space");
        std::fs::create_dir(&bin_dir).unwrap();
        let bin = make_existing_binary(&bin_dir);
        let err = AlignerPreset::BwaMem3.validate(&ref_path, Some(&bin)).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("shell-unsafe"), "got: {msg}");
        assert!(msg.contains("--aligner-bin"), "flag label must be named: {msg}");
    }

    /// A reference path whose interpolated word starts with `-` is rejected
    /// (the aligner would parse it as a flag).
    #[test]
    fn test_check_shell_safe_path_rejects_leading_dash() {
        let err = check_shell_safe_path(Path::new("-genome.fa"), "--ref").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("starts with '-'"), "got: {msg}");
        // A `./`-prefixed sibling of the same name is accepted.
        check_shell_safe_path(Path::new("./-genome.fa"), "--ref").unwrap();
    }

    /// `relay_stderr` retains only the last `ring_size` lines, evicting the
    /// oldest, and drains the whole stream regardless.
    #[test]
    fn test_relay_stderr_ring_cap_evicts_oldest() {
        let data = b"l1\nl2\nl3\nl4\nl5\nl6\nl7\nl8\nl9\nl10\n".to_vec();
        let capture = parking_lot::Mutex::new(StderrCapture::new(3));
        relay_stderr(std::io::Cursor::new(data), &capture);
        assert_eq!(capture.into_inner().lines, ["l8", "l9", "l10"]);
    }

    // ------------------------------------------------------------------
    // --methylation-mode (conversion-aware alignment)
    // ------------------------------------------------------------------

    /// A reference whose plain bwa-mem3 index files exist, plus the
    /// `bwa-mem3 index --meth` files when `with_meth_index`.
    fn reference_with_index(dir: &Path, with_meth_index: bool) -> PathBuf {
        let ref_path = dir.join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        let seed = append_extension(&ref_path, BWA_MEM3_METH_INDEX_SUFFIX);
        for ext in AlignerPreset::BwaMem3.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
            if with_meth_index {
                std::fs::write(append_extension(&seed, ext), b"x").unwrap();
            }
        }
        ref_path
    }

    #[rstest]
    #[case::disabled(
        AlignerPreset::BwaMem3,
        MethylationMode::Disabled,
        "bwa-mem3 mem --bam=0 -p -K 1000"
    )]
    #[case::emseq(
        AlignerPreset::BwaMem3,
        MethylationMode::EmSeq,
        "bwa-mem3 mem --bam=0 --meth --meth-scoring genomic -p -K 1000"
    )]
    #[case::taps(
        AlignerPreset::BwaMem3,
        MethylationMode::Taps,
        "bwa-mem3 mem --bam=0 --meth=taps --meth-scoring neutral -p -K 1000"
    )]
    #[case::bwa_disabled(AlignerPreset::Bwa, MethylationMode::Disabled, "bwa mem -p -K 1000")]
    #[case::bwa_taps_aligns_plain(AlignerPreset::Bwa, MethylationMode::Taps, "bwa mem -p -K 1000")]
    fn build_command_adds_the_meth_flags_for_bwa_mem3(
        #[case] preset: AlignerPreset,
        #[case] methylation: MethylationMode,
        #[case] expected_prefix: &str,
    ) {
        let cmd = preset.build_command(Path::new("/refs/ref.fa"), 4, 1000, None, methylation);
        assert!(cmd.starts_with(expected_prefix), "got: {cmd}");
    }

    /// Disabled needs no `.meth` index with any preset; neither does `bwa` under
    /// TAPS, which it aligns plain.
    #[rstest]
    #[case::disabled_bwa(AlignerPreset::Bwa, MethylationMode::Disabled)]
    #[case::disabled_bwa_mem3(AlignerPreset::BwaMem3, MethylationMode::Disabled)]
    #[case::disabled_inproc(AlignerPreset::BwaMem3InProc, MethylationMode::Disabled)]
    #[case::taps_bwa(AlignerPreset::Bwa, MethylationMode::Taps)]
    fn validate_methylation_needs_no_meth_index(
        #[case] preset: AlignerPreset,
        #[case] methylation: MethylationMode,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), false);
        preset.validate_methylation(&ref_path, methylation).unwrap();
    }

    /// A preset consumes `--methylation-mode` when it aligns conversion-aware or
    /// validation rejects it (`bwa` or `star` under EM-seq); `bwa` and `star`
    /// align TAPS plain.
    #[rstest]
    #[case::bwa_mem3_emseq(AlignerPreset::BwaMem3, MethylationMode::EmSeq, true)]
    #[case::bwa_mem3_taps(AlignerPreset::BwaMem3, MethylationMode::Taps, true)]
    #[case::inproc_taps(AlignerPreset::BwaMem3InProc, MethylationMode::Taps, true)]
    #[case::bwa_mem3_disabled(AlignerPreset::BwaMem3, MethylationMode::Disabled, false)]
    #[case::bwa_emseq_rejected(AlignerPreset::Bwa, MethylationMode::EmSeq, true)]
    #[case::bwa_taps_plain(AlignerPreset::Bwa, MethylationMode::Taps, false)]
    #[case::bwa_disabled(AlignerPreset::Bwa, MethylationMode::Disabled, false)]
    #[case::star_emseq_rejected(AlignerPreset::Star, MethylationMode::EmSeq, true)]
    #[case::star_taps_plain(AlignerPreset::Star, MethylationMode::Taps, false)]
    #[case::star_disabled(AlignerPreset::Star, MethylationMode::Disabled, false)]
    fn consumes_methylation_cases(
        #[case] preset: AlignerPreset,
        #[case] methylation: MethylationMode,
        #[case] expected: bool,
    ) {
        assert_eq!(preset.consumes_methylation(methylation), expected);
    }

    #[test]
    fn validate_methylation_rejects_bwa() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), true);
        let err = AlignerPreset::Bwa
            .validate_methylation(&ref_path, MethylationMode::EmSeq)
            .unwrap_err()
            .to_string();
        assert!(err.contains("`bwa` has none"), "got: {err}");
    }

    /// Both bwa-mem3 presets align EM-seq and TAPS against the `--meth` dual
    /// index, so both chemistries need it; the error names the mode and flag.
    #[rstest]
    fn validate_methylation_requires_the_meth_dual_index(
        #[values(AlignerPreset::BwaMem3, AlignerPreset::BwaMem3InProc)] preset: AlignerPreset,
        #[values((MethylationMode::EmSeq, "em-seq", "--meth,"), (MethylationMode::Taps, "taps", "--meth=taps,"))]
        mode: (MethylationMode, &str, &str),
    ) {
        let (methylation, mode_name, flag) = mode;
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), false);
        let err = preset.validate_methylation(&ref_path, methylation).unwrap_err().to_string();
        assert!(err.contains("ref.fa.meth."), "names the missing .meth file; got: {err}");
        assert!(err.contains("bwa-mem3 index --meth --force"), "gives the fix; got: {err}");
        assert!(err.contains("before 0.18.0"), "says which releases lack --force; got: {err}");
        assert!(err.contains(&format!("--methylation-mode {mode_name} ")), "got: {err}");
        assert!(err.contains(&format!("bwa-mem3 {flag}")), "got: {err}");
        // Only TAPS, whose reads also align plain, points at the plain fallback.
        assert_eq!(
            err.contains("--aligner::preset bwa (which needs a bwa index)"),
            methylation == MethylationMode::Taps,
            "got: {err}"
        );
        assert_eq!(
            err.contains("running plain `bwa-mem3 mem`"),
            methylation == MethylationMode::Taps,
            "got: {err}"
        );

        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), true);
        preset.validate_methylation(&ref_path, methylation).unwrap();
    }

    /// A reference with a `--meth-only` index ([`METH_ONLY_INDEX_FILES`]),
    /// leaving out `omit` when given.
    fn reference_with_meth_only_index(dir: &Path, omit: Option<&str>) -> PathBuf {
        let ref_path = dir.join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in METH_ONLY_INDEX_FILES.iter().filter(|ext| Some(**ext) != omit) {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        ref_path
    }

    /// `bwa-mem3 mem --meth` and `--meth=taps` never read the plain
    /// `<ref>.bwt.2bit.64`, so the subprocess preset accepts a `--meth-only`
    /// index under EM-seq and TAPS.
    #[rstest]
    #[case::emseq(MethylationMode::EmSeq)]
    #[case::taps(MethylationMode::Taps)]
    fn validate_with_methylation_accepts_a_meth_only_index(#[case] methylation: MethylationMode) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_meth_only_index(tmp.path(), None);
        let bin = make_existing_binary(tmp.path());
        AlignerPreset::BwaMem3
            .validate_with_methylation(&ref_path, Some(&bin), methylation)
            .unwrap();
    }

    /// Under EM-seq and TAPS every file of the `--meth-only` layout is still
    /// required; the error names the missing file and suggests the full dual
    /// index (not `--meth-only`, which cannot align plain), with `--force` to
    /// replace one.
    #[rstest]
    fn validate_with_methylation_requires_each_meth_only_file(
        #[values(MethylationMode::EmSeq, MethylationMode::Taps)] methylation: MethylationMode,
        #[values(
            ".amb",
            ".ann",
            ".pac",
            ".meth.amb",
            ".meth.ann",
            ".meth.bwt.2bit.64",
            ".meth.pac"
        )]
        missing: &str,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_meth_only_index(tmp.path(), Some(missing));
        let bin = make_existing_binary(tmp.path());
        let err = AlignerPreset::BwaMem3
            .validate_with_methylation(&ref_path, Some(&bin), methylation)
            .unwrap_err()
            .to_string();
        assert!(err.contains(&format!("ref.fa{missing} ")), "names {missing}; got: {err}");
        assert!(
            err.contains("run `bwa-mem3 index --meth --force ") && err.contains("before 0.18.0"),
            "gives the fix; got: {err}"
        );
        assert!(!err.contains("--meth-only"), "does not suggest --meth-only; got: {err}");
    }

    /// Plain alignment reads the plain index, so a `--meth-only` index (no
    /// `<ref>.bwt.2bit.64`) is still rejected for it, as it is for `bwa` under
    /// TAPS, which aligns plain.
    #[rstest]
    #[case::disabled(AlignerPreset::BwaMem3, MethylationMode::Disabled)]
    #[case::bwa_taps(AlignerPreset::Bwa, MethylationMode::Taps)]
    fn validate_with_methylation_requires_the_plain_index_when_aligning_plain(
        #[case] preset: AlignerPreset,
        #[case] methylation: MethylationMode,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_meth_only_index(tmp.path(), None);
        let bin = make_existing_binary(tmp.path());
        let err = preset
            .validate_with_methylation(&ref_path, Some(&bin), methylation)
            .unwrap_err()
            .to_string();
        assert!(err.contains("required index file not found: "), "got: {err}");
        assert!(!err.contains("--meth"), "plain mode suggests a plain index; got: {err}");
    }

    /// The in-process backend's index loader checks for the plain
    /// `<ref>.bwt.2bit.64` even for a dual index, so it rejects a `--meth-only`
    /// index up front under EM-seq and TAPS, says it needs the full dual index,
    /// and suggests building it (with `--force` to replace the existing files).
    #[rstest]
    #[case::emseq(MethylationMode::EmSeq)]
    #[case::taps(MethylationMode::Taps)]
    fn validate_with_methylation_inproc_rejects_a_meth_only_index(
        #[case] methylation: MethylationMode,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_meth_only_index(tmp.path(), None);
        let err = AlignerPreset::BwaMem3InProc
            .validate_with_methylation(&ref_path, None, methylation)
            .unwrap_err()
            .to_string();
        assert!(err.contains("ref.fa.bwt.2bit.64 "), "names the plain index; got: {err}");
        assert!(
            err.contains("needs the full dual index, including the plain `"),
            "explains why; got: {err}"
        );
        assert!(err.contains("index --meth --force"), "mentions --force; got: {err}");
    }

    /// `bwa` is rejected under EM-seq before its binary or index is checked:
    /// no binary or index would make it bisulfite-aware.
    #[test]
    fn validate_with_methylation_rejects_bwa_first() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("missing.fa");
        let err = AlignerPreset::Bwa
            .validate_with_methylation(&ref_path, None, MethylationMode::EmSeq)
            .unwrap_err()
            .to_string();
        assert!(err.contains("`bwa` has none"), "got: {err}");
    }

    /// The subprocess presets carry `--methylation-mode` into their command:
    /// bwa-mem3 gets `--meth --meth-scoring genomic` for EM-seq and
    /// `--meth=taps --meth-scoring neutral` for TAPS (both need the `.meth`
    /// index); `bwa` aligns TAPS plain and needs none.
    #[rstest]
    #[case::emseq(
        AlignerPreset::BwaMem3,
        MethylationMode::EmSeq,
        true,
        Some("--meth --meth-scoring genomic")
    )]
    #[case::taps(
        AlignerPreset::BwaMem3,
        MethylationMode::Taps,
        true,
        Some("--meth=taps --meth-scoring neutral")
    )]
    #[case::bwa_taps(AlignerPreset::Bwa, MethylationMode::Taps, false, None)]
    fn resolve_preset_adds_the_meth_flags(
        #[case] preset: AlignerPreset,
        #[case] methylation: MethylationMode,
        #[case] with_meth_index: bool,
        #[case] expected_flags: Option<&str>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), with_meth_index);
        if preset == AlignerPreset::Bwa {
            for ext in AlignerPreset::Bwa.index_extensions() {
                std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
            }
        }
        let bin = make_existing_binary(tmp.path());
        let opts = AlignerOptions {
            preset: Some(preset),
            methylation_mode: methylation,
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, Some(&bin)).unwrap();
        let ResolvedBackend::Subprocess { command, .. } = resolved.backend else {
            panic!("a subprocess preset resolves to a subprocess");
        };
        match expected_flags {
            Some(flags) => assert!(command.contains(&format!(" {flags} ")), "got: {command}"),
            None => assert!(!command.contains("--meth"), "got: {command}"),
        }
    }

    /// `--aligner::preset bwa-mem3` resolves against a
    /// `bwa-mem3 index --meth --meth-only` index under EM-seq and TAPS.
    #[rstest]
    #[case::emseq(MethylationMode::EmSeq, "--meth --meth-scoring genomic")]
    #[case::taps(MethylationMode::Taps, "--meth=taps --meth-scoring neutral")]
    fn resolve_preset_accepts_a_meth_only_index(
        #[case] methylation: MethylationMode,
        #[case] flags: &str,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_meth_only_index(tmp.path(), None);
        let bin = make_existing_binary(tmp.path());
        let opts = AlignerOptions {
            preset: Some(AlignerPreset::BwaMem3),
            methylation_mode: methylation,
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, Some(&bin)).unwrap();
        let ResolvedBackend::Subprocess { command, .. } = resolved.backend else {
            panic!("the bwa-mem3 preset resolves to a subprocess");
        };
        assert!(command.contains(&format!(" {flags} ")), "got: {command}");
    }

    /// TAPS with the bwa-mem3 preset aligns with `--meth=taps`, so a reference
    /// without the `.meth` index is rejected at resolve time.
    #[test]
    fn resolve_taps_preset_requires_the_meth_index() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), false);
        let bin = make_existing_binary(tmp.path());
        let opts = AlignerOptions {
            preset: Some(AlignerPreset::BwaMem3),
            methylation_mode: MethylationMode::Taps,
            ..AlignerOptions::default()
        };
        let err = opts.resolve(&ref_path, 4, Some(&bin)).unwrap_err().to_string();
        assert!(err.contains("bwa-mem3 index --meth"), "got: {err}");
    }

    /// `bwa` cannot align bisulfite reads, so `--methylation-mode em-seq` rejects it at
    /// resolve time rather than emitting an unconverted alignment.
    #[test]
    fn resolve_rejects_bwa_under_methylation_mode() {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in AlignerPreset::Bwa.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let bin = make_existing_binary(tmp.path());
        let opts = AlignerOptions {
            preset: Some(AlignerPreset::Bwa),
            methylation_mode: MethylationMode::EmSeq,
            ..AlignerOptions::default()
        };
        let err = opts.resolve(&ref_path, 4, Some(&bin)).unwrap_err().to_string();
        assert!(err.contains("`bwa` has none"), "got: {err}");
    }

    #[rstest]
    #[cfg(feature = "aligner-bwa-mem3")]
    fn resolve_inproc_carries_methylation_mode(
        #[values(MethylationMode::EmSeq, MethylationMode::Taps)] methylation: MethylationMode,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = reference_with_index(tmp.path(), true);
        let opts = AlignerOptions {
            preset: Some(AlignerPreset::BwaMem3InProc),
            methylation_mode: methylation,
            ..AlignerOptions::default()
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        assert!(
            matches!(
                resolved.backend,
                ResolvedBackend::InProcessBwaMem3 { methylation: m, .. } if m == methylation
            ),
            "got: {:?}",
            resolved.backend
        );
    }
}
