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
use std::thread::{self, JoinHandle};
use std::time::Duration;

use anyhow::{Context, Result, bail};
use clap::Args;
use fgumi_cli_macros::multi_options;

// ============================================================================
// AlignerProcess
// ============================================================================

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
    /// The child's stdin pipe, available until [`Self::take_stdin`] is called.
    stdin: Option<ChildStdin>,
    /// The child's stdout pipe, available until [`Self::take_stdout`] is called.
    stdout: Option<ChildStdout>,
    /// Background thread that relays child stderr; returns the captured ring buffer.
    stderr_thread: Option<JoinHandle<Vec<String>>>,
}

/// The mimalloc option environment variables (`MIMALLOC_PURGE_DELAY` and its
/// legacy name) a user can set to choose the purge delay themselves.
const MIMALLOC_PURGE_ENV: [&str; 2] = ["MIMALLOC_PURGE_DELAY", "MIMALLOC_RESET_DELAY"];

/// Whether the user set mimalloc's purge delay in the environment.
pub(crate) fn user_set_mimalloc_purge() -> bool {
    names_set_mimalloc_purge(std::env::vars_os().map(|(name, _)| name))
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
        // after 1 s and costs the aligner page faults on every batch: never
        // purging measured -0.5% wall on the bwa-mem3 CLI alone. A user-set
        // value (either name, any case) is inherited unchanged; aligners
        // without mimalloc ignore the variable. The variable reaches every
        // process in the shell command, so with `--aligner::command` any
        // mimalloc-linked tool in the user's pipeline also keeps its freed
        // pages (higher RSS); set `MIMALLOC_PURGE_DELAY` to opt out.
        if !user_set_mimalloc_purge() {
            cmd.env("MIMALLOC_PURGE_DELAY", "-1");
        }
        let mut child =
            cmd.spawn().with_context(|| format!("failed to spawn aligner command: {command}"))?;

        let stdin = child.stdin.take();
        let stdout = child.stdout.take();
        let child_stderr = child.stderr.take().expect("stderr was configured as piped");

        let stderr_thread = thread::spawn(move || relay_stderr(child_stderr, ring_size));

        Ok(Self { child, stdin, stdout, stderr_thread: Some(stderr_thread) })
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
    /// non-zero status.
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

        let status = self.child.wait().context("failed to wait for aligner process")?;

        // Join the stderr relay thread so its resources are cleaned up.
        // A panic in the relay thread is logged (not silently dropped)
        // so operators can tell apart "child produced no stderr" from
        // "our relay died" in the error message.
        let last_lines = match self.stderr_thread.take().map(JoinHandle::join) {
            Some(Ok(lines)) => lines,
            Some(Err(_)) => {
                log::warn!("aligner stderr relay thread panicked; captured ring may be incomplete");
                Vec::new()
            }
            None => Vec::new(),
        };

        if !status.success() {
            let code = status.code().map_or_else(|| "signal".to_string(), |c| c.to_string());
            let tail = if last_lines.is_empty() {
                "(no stderr captured)".to_string()
            } else {
                last_lines.join("\n")
            };
            bail!("aligner process exited with status {code}. Last stderr lines:\n{tail}");
        }

        Ok(())
    }

    /// Kill the aligner process.
    ///
    /// Sends SIGKILL immediately (`Child::kill()` always sends SIGKILL
    /// on Unix). Waits up to 1 second for the process to exit.
    ///
    /// Returns `true` if the child was reaped within the 1s deadline (or was
    /// already gone), and `false` if it was still running when the deadline
    /// expired (i.e. stuck, e.g. in uninterruptible I/O). Callers MUST honor a
    /// `false` return by NOT issuing a subsequent unbounded [`wait`](Self::wait)
    /// on the child, which would re-introduce the very teardown hang this
    /// bounded kill exists to prevent.
    pub fn kill(&mut self) -> bool {
        let pid = self.child.id();
        let _ = self.child.kill();

        // Wait for the process to actually exit so its resources are cleaned up.
        let deadline = std::time::Instant::now() + Duration::from_secs(1);
        // Tracks whether the loop gave up because the 1s deadline elapsed while
        // the process was still alive. Only that case warrants the leak warning.
        let mut deadline_exceeded = false;
        // Set only when exit is *positively confirmed*: a clean reap
        // (`Ok(Some)`) or an `ECHILD` error (already reaped by a SIGCHLD
        // handler). A deadline timeout or an unclassified `try_wait` errno
        // leaves this `false`, so the caller must NOT follow up with an
        // unbounded `wait()`.
        let mut reaped = false;
        loop {
            match self.child.try_wait() {
                // Reaped cleanly: resources cleaned up, no warning.
                Ok(Some(_)) => {
                    reaped = true;
                    break;
                }
                // Still alive: keep polling until the deadline.
                Ok(None) => {
                    if std::time::Instant::now() >= deadline {
                        deadline_exceeded = true;
                        break;
                    }
                    thread::sleep(Duration::from_millis(50));
                }
                // `try_wait` failed. `ECHILD` (the child was already reaped,
                // e.g. by a SIGCHLD handler) is a benign non-leak that confirms
                // exit. Any *other* errno is unexpected — the process may still
                // be around — so we do NOT report it as reaped; we cannot make
                // progress on a `try_wait` error either way, so stop polling.
                Err(err) => {
                    if is_already_reaped(&err) {
                        log::debug!(
                            "aligner process (pid {pid}) try_wait returned ECHILD after \
                             SIGKILL; already reaped"
                        );
                        reaped = true;
                    } else {
                        log::warn!(
                            "aligner process (pid {pid}) try_wait failed after SIGKILL: {err}; \
                             unable to confirm exit, treating as not reaped"
                        );
                    }
                    break;
                }
            }
        }
        if deadline_exceeded {
            log::warn!(
                "aligner process (pid {pid}) did not exit within 1s of SIGKILL; it may be \
                 stuck (e.g. uninterruptible I/O) and left unreaped"
            );
        }
        // Report reaped only when exit was positively confirmed (clean reap or
        // `ECHILD`). A timeout on a still-running child or an unclassified
        // `try_wait` errno reports `false` so the caller avoids a follow-up
        // unbounded `wait()`.
        reaped
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

impl Drop for AlignerProcess {
    fn drop(&mut self) {
        // Determine whether the child's exit is *confirmed*. Only then is it
        // safe to join the stderr relay thread: that thread blocks reading the
        // child's stderr fd, which stays open while the child lives, so joining
        // a still-alive child would hang teardown forever — the very hang the
        // bounded `kill()` exists to prevent.
        let exit_confirmed = match self.child.try_wait() {
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
            Err(err) => {
                if is_already_reaped(&err) {
                    true
                } else {
                    log::warn!(
                        "aligner process (pid {}) try_wait failed in Drop: {err}; \
                         attempting bounded kill",
                        self.child.id()
                    );
                    self.kill()
                }
            }
        };
        // Join the stderr relay only when exit is confirmed. If the child could
        // not be confirmed dead (kill timed out / unclassified errno), skip the
        // join: leaking a detached relay thread is the lesser evil versus
        // hanging Drop indefinitely on a still-open stderr fd.
        if let Some(handle) = self.stderr_thread.take()
            && exit_confirmed
        {
            let _ = handle.join();
        }
    }
}

/// Read lines from `stderr`, print each to fgumi's stderr, and retain
/// the last `ring_size` lines in a [`VecDeque`] ring buffer.
///
/// Returns the captured lines when the child closes its stderr fd.
fn relay_stderr(stderr: impl std::io::Read, ring_size: usize) -> Vec<String> {
    let mut reader = BufReader::new(stderr);
    let mut ring: VecDeque<String> = VecDeque::with_capacity(ring_size);
    let mut buf: Vec<u8> = Vec::new();

    // Byte-oriented drain: `BufRead::lines()` yields `Err` on a non-UTF-8 line,
    // and breaking on it would STOP draining the aligner's stderr — the aligner
    // then blocks once its stderr pipe fills (~64 KiB) and the whole align
    // pipeline deadlocks. Read raw lines and decode lossily so invalid UTF-8
    // never halts the relay.
    loop {
        buf.clear();
        match reader.read_until(b'\n', &mut buf) {
            // EOF (child closed stderr) or an unrecoverable pipe read error:
            // stop draining either way.
            Ok(0) | Err(_) => break,
            Ok(_) => {}
        }
        // Strip the line terminator to match `lines()` semantics.
        if buf.last() == Some(&b'\n') {
            buf.pop();
            if buf.last() == Some(&b'\r') {
                buf.pop();
            }
        }
        let line = String::from_utf8_lossy(&buf).into_owned();
        eprintln!("{line}");
        if ring_size > 0 {
            if ring.len() == ring_size {
                ring.pop_front();
            }
            ring.push_back(line);
        }
    }

    ring.into_iter().collect()
}

// ============================================================================
// AlignerPreset
// ============================================================================

/// Preset aligner configurations with command-building and index validation.
///
/// Three presets at landing — `bwa-mem3` (the project's primary
/// subprocess aligner), `bwa` (legacy reference), and `bwa-mem3-inproc`
/// (bwa-mem3 linked directly into fgumi, gated by the `aligner-bwa-mem3`
/// build feature; see `AlignerOptions::resolve`). Methylation-aware
/// presets (e.g. `bwameth`, `bwa-mem3 --methylation-mode em-seq`) are a
/// follow-up; EM-seq users today route through `--aligner::command "..."`
/// (free-form mode) instead of a preset.
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
}

impl std::fmt::Display for AlignerPreset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.binary_name())
    }
}

impl AlignerPreset {
    /// Default binary name (without path) for this preset. Also the
    /// kebab-case CLI value (`bwa`, `bwa-mem3`). Internal helper —
    /// callers outside this module should use [`Display`](std::fmt::Display) instead.
    #[must_use]
    pub(crate) fn binary_name(self) -> &'static str {
        match self {
            Self::Bwa => "bwa",
            Self::BwaMem3 => "bwa-mem3",
            Self::BwaMem3InProc => "bwa-mem3-inproc",
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
        }
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
    // `pub(crate)`, not `pub`: this interpolates `reference` / `binary_override`
    // into a string later run via `/bin/bash -c`, and it does not itself validate
    // those paths — it trusts that `validate` (which runs `check_shell_safe_path`)
    // was called first. Restricting the method to in-crate callers keeps a direct
    // external caller from bypassing that validation and injecting shell syntax.
    #[must_use]
    pub(crate) fn build_command(
        self,
        reference: &Path,
        threads: usize,
        chunk_size: u64,
        binary_override: Option<&Path>,
    ) -> String {
        let binary = match binary_override {
            Some(path) => path.display().to_string(),
            None => self.binary_name().to_string(),
        };
        let ref_str = reference.display();
        let output = match self {
            Self::BwaMem3 | Self::BwaMem3InProc => " --bam=0",
            Self::Bwa => "",
        };
        format!("{binary} mem{output} -p -K {chunk_size} -t {threads} {ref_str} /dev/stdin")
    }

    /// Validate that the aligner binary and required index files are present.
    ///
    /// Checks, in order:
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
    /// # Errors
    ///
    /// Returns an error with a fix-it hint if any of the above fail.
    pub fn validate(self, reference: &Path, binary_override: Option<&Path>) -> Result<()> {
        let binary_name = self.binary_name();

        // 1. Reference FASTA presence (separated from the index-loop
        //    error so users see "ref doesn't exist" before "index
        //    missing" when the path is wrong).
        if !reference.is_file() {
            bail!(
                "reference FASTA not found or not a regular file: {} \
                 (preset `{binary_name}` requires the FASTA + its index files)",
                reference.display()
            );
        }

        // 2. Reference path shell-safety (preset mode argv flows through
        //    /bin/bash -c, so paths with spaces, quotes, or shell metas
        //    silently break). The in-process backend loads the index
        //    through the library, so only the `@PG CL:` line constrains it.
        if self.requires_binary() {
            check_shell_safe_path(reference, "--ref")?;
        } else {
            check_header_safe_path(reference, "--ref")?;
        }

        // 3. Binary discovery — skipped for presets that don't need one.
        if self.requires_binary() {
            match binary_override {
                Some(path) => {
                    check_shell_safe_path(path, "--aligner-bin")?;
                    if !path.is_file() {
                        bail!(
                            "--aligner-bin path is not a regular file: {} \
                             (preset `{binary_name}` requires the binary)",
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
                                 (preset `{binary_name}` requires an executable binary; \
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
                             (preset `{binary_name}`; pass --aligner-bin <path> to override)"
                        )
                    })?;
                }
            }
        }

        // 4. Index files for the chosen preset.
        for ext in self.index_extensions() {
            let index_path = append_extension(reference, ext);
            if !index_path.is_file() {
                bail!(
                    "required index file not found: {} (run `{} index {}`)",
                    index_path.display(),
                    self.index_binary_hint(),
                    reference.display()
                );
            }
        }

        Ok(())
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
#[multi_options("aligner", "Aligner Options")]
#[derive(Args, Debug, Clone)]
pub struct AlignerOptions {
    /// Named aligner preset. Mutually exclusive with `--aligner::command`.
    #[arg(long = "preset", value_enum)]
    pub preset: Option<AlignerPreset>,

    /// Free-form aligner command template. Supports `{ref}` and
    /// `{threads}` placeholders. Mutually exclusive with
    /// `--aligner::preset` and every other `--aligner::*` flag.
    #[arg(long = "command")]
    pub command: Option<String>,

    /// Thread count passed to the aligner via `-t` (preset mode only).
    /// Defaults to runall's top-level `--threads` value when unset.
    /// The aligner uses these threads concurrently with the rest of
    /// the AAM pipeline; oversubscription is the OS scheduler's
    /// problem.
    #[arg(long = "threads")]
    pub threads: Option<usize>,

    /// Bases per aligner batch (preset mode `-K` flag; also drives the
    /// Step-1 batch size in both modes).
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
        }
    }
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
    },
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
    /// - Command mode + a preset-only flag (`--aligner-bin` or
    ///   `--aligner::threads`).
    /// - Command mode template missing `{ref}`.
    /// - `--aligner::preset bwa-mem3-inproc` + `--aligner::threads` or
    ///   `--aligner-bin` (the in-process backend has no subprocess and
    ///   shares runall's thread budget).
    /// - `--aligner::preset bwa-mem3-inproc` built without the
    ///   `aligner-bwa-mem3` feature.
    /// - Preset-mode index files / binary missing (delegated to
    ///   [`AlignerPreset::validate`]).
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
        if self.preset.is_some() && self.chunk_size > MAX_PRESET_CHUNK_SIZE {
            bail!(
                "--aligner::chunk-size must be at most {MAX_PRESET_CHUNK_SIZE} (got {}) with \
                 --aligner::preset: bwa parses -K as a 32-bit int",
                self.chunk_size
            );
        }
        // `--aligner::sub-batch-templates` is a hidden, in-process-only knob.
        // Reject it early (before the preset/command dispatch below) so it is
        // rejected uniformly regardless of mode, and regardless of whether
        // this binary was built with `aligner-bwa-mem3`.
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
        match (self.preset, self.command) {
            (None, None) => bail!(
                "--start-from align requires one of `--aligner::preset` \
                 (e.g. bwa-mem3, bwa) or `--aligner::command \"...\"`"
            ),
            (Some(_), Some(_)) => bail!(
                "--aligner::preset and --aligner::command are mutually exclusive; \
                 pass only one"
            ),
            (Some(AlignerPreset::BwaMem3InProc), None) => resolve_inproc(
                reference,
                self.threads,
                aligner_bin,
                self.sub_batch_templates,
                self.dedup_reads,
                self.chunk_size,
            ),
            (Some(preset), None) => {
                // Preset mode: validate indexes + binary, build argv.
                preset.validate(reference, aligner_bin)?;
                // An explicit `--aligner::threads 0` would reach the aligner as
                // `-t 0`. `None` means "default to runall's --threads" (not
                // zero), so only reject an explicit zero.
                if self.threads == Some(0) {
                    bail!("--aligner::threads must be greater than 0 (got 0)");
                }
                let threads = self.threads.unwrap_or(top_threads);
                let command =
                    preset.build_command(reference, threads, self.chunk_size, aligner_bin);
                Ok(ResolvedAligner {
                    // Both subprocess presets run `mem -p -K`, whose smart
                    // pairing splits a pair across a chunk cut (bwa and
                    // bwa-mem3 share `bseq_read`'s even-count cut and
                    // `bseq_classify`).
                    backend: ResolvedBackend::Subprocess { command, accept_mid_pair_split: true },
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
                    backend: ResolvedBackend::Subprocess { command, accept_mid_pair_split: false },
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
        AlignerPreset::BwaMem3InProc.validate(reference, None)?;
        let sub_batch_templates = sub_batch_templates.unwrap_or_else(default_sub_batch_templates);
        Ok(ResolvedAligner {
            backend: ResolvedBackend::InProcessBwaMem3 {
                reference: reference.to_path_buf(),
                sub_batch_templates,
                dedup_reads: dedup_reads != Some(DedupReads::Off),
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

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    use rstest::rstest;

    use super::*;

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

        let ring = relay_stderr(std::io::Cursor::new(data), 8);

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
    /// reaps it well inside the 1s window — and a `true` return is what lets
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

    /// `is_already_reaped` confirms exit only for `ECHILD`. An unrelated errno
    /// (e.g. `EINVAL`) must NOT be treated as a confirmed reap, so `kill()`
    /// reports `false` and the caller avoids a follow-up unbounded `wait()`.
    #[test]
    fn test_is_already_reaped_only_for_echild() {
        let echild = std::io::Error::from_raw_os_error(ECHILD);
        assert!(is_already_reaped(&echild), "ECHILD must count as already reaped");

        // EINVAL (22) stands in for any unrelated errno: it must NOT count.
        let einval = std::io::Error::from_raw_os_error(22);
        assert!(!is_already_reaped(&einval), "an unrelated errno must not count as reaped");

        // A non-OS error (no errno at all) is likewise not a confirmed reap.
        let not_os = std::io::Error::other("synthetic non-os error");
        assert!(!is_already_reaped(&not_os), "a non-OS error must not count as reaped");
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
        let cmd = preset.build_command(reference, threads, chunk_size, None);
        assert_eq!(cmd, expected);
    }

    /// `--aligner-bin` override replaces the bare binary name.
    #[test]
    fn test_build_command_with_binary_override() {
        let reference = Path::new("/data/hg38.fa");
        let override_path = PathBuf::from("/opt/bwa-mem3-2.2.1/bwa-mem3");
        let cmd =
            AlignerPreset::BwaMem3.build_command(reference, 16, 100_000, Some(&override_path));
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
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        let ResolvedBackend::Subprocess { command, .. } = resolved.backend else {
            panic!("command mode must resolve to ResolvedBackend::Subprocess");
        };
        assert!(command.contains(&ref_path.display().to_string()));
        assert!(command.contains("-t 4"));
        assert!(matches!(resolved.mode, ResolvedAlignerMode::Command));
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
    fn resolve_preset_accepts_mid_pair_split(#[case] preset: AlignerPreset) {
        let tmp = tempfile::tempdir().unwrap();
        let ref_path = tmp.path().join("ref.fa");
        std::fs::write(&ref_path, b">chr1\nACGT\n").unwrap();
        for ext in preset.index_extensions() {
            std::fs::write(append_extension(&ref_path, ext), b"x").unwrap();
        }
        let bin = make_existing_binary(tmp.path());
        let opts = AlignerOptions { preset: Some(preset), ..AlignerOptions::default() };
        let resolved = opts.resolve(&ref_path, 4, Some(&bin)).unwrap();
        let ResolvedBackend::Subprocess { command, accept_mid_pair_split } = resolved.backend
        else {
            panic!("a subprocess preset must resolve to ResolvedBackend::Subprocess");
        };
        assert!(accept_mid_pair_split, "preset mode must accept bwa's mid-pair split");
        assert!(command.starts_with(&bin.display().to_string()), "got: {command}");
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

    /// `relay_stderr` with `ring_size = 0` still drains the pipe (so
    /// the child doesn't block on stderr-write) but returns an empty
    /// ring buffer.
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
        let ring = relay_stderr(std::io::Cursor::new(data), 3);
        assert_eq!(ring, vec!["l8".to_string(), "l9".to_string(), "l10".to_string()]);
    }
}
