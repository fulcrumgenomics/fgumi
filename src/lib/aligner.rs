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
//! Used by the `AlignAndMergeStep` in
//! `src/lib/pipeline/steps/align_and_merge.rs`. This module is the
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
        let mut child = Command::new("/bin/bash")
            .args(["-c", command])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .with_context(|| format!("failed to spawn aligner command: {command}"))?;

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
        // the process was still alive. Only that case warrants a warning; a
        // successful reap or a non-leak `try_wait` error (e.g. ECHILD) does not.
        let mut deadline_exceeded = false;
        // `Ok(None)` means still alive — keep polling. Any other result
        // (`Ok(Some)` = reaped, or `Err` such as ECHILD = already reaped) is a
        // non-leak: exit without warning. Only hitting the deadline while still
        // alive sets `deadline_exceeded`.
        while let Ok(None) = self.child.try_wait() {
            if std::time::Instant::now() >= deadline {
                deadline_exceeded = true;
                break;
            }
            thread::sleep(Duration::from_millis(50));
        }
        if deadline_exceeded {
            log::warn!(
                "aligner process (pid {pid}) did not exit within 1s of SIGKILL; it may be \
                 stuck (e.g. uninterruptible I/O) and left unreaped"
            );
        }
        reaped
    }

    /// Return the process ID of the aligner child process.
    #[must_use]
    pub fn pid(&self) -> u32 {
        self.child.id()
    }
}

impl Drop for AlignerProcess {
    fn drop(&mut self) {
        // If the child process is still running when the struct is
        // dropped (e.g. due to a panic or early return), kill it to
        // prevent orphaned processes.
        if let Ok(None) = self.child.try_wait() {
            self.kill();
        }
        // Always join the stderr relay thread to avoid resource leaks.
        if let Some(handle) = self.stderr_thread.take() {
            let _ = handle.join();
        }
    }
}

/// Read lines from `stderr`, print each to fgumi's stderr, and retain
/// the last `ring_size` lines in a [`VecDeque`] ring buffer.
///
/// Returns the captured lines when the child closes its stderr fd.
fn relay_stderr(stderr: impl std::io::Read, ring_size: usize) -> Vec<String> {
    let reader = BufReader::new(stderr);
    let mut ring: VecDeque<String> = VecDeque::with_capacity(ring_size);

    for line in reader.lines() {
        let Ok(line) = line else { break };
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
/// Two presets at landing — `bwa-mem3` (the project's primary aligner)
/// and `bwa` (legacy reference). Methylation-aware presets (e.g.
/// `bwameth`, `bwa-mem3 --methylation-mode em-seq`) are a follow-up;
/// EM-seq users today route through `--aligner::command "..."`
/// (free-form mode) instead of a preset.
///
/// The CLI value for each variant matches the on-disk binary name
/// (`bwa`, `bwa-mem3`). Variant identifiers (`Bwa`, `BwaMem3`) are
/// chosen so `rename_all = "kebab-case"` round-trips through the
/// identifier ↔ binary-name mapping without divergence: `Bwa` ↔
/// `"bwa"`, `BwaMem3` ↔ `"bwa-mem3"`. `Display` emits the same
/// strings so error messages match what the user typed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum AlignerPreset {
    /// Classic BWA aligner (`bwa mem`). CLI value: `bwa`.
    Bwa,
    /// BWA-MEM3 aligner (`bwa-mem3 mem`). CLI value: `bwa-mem3`.
    BwaMem3,
}

impl std::fmt::Display for AlignerPreset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.binary_name())
    }
}

impl AlignerPreset {
    /// Default binary name (without path) for this preset. Also the
    /// kebab-case CLI value (`bwa`, `bwa-mem3`). Internal helper —
    /// callers outside this module should use [`Display`] instead.
    #[must_use]
    pub(crate) fn binary_name(self) -> &'static str {
        match self {
            Self::Bwa => "bwa",
            Self::BwaMem3 => "bwa-mem3",
        }
    }

    /// Index file extensions this preset requires alongside the reference
    /// FASTA. All must be present for `validate` to succeed.
    #[must_use]
    pub(crate) fn index_extensions(self) -> &'static [&'static str] {
        match self {
            // BWA classic: amb / ann / bwt / pac / sa
            Self::Bwa => &[".amb", ".ann", ".bwt", ".pac", ".sa"],
            // BWA-MEM3: amb / ann / bwt.2bit.64 / pac. `bwa-mem3 0.4.0+` pac-fetches
            // the reference from `.pac` on demand and no longer writes `.0123` by
            // default (fg-labs/bwa-mem3#177); `.0123` is now opt-in via `index
            // --emit-unpacked-ref` and `mem` ignores any present, so requiring it
            // here rejects a valid 0.4.0 index. `.bwt.2bit.64` is the canonical
            // index sentinel. Indexes built by older bwa-mem3 still carry these
            // four files, so dropping `.0123` is backward-compatible.
            Self::BwaMem3 => &[".amb", ".ann", ".bwt.2bit.64", ".pac"],
        }
    }

    /// Build the aligner shell command for this preset.
    ///
    /// The command reads FASTQ from `/dev/stdin` and writes SAM to
    /// stdout, allowing it to be connected to fgumi's pipeline via
    /// stdin/stdout pipes.
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
    #[must_use]
    pub fn build_command(
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
        format!("{binary} mem -p -K {chunk_size} -t {threads} {ref_str} /dev/stdin")
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
    ///    two arguments).
    /// 3. The aligner binary is reachable: `--aligner-bin` override
    ///    path (must be an existing file), else `which::which("<binary_name>")`
    ///    on `PATH`. `--aligner-bin` paths are also checked for
    ///    shell-unsafe characters.
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
        //    silently break).
        check_shell_safe_path(reference, "--ref")?;

        // 3. Binary discovery.
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
                            format!("reading metadata for --aligner-bin path: {}", path.display())
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

        // 4. Index files for the chosen preset.
        for ext in self.index_extensions() {
            let index_path = append_extension(reference, ext);
            if !index_path.is_file() {
                bail!(
                    "required index file not found: {} (run `{binary_name} index {}`)",
                    index_path.display(),
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
             \"bwa-mem3 mem -p -K 150000000 -t {{threads}} {{ref}} /dev/stdin\""
        );
    }
    let with_ref = template.replace("{ref}", &reference.display().to_string());
    Ok(with_ref.replace("{threads}", &threads.to_string()))
}

// ============================================================================
// AlignerOptions — runall CLI surface
// ============================================================================

/// Default chunk size in bases per aligner batch (`-K` flag). Matches the
/// `bwa mem -K` examples in the bwa-mem3 documentation; large enough that
/// the aligner sees full batches but small enough that ~2 batches in
/// flight fit in a few GB of RAM.
pub const DEFAULT_ALIGNER_CHUNK_SIZE: u64 = 150_000_000;

/// Per-stage aligner tuning knobs.
///
/// Annotated with `#[multi_options("aligner", "Aligner Options")]` so
/// `runall` exposes each field as `--aligner::<flag>` via the generated
/// `MultiAlignerOptions` companion struct. Mutual-exclusion between
/// `preset` and `command` is enforced inside [`Self::resolve`] rather
/// than at clap-parse time — the macro doesn't rewrite clap's
/// `conflicts_with` field-name references when the prefixed identifiers
/// shift, so we validate logically after `MultiAlignerOptions::validate`.
///
/// Field shapes:
/// - `preset` — `Option<AlignerPreset>` so clap parses it as
///   `--aligner::preset {bwa-mem3|bwa}` (no default; if unset and
///   `command` is also unset, [`Self::resolve`] errors).
/// - `command` — `Option<String>` for the free-form mode.
/// - `threads` — `Option<usize>` so [`Self::resolve`] can default it
///   from `std::thread::available_parallelism()` at runtime (the macro
///   can't represent "all cores" as a const default).
/// - `chunk_size` — `u64` with a const default; drives the `-K` flag
///   in preset mode and the Step-1 batch size in both modes.
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
        Self { preset: None, command: None, threads: None, chunk_size: DEFAULT_ALIGNER_CHUNK_SIZE }
    }
}

/// Result of [`AlignerOptions::resolve`] — a ready-to-spawn aligner
/// invocation plus the parameters the downstream chain needs.
///
/// `pub(crate)` because only the runall AAM dispatch consumes it;
/// promote to `pub` if a cross-crate caller materializes.
///
/// Fields are flagged `#[allow(dead_code)]` because the AAM validation path
/// currently only validates them (via `_resolved` in
/// `validate_align_and_merge`); the wired chain sources its parameters
/// independently.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct ResolvedAligner {
    /// The shell command to spawn (already substituted and validated).
    /// Passed directly to [`AlignerProcess::spawn`].
    pub command: String,
    /// Bases per batch. Step 1 of the AAM chain uses this to size
    /// outgoing batches so the aligner sees full `-K` chunks.
    pub chunk_size: u64,
    /// Resolved thread count (default-substituted for preset mode;
    /// `None` if command mode let the user hardcode their own).
    pub threads: Option<usize>,
    /// Which mode produced this — preset (with which preset) or
    /// command. Used in info-logging.
    pub mode: ResolvedAlignerMode,
}

/// How a [`ResolvedAligner`] was produced. Same dead-code rationale
/// as [`ResolvedAligner`].
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub(crate) enum ResolvedAlignerMode {
    /// `--aligner::preset <preset>`. Carries the chosen preset so the
    /// caller can log it.
    Preset(AlignerPreset),
    /// `--aligner::command "..."`.
    Command,
}

impl AlignerOptions {
    /// Validate the option combination and produce a [`ResolvedAligner`]
    /// ready for [`AlignerProcess::spawn`].
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
    /// - Neither `--aligner::preset` nor `--aligner::command` was set.
    /// - Both were set (mutual exclusion violation).
    /// - Command mode + a preset-only flag (`--aligner-bin` or
    ///   `--aligner::threads`).
    /// - Command mode template missing `{ref}`.
    /// - Preset-mode index files / binary missing (delegated to
    ///   [`AlignerPreset::validate`]).
    pub(crate) fn resolve(
        self,
        reference: &Path,
        top_threads: usize,
        aligner_bin: Option<&Path>,
    ) -> Result<ResolvedAligner> {
        match (self.preset, self.command) {
            (None, None) => bail!(
                "--start-from align requires one of `--aligner::preset` \
                 (e.g. bwa-mem3, bwa) or `--aligner::command \"...\"`"
            ),
            (Some(_), Some(_)) => bail!(
                "--aligner::preset and --aligner::command are mutually exclusive; \
                 pass only one"
            ),
            (Some(preset), None) => {
                // Preset mode: validate indexes + binary, build argv.
                preset.validate(reference, aligner_bin)?;
                let threads = self.threads.unwrap_or(top_threads);
                let command =
                    preset.build_command(reference, threads, self.chunk_size, aligner_bin);
                Ok(ResolvedAligner {
                    command,
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
                // `--aligner::chunk-size` is NOT preset-only — Step 1
                // of the AAM chain uses it for BAM-record batching
                // regardless of mode, and command-mode users may
                // legitimately want to align it with their `-K` choice
                // in the template.
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
                    command,
                    chunk_size: self.chunk_size,
                    threads: None,
                    mode: ResolvedAlignerMode::Command,
                })
            }
        }
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

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    use super::*;

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

    /// Spawn a command that writes to both stdout and stderr; verify
    /// both are captured.
    #[test]
    fn test_stderr_capture() {
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
    /// `AlignAndMergeStep::Drop` safely issue its follow-up `wait()` without
    /// risking a teardown hang.
    #[test]
    fn test_kill_reports_reaped_for_killable_child() {
        let mut proc = AlignerProcess::spawn("sleep 30", 10).expect("spawn should succeed");
        assert!(proc.kill(), "a killable child must be reaped within the deadline");
        // A second kill on an already-gone child still reports reaped, never
        // a spurious timeout.
        assert!(proc.kill(), "killing an already-reaped child must still report reaped");
    }

    /// `bwa mem` command string with default-PATH binary.
    #[test]
    fn test_bwa_mem_command() {
        let reference = Path::new("/ref/genome.fa");
        let cmd = AlignerPreset::Bwa.build_command(reference, 8, 10_000_000, None);
        assert_eq!(cmd, "bwa mem -p -K 10000000 -t 8 /ref/genome.fa /dev/stdin");
    }

    /// `bwa-mem3` command string with default-PATH binary.
    #[test]
    fn test_bwa_mem3_command() {
        let reference = Path::new("/ref/genome.fa");
        let cmd = AlignerPreset::BwaMem3.build_command(reference, 4, 5_000_000, None);
        assert_eq!(cmd, "bwa-mem3 mem -p -K 5000000 -t 4 /ref/genome.fa /dev/stdin");
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
        };
        let resolved = opts.resolve(&ref_path, 4, None).unwrap();
        assert!(resolved.command.contains(&ref_path.display().to_string()));
        assert!(resolved.command.contains("-t 4"));
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
}
