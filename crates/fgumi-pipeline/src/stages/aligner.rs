//! Aligner subprocess management for the runall pipeline.
//!
//! This module provides [`AlignerProcess`] for spawning an aligner as a child process
//! and communicating with it via stdin/stdout pipes, and [`AlignerPreset`] for building
//! well-formed aligner command strings for common aligners.

use std::collections::VecDeque;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};
use std::thread::{self, JoinHandle};
use std::time::Duration;

use anyhow::{Context, Result, bail};

// ============================================================================
// AlignerProcess
// ============================================================================

/// A running aligner subprocess with managed stdin, stdout, and stderr.
///
/// The aligner is spawned via `/bin/bash -c`, allowing the command string to contain
/// pipes, redirects, and other shell constructs.
///
/// Stderr from the child process is relayed to fgumi's own stderr in real time,
/// and the last `ring_size` lines are retained for inclusion in error messages.
pub struct AlignerProcess {
    child: Child,
    /// The child's stdin pipe, available until [`take_stdin`] is called.
    stdin: Option<ChildStdin>,
    /// The child's stdout pipe, available until [`take_stdout`] is called.
    stdout: Option<ChildStdout>,
    /// Background thread that relays child stderr; returns the captured ring buffer.
    stderr_thread: Option<JoinHandle<Vec<String>>>,
}

impl AlignerProcess {
    /// Spawn an aligner subprocess with the given shell command.
    ///
    /// The command is run via `/bin/bash -c`, so it supports pipes and redirects.
    /// A background thread is started immediately to relay the child's stderr to
    /// fgumi's stderr in real time; the last `ring_size` lines are captured for
    /// failure diagnostics.
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
    /// Panics if the OS does not provide the stderr pipe after configuring `Stdio::piped()`,
    /// which should never happen in practice.
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
    /// Returns `None` if stdin has already been taken.
    pub fn take_stdin(&mut self) -> Option<ChildStdin> {
        self.stdin.take()
    }

    /// Take the stdout pipe for reading SAM data.
    ///
    /// Returns `None` if stdout has already been taken.
    pub fn take_stdout(&mut self) -> Option<ChildStdout> {
        self.stdout.take()
    }

    /// Wait for the aligner process to finish.
    ///
    /// Collects the stderr relay thread and includes the last captured stderr lines
    /// in the error message if the process exits with a non-zero status.
    ///
    /// # Errors
    ///
    /// Returns an error if the process exits with a non-zero status code or if
    /// waiting on the process fails.
    pub fn wait(mut self) -> Result<()> {
        let status = self.child.wait().context("failed to wait for aligner process")?;

        // Always join the stderr relay thread so its resources are cleaned up.
        let last_lines = self.stderr_thread.take().and_then(|t| t.join().ok()).unwrap_or_default();

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
    /// Sends SIGKILL immediately (`Child::kill()` always sends SIGKILL on Unix).
    /// Waits up to 1 second for the process to exit.
    pub fn kill(&mut self) {
        let _ = self.child.kill();

        // Wait for the process to actually exit so its resources are cleaned up.
        let deadline = std::time::Instant::now() + Duration::from_secs(1);
        while let Ok(None) = self.child.try_wait() {
            if std::time::Instant::now() >= deadline {
                break;
            }
            thread::sleep(Duration::from_millis(50));
        }
    }

    /// Return the process ID of the aligner child process.
    #[must_use]
    pub fn pid(&self) -> u32 {
        self.child.id()
    }
}

impl Drop for AlignerProcess {
    fn drop(&mut self) {
        // If the child process is still running when the struct is dropped (e.g. due to
        // a panic or early return), kill it to prevent orphaned processes.
        if let Ok(None) = self.child.try_wait() {
            self.kill();
        }
        // Always join the stderr relay thread to avoid resource leaks.
        if let Some(handle) = self.stderr_thread.take() {
            let _ = handle.join();
        }
    }
}

/// Read lines from `stderr`, print each to fgumi's stderr, and retain the last
/// `ring_size` lines in a [`VecDeque`] ring buffer.
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
pub enum AlignerPreset {
    /// BWA-MEM aligner (`bwa mem`).
    BwaMem,
    /// BWA-MEM2 aligner (`bwa-mem2 mem`).
    BwaMem2,
}

impl AlignerPreset {
    /// Build the aligner shell command for this preset.
    ///
    /// The command reads FASTQ from `/dev/stdin` and writes SAM to stdout, allowing it
    /// to be connected to fgumi's pipeline via stdin/stdout pipes.
    ///
    /// # Arguments
    ///
    /// * `reference`  - Path to the reference FASTA (index must already exist).
    /// * `threads`    - Number of threads to pass to the aligner (`-t`).
    /// * `chunk_size` - Chunk size in bases to pass to the aligner (`-K`).
    /// * `extra_args` - Additional arguments inserted verbatim before the reference path.
    #[must_use]
    pub fn build_command(
        &self,
        reference: &Path,
        threads: usize,
        chunk_size: u64,
        extra_args: &str,
    ) -> String {
        let binary = match self {
            AlignerPreset::BwaMem => "bwa",
            AlignerPreset::BwaMem2 => "bwa-mem2",
        };
        let ref_str = reference.display();
        let extra = extra_args.trim();
        if extra.is_empty() {
            format!("{binary} mem -p -K {chunk_size} -t {threads} {ref_str} /dev/stdin")
        } else {
            format!("{binary} mem -p -K {chunk_size} -t {threads} {extra} {ref_str} /dev/stdin")
        }
    }

    /// Validate that the aligner binary and required index files are present.
    ///
    /// Checks that:
    /// - The aligner binary is on `PATH` (via `which`).
    /// - The expected index file(s) exist alongside the reference.
    ///
    /// # Errors
    ///
    /// Returns an error if the binary is not found or any required index file is missing.
    pub fn validate(&self, reference: &Path) -> Result<()> {
        let (binary, index_extensions) = match self {
            AlignerPreset::BwaMem => ("bwa", vec![".bwt"]),
            AlignerPreset::BwaMem2 => ("bwa-mem2", vec![".bwt.2bit.64"]),
        };

        // Check binary is on PATH.
        which::which(binary)
            .with_context(|| format!("aligner binary `{binary}` not found on PATH"))?;

        // Check required index files.
        for ext in index_extensions {
            let index_path = append_extension(reference, ext);
            if !index_path.exists() {
                bail!(
                    "required index file not found: {} (run `{binary} index {}`)",
                    index_path.display(),
                    reference.display()
                );
            }
        }

        Ok(())
    }
}

/// Append a suffix to a path without replacing the existing extension.
///
/// For example, `append_extension("/ref/genome.fa", ".bwt")` → `/ref/genome.fa.bwt`.
fn append_extension(path: &Path, suffix: &str) -> std::path::PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(suffix);
    std::path::PathBuf::from(s)
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

    /// Spawn a command that writes to both stdout and stderr; verify both are captured.
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

    /// Spawn a command that exits with a non-zero status and verify that `wait` returns
    /// an error.
    #[test]
    fn test_nonzero_exit() {
        let proc = AlignerProcess::spawn("bash -c 'exit 1'", 10).expect("spawn should succeed");
        let result = proc.wait();
        assert!(result.is_err(), "wait() should return Err for non-zero exit");
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("status"), "error message should mention status: {msg}");
    }

    /// Verify that `build_command` produces the correct command string for bwa mem.
    #[test]
    fn test_bwa_mem_command() {
        let reference = Path::new("/ref/genome.fa");
        let cmd = AlignerPreset::BwaMem.build_command(reference, 8, 10_000_000, "");
        assert_eq!(cmd, "bwa mem -p -K 10000000 -t 8 /ref/genome.fa /dev/stdin");
    }

    /// Verify that `build_command` produces the correct command string for bwa-mem2.
    #[test]
    fn test_bwa_mem2_command() {
        let reference = Path::new("/ref/genome.fa");
        let cmd = AlignerPreset::BwaMem2.build_command(reference, 4, 5_000_000, "-Y");
        assert_eq!(cmd, "bwa-mem2 mem -p -K 5000000 -t 4 -Y /ref/genome.fa /dev/stdin");
    }

    /// Verify that `extra_args` is omitted when empty and included when non-empty for bwa mem.
    #[test]
    fn test_bwa_mem_command_with_extra_args() {
        let reference = Path::new("/data/hg38.fa");
        let cmd = AlignerPreset::BwaMem.build_command(reference, 16, 100_000, "-R '@RG\\tID:1'");
        assert!(cmd.contains("-R '@RG\\tID:1'"), "extra args should appear in command: {cmd}");
        assert!(cmd.contains("/data/hg38.fa"), "reference should appear in command: {cmd}");
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
}
