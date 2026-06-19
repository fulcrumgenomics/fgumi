//! Simulate a streaming aligner for benchmarking `fgumi runall`'s align stage.
//!
//! `fgumi simulate aligner` is a *fake* aligner subprocess: it reads the
//! interleaved FASTQ that `fgumi runall`'s `AlignAndMerge` writes to its stdin,
//! discards it, and streams a pre-captured aligned BAM back on its stdout. This
//! lets the benchmark measure fgumi's own pipeline overhead without a real
//! aligner's compute dominating every run.
//!
//! It exists to replace the fgumi-benchmarks `replay-aligner.sh` shell script,
//! which read stdin and wrote stdout in lockstep on a single thread and thereby
//! deadlocked `AlignAndMerge`: the writer fills its in-flight gate feeding FASTQ
//! and blocks; the reader stalls waiting for aligner output; the shell replay
//! won't emit because it is mid-drain wanting more stdin the blocked writer
//! cannot supply.
//!
//! The fix is to read and write **concurrently** — exactly what a real streaming
//! aligner does. Two workers that run concurrently and never block each other
//! while streaming:
//!
//!   - a *drainer* that reads stdin to EOF and discards it, and
//!   - an *emitter* that copies the replay BAM verbatim to stdout.
//!
//! Because the replay file already *is* BGZF/BAM on disk — exactly the byte
//! stream `AlignAndMerge`'s reader expects from a BAM-emitting aligner — the
//! emitter is a raw byte copy: the BAM header is emitted first, and every
//! record (including SECONDARY/SUPPLEMENTARY) is passed through in stored order,
//! all for free. `AlignAndMerge` regroups records into templates by queryname,
//! so the replay BAM must be in input (queryname-grouped) order, which is how
//! `bwa mem` (and the benchmark's `bwa_align` rule) emit it.

use crate::commands::command::Command;
use anyhow::{Context, Result};
use clap::Parser;
use log::info;
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::PathBuf;
use std::thread;

/// Simulate a streaming aligner by replaying a pre-captured aligned BAM.
#[derive(Parser, Debug)]
#[command(
    name = "aligner",
    about = "Simulate a streaming aligner by replaying a pre-captured aligned BAM",
    long_about = r#"
Fake aligner subprocess for `fgumi runall --aligner::command`.

Reads the interleaved FASTQ on stdin (discarding it) while concurrently
streaming a pre-captured aligned BAM (`--replay-bam`) to stdout. Reading and
writing run on independent threads so the process never deadlocks the caller's
align stage, the way a single-threaded lockstep replay does.

The replay BAM is streamed verbatim (header first; every record, including
secondary/supplementary, in stored order), so it must be in input /
queryname-grouped order (as `bwa mem` emits it). A positional reference path is
accepted and ignored for drop-in compatibility with real-aligner command
templates that require a `{ref}` token.
"#
)]
pub struct Aligner {
    /// Pre-captured aligned BAM to stream back on stdout.
    #[arg(long = "replay-bam")]
    pub replay_bam: PathBuf,

    /// Positional tokens substituted into the aligner command template
    /// (`{ref}`, `{threads}`, …) — accepted and ignored, so the command is a
    /// drop-in for real-aligner templates regardless of which tokens they fill.
    #[arg(value_name = "IGNORED", trailing_var_arg = true, allow_hyphen_values = true)]
    pub ignored: Vec<String>,
}

impl Command for Aligner {
    fn execute(&self, _command_line: &str) -> Result<()> {
        let replay = File::open(&self.replay_bam)
            .with_context(|| format!("opening replay BAM {}", self.replay_bam.display()))?;
        info!("simulate aligner: replaying {}", self.replay_bam.display());

        // Use the `Stdin`/`Stdout` handles (which are `Send`) rather than their
        // locks (`StdinLock` holds a non-`Send` guard) so the drainer can own
        // stdin on its own thread.
        simulate_aligner(io::stdin(), replay, io::stdout())?;
        Ok(())
    }
}

/// Concurrently drain `stdin` (discarding it) and copy `replay` to `stdout`.
///
/// The drainer and emitter run on separate threads and never wait on each other,
/// so the function cannot deadlock a caller that backpressures one direction
/// while feeding the other. A closed stdout read end (the caller went away) is
/// treated as a clean shutdown — the same way a real aligner exits on SIGPIPE —
/// rather than an error.
///
/// Both arguments are consumed (moved into the worker that owns them) so the
/// underlying handles close when the function returns.
fn simulate_aligner<I, B, O>(mut stdin: I, mut replay: B, mut stdout: O) -> io::Result<()>
where
    I: Read + Send,
    B: Read,
    O: Write,
{
    /// A closed pipe (the caller went away) is a clean shutdown for an aligner,
    /// not an error — map `BrokenPipe` to `Ok`.
    fn allow_broken_pipe(r: io::Result<()>) -> io::Result<()> {
        match r {
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
            other => other,
        }
    }

    thread::scope(|scope| {
        // Drainer: read all of stdin and discard. Runs concurrently with the
        // emitter so the caller's FASTQ writer is never blocked by our output.
        let drainer = scope.spawn(move || io::copy(&mut stdin, &mut io::sink()).map(|_| ()));

        // Emitter (this thread): stream the replay BAM verbatim to stdout.
        let emit =
            allow_broken_pipe(io::copy(&mut replay, &mut stdout).and_then(|_| stdout.flush()));

        // Close stdout so the caller's reader sees EOF and can wind down — even
        // when the emitter errored partway (e.g. a corrupt replay BAM). Without
        // this, the drainer join below would block forever on a caller that is
        // still feeding FASTQ and waiting on output we will never produce.
        drop(stdout);

        // Wait for the caller to finish feeding FASTQ before we exit, so we never
        // SIGPIPE its writer by closing stdin early.
        let drain = allow_broken_pipe(drainer.join().expect("simulate aligner drainer panicked"));

        emit.and(drain)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::io::{Cursor, Read, Write};
    use std::time::{Duration, Instant};

    /// A reader that yields `remaining` bytes then fails — models a corrupt or
    /// truncated replay BAM whose read errors partway through.
    struct ErrAfter {
        remaining: usize,
    }
    impl Read for ErrAfter {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.remaining == 0 {
                return Err(io::Error::other("simulated replay read error"));
            }
            let n = buf.len().min(self.remaining);
            buf[..n].fill(0xEE);
            self.remaining -= n;
            Ok(n)
        }
    }

    #[test]
    fn accepts_and_ignores_extra_positional_tokens() {
        // A `--aligner::command` template may substitute `{ref}` and `{threads}`,
        // producing extra positionals; they must be accepted and ignored.
        let a = Aligner::try_parse_from(["aligner", "--replay-bam", "x.bam", "/ref.fa", "8"])
            .expect("parse with ref + threads tokens");
        assert_eq!(a.replay_bam, PathBuf::from("x.bam"));
        assert_eq!(a.ignored, vec!["/ref.fa".to_string(), "8".to_string()]);
    }

    #[test]
    fn emit_error_surfaces_as_err_without_hanging() {
        // When the emitter errors partway (corrupt replay) against a live caller
        // that is still feeding stdin, the function must close stdout (so the
        // caller winds down) and return the error — not hang on the drainer join.
        let (stdin_r, mut stdin_w) = io::pipe().expect("stdin pipe");
        let (mut stdout_r, stdout_w) = io::pipe().expect("stdout pipe");
        let replay = ErrAfter { remaining: 1024 };

        let sim = thread::spawn(move || simulate_aligner(stdin_r, replay, stdout_w));

        // Caller mimics AlignAndMerge: feed some FASTQ, then read stdout to EOF
        // (released by the emitter's drop(stdout) on error), then close stdin.
        let caller = thread::spawn(move || {
            let _ = stdin_w.write_all(&[0x40u8; 4096]);
            let mut got = Vec::new();
            let _ = stdout_r.read_to_end(&mut got);
            drop(stdin_w);
        });

        let deadline = Instant::now() + Duration::from_secs(10);
        while !(sim.is_finished() && caller.is_finished()) {
            assert!(Instant::now() < deadline, "simulate aligner hung on emit error");
            thread::sleep(Duration::from_millis(20));
        }
        let res = sim.join().expect("sim thread");
        caller.join().expect("caller thread");
        assert!(res.is_err(), "emit error must surface as Err, got {res:?}");
    }

    #[test]
    fn streams_replay_verbatim_and_drains_stdin() {
        let stdin = Cursor::new(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nIIII\n".to_vec());
        let replay = Cursor::new(b"\x1f\x8b\x08replay-bam-bytes-verbatim".to_vec());
        let mut stdout: Vec<u8> = Vec::new();

        simulate_aligner(stdin, replay, &mut stdout).expect("simulate aligner");

        assert_eq!(stdout, b"\x1f\x8b\x08replay-bam-bytes-verbatim");
    }

    #[test]
    fn does_not_deadlock_against_a_single_threaded_write_then_read_caller() {
        // Model the cycle that wedged the shell replay. `AlignAndMerge` is, in
        // effect, a single-threaded caller: it writes FASTQ to the aligner's
        // stdin and reads alignments from its stdout, and it cannot finish
        // writing more input until its reader drains output (its in-flight gate),
        // nor read output the aligner hasn't produced. We reproduce that coupling
        // with ONE caller thread that writes all of stdin, THEN reads all of
        // stdout, against bounded OS pipes — and inputs/outputs both larger than a
        // pipe buffer so each direction blocks until serviced.
        //
        // A correct aligner drains stdin and emits stdout *concurrently*, so the
        // caller's write phase always makes progress (drainer) and its read phase
        // always finds output (emitter). A lockstep aligner that emits before it
        // drains deadlocks here: the caller blocks writing stdin (pipe full, not
        // drained) while the aligner blocks writing stdout (pipe full, not read).
        let (stdin_r, mut stdin_w) = io::pipe().expect("stdin pipe");
        let (mut stdout_r, stdout_w) = io::pipe().expect("stdout pipe");

        // Both larger than any pipe buffer, forcing both directions to block.
        let stdin_bytes = vec![0x40u8; 4 * 1024 * 1024];
        let replay_bytes = vec![0xABu8; 4 * 1024 * 1024];
        let replay = Cursor::new(replay_bytes.clone());

        let sim = thread::spawn(move || simulate_aligner(stdin_r, replay, stdout_w));

        // Caller on a SINGLE thread: write all of stdin, THEN read all of stdout
        // — sequentially, so the two directions are not independently serviced.
        // This is the coupling that gives the test teeth.
        let caller = thread::spawn(move || {
            stdin_w.write_all(&stdin_bytes).expect("write stdin");
            drop(stdin_w); // EOF so the drainer completes
            let mut got = Vec::new();
            stdout_r.read_to_end(&mut got).expect("read stdout");
            got
        });

        // Watchdog on main: a deadlock shows up as the worker threads never
        // finishing. (We can't join-with-timeout in std, so poll is_finished.)
        let deadline = Instant::now() + Duration::from_secs(10);
        while !(sim.is_finished() && caller.is_finished()) {
            assert!(Instant::now() < deadline, "simulate aligner deadlocked");
            thread::sleep(Duration::from_millis(20));
        }
        sim.join().expect("sim thread").expect("simulate aligner result");
        let got = caller.join().expect("caller thread");

        assert_eq!(got, replay_bytes, "stdout must be the verbatim replay");
    }

    #[test]
    fn broken_stdout_pipe_is_a_clean_exit() {
        // Caller drops its stdout read end immediately; the emitter's write
        // should surface as a clean shutdown, not an error.
        let (stdout_r, stdout_w) = io::pipe().expect("stdout pipe");
        drop(stdout_r); // reader gone before we write

        let stdin = Cursor::new(b"@r1\nACGT\n+\nIIII\n".to_vec());
        let replay = Cursor::new(vec![0xCDu8; 1024 * 1024]);

        let res = simulate_aligner(stdin, replay, stdout_w);
        assert!(res.is_ok(), "BrokenPipe on stdout must be a clean exit: {res:?}");
    }
}
