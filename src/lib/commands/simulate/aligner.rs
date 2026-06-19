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
//! aligner does. Two independent workers that never wait on each other:
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

    /// Reference path (accepted and ignored; for `{ref}`-template compatibility).
    #[arg(value_name = "REF")]
    pub reference: Option<PathBuf>,
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

        // Wait for the caller to finish feeding FASTQ before we exit, so we never
        // SIGPIPE its writer by closing stdin early.
        let drain = allow_broken_pipe(drainer.join().expect("simulate aligner drainer panicked"));

        emit.and(drain)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Read, Write};
    use std::time::{Duration, Instant};

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
