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
    ///
    /// Because this field is `trailing_var_arg` + `allow_hyphen_values`, once the
    /// first positional token is seen *every* remaining token — including
    /// flag-like ones such as `-t 8` or even `--replay-bam` itself — is swallowed
    /// into `ignored` rather than parsed. Therefore `--replay-bam` MUST be
    /// supplied before any positional token; placing it after a positional makes
    /// it disappear into `ignored` and the command fails with a missing
    /// `--replay-bam` error.
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
/// The drainer and emitter run concurrently and never wait on each other, so the
/// function cannot deadlock a caller that backpressures one direction while
/// feeding the other.
///
/// The drainer runs on a *detached* thread, joined only on the full-emit success
/// path. On the happy path the full replay BAM (with its in-band BGZF EOF marker)
/// reaches the caller, whose BAM reader stops on its own; we then join the drainer
/// so the caller finishes feeding FASTQ before we exit and we never SIGPIPE its
/// writer.
///
/// The other two outcomes return *immediately without joining the drainer*,
/// because joining it could block forever on stdin that never reaches EOF:
///
///   - A `BrokenPipe` on stdout means the caller closed its read end (it went
///     away) — a clean shutdown, the same way a real aligner exits on SIGPIPE. We
///     return `Ok(())`: there is no reader left to finish feeding FASTQ for, so
///     there is nothing to wait on.
///   - On an emit error the replay stream is truncated (no EOF marker), so the
///     caller's reader would block waiting for bytes that never come — and if it
///     is coupled to its FASTQ writer by an in-flight gate (as `AlignAndMerge`
///     is), the writer wedges too and never closes our stdin. We return the error;
///     process teardown then closes our fd 1, the caller's reader sees EOF and the
///     whole pipeline winds down. We cannot rely on dropping `stdout` to deliver
///     that EOF — dropping `io::Stdout` does not close fd 1 (it is intentionally
///     kept open for the process lifetime).
fn simulate_aligner<I, B, O>(mut stdin: I, mut replay: B, mut stdout: O) -> io::Result<()>
where
    I: Read + Send + 'static,
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

    // Drainer: read all of stdin and discard, on a detached thread that runs
    // concurrently with the emitter so the caller's FASTQ writer is never blocked
    // by our output. Detached (not scoped) so an emit error can return without
    // waiting on a caller that may never close stdin (see the function doc).
    let drainer = thread::spawn(move || io::copy(&mut stdin, &mut io::sink()).map(|_| ()));

    // Emitter (this thread): stream the replay BAM verbatim to stdout.
    //
    // A `BrokenPipe` means the caller closed its stdout read end (it went away):
    // there is no reader left to finish feeding FASTQ for, so we return a clean
    // shutdown *immediately* without joining the drainer — joining it could block
    // forever if the (departed) caller never closes our stdin. Only a full,
    // successful emit falls through to join the drainer below.
    match io::copy(&mut replay, &mut stdout).and_then(|_| stdout.flush()) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => return Ok(()),
        Err(e) => return Err(e),
    }

    // Emit succeeded: the caller has the complete BAM and its reader will stop on
    // the in-band EOF marker. Wait for the caller to finish feeding FASTQ before
    // we exit, so we never SIGPIPE its writer by dropping our stdin read end.
    allow_broken_pipe(drainer.join().expect("simulate aligner drainer panicked"))
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
    fn emit_error_returns_without_joining_a_caller_that_never_closes_stdin() {
        // Models the production deadlock: the emitter errors partway (corrupt
        // replay) while the caller is wedged — its BAM reader blocks on the
        // truncated stream and its in-flight gate keeps its FASTQ writer from ever
        // closing our stdin. The function MUST surface the error promptly without
        // joining the (detached) drainer, which would otherwise block forever on
        // stdin that never reaches EOF. We deliberately keep `stdin_w` open for the
        // whole test, so the drainer can never complete; a join-the-drainer
        // implementation hangs here.
        //
        // This is the case the old `drop(stdout)` could not cover with a real
        // `io::Stdout`: dropping it does not close fd 1, so a join-based version
        // never unblocked. A `PipeWriter` (which does close on drop) would mask
        // that, so this test instead asserts the no-join behavior directly.
        let (stdin_r, stdin_w) = io::pipe().expect("stdin pipe");
        // Hold the stdout read end open so the partial emit does not BrokenPipe;
        // the error must come from the replay read, not a closed reader.
        let (stdout_r, stdout_w) = io::pipe().expect("stdout pipe");
        let replay = ErrAfter { remaining: 1024 };

        let sim = thread::spawn(move || simulate_aligner(stdin_r, replay, stdout_w));

        let deadline = Instant::now() + Duration::from_secs(10);
        while !sim.is_finished() {
            assert!(
                Instant::now() < deadline,
                "simulate aligner hung on the drainer join after an emit error"
            );
            thread::sleep(Duration::from_millis(20));
        }
        let res = sim.join().expect("sim thread");
        assert!(res.is_err(), "emit error must surface as Err, got {res:?}");

        // Cleanup: closing stdin_w lets the detached drainer finish.
        drop(stdin_w);
        drop(stdout_r);
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

    #[test]
    fn broken_stdout_pipe_returns_without_joining_a_caller_that_never_closes_stdin() {
        // The caller's reader is gone (stdout read end dropped) so the emit hits
        // `BrokenPipe`, but the caller never closes our stdin. A correct aligner
        // must NOT wait on the drainer here — there is no reader left to finish
        // feeding FASTQ for, so joining the (detached) drainer would block forever
        // on stdin that never reaches EOF. We hold `stdin_w` open for the whole
        // test, so a join-on-BrokenPipe implementation hangs here.
        let (stdin_r, stdin_w) = io::pipe().expect("stdin pipe");
        let (stdout_r, stdout_w) = io::pipe().expect("stdout pipe");
        drop(stdout_r); // reader gone: the emit will BrokenPipe

        let replay = Cursor::new(vec![0xCDu8; 1024 * 1024]);
        let sim = thread::spawn(move || simulate_aligner(stdin_r, replay, stdout_w));

        let deadline = Instant::now() + Duration::from_secs(10);
        while !sim.is_finished() {
            assert!(
                Instant::now() < deadline,
                "simulate aligner hung on the drainer join after a BrokenPipe stdout"
            );
            thread::sleep(Duration::from_millis(20));
        }
        let res = sim.join().expect("sim thread");
        assert!(res.is_ok(), "BrokenPipe on stdout must be a clean exit: {res:?}");

        // Cleanup: closing stdin_w lets the detached drainer finish.
        drop(stdin_w);
    }

    #[test]
    fn replayed_bytes_are_a_valid_bam_the_readers_accept() {
        // End-to-end: the module's whole purpose is to feed a *BAM stream*
        // (header + records + in-band BGZF EOF) to a real BAM reader. Write a
        // small valid BAM to a temp file, replay it through `simulate_aligner`
        // capturing stdout, then decode the captured bytes with the project's
        // BAM reader and assert the header and record count round-trip. A
        // regression that dropped a flush or truncated the tail would pass the
        // opaque-byte tests above but fail here.
        use fgumi_bam_io::PipelineReaderOpts;
        use fgumi_raw_bam::RawRecord;
        use fgumi_raw_bam::testutil::make_bam_bytes;
        use noodles::sam::Header;

        const N_RECORDS: usize = 5;
        const HEADER_MARKER: &str = "simulate-aligner-e2e-marker";

        // A header with an identifiable comment so we can assert it survives.
        let header = Header::builder().add_comment(HEADER_MARKER).build();

        let replay_path = tempfile::NamedTempFile::new().expect("replay temp").into_temp_path();
        let mut writer =
            fgumi_bam_io::create_raw_bam_writer(&replay_path, &header, 1, 1).expect("writer");
        for i in 0..N_RECORDS {
            let name = format!("read_{i}");
            let record_bytes = make_bam_bytes(
                -1,  // tid (unmapped)
                -1,  // pos
                0x4, // flag (unmapped)
                name.as_bytes(),
                &[], // cigar_ops
                0,   // seq_len
                -1,  // mate_tid
                -1,  // mate_pos
                &[], // aux_data
            );
            writer.write_raw_record(&record_bytes).expect("write record");
        }
        writer.finish().expect("finish writer");

        // Replay the BAM through the aligner, capturing the verbatim stdout.
        let replay = File::open(&replay_path).expect("open replay");
        let stdin = Cursor::new(b"@r1\nACGT\n+\nIIII\n".to_vec());
        let mut captured: Vec<u8> = Vec::new();
        simulate_aligner(stdin, replay, &mut captured).expect("simulate aligner");

        // Decode the captured bytes with the project's BAM reader by routing
        // them back through a temp file (the reader factory takes a path).
        let out_path = tempfile::NamedTempFile::new().expect("out temp").into_temp_path();
        std::fs::write(&out_path, &captured).expect("write captured bytes");

        let (mut reader, got_header) = fgumi_bam_io::create_raw_bam_reader_with_opts(
            &out_path,
            1,
            PipelineReaderOpts::default(),
        )
        .expect("reader accepts replayed bytes");

        // Header round-trips: the marker comment survives the verbatim copy.
        let comments: Vec<String> = got_header.comments().iter().map(|c| c.to_string()).collect();
        assert!(
            comments.iter().any(|c| c == HEADER_MARKER),
            "header comment must round-trip, got comments: {comments:?}"
        );

        // Record count round-trips, and the reader stops cleanly on the in-band
        // BGZF EOF marker (read_record returns 0 rather than erroring).
        let mut record = RawRecord::default();
        let mut n = 0;
        while reader.read_record(&mut record).expect("read replayed record") > 0 {
            n += 1;
        }
        assert_eq!(n, N_RECORDS, "record count must round-trip through the replay");
    }

    #[test]
    fn replay_bam_flag_after_a_positional_is_swallowed_into_ignored() {
        // `ignored` is trailing_var_arg + allow_hyphen_values, so once the first
        // positional is seen every later token — including `--replay-bam` — is
        // captured into `ignored` rather than parsed. Placing `--replay-bam`
        // after a positional therefore makes it disappear and parsing fails with
        // the required-arg error. This pins the ordering constraint documented on
        // the `ignored` field.
        let err = Aligner::try_parse_from(["aligner", "/ref.fa", "--replay-bam", "x.bam"])
            .expect_err("--replay-bam after a positional must fail to parse");
        // clap reports the missing required `--replay-bam`.
        assert_eq!(err.kind(), clap::error::ErrorKind::MissingRequiredArgument);
    }
}
