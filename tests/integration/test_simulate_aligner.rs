//! Integration test for `fgumi simulate aligner`.
//!
//! Exercises the real binary as a subprocess through OS pipes — the exact shape
//! `fgumi runall --aligner::command` uses it in. The tool is format-agnostic (it
//! byte-copies the replay file), so we use arbitrary bytes as the "replay BAM"
//! and arbitrary bytes as the FASTQ on stdin. Both are sized larger than a pipe
//! buffer so the concurrent drain/emit is genuinely exercised: a single-threaded
//! lockstep replay would deadlock here.

use std::io::{Read, Write};
use std::process::{Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

use tempfile::TempDir;

/// Spawn `fgumi simulate aligner`, feed `stdin_bytes`, and return its stdout.
///
/// Writes stdin on a worker thread and reads stdout on the main thread so the
/// two directions cannot deadlock the test itself. A watchdog kills the child if
/// it has not exited within `timeout`, so a regression surfaces as a failed
/// assertion rather than a hung test.
fn run_simulate_aligner(
    replay: &std::path::Path,
    stdin_bytes: Vec<u8>,
    extra_args: &[&str],
) -> Vec<u8> {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.args(["simulate", "aligner", "--replay-bam"])
        .arg(replay)
        .args(extra_args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null());
    let mut child = cmd.spawn().expect("spawn fgumi simulate aligner");

    let mut child_stdin = child.stdin.take().expect("child stdin");
    let writer = thread::spawn(move || {
        // Ignore write errors: if the child exits early the write end breaks,
        // which is a legitimate (BrokenPipe) outcome, not a test failure.
        let _ = child_stdin.write_all(&stdin_bytes);
        drop(child_stdin); // EOF so the child's drainer completes
    });

    let mut stdout = child.stdout.take().expect("child stdout");
    let mut out = Vec::new();
    stdout.read_to_end(&mut out).expect("read child stdout");

    // Watchdog: the child should exit promptly once both pipes are serviced.
    let deadline = Instant::now() + Duration::from_secs(30);
    loop {
        if let Some(status) = child.try_wait().expect("try_wait") {
            assert!(status.success(), "simulate aligner exited non-zero: {status:?}");
            break;
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            panic!("simulate aligner did not exit within 30s (deadlock?)");
        }
        thread::sleep(Duration::from_millis(20));
    }
    writer.join().expect("writer thread");
    out
}

#[test]
fn replays_verbatim_while_draining_stdin_as_a_subprocess() {
    let tmp = TempDir::new().unwrap();
    let replay_path = tmp.path().join("replay.bam");

    // Replay + stdin both > a 64 KiB pipe buffer, so both directions block until
    // serviced — only the concurrent drain/emit keeps the subprocess live.
    let replay_bytes: Vec<u8> = (0..(2u32 * 1024 * 1024)).map(|i| (i % 251) as u8).collect();
    std::fs::write(&replay_path, &replay_bytes).unwrap();
    let stdin_bytes = vec![b'@'; 2 * 1024 * 1024];

    // Pass extra positional tokens like a `--aligner::command` template would
    // substitute for `{ref}` and `{threads}` — they must be accepted and ignored.
    let out = run_simulate_aligner(&replay_path, stdin_bytes, &["/fake/reference.fa", "8"]);

    assert_eq!(out, replay_bytes, "stdout must be the replay file, byte-for-byte");
}

#[test]
fn errors_when_replay_bam_is_missing() {
    let tmp = TempDir::new().unwrap();
    let missing = tmp.path().join("does-not-exist.bam");

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["simulate", "aligner", "--replay-bam"])
        .arg(&missing)
        .arg("/fake/reference.fa")
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("spawn fgumi simulate aligner");

    assert!(!status.success(), "missing replay BAM must be a non-zero exit");
}
