//! Integration test for `fgumi runall`'s Ctrl-C (SIGINT) handling.
//!
//! Spawns `fgumi runall` as a subprocess on a fixture that is large enough to
//! keep the pipeline running for at least a couple hundred milliseconds,
//! sends SIGINT, and asserts that:
//!   1. The process exits non-zero.
//!   2. No final output BAM is left behind.
//!   3. No `.tmp` partial output is left behind (cleaned up by the
//!      `AtomicOutputGuard` drop).
//!   4. The process terminates promptly (within a generous bound).
//!
//! Unix-only: uses `libc::kill` to deliver SIGINT to a specific PID. Windows
//! has different signal semantics; a Windows variant can be added later.

#![cfg(unix)]

use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;

use crate::helpers::bam_generator::create_minimal_header;

/// Build a template-coordinate-sorted BAM with many position groups so that
/// `runall --start-from group` has enough work to still be running when we
/// deliver SIGINT a couple hundred milliseconds after spawn.
///
/// The fixture contains `num_groups` UMI families spread across `num_groups`
/// distinct reference positions (position groups of size 8 each). This is
/// substantially larger than the tiny canonical fixture in
/// `grouped_bam_fixture.rs` (which has only 12 records total).
#[allow(clippy::cast_possible_wrap, clippy::cast_possible_truncation)]
fn build_large_grouped_bam(path: &Path, num_groups: usize) {
    let header = create_minimal_header("chr_test", 1_000_000);
    let file = std::fs::File::create(path).expect("create fixture BAM");
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).expect("write header");

    let sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bp
    let qualities = [30u8; 32];

    for g in 0..num_groups {
        // Spread position groups by 100 bp so each group is distinct. 8 reads
        // per position group; all reads share a unique RX and MI per group.
        let pos = 100 + (g * 100);
        let umi = format!("{:08}", g % 100_000_000);
        let mi = (g as u32 + 1).to_string();
        for read_idx in 0..8 {
            let name = format!("grp{g}_r{read_idx}");
            let rec = RecordBuilder::new()
                .name(&name)
                .sequence(sequence)
                .qualities(&qualities)
                .reference_sequence_id(0)
                .alignment_start(pos)
                .mapping_quality(60)
                .tag("RX", umi.as_str())
                .tag("MI", mi.as_str())
                .build();
            writer.write_alignment_record(&header, &rec).expect("write record");
        }
    }

    writer.try_finish().expect("finish BAM");
}

/// Wait up to `timeout` for `child` to exit. Returns the exit status if it
/// exits in time, otherwise kills the child (SIGKILL) and panics.
fn wait_for_exit(child: &mut std::process::Child, timeout: Duration) -> std::process::ExitStatus {
    let start = Instant::now();
    loop {
        match child.try_wait() {
            Ok(Some(status)) => return status,
            Ok(None) if start.elapsed() > timeout => {
                let _ = child.kill();
                panic!("child did not exit within {timeout:?} after SIGINT; forcibly killed");
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(50)),
            Err(e) => panic!("wait error: {e}"),
        }
    }
}

/// Deliver SIGINT to the given PID. Isolated `#[allow(unsafe_code)]` so the
/// rest of the test file stays in the safe subset enforced by the crate-level
/// `#![deny(unsafe_code)]`.
#[allow(unsafe_code)]
fn send_sigint(pid: u32) {
    // `pid` came from `Child::id()` (u32) but POSIX kill takes pid_t (i32).
    // The child's PID is always a small positive integer in practice, but
    // `cast_signed()` is the lint-clean conversion.
    let pid_i32 = pid.cast_signed();
    // SAFETY: `libc::kill` is a thin wrapper over the POSIX kill(2) syscall.
    // We pass a PID we just obtained from `Child::id()` (a live child we
    // spawned) and SIGINT, a well-defined signal. Failure here just means
    // the child already exited, which the test will detect via try_wait.
    unsafe {
        let _ = libc::kill(pid_i32, libc::SIGINT);
    }
}

#[test]
fn runall_ctrl_c_leaves_no_corrupt_output() {
    let tmp = tempfile::tempdir().expect("tempdir");
    let input = tmp.path().join("input.bam");
    let output = tmp.path().join("output.bam");
    let tmp_output = output.with_extension("bam.tmp");

    // 100k position groups × 8 reads = 800k records. Enough work that the
    // pipeline is still running when SIGINT lands ~200ms after spawn.
    build_large_grouped_bam(&input, 100_000);

    let mut child = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        // Ensure tracing::warn! output is visible on stderr so we can assert on
        // the cancel message. RUST_LOG=info captures both the handler's
        // "cancelling pipeline" warn and the upstream pipeline info logs, so
        // assertion failures include plenty of diagnostic context.
        .env("RUST_LOG", "info")
        .arg("runall")
        .arg("--start-from")
        .arg("group")
        .arg("--threads")
        .arg("1")
        .arg("--filter::min-reads")
        .arg("1")
        .arg("--filter::min-base-quality")
        .arg("0")
        .arg("--input")
        .arg(&input)
        .arg("--output")
        .arg(&output)
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn fgumi");

    let pid = child.id();

    // Give runall time to start (clap parse, validate, install signal
    // handler, open BAM, prime stages). 1000ms is conservative: well above
    // fgumi's cold-start cost even in debug builds, yet well below the time
    // to process 800k records through the group/consensus/filter pipeline.
    std::thread::sleep(Duration::from_millis(1000));

    send_sigint(pid);

    let status = wait_for_exit(&mut child, Duration::from_secs(30));

    // Collect stderr for diagnostics and the cancel-message assertion.
    let mut stderr_bytes = Vec::new();
    if let Some(mut err) = child.stderr.take() {
        use std::io::Read;
        let _ = err.read_to_end(&mut stderr_bytes);
    }
    let stderr = String::from_utf8_lossy(&stderr_bytes);

    assert!(
        !status.success(),
        "runall should exit non-zero on SIGINT; got {status:?}. stderr:\n{stderr}",
    );

    assert!(
        !output.exists(),
        "final output BAM must not exist after cancel: {}. stderr:\n{stderr}",
        output.display()
    );

    assert!(
        !tmp_output.exists(),
        "tmp output BAM must be cleaned up after cancel: {}. stderr:\n{stderr}",
        tmp_output.display()
    );

    // Sanity: the cancel path should log something recognizable. We tolerate
    // either our explicit "cancelling pipeline" log from the signal handler
    // or the "pipeline cancelled by signal" error chain.
    let lower = stderr.to_ascii_lowercase();
    assert!(
        lower.contains("cancel") || lower.contains("interrupt"),
        "expected stderr to mention cancel/interrupt; got:\n{stderr}"
    );
}
