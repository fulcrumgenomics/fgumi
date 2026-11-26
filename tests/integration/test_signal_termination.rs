//! Termination signals (SIGINT, SIGTERM, SIGHUP) cancel the running command,
//! let it clean up, and then kill fgumi *by that signal* (issue #1148).
//!
//! Every test spawns the real binary in its own process group, so a `killpg`
//! on that group delivers exactly what a terminal's Ctrl-C delivers to a
//! foreground job, and asserts on `ExitStatusExt::signal()`, never on a
//! shell-style numeric code.
#![cfg(unix)]

use std::fs::File;
use std::io::Write;
use std::os::unix::process::{CommandExt, ExitStatusExt};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::sync::mpsc;
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use rstest::rstest;
use rustix::process::{Pid, Signal};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
use crate::helpers::chain_spec::records;
use crate::helpers::cli::utf8_path as p;
use crate::test_realign_command::{process_is_running, write_consensus};

/// The notice fgumi logs when it accepts a termination signal.
const NOTICE: &str = "stopping and cleaning up";

#[derive(Debug, Clone, Copy)]
enum Sig {
    Int,
    Term,
    Hup,
}

impl Sig {
    fn signal(self) -> Signal {
        match self {
            Sig::Int => Signal::INT,
            Sig::Term => Signal::TERM,
            Sig::Hup => Signal::HUP,
        }
    }

    fn raw(self) -> i32 {
        self.signal().as_raw()
    }

    fn name(self) -> &'static str {
        match self {
            Sig::Int => "SIGINT",
            Sig::Term => "SIGTERM",
            Sig::Hup => "SIGHUP",
        }
    }
}

/// Whom a signal is sent to.
#[derive(Debug, Clone, Copy)]
enum Target {
    /// fgumi's pid alone (`kill`).
    Process,
    /// fgumi's whole process group (`killpg`), as a terminal's Ctrl-C does.
    Group,
}

/// A spawned fgumi, leading its own process group, with stderr in a file.
/// Dropping it SIGKILLs the group and any recorded aligner pid, so a failing
/// test leaves nothing behind.
struct Fgumi {
    child: Child,
    stderr: PathBuf,
    extra_pids: Vec<PathBuf>,
}

impl Fgumi {
    fn spawn(mut command: Command, dir: &Path) -> Self {
        let stderr = dir.join("fgumi.stderr");
        command
            .env("RUST_LOG", "info")
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(File::create(&stderr).expect("create stderr file"))
            .process_group(0);
        let child = command.spawn().expect("spawn fgumi");
        Self { child, stderr, extra_pids: Vec::new() }
    }

    fn pid(&self) -> Pid {
        Pid::from_child(&self.child)
    }

    fn send(&self, sig: Sig, target: Target) {
        match target {
            Target::Process => rustix::process::kill_process(self.pid(), sig.signal()),
            Target::Group => rustix::process::kill_process_group(self.pid(), sig.signal()),
        }
        .expect("send signal");
    }

    fn stderr(&self) -> String {
        std::fs::read_to_string(&self.stderr).unwrap_or_default()
    }

    /// Wait until `ready` holds, panicking (with stderr) after `within`.
    fn wait_until(&mut self, within: Duration, what: &str, ready: impl Fn(&Self) -> bool) {
        let deadline = Instant::now() + within;
        while !ready(self) {
            if let Some(status) = self.child.try_wait().expect("poll fgumi") {
                panic!("fgumi exited ({status}) before {what}:\n{}", self.stderr());
            }
            assert!(Instant::now() < deadline, "timed out waiting for {what}:\n{}", self.stderr());
            std::thread::sleep(Duration::from_millis(10));
        }
    }

    fn wait_for_notice(&mut self) {
        self.wait_until(Duration::from_secs(10), "the signal notice", |f| {
            f.stderr().contains(NOTICE)
        });
    }

    /// The exit status, or `None` if fgumi is still running after `within`.
    fn exit_within(&mut self, within: Duration) -> Option<ExitStatus> {
        let deadline = Instant::now() + within;
        loop {
            if let Some(status) = self.child.try_wait().expect("poll fgumi") {
                return Some(status);
            }
            if Instant::now() >= deadline {
                return None;
            }
            std::thread::sleep(Duration::from_millis(10));
        }
    }

    /// Assert fgumi dies by `sig` within `within`, returning how long it took.
    fn assert_dies_by(&mut self, sig: Sig, within: Duration) -> Duration {
        let start = Instant::now();
        let status = self
            .exit_within(within)
            .unwrap_or_else(|| panic!("fgumi still running after {within:?}:\n{}", self.stderr()));
        assert_eq!(
            (status.signal(), status.code()),
            (Some(sig.raw()), None),
            "fgumi must die by {} ({status}):\n{}",
            sig.name(),
            self.stderr()
        );
        start.elapsed()
    }

    fn assert_alive_for(&mut self, duration: Duration) {
        if let Some(status) = self.exit_within(duration) {
            panic!("fgumi exited ({status}) but should still be running:\n{}", self.stderr());
        }
    }
}

impl Drop for Fgumi {
    fn drop(&mut self) {
        let _ = rustix::process::kill_process_group(self.pid(), Signal::KILL);
        let _ = self.child.wait();
        for pid_file in &self.extra_pids {
            if let Some(pid) = read_pid(pid_file) {
                let _ = rustix::process::kill_process(pid, Signal::KILL);
            }
        }
    }
}

fn read_pid(path: &Path) -> Option<Pid> {
    std::fs::read_to_string(path).ok()?.trim().parse::<i32>().ok().and_then(Pid::from_raw)
}

/// Whether `pid` is still running. A zombie (an orphan whose exit init has not
/// reaped yet, as in a container whose PID 1 does not reap) counts as stopped.
fn alive(pid: Pid) -> bool {
    process_is_running(pid, Path::new("/bin/ps"))
}

fn fgumi() -> Command {
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
}

// ---------------------------------------------------------------------------
// Subprocess aligner
// ---------------------------------------------------------------------------

/// A stand-in aligner command that waits for its first input byte, then
/// backgrounds a long-lived grandchild, records its pid in `pid_file`, and
/// waits without touching its pipes again. Input arrives only once the pipeline
/// runs (the aligner is spawned while it is built), so a recorded pid means a
/// signal reaches the running pipeline, not its build.
fn idle_aligner(pid_file: &Path, trap: &str) -> String {
    format!(
        "{trap}head -c 1 > /dev/null; sleep 600 & echo $! > '{}'; wait # {{ref}}",
        pid_file.display()
    )
}

/// `fgumi realign` with a stand-in aligner (see [`idle_aligner`]). With
/// `ignore_term` the aligner's group ignores SIGTERM, so stopping it needs
/// fgumi's SIGKILL escalation.
fn realign_with_idle_aligner(dir: &Path, ignore_term: bool, envs: &[(&str, &str)]) -> Fgumi {
    let input = write_consensus(dir, 2, |_| (Some("chr1:1-100:+"), None));
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let pid_file = dir.join("aligner.pid");
    let trap = if ignore_term { "trap '' TERM; " } else { "" };
    let aligner = idle_aligner(&pid_file, trap);
    let mut command = fgumi();
    command
        .args(["realign", "-i", p(&input), "-o", p(&output), "-r", p(&reference)])
        .args(["--aligner::command", &aligner])
        .envs(envs.iter().copied());
    let mut fgumi = Fgumi::spawn(command, dir);
    fgumi.extra_pids.push(pid_file.clone());
    fgumi.wait_until(Duration::from_secs(30), "the aligner to start", |_| {
        read_pid(&pid_file).is_some()
    });
    fgumi
}

/// The aligner's grandchild is gone within `within` (an orphan is reaped by
/// init shortly after it dies).
fn assert_aligner_stopped(dir: &Path, within: Duration) {
    let pid = read_pid(&dir.join("aligner.pid")).expect("aligner pid recorded");
    let deadline = Instant::now() + within;
    while alive(pid) {
        assert!(Instant::now() < deadline, "the aligner's process {pid:?} survived fgumi");
        std::thread::sleep(Duration::from_millis(20));
    }
}

#[rstest]
fn signal_stops_the_aligner_group_and_kills_fgumi_by_that_signal(
    #[values(Sig::Int, Sig::Term, Sig::Hup)] sig: Sig,
    #[values(Target::Process, Target::Group)] target: Target,
    #[values(false, true)] ignore_term: bool,
) {
    let tmp = TempDir::new().expect("tempdir");
    let mut fgumi = realign_with_idle_aligner(tmp.path(), ignore_term, &[]);

    fgumi.send(sig, target);
    fgumi.assert_dies_by(sig, Duration::from_secs(7));
    assert_aligner_stopped(tmp.path(), Duration::from_secs(10));

    let stderr = fgumi.stderr();
    assert!(stderr.contains(&format!("Received {}: {NOTICE}", sig.name())), "{stderr}");
    assert_reported_as_cancellation(&stderr);
}

/// stderr of a pipeline command stopped by a signal: the run is reported as
/// cancelled, not as a failure.
fn assert_reported_as_cancellation(stderr: &str) {
    assert!(stderr.contains("Run cancelled; its output is incomplete"), "{stderr}");
    assert!(!stderr.contains("Run aborted"), "{stderr}");
    assert!(!stderr.contains("Error:"), "a cancellation is not an error:\n{stderr}");
    assert!(!stderr.contains("panicked"), "{stderr}");
}

/// `FGUMI_SIGNAL_GRACE_SECS=0` turns the handling off: the signal's default
/// action kills fgumi at once, with no notice.
#[test]
fn zero_grace_keeps_the_default_action() {
    let tmp = TempDir::new().expect("tempdir");
    let mut fgumi =
        realign_with_idle_aligner(tmp.path(), false, &[("FGUMI_SIGNAL_GRACE_SECS", "0")]);
    fgumi.send(Sig::Term, Target::Process);
    fgumi.assert_dies_by(Sig::Term, Duration::from_secs(2));
    assert!(!fgumi.stderr().contains(NOTICE), "{}", fgumi.stderr());
}

/// A signal inherited as ignored (`nohup`) stays ignored; the others still
/// stop fgumi gracefully.
#[test]
fn inherited_ignored_signal_stays_ignored() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let input = write_consensus(dir, 2, |_| (Some("chr1:1-100:+"), None));
    let reference = create_test_reference(dir);
    let output = dir.join("out.bam");
    let pid_file = dir.join("aligner.pid");
    let aligner = idle_aligner(&pid_file, "");
    let mut command = Command::new("/bin/sh");
    command
        .arg("-c")
        .arg("trap '' HUP; exec \"$0\" \"$@\"")
        .arg(env!("CARGO_BIN_EXE_fgumi"))
        .args(["realign", "-i", p(&input), "-o", p(&output), "-r", p(&reference)])
        .args(["--aligner::command", &aligner]);
    let mut fgumi = Fgumi::spawn(command, dir);
    fgumi.extra_pids.push(pid_file.clone());
    fgumi.wait_until(Duration::from_secs(30), "the aligner to start", |_| {
        read_pid(&pid_file).is_some()
    });

    fgumi.send(Sig::Hup, Target::Process);
    fgumi.assert_alive_for(Duration::from_secs(1));
    assert!(alive(read_pid(&pid_file).expect("pid")), "the aligner must survive SIGHUP");

    fgumi.send(Sig::Term, Target::Process);
    fgumi.assert_dies_by(Sig::Term, Duration::from_secs(7));
    assert_aligner_stopped(dir, Duration::from_secs(10));
}

// ---------------------------------------------------------------------------
// Sort and merge cleanup
// ---------------------------------------------------------------------------

/// The BGZF end-of-file marker every complete BAM ends with.
const BGZF_EOF_LEN: usize = 28;

/// Create a FIFO at `path`.
fn make_fifo(path: &Path) {
    let status = Command::new("mkfifo").arg(path).status().expect("run mkfifo");
    assert!(status.success(), "mkfifo {} failed", path.display());
}

/// A FIFO fed with `bytes` (minus the BGZF EOF marker, so the input can never
/// complete) by a thread that keeps its write end open until told to close
/// it, as an upstream process does until it too is interrupted.
struct HeldFifo {
    close: Option<mpsc::Sender<()>>,
    writer: Option<JoinHandle<()>>,
}

impl HeldFifo {
    /// Create the FIFO at `path` and start feeding it.
    fn new(path: &Path, bytes: Vec<u8>) -> Self {
        make_fifo(path);
        Self::feed(path, bytes)
    }

    /// Start feeding the existing FIFO at `path`.
    fn feed(path: &Path, bytes: Vec<u8>) -> Self {
        let (close, closed) = mpsc::channel::<()>();
        let path = path.to_path_buf();
        let writer = std::thread::spawn(move || {
            let Ok(mut fifo) = File::options().write(true).open(&path) else { return };
            let _ = fifo.write_all(&bytes[..bytes.len() - BGZF_EOF_LEN]);
            let _ = closed.recv_timeout(Duration::from_mins(1));
        });
        Self { close: Some(close), writer: Some(writer) }
    }

    fn close(&mut self) {
        self.close.take();
        if let Some(writer) = self.writer.take() {
            let _ = writer.join();
        }
    }
}

impl Drop for HeldFifo {
    fn drop(&mut self) {
        self.close.take();
        // The writer may still be blocked opening a FIFO nobody read; leave it.
    }
}

fn unsorted_bam(dir: &Path, n: usize) -> PathBuf {
    let path = dir.join("unsorted.bam");
    write_bam(&path, &create_minimal_header("chr1", 1_000_000), &records(n, |i| format!("q{i}")));
    path
}

/// Every regular file under `dir`, recursively.
fn files_under(dir: &Path) -> Vec<PathBuf> {
    let mut found = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).into_iter().flatten().flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else {
                found.push(path);
            }
        }
    }
    found
}

fn entries_in(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .expect("read dir")
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

const SENTINEL: &[u8] = b"an earlier run's output, which an interrupted run must not replace";

/// An interrupted sort removes its spill directory and staged output and
/// leaves the existing destination (and its index) untouched.
#[rstest]
#[case::int(Sig::Int, false)]
#[case::term(Sig::Term, false)]
#[case::hup(Sig::Hup, false)]
#[case::term_with_index(Sig::Term, true)]
fn interrupted_sort_cleans_up_and_keeps_the_destination(
    #[case] sig: Sig,
    #[case] write_index: bool,
) {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let bam = std::fs::read(unsorted_bam(dir, 20_000)).expect("read fixture");
    let fifo_path = dir.join("input.fifo");
    let mut fifo = HeldFifo::new(&fifo_path, bam);
    let out_dir = dir.join("out");
    let spill_dir = dir.join("spill");
    std::fs::create_dir_all(&out_dir).expect("mkdir");
    std::fs::create_dir_all(&spill_dir).expect("mkdir");
    let output = out_dir.join("dest.bam");
    std::fs::write(&output, SENTINEL).expect("write sentinel");
    if write_index {
        std::fs::write(out_dir.join("dest.bam.bai"), SENTINEL).expect("write sentinel index");
    }
    let before = entries_in(&out_dir);

    let mut command = fgumi();
    command
        .args(["sort", "-i", p(&fifo_path), "-o", p(&output), "--tmp-dir", p(&spill_dir)])
        .args(["--order", "coordinate", "--threads", "2"])
        .args(["--max-memory", "64K", "--memory-per-thread", "false"]);
    if write_index {
        command.arg("--write-index");
    }
    let mut fgumi = Fgumi::spawn(command, dir);
    fgumi.wait_until(Duration::from_secs(30), "a spill file", |_| {
        !files_under(&spill_dir).is_empty()
    });

    fgumi.send(sig, Target::Group);
    fgumi.wait_for_notice();
    fifo.close();
    fgumi.assert_dies_by(sig, Duration::from_secs(7));

    assert_eq!(files_under(&spill_dir), Vec::<PathBuf>::new(), "spill files left behind");
    assert_eq!(entries_in(&spill_dir), Vec::<String>::new(), "spill directory left behind");
    assert_eq!(entries_in(&out_dir), before, "staged output left behind");
    assert_eq!(std::fs::read(&output).expect("read dest"), SENTINEL);
    if write_index {
        assert_eq!(std::fs::read(out_dir.join("dest.bam.bai")).expect("read bai"), SENTINEL);
    }
    assert_reported_as_cancellation(&fgumi.stderr());
}

/// An interrupted merge removes its staged output and leaves the existing
/// destination untouched.
#[rstest]
fn interrupted_merge_cleans_up_and_keeps_the_destination(#[values(Sig::Int, Sig::Term)] sig: Sig) {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let unsorted = unsorted_bam(dir, 20_000);
    let sorted = dir.join("sorted.bam");
    let status = fgumi()
        .args(["sort", "-i", p(&unsorted), "-o", p(&sorted), "--order", "coordinate"])
        .stderr(Stdio::null())
        .status()
        .expect("run fgumi sort");
    assert!(status.success(), "fixture sort failed");
    let bytes = std::fs::read(&sorted).expect("read sorted");
    let fifo_path = dir.join("input.fifo");
    let mut fifo = HeldFifo::new(&fifo_path, bytes);
    let out_dir = dir.join("out");
    std::fs::create_dir_all(&out_dir).expect("mkdir");
    let output = out_dir.join("dest.bam");
    std::fs::write(&output, SENTINEL).expect("write sentinel");
    let before = entries_in(&out_dir);

    let mut command = fgumi();
    command
        .args(["merge", p(&sorted), p(&fifo_path), "-o", p(&output)])
        .args(["--order", "coordinate"]);
    let mut fgumi = Fgumi::spawn(command, dir);
    fgumi.wait_until(Duration::from_secs(30), "the staged output", |_| {
        entries_in(&out_dir).len() > before.len()
    });

    fgumi.send(sig, Target::Group);
    fgumi.wait_for_notice();
    fifo.close();
    fgumi.assert_dies_by(sig, Duration::from_secs(7));

    assert_eq!(entries_in(&out_dir), before, "staged output left behind");
    assert_eq!(std::fs::read(&output).expect("read dest"), SENTINEL);
    assert_reported_as_cancellation(&fgumi.stderr());
}

// ---------------------------------------------------------------------------
// Signals with nothing to clean up, or during setup
// ---------------------------------------------------------------------------

/// A command with no cleanup scope (`sort --verify` only reads) dies by the
/// signal at once, with no notice.
#[test]
fn signal_with_nothing_registered_kills_at_once() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let fifo_path = dir.join("input.fifo");
    make_fifo(&fifo_path);
    let mut command = fgumi();
    command.args(["sort", "--verify", "-i", p(&fifo_path), "--order", "coordinate"]);
    let mut fgumi = Fgumi::spawn(command, dir);
    // Blocked opening a FIFO nobody writes.
    fgumi.wait_until(Duration::from_secs(10), "startup", |f| {
        f.stderr().contains("Running fgumi version")
    });
    fgumi.send(Sig::Int, Target::Group);
    fgumi.assert_dies_by(Sig::Int, Duration::from_secs(2));
    assert!(!fgumi.stderr().contains(NOTICE), "{}", fgumi.stderr());
}

/// `fgumi sort` stuck in setup: its output is staged (so a cleanup scope is
/// open) and it is blocked opening a FIFO input nobody writes yet, so an
/// accepted signal cannot make progress until the FIFO is fed. Returns fgumi,
/// the FIFO and the output directory.
fn sort_blocked_in_setup(dir: &Path, grace: &str) -> (Fgumi, PathBuf, PathBuf) {
    let fifo_path = dir.join("input.fifo");
    make_fifo(&fifo_path);
    let out_dir = dir.join("out");
    std::fs::create_dir_all(&out_dir).expect("mkdir");
    let output = out_dir.join("dest.bam");
    let mut command = fgumi();
    command
        .args(["sort", "-i", p(&fifo_path), "-o", p(&output), "--order", "coordinate"])
        .env("FGUMI_SIGNAL_GRACE_SECS", grace);
    let mut fgumi = Fgumi::spawn(command, dir);
    fgumi.wait_until(Duration::from_secs(10), "the staged output", |_| {
        !entries_in(&out_dir).is_empty()
    });
    (fgumi, fifo_path, out_dir)
}

/// The watchdog bounds how long fgumi waits for a cleanup that cannot
/// progress, and kills by the first signal.
#[test]
fn watchdog_terminates_a_stalled_cleanup() {
    let tmp = TempDir::new().expect("tempdir");
    let (mut fgumi, _fifo, _out) = sort_blocked_in_setup(tmp.path(), "1.5");
    fgumi.send(Sig::Hup, Target::Group);
    let elapsed = fgumi.assert_dies_by(Sig::Hup, Duration::from_secs(5));
    assert!(elapsed >= Duration::from_millis(1400), "died before the grace period: {elapsed:?}");
}

/// A second signal, sent after the debounce, kills at once by that signal.
#[test]
fn second_signal_terminates_a_stalled_cleanup() {
    let tmp = TempDir::new().expect("tempdir");
    let (mut fgumi, _fifo, _out) = sort_blocked_in_setup(tmp.path(), "60");
    fgumi.send(Sig::Term, Target::Process);
    fgumi.wait_for_notice();
    fgumi.assert_alive_for(Duration::from_millis(700));
    fgumi.send(Sig::Int, Target::Process);
    fgumi.assert_dies_by(Sig::Int, Duration::from_secs(2));
}

/// A second signal inside the debounce is a duplicate of the first, not a
/// request to stop at once; a later one still kills. Two different signals,
/// because repeats of one signal coalesce before fgumi sees them.
#[test]
fn second_signal_inside_the_debounce_is_ignored() {
    let tmp = TempDir::new().expect("tempdir");
    let (mut fgumi, _fifo, _out) = sort_blocked_in_setup(tmp.path(), "60");
    fgumi.send(Sig::Term, Target::Process);
    fgumi.send(Sig::Int, Target::Process);
    fgumi.assert_alive_for(Duration::from_millis(700));
    fgumi.send(Sig::Hup, Target::Process);
    fgumi.assert_dies_by(Sig::Hup, Duration::from_secs(2));
}

/// A signal accepted while the chain cannot be built yet is held: once setup
/// can proceed, the build notices it, the staged output is removed, and fgumi
/// dies by the signal promptly, well inside the grace period.
#[test]
fn signal_during_setup_cancels_once_setup_can_proceed() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let bam = std::fs::read(unsorted_bam(dir, 2_000)).expect("read fixture");
    let (mut fgumi, fifo_path, out_dir) = sort_blocked_in_setup(dir, "30");
    fgumi.send(Sig::Term, Target::Process);
    fgumi.wait_for_notice();
    let _fifo = HeldFifo::feed(&fifo_path, bam);
    fgumi.assert_dies_by(Sig::Term, Duration::from_secs(7));
    assert_eq!(entries_in(&out_dir), Vec::<String>::new(), "staged output left behind");
    let stderr = fgumi.stderr();
    assert!(!stderr.contains("Error:"), "{stderr}");
}

/// A command that writes its output directly (`group`) and is signalled
/// before its sink exists leaves an existing output untouched: the build
/// notices the request before it creates (and truncates) the output.
#[test]
fn signal_during_setup_leaves_a_direct_output_untouched() {
    let tmp = TempDir::new().expect("tempdir");
    let dir = tmp.path();
    let unsorted = unsorted_bam(dir, 2_000);
    let sorted = dir.join("sorted.bam");
    let status = fgumi()
        .args(["sort", "-i", p(&unsorted), "-o", p(&sorted), "--order", "template-coordinate"])
        .stderr(Stdio::null())
        .status()
        .expect("run fgumi sort");
    assert!(status.success(), "fixture sort failed");
    let bam = std::fs::read(&sorted).expect("read sorted");
    let fifo_path = dir.join("input.fifo");
    make_fifo(&fifo_path);
    let output = dir.join("existing.bam");
    std::fs::write(&output, SENTINEL).expect("write sentinel");

    let mut command = fgumi();
    command
        .args(["group", "-i", p(&fifo_path), "-o", p(&output), "--strategy", "identity"])
        .env("FGUMI_SIGNAL_GRACE_SECS", "30");
    let mut fgumi = Fgumi::spawn(command, dir);
    fgumi.wait_until(Duration::from_secs(10), "startup", |f| {
        f.stderr().contains("Running fgumi version")
    });
    fgumi.send(Sig::Term, Target::Process);
    fgumi.wait_for_notice();
    let _fifo = HeldFifo::feed(&fifo_path, bam);
    fgumi.assert_dies_by(Sig::Term, Duration::from_secs(7));
    assert_eq!(std::fs::read(&output).expect("read output"), SENTINEL, "output was replaced");
}
