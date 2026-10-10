//! Output sinks whose close is checked.
//!
//! Dropping a [`File`] closes its descriptor and discards the result, and
//! `File::flush` is a no-op on Unix. So an error the OS reports only at `close`
//! (`ENOSPC`/`EDQUOT` or `EIO` on an NFS mount, which flushes dirty pages
//! there) is lost, and the command exits zero with a short or corrupt output.
//!
//! [`OutputFile`] owns an output file and finishes it with [`OutputFile::close`],
//! a checked `close(2)`. It deliberately does not sync the file's data to
//! storage first, as fgumi 0.7.0 did not. This differs from htslib (and so
//! samtools), whose `bgzf_close` reaches `fdatasync` through `hflush`. Syncing
//! made the final close wait out the whole write-back: closing a 49 GB
//! coordinate-sorted output took 61.5 s with the sync against 20-28 s without
//! it (c7g.4xlarge, gp3). Durability against a host crash is the caller's job
//! (e.g. `sync` after the command). The cost is that a write-back error the OS
//! reports only after the close is not seen; one reported at the close itself
//! still fails the command.
//!
//! Two outputs outside this module do sync, each for its own reason: the
//! metrics writer's temp-and-rename (`fgumi-metrics`, small files) and
//! `review`'s staged grouped BAM.
//! [`OutputSink`] is the type-erased form, so a writer that may target either a
//! file or stdout can be finished the same way; [`open_output_sink`] opens one
//! for a path, honouring `-` and `/dev/stdout`.

use std::fs::File;
use std::io::{self, BufWriter, Stdout, Write};
use std::path::Path;

use anyhow::Context;
use tempfile::NamedTempFile;

use crate::paths::is_stdout_path;

/// A writer with an explicit, fallible close.
///
/// Call [`OutputSink::close`] once all data has been written. Dropping a sink
/// without closing it releases it but discards any error, so drop only on an
/// abort path where the output is already known to be bad.
pub trait OutputSink: Write + Send {
    /// Flush and release the sink, returning the first error encountered.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing or closing fails.
    fn close(self: Box<Self>) -> io::Result<()>;
}

impl<S: OutputSink + ?Sized> OutputSink for Box<S> {
    fn close(self: Box<Self>) -> io::Result<()> {
        S::close(*self)
    }
}

/// Stdout cannot be closed, so closing it only flushes. Used only
/// where stdout cannot be duplicated (targets without file descriptors); on
/// Unix [`open_output_sink`] returns a checked duplicate instead.
impl OutputSink for Stdout {
    fn close(mut self: Box<Self>) -> io::Result<()> {
        self.flush()
    }
}

/// An output file that is closed with its errors checked. Its data is not
/// synced to storage.
///
/// Created with [`OutputFile::create`] or wrapped around an open [`File`] with
/// [`From`]. Writes go straight to the file (wrap it in a [`BufWriter`] for
/// small writes). Finish with [`OutputFile::close`], or [`close_buffered`] when
/// it sits under a `BufWriter`.
#[derive(Debug)]
pub struct OutputFile {
    file: File,
}

impl OutputFile {
    /// Create (or truncate) the file at `path`.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created.
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        File::create(path).map(Self::from)
    }

    /// Close the file, returning any error `close(2)` reports.
    ///
    /// The data is not synced to storage first, so an error the OS defers past
    /// the close (write-back failing after it) is not seen; one reported at the
    /// close itself (e.g. `ENOSPC`/`EDQUOT` on NFS) is. On Unix an `EINTR` from
    /// `close` is treated as success: the descriptor is released regardless,
    /// and retrying could close a descriptor another thread has since opened.
    ///
    /// # Errors
    ///
    /// Returns an error if the close fails.
    pub fn close(self) -> io::Result<()> {
        close_file(self.file)
    }
}

impl From<File> for OutputFile {
    fn from(file: File) -> Self {
        Self { file }
    }
}

impl Write for OutputFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.file.write(buf)
    }

    fn write_vectored(&mut self, bufs: &[io::IoSlice<'_>]) -> io::Result<usize> {
        self.file.write_vectored(bufs)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.flush()
    }
}

impl OutputSink for OutputFile {
    fn close(self: Box<Self>) -> io::Result<()> {
        OutputFile::close(*self)
    }
}

/// Flush `writer`, then close the sink beneath it.
///
/// # Errors
///
/// Returns an error if the flush, or closing the inner sink, fails.
pub fn close_buffered<W: OutputSink>(writer: BufWriter<W>) -> io::Result<()> {
    let inner = writer.into_inner().map_err(io::IntoInnerError::into_error)?;
    Box::new(inner).close()
}

/// Open an output sink for `path`: stdout for `-` or `/dev/stdout`, otherwise
/// a newly created (or truncated) [`OutputFile`].
///
/// Stdout is a block-buffered duplicate of fd 1 rather than
/// [`std::io::Stdout`], whose `LineWriter` splits every write at the last
/// `\n`; binary output such as BGZF carries `0x0a` at arbitrary offsets, so a
/// large flush would be torn into many small writes. Its close is checked (on
/// NFS every `close(2)` flushes the file's dirty pages, so a `> out.bam`
/// redirect onto a full or failing mount reports the error there). The
/// duplicate is closed, never fd 1 itself.
///
/// # Errors
///
/// Returns an error if the file cannot be created or stdout cannot be
/// duplicated.
pub fn open_output_sink<P: AsRef<Path>>(path: P) -> io::Result<Box<dyn OutputSink>> {
    let path = path.as_ref();
    if is_stdout_path(path) {
        block_buffered_stdout()
    } else {
        Ok(Box::new(OutputFile::create(path)?))
    }
}

/// A block-buffered, close-checked duplicate of this process's stdout (see
/// [`open_output_sink`]).
#[cfg(unix)]
fn block_buffered_stdout() -> io::Result<Box<dyn OutputSink>> {
    use std::os::fd::AsFd;

    let dup = io::stdout().as_fd().try_clone_to_owned()?;
    Ok(Box::new(OutputFile::from(File::from(dup))))
}

/// Stdout itself: targets without file descriptors cannot duplicate it, so the
/// line-buffering concern on the Unix variant cannot be worked around there.
#[cfg(not(unix))]
#[expect(clippy::unnecessary_wraps, reason = "matches the fallible Unix variant")]
fn block_buffered_stdout() -> io::Result<Box<dyn OutputSink>> {
    Ok(Box::new(io::stdout()))
}

/// Finish a temp file and rename it onto `dest`: `close` releases the temp's
/// own handle, and only once that succeeds is the temp renamed.
///
/// Production callers pass `|f| OutputFile::from(f).close()`; `close` is a
/// parameter so tests (here and in `fgumi-sort`) can inject a close failure,
/// which a real file cannot be made to produce portably.
///
/// So a close error fails the write instead of leaving a renamed but
/// incomplete file at `dest`. On any error the temp is removed and `dest` is
/// left untouched. The rename is atomic with respect to a failure of this
/// process only: the data is not synced before it, so after a host crash or
/// power loss the rename can be durable ahead of the data, leaving an
/// incomplete file at `dest`. Re-stamp the temp's mode (`restamp_for_persist`)
/// before calling this.
///
/// # Errors
///
/// Returns an error if `close` fails or the rename fails.
pub fn persist_after_close<F>(tmp: NamedTempFile, dest: &Path, close: F) -> anyhow::Result<()>
where
    F: FnOnce(File) -> io::Result<()>,
{
    let (file, temp_path) = tmp.into_parts();
    // On error `temp_path` is dropped, removing the temp.
    close(file).with_context(|| {
        format!("Failed to close temp file before renaming onto: {}", dest.display())
    })?;
    temp_path
        .persist(dest)
        .map_err(|e| e.error)
        .with_context(|| format!("Failed to rename temp file onto: {}", dest.display()))
}

/// Close `file`, returning the error `close(2)` reports.
#[cfg(unix)]
fn close_file(file: File) -> io::Result<()> {
    close_result(nix::unistd::close(file))
}

/// Close `file`. Targets without file descriptors offer no checked close.
#[cfg(not(unix))]
#[allow(clippy::unnecessary_wraps)]
fn close_file(file: File) -> io::Result<()> {
    drop(file);
    Ok(())
}

/// Map a `close(2)` result to an I/O result, treating `EINTR` as success.
///
/// POSIX leaves the descriptor's state unspecified after `EINTR`, but Linux,
/// macOS, and the BSDs all release it, so the close must not be retried.
#[cfg(unix)]
fn close_result(result: nix::Result<()>) -> io::Result<()> {
    match result {
        Ok(()) | Err(nix::errno::Errno::EINTR) => Ok(()),
        Err(errno) => Err(io::Error::from(errno)),
    }
}

/// Test doubles shared by this crate's tests.
#[cfg(test)]
pub(crate) mod test_support {
    use std::io::{self, Write};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};

    use super::OutputSink;

    /// The error a failing [`CloseProbe`] close returns.
    pub(crate) const CLOSE_ERROR: &str = "input/output error at close";

    /// An in-memory [`OutputSink`] that records its bytes, counts closes, and
    /// can fail its writes or its close (modelling an error reported at close,
    /// such as an NFS flush failing). Clones share state.
    #[derive(Clone, Default)]
    pub(crate) struct CloseProbe {
        pub(crate) written: Arc<Mutex<Vec<u8>>>,
        pub(crate) closes: Arc<AtomicUsize>,
        pub(crate) fail_write: bool,
        pub(crate) fail_close: bool,
    }

    impl CloseProbe {
        /// A probe whose close fails with [`CLOSE_ERROR`].
        pub(crate) fn failing_close() -> Self {
            Self { fail_close: true, ..Self::default() }
        }

        pub(crate) fn bytes(&self) -> Vec<u8> {
            self.written.lock().expect("probe lock").clone()
        }

        pub(crate) fn closes(&self) -> usize {
            self.closes.load(Ordering::SeqCst)
        }
    }

    impl Write for CloseProbe {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.fail_write {
                return Err(io::Error::new(io::ErrorKind::StorageFull, "write failed"));
            }
            self.written.lock().expect("probe lock").extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    impl OutputSink for CloseProbe {
        fn close(self: Box<Self>) -> io::Result<()> {
            self.closes.fetch_add(1, Ordering::SeqCst);
            if self.fail_close {
                return Err(io::Error::other(CLOSE_ERROR));
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::{CLOSE_ERROR, CloseProbe};
    use super::*;
    use std::io::Read;

    #[test]
    fn test_output_file_close_writes_all_bytes() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.txt");
        let mut out = OutputFile::create(&path).unwrap();
        out.write_all(b"hello world").unwrap();
        out.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"hello world");
    }

    /// Closing an output that is not a regular file (a pipe, `/dev/null`) must
    /// succeed.
    #[cfg(unix)]
    #[test]
    fn test_output_file_close_succeeds_on_non_regular_files() {
        let dev_null = std::fs::OpenOptions::new().write(true).open("/dev/null").unwrap();
        let mut out = OutputFile::from(dev_null);
        out.write_all(b"discarded").unwrap();
        out.close().expect("closing /dev/null must succeed");

        let (mut reader, writer) = std::io::pipe().unwrap();
        let mut out = OutputFile::from(File::from(std::os::fd::OwnedFd::from(writer)));
        out.write_all(b"through the pipe").unwrap();
        out.close().expect("closing a pipe must succeed");
        let mut received = Vec::new();
        reader.read_to_end(&mut received).unwrap();
        assert_eq!(received, b"through the pipe");
    }

    /// Whether this test runs alone in its process, as nextest runs every test.
    ///
    /// The `EBADF` tests below close a descriptor out from under its `File`.
    /// Between that close and the checked one, another thread can be handed the
    /// same descriptor number, and the checked close would then close *its*
    /// file. Under the threaded `cargo test` harness they skip rather than risk
    /// closing a sibling test's file.
    #[cfg(unix)]
    fn runs_in_own_process() -> bool {
        let isolated =
            std::env::var("NEXTEST_EXECUTION_MODE").is_ok_and(|mode| mode == "process-per-test");
        if !isolated {
            eprintln!("skipped: closing a raw descriptor needs nextest's process-per-test mode");
        }
        isolated
    }

    /// The close error must surface. The descriptor is closed underneath the
    /// `File`, so the checked close sees `EBADF`; an unchecked drop would hide it.
    #[cfg(unix)]
    #[test]
    fn test_close_file_surfaces_close_error() {
        use std::os::fd::AsRawFd;

        if !runs_in_own_process() {
            return;
        }
        let file = tempfile::tempfile().unwrap();
        nix::unistd::close(file.as_raw_fd()).unwrap();
        let err = close_file(file).expect_err("closing a closed descriptor must fail");
        assert_eq!(err.raw_os_error(), Some(nix::errno::Errno::EBADF as i32), "{err}");
    }

    /// The same failure through the public API: an `OutputFile` whose
    /// descriptor is no longer valid must not close cleanly.
    #[cfg(unix)]
    #[test]
    fn test_output_file_close_surfaces_error() {
        use std::os::fd::AsRawFd;

        if !runs_in_own_process() {
            return;
        }
        let dir = tempfile::tempdir().unwrap();
        let file = File::create(dir.path().join("out.txt")).unwrap();
        let fd = file.as_raw_fd();
        let out = OutputFile::from(file);
        nix::unistd::close(fd).unwrap();
        let err = out.close().expect_err("closing a closed descriptor must fail");
        assert_eq!(err.raw_os_error(), Some(nix::errno::Errno::EBADF as i32), "{err}");
    }

    /// A boxed output file (the stdout duplicate) is flushed and its
    /// descriptor released on close; a pipe stands in for stdout so the bytes
    /// can be read back.
    #[cfg(unix)]
    #[test]
    fn test_boxed_output_file_close_flushes_and_closes() {
        let (mut reader, writer) = std::io::pipe().unwrap();
        let mut sink: Box<dyn OutputSink> =
            Box::new(OutputFile::from(File::from(std::os::fd::OwnedFd::from(writer))));
        sink.write_all(b"to stdout").unwrap();
        sink.close().expect("closing the stdout duplicate must succeed");
        // EOF is seen only once the duplicate (the pipe's last writer) is closed.
        let mut received = Vec::new();
        reader.read_to_end(&mut received).unwrap();
        assert_eq!(received, b"to stdout");
    }

    /// fgumi's outputs must not sync their data to storage on close: on a large
    /// output that waits out the whole write-back (61.5 s against 20-28 s to
    /// close a 49 GB coordinate-sorted output).
    ///
    /// A sync leaves no trace a test can observe without a seam built only for
    /// it, so this scans the source of every crate in the workspace (`src/` and
    /// `crates/`, production and test code alike) for a sync call. The only
    /// syncs allowed are the [`ALLOWED_SYNCS`] exceptions, and each must still
    /// be present, so the list cannot go stale. Outside the workspace (e.g. a
    /// packaged crate) there is no source tree to scan and the test is skipped.
    #[test]
    fn test_workspace_issues_no_unlisted_sync() {
        // Spelled in pieces so this test does not match itself.
        let calls = [
            concat!("sync", "_data"),
            concat!("sync", "_all"),
            concat!("fdata", "sync"),
            concat!("f", "sync"),
            concat!("F_FULL", "FSYNC"),
        ];
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
        if !root.join("crates/fgumi-bam-io/src/output.rs").is_file() {
            eprintln!("skipped: not built inside the fgumi workspace");
            return;
        }
        let mut files = Vec::new();
        collect_rust_files(&root.join("src"), &mut files);
        collect_rust_files(&root.join("crates"), &mut files);
        assert!(
            files.len() > 50,
            "found only {} source files under {}",
            files.len(),
            root.display()
        );

        let mut found: Vec<(String, String)> = Vec::new();
        for file in &files {
            let source = std::fs::read_to_string(file).unwrap();
            let rel = file.strip_prefix(&root).unwrap().to_string_lossy().replace('\\', "/");
            for line in source.lines().filter(|line| !line.trim_start().starts_with("//")) {
                for call in calls.iter().filter(|call| line.contains(**call)) {
                    found.push((rel.clone(), (*call).to_string()));
                }
            }
        }
        found.sort();
        let mut allowed: Vec<(String, String)> =
            ALLOWED_SYNCS.iter().map(|(f, c)| ((*f).to_string(), (*c).to_string())).collect();
        allowed.sort();
        assert_eq!(found, allowed, "sync calls in the workspace differ from ALLOWED_SYNCS");
    }

    /// The syncs fgumi does issue: `(path from the workspace root, call)`.
    const ALLOWED_SYNCS: &[(&str, &str)] = &[
        // Metrics temp-and-rename: small files, synced so a host crash cannot
        // persist the rename ahead of the data.
        ("crates/fgumi-metrics/src/writer.rs", concat!("sync", "_all")),
        // `review`'s staged grouped BAM.
        ("src/lib/commands/review.rs", concat!("sync", "_all")),
    ];

    /// Append every `.rs` file under `dir` to `files`, skipping build output.
    fn collect_rust_files(dir: &Path, files: &mut Vec<std::path::PathBuf>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                if path.file_name().is_some_and(|name| name != "target") {
                    collect_rust_files(&path, files);
                }
            } else if path.extension().is_some_and(|ext| ext == "rs") {
                files.push(path);
            }
        }
    }

    /// Both spellings of stdout open a sink that closes cleanly.
    #[test]
    fn test_open_output_sink_stdout_closes_cleanly() {
        for path in ["-", "/dev/stdout"] {
            let mut sink = open_output_sink(path).unwrap();
            sink.flush().unwrap();
            sink.close().unwrap();
        }
    }

    #[test]
    fn test_open_output_sink_creates_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.txt");
        let mut sink = open_output_sink(&path).unwrap();
        sink.write_all(b"file").unwrap();
        sink.close().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"file");
    }

    #[cfg(unix)]
    #[test]
    fn test_close_result_treats_eintr_as_success() {
        use nix::errno::Errno;

        close_result(Ok(())).unwrap();
        close_result(Err(Errno::EINTR)).expect("EINTR releases the descriptor; not an error");
        let err = close_result(Err(Errno::EIO)).expect_err("EIO must surface");
        assert_eq!(err.raw_os_error(), Some(Errno::EIO as i32));
    }

    #[test]
    fn test_close_buffered_flushes_then_closes() {
        let probe = CloseProbe::default();
        let mut writer = BufWriter::new(probe.clone());
        writer.write_all(b"buffered").unwrap();
        close_buffered(writer).unwrap();
        assert_eq!(probe.bytes(), b"buffered");
        assert_eq!(probe.closes(), 1, "inner sink must be closed");
    }

    #[test]
    fn test_close_buffered_surfaces_close_error() {
        let err = close_buffered(BufWriter::new(CloseProbe::failing_close()))
            .expect_err("close error must surface");
        assert_eq!(err.to_string(), CLOSE_ERROR);
    }

    /// Writing out the buffered bytes can fail; that error must surface, and
    /// the sink must not be closed as though the data were complete.
    #[test]
    fn test_close_buffered_surfaces_buffer_write_error() {
        let probe = CloseProbe { fail_write: true, ..CloseProbe::default() };
        let mut writer = BufWriter::new(probe.clone());
        writer.write_all(b"pending").unwrap();
        let err = close_buffered(writer).expect_err("buffered write error must surface");
        assert_eq!(err.kind(), io::ErrorKind::StorageFull);
        assert_eq!(probe.closes(), 0);
    }

    /// A boxed sink closes the sink inside it, so `Box<dyn OutputSink>` works
    /// wherever a concrete sink does.
    #[test]
    fn test_boxed_sink_delegates_close() {
        let sink: Box<dyn OutputSink> = Box::new(CloseProbe::failing_close());
        let err = close_buffered(BufWriter::new(sink)).expect_err("close error must surface");
        assert_eq!(err.to_string(), CLOSE_ERROR);
    }

    /// The temp is closed before the rename: when `close` runs, `dest` does
    /// not exist yet, and afterwards it holds the temp's bytes.
    #[test]
    fn test_persist_after_close_closes_before_rename() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("out.txt");
        let mut tmp = NamedTempFile::new_in(dir.path()).unwrap();
        tmp.write_all(b"payload").unwrap();
        persist_after_close(tmp, &dest, |file| {
            assert!(!dest.exists(), "the temp must be closed before it is renamed");
            OutputFile::from(file).close()
        })
        .unwrap();
        assert_eq!(std::fs::read(&dest).unwrap(), b"payload");
    }

    /// A failed close fails the write: nothing is renamed onto `dest` (an
    /// existing file there is untouched) and the temp is removed.
    #[test]
    fn test_persist_after_close_failure_leaves_dest_and_removes_temp() {
        let dir = tempfile::tempdir().unwrap();
        let dest = dir.path().join("out.txt");
        std::fs::write(&dest, b"previous").unwrap();
        let mut tmp = NamedTempFile::new_in(dir.path()).unwrap();
        tmp.write_all(b"new").unwrap();
        let temp_path = tmp.path().to_path_buf();
        let err = persist_after_close(tmp, &dest, |_| Err(io::Error::other(CLOSE_ERROR)))
            .expect_err("a close error must fail the write");
        assert!(format!("{err:#}").contains(CLOSE_ERROR), "{err:#}");
        assert!(format!("{err:#}").contains("Failed to close temp file"), "{err:#}");
        assert_eq!(std::fs::read(&dest).unwrap(), b"previous");
        assert!(!temp_path.exists(), "the temp must be removed");
    }
}
