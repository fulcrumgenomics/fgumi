//! Output sinks whose close is checked.
//!
//! Dropping a [`File`] closes its descriptor and discards the result, and
//! `File::flush` is a no-op on Unix. So an error the OS reports only after the
//! last `write` (deferred write-back failing with `EIO`, or `ENOSPC`/`EDQUOT` on
//! an NFS mount, which flushes dirty pages at `close`) is lost, and the command
//! exits zero with a short or corrupt output.
//!
//! [`OutputFile`] owns an output file and finishes it with [`OutputFile::close`]:
//! `sync_data` (regular files only), then a checked `close(2)`. This follows
//! htslib's `bgzf_close`, which syncs through `hflush` and then checks `close`;
//! like htslib, a sync the file or filesystem does not support is skipped
//! rather than failed (see [`OutputFile::close`]).
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
    /// Flush, make the data durable where applicable, and release the sink,
    /// returning the first error encountered.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing, syncing, or closing fails.
    fn close(self: Box<Self>) -> io::Result<()>;
}

impl<S: OutputSink + ?Sized> OutputSink for Box<S> {
    fn close(self: Box<Self>) -> io::Result<()> {
        S::close(*self)
    }
}

/// Stdout cannot be closed or synced, so closing it only flushes. Used only
/// where stdout cannot be duplicated (targets without file descriptors); on
/// Unix [`open_output_sink`] returns a checked duplicate instead.
impl OutputSink for Stdout {
    fn close(mut self: Box<Self>) -> io::Result<()> {
        self.flush()
    }
}

/// An output file that is closed with its errors checked, and synced first
/// unless created with [`OutputFile::unsynced`].
///
/// Created with [`OutputFile::create`] or wrapped around an open [`File`] with
/// [`From`]. Writes go straight to the file (wrap it in a [`BufWriter`] for
/// small writes). Finish with [`OutputFile::close`], or [`close_buffered`] when
/// it sits under a `BufWriter`.
#[derive(Debug)]
pub struct OutputFile {
    file: File,
    sync: bool,
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

    /// Wrap `file` so that [`close`](Self::close) checks the close but does
    /// not sync: for a descriptor whose data is synced elsewhere, or that this
    /// process does not own (a duplicate of stdout, whose redirect target is
    /// the shell's file).
    #[must_use]
    pub fn unsynced(file: File) -> Self {
        Self { file, sync: false }
    }

    /// Sync the file's data to storage if it is a regular file, then close it,
    /// returning any error either step reports.
    ///
    /// Syncing is skipped for pipes, character devices (e.g. `/dev/null`), and
    /// other non-regular files, and for a file created with
    /// [`unsynced`](Self::unsynced). A sync that fails with `EINVAL`,
    /// `ENOTSUP`, `EOPNOTSUPP` or `ENOTTY` (the file or filesystem does not
    /// support it, e.g. `F_FULLFSYNC` on some network mounts on macOS) is
    /// logged and skipped, as htslib does; any other sync error is returned. On Unix an
    /// `EINTR` from `close` is treated as success: the descriptor is released
    /// regardless, and retrying could close a descriptor another thread has
    /// since opened.
    ///
    /// The file is closed even when the sync fails; the sync's error is then
    /// the one returned.
    ///
    /// # Errors
    ///
    /// Returns an error if the sync or the close fails.
    pub fn close(self) -> io::Result<()> {
        if self.sync { sync_then_close(self.file, sync_if_regular) } else { close_file(self.file) }
    }
}

impl From<File> for OutputFile {
    fn from(file: File) -> Self {
        Self { file, sync: true }
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
/// redirect onto a full or failing mount reports the error there) but it is
/// not synced: stdout is usually a pipe, and a redirect target is the shell's
/// file, not ours. The duplicate is closed, never fd 1 itself.
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
    Ok(Box::new(OutputFile::unsynced(File::from(dup))))
}

/// Stdout itself: targets without file descriptors cannot duplicate it, so the
/// line-buffering concern on the Unix variant cannot be worked around there.
#[cfg(not(unix))]
#[expect(clippy::unnecessary_wraps, reason = "matches the fallible Unix variant")]
fn block_buffered_stdout() -> io::Result<Box<dyn OutputSink>> {
    Ok(Box::new(io::stdout()))
}

/// Finish a temp file and atomically rename it onto `dest`: `close` releases
/// the temp's own handle (e.g. `|f| OutputFile::from(f).close()` to sync and
/// close it), and only once that succeeds is the temp renamed.
///
/// So a write-back or close error fails the write instead of leaving a renamed
/// but incomplete file at `dest`. On any error the temp is removed and `dest`
/// is left untouched. Re-stamp the temp's mode (`restamp_for_persist`) before
/// calling this.
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
        format!("Failed to sync/close temp file before renaming onto: {}", dest.display())
    })?;
    temp_path
        .persist(dest)
        .map_err(|e| e.error)
        .with_context(|| format!("Failed to rename temp file onto: {}", dest.display()))
}

/// Run `sync` on `file`, then close it even if the sync failed, returning the
/// sync's error first.
fn sync_then_close<F>(file: File, sync: F) -> io::Result<()>
where
    F: FnOnce(&File) -> io::Result<bool>,
{
    let synced = sync(&file);
    let closed = close_file(file);
    synced.and(closed)
}

/// Sync `file`'s data to storage if it is a regular file. Returns whether a
/// sync was performed.
fn sync_if_regular(file: &File) -> io::Result<bool> {
    sync_regular_with(file, File::sync_data)
}

/// [`sync_if_regular`] with the sync call supplied, so tests can make it fail.
fn sync_regular_with<F>(file: &File, sync: F) -> io::Result<bool>
where
    F: FnOnce(&File) -> io::Result<()>,
{
    if !file.metadata()?.file_type().is_file() {
        return Ok(false);
    }
    match sync(file) {
        Ok(()) => Ok(true),
        Err(e) if sync_unsupported(&e) => {
            log::debug!("output file does not support sync; skipping it: {e}");
            Ok(false)
        }
        Err(e) => Err(e),
    }
}

/// Whether a sync error means the file or filesystem cannot be synced, rather
/// than that the data failed to reach storage.
///
/// Matches htslib's `fd_flush` (`hfile.c`), which `bgzf_close` reaches through
/// `hflush`: it ignores `EINVAL` (e.g. a pipe) and `ENOTSUP` ("operation-not-
/// supported errors (Mac OS X)") from `fdatasync`/`fsync` and fails on any
/// other error. `EOPNOTSUPP` is included because it is a distinct value from
/// `ENOTSUP` on macOS and the BSDs (Linux defines them equal). On macOS std's
/// `sync_data` is `fcntl(F_FULLFSYNC)`, which filesystems without it reject
/// with `ENOTSUP` (e.g. SMB mounts) or, where the filesystem has no handler
/// for the request at all, `ENOTTY`; htslib does not see `ENOTTY` because it
/// calls plain `fsync`, so it is added here.
#[cfg(unix)]
fn sync_unsupported(e: &io::Error) -> bool {
    use nix::errno::Errno;

    let unsupported = [Errno::EINVAL, Errno::ENOTSUP, Errno::EOPNOTSUPP, Errno::ENOTTY];
    e.raw_os_error().is_some_and(|code| unsupported.iter().any(|&errno| errno as i32 == code))
}

/// Whether a sync error means the file cannot be synced. Targets without
/// errno values report it only as [`io::ErrorKind::Unsupported`].
#[cfg(not(unix))]
fn sync_unsupported(e: &io::Error) -> bool {
    e.kind() == io::ErrorKind::Unsupported
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
    /// can fail its writes or its close (modelling write-back, or an NFS
    /// flush, failing only when the file is closed). Clones share state.
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

    #[test]
    fn test_sync_if_regular_syncs_only_regular_files() {
        let tmp = tempfile::tempfile().unwrap();
        assert!(sync_if_regular(&tmp).unwrap(), "a regular file must be synced");

        #[cfg(unix)]
        {
            let dev_null = std::fs::OpenOptions::new().write(true).open("/dev/null").unwrap();
            assert!(!sync_if_regular(&dev_null).unwrap(), "/dev/null must not be synced");

            let (_reader, writer) = std::io::pipe().unwrap();
            let pipe = File::from(std::os::fd::OwnedFd::from(writer));
            assert!(!sync_if_regular(&pipe).unwrap(), "a pipe must not be synced");
        }
    }

    /// The sync errors htslib treats as "sync not supported" skip the sync;
    /// every other sync error (EIO, ENOSPC, EDQUOT, ...) is fatal.
    #[cfg(unix)]
    #[test]
    fn test_sync_unsupported_classifies_errnos() {
        use nix::errno::Errno;

        for errno in [Errno::EINVAL, Errno::ENOTSUP, Errno::EOPNOTSUPP, Errno::ENOTTY] {
            let err = io::Error::from_raw_os_error(errno as i32);
            assert!(sync_unsupported(&err), "{errno} means sync is unsupported");
        }
        for errno in [Errno::EIO, Errno::ENOSPC, Errno::EDQUOT, Errno::EBADF] {
            let err = io::Error::from_raw_os_error(errno as i32);
            assert!(!sync_unsupported(&err), "{errno} must stay fatal");
        }
        assert!(!sync_unsupported(&io::Error::other("no errno")));
    }

    /// A regular file whose sync is unsupported is reported as not synced, and
    /// closing it still succeeds with its bytes on disk; a real sync failure on
    /// a regular file still fails the close.
    #[cfg(unix)]
    #[test]
    fn test_unsupported_sync_is_skipped_and_close_still_checked() {
        use nix::errno::Errno;

        let dir = tempfile::tempdir().unwrap();
        for errno in [Errno::EINVAL, Errno::ENOTSUP, Errno::EOPNOTSUPP, Errno::ENOTTY] {
            let path = dir.path().join(format!("{errno}.txt"));
            let mut file = File::create(&path).unwrap();
            file.write_all(b"data").unwrap();
            let fail = |_: &File| Err(io::Error::from_raw_os_error(errno as i32));
            assert!(!sync_regular_with(&file, fail).unwrap(), "{errno}: not synced");
            sync_then_close(file, |f| sync_regular_with(f, fail))
                .unwrap_or_else(|e| panic!("{errno}: an unsupported sync must not fail: {e}"));
            assert_eq!(std::fs::read(&path).unwrap(), b"data");
        }

        let file = tempfile::tempfile().unwrap();
        let eio = |_: &File| Err(io::Error::from_raw_os_error(Errno::EIO as i32));
        let err = sync_then_close(file, |f| sync_regular_with(f, eio))
            .expect_err("EIO from sync must surface");
        assert_eq!(err.raw_os_error(), Some(Errno::EIO as i32));
    }

    /// Closing an output that is not a regular file must still succeed: syncing
    /// a pipe or `/dev/null` fails (`EINVAL`, or `ENOTSUP`/`ENOTTY` on macOS).
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

    /// The same failure through the public API, synced or not: an
    /// `OutputFile` whose descriptor is no longer valid must not close cleanly.
    #[cfg(unix)]
    #[test]
    fn test_output_file_close_surfaces_error() {
        use std::os::fd::AsRawFd;

        if !runs_in_own_process() {
            return;
        }
        let dir = tempfile::tempdir().unwrap();
        for (name, wrap) in [
            ("synced", OutputFile::from as fn(File) -> OutputFile),
            ("unsynced", OutputFile::unsynced),
        ] {
            let file = File::create(dir.path().join(name)).unwrap();
            let fd = file.as_raw_fd();
            let out = wrap(file);
            nix::unistd::close(fd).unwrap();
            let err = out.close().expect_err("closing a closed descriptor must fail");
            assert_eq!(err.raw_os_error(), Some(nix::errno::Errno::EBADF as i32), "{name}: {err}");
        }
    }

    /// A failed sync is reported even though the close that follows succeeds,
    /// and the file is still closed (its bytes are on disk).
    #[test]
    fn test_sync_error_surfaces_when_close_succeeds() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.txt");
        let mut file = File::create(&path).unwrap();
        file.write_all(b"data").unwrap();
        let err = sync_then_close(file, |_| Err(io::Error::from_raw_os_error(5)))
            .expect_err("a sync error must surface");
        assert_eq!(err.raw_os_error(), Some(5));
        assert_eq!(std::fs::read(&path).unwrap(), b"data");
    }

    /// An unsynced file (the stdout duplicate) is flushed and its descriptor
    /// released on close; a pipe stands in for stdout so the bytes can be read
    /// back.
    #[cfg(unix)]
    #[test]
    fn test_unsynced_close_flushes_and_closes() {
        let (mut reader, writer) = std::io::pipe().unwrap();
        let mut sink: Box<dyn OutputSink> =
            Box::new(OutputFile::unsynced(File::from(std::os::fd::OwnedFd::from(writer))));
        sink.write_all(b"to stdout").unwrap();
        sink.close().expect("closing the stdout duplicate must succeed");
        // EOF is seen only once the duplicate (the pipe's last writer) is closed.
        let mut received = Vec::new();
        reader.read_to_end(&mut received).unwrap();
        assert_eq!(received, b"to stdout");
    }

    /// An unsynced file is closed without syncing it.
    #[test]
    fn test_unsynced_output_file_skips_sync() {
        let out = OutputFile::unsynced(tempfile::tempfile().unwrap());
        assert!(!out.sync);
        out.close().unwrap();
        assert!(OutputFile::from(tempfile::tempfile().unwrap()).sync);
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
        assert!(format!("{err:#}").contains("sync/close"), "{err:#}");
        assert_eq!(std::fs::read(&dest).unwrap(), b"previous");
        assert!(!temp_path.exists(), "the temp must be removed");
    }
}
