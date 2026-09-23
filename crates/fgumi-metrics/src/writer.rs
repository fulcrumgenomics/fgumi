//! Utilities for writing metrics files.
//!
//! This module provides convenience functions for writing metrics to TSV files
//! with consistent error handling.

use anyhow::{Context, Result};
use fgoxide::io::{DelimFile, Io};
use flate2::Compression;
use flate2::write::GzEncoder;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::Metric;

/// Write metrics to a TSV file with consistent error handling.
///
/// Writes the same tab-delimited format as `DelimFile::write_tsv`, with consistent
/// error messages across all commands, and — unlike a bare `DelimFile::write_tsv` —
/// **always writes the column header**, even for an empty slice.
///
/// The underlying csv writer emits headers lazily (on the first record), so writing
/// a zero-row slice through it produces a 0-byte file. fgbio writes the header
/// eagerly in its `Metric` writer, so fgbio's `Metric.read` rejects such a file with
/// "No header found". Emitting a header-only file keeps an empty metrics output
/// re-readable by both fgbio and fgumi (it parses back to zero rows).
///
/// # Arguments
/// * `path` - Path to the output TSV file
/// * `metrics` - The metrics to write (must implement `Serialize` and `Default`)
/// * `description` - Human-readable description of the metrics for error messages
///
/// # Errors
/// Returns an error if the file cannot be created or written to
///
/// # Example
/// ```no_run
/// use fgumi_metrics::writer::write_metrics;
/// use serde::Serialize;
/// use std::path::Path;
///
/// #[derive(Serialize, Default)]
/// struct MyMetrics {
///     count: usize,
///     value: f64,
/// }
///
/// let metrics = vec![
///     MyMetrics { count: 10, value: 1.5 },
///     MyMetrics { count: 20, value: 2.5 },
/// ];
///
/// write_metrics(Path::new("metrics.txt"), &metrics, "processing").unwrap();
/// ```
pub fn write_metrics<P: AsRef<Path>, T: Serialize + Default>(
    path: P,
    metrics: &[T],
    description: &str,
) -> Result<()> {
    let path_ref = path.as_ref();
    write_metrics_atomic::<_, T>(path_ref, metrics)
        .with_context(|| format!("Failed to write {} metrics: {}", description, path_ref.display()))
}

/// Serializes `metrics` to `path`, always emitting the column header.
///
/// Regular-file destinations are written atomically: the rows go to a temporary file in
/// the destination directory, which is fsynced and then renamed into place. A crash
/// mid-write therefore never leaves the target holding a partial row set; it is either
/// absent or complete. The rows are written through the temp file's own handle, never by
/// re-opening its path. A `.gz`/`.bgz` destination is written gzip-compressed (compression
/// is chosen by the destination's extension). A symlinked destination is resolved first, so
/// the file it points at is replaced and the link itself is preserved.
///
/// Non-regular destinations (`/dev/stdout`, a FIFO, a `>(...)` process substitution) cannot
/// be the target of a rename, so they are written in place. So are file-descriptor aliases
/// (`/dev/stdout`, `/dev/fd/N`, `/proc/<pid>/fd/N`, reached directly or through a symlink)
/// even when they resolve to a regular file, e.g. `--metrics /dev/stdout >> run.log`: renaming over `run.log` would
/// discard its earlier contents and strand later writes on the unlinked inode. In-place
/// writes open in append mode, never truncating: on Linux, opening `/dev/stdout` re-opens
/// the file behind it with the caller's flags, so a truncating open would wipe `run.log`.
///
/// An empty slice yields a header-only file: the csv writer emits its header lazily (on the
/// first record), so a zero-row slice would otherwise produce a 0-byte file, which fgbio's
/// `Metric.read` rejects with "No header found".
///
/// The temp file is created via [`create_new_file`] (through `Builder::make_in`) rather than
/// the default `NamedTempFile`, so it inherits the umask-applied `0o666 & !umask` mode.
/// Otherwise `NamedTempFile`'s hard-coded owner-only `0o600` would make every metrics file
/// owner-only.
fn write_metrics_atomic<P: AsRef<Path>, T: Serialize + Default>(
    path: P,
    metrics: &[T],
) -> Result<()> {
    let path = path.as_ref();
    if is_non_regular_file(path) || routes_through_kernel_path(path) {
        let file = std::fs::OpenOptions::new().append(true).create(true).open(path)?;
        return write_metrics_to(file, Io::is_gzip_path(path), metrics);
    }
    let target = resolve_symlink(path);
    // Rename is only atomic within a filesystem, so create the temp file next to the target.
    let dir = target.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
    let tmp = tempfile::Builder::new().make_in(dir, create_new_file)?;
    write_metrics_to(tmp.as_file().try_clone()?, Io::is_gzip_path(&target), metrics)?;
    // Flush the rows to disk before the rename; otherwise a system crash can persist the
    // rename ahead of the data and leave a truncated (yet parseable) file at the target.
    tmp.as_file().sync_all()?;
    tmp.persist(&target)?;
    Ok(())
}

/// Creates `path` exclusively (`O_CREAT | O_EXCL`) with the default `0o666 & !umask` mode.
///
/// Unlike `File::create`, this fails with `AlreadyExists` when `path` exists — as a regular
/// file or as a symlink, which is never followed — so `Builder::make_in` retries with a fresh
/// name instead of truncating or writing through a path another process created first.
fn create_new_file(path: &Path) -> std::io::Result<std::fs::File> {
    std::fs::OpenOptions::new().write(true).create_new(true).open(path)
}

/// gzip level for `.gz`/`.bgz` metrics, matching `fgoxide::io::Io::default()`.
const GZIP_COMPRESSION_LEVEL: u32 = 5;

/// Writes `metrics` (or, when empty, just the header row) to `file`, gzip-compressed when
/// `gzip` is set. Every buffer and the gzip trailer are flushed explicitly, so a failed
/// final write is reported rather than lost in a `Drop`.
fn write_metrics_to<T: Serialize + Default>(file: File, gzip: bool, metrics: &[T]) -> Result<()> {
    let out = BufWriter::new(file);
    let out = if gzip {
        let encoder = GzEncoder::new(out, Compression::new(GZIP_COMPRESSION_LEVEL));
        write_rows(encoder, metrics)?.finish()?
    } else {
        write_rows(out, metrics)?
    };
    out.into_inner().map_err(std::io::IntoInnerError::into_error)?;
    Ok(())
}

/// Serializes `metrics` (or, when empty, just the header row) to `sink` and returns it.
fn write_rows<W: Write, T: Serialize + Default>(mut sink: W, metrics: &[T]) -> Result<W> {
    if metrics.is_empty() {
        writeln!(sink, "{}", header_of::<T>()?)?;
        return Ok(sink);
    }
    let mut writer = tsv_writer(sink);
    for metric in metrics {
        writer.serialize(metric)?;
    }
    writer.into_inner().map_err(|e| anyhow::anyhow!("failed to flush TSV writer: {}", e.error()))
}

/// A csv writer configured like `DelimFile::write_tsv`: tab-delimited, a header row, and
/// quoting only where necessary.
fn tsv_writer<W: Write>(sink: W) -> csv::Writer<W> {
    csv::WriterBuilder::new()
        .delimiter(b'\t')
        .quote_style(csv::QuoteStyle::Necessary)
        .from_writer(sink)
}

/// The TSV header row for `T`, obtained by serializing one default row in memory (the csv
/// writer only knows the column names once it sees a record).
fn header_of<T: Serialize + Default>() -> Result<String> {
    let mut writer = tsv_writer(Vec::new());
    writer.serialize(T::default())?;
    let bytes = writer
        .into_inner()
        .map_err(|e| anyhow::anyhow!("failed to flush TSV writer: {}", e.error()))?;
    Ok(String::from_utf8(bytes)?.lines().next().unwrap_or_default().to_string())
}

/// True when `path` exists and (following symlinks) is neither a regular file nor a
/// directory, e.g. a character device, FIFO, or `/dev/stdout`.
fn is_non_regular_file(path: &Path) -> bool {
    std::fs::metadata(path).is_ok_and(|m| !m.file_type().is_file() && !m.file_type().is_dir())
}

/// True for file-descriptor aliases — `/dev/stdout`, `/dev/stderr`, `/dev/stdin`,
/// `/dev/fd/N` and `/proc/<pid>/fd/N` — which must be opened, never renamed over, even when
/// they currently resolve to a regular file.
///
/// Other paths under `/dev` or `/proc` are not matched: a regular file there (e.g. on the
/// `/dev/shm` tmpfs) is replaced atomically like any other, and a real device or FIFO is
/// already written in place via [`is_non_regular_file`].
fn is_kernel_special_path(path: &Path) -> bool {
    let is_fd_dir = |dir: &Path| {
        dir == Path::new("/dev/fd")
            || (dir.starts_with("/proc") && dir.file_name().is_some_and(|name| name == "fd"))
    };
    matches!(path.to_str(), Some("/dev/stdout" | "/dev/stderr" | "/dev/stdin"))
        || path.parent().is_some_and(is_fd_dir)
}

/// True when `path`, or any hop of the symlink chain starting at it, is a kernel special
/// path (see [`is_kernel_special_path`]). The walk is bounded like [`resolve_symlink`]'s.
fn routes_through_kernel_path(path: &Path) -> bool {
    let mut current = path.to_path_buf();
    for _ in 0..MAX_SYMLINK_HOPS {
        if is_kernel_special_path(&current) {
            return true;
        }
        let Ok(next) = std::fs::read_link(&current) else { return false };
        current = current.parent().unwrap_or(Path::new(".")).join(next);
    }
    false
}

/// Resolves a symlinked `path` to the file it points at, so the atomic rename replaces
/// that file rather than the link. Non-symlinks, unreadable links, and link loops are
/// returned as-is.
///
/// A dangling link (its target does not exist yet) cannot be canonicalized, so it is
/// followed hop by hop with `read_link` instead; relative targets resolve against the
/// link's own directory. The write then creates the missing target, as `File::create` would.
fn resolve_symlink(path: &Path) -> std::path::PathBuf {
    let is_symlink =
        |p: &Path| std::fs::symlink_metadata(p).is_ok_and(|m| m.file_type().is_symlink());
    if !is_symlink(path) {
        return path.to_path_buf();
    }
    if let Ok(resolved) = std::fs::canonicalize(path) {
        return resolved;
    }
    let mut current = path.to_path_buf();
    for _ in 0..MAX_SYMLINK_HOPS {
        if !is_symlink(&current) {
            return current;
        }
        let Ok(next) = std::fs::read_link(&current) else { break };
        current = current.parent().unwrap_or(Path::new(".")).join(next);
    }
    path.to_path_buf()
}

/// Matches Linux's `MAXSYMLINKS`; bounds symlink walks if the links form a loop.
const MAX_SYMLINK_HOPS: usize = 40;

/// Write metrics implementing the Metric trait to a TSV file.
///
/// This version uses the metric's own name for error messages, providing
/// a more concise API when the metrics type is known at compile time.
///
/// # Arguments
/// * `path` - Path to the output TSV file
/// * `metrics` - The metrics to write (must implement Metric)
///
/// # Errors
/// Returns an error if the file cannot be created or written to
///
/// # Example
/// ```no_run
/// use fgumi_metrics::writer::write_metrics_auto;
/// use fgumi_metrics::consensus::ConsensusMetrics;
/// use std::path::Path;
///
/// let metrics = vec![ConsensusMetrics::default()];
/// write_metrics_auto(Path::new("metrics.txt"), &metrics).unwrap();
/// ```
pub fn write_metrics_auto<P: AsRef<Path>, T: Metric>(path: P, metrics: &[T]) -> Result<()> {
    write_metrics(path, metrics, T::metric_name())
}

/// Read metrics from a TSV file with consistent error handling.
///
/// # Arguments
/// * `path` - Path to the TSV file
/// * `description` - Human-readable description for error messages
///
/// # Errors
/// Returns an error if the file cannot be read or parsed
pub fn read_metrics<P: AsRef<Path>, T: for<'de> Deserialize<'de>>(
    path: P,
    description: &str,
) -> Result<Vec<T>> {
    let path_ref = path.as_ref();
    DelimFile::default()
        .read_tsv(path_ref)
        .with_context(|| format!("Failed to read {} metrics: {}", description, path_ref.display()))
}

/// Read metrics implementing the Metric trait from a TSV file.
///
/// Uses the metric's own name for error messages.
///
/// # Arguments
/// * `path` - Path to the TSV file
///
/// # Errors
/// Returns an error if the file cannot be read or parsed
pub fn read_metrics_auto<P: AsRef<Path>, T: Metric>(path: P) -> Result<Vec<T>> {
    read_metrics(path, T::metric_name())
}

/// Test-only guard: writing `rows`, reading them back, and writing them again
/// must produce byte-identical output, so every value survives a trip through the
/// metrics TSV reader (a field the reader drops, zeroes or truncates re-serializes
/// differently). No `PartialEq` is needed on the struct.
///
/// Pass rows with non-default values (non-zero counts, fractional and non-finite
/// floats, `Some` options): an all-default row can only prove that zeros survive.
/// This does not check which *tokens* are written (e.g. `Infinity` vs `inf`); the
/// workspace float-encoding contract test guards that.
#[cfg(test)]
pub(crate) fn assert_roundtrip_stable<T: Metric>(rows: &[T]) {
    use tempfile::NamedTempFile;

    let first = NamedTempFile::new().expect("temp file");
    write_metrics_auto(first.path(), rows).expect("write first");
    let back: Vec<T> = read_metrics_auto(first.path()).expect("read back");
    let second = NamedTempFile::new().expect("temp file");
    write_metrics_auto(second.path(), &back).expect("write second");

    assert_eq!(
        std::fs::read_to_string(first.path()).expect("read first"),
        std::fs::read_to_string(second.path()).expect("read second"),
        "round-trip not stable for metric {}",
        T::metric_name(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use serde::Deserialize;
    use std::fs;
    use tempfile::NamedTempFile;

    #[derive(Debug, Serialize, Deserialize, PartialEq, Clone, Default)]
    struct TestMetrics {
        name: String,
        count: usize,
        value: f64,
    }

    impl Metric for TestMetrics {
        fn metric_name() -> &'static str {
            "test"
        }
    }

    #[test]
    fn test_write_metrics_success() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let metrics = vec![
            TestMetrics { name: "test1".to_string(), count: 10, value: 1.5 },
            TestMetrics { name: "test2".to_string(), count: 20, value: 2.5 },
        ];

        write_metrics(temp_file.path(), &metrics, "test")?;

        // Verify the file was written
        let content = fs::read_to_string(temp_file.path())?;
        assert!(content.contains("name"));
        assert!(content.contains("count"));
        assert!(content.contains("value"));
        assert!(content.contains("test1"));
        assert!(content.contains("test2"));

        Ok(())
    }

    #[test]
    fn test_write_metrics_invalid_path() {
        let metrics = vec![TestMetrics { name: "test".to_string(), count: 10, value: 1.5 }];

        let result = write_metrics("/invalid/path/metrics.txt", &metrics, "test");
        assert!(result.is_err());
        if let Err(e) = result {
            let err_msg = e.to_string();
            assert!(err_msg.contains("Failed to write test metrics"));
        }
    }

    #[test]
    fn test_write_metrics_empty_is_header_only_and_round_trips() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let metrics: Vec<TestMetrics> = vec![];

        write_metrics(temp_file.path(), &metrics, "empty")?;

        // An empty metrics slice must still produce the column header (fgbio's
        // Metric.read rejects a 0-byte file with "No header found"), and no data rows.
        let content = fs::read_to_string(temp_file.path())?;
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 1, "expected header-only file, got: {content:?}");
        assert_eq!(lines[0], "name\tcount\tvalue");

        // And it round-trips back to zero rows.
        let read_back: Vec<TestMetrics> = read_metrics(temp_file.path(), "empty")?;
        assert!(read_back.is_empty());

        Ok(())
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::populated(vec![TestMetrics { name: "a".to_string(), count: 1, value: 1.0 }])]
    fn test_write_metrics_is_atomic_leaves_no_temp_files(
        #[case] metrics: Vec<TestMetrics>,
    ) -> Result<()> {
        // Both the empty (header-only) and populated write paths go through a temp file that
        // must be atomically renamed into place (and never left behind): after a successful
        // write the destination directory should hold exactly the target file, with no stray
        // scratch files. This pins the atomicity of the populated path, not just the empty one.
        let dir = tempfile::tempdir()?;
        let target = dir.path().join("metrics.txt");

        write_metrics(&target, &metrics, "atomic")?;

        let entries: Vec<_> = fs::read_dir(dir.path())?.collect::<std::io::Result<_>>()?;
        assert_eq!(
            entries.len(),
            1,
            "expected only the target file, found: {:?}",
            entries.iter().map(std::fs::DirEntry::path).collect::<Vec<_>>()
        );
        assert_eq!(entries[0].path(), target);

        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn test_write_metrics_empty_matches_populated_file_mode() -> Result<()> {
        // An empty (header-only) metrics file must not be more restrictive than a populated
        // one. The populated path writes directly via `File::create` (umask-applied mode),
        // while the header-only path routes through a temp file; if that temp is left at
        // `NamedTempFile`'s hard-coded 0600 the two outputs diverge. Assert the modes match.
        use std::os::unix::fs::PermissionsExt;

        let dir = tempfile::tempdir()?;

        let populated_path = dir.path().join("populated.txt");
        write_metrics(
            &populated_path,
            &[TestMetrics { name: "x".to_string(), count: 1, value: 1.0 }],
            "populated",
        )?;
        let populated_mode = fs::metadata(&populated_path)?.permissions().mode() & 0o777;

        let empty_path = dir.path().join("empty.txt");
        let empty: Vec<TestMetrics> = vec![];
        write_metrics(&empty_path, &empty, "empty")?;
        let empty_mode = fs::metadata(&empty_path)?.permissions().mode() & 0o777;

        assert_eq!(
            empty_mode, populated_mode,
            "empty metrics file mode {empty_mode:o} should match populated {populated_mode:o}"
        );

        Ok(())
    }

    #[test]
    fn test_roundtrip_tsv() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let original_metrics = vec![
            TestMetrics { name: "first".to_string(), count: 100, value: 12.34 },
            TestMetrics { name: "second".to_string(), count: 200, value: 56.78 },
        ];

        // Write metrics
        write_metrics(temp_file.path(), &original_metrics, "roundtrip")?;

        // Read them back using DelimFile directly
        let read_metrics: Vec<TestMetrics> = DelimFile::default().read_tsv(temp_file.path())?;

        // Verify they match
        assert_eq!(original_metrics.len(), read_metrics.len());
        for (orig, read) in original_metrics.iter().zip(read_metrics.iter()) {
            assert_eq!(orig, read);
        }

        Ok(())
    }

    /// A `.gz` destination must hold real gzip data that the metrics reader (which picks
    /// decompression by extension) can read back.
    #[test]
    fn test_write_metrics_gz_destination_is_gzip_and_round_trips() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("metrics.txt.gz");
        let rows = vec![TestMetrics { name: "a".to_string(), count: 1, value: 0.5 }];
        write_metrics(&path, &rows, "gz")?;
        assert_eq!(&fs::read(&path)?[..2], &[0x1f, 0x8b], "expected gzip magic bytes");
        let back: Vec<TestMetrics> = read_metrics(&path, "gz")?;
        assert_eq!(back, rows);
        Ok(())
    }

    /// An empty write to a `.gz` destination must be a gzipped header-only file.
    #[test]
    fn test_write_metrics_empty_gz_destination_is_gzipped_header() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("metrics.txt.gz");
        write_metrics(&path, &Vec::<TestMetrics>::new(), "gz")?;
        assert_eq!(&fs::read(&path)?[..2], &[0x1f, 0x8b], "expected gzip magic bytes");
        let lines = fgoxide::io::Io::default().read_lines(&path)?;
        assert_eq!(lines, vec!["name\tcount\tvalue".to_string()]);
        let back: Vec<TestMetrics> = read_metrics(&path, "gz")?;
        assert!(back.is_empty());
        Ok(())
    }

    /// A non-regular destination (character device, FIFO, `/dev/stdout`) is written in
    /// place: it cannot be the target of a temp-file rename.
    #[cfg(unix)]
    #[test]
    fn test_write_metrics_to_character_device_writes_in_place() -> Result<()> {
        let rows = vec![TestMetrics { name: "a".to_string(), count: 1, value: 0.5 }];
        write_metrics("/dev/null", &rows, "devnull")?;
        write_metrics("/dev/null", &Vec::<TestMetrics>::new(), "devnull")?;
        Ok(())
    }

    /// A kernel special path (`/dev/fd/N`, `/dev/stdout`) that resolves to a regular file —
    /// e.g. `--metrics /dev/stdout >> run.log` — is written in place, never renamed over, and
    /// never truncated: the log's earlier contents survive and the rows follow them. Renaming
    /// would replace the file the descriptor refers to, leaving the open handle on an
    /// unlinked inode, so the rows are read back through the held handle. On Linux, opening
    /// `/dev/fd/N` re-opens the file with the caller's flags, so a truncating open would wipe
    /// the log. With `via_symlink`, the destination is a user symlink pointing at the special
    /// path.
    #[cfg(unix)]
    #[rstest]
    fn test_write_metrics_to_fd_path_backed_by_regular_file_writes_in_place(
        #[values(false, true)] via_symlink: bool,
    ) -> Result<()> {
        use std::io::{Read, Seek, SeekFrom};
        use std::os::fd::AsRawFd;

        let dir = tempfile::tempdir()?;
        let backing = dir.path().join("run.log");
        fs::write(&backing, "earlier\n")?;
        // The shell's `>> run.log`: an append-mode descriptor on a non-empty file.
        let mut handle = fs::OpenOptions::new().read(true).append(true).open(&backing)?;
        let fd_path = std::path::PathBuf::from(format!("/dev/fd/{}", handle.as_raw_fd()));
        let dest = if via_symlink {
            let link = dir.path().join("link.txt");
            std::os::unix::fs::symlink(&fd_path, &link)?;
            link
        } else {
            fd_path
        };

        let rows = vec![TestMetrics { name: "a".to_string(), count: 1, value: 0.5 }];
        write_metrics(&dest, &rows, "fd")?;

        let mut content = String::new();
        handle.seek(SeekFrom::Start(0))?;
        handle.read_to_string(&mut content)?;
        assert_eq!(content, "earlier\nname\tcount\tvalue\na\t1\t0.5\n");
        Ok(())
    }

    /// Only descriptor aliases are kernel special paths; ordinary files that merely live
    /// under `/dev` or `/proc` (e.g. tmpfs at `/dev/shm`) are not, so they keep the atomic
    /// replace path. Real devices and FIFOs are caught separately by `is_non_regular_file`.
    #[rstest]
    #[case::dev_stdout("/dev/stdout", true)]
    #[case::dev_stderr("/dev/stderr", true)]
    #[case::dev_stdin("/dev/stdin", true)]
    #[case::dev_fd("/dev/fd/3", true)]
    #[case::proc_self_fd("/proc/self/fd/1", true)]
    #[case::proc_pid_fd("/proc/1234/fd/4", true)]
    #[case::proc_thread_self_fd("/proc/thread-self/fd/2", true)]
    #[case::dev_shm_file("/dev/shm/run1/metrics.txt", false)]
    #[case::dev_fd_dir("/dev/fd", false)]
    #[case::dev_null("/dev/null", false)]
    #[case::proc_status("/proc/self/status", false)]
    #[case::fd_outside_proc("/tmp/fd/3", false)]
    #[case::relative("metrics.txt", false)]
    fn test_is_kernel_special_path(#[case] path: &str, #[case] expected: bool) {
        assert_eq!(is_kernel_special_path(Path::new(path)), expected);
    }

    /// A regular file under `/dev/shm` is replaced atomically like any other regular file:
    /// writing it twice leaves exactly the second write, not an appended second header.
    #[cfg(target_os = "linux")]
    #[test]
    fn test_write_metrics_to_dev_shm_regular_file_replaces() -> Result<()> {
        let dir = tempfile::tempdir_in("/dev/shm")?;
        let path = dir.path().join("metrics.txt");
        let first = vec![TestMetrics { name: "a".to_string(), count: 1, value: 0.5 }];
        let second = vec![TestMetrics { name: "b".to_string(), count: 2, value: 1.5 }];
        write_metrics(&path, &first, "shm")?;
        write_metrics(&path, &second, "shm")?;
        assert_eq!(fs::read_to_string(&path)?, "name\tcount\tvalue\nb\t2\t1.5\n");
        Ok(())
    }

    /// The temp-file opener must refuse a path that already exists — a regular file or a
    /// symlink — rather than truncating or following it, so `Builder::make_in` retries with
    /// a fresh name instead of writing through a path another process substituted.
    #[cfg(unix)]
    #[rstest]
    fn test_create_new_file_refuses_existing_path(
        #[values(false, true)] as_symlink: bool,
    ) -> Result<()> {
        let dir = tempfile::tempdir()?;
        let victim = dir.path().join("victim.txt");
        fs::write(&victim, "keep")?;
        let candidate = if as_symlink {
            let link = dir.path().join("tmp.txt");
            std::os::unix::fs::symlink(&victim, &link)?;
            link
        } else {
            victim.clone()
        };
        let err = create_new_file(&candidate).expect_err("opened an existing path");
        assert_eq!(err.kind(), std::io::ErrorKind::AlreadyExists);
        assert_eq!(fs::read_to_string(&victim)?, "keep");
        Ok(())
    }

    /// Writing through a symlink updates the file it points at and leaves the link intact.
    #[cfg(unix)]
    #[test]
    fn test_write_metrics_through_symlink_keeps_link() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let real = dir.path().join("real.txt");
        let link = dir.path().join("link.txt");
        fs::write(&real, "stale")?;
        std::os::unix::fs::symlink(&real, &link)?;
        let rows = vec![TestMetrics { name: "a".to_string(), count: 1, value: 0.5 }];
        write_metrics(&link, &rows, "symlink")?;
        assert!(fs::symlink_metadata(&link)?.file_type().is_symlink(), "link was replaced");
        let back: Vec<TestMetrics> = read_metrics(&real, "symlink")?;
        assert_eq!(back, rows);
        Ok(())
    }

    /// Writing through a dangling symlink (its target does not exist yet) creates the target
    /// and leaves every link in the chain intact, as `File::create` would. `links` is the
    /// chain `link.txt -> links[0] -> links[1] -> ...`; the last entry names the missing
    /// target, and absolute entries are written with `{dir}` standing for the temp dir.
    #[cfg(unix)]
    #[rstest]
    #[case::relative_target(&["real.txt"], "real.txt")]
    #[case::absolute_target(&["{dir}/real.txt"], "real.txt")]
    #[case::relative_target_in_subdir(&["sub/real.txt"], "sub/real.txt")]
    #[case::chain_of_dangling_links(&["mid.txt", "real.txt"], "real.txt")]
    fn test_write_metrics_through_dangling_symlink_creates_target(
        #[case] links: &[&str],
        #[case] expected_target: &str,
    ) -> Result<()> {
        let dir = tempfile::tempdir()?;
        fs::create_dir(dir.path().join("sub"))?;
        let dir_str = dir.path().to_string_lossy();
        let mut next_link = dir.path().join("link.txt");
        let mut link_paths = Vec::new();
        for target in links {
            let target = target.replace("{dir}", &dir_str);
            std::os::unix::fs::symlink(&target, &next_link)?;
            link_paths.push(next_link.clone());
            next_link = dir.path().join(target);
        }
        let rows = vec![TestMetrics { name: "a".to_string(), count: 1, value: 0.5 }];
        write_metrics(dir.path().join("link.txt"), &rows, "dangling")?;
        for link in &link_paths {
            assert!(
                fs::symlink_metadata(link)?.file_type().is_symlink(),
                "link {} was replaced",
                link.display()
            );
        }
        let back: Vec<TestMetrics> = read_metrics(dir.path().join(expected_target), "dangling")?;
        assert_eq!(back, rows);
        Ok(())
    }

    #[test]
    fn test_write_metrics_auto() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let metrics = vec![TestMetrics { name: "auto".to_string(), count: 42, value: 99.5 }];

        write_metrics_auto(temp_file.path(), &metrics)?;

        let content = fs::read_to_string(temp_file.path())?;
        assert!(content.contains("auto"));
        assert!(content.contains("42"));

        Ok(())
    }

    #[test]
    fn test_read_metrics_roundtrip() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let original = vec![
            TestMetrics { name: "a".to_string(), count: 1, value: 1.1 },
            TestMetrics { name: "b".to_string(), count: 2, value: 2.2 },
        ];

        write_metrics(temp_file.path(), &original, "test")?;
        let read_back: Vec<TestMetrics> = read_metrics(temp_file.path(), "test")?;

        assert_eq!(original, read_back);

        Ok(())
    }

    #[test]
    fn test_read_metrics_invalid_path() {
        let result: Result<Vec<TestMetrics>> = read_metrics("/nonexistent/path/file.tsv", "test");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Failed to read test metrics"));
    }

    #[test]
    fn test_read_metrics_auto_roundtrip() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let original = vec![TestMetrics { name: "auto".to_string(), count: 42, value: 1.5 }];

        write_metrics_auto(temp_file.path(), &original)?;
        let read_back: Vec<TestMetrics> = read_metrics_auto(temp_file.path())?;

        assert_eq!(original, read_back);

        Ok(())
    }
}
