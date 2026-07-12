//! Utilities for writing metrics files.
//!
//! This module provides convenience functions for writing metrics to TSV files
//! with consistent error handling.

use anyhow::{Context, Result};
use fgoxide::io::DelimFile;
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::Metric;

/// Write metrics to a TSV file with consistent error handling.
///
/// This is a convenience wrapper around `DelimFile::write_tsv` that provides
/// consistent error messages across all commands, and — unlike a bare
/// `DelimFile::write_tsv` — **always writes the column header**, even for an empty
/// slice.
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

/// Serializes `metrics` to `path` atomically, always emitting the column header.
///
/// Both the populated and empty cases write to a temporary file in the destination directory
/// and are then atomically renamed into place. A crash mid-write therefore never leaves the
/// target holding a partially-written row set — or, for the empty case, the bogus all-default
/// data row used to materialize the header (which would look like a real collected metric);
/// the target is either absent or holds the complete, correct content.
///
/// An empty slice yields a header-only file: the csv writer emits its header lazily (on the
/// first record), so a zero-row slice would otherwise produce a 0-byte file, which fgbio's
/// `Metric.read` rejects with "No header found". We serialize a single `T::default()` row to
/// materialize the header, then truncate to just that header line (it parses back to zero
/// rows), so an empty metrics output stays re-readable by both fgbio and fgumi.
///
/// The temp file is created via `File::create` (through `Builder::make_in`) rather than the
/// default `NamedTempFile`, so it inherits the umask-applied `0o666 & !umask` mode. Otherwise
/// `NamedTempFile`'s hard-coded owner-only `0o600` would make every metrics file owner-only.
fn write_metrics_atomic<P: AsRef<Path>, T: Serialize + Default>(
    path: P,
    metrics: &[T],
) -> Result<()> {
    let path_ref = path.as_ref();
    // Rename is only atomic within a filesystem, so create the temp file next to the target.
    let dir = path_ref.parent().unwrap_or_else(|| Path::new("."));
    let tmp = tempfile::Builder::new().make_in(dir, |p| std::fs::File::create(p))?;
    if metrics.is_empty() {
        DelimFile::default().write_tsv(tmp.path(), [T::default()])?;
        let content = std::fs::read_to_string(tmp.path())?;
        let header = content.lines().next().unwrap_or_default();
        std::fs::write(tmp.path(), format!("{header}\n"))?;
    } else {
        DelimFile::default().write_tsv(tmp.path(), metrics)?;
    }
    tmp.persist(path_ref)?;
    Ok(())
}

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
