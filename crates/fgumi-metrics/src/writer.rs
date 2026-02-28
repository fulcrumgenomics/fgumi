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
/// consistent error messages across all commands.
///
/// # Arguments
/// * `path` - Path to the output TSV file
/// * `metrics` - The metrics to write (must implement Serialize)
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
/// #[derive(Serialize)]
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
pub fn write_metrics<P: AsRef<Path>, T: Serialize>(
    path: P,
    metrics: &[T],
    description: &str,
) -> Result<()> {
    let path_ref = path.as_ref();
    DelimFile::default()
        .write_tsv(path_ref, metrics)
        .with_context(|| format!("Failed to write {} metrics: {}", description, path_ref.display()))
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
    fn test_write_metrics_empty() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let metrics: Vec<TestMetrics> = vec![];

        write_metrics(temp_file.path(), &metrics, "empty")?;

        // Verify the file was created (even if empty/header-only)
        assert!(temp_file.path().exists());

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
        let result: Result<Vec<TestMetrics>> =
            read_metrics("/nonexistent/path/file.tsv", "test");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Failed to read test metrics"));
    }

    #[test]
    fn test_read_metrics_auto_roundtrip() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let original = vec![TestMetrics { name: "auto".to_string(), count: 42, value: 3.14 }];

        write_metrics_auto(temp_file.path(), &original)?;
        let read_back: Vec<TestMetrics> = read_metrics_auto(temp_file.path())?;

        assert_eq!(original, read_back);

        Ok(())
    }
}
