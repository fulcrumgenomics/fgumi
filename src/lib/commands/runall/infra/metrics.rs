//! `--metrics` TSV writer for `fgumi runall`.
//!
//! Chain runners that have wired metrics support write per-stage
//! metric rows to the `--metrics` TSV via [`MetricsWriter`]. PR 1
//! ships the writer surface; the only chain runner registered in
//! PR 1 ([`crate::commands::runall::chains::NotImplementedRunner`])
//! ignores `--metrics` because it never gets to a working chain.
//!
//! TSV schema:
//!
//! ```text
//! stage<TAB>wall_time_secs<TAB>records_in<TAB>records_out
//! ```
//!
//! Suppressed `dead_code`: PR 1 ships the writer surface but no
//! production code path constructs a row yet (the only registered
//! runner is the [`crate::commands::runall::chains::NotImplementedRunner`]
//! fallback, which ignores `--metrics`). Real chain runners in PR 2/3
//! consume these items.

#![allow(dead_code)]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

/// Per-stage row appended to the `--metrics` TSV.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MetricsRow {
    /// Short stage name (e.g. `"extract"`, `"correct"`).
    pub stage: String,
    /// Wall-clock time in seconds for this stage.
    pub wall_time_secs: f64,
    /// Records consumed by this stage.
    pub records_in: u64,
    /// Records emitted by this stage.
    pub records_out: u64,
}

/// Writer that buffers and emits `MetricsRow` values as TSV.
///
/// Construct via [`MetricsWriter::create`] (writes a header row
/// immediately) and feed rows through [`MetricsWriter::write_row`].
/// Drop the writer to flush.
pub(crate) struct MetricsWriter {
    /// Buffered output sink.
    inner: BufWriter<File>,
}

impl MetricsWriter {
    /// Create the metrics file at `path` and write the header row.
    pub(crate) fn create(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create metrics file {}", path.display()))?;
        let mut inner = BufWriter::new(file);
        writeln!(inner, "stage\twall_time_secs\trecords_in\trecords_out")
            .with_context(|| format!("Failed to write metrics header to {}", path.display()))?;
        Ok(Self { inner })
    }

    /// Append a row of stage metrics to the TSV.
    pub(crate) fn write_row(&mut self, row: &MetricsRow) -> Result<()> {
        writeln!(
            self.inner,
            "{}\t{:.6}\t{}\t{}",
            row.stage, row.wall_time_secs, row.records_in, row.records_out,
        )
        .context("Failed to write metrics row")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn writer_emits_header_and_rows() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("metrics.tsv");

        {
            let mut writer = MetricsWriter::create(&path).expect("create");
            writer
                .write_row(&MetricsRow {
                    stage: "extract".into(),
                    wall_time_secs: 1.5,
                    records_in: 10_000,
                    records_out: 10_000,
                })
                .unwrap();
            writer
                .write_row(&MetricsRow {
                    stage: "correct".into(),
                    wall_time_secs: 2.25,
                    records_in: 10_000,
                    records_out: 9_950,
                })
                .unwrap();
        }

        let contents = fs::read_to_string(&path).unwrap();
        let mut lines = contents.lines();
        assert_eq!(lines.next().unwrap(), "stage\twall_time_secs\trecords_in\trecords_out");
        assert_eq!(lines.next().unwrap(), "extract\t1.500000\t10000\t10000");
        assert_eq!(lines.next().unwrap(), "correct\t2.250000\t10000\t9950");
        assert!(lines.next().is_none());
    }

    #[test]
    fn writer_create_errors_when_directory_missing() {
        let path = std::path::Path::new("/nonexistent/dir/metrics.tsv");
        let err = match MetricsWriter::create(path) {
            Ok(_) => panic!("expected error for missing parent directory"),
            Err(e) => e,
        };
        let msg = format!("{err:#}");
        assert!(msg.contains("Failed to create metrics file"), "got: {msg}");
    }
}
