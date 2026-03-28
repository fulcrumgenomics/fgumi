//! Pipeline stage metrics and atomic output guard.

use std::path::PathBuf;

use anyhow::{Context, Result};

/// Drop guard that removes the temporary output file on error or cancellation.
///
/// Call `disarm()` after a successful rename to prevent cleanup.
pub(super) struct AtomicOutputGuard {
    pub(super) tmp_path: PathBuf,
}

impl AtomicOutputGuard {
    /// Disarm the guard so the temp file is not deleted on drop.
    pub(super) fn disarm(self) {
        std::mem::forget(self);
    }
}

impl Drop for AtomicOutputGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.tmp_path);
    }
}

/// Metrics for a single pipeline stage.
pub(super) struct StageMetric {
    pub(super) stage: String,
    pub(super) wall_time_secs: f64,
    pub(super) records_in: u64,
    pub(super) records_out: u64,
}

/// Collection of stage metrics for the entire pipeline run.
pub(super) struct PipelineMetrics {
    pub(super) stages: Vec<StageMetric>,
}

impl PipelineMetrics {
    pub(super) fn new() -> Self {
        Self { stages: Vec::new() }
    }

    pub(super) fn add(
        &mut self,
        stage: &str,
        wall_time_secs: f64,
        records_in: u64,
        records_out: u64,
    ) {
        self.stages.push(StageMetric {
            stage: stage.to_string(),
            wall_time_secs,
            records_in,
            records_out,
        });
    }

    /// Write metrics to a TSV file.
    pub(super) fn write_to_file(&self, path: &std::path::Path) -> Result<()> {
        use std::io::Write;
        let mut file = std::fs::File::create(path)
            .with_context(|| format!("Failed to create metrics file: {}", path.display()))?;
        writeln!(file, "stage\twall_time_secs\trecords_in\trecords_out")?;
        let mut total_time = 0.0;
        let mut first_in = 0u64;
        let mut last_out = 0u64;
        for (i, m) in self.stages.iter().enumerate() {
            writeln!(
                file,
                "{}\t{:.1}\t{}\t{}",
                m.stage, m.wall_time_secs, m.records_in, m.records_out
            )?;
            total_time += m.wall_time_secs;
            if i == 0 {
                first_in = m.records_in;
            }
            last_out = m.records_out;
        }
        writeln!(file, "total\t{total_time:.1}\t{first_in}\t{last_out}")?;
        Ok(())
    }
}
