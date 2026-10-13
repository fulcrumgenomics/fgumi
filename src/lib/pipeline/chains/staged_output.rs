//! A sort chain's output staged beside its destination, published only when
//! the run succeeded and wrote every record it read.
//!
//! `fgumi sort`, `fgumi_lib::sort::sort_bam` and `fgumi simulate` all write
//! through it: the chain's sink writes a temp
//! ([`fgumi_sort::MergeOutputTarget`]: stdout streamed directly, a symlinked
//! destination followed, a non-regular destination refused, the file's mode
//! re-stamped on publish), and the BAI the sink derives beside it. A failed
//! run, or one whose record counts disagree, publishes nothing: the
//! destination (and its index) keep whatever they held.
//!
//! Publishing an indexed output (`--write-index`) takes three steps: the
//! destination's old BAI is removed, the BAM is renamed into place, then the
//! BAI. A failure (or a crash) between them leaves a BAM with no index, never a
//! BAM beside an index of another BAM, whose offsets would answer region
//! queries with the wrong records and no error. The price is that a failure to
//! publish the BAM, after the old index was removed, leaves the old BAM without
//! its index.
//!
//! That guarantee holds for indexed runs only. A run without an index publishes
//! the BAM and does not touch an existing sidecar, as `samtools sort` does: an
//! earlier run's `.bai` beside the destination is left in place, and it is
//! stale once the new BAM replaces the one it indexed.

use std::path::{Path, PathBuf};

use anyhow::{Result, bail, ensure};

use super::SinkSpec;

/// A staged sort output (see the module docs).
pub(crate) struct StagedSortOutput {
    target: fgumi_sort::MergeOutputTarget,
    output: PathBuf,
    /// The staged BAI the sink writes, removed unless published.
    staged_index: Option<RemoveOnDrop>,
}

/// A file removed on drop unless disarmed.
struct RemoveOnDrop(Option<PathBuf>);

impl Drop for RemoveOnDrop {
    fn drop(&mut self) {
        if let Some(path) = self.0.take() {
            let _ = std::fs::remove_file(path);
        }
    }
}

impl StagedSortOutput {
    /// Stage `command`'s output for `output`, with its BAI when `write_index`.
    ///
    /// # Errors
    ///
    /// An index requested for stdout, or a destination that cannot be staged.
    pub(crate) fn create(output: &Path, write_index: bool, command: &str) -> Result<Self> {
        if write_index && fgumi_bam_io::is_stdout_path(output) {
            bail!("a BAI index needs a file output, not stdout");
        }
        let target = fgumi_sort::MergeOutputTarget::create_for(output, command)?;
        let staged_index =
            write_index.then(|| RemoveOnDrop(Some(fgumi_bam_io::bai_sidecar_path(target.path()))));
        Ok(Self { target, output: output.to_path_buf(), staged_index })
    }

    /// The sink the chain writes: the staged BAM, with its BAI when indexing.
    pub(crate) fn sink(&self) -> SinkSpec {
        let path = self.target.path().to_path_buf();
        if self.staged_index.is_some() { SinkSpec::BamWithIndex(path) } else { SinkSpec::Bam(path) }
    }

    /// The staged BAM's path.
    pub(crate) fn path(&self) -> &Path {
        self.target.path()
    }

    /// Publish the output after a successful run whose `stats` record every
    /// record read as written (the sort's permutation guard); otherwise
    /// nothing is published.
    ///
    /// # Errors
    ///
    /// The counts disagree, the destination's old index cannot be removed, or
    /// publishing the BAM or its index failed (see the module docs for what
    /// each failure leaves behind).
    pub(crate) fn finish(mut self, stats: &fgumi_sort::SortStats) -> Result<()> {
        ensure_all_written(stats, &self.output)?;
        let staged_index = self.staged_index.take();
        if staged_index.is_some() {
            // This run writes an index: never let the new BAM sit beside the
            // old BAM's index. A run without one leaves any sidecar alone (see
            // the module docs).
            let dest = fgumi_bam_io::bai_sidecar_path(&self.output);
            match std::fs::remove_file(&dest) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => bail!("failed to remove the old index {}: {e}", dest.display()),
            }
        }
        self.target.persist()?;
        if let Some(mut index) = staged_index {
            let staged = index.0.take().expect("armed until published");
            let dest = fgumi_bam_io::bai_sidecar_path(&self.output);
            if let Err(e) = std::fs::rename(&staged, &dest) {
                let _ = std::fs::remove_file(&staged);
                bail!("failed to publish the index {}: {e}", dest.display());
            }
        }
        Ok(())
    }
}

/// The sort's permutation guard: every record read was written.
fn ensure_all_written(stats: &fgumi_sort::SortStats, output: &Path) -> Result<()> {
    ensure!(
        stats.total_records == stats.output_records,
        "sort lost records: read {} but wrote {} (differ by {}); nothing was written to {}",
        stats.total_records,
        stats.output_records,
        stats.total_records.abs_diff(stats.output_records),
        output.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A count mismatch publishes nothing: the staged BAM and BAI are removed
    /// and an existing destination keeps its bytes. Catches: a short output
    /// published (the retired engine's `enforce_record_count`).
    #[test]
    fn a_record_count_mismatch_publishes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.bam");
        std::fs::write(&out, b"previous").unwrap();
        let staged = StagedSortOutput::create(&out, true, "sort").unwrap();
        let (bam, bai) =
            (staged.path().to_path_buf(), fgumi_bam_io::bai_sidecar_path(staged.path()));
        std::fs::write(&bam, b"short").unwrap();
        std::fs::write(&bai, b"index").unwrap();
        let stats = fgumi_sort::SortStats { total_records: 10, output_records: 9, runs_written: 0 };
        let err = staged.finish(&stats).unwrap_err();
        assert!(err.to_string().contains("read 10 but wrote 9"), "{err:#}");
        assert_eq!(std::fs::read(&out).unwrap(), b"previous", "the destination is untouched");
        assert!(!bam.exists() && !bai.exists(), "the staged files are removed");
        assert!(!fgumi_bam_io::bai_sidecar_path(&out).exists());
    }

    /// A good run publishes the BAM and its index under the destination's
    /// names. Catches: the index left at the staged name.
    #[test]
    fn a_complete_run_publishes_the_bam_and_its_index() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.bam");
        let staged = StagedSortOutput::create(&out, true, "sort").unwrap();
        assert!(matches!(staged.sink(), SinkSpec::BamWithIndex(_)));
        std::fs::write(staged.path(), b"bam").unwrap();
        std::fs::write(fgumi_bam_io::bai_sidecar_path(staged.path()), b"bai").unwrap();
        let stats = fgumi_sort::SortStats { total_records: 3, output_records: 3, runs_written: 0 };
        staged.finish(&stats).unwrap();
        assert_eq!(std::fs::read(&out).unwrap(), b"bam");
        assert_eq!(std::fs::read(fgumi_bam_io::bai_sidecar_path(&out)).unwrap(), b"bai");
        let leftovers: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(Result::ok)
            .filter(|e| e.file_name().to_string_lossy().starts_with(".fgumi-"))
            .collect();
        assert!(leftovers.is_empty(), "{leftovers:?}");
    }

    /// An index that fails to publish leaves the new BAM with no index, not
    /// beside the previous run's index (whose offsets point into the old
    /// BAM). The staged BAI is missing here, so its rename fails.
    #[test]
    fn a_failed_index_publish_leaves_no_stale_index() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.bam");
        let dest_bai = fgumi_bam_io::bai_sidecar_path(&out);
        std::fs::write(&out, b"old bam").unwrap();
        std::fs::write(&dest_bai, b"old bai").unwrap();
        let staged = StagedSortOutput::create(&out, true, "sort").unwrap();
        std::fs::write(staged.path(), b"new bam").unwrap();
        let stats = fgumi_sort::SortStats { total_records: 1, output_records: 1, runs_written: 0 };
        let err = staged.finish(&stats).unwrap_err();
        assert!(err.to_string().contains("failed to publish the index"), "{err:#}");
        assert_eq!(std::fs::read(&out).unwrap(), b"new bam");
        assert!(!dest_bai.exists(), "the old BAM's index must not survive beside the new BAM");
    }

    /// An old index that cannot be removed (here a non-empty directory at
    /// its path) fails the publish before the BAM is touched; a missing old
    /// index is not an error (`a_complete_run_publishes_the_bam_and_its_index`).
    #[test]
    fn an_unremovable_old_index_publishes_nothing() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.bam");
        let dest_bai = fgumi_bam_io::bai_sidecar_path(&out);
        std::fs::write(&out, b"old bam").unwrap();
        std::fs::create_dir(&dest_bai).unwrap();
        std::fs::write(dest_bai.join("keep"), b"x").unwrap();
        let staged = StagedSortOutput::create(&out, true, "sort").unwrap();
        std::fs::write(staged.path(), b"new bam").unwrap();
        std::fs::write(fgumi_bam_io::bai_sidecar_path(staged.path()), b"new bai").unwrap();
        let stats = fgumi_sort::SortStats { total_records: 1, output_records: 1, runs_written: 0 };
        let err = staged.finish(&stats).unwrap_err();
        assert!(err.to_string().contains("failed to remove the old index"), "{err:#}");
        assert_eq!(std::fs::read(&out).unwrap(), b"old bam", "the BAM is not replaced");
    }
}
