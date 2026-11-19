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
use super::sink_spec::IndexKind;

/// A staged sort output (see the module docs).
pub(crate) struct StagedSortOutput {
    target: fgumi_sort::MergeOutputTarget,
    output: PathBuf,
    /// The staged index the sink writes (`.bai` or `.crai`), removed unless
    /// published; `None` when not indexing.
    staged_index: RemoveOnDrop,
}

/// A file (if any) removed on drop unless disarmed.
struct RemoveOnDrop(Option<PathBuf>);

impl RemoveOnDrop {
    /// Stop guarding the file, returning its path.
    fn disarm(&mut self) -> Option<PathBuf> {
        self.0.take()
    }
}

impl Drop for RemoveOnDrop {
    fn drop(&mut self) {
        if let Some(path) = self.0.take() {
            let _ = std::fs::remove_file(path);
        }
    }
}

impl StagedSortOutput {
    /// Stage `command`'s output for `output`, with its index (`.bai`, or
    /// `.crai` for a `*.cram` output) when `write_index`.
    ///
    /// # Errors
    ///
    /// An index requested for stdout, or a destination that cannot be staged.
    pub(crate) fn create(output: &Path, write_index: bool, command: &str) -> Result<Self> {
        if write_index && fgumi_bam_io::is_stdout_path(output) {
            let kind = IndexKind::for_path(output).name();
            bail!("a {kind} index needs a file output, not stdout");
        }
        let target = fgumi_sort::MergeOutputTarget::create_for(output, command)?;
        let staged_index =
            RemoveOnDrop(write_index.then(|| IndexKind::for_path(output).sidecar(target.path())));
        Ok(Self { target, output: output.to_path_buf(), staged_index })
    }

    /// Protect the input from the final sidecar, not just the staged one.
    pub(crate) fn validate_input(&self, input: &Path) -> Result<()> {
        if self.staged_index.0.is_some() {
            crate::commands::common::reject_writes_aliasing_inputs(
                &[(input, "--input")],
                &[(&destination_index(&self.output), "--write-index sidecar")],
            )?;
        }
        Ok(())
    }

    /// The sink the chain writes: the staged output, with its index when
    /// indexing, in the format the destination's name selects.
    pub(crate) fn sink(&self) -> SinkSpec {
        let path = self.target.path().to_path_buf();
        SinkSpec::alignment(path, &self.output, self.staged_index.0.is_some())
    }

    /// The staged output's path.
    #[cfg(test)]
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
        // The staged index stays armed in `self` (removed on any early return)
        // until it is renamed into place.
        if self.staged_index.0.is_some() {
            // This run writes an index: never let the new BAM sit beside the
            // old BAM's index. A run without one leaves any sidecar alone (see
            // the module docs).
            let dest = destination_index(&self.output);
            if let fgumi_sort::MergeOutputTarget::Temp { dest: primary, .. } = &self.target {
                ensure!(
                    !crate::commands::common::paths_refer_to_same_file(primary, &dest),
                    "the index {} and resolved output {} are the same file; choose a different path",
                    dest.display(),
                    primary.display()
                );
            }
            match std::fs::remove_file(&dest) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => bail!("failed to remove the old index {}: {e}", dest.display()),
            }
        }
        self.target.persist()?;
        if let Some(staged) = self.staged_index.disarm() {
            let dest = destination_index(&self.output);
            if let Err(e) = std::fs::rename(&staged, &dest) {
                let _ = std::fs::remove_file(&staged);
                bail!("failed to publish the index {}: {e}", dest.display());
            }
        }
        Ok(())
    }
}

/// The destination's index sidecar: `.crai` beside a CRAM, `.bai` beside a BAM.
fn destination_index(output: &Path) -> PathBuf {
    IndexKind::for_path(output).sidecar(output)
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
    use rstest::rstest;

    /// The sink's format, index and CRAM file ID follow the destination name
    /// and `--write-index`, while it writes to the staged path. Catches: the
    /// format or the CRAM file ID taken from the staged temp's name.
    #[rstest]
    #[case::bam("out.bam", false, false, false, None)]
    #[case::bam_indexed("out.bam", true, false, true, Some("out.bam.bai"))]
    #[case::cram("out.cram", false, true, false, None)]
    #[case::cram_indexed("out.cram", true, true, true, Some("out.cram.crai"))]
    fn the_sink_follows_the_destination_name(
        #[case] name: &str,
        #[case] write_index: bool,
        #[case] is_cram: bool,
        #[case] has_index: bool,
        #[case] dest_index: Option<&str>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join(name);
        let staged = StagedSortOutput::create(&out, write_index, "sort").unwrap();
        let sink = staged.sink();
        assert_eq!(sink.is_cram(), is_cram);
        assert_eq!(sink.has_index(), has_index);
        assert_eq!(sink.cram_name(), is_cram.then_some(out.as_path()));
        assert_eq!(sink.path().map(PathBuf::as_path), Some(staged.path()), "writes the stage");
        assert_ne!(staged.path(), out.as_path());
        let staged_index = staged.staged_index.0.clone();
        let sink_index = sink.index().map(|(_, sidecar)| sidecar);
        assert_eq!(staged_index, sink_index, "the staged index is the sink's");
        assert_eq!(
            write_index.then(|| destination_index(&out)),
            dest_index.map(|n| dir.path().join(n))
        );
    }

    /// Publishing must not replace the primary file with its own sidecar.
    #[cfg(unix)]
    #[rstest]
    #[case::bam("out.bam")]
    #[case::cram("out.cram")]
    fn a_sidecar_aliasing_the_resolved_destination_is_rejected(#[case] name: &str) {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join(name);
        let index = destination_index(&out);
        std::fs::write(&index, b"previous").unwrap();
        std::os::unix::fs::symlink(index.file_name().unwrap(), &out).unwrap();
        let staged = StagedSortOutput::create(&out, true, "sort").unwrap();
        std::fs::write(staged.path(), b"alignment").unwrap();
        std::fs::write(staged.staged_index.0.as_ref().unwrap(), b"index").unwrap();
        let stats = fgumi_sort::SortStats { total_records: 1, output_records: 1, runs_written: 0 };
        let err = staged.finish(&stats).unwrap_err();
        assert!(err.to_string().contains("same file"), "{err:#}");
        assert_eq!(std::fs::read(&index).unwrap(), b"previous");
        assert!(out.is_symlink());
    }

    /// A symlink destination remains supported when its sidecar is distinct.
    #[cfg(unix)]
    #[rstest]
    #[case::bam("out.bam")]
    #[case::cram("out.cram")]
    fn a_distinct_symlink_destination_is_published(#[case] name: &str) {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join(name);
        let target = dir.path().join("alignment");
        std::fs::write(&target, b"previous").unwrap();
        std::os::unix::fs::symlink("alignment", &out).unwrap();
        let staged = StagedSortOutput::create(&out, true, "sort").unwrap();
        std::fs::write(staged.path(), b"alignment").unwrap();
        std::fs::write(staged.staged_index.0.as_ref().unwrap(), b"index").unwrap();
        let stats = fgumi_sort::SortStats { total_records: 1, output_records: 1, runs_written: 0 };
        staged.finish(&stats).unwrap();
        assert_eq!(std::fs::read(&target).unwrap(), b"alignment");
        assert_eq!(std::fs::read(destination_index(&out)).unwrap(), b"index");
        assert!(out.is_symlink());
    }

    /// Library sorting shares the command's final-sidecar input protection.
    #[rstest]
    #[case::bam("out.bam")]
    #[cfg_attr(feature = "cram", case::cram("out.cram"))]
    fn library_sort_rejects_input_at_the_final_index(#[case] name: &str) {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().join(name);
        let input = destination_index(&output);
        fgumi_bam_io::create_raw_bam_writer(&input, &noodles::sam::Header::default(), 1, 1)
            .unwrap()
            .finish()
            .unwrap();
        let before = std::fs::read(&input).unwrap();
        let config = crate::sort::SortConfig {
            input: input.clone(),
            output: output.clone(),
            order: crate::commands::sort::SortOrderArg::Coordinate,
            threads: 1,
            max_memory: Some(64 << 20),
            tmp_dirs: vec![],
            compression: 1,
            write_index: true,
            key_types: fgumi_sort::KeyTypesSpec::Auto,
        };
        let err = crate::sort::sort_bam(&config).unwrap_err();
        assert!(err.to_string().contains("same file"), "{err:#}");
        assert_eq!(std::fs::read(input).unwrap(), before);
        assert!(!output.exists());
    }

    /// An index for stdout is refused, naming the index the name selects.
    #[test]
    fn an_index_for_stdout_is_refused() {
        let err = StagedSortOutput::create(Path::new("-"), true, "sort").err().unwrap();
        assert_eq!(err.to_string(), "a BAI index needs a file output, not stdout");
    }

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
