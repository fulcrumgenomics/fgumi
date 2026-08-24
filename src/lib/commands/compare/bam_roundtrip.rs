//! `fgumi compare bam-roundtrip` — Phase 3 bench-prep gate.
//!
//! Runs the new typed-step DSL pipeline (`ReadBgzfBlocks → BgzfDecompress
//! → FindBamBoundaries → DecodeRecords → GroupBam → SerializeBamRecords
//! → BgzfCompress → WriteBgzfFile`) on an input BAM and writes the output
//! to a temp file, then compares record-count parity to the input.
//!
//! Phase 3's bench gate is **partial**: it verifies the chain processes
//! and re-serializes records correctly, but does not assert byte-for-byte
//! equivalence of header bytes (which the `WriteBgzfFile` rewrites with a
//! fresh BGZF compressor). Phase 4's `--preset correct` is the original
//! design's full equivalence gate.

use anyhow::{Context, Result, anyhow};
use clap::Parser;
use std::path::{Path, PathBuf};

use crate::commands::command::Command;
use crate::pipeline::steps::{RoundtripConfig, run_bam_roundtrip};

/// Run the new pipeline as a no-op chain on an input BAM and assert
/// record-count parity.
#[derive(Debug, Parser)]
#[command(
    name = "bam-roundtrip",
    about = "Run the new pipeline on a BAM and verify record-count parity"
)]
pub struct CompareBamRoundtrip {
    /// Input BAM file.
    #[arg(index = 1)]
    pub input: PathBuf,

    /// Optional output BAM path. If unset, output goes to a tempfile and
    /// is discarded after the comparison.
    #[arg(long = "output")]
    pub output: Option<PathBuf>,

    /// Worker threads (1–1024). The chain runs its BAM reader and writer on
    /// separate worker threads when 2+ are given; a single thread falls back to
    /// the fused single-thread path. The upper bound is a guard against a
    /// pathological value allocating that many worker slots and OS threads —
    /// `PipelineConfig` propagates the count unclamped, so it is bounded here.
    #[arg(
        short = 't',
        long = "threads",
        default_value = "4",
        value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..=1024)
    )]
    pub threads: usize,

    /// BGZF compression level for the output (0–12).
    #[arg(
        long = "compression-level",
        default_value = "1",
        value_parser = clap::value_parser!(u32).range(0..=12)
    )]
    pub compression_level: u32,
}

impl Command for CompareBamRoundtrip {
    fn execute(&self, _command_line: &str) -> Result<()> {
        // Fail early with a clear message if the input is missing, matching the
        // sibling compare subcommands rather than surfacing a raw open error.
        crate::validation::validate_file_exists(&self.input, "Input BAM")?;

        // Guard against pointing --output at the input BAM, which would
        // overwrite the file while it is still being read.
        if let Some(p) = &self.output
            && paths_refer_to_same_file(&self.input, p)
        {
            return Err(anyhow!(
                "--output {} would clobber the input BAM {} (choose a different output path)",
                p.display(),
                self.input.display()
            ));
        }

        // Always write to a temporary staging file, then move it into place only
        // after the parity check passes. This keeps --output atomic: a failed run
        // (or a parity mismatch) never leaves a truncated or unverified BAM at the
        // requested path — the staging file is discarded on any early return. When
        // --output is given, the staging file is created in that file's own
        // directory so the final persist is a same-filesystem rename, not a
        // cross-device copy.
        let staging = match &self.output {
            Some(p) => {
                let dir = p
                    .parent()
                    .filter(|d| !d.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                tempfile::NamedTempFile::new_in(dir).context("create staging tempfile")?
            }
            None => tempfile::NamedTempFile::new().context("create tempfile")?,
        };
        let staging_path = staging.path().to_path_buf();

        let cfg = RoundtripConfig::auto_tuned(self.threads)
            .with_compression_level(self.compression_level);

        log::info!(
            "bam-roundtrip: input={} output={} threads={}",
            self.input.display(),
            self.output.as_deref().unwrap_or(&staging_path).display(),
            self.threads
        );
        run_bam_roundtrip(&self.input, &staging_path, cfg)
            .with_context(|| format!("run_bam_roundtrip on {}", self.input.display()))?;

        // Record-count parity check against the staging file.
        let input_n = count_records(&self.input)?;
        let output_n = count_records(&staging_path)?;
        if input_n != output_n {
            return Err(anyhow!(
                "bam-roundtrip record count mismatch: input={input_n}, output={output_n}"
            ));
        }

        // Validation passed: move the staging file into place (if requested).
        if let Some(p) = &self.output {
            staging.persist(p).map_err(|e| {
                anyhow!("failed to move validated output into {}: {}", p.display(), e.error)
            })?;
        }
        log::info!("bam-roundtrip: PASS — {input_n} records round-tripped");
        Ok(())
    }
}

/// Returns `true` if `a` and `b` resolve to the same on-disk file.
///
/// When both files already exist, compares filesystem identity (device + inode
/// on Unix) so that *hard links* to the same file are caught as well as symlinks
/// — a hard link shares the input's inode but has a distinct canonical path that
/// path comparison alone would miss, letting `--output` overwrite the input BAM
/// mid-run. Falls back to canonical-path comparison, then to a literal path
/// comparison, when either file does not exist yet (e.g. a fresh `--output`
/// whose parent cannot be resolved) so an exact-string match is still rejected.
fn paths_refer_to_same_file(a: &std::path::Path, b: &std::path::Path) -> bool {
    // Prefer filesystem identity when both paths exist: `std::fs::metadata`
    // follows symlinks, so this also collapses symlink aliases onto their target,
    // and unlike a canonical-path compare it detects hard links (same dev+inode,
    // different path).
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(ma), Ok(mb)) = (std::fs::metadata(a), std::fs::metadata(b)) {
            return ma.dev() == mb.dev() && ma.ino() == mb.ino();
        }
    }
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(ca), Ok(cb)) => ca == cb,
        // If either side can't be canonicalized (e.g. `b` doesn't exist yet),
        // fall back to a direct comparison so an exact-string match is still
        // rejected.
        _ => a == b,
    }
}

fn count_records(path: &std::path::Path) -> Result<u64> {
    use fgumi_bam_io::PipelineReaderOpts;
    use fgumi_raw_bam::RawRecord;

    let (mut reader, _hdr) =
        fgumi_bam_io::create_raw_bam_reader_with_opts(path, 1, PipelineReaderOpts::default())
            .with_context(|| format!("open BAM {}", path.display()))?;
    let mut record = RawRecord::default();
    let mut n: u64 = 0;
    loop {
        match reader.read_record(&mut record).context("read_record")? {
            0 => break,
            _ => n += 1,
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
    use noodles::sam::Header;

    /// Write a small valid BAM of `n` records to a temp file for the command to
    /// consume.
    fn write_input_bam(n: usize) -> tempfile::NamedTempFile {
        use noodles::sam::alignment::io::Write as _;
        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let header = Header::default();
        let mut writer =
            noodles::bam::io::Writer::new(std::fs::File::create(tmp.path()).expect("create BAM"));
        writer.write_header(&header).expect("write header");
        for i in 0..n {
            let record: RawRecord = SamBuilder::new()
                .read_name(format!("read{i:07}").as_bytes())
                .sequence(&[b'A'; 60])
                .qualities(&[30; 60])
                .flags(flags::FIRST_SEGMENT)
                .build();
            let buf = fgumi_raw_bam::raw_record_to_record_buf(&record, &header)
                .expect("raw_record_to_record_buf");
            writer.write_alignment_record(&header, &buf).expect("write record");
        }
        writer.try_finish().expect("finish BAM");
        tmp
    }

    /// Happy path: a real BAM round-trips to an explicit `--output` and the
    /// parity check inside `execute` passes.
    #[test]
    fn execute_roundtrips_to_explicit_output() {
        let input = write_input_bam(1500);
        let output = tempfile::NamedTempFile::new().expect("create output path");
        let cmd = CompareBamRoundtrip {
            input: input.path().to_path_buf(),
            output: Some(output.path().to_path_buf()),
            threads: 4,
            compression_level: 1,
        };
        cmd.execute("compare bam-roundtrip").expect("round-trip should pass parity");
        assert_eq!(count_records(output.path()).expect("count output"), 1500);
    }

    /// `--output` is written atomically: a run that fails before the parity check
    /// passes must not leave a partial or unverified BAM at the requested path.
    /// Here the output's directory does not exist, so staging fails early — and
    /// the requested `--output` path must not exist afterward.
    #[test]
    fn execute_does_not_create_output_on_early_failure() {
        let input = write_input_bam(10);
        let dir = tempfile::tempdir().expect("temp dir");
        let missing_output = dir.path().join("no_such_subdir").join("out.bam");
        let cmd = CompareBamRoundtrip {
            input: input.path().to_path_buf(),
            output: Some(missing_output.clone()),
            threads: 4,
            compression_level: 1,
        };
        cmd.execute("compare bam-roundtrip").expect_err("staging into a missing dir must error");
        assert!(!missing_output.exists(), "no output must be left behind on failure");
    }

    /// `--threads` is range-checked at parse time (1..=1024) so a pathological
    /// value is rejected by clap before it reaches `PipelineConfig`, which
    /// propagates the worker count unclamped.
    #[test]
    fn rejects_out_of_range_threads() {
        for ok in ["1", "4", "1024"] {
            assert!(
                CompareBamRoundtrip::try_parse_from(["bam-roundtrip", "in.bam", "-t", ok]).is_ok(),
                "threads {ok} should parse"
            );
        }
        for bad in ["0", "1025"] {
            assert!(
                CompareBamRoundtrip::try_parse_from(["bam-roundtrip", "in.bam", "-t", bad])
                    .is_err(),
                "threads {bad} must be rejected at parse time"
            );
        }
    }

    /// A single worker is a valid configuration: the linear read->...->write
    /// chain is fusible, so `threads=1` runs on the fused single-thread path and
    /// still passes the parity check.
    #[test]
    fn execute_succeeds_at_one_thread() {
        let input = write_input_bam(500);
        let output = tempfile::NamedTempFile::new().expect("create output path");
        let cmd = CompareBamRoundtrip {
            input: input.path().to_path_buf(),
            output: Some(output.path().to_path_buf()),
            threads: 1,
            compression_level: 1,
        };
        cmd.execute("compare bam-roundtrip").expect("threads=1 should round-trip");
        assert_eq!(count_records(output.path()).expect("count output"), 500);
    }

    /// A missing input file is rejected up front with a clear message, matching
    /// the sibling compare subcommands rather than surfacing a raw open error.
    #[test]
    fn execute_rejects_missing_input() {
        let dir = tempfile::tempdir().expect("temp dir");
        let cmd = CompareBamRoundtrip {
            input: dir.path().join("does_not_exist.bam"),
            output: None,
            threads: 4,
            compression_level: 1,
        };
        cmd.execute("compare bam-roundtrip").expect_err("missing input must error");
    }

    /// Pointing `--output` at the input BAM must be rejected: the chain reads the
    /// input while writing the output, so clobbering it would corrupt the run.
    #[test]
    fn execute_rejects_output_clobbering_input() {
        let input = write_input_bam(1);
        let cmd = CompareBamRoundtrip {
            input: input.path().to_path_buf(),
            output: Some(input.path().to_path_buf()),
            threads: 4,
            compression_level: 1,
        };
        let err = cmd.execute("compare bam-roundtrip").expect_err("clobber must error");
        assert!(err.to_string().contains("clobber"), "unexpected error: {err}");
    }

    /// `--compression-level` is range-checked at parse time so an out-of-range
    /// value is rejected by clap before any work runs — the value otherwise
    /// reaches `InlineBgzfCompressor::new`, which panics for anything above 12.
    #[test]
    fn rejects_out_of_range_compression_level() {
        // In range (0..=12) parses.
        for level in ["0", "1", "12"] {
            assert!(
                CompareBamRoundtrip::try_parse_from([
                    "bam-roundtrip",
                    "in.bam",
                    "--compression-level",
                    level,
                ])
                .is_ok(),
                "level {level} should parse"
            );
        }
        // Out of range is rejected by the value_parser, not by a later panic.
        assert!(
            CompareBamRoundtrip::try_parse_from([
                "bam-roundtrip",
                "in.bam",
                "--compression-level",
                "13",
            ])
            .is_err(),
            "compression level 13 must be rejected at parse time"
        );
    }

    /// The same-file guard canonicalizes, so a genuinely distinct path spelling
    /// that resolves to one file is caught, while two distinct files are not.
    /// The symlink case is what a naive literal `a == b` implementation would
    /// miss — `Path` equality normalizes away `.` but not a symlink indirection.
    #[test]
    fn same_file_guard_matches_distinct_spellings() {
        let f = tempfile::NamedTempFile::new().expect("temp file");
        let abs = f.path().to_path_buf();
        assert!(paths_refer_to_same_file(&abs, &abs), "identical paths must match");

        let other = tempfile::NamedTempFile::new().expect("temp file 2");
        assert!(!paths_refer_to_same_file(&abs, other.path()), "distinct files must not match");

        // A symlink is a distinct path spelling that only canonicalization
        // collapses onto its target; a literal path comparison would not.
        #[cfg(unix)]
        {
            let link = abs.with_extension("link");
            std::os::unix::fs::symlink(&abs, &link).expect("create symlink");
            assert_ne!(link, abs, "the symlink path differs literally from its target");
            assert!(
                paths_refer_to_same_file(&abs, &link),
                "a symlink to the input must be caught by canonicalization"
            );
            std::fs::remove_file(&link).ok();
        }
    }

    /// A hard link shares the input's inode but has a genuinely distinct
    /// canonical path, so a path-only comparison misses it. Filesystem-identity
    /// matching (device + inode) is what catches an `--output` hard-linked to the
    /// input, which would otherwise overwrite the source file mid-run.
    #[cfg(unix)]
    #[test]
    fn same_file_guard_matches_hard_link() {
        let f = tempfile::NamedTempFile::new().expect("temp file");
        let target = f.path().to_path_buf();
        let link = target.with_extension("hardlink");
        std::fs::hard_link(&target, &link).expect("create hard link");

        // Both paths canonicalize to themselves — the hard link is not a symlink,
        // so canonicalization does not collapse the two spellings.
        assert_ne!(
            std::fs::canonicalize(&target).expect("canon target"),
            std::fs::canonicalize(&link).expect("canon link"),
            "a hard link has a distinct canonical path from its target"
        );
        assert!(
            paths_refer_to_same_file(&target, &link),
            "a hard link to the input must be caught by filesystem identity"
        );
        std::fs::remove_file(&link).ok();
    }
}
