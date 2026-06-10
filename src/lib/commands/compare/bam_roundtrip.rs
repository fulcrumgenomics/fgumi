//! `fgumi compare bam-roundtrip` — Phase 3 bench-prep gate.
//!
//! Runs the new typed-step DSL pipeline (`ReadBgzfBlocks → BgzfDecompress
//! → FindBamBoundaries → ParseBamRecords → GroupBam → SerializeBamRecords
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
use std::path::PathBuf;

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

    /// Worker threads (must be ≥ 2; chain has 2 `Exclusive` steps).
    #[arg(short = 't', long = "threads", default_value = "4")]
    pub threads: usize,

    /// BGZF compression level for the output (1–12).
    #[arg(long = "compression-level", default_value = "1")]
    pub compression_level: u32,
}

impl Command for CompareBamRoundtrip {
    fn execute(&self, _command_line: &str) -> Result<()> {
        if self.threads < 2 {
            return Err(anyhow!(
                "bam-roundtrip needs --threads ≥ 2 (chain has Exclusive read + write)"
            ));
        }

        // Decide where the output lands.
        let (output_path, _temp) = if let Some(p) = &self.output {
            // Guard against pointing --output at the input BAM, which would
            // overwrite the file while it is still being read.
            if paths_refer_to_same_file(&self.input, p) {
                return Err(anyhow!(
                    "--output {} would clobber the input BAM {} (choose a different output path)",
                    p.display(),
                    self.input.display()
                ));
            }
            (p.clone(), None)
        } else {
            let f = tempfile::NamedTempFile::new().context("tempfile")?;
            let p = f.path().to_path_buf();
            (p, Some(f))
        };

        let cfg = RoundtripConfig::auto_tuned(self.threads)
            .with_compression_level(self.compression_level);

        log::info!(
            "bam-roundtrip: input={} output={} threads={}",
            self.input.display(),
            output_path.display(),
            self.threads
        );
        run_bam_roundtrip(&self.input, &output_path, cfg)
            .with_context(|| format!("run_bam_roundtrip on {}", self.input.display()))?;

        // Record-count parity check.
        let input_n = count_records(&self.input)?;
        let output_n = count_records(&output_path)?;
        if input_n != output_n {
            return Err(anyhow!(
                "bam-roundtrip record count mismatch: input={input_n}, output={output_n}"
            ));
        }
        log::info!("bam-roundtrip: PASS — {input_n} records round-tripped");
        Ok(())
    }
}

/// Returns `true` if `a` and `b` resolve to the same on-disk file.
///
/// Canonicalizes both paths so that distinct spellings of the same file
/// (relative vs. absolute, symlinks, `./`) are caught. Falls back to a literal
/// path comparison when canonicalization fails (e.g. the output does not exist
/// yet and its parent cannot be resolved).
fn paths_refer_to_same_file(a: &std::path::Path, b: &std::path::Path) -> bool {
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
