//! End-to-end no-op BAM round-trip helper. Wires the full block-level
//! BAM chain (`ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries →
//! DecodeRecords → GroupBam → SerializeBamRecords → BgzfCompress →
//! WriteBgzfFile`) and runs it against an input BAM, producing an output
//! BAM. Convenience for the `fgumi compare bam-roundtrip <bam>`
//! gate (Phase 3 bench-prep) and for any caller that wants to drive the
//! framework end-to-end on a real BAM.

use std::io;
use std::path::Path;

use fgumi_bam_io::PipelineReaderOpts;

use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::steps::bgzf::compress::BgzfCompress;
use crate::pipeline::steps::bgzf::decompress::BgzfDecompress;
use crate::pipeline::steps::boundaries::bam::FindBamBoundaries;
use crate::pipeline::steps::group::bam::GroupBam;
use crate::pipeline::steps::parse::decode::DecodeRecords;
use crate::pipeline::steps::serialize::SerializeBamRecords;
use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;
use crate::pipeline::steps::source::read_bam::read_bam;
use crate::pipeline::steps::tuning::BamPipelineTuning;
use fgumi_bam_io::GroupKeyConfig;

/// Configuration for [`run_bam_roundtrip`]. Backed by
/// [`BamPipelineTuning`] which mirrors legacy `auto_tuned` defaults
/// (`blocks_per_batch` 16-64 by thread count, `template_batch_size=500`,
/// etc.).
///
/// This is intentionally a thin newtype over [`BamPipelineTuning`] rather than a
/// bare type alias: it is the seam where round-trip-only knobs would attach
/// (e.g. a future `--verify` mode) without widening the shared tuning struct.
/// Until such a knob exists its two methods delegate straight through.
#[derive(Debug, Clone, Copy)]
pub struct RoundtripConfig {
    pub tuning: BamPipelineTuning,
}

impl RoundtripConfig {
    /// Auto-tune for `threads` workers. Equivalent to
    /// `RoundtripConfig { tuning: BamPipelineTuning::auto_tuned(threads) }`.
    #[must_use]
    pub fn auto_tuned(threads: usize) -> Self {
        Self { tuning: BamPipelineTuning::auto_tuned(threads) }
    }

    /// Override the BGZF compression level on the underlying tuning.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.tuning = self.tuning.with_compression_level(level);
        self
    }
}

impl Default for RoundtripConfig {
    fn default() -> Self {
        Self::auto_tuned(4)
    }
}

/// Run the full block-level BAM chain on `input` and write the output to
/// `output`. The output holds the same *set* of records as the input (same
/// count, same per-record bytes) modulo header `@PG` (Phase 3 has no `@PG`
/// mutation).
///
/// Record order is **not** byte-for-byte preserved for multi-record templates:
/// the chain routes through `GroupBam`, whose `Template::from_records`
/// canonicalizes within-template order (R1 before R2; supplementary/secondary
/// reads regrouped to match fgbio). Templates already in canonical order round
/// trip unchanged, but e.g. an input `[R1-primary, R1-supplementary, R2-primary]`
/// comes back as `[R1-primary, R2-primary, R1-supplementary]`. The CLI gate over
/// this function checks record *count* only, so it does not observe that
/// reordering.
///
/// # Errors
///
/// Returns I/O errors from file open / read / write, BAM parsing errors,
/// or `BuildError` propagated as `io::Error::other` if the pipeline graph
/// is malformed. Returns an error (rather than panicking) if the configured
/// compression level is outside `0..=12` or the thread count is outside
/// `1..=1024`.
///
/// **Output is not written atomically.** `WriteBgzfFile::new` creates `output`
/// and writes the BGZF header before the chain runs, so an error returned from
/// `build`/`run` (a pipeline or write failure) can leave a truncated or
/// partial file at `output`. The `compare bam-roundtrip` CLI stages its output
/// through a `NamedTempFile` and renames on success, so it never exposes a
/// partial file; direct library callers that need that guarantee must stage and
/// rename themselves.
pub fn run_bam_roundtrip(input: &Path, output: &Path, cfg: RoundtripConfig) -> io::Result<()> {
    let t = cfg.tuning;
    // Validate the compression level and thread count up front — before the
    // output file is created — so an out-of-range value returns an error
    // instead of panicking deep inside `InlineBgzfCompressor::new` (compression)
    // or the worker-thread spawn loop in `Pipeline::run` (threads). The CLI
    // guards both with clap `value_parser` ranges, but this function is public
    // API and a direct caller must get an error, not a panic. `auto_tuned`
    // floors `threads` to 1, so this catches a `RoundtripConfig` whose `tuning`
    // was built with a raw 0 or an unbounded count.
    if t.compression_level > 12 {
        return Err(io::Error::other(format!(
            "compression level {} is out of range (expected 0..=12)",
            t.compression_level
        )));
    }
    if !(1..=1024).contains(&t.threads) {
        return Err(io::Error::other(format!(
            "thread count {} is out of range (expected 1..=1024)",
            t.threads
        )));
    }
    let opts = PipelineReaderOpts::default();
    let (read_step, header) = read_bam(input, opts, t.blocks_per_batch, t.per_step_byte_limit)?;

    let builder = Pipeline::builder();
    builder
        .chain(read_step)
        .chain(BgzfDecompress::new(t.per_step_byte_limit))
        .chain(FindBamBoundaries::new(t.per_step_byte_limit))
        // DecodeRecords does parse + per-record GroupKey computation
        // (CIGAR walk + tag scan + library/cell-barcode hashing) in one
        // parallel pass. Matches the legacy pipeline's combined Decode
        // step; drops the intermediate Parse → Decode queue.
        // `GroupKeyConfig::default()` (empty library index, no cell tag) is used
        // deliberately: this round-trip validates byte-equivalence only and does
        // not exercise UMI/library/cell grouping (`GroupBam` batches by qname
        // regardless), so the computed group key does not affect output bytes.
        .chain(DecodeRecords::new(GroupKeyConfig::default(), t.per_step_byte_limit))
        .chain(GroupBam::new(t.template_batch_size, t.per_step_byte_limit))
        .chain(SerializeBamRecords::new(t.per_step_byte_limit))
        .chain(BgzfCompress::new(t.compression_level, t.per_step_byte_limit, false))
        .chain(
            WriteBgzfFile::new(output, &header, t.compression_level)
                .map_err(|e| io::Error::other(format!("WriteBgzfFile::new: {e}")))?,
        )
        .into_sink_marker();
    let pipeline =
        builder.build().map_err(|e| io::Error::other(format!("Pipeline::build: {e:?}")))?;
    pipeline
        .run(PipelineConfig { threads: t.threads, ..Default::default() })
        .map_err(|e| io::Error::other(format!("Pipeline::run: {e:?}")))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use fgumi_bam_io::PipelineReaderOpts;
    use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
    use noodles::sam::Header;
    use rstest::rstest;

    /// Write `records` to a temp BAM, in order, and return the temp file (kept
    /// alive by the caller so the on-disk file persists for the round-trip).
    fn write_bam(records: &[RawRecord]) -> tempfile::NamedTempFile {
        use noodles::sam::alignment::io::Write as _;
        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let header = Header::default();
        let mut writer =
            noodles::bam::io::Writer::new(std::fs::File::create(tmp.path()).expect("create BAM"));
        writer.write_header(&header).expect("write header");
        for record in records {
            let buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
                .expect("raw_record_to_record_buf");
            writer.write_alignment_record(&header, &buf).expect("write record");
        }
        writer.try_finish().expect("finish BAM");
        tmp
    }

    /// Read every record body from a BAM, in file order, as raw bytes. Two BAMs
    /// hold the same records in the same order iff these vectors are equal:
    /// `RawRecord` derefs to its record-body `[u8]`, so byte equality is record
    /// identity (name, flags, SEQ, QUAL, CIGAR, tags — everything but the
    /// header, which the round-trip legitimately re-serializes).
    fn record_bodies(path: &Path) -> Vec<Vec<u8>> {
        let (mut reader, _hdr) =
            fgumi_bam_io::create_raw_bam_reader_with_opts(path, 1, PipelineReaderOpts::default())
                .expect("open BAM");
        let mut record = RawRecord::default();
        let mut bodies = Vec::new();
        while reader.read_record(&mut record).expect("read_record") != 0 {
            bodies.push(record.to_vec());
        }
        bodies
    }

    /// A heterogeneous input that exercises the paths a count-only check misses:
    /// enough singleton records to span multiple BGZF blocks and template
    /// batches (`template_batch_size` defaults to 500); a multi-record template
    /// (paired R1/R2 sharing one QNAME, so a drop/duplicate *within* a template
    /// is visible); and degenerate records (unmapped, empty SEQ) that a
    /// fixed-shape input would never cover.
    fn varied_records() -> Vec<RawRecord> {
        let mut records = Vec::new();
        // Bulk singletons across several blocks/batches.
        for i in 0..1500 {
            records.push(
                SamBuilder::new()
                    .read_name(format!("read{i:07}").as_bytes())
                    .sequence(&[b'A'; 60])
                    .qualities(&[30; 60])
                    .flags(flags::FIRST_SEGMENT)
                    .build(),
            );
        }
        // Multi-record template: two mates under one QNAME, emitted adjacently.
        records.push(
            SamBuilder::new()
                .read_name(b"pair_template")
                .sequence(&[b'C'; 40])
                .qualities(&[35; 40])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT)
                .build(),
        );
        records.push(
            SamBuilder::new()
                .read_name(b"pair_template")
                .sequence(&[b'G'; 40])
                .qualities(&[35; 40])
                .flags(flags::PAIRED | flags::LAST_SEGMENT)
                .build(),
        );
        // Degenerate: unmapped, empty SEQ/QUAL.
        records.push(SamBuilder::new().read_name(b"empty_seq").flags(flags::UNMAPPED).build());
        records
    }

    /// The core contract of the no-op round-trip: the output holds exactly the
    /// input records, byte-for-byte, in the same order — not merely the same
    /// count. Runs at 1, 2, and 4 threads: threads=1 takes the fused
    /// single-thread path, threads>=2 the scheduled pool where the reader/writer
    /// affinity split takes effect, so all three must agree.
    #[rstest]
    fn roundtrip_preserves_records_identically(#[values(1, 2, 4)] threads: usize) {
        let input = write_bam(&varied_records());
        let output = tempfile::NamedTempFile::new().expect("create output BAM");

        run_bam_roundtrip(input.path(), output.path(), RoundtripConfig::auto_tuned(threads))
            .unwrap_or_else(|e| panic!("run_bam_roundtrip at threads={threads}: {e}"));

        assert_eq!(
            record_bodies(output.path()),
            record_bodies(input.path()),
            "round-trip at threads={threads} did not reproduce the input records identically"
        );
    }

    /// An empty BAM (header only, zero records) must round-trip to an empty BAM,
    /// not error and not fabricate records — the drain path has to flush a chain
    /// that never sees an input record.
    #[test]
    fn roundtrip_handles_empty_input() {
        let input = write_bam(&[]);
        let output = tempfile::NamedTempFile::new().expect("create output BAM");

        run_bam_roundtrip(input.path(), output.path(), RoundtripConfig::auto_tuned(4))
            .expect("run_bam_roundtrip on empty input");

        assert!(record_bodies(output.path()).is_empty(), "empty input must yield empty output");
    }

    /// A compression level above 12 returns an error from the library API
    /// instead of panicking deep in `InlineBgzfCompressor::new`, and does so
    /// without creating the output file.
    #[test]
    fn out_of_range_compression_level_errors_without_touching_output() {
        let input = write_bam(&varied_records());
        let out_dir = tempfile::tempdir().expect("temp dir");
        let output = out_dir.path().join("should_not_exist.bam");

        let cfg = RoundtripConfig::auto_tuned(4).with_compression_level(13);
        let err = run_bam_roundtrip(input.path(), &output, cfg)
            .expect_err("compression level 13 must error, not panic");
        assert!(err.to_string().contains("compression level"), "unexpected error: {err}");
        assert!(!output.exists(), "no output file should be created on a config error");
    }

    /// A thread count outside `1..=1024` returns an error from the library API
    /// instead of attempting an unbounded worker-thread spawn (which panics on
    /// the first failed `thread::Builder::spawn`), and does so without creating
    /// the output file. `auto_tuned` floors to 1 but does not cap, so a large
    /// count reaches this guard; the CLI caps it at parse time, but a direct
    /// library caller must get an error, not a panic.
    #[test]
    fn out_of_range_thread_count_errors_without_touching_output() {
        let input = write_bam(&varied_records());
        let out_dir = tempfile::tempdir().expect("temp dir");
        let output = out_dir.path().join("should_not_exist.bam");

        let cfg = RoundtripConfig::auto_tuned(2000);
        let err = run_bam_roundtrip(input.path(), &output, cfg)
            .expect_err("thread count 2000 must error, not spawn 2000 threads");
        assert!(err.to_string().contains("thread count"), "unexpected error: {err}");
        assert!(!output.exists(), "no output file should be created on a config error");
    }

    /// A multi-record template whose input order is non-canonical round-trips
    /// with the same record *set* but reordered: the chain routes through
    /// `GroupBam`, which canonicalizes within-template order (R1 before R2).
    /// This pins the behavior `run_bam_roundtrip`'s doc calls out — the
    /// round-trip preserves records, not per-template order — which the
    /// canonical-order inputs in `varied_records` cannot exercise.
    #[test]
    fn non_canonical_template_order_is_canonicalized() {
        // Emit the mates R2-then-R1 (non-canonical); GroupBam will swap them.
        let r2 = SamBuilder::new()
            .read_name(b"pair")
            .sequence(&[b'G'; 40])
            .qualities(&[35; 40])
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .build();
        let r1 = SamBuilder::new()
            .read_name(b"pair")
            .sequence(&[b'C'; 40])
            .qualities(&[35; 40])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .build();
        let input = write_bam(&[r2, r1]);
        let output = tempfile::NamedTempFile::new().expect("create output BAM");

        run_bam_roundtrip(input.path(), output.path(), RoundtripConfig::auto_tuned(1))
            .expect("run_bam_roundtrip");

        let in_bodies = record_bodies(input.path());
        let out_bodies = record_bodies(output.path());

        // The record set is preserved (same count, same per-record bytes)...
        let (mut in_sorted, mut out_sorted) = (in_bodies.clone(), out_bodies.clone());
        in_sorted.sort();
        out_sorted.sort();
        assert_eq!(out_sorted, in_sorted, "the record set must be preserved");
        // ...but the non-canonical [R2, R1] input is emitted canonicalized, so
        // the file order changes — which the count-only CLI gate cannot see.
        assert_ne!(out_bodies, in_bodies, "non-canonical template order must be canonicalized");
    }
}
