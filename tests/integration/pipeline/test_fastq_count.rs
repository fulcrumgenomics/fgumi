//! End-to-end FASTQ record-counting integration tests for `pipeline`.
//!
//! Each test writes a small FASTQ file (plain, gzip, or BGZF), wires up the
//! full FASTQ input pipeline chain
//! (`FastqFileRead` -> [`BgzfDecompress` ->] [`FastqBlockParse` +
//! `FastqBlockMerge` | `FastqParse` + `FastqPair`] -> `CountStage` ->
//! `CountSink`), runs it, and asserts the observed record count matches what
//! was written.
//!
//! Paired-end tests expect the count to equal the per-file record count (not
//! the total), because each pair of records zips into a single template.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use fgumi_lib::runall::engine::cancel::CancelToken;
use fgumi_lib::runall::engine::driver::{PipelineConfig, StageEntry, run_multi_stage};
use fgumi_lib::runall::engine::sink::count::CountSink;
use fgumi_lib::runall::engine::source::fastq::{FastqFileRead, FastqFormat, detect_fastq_format};
use fgumi_lib::runall::engine::stage_source::{ErasedStageSource, StageSource};
use fgumi_lib::runall::engine::stages::{
    BgzfDecompress, CountStage, FastqBlockMerge, FastqBlockParse, FastqPair, FastqParse,
};

// ============================================================================
// Helpers: write test FASTQ files in each supported format
// ============================================================================

fn write_plain_fastq(path: &Path, records: usize) {
    let mut f = std::fs::File::create(path).unwrap();
    for i in 0..records {
        writeln!(f, "@read{i}").unwrap();
        writeln!(f, "ACGTACGTACGT").unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "IIIIIIIIIIII").unwrap();
    }
}

fn write_gzip_fastq(path: &Path, records: usize) {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    let f = std::fs::File::create(path).unwrap();
    let mut enc = GzEncoder::new(f, Compression::default());
    for i in 0..records {
        writeln!(enc, "@read{i}").unwrap();
        writeln!(enc, "ACGTACGTACGT").unwrap();
        writeln!(enc, "+").unwrap();
        writeln!(enc, "IIIIIIIIIIII").unwrap();
    }
    enc.finish().unwrap();
}

fn write_bgzf_fastq(path: &Path, records: usize) {
    use fgumi_bgzf::InlineBgzfCompressor;
    let mut compressor = InlineBgzfCompressor::new(6);
    for i in 0..records {
        let record = format!("@read{i}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n");
        compressor.write_all(record.as_bytes()).unwrap();
    }
    compressor.flush().unwrap();
    let blocks = compressor.take_blocks();
    let mut f = std::fs::File::create(path).unwrap();
    for block in blocks {
        f.write_all(&block.data).unwrap();
    }
    f.flush().unwrap();
}

// ============================================================================
// Pipeline runner
// ============================================================================

fn run_count_pipeline(paths: Vec<PathBuf>) -> u64 {
    let counter = Arc::new(AtomicU64::new(0));
    let num_streams = paths.len();
    let fmt = detect_fastq_format(&paths[0]).unwrap();
    for p in &paths[1..] {
        assert_eq!(
            detect_fastq_format(p).unwrap(),
            fmt,
            "all inputs must be same format for this test harness"
        );
    }

    let source = Box::new(FastqFileRead::new(paths).unwrap());

    let stages: Vec<StageEntry> = match fmt {
        FastqFormat::Bgzf => vec![
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(
                BgzfDecompress::new(),
            ))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(FastqBlockParse))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(
                FastqBlockMerge::new(num_streams),
            ))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(CountStage::new(
                counter.clone(),
            )))),
        ],
        FastqFormat::Gzip | FastqFormat::Plain => vec![
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(FastqParse))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(FastqPair::new(
                num_streams,
            )))),
            StageEntry::Pool(ErasedStageSource::from_source(StageSource::Clone(CountStage::new(
                counter.clone(),
            )))),
        ],
    };

    let sink = Box::new(CountSink);

    run_multi_stage(
        source,
        stages,
        sink,
        &CancelToken::new(),
        PipelineConfig {
            worker_threads: 4,
            queue_capacity: 64,
            queue_memory_limit: 4 * 1024 * 1024,
            global_memory_limit: 64 * 1024 * 1024,
        },
    )
    .unwrap();

    counter.load(Ordering::Relaxed)
}

// ============================================================================
// Tests
// ============================================================================

#[test]
fn test_count_single_plain_fastq() {
    let tmp = tempfile::tempdir().unwrap();
    let p = tmp.path().join("test.fastq");
    write_plain_fastq(&p, 1000);
    assert_eq!(run_count_pipeline(vec![p]), 1000);
}

#[test]
fn test_count_single_gzip_fastq() {
    let tmp = tempfile::tempdir().unwrap();
    let p = tmp.path().join("test.fastq.gz");
    write_gzip_fastq(&p, 1000);
    assert_eq!(run_count_pipeline(vec![p]), 1000);
}

#[test]
fn test_count_single_bgzf_fastq() {
    let tmp = tempfile::tempdir().unwrap();
    let p = tmp.path().join("test.fastq.bgz");
    write_bgzf_fastq(&p, 1000);
    assert_eq!(run_count_pipeline(vec![p]), 1000);
}

#[test]
fn test_count_paired_plain_fastq() {
    let tmp = tempfile::tempdir().unwrap();
    let p1 = tmp.path().join("R1.fastq");
    let p2 = tmp.path().join("R2.fastq");
    write_plain_fastq(&p1, 500);
    write_plain_fastq(&p2, 500);
    assert_eq!(run_count_pipeline(vec![p1, p2]), 500);
}

#[test]
fn test_count_paired_gzip_fastq() {
    let tmp = tempfile::tempdir().unwrap();
    let p1 = tmp.path().join("R1.fastq.gz");
    let p2 = tmp.path().join("R2.fastq.gz");
    write_gzip_fastq(&p1, 500);
    write_gzip_fastq(&p2, 500);
    assert_eq!(run_count_pipeline(vec![p1, p2]), 500);
}

#[test]
fn test_count_paired_bgzf_fastq() {
    let tmp = tempfile::tempdir().unwrap();
    let p1 = tmp.path().join("R1.fastq.bgz");
    let p2 = tmp.path().join("R2.fastq.bgz");
    write_bgzf_fastq(&p1, 500);
    write_bgzf_fastq(&p2, 500);
    assert_eq!(run_count_pipeline(vec![p1, p2]), 500);
}
