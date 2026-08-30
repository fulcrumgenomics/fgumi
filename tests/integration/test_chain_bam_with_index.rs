//! Chain-direct integration test for the inline BAI indexer wired onto
//! `SinkSpec::BamWithIndex` (Task 6).
//!
//! Builds a `ChainSpec` directly — no command-level driver (`Sort::execute`)
//! involved — with a `[Stage::Sort]` terminal chain and `SinkSpec::BamWithIndex`,
//! runs it via `build_for(spec)?.run()` over a programmatically generated
//! unsorted BAM, and asserts the `.bai` sidecar the *inline* indexer writes
//! (as part of the arena sink's drained-finish, not a post-pipeline re-read
//! hook) exists and parses. `IndexBamFinalizeHook` no longer exists, so a
//! regression that resurrected the re-read hook would show up as a compile
//! failure elsewhere, not here — this test's job is to prove the sink itself
//! produces a working index without it.

use std::process::Command;

use fgumi_lib::commands::common::{
    CompressionOptions, MaxTempFiles, MemoryLimit, MemoryReserve, QueueMemoryOptions,
    SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::sort::{SortOptions, SortOrderArg};
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use fgumi_raw_bam::SamBuilder;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};

/// `n` mapped, single-end records on one reference, placed at *descending*
/// positions so the input is deliberately unsorted — the coordinate sort
/// stage has real work to do before the sink writes (and indexes) its output.
fn unsorted_records(n: usize) -> Vec<fgumi_raw_bam::RawRecord> {
    (0..n)
        .map(|i| {
            let pos = i32::try_from((n - i) * 100).expect("pos fits i32");
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            b.build()
        })
        .collect()
}

fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

#[test]
fn bam_with_index_produces_inline_sidecar_not_reread() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    let header = create_minimal_header("chr1", 1_000_000);
    let records = unsorted_records(500);
    write_bam(&input_bam, &header, &records);

    let sort_options = SortOptions {
        order: SortOrderArg::Coordinate,
        key_types: None,
        max_memory: MemoryLimit::Fixed(64 * 1024 * 1024),
        memory_reserve: MemoryReserve::Auto,
        memory_per_thread: true,
        tmp_dirs: Vec::new(),
        sort_threads: None,
        merge_threads: None,
        temp_compression: 1,
        temp_codec: fgumi_sort::SpillCodec::default(),
        max_temp_files: MaxTempFiles::Auto,
        block_batch: 4,
        file_granularity: false,
    };

    let spec = ChainSpec {
        stages: vec![Stage::Sort],
        source: SourceSpec::Bam(input_bam.clone()),
        sink: SinkSpec::BamWithIndex(output_bam.clone()),
        stage_opts: StageOptionsBag { sort: Some(sort_options), ..Default::default() },
        threading: ThreadingOptions { threads: None },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        verify_crc: true,
        command_line: "fgumi sort --write-index".to_string(),
    };

    build_for(spec)
        .expect("build_for should accept a Sort-terminal BamWithIndex chain")
        .run()
        .expect(
            "running the chain should write the coordinate-sorted BAM and its inline .bai sidecar",
        );

    assert!(output_bam.exists(), "sorted output BAM was not written");

    let bai_path = fgumi_bam_io::bai_sidecar_path(&output_bam);
    assert!(bai_path.exists(), "inline indexer did not write a .bai sidecar at {bai_path:?}");

    // The sidecar must parse as a well-formed BAI: one reference (matching the
    // header), covering the records the sort just wrote.
    let index = noodles::bam::bai::fs::read(&bai_path)
        .expect("inline .bai must parse as a valid BAI index");
    assert_eq!(
        index.reference_sequences().len(),
        1,
        "expected one reference sequence's worth of bins in the inline .bai"
    );

    // Full samtools-equivalence (idxstats / region-query parity) is Task 7's
    // job; here a structural parse plus a basic region query is enough to
    // confirm the inline path produces a *usable* index, not just bytes that
    // happen to exist.
    if samtools_available() {
        let status = Command::new("samtools")
            .args(["quickcheck", output_bam.to_str().unwrap()])
            .status()
            .expect("run samtools quickcheck");
        assert!(status.success(), "samtools quickcheck failed on inline-indexed output BAM");

        let out = Command::new("samtools")
            .args(["view", "-c", output_bam.to_str().unwrap(), "chr1:1-50000"])
            .output()
            .expect("run samtools view region query");
        assert!(
            out.status.success(),
            "samtools view region query against the inline .bai failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let count: usize = String::from_utf8_lossy(&out.stdout).trim().parse().unwrap_or(0);
        assert!(count > 0, "expected at least one record in chr1:1-50000 via the inline .bai");
    } else {
        eprintln!("skipping samtools region-query check: samtools not available in PATH");
    }
}
