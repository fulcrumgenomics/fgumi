//! Shared fixtures for the in-process sort-chain tests: a `ChainSpec` for a
//! BAM-to-BAM chain, default `SortOptions`, tie-heavy unsorted records and a
//! reader for the counts `fgumi sort` logs.

use std::path::Path;

use fgumi_lib::commands::common::{
    CompressionOptions, MaxTempFiles, MemoryLimit, MemoryReserve, QueueMemoryOptions,
    SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::group::GroupOptions;
use fgumi_lib::commands::sort::{SortOptions, SortOrderArg};
use fgumi_lib::pipeline::chains::{ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};

/// Coordinate-order `SortOptions` with a 64 MiB per-thread budget and the
/// given phase-thread overrides.
pub fn sort_options(sort_threads: Option<usize>, merge_threads: Option<usize>) -> SortOptions {
    SortOptions {
        order: SortOrderArg::Coordinate,
        key_types: None,
        max_memory: MemoryLimit::Fixed(64 * 1024 * 1024),
        memory_reserve: MemoryReserve::Auto,
        memory_per_thread: true,
        tmp_dirs: Vec::new(),
        sort_threads,
        merge_threads,
        temp_compression: 1,
        temp_codec: fgumi_sort::SpillCodec::default(),
        max_temp_files: MaxTempFiles::Auto,
        sort_stats: false,
    }
}

/// A chain spec running `stages` from `source` into a BAM at `output`.
pub fn spec_with(
    stages: Vec<Stage>,
    source: SourceSpec,
    output: &Path,
    threads: Option<usize>,
    stage_opts: StageOptionsBag,
) -> ChainSpec {
    ChainSpec {
        stages,
        source,
        sink: SinkSpec::Bam(output.to_path_buf()),
        stage_opts,
        threading: ThreadingOptions { threads },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: true,
        command_line: "fgumi sort <chain test>".to_string(),
    }
}

/// A BAM-to-BAM chain spec over `stages` with `sort` (and default group
/// options when `stages` has a group stage).
pub fn spec(
    stages: Vec<Stage>,
    input: &Path,
    output: &Path,
    threads: usize,
    sort: SortOptions,
) -> ChainSpec {
    let group = stages.contains(&Stage::Group).then(GroupOptions::default);
    let bag = StageOptionsBag { sort: Some(sort), group, ..Default::default() };
    spec_with(stages, SourceSpec::Bam(input.to_path_buf()), output, Some(threads), bag)
}

/// Unsorted single-end records with a unique name per record (for commands such
/// as `group` that reject two primary reads sharing a name); see [`records`].
pub fn uniquely_named_records(n: usize) -> Vec<RawRecord> {
    records(n, |i| format!("q{i}"))
}

/// `n` unsorted single-end records engineered to tie under every sort order:
/// a small pool of positions (coordinate ties), names from `name_of`, and a
/// unique `RX` per record so two tied records are still distinguishable in the
/// output bytes. Every 17th record is unmapped (no reference, no position), so
/// all of those tie on the coordinate sort's unmapped sentinel key.
pub fn records<F: Fn(usize) -> String>(n: usize, name_of: F) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let name = name_of(i);
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes()).sequence(b"ACGTACGTAC").qualities(&[30u8; 10]);
            if i % 17 == 0 {
                b.ref_id(-1).pos(-1).mapq(0).flags(flags::UNMAPPED);
            } else {
                let pos = i32::try_from(1 + (i * 7919 % 500) * 10).expect("pos fits i32");
                let strand = if i % 2 == 0 { 0 } else { flags::REVERSE };
                b.ref_id(0).pos(pos).mapq(60).flags(strand).cigar_ops(&[10u32 << 4]); // 10M
            }
            b.add_string_tag(SamTag::RX, format!("u{i}").as_bytes());
            b.build()
        })
        .collect()
}

/// The integer after `label` on the first stderr line containing it.
pub fn logged_count(stderr: &str, label: &str) -> Option<u64> {
    let line = stderr.lines().find(|l| l.contains(label))?;
    let rest = &line[line.find(label)? + label.len()..];
    rest.split_whitespace().next()?.parse().ok()
}
