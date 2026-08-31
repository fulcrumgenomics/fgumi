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
//!
//! ## Task 7 additions: the arena-chain correctness gate
//!
//! `fgumi sort --write-index` now routes through this same chain builder (via
//! `SinkSpec::BamWithIndex`, constructed in `Sort::execute`) — see
//! `task-7-report.md` for the investigation that predates the cutover. This
//! test still builds the `ChainSpec` directly rather than going through the
//! CLI, so it stays a chain-level gate independent of `Sort::execute`'s own
//! argument handling. The three tests below are therefore the correctness
//! gate for the arena sink's inline `.bai`, run the same chain-direct way as
//! `bam_with_index_produces_inline_sidecar_not_reread` above, but with the
//! full samtools-equivalence and byte-identity assertions Task 7 specifies:
//!
//! - [`arena_bam_with_index_idxstats_equivalent_to_samtools`][]: region-query
//!   and `idxstats` parity against samtools' own index, over placed, placed-
//!   but-unmapped, and truly-unplaced reads across multiple references.
//! - [`arena_bam_with_index_multi_block_and_straddle`][]: the same parity
//!   check forced through many physical BGZF blocks and one record that
//!   straddles a physical-block boundary by itself.
//! - [`arena_bam_with_index_detached_writer_is_deterministic_and_correct`][]:
//!   Task 7's third requirement (Detached-vs-Serial byte-identity) is only
//!   *partially* achievable here — see that test's doc comment for why, and
//!   what is asserted instead.

use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::Command;

use fgumi_lib::commands::common::{
    CompressionOptions, MaxTempFiles, MemoryLimit, MemoryReserve, QueueMemoryOptions,
    SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::sort::{SortOptions, SortOrderArg};
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
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
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
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

// =============================================================================
// Task 7: arena-chain correctness gate — shared helpers
// =============================================================================

/// A header with `refs.len()` references, each `(name, length)`.
fn multi_ref_header(refs: &[(&str, usize)]) -> Header {
    let mut builder = Header::builder();
    for &(name, len) in refs {
        let seq = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(len).expect("reference length must be non-zero"),
        );
        builder = builder.add_reference_sequence(bstr::BString::from(name), seq);
    }
    builder.build()
}

/// A deterministic scattered (i.e. deliberately unsorted) position in
/// `[1, modulus]`, salted per caller so different record categories don't
/// collide on the same positions.
fn scattered_pos(i: usize, salt: usize, modulus: usize) -> i32 {
    i32::try_from(1 + (i * 97 + salt) % modulus).expect("scattered position fits i32")
}

/// An unsorted input spanning `refs.len()` references, covering every record
/// shape the inline `.bai` must position-bin the same way samtools does:
/// normally placed mapped reads, placed-but-unmapped reads (`FUNMAP` set but
/// carrying a real ref+pos — e.g. the unmapped mate of a mapped read; both
/// fgumi and samtools bin these by position), and a trailing block of truly
/// unplaced reads (no reference at all — samtools' `n_no_coor`).
fn diverse_multiref_records(
    mapped_per_ref: usize,
    num_refs: usize,
    ref_len: usize,
) -> Vec<RawRecord> {
    let mut records = Vec::new();
    let mut b = SamBuilder::new();
    let modulus = ref_len.saturating_sub(100).max(1);

    for ref_idx in 0..num_refs {
        for i in 0..mapped_per_ref {
            let pos = scattered_pos(i, ref_idx * 31 + 13, modulus);
            b.clear();
            b.read_name(format!("m{ref_idx}_{i}").as_bytes())
                .ref_id(i32::try_from(ref_idx).expect("ref_id fits i32"))
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            records.push(b.build());
        }
    }

    // Placed-but-unmapped: FUNMAP set, real ref+pos.
    for ref_idx in 0..num_refs {
        for i in 0..(mapped_per_ref / 4).max(1) {
            let pos = scattered_pos(i, ref_idx * 7 + 5, modulus);
            b.clear();
            b.read_name(format!("pu{ref_idx}_{i}").as_bytes())
                .ref_id(i32::try_from(ref_idx).expect("ref_id fits i32"))
                .pos(pos)
                .mapq(0)
                .flags(flags::UNMAPPED)
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            records.push(b.build());
        }
    }

    // Trailing truly-unplaced region: no reference, no position.
    for i in 0..(mapped_per_ref / 8).max(50) {
        b.clear();
        b.read_name(format!("u{i}").as_bytes())
            .ref_id(-1)
            .pos(-1)
            .mapq(0)
            .flags(flags::UNMAPPED)
            .sequence(b"ACGT")
            .qualities(&[30u8; 4]);
        records.push(b.build());
    }

    records
}

/// `mapped_per_ref` scattered mapped reads on each of `num_refs` references —
/// enough, at default BGZF compression, to span many 65280-byte physical
/// blocks (Task 7 requirement 2). No unmapped reads: this generator's only
/// job is bulk plus one oversized record (below); category coverage is
/// `diverse_multiref_records`'s job.
fn many_mapped_records(mapped_per_ref: usize, num_refs: usize, ref_len: usize) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(mapped_per_ref * num_refs);
    let mut b = SamBuilder::new();
    let modulus = ref_len.saturating_sub(100).max(1);
    for ref_idx in 0..num_refs {
        for i in 0..mapped_per_ref {
            let pos = scattered_pos(i, ref_idx * 31 + 13, modulus);
            b.clear();
            b.read_name(format!("m{ref_idx}_{i}").as_bytes())
                .ref_id(i32::try_from(ref_idx).expect("ref_id fits i32"))
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            records.push(b.build());
        }
    }
    records
}

/// One record whose SEQ/QUAL/CIGAR are `len` bases long — well over one BGZF
/// physical block (`fgumi_bgzf::BGZF_MAX_BLOCK_SIZE` is 65280 bytes) — so the
/// record's own bytes straddle a physical-block boundary in the compressed
/// output (Task 7 requirement 2's "record spans physical blocks" case).
fn oversized_record(ref_id: i32, pos: i32, len: usize) -> RawRecord {
    let bases = vec![b'A'; len];
    let quals = vec![30u8; len];
    let mut b = SamBuilder::new();
    b.read_name(b"oversized")
        .ref_id(ref_id)
        .pos(pos)
        .mapq(60)
        .flags(0)
        .cigar_ops(&[u32::try_from(len).expect("record length fits u32 CIGAR op") << 4])
        .sequence(&bases)
        .qualities(&quals);
    b.build()
}

/// A `SortOptions` covering every field `add_sort` reads, matching
/// `Sort::to_sort_options()`'s field set (mirrors Task 6's test). Only
/// `max_memory` and `tmp_dirs` vary per test below.
fn sort_options(max_memory: MemoryLimit, tmp_dirs: Vec<PathBuf>) -> SortOptions {
    SortOptions {
        order: SortOrderArg::Coordinate,
        key_types: None,
        max_memory,
        memory_reserve: MemoryReserve::Auto,
        memory_per_thread: true,
        tmp_dirs,
        sort_threads: None,
        merge_threads: None,
        temp_compression: 1,
        temp_codec: fgumi_sort::SpillCodec::default(),
        max_temp_files: MaxTempFiles::Auto,
        block_batch: 4,
        file_granularity: false,
    }
}

/// A `[Stage::Sort]`-terminal `SinkSpec::BamWithIndex` chain spec — the exact
/// shape `SinkSpec::BamWithIndex` requires (`validate.rs` Rule 3) and the only
/// shape from which `ChainBuilder::add_sort` sets `detached_writer = true`
/// (`is_standalone_sort`, i.e. `stages == [Stage::Sort]`). See
/// [`arena_bam_with_index_detached_writer_is_deterministic_and_correct`]'s doc
/// comment for why that ties every chain-direct `BamWithIndex` test to the
/// Detached writer.
fn bam_with_index_spec(
    input: PathBuf,
    output: PathBuf,
    sort_options: SortOptions,
    threads: Option<usize>,
) -> ChainSpec {
    ChainSpec {
        stages: vec![Stage::Sort],
        source: SourceSpec::Bam(input),
        sink: SinkSpec::BamWithIndex(output),
        stage_opts: StageOptionsBag { sort: Some(sort_options), ..Default::default() },
        threading: ThreadingOptions { threads },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: true,
        command_line: "fgumi sort --write-index".to_string(),
    }
}

/// `samtools view -c [extra] <bam> <region>` region query via the `.bai`.
/// `extra` carries filter flags (e.g. `-F 4`).
///
/// Returns the records rather than a count: a matching count alone cannot
/// distinguish "returned the right records" from "returned the right
/// *number* of the wrong records," which a mis-recovered virtual offset or a
/// bin pointing at a neighbouring chunk would produce identically. Mirrors
/// `tests/integration/test_sort_write_index.rs`'s helper of the same name.
fn region_records(bam: &Path, region: &str, extra: &[&str]) -> String {
    let out = Command::new("samtools")
        .arg("view")
        .args(extra)
        .args([bam.to_str().unwrap(), region])
        .output()
        .expect("run samtools view");
    assert!(
        out.status.success(),
        "samtools view {region} failed on {}: {}",
        bam.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// Assert two region-query outputs hold the identical records, reporting the
/// first divergence compactly (these regions can return thousands of records,
/// so a whole-output `assert_eq!` diff would be unreadable).
fn assert_same_records(region: &str, via_fgumi: &str, via_samtools: &str) {
    if via_fgumi == via_samtools {
        return;
    }
    let fgumi_lines: Vec<&str> = via_fgumi.lines().collect();
    let samtools_lines: Vec<&str> = via_samtools.lines().collect();
    let detail = match fgumi_lines.iter().zip(samtools_lines.iter()).position(|(a, b)| a != b) {
        Some(i) => format!(
            "first differing record at index {i}:\n  via fgumi:    {}\n  via samtools: {}",
            fgumi_lines[i], samtools_lines[i]
        ),
        None => "one output is a strict prefix of the other".to_string(),
    };
    panic!(
        "region {region}: records retrieved via fgumi's inline .bai differ from samtools' .bai \
         ({} vs {} records); {detail}",
        fgumi_lines.len(),
        samtools_lines.len()
    );
}

/// `samtools idxstats <bam>` — per-reference mapped/unmapped counts (plus the
/// trailing no-coordinate count) from the `.bai`.
fn idxstats(bam: &Path) -> String {
    let out = Command::new("samtools")
        .args(["idxstats", bam.to_str().unwrap()])
        .output()
        .expect("run samtools idxstats");
    assert!(
        out.status.success(),
        "samtools idxstats failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// Parse the inline `.bai` sidecar as a well-formed BAI and assert it carries
/// one reference sequence's worth of bins per header reference.
///
/// Called unconditionally — before any `samtools_available()` gate — so the
/// chain, the sidecar write, and the BAI parser are exercised on every run,
/// even the samtools-less CI job that only runs `#[ignore]`d tests. The
/// samtools parity assertions each test adds on top are what the gate skips.
fn assert_inline_bai_valid(bai_path: &Path, expected_refs: usize) {
    let index =
        noodles::bam::bai::fs::read(bai_path).expect("inline .bai must parse as a valid BAI index");
    assert_eq!(
        index.reference_sequences().len(),
        expected_refs,
        "inline .bai reference-sequence count must match the header",
    );
}

/// Task 7 requirement 1: the inline `.bai` `SinkSpec::BamWithIndex` writes
/// must be idxstats/region-query equivalent to one samtools builds for the
/// identical output bytes, across the record shapes samtools and fgumi must
/// bin the same way — placed mapped reads, placed-but-unmapped reads, and a
/// trailing truly-unplaced region — spread over three references.
#[test]
fn arena_bam_with_index_idxstats_equivalent_to_samtools() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    let refs = [("chr1", 200_000usize), ("chr2", 200_000usize), ("chr3", 200_000usize)];
    let header = multi_ref_header(&refs);
    let records = diverse_multiref_records(800, refs.len(), 200_000);
    write_bam(&input_bam, &header, &records);

    let spec = bam_with_index_spec(
        input_bam,
        output_bam.clone(),
        sort_options(MemoryLimit::Fixed(64 * 1024 * 1024), Vec::new()),
        None,
    );

    build_for(spec)
        .expect("build_for should accept a Sort-terminal BamWithIndex chain")
        .run()
        .expect("running the chain should write the sorted BAM and its inline .bai");

    assert!(output_bam.exists(), "sorted output BAM was not written");
    let fgumi_bai = fgumi_bam_io::bai_sidecar_path(&output_bam);
    assert!(fgumi_bai.exists(), "inline indexer did not write a .bai sidecar at {fgumi_bai:?}");
    // Exercise the chain, sidecar, and parser regardless of samtools.
    assert_inline_bai_valid(&fgumi_bai, refs.len());

    if !samtools_available() {
        eprintln!(
            "arena_bam_with_index_idxstats_equivalent_to_samtools: samtools not on PATH; \
             validated the inline .bai only, skipped the samtools parity assertions"
        );
        return;
    }

    // A samtools-built index over an independent copy of the identical output
    // bytes: this isolates the two indexes (each colocated with its own copy
    // of the BAM) rather than overwriting fgumi's inline sidecar in place.
    let reference = dir.path().join("reference.bam");
    std::fs::copy(&output_bam, &reference).expect("copy output bam for reference indexing");
    let status = Command::new("samtools")
        .args(["index", reference.to_str().unwrap()])
        .status()
        .expect("run samtools index");
    assert!(status.success(), "samtools index failed");

    let regions = [
        "chr1",
        "chr2",
        "chr3",
        "chr1:1-100000",
        "chr2:50000-150000",
        "chr3:150000-200000",
        "chr1:199000-200000",
    ];
    let mut total = 0usize;
    for r in regions {
        let via_fgumi = region_records(&output_bam, r, &[]);
        let via_samtools = region_records(&reference, r, &[]);
        assert_same_records(r, &via_fgumi, &via_samtools);
        total += via_fgumi.lines().count();
    }
    assert!(total > 0, "expected non-empty region queries");

    assert_eq!(
        idxstats(&output_bam),
        idxstats(&reference),
        "idxstats via the inline .bai must match samtools' idxstats over the identical bytes"
    );
}

/// Task 7 requirement 2: force many 65280-byte physical BGZF blocks (a small
/// in-memory budget forces many spill runs and a multi-block k-way merge)
/// plus one record whose own bytes are larger than one physical block, and
/// assert the inline `.bai` still answers region queries identically to
/// samtools — proving the join between the merge's recovered virtual offsets
/// and `BgzfCompress`'s per-record manifest survives both a batch boundary
/// that doesn't line up with a physical block AND a record straddling a
/// physical block by itself.
#[test]
fn arena_bam_with_index_multi_block_and_straddle() {
    let dir = TempDir::new().expect("create temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    let refs = [("chr1", 2_000_000usize), ("chr2", 2_000_000usize)];
    let header = multi_ref_header(&refs);
    let mut records = many_mapped_records(12_000, refs.len(), 2_000_000);
    records.push(oversized_record(0, 100_000, 70_000)); // > BGZF_MAX_BLOCK_SIZE (65280)
    write_bam(&input_bam, &header, &records);

    // A 64 KiB in-memory budget over ~24k records forces many real spill
    // files, so both the spill k-way merge and the sink's compression cross
    // physical-block boundaries repeatedly. Spill temp files go to the test's
    // own disk-backed TempDir (`--tmp-dir`), not the default tmp location.
    let spec = bam_with_index_spec(
        input_bam,
        output_bam.clone(),
        sort_options(MemoryLimit::Fixed(64 * 1024), vec![dir.path().to_path_buf()]),
        Some(4),
    );

    build_for(spec)
        .expect("build_for should accept a Sort-terminal BamWithIndex chain")
        .run()
        .expect("running the chain should write the sorted BAM and its inline .bai");

    assert!(output_bam.exists(), "sorted output BAM was not written");
    let fgumi_bai = fgumi_bam_io::bai_sidecar_path(&output_bam);
    assert!(fgumi_bai.exists(), "inline indexer did not write a .bai sidecar at {fgumi_bai:?}");
    // Exercise the chain, sidecar, and parser regardless of samtools.
    assert_inline_bai_valid(&fgumi_bai, refs.len());

    if !samtools_available() {
        eprintln!(
            "arena_bam_with_index_multi_block_and_straddle: samtools not on PATH; validated the \
             inline .bai only, skipped the samtools parity assertions"
        );
        return;
    }

    let reference = dir.path().join("reference.bam");
    std::fs::copy(&output_bam, &reference).expect("copy output bam for reference indexing");
    let status = Command::new("samtools")
        .args(["index", reference.to_str().unwrap()])
        .status()
        .expect("run samtools index");
    assert!(status.success(), "samtools index failed");

    let regions = [
        "chr1",
        "chr2",
        "chr1:1-500000",
        "chr1:90000-180000", // spans the oversized record
        "chr1:1500000-2000000",
        "chr2:1-2000000",
    ];
    let mut total = 0usize;
    for r in regions {
        let via_fgumi = region_records(&output_bam, r, &[]);
        let via_samtools = region_records(&reference, r, &[]);
        assert_same_records(r, &via_fgumi, &via_samtools);
        total += via_fgumi.lines().count();
    }
    assert!(total > 0, "expected non-empty region queries");

    assert_eq!(
        idxstats(&output_bam),
        idxstats(&reference),
        "idxstats via the inline .bai must match samtools' idxstats over the identical bytes"
    );
}

/// Task 7 requirement 3: Detached-vs-Serial writer byte-identity.
///
/// **No seam exists to force the Serial (non-detached) writer from a
/// chain-direct test, so this asserts the achievable subset instead of
/// faking the comparison — per the ruling's explicit fallback.**
///
/// `ChainBuilder::detached_writer` (private,
/// `src/lib/pipeline/chains/builder.rs`) is set to `true` only when
/// `self.spec.is_sort_terminal()` — i.e. `spec.stages == [Stage::Sort]`
/// exactly — and the Sort stage sits at the chain's terminal position
/// (`add_sort`, guarded by `is_standalone_sort`). `SinkSpec::BamWithIndex`
/// only requires `stages.last() == Some(Stage::Sort)` (`validate.rs` Rule 3),
/// which is *weaker* than `is_sort_terminal()` — so in principle a
/// `BamWithIndex` chain with a leading non-Sort stage (e.g.
/// `[Stage::Correct, Stage::Sort]`) would leave `detached_writer` at its
/// default `false` (Serial). But inserting a real upstream stage changes what
/// is being compared: `Correct` (or any other stage) decodes and rewrites
/// every record before Sort ever sees it, so a byte-identity check between
/// that chain's output and a plain `[Stage::Sort]` chain's output would be
/// confounded by the upstream stage's own re-encoding, not isolate the
/// writer-threading variable Task 7 asks about. No `Stage` is a verified
/// byte-preserving no-op, and inventing a test-only production seam to flip
/// `detached_writer` directly from outside `builder.rs` was out of scope for
/// a correctness-gate task (see `task-7-report.md` for the full reasoning).
///
/// What IS asserted here: every `SinkSpec::BamWithIndex` chain built
/// chain-direct takes the Detached writer, so this runs that (only reachable)
/// configuration twice over the same input and asserts (a) the two runs are
/// byte-identical to each other — the Detached writer is deterministic, not
/// just "correct once" — and (b) the produced `.bai` is idxstats-equivalent
/// to samtools, re-confirming correctness for this configuration
/// independently of the two tests above. This is confirmatory, not
/// load-bearing, for the Detached/Serial claim itself: Task 5's review
/// already reasoned that `with_detached()` only relocates which thread drives
/// the sink without touching the manifest/join math that produces the index
/// bytes, so Serial's BAM/`.bai` bytes are expected — not proven here — to be
/// identical to Detached's.
#[test]
fn arena_bam_with_index_detached_writer_is_deterministic_and_correct() {
    let dir = TempDir::new().expect("create temp dir");
    let refs = [("chr1", 100_000usize), ("chr2", 100_000usize)];
    let header = multi_ref_header(&refs);
    let records = diverse_multiref_records(300, refs.len(), 100_000);
    let input_bam = dir.path().join("in.bam");
    write_bam(&input_bam, &header, &records);

    let run = |name: &str| -> PathBuf {
        let output_bam = dir.path().join(name);
        let spec = bam_with_index_spec(
            input_bam.clone(),
            output_bam.clone(),
            sort_options(MemoryLimit::Fixed(64 * 1024 * 1024), Vec::new()),
            None,
        );
        build_for(spec)
            .expect("build_for should accept a Sort-terminal BamWithIndex chain")
            .run()
            .expect("running the chain should write the sorted BAM and its inline .bai");
        output_bam
    };

    let bam_run_a = run("run_a.bam");
    let bam_run_b = run("run_b.bam");
    let index_run_a = fgumi_bam_io::bai_sidecar_path(&bam_run_a);
    let index_run_b = fgumi_bam_io::bai_sidecar_path(&bam_run_b);

    assert_eq!(
        std::fs::read(&bam_run_a).expect("read run_a.bam"),
        std::fs::read(&bam_run_b).expect("read run_b.bam"),
        "two runs of the arena chain's Detached-writer sink over the same input must produce \
         byte-identical BAM output"
    );
    assert_eq!(
        std::fs::read(&index_run_a).expect("read run_a.bam.bai"),
        std::fs::read(&index_run_b).expect("read run_b.bam.bai"),
        "two runs of the arena chain's Detached-writer sink over the same input must produce a \
         byte-identical inline .bai"
    );
    // The deterministic byte-identity checks above need no samtools; the sidecar
    // parse also runs regardless of it.
    assert_inline_bai_valid(&index_run_a, refs.len());

    if !samtools_available() {
        eprintln!(
            "arena_bam_with_index_detached_writer_is_deterministic_and_correct: samtools not on \
             PATH; validated determinism and the inline .bai only, skipped the samtools parity \
             assertions"
        );
        return;
    }

    let reference = dir.path().join("reference.bam");
    std::fs::copy(&bam_run_a, &reference).expect("copy output bam for reference indexing");
    let status = Command::new("samtools")
        .args(["index", reference.to_str().unwrap()])
        .status()
        .expect("run samtools index");
    assert!(status.success(), "samtools index failed");
    assert_eq!(
        idxstats(&bam_run_a),
        idxstats(&reference),
        "idxstats via the inline .bai must match samtools' idxstats over the identical bytes"
    );
}
