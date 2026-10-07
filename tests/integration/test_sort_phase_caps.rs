//! The thread-flag contract at the chain level: an explicit `--sort-threads` /
//! `--merge-threads` becomes one shared admission cap per phase, visible on
//! `Pipeline::dag()`, bounded by the pool the chain runs (`--threads`, or a
//! zipper / aligner floor above it), and absent from steps the sort does not
//! own. An unset flag creates no cap at all.

use std::ffi::OsStr;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use fgumi_lib::commands::common::{
    CompressionOptions, MaxTempFiles, MemoryLimit, MemoryReserve, QueueMemoryOptions,
    SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::group::GroupOptions;
use fgumi_lib::commands::sort::{SortOptions, SortOrderArg};
use fgumi_lib::commands::zipper::ZipperOptions;
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_test_reference, write_bam};
use crate::helpers::cutover::decompressed_records_without_pg;

/// `n` records (whole FR pairs) in reverse order of every sort key: template
/// names count down (unsorted for queryname), positions count down (unsorted
/// for coordinate and template-coordinate). Every three templates share one
/// position and strand, so equal-key ties — where a thread count could
/// conceivably change record order — are part of the byte-identity check.
fn unsorted_records(n: usize) -> Vec<RawRecord> {
    let templates = n / 2;
    let mut records = Vec::with_capacity(templates * 2);
    for t in 0..templates {
        let id = templates - t;
        let pos = i32::try_from((id / 3 + 1) * 100).expect("pos fits i32");
        let mate_pos = pos + 50;
        for first in [true, false] {
            let (own, mate, strand) = if first {
                (pos, mate_pos, flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            } else {
                (mate_pos, pos, flags::LAST_SEGMENT | flags::REVERSE)
            };
            let tlen = if first { 54 } else { -54 };
            let mut b = SamBuilder::new();
            b.read_name(format!("read{id}").as_bytes())
                .ref_id(0)
                .pos(own)
                .mapq(60)
                .flags(flags::PAIRED | flags::PROPER_PAIR | strand)
                .mate_ref_id(0)
                .mate_pos(mate)
                .template_length(tlen)
                .cigar_ops(&[4u32 << 4])
                .sequence(b"ACGT")
                .qualities(&[30u8; 4])
                .add_string_tag(SamTag::MC, b"4M");
            records.push(b.build());
        }
    }
    records
}

fn sort_options(sort_threads: Option<usize>, merge_threads: Option<usize>) -> SortOptions {
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
        block_batch: 4,
        file_granularity: false,
        sort_stats: false,
    }
}

fn spec_with(
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
        command_line: "fgumi sort <phase-cap test>".to_string(),
    }
}

fn spec(
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

/// The `cap=` token on the DAG line of `step`, if any.
fn cap_of(dag: &str, step: &str) -> Option<String> {
    let line = dag
        .lines()
        .find(|l| l.split_whitespace().nth(1) == Some(step))
        .unwrap_or_else(|| panic!("no DAG line for {step}:\n{dag}"));
    line.split_whitespace().find_map(|tok| tok.strip_prefix("cap=").map(str::to_owned))
}

/// Every WARN line about the per-phase thread flags — matched on the level and
/// the flag name, not the wording, so a reworded warning cannot make a
/// "no warning" check pass vacuously.
fn thread_flag_warnings(stderr: &str) -> Vec<&str> {
    stderr
        .lines()
        .filter(|l| {
            l.contains(" WARN ") && (l.contains("sort-threads") || l.contains("merge-threads"))
        })
        .collect()
}

/// Standalone sort: with an explicit `--sort-threads`, every phase-1 step —
/// the pool steps and the per-run sort (`FindBoundariesAndSort`) — shares
/// `sort-phase1(min(sort,threads))`; with an explicit `--merge-threads`, every
/// phase-2 step — the pool steps and the fast-path gather (`SortMerge`) —
/// shares `sort-phase2(min(merge,threads))`. An unset flag caps nothing in its
/// phase. Capping an unset phase, or not handing the cap to the per-run sort or
/// the gather, fails a case here. Whether those two then actually take the
/// whole cap is pinned by the `fgumi-pipeline-io` unit tests
/// (`capped_parallel_seal_holds_the_whole_cap_while_it_sorts`,
/// `capped_parallel_fast_path_gather_holds_the_whole_phase2_cap`).
#[rstest]
#[case::defaults(16, None, None, None, None)]
#[case::sort_below(16, Some(4), None, Some("sort-phase1(4)"), None)]
#[case::merge_below(16, None, Some(2), None, Some("sort-phase2(2)"))]
#[case::sort_above_clamps(16, Some(32), None, Some("sort-phase1(16)"), None)]
#[case::zero_is_one(4, Some(0), Some(0), Some("sort-phase1(1)"), Some("sort-phase2(1)"))]
fn standalone_sort_caps_every_step_of_an_explicit_phase(
    #[case] threads: usize,
    #[case] sort_threads: Option<usize>,
    #[case] merge_threads: Option<usize>,
    #[case] phase1: Option<&str>,
    #[case] phase2: Option<&str>,
) {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &unsorted_records(100));

    let built = build_for(spec(
        vec![Stage::Sort],
        &input,
        &output,
        threads,
        sort_options(sort_threads, merge_threads),
    ))
    .expect("standalone sort chain builds");
    let dag = built.pipeline.dag();

    for step in ["InflateToArena", "FindBoundariesAndSort", "SpillBlockCompress"] {
        assert_eq!(cap_of(&dag, step).as_deref(), phase1, "{step}:\n{dag}");
    }
    for step in ["SortSpillDecompress", "SortMerge", "BgzfCompress"] {
        assert_eq!(cap_of(&dag, step).as_deref(), phase2, "{step}:\n{dag}");
    }
    if phase1.is_none() && phase2.is_none() {
        assert!(built.pipeline.phase_caps().is_empty(), "unset flags create no cap:\n{dag}");
    }
}

/// The cap belongs to the phase, not the step: in a built standalone sort every
/// phase-1 step — including the per-run sort, `FindBoundariesAndSort` — holds
/// the SAME counter, as does every phase-2 step including the fast-path gather
/// in `SortMerge`, and the phases' counters are distinct. Identity is checked
/// by address and by behaviour — exhausting the counter through one step's cap
/// makes the other step's cap refuse.
#[test]
fn standalone_sort_steps_of_one_phase_share_one_counter() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &unsorted_records(100));

    let built =
        build_for(spec(vec![Stage::Sort], &input, &output, 8, sort_options(Some(2), Some(3))))
            .expect("standalone sort chain builds");
    let caps = built.pipeline.phase_caps();
    let cap = |step: &str| {
        caps.iter()
            .find(|(name, _)| *name == step)
            .map_or_else(|| panic!("{step} has no phase cap"), |(_, cap)| *cap)
    };
    let (inflate, spill) = (cap("InflateToArena"), cap("SpillBlockCompress"));
    let (decompress, output_compress) = (cap("SortSpillDecompress"), cap("BgzfCompress"));

    assert!(std::ptr::eq(inflate, spill), "phase-1 steps must share one counter");
    assert!(std::ptr::eq(decompress, output_compress), "phase-2 steps must share one counter");
    assert!(
        std::ptr::eq(cap("FindBoundariesAndSort"), inflate),
        "the per-run sort must hold the phase-1 pool steps' counter"
    );
    assert!(
        std::ptr::eq(cap("SortMerge"), decompress),
        "the fast-path gather must hold the phase-2 pool steps' counter"
    );
    assert!(!std::ptr::eq(inflate, decompress), "the phases have separate counters");

    let held: Vec<_> = (0..2).map(|_| inflate.try_acquire().expect("phase-1 permit")).collect();
    assert!(spill.try_acquire().is_none(), "spill compress shares inflate's exhausted counter");
    assert!(decompress.try_acquire().is_some(), "phase 2 is unaffected by phase 1");
    drop(held);
}

/// A fused chain's sort is intermediate: its spill steps are still capped, but
/// the terminal `BgzfCompress` compresses the *group* output and must stay
/// uncapped — `detached_writer` is false, so `add_sink` never receives a cap,
/// and an uncapped step shows no `cap=` token at all.
#[test]
fn runall_shaped_chain_leaves_terminal_compress_uncapped() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &unsorted_records(100));

    let mut sort = sort_options(Some(2), None);
    sort.order = SortOrderArg::TemplateCoordinate;
    let built = build_for(spec(vec![Stage::Sort, Stage::Group], &input, &output, 8, sort))
        .expect("fused sort+group chain builds");
    let dag = built.pipeline.dag();

    assert_eq!(cap_of(&dag, "InflateToArena").as_deref(), Some("sort-phase1(2)"), "{dag}");
    assert_eq!(cap_of(&dag, "SortSpillDecompress"), None, "--merge-threads unset: {dag}");
    assert_eq!(cap_of(&dag, "BgzfCompress"), None, "{dag}");
}

/// Matched unmapped + mapped BAMs and a reference for a zipper start: four
/// templates, enough to build and run `zipper → sort`.
struct ZipperInputs {
    unmapped: PathBuf,
    mapped: PathBuf,
    reference: PathBuf,
}

fn zipper_inputs(dir: &Path) -> ZipperInputs {
    let unmapped = dir.join("unmapped.bam");
    let mapped = dir.join("mapped.bam");
    let reference = create_test_reference(dir);
    let names = ["read1", "read2", "read3", "read4"];
    let unmapped_records: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::UNMAPPED)
                .add_string_tag(SamTag::RX, format!("AACC{i}").as_bytes());
            b.build()
        })
        .collect();
    write_bam(&unmapped, &noodles::sam::Header::default(), &unmapped_records);
    let mapped_records: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(i32::try_from(400 - 50 * i).expect("pos fits i32"))
                .mapq(60)
                .flags(0)
                .cigar_ops(&[8u32 << 4])
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8]);
            b.build()
        })
        .collect();
    write_bam(&mapped, &create_minimal_header("chr1", 10_000), &mapped_records);
    ZipperInputs { unmapped, mapped, reference }
}

/// Zipper raises the pool to at least 4 workers, so an explicit sort cap
/// follows that pool, not the raw `--threads`: `--sort::sort-threads 4` with
/// `--threads 2` is honoured, not clamped. Unset, no step is capped. The
/// zipper's template output reaches the sort through `TemplatesToRecordBatch`,
/// which is sort ingest and so shares the phase-1 cap.
#[rstest]
#[case::threads_unset(None, None, None)]
#[case::threads_two_sort_four(Some(2), Some(4), Some("sort-phase1(4)"))]
#[case::threads_two_sort_two(Some(2), Some(2), Some("sort-phase1(2)"))]
fn zipper_floor_widens_the_sort_caps(
    #[case] threads: Option<usize>,
    #[case] sort_threads: Option<usize>,
    #[case] phase1: Option<&str>,
) {
    let dir = TempDir::new().expect("tempdir");
    let inputs = zipper_inputs(dir.path());
    let output = dir.path().join("out.bam");
    let mut sort = sort_options(sort_threads, None);
    sort.order = SortOrderArg::TemplateCoordinate;
    let bag = StageOptionsBag {
        zipper: Some(ZipperOptions::default()),
        sort: Some(sort),
        ..Default::default()
    };
    let source = SourceSpec::PairedBams {
        unmapped: inputs.unmapped,
        mapped: inputs.mapped,
        reference: inputs.reference,
    };
    let built =
        build_for(spec_with(vec![Stage::Zipper, Stage::Sort], source, &output, threads, bag))
            .expect("zipper→sort chain builds");
    let dag = built.pipeline.dag();
    assert_eq!(cap_of(&dag, "TemplatesToRecordBatch").as_deref(), phase1, "{dag}");
    assert_eq!(cap_of(&dag, "SpillBlockCompress").as_deref(), phase1, "{dag}");
    if phase1.is_some() {
        let caps = built.pipeline.phase_caps();
        let cap = |step: &str| {
            caps.iter()
                .find(|(name, _)| *name == step)
                .map_or_else(|| panic!("{step} has no phase cap"), |(_, cap)| *cap)
        };
        let spill = cap("SpillBlockCompress");
        for step in ["TemplatesToRecordBatch", "SortBuffer"] {
            assert!(std::ptr::eq(cap(step), spill), "{step} must share the phase-1 counter");
        }
    }
    assert_eq!(cap_of(&dag, "SortSpillDecompress"), None, "--merge-threads unset: {dag}");
    // Sort is the last stage but not the only one, so the terminal output
    // compressor is not phase-2 work (matches `--sort::merge-threads` help).
    assert_eq!(cap_of(&dag, "BgzfCompress"), None, "{dag}");
    built.run().expect("zipper→sort chain runs");
}

/// End to end through runall: a zipper→sort run at `--threads 2` with
/// `--sort::sort-threads 4` warns about nothing (the pool has 4 workers), the
/// debug split shows phase 1 at 4, and the memory budget multiplier is the
/// same as `origin/main`'s: 1 with nothing set (`--threads` unset), 4 for the
/// explicit `--sort::sort-threads 4`.
#[rstest]
#[case::threads_unset(&[], 1)]
#[case::threads_two_sort_four(&["--threads", "2", "--sort::sort-threads", "4"], 4)]
fn runall_zipper_sort_caps_follow_the_pool_without_a_false_warning(
    #[case] extra: &[&str],
    #[case] budget_threads: usize,
) {
    let dir = TempDir::new().expect("tempdir");
    let inputs = zipper_inputs(dir.path());
    let output = dir.path().join("out.bam");
    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "debug")
        .args(["runall", "--start-from", "zipper", "--stop-after", "sort"])
        .arg("-i")
        .arg(&inputs.mapped)
        .arg("--unmapped")
        .arg(&inputs.unmapped)
        .arg("--ref")
        .arg(&inputs.reference)
        .arg("-o")
        .arg(&output)
        .args(extra)
        .output()
        .expect("run fgumi runall");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(result.status.success(), "runall zipper→sort failed:\n{stderr}");
    assert_eq!(thread_flag_warnings(&stderr), Vec::<&str>::new(), "no warning in a 4-worker pool");
    let split = stderr
        .lines()
        .find(|l| l.contains("per-phase thread split"))
        .unwrap_or_else(|| panic!("no per-phase split logged:\n{stderr}"));
    assert!(split.contains("(phase1=4, phase2=4) within a 4-worker pool"), "{split}");
    assert!(split.contains(&format!("memory budget x{budget_threads}")), "{split}");
}

/// Run `fgumi sort` as a subprocess under a tiny budget (so phase 1 spills and
/// phase 2 merges — both caps do real work) and return (output, stderr).
fn sort_subprocess(
    input: &Path,
    dir: &Path,
    name: &str,
    order: &str,
    extra: &[&str],
) -> (PathBuf, String) {
    let output = dir.join(format!("{name}.bam"));
    let tmp = dir.join(format!("{name}-tmp"));
    std::fs::create_dir_all(&tmp).expect("tmp dir");
    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args([
            OsStr::new("sort"),
            OsStr::new("-i"),
            input.as_os_str(),
            OsStr::new("-o"),
            output.as_os_str(),
        ])
        .args([
            "--order",
            order,
            "--max-memory",
            "64K",
            "--memory-per-thread",
            "false",
            "--temp-compression",
            "1",
        ])
        .arg("--tmp-dir")
        .arg(&tmp)
        .args(extra)
        .output()
        .expect("run fgumi sort");
    let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(result.status.success(), "fgumi sort {name} failed:\n{stderr}");
    // A real k-way merge, not one run passed through: the merge must read
    // more than one source.
    let sources: usize = stderr
        .lines()
        .find_map(|l| l.split_once("Merge sources: ").map(|(_, n)| n.trim().to_owned()))
        .unwrap_or_else(|| panic!("{name}: no `Merge sources:` line:\n{stderr}"))
        .parse()
        .expect("merge source count");
    assert!(sources > 1, "{name} must merge several runs so both phases work:\n{stderr}");
    (output, stderr)
}

/// Capped and clamped runs are byte-identical to the plain `--threads 4` run
/// (ties included), the clamped runs warn once per flag, and the summary
/// reports each phase cap's high-water mark within its limit. The peaks are a
/// sanity check of the reporting, not proof that each step honours its cap —
/// the per-step `full_cap_refuses_*` unit tests pin that, and the
/// `fgumi-pipeline-io` whole-cap unit tests pin that the per-run sort and the
/// gather hold the whole cap while they run.
#[rstest]
#[case::caps_of_one(&["--threads", "4", "--sort-threads", "1", "--merge-threads", "1"], 1, 1, 0)]
#[case::clamped_above_threads(&["--threads", "2", "--sort-threads", "8", "--merge-threads", "8"], 2, 2, 2)]
#[case::threads_one_clamped(&["--threads", "1", "--sort-threads", "4", "--merge-threads", "4"], 1, 1, 2)]
#[case::sort_cap_only(&["--threads", "4", "--sort-threads", "2"], 2, 0, 0)]
fn capped_and_clamped_sorts_are_byte_identical_and_warn(
    #[case] extra: &[&str],
    #[case] phase1_cap: usize,
    #[case] phase2_cap: usize,
    #[case] expected_warnings: usize,
) {
    for order in ["coordinate", "template-coordinate", "queryname::natural"] {
        let dir = TempDir::new().expect("tempdir");
        let input = dir.path().join("in.bam");
        write_bam(&input, &create_minimal_header("chr1", 1_000_000), &unsorted_records(20_000));

        let (reference, _) =
            sort_subprocess(&input, dir.path(), "reference", order, &["--threads", "4"]);
        let (capped, stderr) = sort_subprocess(&input, dir.path(), "capped", order, extra);

        assert_eq!(
            decompressed_records_without_pg(&capped),
            decompressed_records_without_pg(&reference),
            "{order}: capped run must be byte-identical to --threads 4"
        );
        let warnings = thread_flag_warnings(&stderr).len();
        assert_eq!(warnings, expected_warnings, "{order}: clamp warnings\n{stderr}");

        let caps_line = stderr
            .lines()
            .find(|l| l.contains("Phase caps"))
            .unwrap_or_else(|| panic!("{order}: summary must report phase caps:\n{stderr}"));
        // Each explicit phase is reported with its configured limit and a
        // peak within it; an unset phase (limit 0 here) has no cap to report.
        for (phase, limit) in [("sort-phase1", phase1_cap), ("sort-phase2", phase2_cap)] {
            let Some(after) = caps_line.split(phase).nth(1) else {
                assert_eq!(limit, 0, "{order}: {phase} missing: {caps_line}");
                continue;
            };
            assert_ne!(limit, 0, "{order}: {phase} has no flag but was reported: {caps_line}");
            let peak: usize =
                after.split_whitespace().nth(1).expect("peak").parse().expect("peak is a number");
            assert!((1..=limit).contains(&peak), "{order}: {caps_line}");
            assert!(caps_line.contains(&format!("{phase} peak {peak} of {limit},")), "{caps_line}");
        }
    }
}

/// A SAM-sourced sort (the `bwa mem | fgumi sort -i -` case) ingests through
/// Parallel `ParseSamChunk` and `DecodedRecordBatchToRecordBatch`, then sorts
/// in `SortBuffer`: all of it is phase-1 work, so every one of those steps
/// shares the single phase-1 counter with spill compression.
#[test]
fn sam_sourced_sort_caps_every_ingest_step_with_the_phase1_counter() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.sam");
    let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n");
    for pos in [500, 100, 400, 200, 300] {
        writeln!(sam, "q{pos}\t0\tchr1\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII").expect("write SAM");
    }
    std::fs::write(&input, sam).expect("write SAM");
    let output = dir.path().join("out.bam");

    let built = build_for(spec(vec![Stage::Sort], &input, &output, 8, sort_options(Some(3), None)))
        .expect("SAM-sourced sort chain builds");
    let dag = built.pipeline.dag();
    for step in
        ["ParseSamChunk", "DecodedRecordBatchToRecordBatch", "SortBuffer", "SpillBlockCompress"]
    {
        assert_eq!(cap_of(&dag, step).as_deref(), Some("sort-phase1(3)"), "{step}:\n{dag}");
    }
    let caps = built.pipeline.phase_caps();
    let cap = |step: &str| {
        caps.iter()
            .find(|(name, _)| *name == step)
            .map_or_else(|| panic!("{step} has no phase cap"), |(_, cap)| *cap)
    };
    let parse = cap("ParseSamChunk");
    for step in ["DecodedRecordBatchToRecordBatch", "SortBuffer", "SpillBlockCompress"] {
        assert!(std::ptr::eq(parse, cap(step)), "{step} must share the parse step's counter");
    }
    built.run().expect("SAM-sourced sort runs");
}

/// The record-buffer budget the BUILT chain uses follows `--sort-threads`: an
/// input that fits a 4-thread budget sorts without spilling at `--threads 4`,
/// but spills once `--sort-threads 1` shrinks the budget to one thread's worth
/// — and the output is byte-identical either way.
#[test]
fn sort_threads_below_threads_shrinks_the_record_buffer_in_the_chain() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &unsorted_records(20_000));
    let run = |name: &str, extra: &[&str]| {
        let output = dir.path().join(format!("{name}.bam"));
        let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .env("RUST_LOG", "info")
            .args([OsStr::new("sort"), OsStr::new("-i"), input.as_os_str()])
            .args([OsStr::new("-o"), output.as_os_str()])
            .args(["--order", "coordinate", "--max-memory", "1M", "--threads", "4"])
            .arg("--tmp-dir")
            .arg(dir.path())
            .args(extra)
            .output()
            .expect("run fgumi sort");
        let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
        assert!(result.status.success(), "fgumi sort {name} failed:\n{stderr}");
        // No spill logs no `Spill runs:` line at all.
        let spills: usize = stderr
            .lines()
            .find_map(|l| {
                l.split_once("Spill runs: ").and_then(|(_, n)| n.split_whitespace().next())
            })
            .map_or(0, |n| n.parse().expect("spill count"));
        (output, spills, stderr)
    };
    let (wide, wide_spills, wide_log) = run("wide", &[]);
    let (narrow, narrow_spills, narrow_log) = run("narrow", &["--sort-threads", "1"]);
    assert_eq!(wide_spills, 0, "a 4 x 1M budget holds the input:\n{wide_log}");
    assert!(narrow_spills >= 1, "a 1 x 1M budget must spill:\n{narrow_log}");
    assert_eq!(
        decompressed_records_without_pg(&narrow),
        decompressed_records_without_pg(&wide),
        "the smaller budget changes spilling, not output"
    );
}
