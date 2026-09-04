//! End-to-end tests for the clip command.
//!
//! These tests invoke `Clip::execute()` in-process and validate:
//! 1. Basic read clipping
//! 2. Fixed clipping (5' and 3' ends)
//! 3. Metrics output

use clap::Parser;
use fgumi_lib::commands::clip::Clip;
use fgumi_lib::commands::command::Command;
use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use rstest::rstest;
use std::fs;
use std::path::Path;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_coordinate_sorted_header, create_minimal_header, create_test_reference, to_record_buf,
};

/// Create a BAM with paired reads using the shared minimal (query-grouped) header.
fn create_paired_bam(path: &Path, pairs: Vec<(RawRecord, RawRecord)>) {
    let header = create_minimal_header("chr1", 10000);
    write_paired_bam_with_header(path, &header, pairs);
}

/// Write paired reads under a caller-supplied header (used to exercise
/// header-less and coordinate-sorted input).
fn write_paired_bam_with_header(
    path: &Path,
    header: &noodles::sam::Header,
    pairs: Vec<(RawRecord, RawRecord)>,
) {
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    for (r1, r2) in pairs {
        writer.write_alignment_record(header, &to_record_buf(&r1)).expect("Failed to write R1");
        writer.write_alignment_record(header, &to_record_buf(&r2)).expect("Failed to write R2");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

// ============================================================================
// --check-crc / --no-check-crc on the chain (`--threads N`) path (#800)
//
// clip always decodes its input through the unified pipeline (the chain is its
// only execution path). These prove that path honors the flag: `--no-check-crc`
// accepts a corrupted block, while the default and `--check-crc` reject it.
// ============================================================================

/// Flip a byte in the last BGZF block's CRC32 footer, so decoding that block
/// fails only when CRC verification is on. The input must span at least two
/// blocks so the corrupted block is not the header's block (which always
/// verifies during the header parse).
fn corrupt_last_block_crc(path: &Path) {
    use fgumi_lib::bgzf_reader::{BGZF_FOOTER_SIZE, RawBgzfBlock, read_raw_blocks};
    let mut bytes = fs::read(path).expect("read bam for corruption");
    let mut cursor: &[u8] = &bytes;
    let blocks = read_raw_blocks(&mut cursor, 100_000).expect("read bgzf blocks from test bam");
    assert!(
        blocks.len() >= 2,
        "test input must span >= 2 BGZF blocks so the corrupted block isn't the header's; got {}",
        blocks.len()
    );
    let offset: usize = blocks[..blocks.len() - 1].iter().map(RawBgzfBlock::len).sum();
    let last = blocks.last().expect("checked len >= 2 above");
    let crc_off = offset + last.len() - BGZF_FOOTER_SIZE;
    bytes[crc_off] ^= 0x01;
    fs::write(path, bytes).expect("write corrupted bam");
}

/// Write enough query-grouped read pairs to span more than one BGZF block, so
/// the corrupted last block is a record block rather than the header's.
fn create_multiblock_clip_bam(path: &Path) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for i in 0..3000u32 {
        let name = format!("read{i}");
        let mut r1 = SamBuilder::new();
        r1.read_name(name.as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        let mut r2 = SamBuilder::new();
        r2.read_name(name.as_bytes())
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        writer.write_alignment_record(&header, &to_record_buf(&r1.build())).expect("R1");
        writer.write_alignment_record(&header, &to_record_buf(&r2.build())).expect("R2");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Build clip's `--threads N` argv with a clipping option set (clip requires at
/// least one), appending any extra flags (e.g. `--no-check-crc`).
fn clip_threads_args<'a>(
    input: &'a str,
    output: &'a str,
    reference: &'a str,
    extra: &[&'a str],
) -> Vec<&'a str> {
    let mut args = vec![
        "clip",
        "--input",
        input,
        "--output",
        output,
        "--ref",
        reference,
        "--clip-overlapping-reads",
        "--threads",
        "2",
    ];
    args.extend_from_slice(extra);
    args
}

/// `--no-check-crc` must let clip's `--threads N` pipeline decode accept a
/// corrupted BGZF CRC32 and complete (#800).
#[test]
fn test_clip_threads_no_check_crc_accepts_corrupted_crc() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_multiblock_clip_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = Clip::try_parse_from(clip_threads_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &["--no-check-crc"],
    ))
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip")
        .expect("--no-check-crc must accept a corrupted BGZF CRC32 under --threads");
    assert!(output_bam.exists(), "Output BAM not created");
}

/// Default (verify-on for file input) must reject the same corrupted CRC32 under
/// `--threads N`.
#[test]
fn test_clip_threads_rejects_corrupted_crc_by_default() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_multiblock_clip_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = Clip::try_parse_from(clip_threads_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &[],
    ))
    .expect("failed to parse clip args");
    let err = cmd
        .execute("fgumi clip")
        .expect_err("default (verify-on for file input) must reject a corrupted BGZF CRC32");
    assert!(
        format!("{err:#}").to_uppercase().contains("CRC32"),
        "error should mention CRC32: {err:#}"
    );
}

/// `--check-crc` must also reject the corrupted CRC32 under `--threads N`.
#[test]
fn test_clip_threads_check_crc_rejects_corrupted_crc() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());
    create_multiblock_clip_bam(&input_bam);
    corrupt_last_block_crc(&input_bam);

    let cmd = Clip::try_parse_from(clip_threads_args(
        input_bam.to_str().unwrap(),
        output_bam.to_str().unwrap(),
        ref_path.to_str().unwrap(),
        &["--check-crc"],
    ))
    .expect("failed to parse clip args");
    let err =
        cmd.execute("fgumi clip").expect_err("--check-crc must reject a corrupted BGZF CRC32");
    assert!(
        format!("{err:#}").to_uppercase().contains("CRC32"),
        "error should mention CRC32: {err:#}"
    );
}

/// CLIP3-05: clip must reject coordinate-sorted (non-query-grouped) input,
/// matching fgbio's `Bams.requireQueryGrouped`. On coordinate-sorted input mates
/// scatter, so pair clip / overlap / mate-fix silently no-op — hard-fail instead.
///
/// Exercised at both worker counts (no `--threads`, i.e. the chain at a single
/// worker, and `--threads N`) since `require_query_grouped` fires on the chain
/// regardless of worker count.
#[rstest]
#[case::single_worker(&[])]
#[case::threaded(&["--threads", "2"])]
fn test_clip_rejects_coordinate_sorted_input(#[case] extra_args: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };
    let header = create_coordinate_sorted_header("chr1", 10000);
    write_paired_bam_with_header(&input_bam, &header, vec![(r1, r2)]);

    let mut args = vec![
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--ref",
        ref_path.to_str().unwrap(),
        "--clip-overlapping-reads",
    ];
    args.extend_from_slice(extra_args);
    let cmd = Clip::try_parse_from(args).expect("failed to parse clip args");

    let err = cmd.execute("fgumi clip").expect_err("must reject coordinate-sorted input");
    assert!(
        err.to_string().contains("queryname sorted or query grouped"),
        "unexpected error message: {err}"
    );
}

/// Builds a query-grouped BAM covering the shapes `--metrics` must count correctly:
/// a lone fragment, a plain overlapping pair (no pre-existing clipping), a pair with
/// pre-existing SOFT clipping on R1, and a pair with pre-existing HARD clipping on R2.
/// Every `ClipCounts` field the chain `--metrics` path can populate (`prior`,
/// `five_prime`/`three_prime` fixed clipping, `overlapping`) ends up non-zero when run
/// with `--clip-overlapping-reads --read-one-five-prime 1 --read-two-three-prime 1`.
///
/// The four-shape block is repeated `repeats` times (unique QNAME per repeat, same
/// positions -- clip only ever compares records *within* one template, so distinct
/// templates sharing coordinates is harmless) so callers can force the chain
/// `--threads N` path across enough `GroupBam` batches to exercise real cross-worker
/// `PerThreadAccumulator` sharding -- see
/// `test_clip_metrics_match_across_threading_modes` for why that matters.
fn build_clip_metrics_parity_bam(path: &Path, repeats: usize) {
    let mut records: Vec<RawRecord> = Vec::new();
    push_clip_metrics_parity_records(&mut records, repeats);
    create_bam_from_records(path, &records);
}

/// Appends `repeats` copies of the four-shape block to `records`. Factored out of
/// [`build_clip_metrics_parity_bam`] so `test_clip_metrics_not_written_on_chain_failure`
/// can append a malformed template after a run of otherwise-valid ones, in the same
/// query-grouped record stream.
fn push_clip_metrics_parity_records(records: &mut Vec<RawRecord>, repeats: usize) {
    for i in 0..repeats {
        // A lone fragment (unpaired primary), no pre-existing clipping.
        {
            let mut b = SamBuilder::new();
            b.read_name(format!("frag1_{i:05}").as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(0)
                .ref_id(0)
                .pos(499)
                .mapq(60)
                .cigar_ops(&[8 << 4]); // 8M
            records.push(b.build());
        }

        // A proper, overlapping pair with no pre-existing clipping.
        {
            let mut r1 = SamBuilder::new();
            r1.read_name(format!("pair1_{i:05}").as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(103)
                .template_length(12);
            records.push(r1.build());

            let mut r2 = SamBuilder::new();
            r2.read_name(format!("pair1_{i:05}").as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(0)
                .pos(103)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(99)
                .template_length(-12);
            records.push(r2.build());
        }

        // A pair with pre-existing SOFT clipping on R1 (2S6M): non-zero `prior`
        // base-clip count for read_one.
        {
            let mut r1 = SamBuilder::new();
            r1.read_name(format!("pair2_{i:05}").as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
                .ref_id(0)
                .pos(199)
                .mapq(60)
                .cigar_ops(&[(2 << 4) | 4, 6 << 4]) // 2S6M
                .mate_ref_id(0)
                .mate_pos(203)
                .template_length(12);
            records.push(r1.build());

            let mut r2 = SamBuilder::new();
            r2.read_name(format!("pair2_{i:05}").as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(0)
                .pos(203)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(199)
                .template_length(-12);
            records.push(r2.build());
        }

        // A pair with pre-existing HARD clipping on R2 (6M2H): non-zero `prior`
        // base-clip count for read_two. Hard-clipped bases are absent from SEQ/QUAL,
        // so R2's sequence is only 6 bases long.
        {
            let mut r1 = SamBuilder::new();
            r1.read_name(format!("pair3_{i:05}").as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
                .ref_id(0)
                .pos(299)
                .mapq(60)
                .cigar_ops(&[8 << 4]) // 8M
                .mate_ref_id(0)
                .mate_pos(303)
                .template_length(10);
            records.push(r1.build());

            let mut r2 = SamBuilder::new();
            r2.read_name(format!("pair3_{i:05}").as_bytes())
                .sequence(b"ACGTAC")
                .qualities(&[30; 6])
                .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                .ref_id(0)
                .pos(303)
                .mapq(60)
                .cigar_ops(&[6 << 4, (2 << 4) | 5]) // 6M2H
                .mate_ref_id(0)
                .mate_pos(299)
                .template_length(-10);
            records.push(r2.build());
        }
    }
}

/// Number of times the 4-shape block (`build_clip_metrics_parity_bam`) is repeated
/// in `test_clip_metrics_match_across_threading_modes`'s fixture: `4 * REPEATS`
/// templates total.
///
/// This must be large enough to force the chain `--threads N` path across *many*
/// `GroupBam` batches, not just one -- `GroupBam` caps every emitted batch at
/// `BamPipelineTuning::template_batch_size` (fixed at 500 templates regardless of
/// `--threads`; see `crates/.../pipeline/steps/tuning.rs`), and each `ClipTemplates`
/// worker call (one call per batch) runs `with_slot` on whichever `PerThreadAccumulator`
/// slot its calling OS thread owns. A too-small fixture (one batch) would let every
/// batch land on a single worker/slot even at `--threads 4`, so a
/// `ClipMetricsFinalizeHook` bug that reduced only a subset of slots (e.g. summed only
/// slot 0) would still merge correctly by accident and this test would not catch it --
/// exactly the gap `parallel_threads_populate_and_merge_distinct_clip_metrics_slots`
/// (`src/lib/per_thread_accumulator.rs`) closes deterministically at the
/// `PerThreadAccumulator` level. Here, `4 * REPEATS = 8000` templates make `16`
/// `GroupBam` batches (`8000 / 500`) available to the pipeline's pull-based `Parallel`
/// worker scheduler (every idle worker thread can claim the next available batch --
/// see `crates/fgumi-pipeline-core/src/runtime/pool.rs`'s "each worker has its own
/// clone" note) -- 4x the `--threads 4` worker count, so multiple distinct worker
/// threads claiming at least one batch each (and therefore populating multiple
/// distinct accumulator slots) is the expected, reliably-reproduced outcome, not a
/// lucky scheduling accident.
const REPEATS: usize = 2000;

/// The chain `--metrics` output must be independent of the worker count: a
/// no-`--threads` run (the chain at a single worker) and `--threads N` runs must
/// produce a byte-identical metrics TSV.
///
/// Detailed metrics are collected via a per-thread `ClippingMetricsCollection`
/// accumulator (one slot per worker), reduced by summing raw `fragment`/`read_one`/
/// `read_two` counters across slots and then calling `finalize()` once -- pure integer
/// addition is commutative and associative, so the reduction is order-independent and
/// the TSV is byte-identical regardless of thread count or how work was sharded. This
/// is the load-bearing worker-count-independence test for `--metrics`. (Byte-parity
/// against the pre-removal single-threaded binary lives in
/// `test_clip_cutover_parity.rs`.)
///
/// The fixture is sized (`REPEATS`, see its doc comment) to force real cross-worker
/// `PerThreadAccumulator` sharding at `--threads 2`/`4`, not just a single-slot
/// accumulation -- so a reduction bug that dropped a worker's slot would change the
/// merged totals and fail the byte-identity assertion below, rather than passing by
/// accident because everything happened to land in one slot.
#[rstest]
#[case::threads_1(1)]
#[case::threads_2(2)]
#[case::threads_4(4)]
fn test_clip_metrics_match_across_threading_modes(#[case] threads: usize) {
    use fgumi_lib::metrics::{ClippingMetrics, ReadType};

    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    build_clip_metrics_parity_bam(&input_bam, REPEATS);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let chain_out = temp_dir.path().join("chain.bam");
    let oracle_metrics = temp_dir.path().join("oracle.metrics.txt");
    let chain_metrics = temp_dir.path().join("chain.metrics.txt");

    let opts =
        ["--clip-overlapping-reads", "--read-one-five-prime", "1", "--read-two-three-prime", "1"];

    // Reference: no --threads (the chain at a single worker).
    {
        let mut args = vec![
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            oracle_out.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--metrics",
            oracle_metrics.to_str().unwrap(),
        ];
        args.extend_from_slice(&opts);
        Clip::try_parse_from(args)
            .expect("parse single-worker args")
            .execute("fgumi clip")
            .expect("single-worker clip failed");
    }

    // Chain: --threads N with --metrics. Must succeed and write the metrics file.
    {
        let threads_arg = threads.to_string();
        let mut args = vec![
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            chain_out.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--metrics",
            chain_metrics.to_str().unwrap(),
            "--threads",
            &threads_arg,
        ];
        args.extend_from_slice(&opts);
        Clip::try_parse_from(args)
            .expect("parse chain args")
            .execute("fgumi clip")
            .expect("chain clip with --metrics must now succeed under --threads");
    }
    assert!(chain_metrics.exists(), "chain --metrics file must be written");

    // Non-vacuous: every ClipCounts field the fixture exercises (prior soft/hard
    // clip, fixed 5'/3' clipping, overlap clipping) must show up as a real,
    // non-zero count in the oracle -- not an all-zero table.
    let oracle_rows =
        fgumi_lib::metrics::read_metrics::<_, ClippingMetrics>(&oracle_metrics, "clipping")
            .expect("parse oracle metrics");
    let by_type = |t: ReadType| {
        oracle_rows
            .iter()
            .find(|m| m.read_type == t)
            .unwrap_or_else(|| panic!("missing {t:?} row in metrics TSV"))
    };
    assert_eq!(by_type(ReadType::Fragment).reads, REPEATS, "one fragment read per repeat");
    assert!(
        by_type(ReadType::ReadOne).bases_clipped_pre > 0,
        "pair2's pre-existing R1 soft clip must count as `prior`"
    );
    assert!(
        by_type(ReadType::ReadTwo).bases_clipped_pre > 0,
        "pair3's pre-existing R2 hard clip must count as `prior`"
    );
    assert!(
        by_type(ReadType::ReadOne).bases_clipped_five_prime > 0
            || by_type(ReadType::ReadTwo).bases_clipped_five_prime > 0,
        "--read-one-five-prime 1 must clip at least one read"
    );
    assert!(
        by_type(ReadType::ReadTwo).bases_clipped_three_prime > 0
            || by_type(ReadType::ReadOne).bases_clipped_three_prime > 0,
        "--read-two-three-prime 1 must clip at least one read"
    );
    assert!(
        by_type(ReadType::Pair).bases_clipped_overlapping > 0,
        "pair1's overlap must be clipped and rolled up into Pair"
    );

    // The chain (--threads N) metrics TSV must be byte-identical to the oracle's.
    let oracle_content = fs::read_to_string(&oracle_metrics).expect("read oracle metrics");
    let chain_content = fs::read_to_string(&chain_metrics).expect("read chain metrics");
    assert_eq!(
        oracle_content, chain_content,
        "metrics TSV must be byte-identical between the oracle and chain (--threads {threads}) paths",
    );
}

/// Builds a query-grouped BAM of `valid_repeats` copies of the four-shape parity
/// block, followed by one malformed template: two primary (non-secondary,
/// non-supplementary) first-of-pair reads sharing one QNAME. The chain's `GroupBam`
/// step groups records into `Template`s via the shared grouper (`src/lib/template.rs`),
/// which rejects that shape ("Multiple non-secondary, non-supplemental R1 records") --
/// the same invariant `find_primary_pair_indices` in `commands::clip` also enforces,
/// defensively, on the `RawRecord` path. Placed last, so grouping fails only after every
/// valid template ahead of it has already been grouped into complete batches and handed
/// to the downstream `ClipTemplates` workers -- `clip_templates_in_batch` commits each
/// template's counts into its calling worker's live `PerThreadAccumulator` slot as it
/// loops, not just at batch end, so by the time `GroupBam` hits the malformed record and
/// fails the run, real (non-zero) metrics have very likely already been written into at
/// least one slot.
fn build_clip_metrics_failure_bam(path: &Path, valid_repeats: usize) {
    let mut records: Vec<RawRecord> = Vec::new();
    push_clip_metrics_parity_records(&mut records, valid_repeats);

    for i in 0..2 {
        let mut b = SamBuilder::new();
        b.read_name(b"malformed1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(900 + i * 20)
            .mapq(60)
            .cigar_ops(&[8 << 4]); // 8M
        records.push(b.build());
    }

    create_bam_from_records(path, &records);
}

/// A failed/partial chain `--threads N` run must leave NO `--metrics` file on disk at
/// all, not a partial/stale one.
///
/// `ClipMetricsFinalizeHook` is registered on `finalize_on_success`
/// (`src/lib/pipeline/chains/builder.rs::add_clip`), which `drain_finalize`
/// (`src/lib/pipeline/chains/finalize.rs`) runs only once `Pipeline::run` AND every
/// always-run `finalize` hook have succeeded -- see
/// `drain_finalize_skips_success_hooks_when_pipeline_fails` for the generic mechanism
/// this test exercises end-to-end through the real `clip` command. The fixture
/// (`build_clip_metrics_failure_bam`) forces a genuine mid-run pipeline failure via a
/// malformed template appended after many valid ones, so real metrics accumulation
/// happens in at least one `PerThreadAccumulator` slot before the failure -- proving the
/// hook actually discards in-progress state rather than merely never having any to write.
#[test]
fn test_clip_metrics_not_written_on_chain_failure() {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");
    let ref_path = create_test_reference(temp_dir.path());
    // REPEATS (many `GroupBam` batches' worth) valid templates ahead of the malformed
    // one, so grouping fails only after real work has already reached `ClipTemplates`.
    build_clip_metrics_failure_bam(&input_bam, REPEATS);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--ref",
        ref_path.to_str().unwrap(),
        "--clip-overlapping-reads",
        "--read-one-five-prime",
        "1",
        "--read-two-three-prime",
        "1",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--threads",
        "4",
    ])
    .expect("failed to parse clip args");

    let err = cmd.execute("fgumi clip").expect_err("malformed template must fail the chain run");
    assert!(
        err.to_string().contains("Multiple non-secondary, non-supplemental R1"),
        "unexpected error message: {err}"
    );
    assert!(
        !metrics_path.exists(),
        "a failed chain run must not write a --metrics file, even a partial one"
    );
}

/// A header carrying only `@SQ` (no `@HD` line) -- mirrors header-less input.
fn headerless_header(ref_name: &str, ref_len: usize) -> noodles::sam::Header {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;

    let header = noodles::sam::Header::builder()
        .add_reference_sequence(
            ref_name,
            Map::<ReferenceSequence>::new(
                std::num::NonZero::new(ref_len).expect("non-zero reference length"),
            ),
        )
        .build();
    assert!(header.header().is_none(), "precondition: header must lack @HD");
    header
}

/// Test basic clip command with fixed-end clipping.
#[test]
fn test_clip_command_basic() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Create a paired-end read
    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };

    create_paired_bam(&input_bam, vec![(r1, r2)]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--read-one-five-prime",
        "1",
        "--read-one-three-prime",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");
    assert!(output_bam.exists(), "Output BAM not created");

    // Verify output has the reads
    let mut reader = bam::io::Reader::new(fs::File::open(&output_bam).unwrap());
    let _header = reader.read_header().unwrap();
    let count = reader.records().count();
    assert_eq!(count, 2, "Should have both reads in output");
}

/// Build a query-grouped BAM of `count` overlapping paired templates. Each
/// template is an 8M R1 at pos 99 and an 8M reverse R2 at pos 103 (mates overlap
/// in [103,106]), so `--clip-overlapping-reads` plus fixed-end clipping both do
/// real work. Read names are `tmpl{i}`, with R1/R2 kept adjacent (query-grouped),
/// as clip requires.
fn build_clip_parity_bam(path: &Path, count: usize) {
    let pairs: Vec<(RawRecord, RawRecord)> = (0..count)
        .map(|i| {
            let name = format!("tmpl{i}");
            let r1 = {
                let mut b = SamBuilder::new();
                b.read_name(name.as_bytes())
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
                    .ref_id(0)
                    .pos(99)
                    .mapq(60)
                    .cigar_ops(&[8 << 4]) // 8M
                    .mate_ref_id(0)
                    .mate_pos(103)
                    .template_length(12);
                b.build()
            };
            let r2 = {
                let mut b = SamBuilder::new();
                b.read_name(name.as_bytes())
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8])
                    .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
                    .ref_id(0)
                    .pos(103)
                    .mapq(60)
                    .cigar_ops(&[8 << 4]) // 8M
                    .mate_ref_id(0)
                    .mate_pos(99)
                    .template_length(-12);
                b.build()
            };
            (r1, r2)
        })
        .collect();
    create_paired_bam(path, pairs);
}

/// Worker-count independence: clip always runs on the declarative chain
/// (`ChainSpec::single_stage(Stage::Clip, ...)` → `build_for(spec)?.run()`), so a
/// no-`--threads` run (the chain at a single worker) and `--threads N` runs must
/// produce the same records AND the same normalized header. `read_bam_output`
/// normalizes only the `@PG` `CL` field, which legitimately differs by `--threads`,
/// so a dropped `@SQ`/`@RG`/`@HD`/`@CO` line or any record divergence fails the
/// test. Parameterized over thread counts to exercise single- and multi-worker
/// scheduling. (Byte-parity against the pre-removal single-threaded binary lives in
/// `test_clip_cutover_parity.rs`.)
#[rstest]
#[case::threads_1(1)]
#[case::threads_2(2)]
#[case::threads_4(4)]
fn test_clip_chain_matches_single_threaded(#[case] threads: usize) {
    let temp_dir = TempDir::new().expect("temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    let ref_path = create_test_reference(temp_dir.path());
    build_clip_parity_bam(&input_bam, 8);

    let oracle_out = temp_dir.path().join("oracle.bam");
    let chain_out = temp_dir.path().join("chain.bam");

    let opts =
        ["--clip-overlapping-reads", "--read-one-five-prime", "1", "--read-two-three-prime", "1"];

    // Reference: no --threads (the chain at a single worker).
    {
        let mut args = vec![
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            oracle_out.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
        ];
        args.extend_from_slice(&opts);
        Clip::try_parse_from(args)
            .expect("parse single-worker args")
            .execute("fgumi clip")
            .expect("single-worker clip failed");
    }

    // Chain: --threads N (declarative chain builder).
    {
        let threads_arg = threads.to_string();
        let mut args = vec![
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            chain_out.to_str().unwrap(),
            "--ref",
            ref_path.to_str().unwrap(),
            "--threads",
            &threads_arg,
        ];
        args.extend_from_slice(&opts);
        Clip::try_parse_from(args)
            .expect("parse chain args")
            .execute("fgumi clip")
            .expect("chain clip failed");
    }

    let (oracle_header, oracle_records) = crate::helpers::read_bam_output(&oracle_out);
    let (chain_header, chain_records) = crate::helpers::read_bam_output(&chain_out);

    // Non-vacuous guards: all 16 records (8 templates x 2) survive, and clipping
    // actually ran -- at least one output CIGAR differs from the input's plain 8M
    // (robust to soft/hard clipping mode). Without this a two-empty-output run
    // would pass the equality check vacuously.
    assert_eq!(oracle_records.len(), 16, "oracle should keep all 16 records");
    assert!(
        oracle_records.iter().any(|r| cigar_ops(r) != vec![(CigarKind::Match, 8)]),
        "fixture must actually clip (expected an output CIGAR other than 8M)",
    );

    assert_eq!(
        oracle_records, chain_records,
        "chain (--threads {threads}) output must match the single-worker run record-for-record",
    );
    assert_eq!(
        oracle_header, chain_header,
        "chain and oracle output headers must match (with @PG CL normalized)",
    );
}

/// CLIP3-05: clip must reject header-less input. A header-less BAM synthesizes
/// `@HD VN:1.6 SO:unsorted` (via `ensure_hd_record`), which is neither queryname
/// sorted nor query grouped, so `require_query_grouped` rejects it — matching
/// fgbio's `Bams.requireQueryGrouped`. The `@HD` synthesis itself is still exercised
/// end-to-end by `correct`/`review` (which have no query-grouped guard).
///
/// Exercised at both worker counts (no `--threads`, i.e. the chain at a single
/// worker, and `--threads N`): the chain runs `ensure_hd_record` before
/// `require_query_grouped`, so header-less input is rejected the same way
/// regardless of worker count.
#[rstest]
#[case::single_worker(&[])]
#[case::threaded(&["--threads", "2"])]
fn test_clip_rejects_headerless_input(#[case] extra_args: &[&str]) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(103)
            .template_length(12);
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(103)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-12);
        b.build()
    };

    let header = headerless_header("chr1", 10000);
    write_paired_bam_with_header(&input_bam, &header, vec![(r1, r2)]);

    let mut args = vec![
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--read-one-five-prime",
        "1",
        "--read-one-three-prime",
        "1",
        "--compression-level",
        "1",
    ];
    args.extend_from_slice(extra_args);
    let cmd = Clip::try_parse_from(args).expect("failed to parse clip args");

    let err = cmd.execute("fgumi clip").expect_err("must reject header-less input");
    assert!(
        err.to_string().contains("queryname sorted or query grouped"),
        "unexpected error message: {err}"
    );
}

/// Test clip command with metrics output.
#[test]
fn test_clip_command_with_metrics() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let metrics_path = temp_dir.path().join("metrics.txt");
    let ref_path = create_test_reference(temp_dir.path());

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(108);
        b.build()
    };

    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE)
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .cigar_ops(&[8 << 4]) // 8M
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-108);
        b.build()
    };

    create_paired_bam(&input_bam, vec![(r1, r2)]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--read-one-five-prime",
        "2",
        "--metrics",
        metrics_path.to_str().unwrap(),
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command with metrics failed");
    assert!(metrics_path.exists(), "Metrics file not created");
}

/// Write a BAM from an explicit list of records (any flags), using the shared minimal header.
fn create_bam_from_records(path: &Path, records: &[RawRecord]) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for r in records {
        writer.write_alignment_record(&header, &to_record_buf(r)).expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// One decoded output record: `(name, cigar-ops, raw-flags)`.
type OutputRecord = (String, Vec<(CigarKind, usize)>, u16);

/// Read output records so tests can assert the *actual* clip state (CIGAR changed /
/// unchanged), not merely the record count.
fn read_output_records(path: &Path) -> Vec<OutputRecord> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader
        .record_bufs(&header)
        .map(|r| {
            let rec = r.expect("Failed to read output record");
            let name = rec
                .name()
                .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                .unwrap_or_default();
            let ops: Vec<(CigarKind, usize)> =
                rec.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect();
            (name, ops, u16::from(rec.flags()))
        })
        .collect()
}

/// Read output records as full noodles `RecordBuf`s so tests can assert mate metadata
/// (mate ref/pos/strand, MC, MQ, TLEN) — not just CIGAR — after clipping.
fn read_output_record_bufs(path: &Path) -> Vec<RecordBuf> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).map(|r| r.expect("Failed to read output record")).collect()
}

/// The record's CIGAR as `(kind, len)` ops (e.g. `[(Match, 98), (HardClip, 2)]`).
fn cigar_ops(rec: &RecordBuf) -> Vec<(CigarKind, usize)> {
    rec.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect()
}

/// A record's identity for exact match/count assertions:
/// `(name, flags, 1-based pos, 1-based mate-pos, CIGAR ops, SEQ bytes, QUAL scores)`.
///
/// SEQ and QUAL are included so a hard-clip upgrade that rewrites the CIGAR (e.g. `10S40M`
/// -> `10H40M`) but fails to drop the hard-clipped bases from the payload — leaving 50
/// SEQ/QUAL bytes under a 40-base-consuming CIGAR — is rejected, not silently accepted.
type RecordIdentity = (String, u16, usize, usize, Vec<(CigarKind, usize)>, Vec<u8>, Vec<u8>);

/// Extracts the [`RecordIdentity`] tuple from a mapped, mate-mapped record.
fn record_identity(rec: &RecordBuf) -> RecordIdentity {
    let name =
        rec.name().map(|n| String::from_utf8_lossy(n.as_ref()).into_owned()).unwrap_or_default();
    (
        name,
        u16::from(rec.flags()),
        usize::from(rec.alignment_start().expect("mapped record has alignment start")),
        usize::from(rec.mate_alignment_start().expect("record has mate alignment start")),
        cigar_ops(rec),
        rec.sequence().as_ref().to_vec(),
        rec.quality_scores().as_ref().to_vec(),
    )
}

/// The record's mate-CIGAR (`MC`) tag as a string, if present.
fn mate_cigar(rec: &RecordBuf) -> Option<String> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&Tag::MATE_CIGAR) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

/// The record's mate-mapping-quality (`MQ`) tag as an integer, if present.
fn mate_mapq(rec: &RecordBuf) -> Option<i64> {
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    rec.data().get(&Tag::MATE_MAPPING_QUALITY).and_then(Value::as_int)
}

/// Build a mapped read with the given flags/position/CIGAR-length (all `M`).
fn mapped_read(name: &[u8], flags: u16, pos: i32, match_len: u32, mapq: u8) -> RawRecord {
    let len = match_len as usize;
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(&vec![b'A'; len])
        .qualities(&vec![30; len])
        .flags(flags)
        .ref_id(0)
        .pos(pos)
        .mapq(mapq)
        .cigar_ops(&[match_len << 4]) // <len>M
        .mate_ref_id(0)
        .mate_pos(pos);
    b.build()
}

/// Build a mapped read with an explicit CIGAR (raw BAM `(len << 4) | op` codes) and a
/// matching-length sequence/qualities. `mate_pos` sets the mate's position (PNEXT); pass the
/// mate read's start, not the read's own. `read_len` must equal the CIGAR's read-consuming length.
fn read_with_cigar(
    name: &[u8],
    flags: u16,
    pos: i32,
    mate_pos: i32,
    cigar_ops: &[u32],
    read_len: usize,
    mapq: u8,
) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(&vec![b'A'; read_len])
        .qualities(&vec![30; read_len])
        .flags(flags)
        .ref_id(0)
        .pos(pos)
        .mapq(mapq)
        .cigar_ops(cigar_ops)
        .mate_ref_id(0)
        .mate_pos(mate_pos);
    b.build()
}

/// Like [`read_with_cigar`] but also sets the `MC` (mate CIGAR) tag and derives the read length
/// from `cigar_ops` (rather than taking it as a separate, easily-mismatched argument). The
/// past-mate-clip geometries below don't strictly require `MC` (`clip_extending_past_mate_ends`
/// reads the mate's own CIGAR from the mate record in hand -- see
/// `fgumi_raw_bam::num_bases_extending_past_mate_vs_mate_raw`), but real query-grouped input
/// from an aligner carries it, so wiring it here keeps these fixtures representative.
fn read_with_cigar_mc(
    name: &[u8],
    flags: u16,
    pos: i32,
    mate_pos: i32,
    cigar_ops: &[u32],
    mapq: u8,
    mate_cigar: &[u8],
) -> RawRecord {
    let read_len = fgumi_raw_bam::query_length_from_cigar(cigar_ops);
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(&vec![b'A'; read_len])
        .qualities(&vec![30; read_len])
        .flags(flags)
        .ref_id(0)
        .pos(pos)
        .mapq(mapq)
        .cigar_ops(cigar_ops)
        .mate_ref_id(0)
        .mate_pos(mate_pos)
        .add_string_tag(SamTag::MC, mate_cigar);
    b.build()
}

/// Runs `clip --upgrade-clipping` (Hard mode) over a template whose only clipped read is a
/// supplementary alignment (`10S40M`), and asserts the supplementary's leading soft clip is
/// upgraded to hard. fgbio `ClipBam` upgrades `template.allReads` (`ClipBam.scala:123`), not just
/// the primary R1/R2, so the supplementary must be upgraded too. `threads` selects a
/// single worker (`None`) vs `--threads` — the fix must hold at both worker counts.
fn run_supplementary_upgrade_case(threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Primary pair, both 50M (no clipping to upgrade — they must stay 50M).
    let r1 = mapped_read(b"t", flags::PAIRED | flags::FIRST_SEGMENT, 99, 50, 60);
    let r2 = mapped_read(b"t", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 99, 50, 60);
    // R1 supplementary with a leading soft clip: 10S40M. The only read carrying clipping.
    let supp = read_with_cigar(
        b"t",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY,
        199,
        99, // mate = primary R2 at position 99, not the supplementary's own position
        &[(10u32 << 4) | 4, 40u32 << 4], // 10S40M
        50,
        60,
    );
    create_bam_from_records(&input_bam, &[r1, r2, supp]);

    let mut args = vec![
        "clip".to_string(),
        "--input".to_string(),
        input_bam.to_str().unwrap().to_string(),
        "--output".to_string(),
        output_bam.to_str().unwrap().to_string(),
        "--reference".to_string(),
        ref_path.to_str().unwrap().to_string(),
        "--clipping-mode".to_string(),
        "hard".to_string(),
        "--upgrade-clipping".to_string(),
        "true".to_string(),
        "--compression-level".to_string(),
        "1".to_string(),
    ];
    if let Some(t) = threads {
        args.push("--threads".to_string());
        args.push(t.to_string());
    }
    let cmd = Clip::try_parse_from(&args).expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    let recs = read_output_record_bufs(&output_bam);
    assert_eq!(recs.len(), 3, "all three reads retained");

    // Assert the exact identity of all three expected records — (name, flags, 1-based pos,
    // 1-based mate_pos, CIGAR) — each present exactly once. Count-only checks would pass even
    // with two copies of the same primary, so this pins down every record.
    //
    // Values are hand-derived from the fixed inputs and the post-clip mate-info fixing:
    //   * Only `--upgrade-clipping` (Hard) is on — no fixed/overlap/extension clipping — so both
    //     primaries stay 50M and the supplementary's leading 10S is upgraded to 10H.
    //   * `set_mate_info_raw` sets MATE_REVERSE on primary R1 (its mate R2 is reverse) and, via
    //     `fix_supplemental_mate_info`, on the supplementary (its mate primary R2 is reverse); it
    //     also points every read's mate at primary position 99 (0-based) => 100 (1-based).
    // Fixtures build reads with SEQ = all `A` and QUAL = all 30 (`read_with_cigar` /
    // `mapped_read`), so the unchanged 50M primaries keep 50 bytes each and the hard-clipped
    // supplementary must drop its leading 10 to 40 bytes.
    let expected: [RecordIdentity; 3] = [
        // Primary R1: forward, unchanged 50M; MATE_REVERSE now set (mate R2 is reverse).
        (
            "t".to_string(),
            flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
            100,
            100,
            vec![(CigarKind::Match, 50)],
            vec![b'A'; 50],
            vec![30; 50],
        ),
        // Primary R2: reverse, unchanged 50M.
        (
            "t".to_string(),
            flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE,
            100,
            100,
            vec![(CigarKind::Match, 50)],
            vec![b'A'; 50],
            vec![30; 50],
        ),
        // Supplementary R1: leading 10S upgraded to 10H (fgbio ClipBam.scala:123 upgrades
        // template.allReads); MATE_REVERSE set from the primary R2 mate. The hard clip drops the
        // leading 10 bases, so SEQ/QUAL shrink 50 -> 40.
        (
            "t".to_string(),
            flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY | flags::MATE_REVERSE,
            200,
            100,
            vec![(CigarKind::HardClip, 10), (CigarKind::Match, 40)],
            vec![b'A'; 40],
            vec![30; 40],
        ),
    ];
    let actual: Vec<RecordIdentity> = recs.iter().map(record_identity).collect();
    for want in &expected {
        let count = actual.iter().filter(|got| *got == want).count();
        assert_eq!(
            count, 1,
            "expected exactly one record matching {want:?}; got {count} in {actual:?}"
        );
    }
}

/// `--upgrade-clipping` upgrades a supplementary alignment's clipping at both a single
/// worker (`None`) and `--threads` (`Some("2")`).
#[rstest]
#[case::single_worker(None)]
#[case::threaded(Some("2"))]
fn test_upgrade_clipping_upgrades_supplementary(#[case] threads: Option<&str>) {
    run_supplementary_upgrade_case(threads);
}

/// `--threads` mode routes through the declarative chain builder (`execute_chain` → `add_clip` →
/// `build_clip_process_step`, which delegates to `ClipParams::clip_template`); exercise the
/// `(Some, Some)` primary-pair
/// branch with fixed R1/R2 clipping, overlap clipping, mate-extension clipping and clip upgrading
/// all enabled. R1 (fwd, 150M) and R2 (rev, 100M) fully overlap, so overlap clipping engages.
#[test]
fn test_clip_command_threads_mode_primary_pair_all_options() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let mut r1 =
        mapped_read(b"p", flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE, 99, 150, 60);
    let mut r2 =
        mapped_read(b"p", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 99, 100, 60);
    // Set a realistic FR insert size on the pair. Overlap and mate-extension clipping gate on the
    // symmetric `is_primary_fr_pair_raw` (the reverse read's CIGAR-derived arm), so FR detection is
    // TLEN-independent; the values are set only because a real FR pair carries them.
    fgumi_raw_bam::set_template_length(r1.as_mut_vec(), 150);
    fgumi_raw_bam::set_template_length(r2.as_mut_vec(), -150);
    create_bam_from_records(&input_bam, &[r1, r2]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-one-five-prime",
        "1",
        "--read-one-three-prime",
        "1",
        "--read-two-five-prime",
        "1",
        "--read-two-three-prime",
        "1",
        "--clip-overlapping-reads",
        "true",
        "--clip-bases-past-mate",
        "true",
        "--upgrade-clipping",
        "true",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    // Exact post-clip CIGARs, derived independently from the clip pipeline (Hard mode),
    // not just "some hard clip exists". Inputs are fixed: R1 fwd 150M@100 (1-based),
    // R2 rev 100M@100.
    //
    //   R1: fixed 5'/3' = 1 each  -> 1H148M1H, ref 101..248
    //       overlap clip trims R1's 3' end down to the pair midpoint
    //       (midpoint(r1_start=101, r2_end=198) = 149), clipping ref 150..248 (99 M bases)
    //       -> 1H49M100H, ref 101..149
    //   R2: fixed 5'/3' = 1 each  -> 1H98M1H, ref 101..198
    //       overlap clip trims R2's 5' (low-coord) end up to midpoint+1 = 150,
    //       clipping ref 101..149 (49 M bases) -> 50H49M1H, ref 150..198
    //
    // Mate-extension clipping then finds nothing past the mate (the reads now meet at the
    // midpoint), so it is a no-op here — but its code path still runs. Totals check out:
    // R1 = 1+49+100 = 150 read bases, R2 = 50+49+1 = 100 read bases.
    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 2, "both reads retained");
    let expected_r1 =
        vec![(CigarKind::HardClip, 1), (CigarKind::Match, 49), (CigarKind::HardClip, 100)];
    let expected_r2 =
        vec![(CigarKind::HardClip, 50), (CigarKind::Match, 49), (CigarKind::HardClip, 1)];
    for (name, ops, flag_bits) in &recs {
        assert_eq!(name, "p", "unexpected read name");
        let is_r1 = flag_bits & flags::FIRST_SEGMENT != 0;
        let (label, expected) = if is_r1 { ("R1", &expected_r1) } else { ("R2", &expected_r2) };
        assert_eq!(ops, expected, "{label} exact post-clip CIGAR mismatch; cigar={ops:?}");
    }
}

/// `--threads` mode with a lone unpaired fragment exercises the `(Some, None)` fragment branch.
#[test]
fn test_clip_command_threads_mode_fragment() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // No PAIRED flag: find_primary_pair_indices resolves R1 only.
    let frag = mapped_read(b"f", 0, 99, 100, 60);
    create_bam_from_records(&input_bam, &[frag]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-one-five-prime",
        "2",
        "--read-one-three-prime",
        "2",
        "--upgrade-clipping",
        "true",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    // The lone forward-strand fragment is hard-clipped 2bp at each end (read-one 5'/3' = 2,
    // Hard mode) => 2H96M2H, proving the (Some, None) branch actually clipped rather than
    // just passing the record through.
    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 1, "fragment retained");
    assert_eq!(
        recs[0].1,
        vec![(CigarKind::HardClip, 2), (CigarKind::Match, 96), (CigarKind::HardClip, 2)],
        "fragment must be hard-clipped 2bp at each end (2H96M2H); cigar={:?}",
        recs[0].1
    );
}

/// `--threads` mode with a template that has no primary read (secondary only) exercises the
/// `(None, None)` no-op branch — the record passes through untouched.
#[test]
fn test_clip_command_threads_mode_secondary_only() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    let secondary =
        mapped_read(b"s", flags::PAIRED | flags::FIRST_SEGMENT | flags::SECONDARY, 99, 100, 60);
    create_bam_from_records(&input_bam, &[secondary]);

    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-one-five-prime",
        "1",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    // The (None, None) branch is a true no-op: the secondary read passes through with its
    // CIGAR unchanged (still 100M, not clipped) and its SECONDARY flag preserved.
    let recs = read_output_records(&output_bam);
    assert_eq!(recs.len(), 1, "secondary read passed through");
    let (_name, ops, flag_bits) = &recs[0];
    assert_eq!(
        ops,
        &vec![(CigarKind::Match, 100)],
        "secondary read must be unchanged (100M no-op); cigar={ops:?}"
    );
    assert_ne!(flag_bits & flags::SECONDARY, 0, "SECONDARY flag must be preserved");
}

/// `--threads` mode with a chimeric/split template — primary R1 + primary R2 + a supplementary
/// R1 — exercises the PR's headline fix end-to-end: after the primary pair is clipped,
/// `fix_supplemental_mate_info` must repair the supplementary alignment's mate metadata to point
/// at the *post-clip* primary R2. The `clip.rs` unit test covers the helper in isolation; this
/// verifies the `--threads` chain-builder wiring (`execute_chain` → `build_clip_process_step`:
/// correct r1/r2 indices, and that the supplemental
/// mate snapshot is taken *after* the primary is clipped, not before). A wiring bug (wrong
/// index, or records mutated out of order) would be invisible to every other test in this file.
#[test]
fn test_clip_command_threads_mode_supplementary_mate_repair() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // Chimeric template: a primary FR pair (R1 fwd 100M @100, R2 rev 100M @300) plus a
    // supplementary R1 mapped far away (50M @5000). `mapped_read` seeds each record's mate fields
    // to its *own* position with no MC/MQ, so the supplementary's input mate info is deliberately
    // stale (points at 5000, not the primary R2) — only a correct repair yields the expected
    // primary-R2-derived values.
    let r1 = mapped_read(
        b"chim",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
        99,
        100,
        60,
    );
    let r2 =
        mapped_read(b"chim", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 299, 100, 40);
    let supp = mapped_read(
        b"chim",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY,
        4999,
        50,
        30,
    );
    create_bam_from_records(&input_bam, &[r1, r2, supp]);

    // Clip only read-two's 5' end by 2bp. R2 is reverse, so its 5' end is the high-coordinate
    // (right) end: 100M @300 -> 98M2H, alignment start unchanged at 300 (1-based). The
    // supplementary R1 is never touched by the clip loop (only the primary pair is), so it keeps
    // 50M — but its mate fields must be rewritten to the post-clip primary R2.
    let cmd = Clip::try_parse_from([
        "clip",
        "--input",
        input_bam.to_str().unwrap(),
        "--output",
        output_bam.to_str().unwrap(),
        "--reference",
        ref_path.to_str().unwrap(),
        "--threads",
        "2",
        "--read-two-five-prime",
        "2",
        "--compression-level",
        "1",
    ])
    .expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    let recs = read_output_record_bufs(&output_bam);
    assert_eq!(recs.len(), 3, "all three records retained");

    let flag_bits = |r: &RecordBuf| u16::from(r.flags());
    let supp = recs
        .iter()
        .find(|r| flag_bits(r) & flags::SUPPLEMENTARY != 0)
        .expect("supplementary record present");
    let r2 = recs
        .iter()
        .find(|r| {
            let f = flag_bits(r);
            f & flags::LAST_SEGMENT != 0
                && f & flags::SUPPLEMENTARY == 0
                && f & flags::SECONDARY == 0
        })
        .expect("primary R2 present");

    // Independent oracle for the clip itself: primary R2 (reverse, 5' clip of 2) -> 98M2H @300,
    // start unchanged, MAPQ 40. Derived by hand from the fixed inputs, not from the code under test.
    assert_eq!(
        cigar_ops(r2),
        vec![(CigarKind::Match, 98), (CigarKind::HardClip, 2)],
        "primary R2 must be 98M2H after a 5' clip of 2"
    );
    assert_eq!(usize::from(r2.alignment_start().unwrap()), 300, "R2 alignment start unchanged");
    assert_eq!(r2.mapping_quality().map(u8::from), Some(40), "R2 MAPQ unchanged");
    assert_ne!(flag_bits(r2) & flags::REVERSE, 0, "primary R2 is reverse");

    // The supplementary alignment itself is not clipped (only the primary pair is).
    assert_eq!(cigar_ops(supp), vec![(CigarKind::Match, 50)], "supplementary unchanged (50M)");

    // Headline fix: the supplementary's mate metadata now points at the post-clip primary R2.
    // Absolute values are hand-derived from the inputs (mate ref 0, mate pos 300, MC 98M2H, MQ 40,
    // mate-reverse set) — so this is an independent oracle, not a self-consistency check.
    assert_eq!(supp.mate_reference_sequence_id(), Some(0), "mate ref = primary R2 ref");
    assert_eq!(
        usize::from(supp.mate_alignment_start().unwrap()),
        300,
        "mate pos = primary R2 start"
    );
    assert_ne!(
        flag_bits(supp) & flags::MATE_REVERSE,
        0,
        "mate-reverse must be set (primary R2 is reverse)"
    );
    assert_eq!(mate_cigar(supp).as_deref(), Some("98M2H"), "MC = primary R2 post-clip CIGAR");
    assert_eq!(mate_mapq(supp), Some(40), "MQ = primary R2 MAPQ");
    // TLEN is computed from the supplementary's own alignment against the post-clip primary R2,
    // not copied from R2's TLEN (see issue #673). Hand-derived from the inputs:
    //   R1 primary  100..199 forward      -> 5' = 100
    //   R2 primary  300..397 reverse      -> 5' = 397   => primaries are +298 / -298
    //   R1 supp    5000..5049 forward     -> 5' = 5000
    // The supplementary is the rightmost segment, so its TLEN is negative: 397 - 5000 - 1.
    assert_eq!(r2.template_length(), -298, "primary R2 TLEN");
    assert_eq!(
        supp.template_length(),
        -4604,
        "supplementary TLEN is computed against primary R2, not copied from it"
    );

    // Cross-check: the repaired mate fields agree with the actual primary R2 record, proving the
    // snapshot was taken from the correct (post-clip) primary and not a stale/wrong index.
    assert_eq!(supp.mate_reference_sequence_id(), r2.reference_sequence_id());
    assert_eq!(supp.mate_alignment_start(), r2.alignment_start());
}

/// Clip always runs on the declarative chain (via `build_clip_process_step` →
/// `ClipParams::clip_template`), so its output must be independent of the worker count: a
/// no-`--threads` run (the chain at a single worker) and a `--threads N` run must yield
/// byte-identical output. This pins that invariant end-to-end: if the two configurations
/// ever diverge (e.g. a future edit to the batching or reduction), this comparison fails.
#[test]
fn test_clip_command_single_and_multi_threaded_outputs_match() {
    let temp_dir = TempDir::new().unwrap();
    let ref_path = create_test_reference(temp_dir.path());
    let input_bam = temp_dir.path().join("input.bam");

    // A fully-overlapping FR pair (exercises fixed + overlap + mate-extension clipping and mate
    // repair) plus a lone fragment (exercises the fragment branch). Non-zero TLEN is required for
    // FR detection to engage overlap/mate-extension clipping.
    let mut r1 = mapped_read(
        b"pair",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
        99,
        150,
        60,
    );
    let mut r2 =
        mapped_read(b"pair", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 99, 100, 60);
    fgumi_raw_bam::set_template_length(r1.as_mut_vec(), 150);
    fgumi_raw_bam::set_template_length(r2.as_mut_vec(), -150);
    let frag = mapped_read(b"frag", 0, 149, 80, 60);
    create_bam_from_records(&input_bam, &[r1, r2, frag]);

    let run_clip = |output: &Path, threads: Option<&str>| {
        let mut args = vec![
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--read-one-five-prime",
            "1",
            "--read-one-three-prime",
            "1",
            "--read-two-five-prime",
            "1",
            "--read-two-three-prime",
            "1",
            "--clip-overlapping-reads",
            "true",
            "--clip-bases-past-mate",
            "true",
            "--compression-level",
            "1",
        ];
        if let Some(t) = threads {
            args.push("--threads");
            args.push(t);
        }
        Clip::try_parse_from(args)
            .expect("failed to parse clip args")
            .execute("fgumi clip")
            .unwrap();
    };

    let single_out = temp_dir.path().join("single.bam");
    let multi_out = temp_dir.path().join("multi.bam");
    run_clip(&single_out, None);
    run_clip(&multi_out, Some("2"));

    // Output order is the input order in both modes, so the full record sequences must be equal.
    let single_records = read_output_record_bufs(&single_out);
    let multi_records = read_output_record_bufs(&multi_out);
    assert_eq!(
        single_records, multi_records,
        "single-worker and multi-worker clip output must be identical"
    );

    // Independent oracle: the cross-mode equality above only proves the two paths agree — a
    // regression that broke (or silently no-op'd) clipping identically in both would still pass it.
    // This asserts clipping actually happened and matches what fgbio produces for this input, so the
    // test fails when clipping stops. `single_records == multi_records`, so checking one covers both.
    assert_expected_clipping_applied(&single_records);
}

/// Independent (non-self-consistency) oracle for the
/// [`test_clip_command_single_and_multi_threaded_outputs_match`] input: an overlapping FR pair
/// (R1 150M, R2 100M) plus a lone fragment (80M), clipped 1 base at each end in the default Hard
/// mode. Asserts the emitted records actually differ from the unclipped input the way fgbio would
/// clip them, so a regression that turns clipping into a no-op is caught rather than passing a
/// path-vs-path comparison.
fn assert_expected_clipping_applied(records: &[RecordBuf]) {
    let named = |name: &[u8]| -> &RecordBuf {
        records
            .iter()
            .find(|r| r.name().is_some_and(|n| AsRef::<[u8]>::as_ref(n) == name))
            .unwrap_or_else(|| panic!("missing record {}", String::from_utf8_lossy(name)))
    };

    // Default clipping mode is Hard, so the lone fragment loses exactly its 1 five-prime + 1
    // three-prime base: input 80M at 1-based 150 becomes 1H78M1H, and the 5' hard clip advances the
    // alignment start by one (150 -> 151). Hard-clipped bases are dropped from the payload (80 -> 78).
    let frag = named(b"frag");
    assert_eq!(
        cigar_ops(frag),
        vec![(CigarKind::HardClip, 1), (CigarKind::Match, 78), (CigarKind::HardClip, 1)],
        "fragment must lose 1 base at each end (80M -> 1H78M1H)"
    );
    assert_eq!(
        usize::from(frag.alignment_start().expect("fragment is mapped")),
        151,
        "fragment 5' hard clip advances the alignment start by one (150 -> 151)"
    );
    assert_eq!(
        frag.sequence().as_ref().len(),
        78,
        "hard-clipped bases must be dropped from the fragment payload"
    );

    // Both reads of the "pair" template are clipped away from their unclipped 150M / 100M inputs.
    let pair_r1 = records
        .iter()
        .find(|r| {
            r.name().is_some_and(|n| AsRef::<[u8]>::as_ref(n) == b"pair")
                && r.flags().is_first_segment()
        })
        .expect("missing pair R1");
    let pair_r2 = records
        .iter()
        .find(|r| {
            r.name().is_some_and(|n| AsRef::<[u8]>::as_ref(n) == b"pair")
                && r.flags().is_last_segment()
        })
        .expect("missing pair R2");
    assert_ne!(
        cigar_ops(pair_r1),
        vec![(CigarKind::Match, 150)],
        "pair R1 must be clipped, not passed through as 150M"
    );
    assert_ne!(
        cigar_ops(pair_r2),
        vec![(CigarKind::Match, 100)],
        "pair R2 must be clipped, not passed through as 100M"
    );
}

/// A lone primary R2 (second-of-pair with no first-of-pair mate in the template) is a malformed
/// template that fgbio `ClipBam` deliberately passes through unclipped (`case _ => ()`,
/// `ClipBam.scala:133`). `clip_template` mirrors that with its explicit `(None, _)` arm, so such a
/// read must emerge byte-for-byte unchanged under both threading modes: clipping thresholds that
/// would trim a fragment or a pair must not touch it. This pins the fgbio-parity behavior so a
/// future edit that starts clipping (or rejecting) lone R2s fails loudly.
#[test]
fn test_clip_command_lone_r2_primary_passed_through_unclipped() {
    let temp_dir = TempDir::new().unwrap();
    let ref_path = create_test_reference(temp_dir.path());
    let input_bam = temp_dir.path().join("input.bam");

    // A single paired, second-of-pair primary read with no R1 anywhere in the template.
    let r2 =
        mapped_read(b"orphan", flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE, 199, 100, 60);
    create_bam_from_records(&input_bam, &[r2]);

    let run_clip = |output: &Path, threads: Option<&str>| {
        // Aggressive thresholds that would visibly clip a fragment or pair — a lone R2 must ignore
        // them entirely.
        let mut args = vec![
            "clip",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--reference",
            ref_path.to_str().unwrap(),
            "--read-one-five-prime",
            "5",
            "--read-one-three-prime",
            "5",
            "--read-two-five-prime",
            "5",
            "--read-two-three-prime",
            "5",
            "--clip-overlapping-reads",
            "true",
            "--clip-bases-past-mate",
            "true",
            "--compression-level",
            "1",
        ];
        if let Some(t) = threads {
            args.push("--threads");
            args.push(t);
        }
        Clip::try_parse_from(args)
            .expect("failed to parse clip args")
            .execute("fgumi clip")
            .unwrap();
    };

    let single_out = temp_dir.path().join("single.bam");
    let multi_out = temp_dir.path().join("multi.bam");
    run_clip(&single_out, None);
    run_clip(&multi_out, Some("2"));

    let single_records = read_output_record_bufs(&single_out);
    let multi_records = read_output_record_bufs(&multi_out);

    // Both modes emit the single input record with no clipping applied (still 100M at 1-based 200).
    for (label, records) in [("single-worker", &single_records), ("multi-worker", &multi_records)] {
        assert_eq!(records.len(), 1, "{label}: expected exactly one output record");
        assert_eq!(
            cigar_ops(&records[0]),
            vec![(CigarKind::Match, 100)],
            "{label}: lone R2 must be passed through unclipped (still 100M)"
        );
        assert_eq!(
            usize::from(records[0].alignment_start().expect("record is mapped")),
            200,
            "{label}: lone R2 alignment start must be unchanged"
        );
    }
    assert_eq!(
        single_records, multi_records,
        "both threading modes agree on the lone-R2 passthrough"
    );
}

// ===================================================================
// #760: query-space past-mate clipping, end-to-end through the LIVE `fgumi clip` command.
//
// The unit-level fix lives in `RawRecordClipper::clip_extending_past_mate_ends`
// (`crates/fgumi-sam/src/clipper.rs`), exercised in isolation by
// `raw_clip_extending_past_mate_ends_query_space`. These tests drive the real command
// (`Clip::execute`) instead, so they prove the fix actually takes effect on the path a user
// runs (`--clip-bases-past-mate`), and pin the final CIGARs the whole pipeline produces --
// something the isolated function's unit tests, which only assert the returned clip-count
// tuple, cannot reach.
// ===================================================================

/// fgbio#1172 knock-on geometry through the LIVE `fgumi clip` command: r1 carries a deletion
/// right at the mate's un-soft-clipped end (`2S124M1D3M`@101/+), r2 is `115M14S`@97/-. The OLD
/// reference-space `clip_extending_past_mate_ends` mapped the deletion boundary to "no read
/// position" and clipped r1's *entire* read, unmapping it -- and because it unmapped r1 before
/// reaching r2, r2's mate was left unclipped too. The query-space fix (`RawRecordClipper`, Task
/// 4R) clips both reads by a small, finite amount instead. `fgumi clip`'s own default clipping
/// mode is Hard, so Soft is forced explicitly here to match the unit test
/// (`deletion_knock_on`) this test's expected CIGARs are hand-derived from.
#[rstest]
#[case::single_worker(None)]
#[case::threaded(Some("2"))]
fn test_clip_command_past_mate_knock_on_regression(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // r1 (forward): 2S124M1D3M @ 1-based 101 (0-based 100).
    let r1 = read_with_cigar_mc(
        b"knockon",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
        100,
        96,
        &[(2u32 << 4) | 4, 124u32 << 4, (1u32 << 4) | 2, 3u32 << 4],
        60,
        b"115M14S",
    );
    // r2 (reverse): 115M14S @ 1-based 97 (0-based 96).
    let r2 = read_with_cigar_mc(
        b"knockon",
        flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE,
        96,
        100,
        &[115u32 << 4, (14u32 << 4) | 4],
        60,
        b"2S124M1D3M",
    );
    create_bam_from_records(&input_bam, &[r1, r2]);

    let mut args = vec![
        "clip".to_string(),
        "--input".to_string(),
        input_bam.to_str().unwrap().to_string(),
        "--output".to_string(),
        output_bam.to_str().unwrap().to_string(),
        "--reference".to_string(),
        ref_path.to_str().unwrap().to_string(),
        "--clipping-mode".to_string(),
        "soft".to_string(),
        "--clip-bases-past-mate".to_string(),
        "true".to_string(),
        "--compression-level".to_string(),
        "1".to_string(),
    ];
    if let Some(t) = threads {
        args.push("--threads".to_string());
        args.push(t.to_string());
    }
    let cmd = Clip::try_parse_from(&args).expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    let recs = read_output_record_bufs(&output_bam);
    assert_eq!(recs.len(), 2, "both reads retained");

    let r1_out = recs.iter().find(|r| r.flags().is_first_segment()).expect("r1 present");
    let r2_out = recs.iter().find(|r| r.flags().is_last_segment()).expect("r2 present");

    // Core regression: neither read was unmapped by the past-mate clip (OLD code unmapped r1).
    assert!(!r1_out.flags().is_unmapped(), "r1 must remain mapped");
    assert!(!r2_out.flags().is_unmapped(), "r2 must remain mapped");
    assert!(r1_out.alignment_start().is_some(), "r1 must have an alignment start");
    assert!(r2_out.alignment_start().is_some(), "r2 must have an alignment start");

    // r1 clipped a small, finite amount (2 query bases at the 3' end): 2S124M1D3M -> 2S124M1D1M2S.
    let r1_cigar = cigar_ops(r1_out);
    assert_eq!(
        r1_cigar,
        vec![
            (CigarKind::SoftClip, 2),
            (CigarKind::Match, 124),
            (CigarKind::Deletion, 1),
            (CigarKind::Match, 1),
            (CigarKind::SoftClip, 2),
        ],
        "r1 CIGAR mismatch; got {r1_cigar:?}"
    );
    assert_eq!(
        usize::from(r1_out.alignment_start().unwrap()),
        101,
        "r1 alignment start unchanged (only its trailing end was clipped)"
    );

    // r2 also clipped (OLD code skipped the mate entirely after unmapping r1): 2 bases at its 3'
    // end, which for a reverse-strand read is the leading (low-coordinate) side, so
    // 115M14S -> 2S113M14S and the alignment start advances by 2 (97 -> 99, 1-based).
    let r2_cigar = cigar_ops(r2_out);
    assert_eq!(
        r2_cigar,
        vec![(CigarKind::SoftClip, 2), (CigarKind::Match, 113), (CigarKind::SoftClip, 14)],
        "r2 CIGAR mismatch; got {r2_cigar:?}"
    );
    assert_eq!(
        usize::from(r2_out.alignment_start().unwrap()),
        99,
        "r2 alignment start advances by the newly added leading clip (97 -> 99)"
    );
}

/// Machine-verifies the two Hard-mode past-mate-clip full-pipeline CIGARs that the
/// `RawRecordClipper` unit test `raw_clip_extending_past_mate_ends_query_space`
/// (`crates/fgumi-sam/src/clipper.rs`) only pins as a `(bases_r1, bases_r2)` return tuple -- the
/// final CIGAR is produced by the whole `fgumi clip` pipeline (which maps that tuple through
/// `clip_3_prime_end_of_read_raw` and, where the request lands inside an already-soft-clipped
/// run, `upgrade_clipping_raw`'s partial soft-to-hard upgrade), not by the isolated function, so
/// this is the only place those CIGARs are checked end-to-end. Both templates share one BAM /
/// one Hard-mode invocation (neither case needs any other clipping option):
///
/// - `ins`: r1 `70M10I23M47S`@100/+, r2 `50S70M30S`@100/-. Unit case `insertion_before_mate_end`
///   pins the return tuple `(3, 0)`. r1's past-mate request (50 query bases) exceeds its
///   existing 47-base trailing soft clip by exactly 3, so `clip_end_of_read_raw` shrinks the
///   alignment by those 3 bases (`23M` -> `20M`) and, in Hard mode, merges them with the
///   existing 47 soft-clipped bases into one 50-base hard-clip run: `70M10I20M50H`. r2's
///   past-mate request is also 50 query bases (Table B's `count_insertion` pins both sides at
///   50), which lands exactly on its existing 50-base leading soft clip -- not past it -- so
///   `clip_start_of_read_raw` takes the upgrade-only branch (contributing `0` newly-*aligned*
///   bases, hence the tuple's `0`) and upgrades the entire existing run to hard:
///   `50S70M30S` -> `50H70M30S`.
/// - `disjoint`: r1 `20M80S`@1000/+, r2 `80S20M`@1020/-. Unit case `disjoint_no_shared_pos` pins
///   `(0, 0)`: the 60-base past-mate request lands entirely inside each read's existing 80-base
///   soft clip, so `clip_end_of_read_raw` takes the upgrade-only branch (no alignment
///   shrinkage, hence the `0` in the tuple) and converts exactly 60 of the 80 existing
///   soft-clipped bases to hard, leaving the other 20 still soft: r1 -> `20M20S60H`, r2
///   (symmetric, on its leading/low-coordinate side) -> `60H20S20M`.
#[rstest]
#[case::single_worker(None)]
#[case::threaded(Some("2"))]
fn test_clip_command_past_mate_hard_mode_cigars(#[case] threads: Option<&str>) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let ref_path = create_test_reference(temp_dir.path());

    // "ins": r1 70M10I23M47S @ 1-based 100 (0-based 99), r2 50S70M30S @ same start.
    let ins_r1 = read_with_cigar_mc(
        b"ins",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
        99,
        99,
        &[70u32 << 4, (10u32 << 4) | 1, 23u32 << 4, (47u32 << 4) | 4],
        60,
        b"50S70M30S",
    );
    let ins_r2 = read_with_cigar_mc(
        b"ins",
        flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE,
        99,
        99,
        &[(50u32 << 4) | 4, 70u32 << 4, (30u32 << 4) | 4],
        60,
        b"70M10I23M47S",
    );

    // "disjoint": r1 20M80S @ 1-based 1000 (0-based 999), r2 80S20M @ 1-based 1020 (0-based 1019).
    let disjoint_r1 = read_with_cigar_mc(
        b"disjoint",
        flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE,
        999,
        1019,
        &[20u32 << 4, (80u32 << 4) | 4],
        60,
        b"80S20M",
    );
    let disjoint_r2 = read_with_cigar_mc(
        b"disjoint",
        flags::PAIRED | flags::LAST_SEGMENT | flags::REVERSE,
        1019,
        999,
        &[(80u32 << 4) | 4, 20u32 << 4],
        60,
        b"20M80S",
    );

    create_bam_from_records(&input_bam, &[ins_r1, ins_r2, disjoint_r1, disjoint_r2]);

    let mut args = vec![
        "clip".to_string(),
        "--input".to_string(),
        input_bam.to_str().unwrap().to_string(),
        "--output".to_string(),
        output_bam.to_str().unwrap().to_string(),
        "--reference".to_string(),
        ref_path.to_str().unwrap().to_string(),
        "--clipping-mode".to_string(),
        "hard".to_string(),
        "--clip-bases-past-mate".to_string(),
        "true".to_string(),
        "--compression-level".to_string(),
        "1".to_string(),
    ];
    if let Some(t) = threads {
        args.push("--threads".to_string());
        args.push(t.to_string());
    }
    let cmd = Clip::try_parse_from(&args).expect("failed to parse clip args");
    cmd.execute("fgumi clip").expect("Clip command failed");

    let recs = read_output_record_bufs(&output_bam);
    assert_eq!(recs.len(), 4, "all four reads retained");

    let named_segment = |name: &[u8], first: bool| -> &RecordBuf {
        recs.iter()
            .find(|r| {
                r.name().is_some_and(|n| AsRef::<[u8]>::as_ref(n) == name)
                    && r.flags().is_first_segment() == first
            })
            .unwrap_or_else(|| {
                panic!(
                    "missing {}/{}",
                    String::from_utf8_lossy(name),
                    if first { "R1" } else { "R2" }
                )
            })
    };

    let ins_r1_out = named_segment(b"ins", true);
    let ins_r2_out = named_segment(b"ins", false);
    let disjoint_r1_out = named_segment(b"disjoint", true);
    let disjoint_r2_out = named_segment(b"disjoint", false);

    let ins_r1_cigar = cigar_ops(ins_r1_out);
    assert_eq!(
        ins_r1_cigar,
        vec![
            (CigarKind::Match, 70),
            (CigarKind::Insertion, 10),
            (CigarKind::Match, 20),
            (CigarKind::HardClip, 50),
        ],
        "ins r1 CIGAR mismatch; got {ins_r1_cigar:?}"
    );
    let ins_r2_cigar = cigar_ops(ins_r2_out);
    assert_eq!(
        ins_r2_cigar,
        vec![(CigarKind::HardClip, 50), (CigarKind::Match, 70), (CigarKind::SoftClip, 30)],
        "ins r2 CIGAR mismatch; got {ins_r2_cigar:?}"
    );

    let disjoint_r1_cigar = cigar_ops(disjoint_r1_out);
    assert_eq!(
        disjoint_r1_cigar,
        vec![(CigarKind::Match, 20), (CigarKind::SoftClip, 20), (CigarKind::HardClip, 60)],
        "disjoint r1 CIGAR mismatch; got {disjoint_r1_cigar:?}"
    );
    let disjoint_r2_cigar = cigar_ops(disjoint_r2_out);
    assert_eq!(
        disjoint_r2_cigar,
        vec![(CigarKind::HardClip, 60), (CigarKind::SoftClip, 20), (CigarKind::Match, 20)],
        "disjoint r2 CIGAR mismatch; got {disjoint_r2_cigar:?}"
    );
}
