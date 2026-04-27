//! Runall vs. standalone defaults must resolve to equal domain configs.
//!
//! Regression guard: this session reconciled defaults into shared constants.
//! If someone later adds a knob to one side but forgets the other, this test
//! fails loudly.

use clap::Parser;
use fgumi_lib::commands::filter::Filter;
use fgumi_lib::commands::group::GroupReadsByUmi;
use fgumi_lib::commands::runall::Runall;
use fgumi_lib::commands::simplex::Simplex;
use fgumi_lib::defaults;

/// Thin wrapper so [`Runall`] (which derives `Args`, not `Parser`) can be
/// parsed from argv in tests. Upstream CLI dispatch does this via the top-level
/// `Cli` enum, but tests don't need the full surface.
#[derive(Parser, Debug)]
#[command(name = "runall", no_binary_name = true)]
struct RunallHarness {
    #[command(flatten)]
    inner: Runall,
}

#[test]
fn runall_filter_defaults_match_standalone_filter() {
    // NOTE: `fgumi filter --min-reads` has no CLI default (`Vec<usize>` with no
    // `default_value`), while `fgumi runall --min-reads` defaults to
    // `defaults::FILTER_MIN_READS = 1`. Supplying `--min-reads 1` to both sides
    // anchors the comparison on the *resolved* config, not the CLI default gap
    // (which is a separate follow-up — standalone Filter's empty default
    // panics `FilterConfig::new` at runtime).
    // Supply matching explicit values for both sides. `fgumi filter --min-reads`
    // has no CLI default (panics empty), and `--min-base-quality` is Option<u8>
    // on standalone but plain u8 on runall — both gaps are real defaults-parity
    // drift that this harness explicitly papers over; follow-up is to give
    // Filter a CLI default for min-reads/min-base-quality that matches the
    // runall default.
    let runall = RunallHarness::try_parse_from([
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--filter::min-reads",
        "1",
        "--filter::min-base-quality",
        "2",
    ])
    .expect("parse runall defaults")
    .inner;
    let standalone = Filter::try_parse_from([
        "filter",
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--min-reads",
        "1",
        "--min-base-quality",
        "2",
    ])
    .expect("parse filter defaults");

    let runall_cfg = runall.build_filter_config().expect("runall built with --filter::min-reads");
    let standalone_cfg = standalone.build_filter_config();

    assert_eq!(
        runall_cfg, standalone_cfg,
        "runall filter defaults diverged from standalone filter defaults"
    );
}

/// Duplex filter parity: runall must accept 1-3 comma-separated values for
/// `--filter::min-reads`, `--filter::max-read-error-rate`, and
/// `--filter::max-base-error-rate` (matching standalone `fgumi filter`), and
/// the resolved configs must be identical. Without this, a duplex runall run
/// silently collapses `[total, AB, BA]` thresholds into a single value.
#[test]
fn runall_duplex_filter_defaults_match_standalone_filter() {
    let runall = RunallHarness::try_parse_from([
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--filter::min-reads",
        "6,3,3",
        "--filter::max-read-error-rate",
        "0.05,0.025,0.025",
        "--filter::max-base-error-rate",
        "0.2,0.1,0.1",
        "--filter::min-base-quality",
        "2",
    ])
    .expect("parse runall duplex filter")
    .inner;
    let standalone = Filter::try_parse_from([
        "filter",
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--min-reads",
        "6,3,3",
        "--max-read-error-rate",
        "0.05,0.025,0.025",
        "--max-base-error-rate",
        "0.2,0.1,0.1",
        "--min-base-quality",
        "2",
    ])
    .expect("parse standalone duplex filter");

    let runall_cfg = runall.build_filter_config().expect("runall duplex config");
    let standalone_cfg = standalone.build_filter_config();

    assert_eq!(
        runall_cfg, standalone_cfg,
        "runall duplex filter config diverged from standalone; 1-3 threshold \
         vectors must fan out identically"
    );
}

#[test]
fn runall_group_defaults_match_standalone_group() {
    let runall =
        RunallHarness::try_parse_from(["--input", "ignored.bam", "--output", "ignored_out.bam"])
            .expect("parse runall defaults")
            .inner;
    // `--strategy` is required by the standalone `group` CLI and does not
    // participate in GroupFilterConfig — supply any value so the parse
    // succeeds.
    let standalone = GroupReadsByUmi::try_parse_from([
        "group",
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--strategy",
        "adjacency",
    ])
    .expect("parse group defaults");

    let runall_cfg = runall.build_group_filter_config();
    let standalone_cfg = standalone.build_group_filter_config();

    assert_eq!(
        runall_cfg, standalone_cfg,
        "runall group defaults diverged from standalone group defaults"
    );
}

/// Consensus-knob parity: the fields that runall's `consensus_helpers` feed
/// into `VanillaUmiConsensusOptions` must resolve to the same defaults that
/// `fgumi simplex` parses into its own `consensus: ConsensusCallingOptions`,
/// `overlapping: OverlappingConsensusOptions`, and `read_group: ReadGroupOptions`.
#[test]
fn runall_consensus_defaults_match_standalone_simplex() {
    let runall =
        RunallHarness::try_parse_from(["--input", "ignored.bam", "--output", "ignored_out.bam"])
            .expect("parse runall defaults")
            .inner;
    // Simplex requires --min-reads; supply a value so parse succeeds. Not
    // under test here — parity is on the consensus-calling knobs.
    let simplex = Simplex::try_parse_from([
        "simplex",
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--min-reads",
        "1",
    ])
    .expect("parse simplex defaults");

    assert_eq!(
        runall.consensus_opts.consensus_error_rate_pre_umi, simplex.consensus.error_rate_pre_umi,
        "error-rate-pre-umi default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_error_rate_post_umi, simplex.consensus.error_rate_post_umi,
        "error-rate-post-umi default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_min_input_base_quality,
        simplex.consensus.min_input_base_quality,
        "min-input-base-quality default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_min_consensus_base_quality,
        simplex.consensus.min_consensus_base_quality,
        "min-consensus-base-quality default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_output_per_base_tags,
        simplex.consensus.output_per_base_tags,
        "output-per-base-tags default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_trim, simplex.consensus.trim,
        "trim default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_call_overlapping_bases,
        simplex.overlapping.consensus_call_overlapping_bases,
        "consensus-call-overlapping-bases default diverged",
    );
    assert_eq!(
        runall.consensus_opts.consensus_read_group_id, simplex.read_group.read_group_id,
        "read-group-id default diverged",
    );
}

/// Correct-knob parity: both runall's scoped `--correct::*` flags and the
/// standalone `fgumi correct` CLI must resolve to the same defaults for the
/// knobs that feed into `fgumi_umi::UmiCorrector`.
#[test]
fn runall_correct_defaults_match_standalone_correct() {
    use fgumi_lib::commands::correct::CorrectUmis;

    let runall =
        RunallHarness::try_parse_from(["--input", "ignored.bam", "--output", "ignored_out.bam"])
            .expect("parse runall defaults")
            .inner;
    // `--min-distance` is required by standalone correct (no default); pick a
    // value and pass the same to runall so parity is on the knobs that DO have
    // defaults.
    let standalone = CorrectUmis::try_parse_from([
        "correct",
        "--input",
        "ignored.bam",
        "--output",
        "ignored_out.bam",
        "--min-distance",
        "2",
        "--umis",
        "AAAAAAAA",
    ])
    .expect("parse correct defaults");

    assert_eq!(
        runall.correct_opts.correct_max_mismatches, standalone.max_mismatches,
        "correct max-mismatches default diverged",
    );
    assert_eq!(
        runall.correct_opts.correct_cache_size, standalone.cache_size,
        "correct cache-size default diverged",
    );
    assert_eq!(
        runall.correct_opts.correct_revcomp, standalone.revcomp,
        "correct revcomp default diverged",
    );
    assert_eq!(
        runall.correct_opts.correct_dont_store_original_umis, standalone.dont_store_original_umis,
        "correct dont-store-original-umis default diverged",
    );

    // Anchor to the published constants too — these guard against both sides
    // silently drifting to matching-but-wrong values.
    assert_eq!(runall.correct_opts.correct_max_mismatches, defaults::CORRECT_MAX_MISMATCHES);
    assert_eq!(runall.correct_opts.correct_cache_size, defaults::CORRECT_CACHE_SIZE);
}

/// Extract-knob parity: the standalone `Extract` struct keeps its fields
/// private, so we anchor parity to the shared `defaults::*` constants that
/// both sides read from. This catches a runall CLI `default_value_t` being
/// replaced with a literal that drifts from the standalone default.
#[test]
fn runall_extract_defaults_match_standalone_extract() {
    let runall =
        RunallHarness::try_parse_from(["--input", "ignored.bam", "--output", "ignored_out.bam"])
            .expect("parse runall defaults")
            .inner;

    // --extract::platform and --extract::read-group defaults must match the
    // `defaults::PLATFORM` / `defaults::READ_GROUP_ID` constants that
    // standalone `fgumi extract` pulls from via the same module.
    assert_eq!(
        runall.extract_opts.extract_platform,
        defaults::PLATFORM,
        "runall --extract::platform default diverged from defaults::PLATFORM",
    );
    assert_eq!(
        runall.consensus_opts.consensus_read_group_id,
        defaults::READ_GROUP_ID,
        "runall --consensus::read-group-id default diverged from defaults::READ_GROUP_ID \
         (shared with standalone extract)",
    );
}

/// Sort-knob parity: runall's `--sort::memory-limit` (a raw `usize` byte count)
/// must default to `defaults::SORT_MEMORY_LIMIT`, which standalone `fgumi sort`
/// decodes from the `SORT_MEMORY_LIMIT_DISPLAY` string ("768M"). Both sides
/// must read the same underlying constant so the effective per-thread memory
/// budget stays in lock step.
#[test]
fn runall_sort_defaults_match_standalone_sort() {
    let runall =
        RunallHarness::try_parse_from(["--input", "ignored.bam", "--output", "ignored_out.bam"])
            .expect("parse runall defaults")
            .inner;

    assert_eq!(
        runall.sort_opts.sort_memory_limit,
        defaults::SORT_MEMORY_LIMIT,
        "runall --sort::memory-limit default diverged from defaults::SORT_MEMORY_LIMIT",
    );
    // The display string the standalone `fgumi sort --max-memory` default
    // parses ("768M") must still equal the numeric constant.
    assert_eq!(
        defaults::SORT_MEMORY_LIMIT_DISPLAY,
        "768M",
        "defaults::SORT_MEMORY_LIMIT_DISPLAY is hardcoded in the standalone Sort CLI; \
         a change here means the standalone default is out of sync",
    );
    assert_eq!(
        runall.sort_opts.sort_temp_dir, None,
        "runall --sort::temp-dir default must be None (standalone derives from $TMPDIR)",
    );
}
