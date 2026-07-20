//! Build-level chain-topology tests (finding S5b1-010).
//!
//! These tests call [`fgumi_lib::pipeline::chains::build_for`] directly and
//! assert it returns `Ok` for fused stage permutations — i.e. the chain
//! *constructs and wires* without a build-time error (missing per-stage
//! options, unwired output branch, sort-order rejection, …). They are
//! deliberately **build-only**: command-mode `--aligner::command` lets a
//! `[Correct, AlignAndMerge]` chain build without spawning an aligner.
//!
//! Scope note: a wrong intermediate→terminal *handle type* is NOT caught here.
//! `build_for` erases handles and builds a mis-typed chain successfully; the
//! "input handle downcast failed — chain topology invariant" check is a
//! **run-time** panic raised at the first dispatch in `TypedStep::resolve_input`
//! (see `fgumi_pipeline_core::builder`). That run-time invariant is exercised
//! end-to-end by the fused mock-aligner runs in `test_runall_parity.rs` (e.g.
//! the `correct → … → align → …` parity tests, which drive the same
//! `Correct` emits `BamTemplateBatch` → `Align` consumes `BamTemplateBatch`
//! (skipping `GroupByQueryname`) hand-off through actual dispatch). What these
//! build-only tests add is cheap, aligner-free confirmation that the chain
//! builder accepts the permutation in the first place.

#![allow(clippy::needless_pass_by_value)]

use std::fs;
use std::path::{Path, PathBuf};

use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_minimal_header, create_paired_umi_family_at, create_test_reference, to_record_buf,
};

/// Writes a tiny template-coordinate-headed, paired-end, UMI-tagged BAM to
/// `dir/input.bam` and returns its path.
///
/// The header carries `SS:template-coordinate` (via [`create_minimal_header`])
/// so the `Group` stage's sort-order check accepts it, and the records carry
/// `RX` + `MC` tags so both `Correct` (reads `RX`) and `Group` (needs `MC` for
/// paired-end template spans) build cleanly. `ChainBuilder::new` opens this
/// file to read the header, so it must be a real on-disk BAM.
fn write_input_bam(dir: &Path) -> PathBuf {
    let path = dir.join("input.bam");
    let header = create_minimal_header("chr1", 10000);
    let mut writer = bam::io::Writer::new(fs::File::create(&path).expect("create input BAM"));
    writer.write_header(&header).expect("write header");
    // Two paired-UMI families at distinct positions — enough to exercise the
    // chain topology without needing many records (these tests only build, they
    // do not run).
    let mut records = Vec::new();
    records.extend(create_paired_umi_family_at(
        "ACGTACGT",
        2,
        "fam_a",
        "ACGTACGTACGT",
        "TGCATGCATGCA",
        30,
        100,
    ));
    records.extend(create_paired_umi_family_at(
        "TGCATGCA",
        2,
        "fam_b",
        "ACGTACGTACGT",
        "TGCATGCATGCA",
        30,
        500,
    ));
    for raw in &records {
        writer.write_alignment_record(&header, &to_record_buf(raw)).expect("write record");
    }
    writer.try_finish().expect("finish input BAM");
    path
}

/// Builds a complete [`ChainSpec`] with the given stages, source BAM, and
/// pre-filled option bag. Mirrors `runall`'s `ChainSpec` construction
/// (see `commands::runall::execute`) with default threading/compression/etc.
fn spec_for(stages: Vec<Stage>, source_bam: &Path, stage_opts: StageOptionsBag) -> ChainSpec {
    use fgumi_lib::commands::common::{
        CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    };
    ChainSpec {
        stages,
        source: SourceSpec::Bam(source_bam.to_path_buf()),
        sink: SinkSpec::Bam(PathBuf::from("/dev/null")),
        stage_opts,
        threading: ThreadingOptions { threads: Some(1) },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        command_line: "fgumi test-chain-build".to_string(),
    }
}

// ─────────────────────── option-bag constructors ───────────────────────
//
// Each mirrors the corresponding arm of `runall::build_stage_options_bag`:
// fill the required fields and the `#[arg(skip)]` slots the chain builder
// reads, leaving everything else at its `Default`.

/// `Stage::Sort` options: template-coordinate order (the only order Group
/// accepts downstream), matching runall's forced `SortOrderArg`.
fn sort_opts() -> fgumi_lib::commands::sort::SortOptions {
    use fgumi_lib::commands::sort::{SortOptions, SortOrderArg};
    SortOptions { order: SortOrderArg::TemplateCoordinate, ..Default::default() }
}

/// `Stage::Group` options for the given strategy, with the `#[arg(skip)]`
/// `effective_*` slots populated exactly as `runall` does.
fn group_opts(strategy: fgumi_lib::assigner::Strategy) -> fgumi_lib::commands::group::GroupOptions {
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::group::GroupOptions;
    let edits = u32::from(!matches!(strategy, Strategy::Identity));
    GroupOptions {
        strategy,
        edits,
        effective_strategy: strategy,
        effective_edits: edits,
        ..Default::default()
    }
}

/// `Stage::Correct` options with an inline UMI whitelist (`-u`), so no
/// whitelist file is needed; `add_correct` reads the UMIs at build time.
fn correct_opts() -> fgumi_lib::commands::correct::CorrectOptions {
    use fgumi_lib::commands::correct::CorrectOptions;
    CorrectOptions {
        umis: vec!["ACGTACGT".to_string(), "TGCATGCA".to_string()],
        ..Default::default()
    }
}

// ─────────────────────────────── tests ───────────────────────────────

/// `[Sort, Group]` — the simplest two-stage fusion. Sort emits the
/// intermediate record stream that Group consumes; pins that hand-off.
#[test]
fn chain_build_sort_to_group() {
    use fgumi_lib::assigner::Strategy;
    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    let bag = StageOptionsBag {
        sort: Some(sort_opts()),
        group: Some(group_opts(Strategy::Adjacency)),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Sort, Stage::Group], &input, bag);
    build_for(spec).expect("build_for([Sort, Group]) should succeed");
}

/// `[Correct, AlignAndMerge]` — the only permutation with no *build-level*
/// coverage elsewhere (parity tests in `test_runall_parity.rs` already exercise
/// the same `Correct → Align` runtime hand-off via fused mock-aligner runs;
/// what is missing is a build-only check). It builds without an aligner:
/// command-mode `--aligner::command` only substitutes `{ref}` into a template,
/// and the reference needs only a `.dict` sibling (provided by
/// `create_test_reference`), not real index files. This pins the `build_for`
/// wiring of the `Correct` (emits `BamTemplateBatch`) → `Align` (consumes
/// `BamTemplateBatch`, skipping `GroupByQueryname`) topology hand-off — exactly
/// the `chain_tail_kind` branch a regression would break at construction time.
#[test]
fn chain_build_correct_to_align() {
    use fgumi_lib::aligner::AlignerOptions;
    use fgumi_lib::pipeline::chains::options_bag::AlignOptions;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());
    let reference = create_test_reference(tmp.path());

    // Command mode: a do-nothing template that contains the required `{ref}`
    // placeholder. `AlignerOptions::resolve` substitutes it and never spawns
    // a process at build time, so no aligner binary is required.
    let bag = StageOptionsBag {
        correct: Some(correct_opts()),
        aligner: Some(AlignOptions {
            aligner: AlignerOptions {
                command: Some("cat {ref} >/dev/null; cat".to_string()),
                ..Default::default()
            },
            reference,
            aligner_bin: None,
        }),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Correct, Stage::Align], &input, bag);
    build_for(spec).expect(
        "build_for([Correct, Align]) should construct the chain — this covers chain \
         construction/topology wiring only; the typed input-handle downcast is a \
         run-time concern checked by `Pipeline::run` (`TypedStep::resolve_input`), not here",
    );
}

// Consensus-terminal fusions. Gated on `consensus` exactly like the bag
// slots (`StageOptionsBag::{simplex,duplex,codec}`) so the
// `--no-default-features` build still compiles.

/// `[Group, Simplex]` — Group's grouped record stream feeds the simplex
/// consensus caller (the terminal stage). Pins that hand-off at build time.
#[cfg(feature = "consensus")]
#[test]
fn chain_build_group_to_simplex() {
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::simplex::SimplexOptions;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    let bag = StageOptionsBag {
        group: Some(group_opts(Strategy::Adjacency)),
        simplex: Some(SimplexOptions::default()),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Group, Stage::Simplex], &input, bag);
    build_for(spec).expect("build_for([Group, Simplex]) should succeed");
}

/// V5 (audit): `runall` builds the chain via `build_for` and never calls
/// `Simplex::execute`, so the `--min-reads` / `--max-reads` range check the
/// standalone command enforces must be re-applied in `add_simplex`. This pins
/// the parity contract: `build_for` must *accept* exactly the read-bound
/// configurations the standalone `fgumi simplex` command accepts and *reject*
/// exactly the degenerate ones it rejects.
///
/// The `expected` column is an explicit expected-validity oracle, authored
/// independently of the validator rather than derived from it: `None` means the
/// configuration is valid and the chain must build; `Some(msg)` means it is
/// degenerate and `build_for` must reject it with an error containing `msg`.
/// Exercising both the accept and reject decisions — not just echoing the
/// validator's own error text — means a shared regression cannot make the test
/// pass: a validator gutted to accept everything breaks the `Some(..)` cases,
/// and one that rejects everything breaks the `None` cases.
#[cfg(feature = "consensus")]
#[rstest::rstest]
#[case::default_bounds(1, None, None)]
#[case::min_one_max_one(1, Some(1), None)]
#[case::min_two_max_five(2, Some(5), None)]
#[case::min_reads_zero(0, None, Some("--min-reads must be >= 1"))]
#[case::max_below_min(5, Some(2), Some("--max-reads (2) must be >= --min-reads (5)"))]
fn chain_build_simplex_enforces_read_bounds(
    #[case] min_reads: usize,
    #[case] max_reads: Option<usize>,
    #[case] expected: Option<&str>,
) {
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::simplex::SimplexOptions;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    let build = || {
        let simplex = SimplexOptions { min_reads, max_reads, ..SimplexOptions::default() };
        let bag = StageOptionsBag {
            group: Some(group_opts(Strategy::Adjacency)),
            simplex: Some(simplex),
            ..Default::default()
        };
        build_for(spec_for(vec![Stage::Group, Stage::Simplex], &input, bag))
    };

    match expected {
        None => {
            build().unwrap_or_else(|err| {
                panic!(
                    "build_for([Group, Simplex]) with min_reads={min_reads} max_reads={max_reads:?} must be accepted, got: {err}"
                )
            });
        }
        Some(msg) => {
            let Err(err) = build() else {
                panic!(
                    "build_for([Group, Simplex]) with min_reads={min_reads} max_reads={max_reads:?} must be rejected"
                );
            };
            assert!(
                err.to_string().contains(msg),
                "expected simplex read-bounds error containing {msg:?}, got: {err}"
            );
        }
    }
}

/// `[Group, Duplex]` — duplex requires the `Paired` grouping strategy (MIs
/// carry `/A`/`/B` suffixes), enforced by the cross-stage validator. Pins
/// the Group→Duplex hand-off under that constraint.
#[cfg(feature = "consensus")]
#[test]
fn chain_build_group_to_duplex() {
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::duplex::DuplexOptions;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    let bag = StageOptionsBag {
        group: Some(group_opts(Strategy::Paired)),
        duplex: Some(DuplexOptions::default()),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Group, Stage::Duplex], &input, bag);
    build_for(spec).expect("build_for([Group, Duplex]) should succeed");
}

/// `[Group, Codec]` — Group's grouped stream feeds the CODEC consensus
/// caller. Pins the Group→Codec hand-off at build time.
#[cfg(feature = "consensus")]
#[test]
fn chain_build_group_to_codec() {
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::codec::CodecOptions;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    let bag = StageOptionsBag {
        group: Some(group_opts(Strategy::Adjacency)),
        codec: Some(CodecOptions::default()),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Group, Stage::Codec], &input, bag);
    build_for(spec).expect("build_for([Group, Codec]) should succeed");
}

/// `[Group, Clip]` — an ordering-valid but type-incompatible chain: Group emits
/// grouped templates (`BamTemplateBatch`) while Clip consumes a record stream
/// (`DecodedRecordBatch`, which it groups internally). The stage-ordering
/// validators pass, but the type-erased pipeline would panic at dispatch. The
/// `chain_tail_kind` guard in `add_clip` must reject it cleanly at build time.
#[test]
fn chain_build_group_to_clip_rejected() {
    use clap::Parser;
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::clip::Clip;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    // Minimal parse only — the guard fires before Clip's own options are read.
    let clip =
        Clip::try_parse_from(["clip", "-i", "/dev/null", "-o", "/dev/null", "-r", "/dev/null"])
            .expect("parse minimal clip");

    let bag = StageOptionsBag {
        group: Some(group_opts(Strategy::Adjacency)),
        clip: Some(clip),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Group, Stage::Clip], &input, bag);
    let Err(err) = build_for(spec) else {
        panic!("Group -> Clip must be rejected at build time");
    };
    assert!(
        err.to_string().contains("Stage::Clip requires a record-stream input (DecodedRecordBatch)"),
        "expected the clip tail-kind guard to fire, got: {err}"
    );
}

/// Filter is the same class of record-stream consumer as Clip/Dedup, so a `Group -> Filter`
/// chain (grouped-template tail feeding a record-stream stage) must also be rejected at build
/// time rather than building and panicking at dispatch.
#[test]
fn chain_build_group_to_filter_rejected() {
    use fgumi_lib::assigner::Strategy;
    use fgumi_lib::commands::filter::FilterOptions;

    let tmp = TempDir::new().unwrap();
    let input = write_input_bam(tmp.path());

    // Default opts — the guard fires before Filter's own options are read.
    let bag = StageOptionsBag {
        group: Some(group_opts(Strategy::Adjacency)),
        filter: Some(FilterOptions::default()),
        ..Default::default()
    };

    let spec = spec_for(vec![Stage::Group, Stage::Filter], &input, bag);
    let Err(err) = build_for(spec) else {
        panic!("Group -> Filter must be rejected at build time");
    };
    assert!(
        err.to_string()
            .contains("Stage::Filter requires a record-stream input (DecodedRecordBatch)"),
        "expected the filter tail-kind guard to fire, got: {err}"
    );
}
