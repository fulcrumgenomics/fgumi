//! Structural chain-shape tests for inline consensus metrics (Task 11, Part F4).
//!
//! Spec §7.1's zero-overhead-when-off guarantee is a *chain-build-time* property:
//! when a consensus stage's `--metrics` field is `None`, the built step graph
//! must contain **no** `MetricsCollectorStep` (and no extra output branch /
//! reorder stage feeding one); when it is `Some`, exactly one collector step is
//! wired onto the extra branch. These tests assert that directly on the built
//! `Pipeline`'s `dag()` rendering, rather than through a runtime proxy — the
//! metrics-off path is literally a smaller step graph, built from the existing
//! unmodified step-factory functions.
//!
//! (The black-box output-file companion check — that `--metrics` absent vs.
//! present changes which files land on disk — lives with the parity suite in
//! Task 12.)

use fgumi_lib::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::simplex::SimplexOptions;
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use fgumi_raw_bam::SamBuilder;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};

/// A couple of minimal mapped records. Content is irrelevant here: `build_for`
/// reads only the header (for the consensus sort-order check and reference
/// sequences) and constructs — but does not run — the pipeline, so no records
/// are ever processed.
fn minimal_records() -> Vec<fgumi_raw_bam::RawRecord> {
    (0..2)
        .map(|i| {
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .ref_id(0)
                .pos(100 + i)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            b.add_string_tag(fgumi_lib::sam::SamTag::RX, b"ACGT");
            b.add_string_tag(fgumi_lib::sam::SamTag::MI, format!("{i}").as_bytes());
            b.build()
        })
        .collect()
}

/// Build a standalone `[Stage::Simplex]` chain over a tiny template-coordinate
/// BAM and return its `Pipeline::dag()` string. `metrics` toggles the inline
/// `--metrics` field on the simplex stage.
fn simplex_chain_dag(metrics: Option<std::path::PathBuf>) -> String {
    let dir = TempDir::new().expect("temp dir");
    let input_bam = dir.path().join("in.bam");
    let output_bam = dir.path().join("out.bam");

    let header = create_minimal_header("chr1", 1_000_000);
    write_bam(&input_bam, &header, &minimal_records());

    let simplex = SimplexOptions { metrics, ..Default::default() };

    let spec = ChainSpec {
        stages: vec![Stage::Simplex],
        source: SourceSpec::Bam(input_bam),
        sink: SinkSpec::Bam(output_bam),
        stage_opts: StageOptionsBag { simplex: Some(simplex), ..Default::default() },
        threading: ThreadingOptions { threads: None },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: true,
        command_line: "fgumi simplex".to_string(),
    };

    build_for(spec).expect("build_for should accept a standalone Simplex chain").pipeline.dag()
}

#[test]
fn metrics_off_chain_has_no_metrics_collector_step() {
    let dag = simplex_chain_dag(None);
    assert!(
        !dag.contains("MetricsCollectorStep"),
        "a --metrics-absent simplex chain must not wire any MetricsCollectorStep:\n{dag}"
    );
}

#[test]
fn metrics_on_chain_wires_exactly_one_metrics_collector_step() {
    let dir = TempDir::new().expect("temp dir");
    let prefix = dir.path().join("metrics_prefix");
    let dag = simplex_chain_dag(Some(prefix));
    // Count step-definition lines (`  [N] MetricsCollectorStep ...`), not raw
    // substring hits — `dag()` also names the step as the *consumer* on the
    // reorder branch feeding it (`.0: ... → MetricsCollectorStep`).
    let collector_steps = dag
        .lines()
        .filter(|l| l.trim_start().starts_with('[') && l.contains("MetricsCollectorStep"))
        .count();
    assert_eq!(
        collector_steps, 1,
        "a --metrics-present simplex chain must wire exactly one MetricsCollectorStep step:\n{dag}"
    );
}

#[test]
fn metrics_on_chain_adds_steps_over_metrics_off() {
    let dir = TempDir::new().expect("temp dir");
    let prefix = dir.path().join("metrics_prefix");
    let dag_off = simplex_chain_dag(None);
    let dag_on = simplex_chain_dag(Some(prefix));

    // The metrics-on chain adds the extra CoordinateGroupFragment branch, its
    // auto-inserted ReorderStage, and the terminal MetricsCollectorStep — so it
    // has strictly more steps than the metrics-off chain. This is spec §7.1's
    // "zero steps added when off" stated as a count.
    let steps_off = dag_off.lines().filter(|l| l.trim_start().starts_with('[')).count();
    let steps_on = dag_on.lines().filter(|l| l.trim_start().starts_with('[')).count();
    assert!(
        steps_on > steps_off,
        "metrics-on chain ({steps_on} steps) must have more steps than metrics-off \
         ({steps_off} steps)\n--- off ---\n{dag_off}\n--- on ---\n{dag_on}"
    );
}
