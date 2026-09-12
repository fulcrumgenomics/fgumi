//! Structural chain-shape tests for inline consensus metrics (Task 11, Part F4;
//! inverted for the parallel T2 producer, Task 3).
//!
//! Spec §7.1's zero-overhead-when-off guarantee is a *chain-build-time* property:
//! when a consensus stage's `--metrics` field is `None`, the built step graph
//! must not gain any extra step, output branch, or reorder stage over the
//! metrics-off shape. As of Task 3, this now holds for `--metrics` `Some` too:
//! the standalone (T2) simplex path records metrics inline in the consensus
//! worker body via a per-thread `ConsensusMetricsSlot` accumulator instead of
//! fanning out to a separate serial collector step, so the metrics-on chain has
//! the SAME step count/shape as metrics-off. These tests assert that directly
//! on the built `Pipeline`'s `dag()` rendering, rather than through a runtime
//! proxy — the metrics-on and metrics-off chains are literally the same step
//! graph, built from the existing unmodified step-factory functions.
//!
//! (An earlier revision of this file asserted the built DAG never rendered the
//! string `"MetricsCollectorStep"` — the name of a now-deleted serial-collector
//! step type from a prior design. Once that type was removed workspace-wide the
//! assertion became vacuous: the string can no longer appear regardless of
//! whether the chain is actually shaped correctly, since no code path emits it.
//! The tests below assert on the DAG's actual step types/ordering instead, so
//! they fail if the metrics-on chain ever again grows an extra step, whatever
//! that step happens to be named.)
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

/// Extract the step *type name* token from each `dag()` step-definition line
/// (`  [{idx}] {name:<24} {kind:?} sticky=... branches=...`), in step order.
/// This deliberately ignores per-step numbering/queue/ordering detail and
/// compares only the sequence of step types actually wired into the chain.
fn dag_step_names(dag: &str) -> Vec<&str> {
    dag.lines()
        .filter(|l| l.trim_start().starts_with('['))
        .map(|l| l.split_whitespace().nth(1).unwrap_or(""))
        .collect()
}

#[test]
fn metrics_on_chain_has_same_step_shape_as_off() {
    let dir = TempDir::new().expect("temp dir");
    let prefix = dir.path().join("metrics_prefix");
    let dag_off = simplex_chain_dag(None);
    let dag_on = simplex_chain_dag(Some(prefix));

    // Stronger than a bare step *count* match (see the test below): this
    // compares the sequence of step *types*, in order, between the
    // metrics-on and metrics-off chains. Metrics recording moved inline into
    // the consensus worker body (a per-thread `ConsensusMetricsSlot`
    // accumulator) rather than fanning out to a separate collector step, so
    // toggling `--metrics` must not add, remove, or substitute any step —
    // the two chains must wire the identical step-type sequence.
    let names_off = dag_step_names(&dag_off);
    let names_on = dag_step_names(&dag_on);
    assert_eq!(
        names_on, names_off,
        "metrics-on chain must wire the SAME step types in the SAME order as metrics-off \
         (no extra or substituted step)\n--- off ---\n{dag_off}\n--- on ---\n{dag_on}"
    );
}

#[test]
fn metrics_on_chain_has_same_step_count_as_off() {
    let dir = TempDir::new().expect("temp dir");
    let prefix = dir.path().join("metrics_prefix");
    let dag_off = simplex_chain_dag(None);
    let dag_on = simplex_chain_dag(Some(prefix));

    // The headline zero-overhead statement (Task 3): metrics collection now
    // happens inside the existing consensus worker, adding no extra branch and
    // no auto-inserted reorder or collector step — so the metrics-on chain has
    // the SAME step count as metrics-off. This is a weaker check than
    // `metrics_on_chain_has_same_step_shape_as_off` above (a count match alone
    // wouldn't catch a step *substitution*), kept because it is the simplest
    // direct statement of the guarantee and fails independently of how `dag()`
    // renders individual step names.
    let steps_off = dag_off.lines().filter(|l| l.trim_start().starts_with('[')).count();
    let steps_on = dag_on.lines().filter(|l| l.trim_start().starts_with('[')).count();
    assert_eq!(
        steps_on, steps_off,
        "metrics-on chain ({steps_on} steps) must have the SAME step count as metrics-off \
         ({steps_off} steps)\n--- off ---\n{dag_off}\n--- on ---\n{dag_on}"
    );
}
