//! `--max-temp-files` consolidation runs through the chain's consolidation tail
//! and wakes `SpillWrite` directly: no completion waits on a park timer, and
//! the tail is built only when the limit can trigger.

use std::path::{Path, PathBuf};
use std::time::Duration;

use fgumi_lib::commands::common::{MaxTempFiles, MemoryLimit};
use fgumi_lib::commands::sort::SortOrderArg;
use fgumi_lib::pipeline::chains::{Stage, build_for};
use fgumi_pipeline_core::runtime::{TestBackoff, TestBackoffTarget};
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, write_bam};
use crate::helpers::chain_spec::{sort_options, spec, uniquely_named_records};
use crate::helpers::cutover::decompressed_records_without_pg;

/// The completion wake reaches `SpillWrite` on the pinned writer worker. In a
/// runall-shaped intermediate sort (Sort → Group), `SpillWrite` is a pool Serial
/// writer step; every worker and driver timer is 10 s, so only
/// `ConsolidationRegistry::complete`'s unpark can resume it after each job.
/// Catches: the registry without its waiter (each job then costs a 10 s timer).
#[test]
fn intermediate_sort_consolidation_completes_on_an_unpark() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &uniquely_named_records(20_000));
    let reference =
        run_intermediate(&input, &dir.path().join("reference.bam"), MaxTempFiles::Fixed(64), false);
    // A detached thread (not a scope): on timeout the test fails at once
    // instead of joining a run that sits out 10 s timers.
    let (tx, rx) = std::sync::mpsc::channel();
    let (input2, bounded) = (input.clone(), dir.path().join("bounded.bam"));
    std::thread::spawn(move || {
        tx.send(run_intermediate(&input2, &bounded, MaxTempFiles::Fixed(2), true)).unwrap();
    });
    let bounded = match rx.recv_timeout(Duration::from_secs(5)) {
        Ok(bounded) => bounded,
        Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
            panic!("consolidation waited on a 10 s timer")
        }
        Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
            panic!("the bounded run panicked; see its thread's output above")
        }
    };
    assert!(
        decompressed_records_without_pg(&bounded) == decompressed_records_without_pg(&reference)
    );
}

/// Sort → Group from `input` into `output` at 4 threads, template-coordinate,
/// 64 KiB budget; `slow` puts every worker and driver on a 10 s idle timer.
fn run_intermediate(input: &Path, output: &Path, limit: MaxTempFiles, slow: bool) -> PathBuf {
    let mut sort = sort_options(None, None);
    sort.order = SortOrderArg::TemplateCoordinate;
    sort.max_memory = MemoryLimit::Fixed(64 * 1024);
    sort.memory_per_thread = false;
    sort.max_temp_files = limit;
    let mut built =
        build_for(spec(vec![Stage::Sort, Stage::Group], input, output, 4, sort)).expect("builds");
    let dag = built.pipeline.dag();
    assert!(
        dag.lines().any(|l| l.contains("SpillWrite") && l.contains("Serial")),
        "SpillWrite must be the pool writer step:\n{dag}"
    );
    if slow {
        built.config.test_backoff = vec![
            TestBackoff { target: TestBackoffTarget::AllWorkers, us: 10_000_000 },
            TestBackoff { target: TestBackoffTarget::AllDrivers, us: 10_000_000 },
        ];
    }
    built.run().expect("runs");
    output.to_path_buf()
}

/// Design issue 8: the consolidation tail is built exactly when the limit can
/// trigger (`consolidation_can_trigger`: a merge needs two live runs, so a
/// limit below 2 cannot). The CLI rejects such a limit, but `SortOptions` (the
/// chain's library surface) carries it, so 1 is the boundary that omits the
/// tail. Any limit of 2 or more builds it, however large: the run count is not
/// known when the chain is built, so a limit of 100,000 can still be reached.
/// Catches: the tail, second demand and branch 2 built unconditionally, or
/// omitted for a large limit.
#[test]
fn the_consolidation_tail_is_built_exactly_when_the_limit_can_trigger() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    write_bam(&input, &create_minimal_header("chr1", 1_000_000), &uniquely_named_records(100));
    let dag_for = |limit| {
        let mut sort = sort_options(None, None);
        sort.max_temp_files = limit;
        build_for(spec(vec![Stage::Sort], &input, &dir.path().join("o.bam"), 4, sort))
            .unwrap()
            .pipeline
            .dag()
    };
    let tail = ["ConsolidateRuns", "SpillBlockCompressConsolidation", "WriteConsolidatedRun"];
    let without = dag_for(MaxTempFiles::Fixed(1));
    for step in tail {
        assert!(!without.contains(step), "{step} present at limit 1:\n{without}");
    }
    for limit in [2, 100_000] {
        let with = dag_for(MaxTempFiles::Fixed(limit));
        for step in tail {
            assert!(with.contains(step), "{step} missing at limit {limit}:\n{with}");
        }
    }
}
