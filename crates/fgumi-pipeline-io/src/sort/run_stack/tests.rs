use std::collections::HashMap;

use rstest::rstest;

use super::*;

/// Bytes each simulated spill run occupies on disk.
const RUN_BYTES: u64 = 1000;

fn run(file_id: u32, bytes: u64) -> SpillRun {
    SpillRun {
        file_id,
        path: PathBuf::from(format!("chunk_{file_id}")),
        bytes,
        records_ingested_so_far: u64::from(file_id),
    }
}

/// Outcome of pushing `sizes.len()` runs through a stack with the given limit,
/// settling it after every push exactly as `SpillWrite` does.
struct Simulation {
    /// Surviving runs, each as the ordered list of original runs it contains.
    survivors: Vec<Vec<u32>>,
    /// Largest number of live runs observed once a push had been settled.
    peak_settled: usize,
    /// Width of every consolidation, in order.
    widths: Vec<usize>,
    /// Bytes read by consolidations, summed.
    bytes_rewritten: u64,
}

fn simulate(max_temp_files: usize, sizes: &[u64]) -> Simulation {
    let mut stack = RunStack::new(max_temp_files);
    // A merged run's path names the originals it holds.
    let mut contents: HashMap<PathBuf, Vec<u32>> = HashMap::new();
    let (mut peak_settled, mut widths, mut bytes_rewritten) = (0, Vec::new(), 0);
    for (id, &bytes) in sizes.iter().enumerate() {
        let r = run(u32::try_from(id).unwrap(), bytes);
        contents.insert(r.path.clone(), vec![r.file_id]);
        stack.push(r);
        while let Some(range) = stack.next_merge() {
            assert!(range.len() >= 2, "a merge must reduce the run count");
            let merged_path = PathBuf::from(format!("merged_{}", widths.len()));
            widths.push(range.len());
            let inputs = &stack.runs()[range.clone()];
            let merged_bytes: u64 = inputs.iter().map(|r| r.bytes).sum();
            bytes_rewritten += merged_bytes;
            let originals: Vec<u32> =
                inputs.iter().flat_map(|r| contents[&r.path].clone()).collect();
            contents.insert(merged_path.clone(), originals);
            stack.complete_merge(range, merged_path, merged_bytes);
        }
        peak_settled = peak_settled.max(stack.runs().len());
    }
    let survivors = stack.into_runs().iter().map(|r| contents[&r.path].clone()).collect();
    Simulation { survivors, peak_settled, widths, bytes_rewritten }
}

/// Bytes the legacy `RawExternalSorter` policy rewrites for the same runs: each time the
/// live count reaches the limit, merge the oldest `limit / 2` runs and put the
/// result back at the front.
fn legacy_policy_bytes(max_temp_files: usize, sizes: &[u64]) -> u64 {
    let take = (max_temp_files / 2).max(2);
    let mut live: Vec<u64> = Vec::new();
    let mut rewritten = 0;
    for &bytes in sizes {
        live.push(bytes);
        if live.len() >= max_temp_files {
            let merged: u64 = live.drain(..take).sum();
            rewritten += merged;
            live.insert(0, merged);
        }
    }
    rewritten
}

fn uniform(runs: usize) -> Vec<u64> {
    vec![RUN_BYTES; runs]
}

#[rstest]
fn settled_stack_respects_the_limit_and_preserves_run_order(
    #[values(2, 3, 4, 5, 8, 64, 1024)] max_temp_files: usize,
    #[values(0, 1, 2, 3, 31, 63, 64, 65, 1024, 1025, 3000)] total_runs: usize,
) {
    let sim = simulate(max_temp_files, &uniform(total_runs));

    assert!(
        sim.peak_settled < max_temp_files || total_runs < 2,
        "L={max_temp_files} R={total_runs}: {} live runs after settling",
        sim.peak_settled
    );
    let fan_in = (max_temp_files / 2).clamp(2, MAX_FAN_IN);
    assert!(sim.widths.iter().all(|&w| (2..=fan_in).contains(&w)), "widths {:?}", sim.widths);
    // Every original run survives exactly once, in the original order — the
    // property that keeps the final merge's tie-break (and so the output)
    // unchanged.
    let flattened: Vec<u32> = sim.survivors.concat();
    assert_eq!(flattened, (0..u32::try_from(total_runs).unwrap()).collect::<Vec<_>>());
}

/// The limit bounds open files; it is not a target. Below it nothing is merged.
#[rstest]
#[case::default_limit(1024, 1000)]
#[case::benchmark_limit(64, 63)]
fn under_the_limit_nothing_is_merged(#[case] max_temp_files: usize, #[case] total_runs: usize) {
    let sim = simulate(max_temp_files, &uniform(total_runs));
    assert!(sim.widths.is_empty(), "merged {:?} under the limit", sim.widths);
    assert_eq!(sim.survivors.len(), total_runs);
}

/// Just past the limit — one consolidation needed — the policy does exactly what
/// the legacy engine does: a single merge of `L / 2` runs. This is the regime the
/// spill-consolidation benchmarks measure, so their numbers stay comparable.
#[test]
fn one_pass_past_the_limit_matches_the_legacy_engine() {
    let sizes = uniform(87);
    let sim = simulate(64, &sizes);
    assert_eq!(sim.widths, vec![32]);
    assert_eq!(sim.bytes_rewritten, legacy_policy_bytes(64, &sizes));
}

/// When the limit is hit repeatedly, merging the cheapest window rewrites far
/// less than the legacy oldest-half policy, whose first merged run absorbed the
/// oldest records on every later pass.
#[rstest]
#[case(64, 350)]
#[case(64, 700)]
#[case(64, 5000)]
#[case(16, 264)]
#[case(8, 264)]
#[case(4, 264)]
fn repeated_consolidation_rewrites_less_than_the_legacy_policy(
    #[case] max_temp_files: usize,
    #[case] total_runs: usize,
) {
    let sizes = uniform(total_runs);
    let ours = simulate(max_temp_files, &sizes).bytes_rewritten;
    let legacy = legacy_policy_bytes(max_temp_files, &sizes);
    assert!(ours * 2 < legacy, "L={max_temp_files} R={total_runs}: {ours} vs legacy {legacy}");
}

/// Sizes, not counts, drive the choice: a window of small runs is merged before
/// a window of the same width holding a large one.
#[test]
fn the_cheapest_window_is_merged_first() {
    let mut stack = RunStack::new(4);
    for (id, bytes) in [(0, 10), (1, 1_000_000), (2, 10), (3, 10)] {
        stack.push(run(id, bytes));
    }
    // Windows of 2: [0,1] and [1,2] include the large run; [2,3] costs 20 for
    // one slot. Windows of 3 each include the large run.
    assert_eq!(stack.next_merge(), Some(2..4));
}

#[test]
fn a_limit_below_two_disables_consolidation() {
    for limit in [0, 1] {
        let sim = simulate(limit, &uniform(100));
        assert!(sim.widths.is_empty(), "limit {limit}");
        assert_eq!(sim.survivors.len(), 100);
    }
}

#[test]
fn merged_run_takes_the_position_and_lowest_file_id_of_its_range() {
    let mut stack = RunStack::new(4);
    for id in [3, 5, 8] {
        stack.push(run(id, RUN_BYTES));
    }
    let absorbed = stack.complete_merge(1..3, PathBuf::from("merged"), 1999);
    assert_eq!(absorbed.iter().map(|r| r.file_id).collect::<Vec<_>>(), vec![5, 8]);
    let runs = stack.runs();
    assert_eq!(runs.len(), 2);
    assert_eq!(runs[0].file_id, 3);
    assert_eq!((runs[1].file_id, runs[1].bytes), (5, 1999));
    assert_eq!(runs[1].path, PathBuf::from("merged"));
    assert_eq!(runs[1].records_ingested_so_far, 8);
}
