# Contributing to fgumi

## Development Setup

### Prerequisites

- Rust (stable toolchain)
- cargo-nextest (for running tests)

### Install Git Hooks

We use pre-commit hooks to ensure code quality. Install them after cloning:

```bash
./scripts/install-hooks.sh
```

This installs hooks that run before each commit:
- `cargo ci-fmt` - Check code formatting
- `cargo ci-lint` - Run clippy lints

### Running Checks Manually

```bash
# Format check (fails if formatting differs)
cargo ci-fmt

# Lint check (fails on any warnings)
cargo ci-lint

# Run all tests
cargo ci-test
```

### Pre-Commit Hook Options

**Run tests in pre-commit hook:**
```bash
FGUMI_PRECOMMIT_TEST=1 git commit -m "message"
```

**Bypass hooks (use sparingly):**
```bash
git commit --no-verify -m "message"
```

## Code Style

- Run `cargo fmt` before committing
- Fix all clippy warnings
- Add backticks around identifiers in doc comments (e.g., `` `read_name` ``)

## Testing

All new features should include tests. Run the full test suite with:

```bash
cargo ci-test
```

## Pull Requests

1. Ensure all CI checks pass (`cargo ci-fmt`, `cargo ci-lint`, `cargo ci-test`)
2. Keep PRs focused and reasonably sized (250-1000 LOC ideal)
3. Include tests for new functionality

## Extending fgumi

The following recipes walk through the three most common ways contributors
extend fgumi. They assume you have already read the relevant module-level docs
(linked from each recipe). Each recipe points at a reference implementation to
copy-and-adapt rather than writing from scratch.

### Adding a new pipeline stage

Pipeline stages are the units of work the runall engine schedules. There are
two trait choices — see the decision tree in `src/lib/pipeline/stage.rs` (and
`src/lib/pipeline/special_stage.rs`) — but most new stages should implement
`Stage`. Pick `SpecialStage` only if your stage emits 0-or-many outputs per
input, has barrier semantics, or owns subprocess / worker threads.

1. **Pick input and output batch types.** Look at existing stages in
   `src/lib/pipeline/stages/` for inspiration. Common batch types live in
   `src/lib/pipeline/fastq_types.rs`, `grouping_types.rs`, and
   `output_types.rs` (e.g. `SerializedBatch` for bytes-of-BAM, `MiGroupBatch`
   for reads grouped by MI tag, `PositionBatch` for 5'-position batches).

2. **Create the stage file** at `src/lib/pipeline/stages/my_stage.rs`. Use
   `src/lib/pipeline/stages/count.rs` as the minimal template — it has the
   full `## Input / ## Output / ## Ordering guarantees / ## Memory model /
   ## Determinism` doc block every stage is expected to carry. Skeleton:

   ```rust
   use crate::pipeline::stage::{Parallelism, Stage};

   #[derive(Clone)] // required for Parallelism::Parallel unless you use stage_factory
   pub struct MyStage { /* fields */ }

   impl Stage for MyStage {
       type Input = /* batch type in */;
       type Output = /* batch type out */;

       fn process(
           &mut self,
           input: Self::Input,
           emit: &mut dyn FnMut(Self::Output),
       ) -> anyhow::Result<()> {
           // transform, call `emit(output)` 0 or 1 times
           Ok(())
       }

       fn parallelism(&self) -> Parallelism { Parallelism::Parallel }

       fn output_memory_estimate(&self, output: &Self::Output) -> usize {
           // heap bytes owned by `output`, not size_of::<Output>()
           output.bytes_owned()
       }

       fn name(&self) -> &'static str { "MyStage" }
   }
   ```

3. **Register the module.** Add `pub mod my_stage;` to
   `src/lib/pipeline/stages/mod.rs`.

4. **Wire it into the plan.** Decide where the stage sits in the
   `--start-from` / `--stop-after` chain. Add a new `StageSpec::MyStage {
   ... }` variant in `src/lib/commands/runall/plan.rs`, push it from the right
   place in `src/lib/commands/runall/planner.rs`, and instantiate it in
   `src/lib/commands/runall/runner.rs`.

5. **Test.** Use the stage harness in
   `tests/integration/helpers/stage_harness.rs`:

   ```rust
   use crate::helpers::stage_harness::run_stage;

   let outputs = run_stage(MyStage::new(/* ... */), vec![input_batch]);
   assert_eq!(outputs.len(), 1);
   ```

   If the stage affects end-to-end output, add a byte-identity gate under
   `tests/integration/runall/` alongside the existing `test_sort.rs`,
   `test_zipper.rs`, etc.

### Adding a new consensus caller

Consensus callers live in the `fgumi-consensus` crate. The `ConsensusCaller`
trait is documented at the top of `crates/fgumi-consensus/src/caller.rs`; use
`crates/fgumi-consensus/src/vanilla_caller.rs` (`VanillaUmiConsensusCaller`)
as the reference implementation.

1. **Create the caller file** at
   `crates/fgumi-consensus/src/my_caller.rs` and implement
   `fgumi_consensus::caller::ConsensusCaller`:

   ```rust
   use crate::caller::{ConsensusCaller, ConsensusCallingStats, ConsensusOutput};
   use crate::error::Result;

   pub struct MyConsensusCaller { stats: ConsensusCallingStats, /* ... */ }

   impl ConsensusCaller for MyConsensusCaller {
       fn consensus_reads(&mut self, records: Vec<Vec<u8>>) -> Result<ConsensusOutput> {
           self.stats.record_input(records.len());
           // filter, call consensus, write raw BAM records into output
           Ok(ConsensusOutput::new())
       }
       fn total_reads(&self) -> usize { self.stats.total_reads }
       fn total_filtered(&self) -> usize { self.stats.filtered_reads }
       fn consensus_reads_constructed(&self) -> usize { self.stats.consensus_reads }
       fn statistics(&self) -> ConsensusCallingStats { self.stats.clone() }
       fn log_statistics(&self) { /* use log_consensus_statistics helper */ }
       fn clear(&mut self) { self.stats = ConsensusCallingStats::new(); /* reset rng etc */ }
   }
   ```

   Each `Vec<u8>` is a raw BAM record **without** the 4-byte `block_size`
   prefix. The returned `ConsensusOutput::data` must be raw BAM records
   **with** `block_size` prefixes so BGZF compression can consume it
   directly. `clear()` must wipe all per-group transient state (stats,
   rejection buffer, RNG cursor) — otherwise output becomes input-order
   dependent.

2. **Register the module and public re-exports** in
   `crates/fgumi-consensus/src/lib.rs` and the matching re-exports in
   `src/lib/consensus/mod.rs`.

3. **Add a `ConsensusMode` variant** in
   `src/lib/commands/runall/options.rs` and a
   matching match arm in `ConsensusMode::build_factory`
   (`src/lib/commands/runall/consensus_helpers.rs`). The factory closure
   returns a fresh caller per pool worker — probe-build once at factory
   construction time to surface constructor errors eagerly, following the
   `Simplex` / `Duplex` branches as a template.

4. **Add integration tests** under `tests/integration/runall/` — parallel to
   `test_consensus_simplex.rs`, `test_consensus_duplex.rs`,
   `test_consensus_codec.rs`.

5. **Update the runall mdbook chapter** at `docs/src/guide/runall.md` to
   mention the new `--consensus-mode` value.

### Adding a new UMI assigner

UMI assigners live in the `fgumi-umi` crate. The `UmiAssigner` trait is
defined in `crates/fgumi-umi/src/assigner.rs`; use `AdjacencyUmiAssigner` in
the same file as the reference implementation (it is the most feature-complete
and multi-threaded).

1. **Add the struct and `impl UmiAssigner`** in
   `crates/fgumi-umi/src/assigner.rs` (or a new submodule you register from
   `lib.rs`):

   ```rust
   pub struct MyAssigner { /* ... */ }

   impl UmiAssigner for MyAssigner {
       fn assign(&self, raw_umis: &[Umi]) -> Vec<MoleculeId> {
           // return one MoleculeId per input UMI, same order as input
           unimplemented!()
       }
       // optional: override `same_umi(&self, a, b)` if equality is non-trivial
       // (e.g. paired UMIs where A-B == B-A)
   }
   ```

   Contract: `result.len() == raw_umis.len()` and `result[i]` is the molecule
   ID for `raw_umis[i]`. UMIs grouped together share the same `MoleculeId`.
   Implementations must be `Send + Sync` for parallel processing.

2. **Add a `Strategy` variant** in `crates/fgumi-umi/src/assigner.rs`
   alongside `Identity`, `Edit`, `Adjacency`, `Paired`.

3. **Dispatch from `Strategy::new_assigner_full`** in the same file — add a
   match arm that constructs your assigner. `new_assigner` and
   `new_assigner_with_threads` delegate to `_full`, so one arm covers all
   three entry points.

4. **Test.** Add unit tests alongside the implementation in `assigner.rs`,
   and a property test under `tests/integration/proptest/umi_adjacency.rs`
   (pattern: generate random UMI populations, assert invariants like "same
   UMI → same molecule ID" and "grouping is a partition").

5. **Update the grouping mdbook chapter** (under `docs/src/guide/`) to
   mention the new `--strategy` value if the assigner is user-visible on the
   CLI.
