# Deterministic `MoleculeId` numbering across runs

> Status: implemented
> See "History" at the bottom of this document for the original proposal
> branch and superseded PR.

## Problem

Two runs of `fgumi group --strategy paired --threads N` (and `fgumi dedup`) on
the same input BAM produce identical molecule groupings but **mismatched
integer `MI:Z` values** on roughly 10–15% of records (100% on some
position-group layouts). The non-determinism propagates into every downstream
stage that derives QNAMEs from the input `MI:Z` (`fgumi simplex`,
`fgumi duplex`, `fgumi codec` — all of which build the read name as
`<prefix>:<input_mi_base>`), and breaks `fgumi compare bams` correspondence
checks because two end-to-end runs of the pipeline disagree on which `MI:Z`
integer corresponds to which molecule.

The single-threaded fast path is already deterministic. The bug only appears
under `--threads N>1`.

### Root cause

The unified BAM pipeline's `Q5` (Process step output → Serialize step input)
is a `crossbeam_queue::ArrayQueue<(u64, Vec<P>)>`. Worker threads run
`process_fn` in parallel and push completions onto Q5 in **arrival order**
(i.e. completion order), not in pipeline-serial order. The Serialize step
then advances a shared `AtomicU64::fetch_add` (`next_mi_base`) inside
`serialize_fn`, observing batches in arrival order. Since each run sees a
different completion order, the cumulative `MoleculeId` offset assigned to
each batch varies across runs even though the per-batch contents are
identical.

`group` and `dedup` are the **only** BAM pipeline commands that mint new
sequential identifiers from shared cumulative state inside `serialize_fn`.
Every other command (`extract`, `clip`, `filter`, `correct`, `simplex`,
`duplex`, `codec`, `merge`) emits per-batch byte content that is fully
determined by that batch's input — they rely on the existing `write_reorder`
buffer at the Write step to reassemble batches into input order on disk, and
do not need any pre-Serialize ordering.

## Goals

1. Make `fgumi group` and `fgumi dedup` produce byte-identical outputs across
   runs of the same input under `--threads N>1` (`MI:Z` integers stable).
2. Do not change the byte-level output of any other command.
3. Do not penalize the throughput of commands that don't need the new
   ordering (no global single-threading of Process or Serialize).
4. Keep `serialize_fn` closures synchronization-free in `group` and `dedup`
   — all ordering logic should live in dedicated pipeline scaffolding.

## Non-goals

* Reordering records within the output file. The existing `write_reorder`
  already ensures records are written in input order; we only need
  consistent integer assignment per record.
* Eliminating `unsafe` or refactoring the unified pipeline beyond what is
  required for this fix.

## Design

### Overview

Add an **opt-in serial-ordered MI Assign stage** that rewrites per-template
`MoleculeId`s from local (per-position-group) values to global cumulative
values before any byte serialization touches them. Conceptually the stage
sits between Process and Serialize:

```text
Read → Decompress → Decode → Group → Process → [MI Assign]? → Serialize → Compress → Write
                              ↑seq    ↑par    ↑opt-in seq     ↑par         ↑par      ↑seq
```

In the shipped implementation the stage is **not** a standalone worker step
with its own intermediate `ArrayQueue`. Instead it is a serial-ordered drain
inside the Serialize step's input path, gated on `mi_assign_fn` being set
(see "Pipeline scaffolding" below for why).

* Group is already sequential (the grouper is stateful) and assigns one
  monotonically increasing `batch_serial` per Q4 push — this is `X` from the
  ordinal scheme.
* Process runs in parallel; each batch is `Vec<P>` and the index of an item
  inside that vector is `Y = idx_in_batch`. Item ordering within a batch is
  deterministic because `try_step_process` iterates the batch with a
  sequential `for` loop, so `(X, Y)` uniquely identifies each
  `ProcessedPositionGroup` in the run and is independent of completion
  timing.
* MI Assign drains Q5 into a serial-keyed reorder buffer and only releases
  the next expected `batch_serial`. Within a batch it walks position groups
  in `idx_in_batch` order, advances a single cumulative offset, and rewrites
  each template's `MoleculeId` from local to global before the batch
  proceeds into `serialize_fn`.
* Serialize stays parallel and reads each template's pre-rewritten MI
  verbatim — no shared state, no synchronization.

### Why opt-in

Only `group` and `dedup` need this stage. Wiring it in unconditionally would:

1. Force a single-threaded coordination point on every other BAM-pipeline
   command, costing throughput for no benefit.
2. Break the existing `test_pipeline_concurrency` invariants, which rely on
   the current Q5 semantics (FIFO with completion-order pop). The closed
   PR demonstrated this — its global Q5 swap had to add deadlock-avoidance
   logic just to keep those tests passing.

The stage is therefore selected at pipeline-construction time by the
command, not enabled by default.

### Hook API

Extend `PipelineFunctions<G, P>` with one new optional field:

```rust
pub struct PipelineFunctions<G: Send, P: Send> {
    pub process_fn: Box<dyn Fn(G) -> io::Result<P> + Send + Sync>,
    pub serialize_fn: Box<dyn Fn(P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync>,
    pub secondary_serialize_fn:
        Option<Box<dyn Fn(&P, &mut Vec<u8>) -> io::Result<u64> + Send + Sync>>,

    /// **NEW.** Optional serial-ordered pre-serialize hook. When set, the
    /// pipeline routes each Process-step batch through a single-threaded
    /// reorder stage that calls this hook with the batch's
    /// `(batch_serial, idx_in_batch)` pair for each item, in input order,
    /// before the batch is handed to `serialize_fn`. The hook is permitted
    /// to mutate the item.
    ///
    /// Used by `fgumi group` and `fgumi dedup` to rewrite per-template
    /// `MoleculeId`s from local to global before Serialize runs.
    pub mi_assign_fn: Option<Box<
        dyn Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync,
    >>,
}

#[derive(Clone, Copy, Debug)]
pub struct BatchOrdinal {
    /// Group-step serial of the enclosing batch (X).
    pub batch_serial: u64,
    /// Zero-based index of this item within the batch (Y).
    pub idx_in_batch: usize,
    /// Total number of items in the batch (so the hook can detect the last
    /// one and finalize per-batch state if needed).
    pub batch_len: usize,
}
```

* `process_fn` and `serialize_fn` keep their existing signatures.
* Commands that do not install `mi_assign_fn` see exactly the same pipeline
  as today.
* The hook is `Fn` (not `FnMut`) so the pipeline scaffolding can hold it
  behind an `Arc` and let multiple workers reach the MI Assign stage even
  though only one of them holds the reorder lock at a time.

### Pipeline scaffolding

The MI Assign stage is implemented as a serial-ordered drain that piggybacks
on the existing Serialize step rather than as a standalone worker step with
its own intermediate `ArrayQueue`. This keeps the change isolated to
`try_step_serialize`'s input path and avoids touching the scheduler /
deadlock detector with a new `PipelineStep` variant.

In `src/lib/unified_pipeline/bam.rs`:

```rust
enum MiAssignPopOutcome<P> {
    Ready { serial: u64, batch: Vec<P> },
    Stalled,        // reorder buffer waiting on an earlier serial
    Empty,          // Q5 truly empty
}

fn try_pop_mi_assigned<G, P>(
    state: &BamPipelineState<G, P>,
    hook: &(dyn Fn(BatchOrdinal, &mut P) -> io::Result<()> + Send + Sync),
) -> MiAssignPopOutcome<P> {
    // Called only when `mi_assign_fn` is installed. Holds the
    // `mi_assign_reorder` mutex for the duration of (drain Q5 → buffer
    // insert → next-serial pop → in-place hook walk):
    //
    //   1. Lock `state.output.mi_assign_reorder` (Mutex<ReorderBuffer<Vec<P>>>).
    //   2. Drain whatever is currently in Q5 (`state.output.processed`)
    //      into the reorder buffer, keyed by `batch_serial`.
    //   3. Try to pop the next expected `batch_serial` from the reorder
    //      buffer. If it isn't buffered yet, return `Stalled` (we're
    //      waiting on an earlier serial). If Q5 was also empty, return
    //      `Empty` so the caller can record a true Q5-empty stall.
    //   4. For each item in the popped batch, call
    //      `hook(BatchOrdinal { batch_serial, idx_in_batch, batch_len }, item)`
    //      in `idx_in_batch` order — the hook mutates each item in place.
    //   5. Return `Ready { serial, batch }`; the caller hands `batch` to
    //      `serialize_fn` exactly as if it had popped directly from Q5.
}
```

Wire-up actually shipped:

* `OutputPipelineQueues` gains one new field that is allocated
  unconditionally but only touched when `mi_assign_fn` is installed:
  * `mi_assign_reorder: Mutex<ReorderBuffer<Vec<P>>>` — the serial-keyed
    reorder buffer that Q5 drains into before the hook fires.
* `OutputPipelineQueues::are_queues_empty` accounts for the new in-flight
  buffer (a non-empty `mi_assign_reorder` means work is parked, not done).
* `try_step_serialize` checks `fns.mi_assign_fn`: if set, it pops via
  `try_pop_mi_assigned` and only records a Q5-empty stall on
  `MiAssignPopOutcome::Empty` (not on `Stalled`, where the reorder buffer
  is just waiting on an earlier serial). If `mi_assign_fn` is `None`,
  `try_step_serialize` keeps its existing behavior of popping directly
  from Q5.
* No new `PipelineStep` variant is added. The scheduler and deadlock
  detector continue to see the work as part of `PipelineStep::Serialize`,
  so existing `test_pipeline_concurrency` invariants are unchanged.
* Heap accounting stays charged to `processed_heap_bytes` while batches
  are parked in `mi_assign_reorder`, so `is_q5_memory_high()` keeps the
  upstream Process step backpressured.
* The single-threaded fast path (`run_bam_pipeline_single_threaded`) calls
  the hook inline via `run_mi_assign_single_threaded` between `process_fn`
  and `serialize_fn`; no reorder buffer is needed because batches arrive
  in serial order by construction.

The doc's earlier draft proposed a separate `try_step_mi_assign` worker
plus a `Q5_5: ArrayQueue<(u64, Vec<P>)>` intermediate queue and a
`PipelineStep::MiAssign` scheduler variant. The shipped implementation
folds all three into the serialize-time drain above; see "Discarded
alternatives" for why.

### `MoleculeId` rewrite

`ProcessedPositionGroup` already carries everything the hook needs:

* `templates: Vec<Template>` where `Template.mi: MoleculeId` holds local
  IDs (`Single(0..N-1)`, `PairedA(0..N-1)`, `PairedB(0..N-1)`).
* `distinct_mi_count: u64` — the number of distinct local IDs in the group.

The hook's body for `group`:

```rust
let base = mi_cumulative
    .fetch_update(
        std::sync::atomic::Ordering::Relaxed,
        std::sync::atomic::Ordering::Relaxed,
        |current| current.checked_add(processed.distinct_mi_count),
    )
    .map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "MoleculeId offset overflow: cumulative MI counter exceeded u64::MAX",
        )
    })?;
for template in &mut processed.templates {
    template.mi = template.mi.with_offset(base);
}
Ok(())
```

`fetch_update` + `checked_add` so wraparound is detected and surfaced
rather than silently reusing MI integers (which would defeat
`MoleculeId::with_offset`'s own overflow check on the per-template add).
The `AtomicU64` is safe here even though the hook closure is `Fn`, because
the MI Assign stage runs single-threaded — only one worker ever calls the
hook at a time. (We could just as easily use a `Cell<u64>` plus
`std::sync::Mutex`; using `Atomic` keeps the closure trait-bound simple.)

`dedup` is identical except `distinct_mi_count` is computed from
`processed.family_sizes.values().sum()` (matching its existing logic).

`MoleculeId::with_offset(offset: u64) -> MoleculeId` is a small new helper
on the existing enum:

```rust
impl MoleculeId {
    pub fn with_offset(self, offset: u64) -> Self {
        match self {
            MoleculeId::None => MoleculeId::None,
            MoleculeId::Single(id) => {
                MoleculeId::Single(id.checked_add(offset).expect("MoleculeId offset overflow"))
            }
            MoleculeId::PairedA(id) => {
                MoleculeId::PairedA(id.checked_add(offset).expect("MoleculeId offset overflow"))
            }
            MoleculeId::PairedB(id) => {
                MoleculeId::PairedB(id.checked_add(offset).expect("MoleculeId offset overflow"))
            }
        }
    }
}
```

After the hook runs, `serialize_fn` in `group` and `dedup` removes its old
`next_mi_base.fetch_add(...)` call entirely and just reads
`template.mi.write_with_offset(0, ...)` (or equivalent) — the rewrite has
already turned local IDs into global ones.

## Implementation plan

1. **MoleculeId helper.** Add `MoleculeId::with_offset` in
   `crates/fgumi-umi/src/lib.rs` plus unit tests covering all four variants.
   _(Smallest, most isolated change — landing first lets later steps keep
   their diffs focused.)_
2. **`BatchOrdinal` + `mi_assign_fn` field on `PipelineFunctions`.** Default
   `mi_assign_fn` to `None` everywhere; no behavior change.
3. **`run_bam_pipeline_single_threaded`.** When `mi_assign_fn` is set, call
   it inline between `process_fn` and `serialize_fn` via
   `run_mi_assign_single_threaded` with `idx_in_batch=0, batch_len=1`. Add
   a unit test asserting the hook fires exactly once per group with the
   right ordinal.
4. **Multi-threaded path: serialize-time MI Assign drain.** Add a single
   `mi_assign_reorder: Mutex<ReorderBuffer<Vec<P>>>` field on
   `OutputPipelineQueues`. Implement `try_pop_mi_assigned` returning
   `MiAssignPopOutcome::{Ready, Stalled, Empty}` and route
   `try_step_serialize`'s input pop through it when `mi_assign_fn` is
   installed. Update `are_queues_empty` to treat the reorder buffer as
   in-flight work. When `mi_assign_fn` is `None`, `try_step_serialize`
   keeps its existing direct `output.processed.pop()` (no behavior change
   for any other command). No new `PipelineStep` variant or intermediate
   `ArrayQueue` is needed.
5. **Switch `group` and `dedup` to the hook.** Each command:
   * stops doing `next_mi_base.fetch_add(...)` inside its `serialize_fn`,
   * installs an `mi_assign_fn` that captures an `Arc<AtomicU64>` cumulative
     offset and rewrites templates via `MoleculeId::with_offset`,
   * `serialize_fn` reads pre-rewritten MIs verbatim.
6. **Regression test.** `tests/integration/test_group_determinism.rs`
   already exists on this branch (six-run determinism check for
   `group --threads 1`, `group --threads 4`, `dedup --threads 4`). It
   reproduces the bug deterministically against `main` and must pass after
   step 5.

Each step is independently reviewable; steps 1–4 land plumbing without
changing behavior, step 5 makes the actual fix, step 6 locks it down.

## Test plan

* Unit: `MoleculeId::with_offset` per-variant tests.
* Unit: a tiny `try_pop_mi_assigned` test using a `PipelineFunctions`
  populated with mock `process_fn`/`serialize_fn`/`mi_assign_fn` that
  asserts:
  * the hook is invoked exactly once per `(batch_serial, idx_in_batch)`
    pair for every batch produced by Process,
  * the hook always observes `batch_serial`s in monotonically increasing
    order even when batches arrive out of order at Q5,
  * the batch reaches `serialize_fn` after the hook has mutated it
    (verify by capturing pre/post snapshots in the mock).
* Concurrency: extend `tests/integration/test_pipeline_concurrency.rs` with
  a stress test that runs `fgumi group --threads 8` on a large synthetic
  BAM and asserts the MI Assign stage neither deadlocks nor reorders the
  output relative to input.
* Integration: `tests/integration/test_group_determinism.rs` (already on
  this branch) — `group --threads 1`, `group --threads 4`,
  `dedup --threads 4`, six runs each, per-`(QNAME, flags)` `MI:Z` parity.
* Real-data smoke (out of CI): `fgumi zipper → sort → group --threads 8` on
  SRR6109255 (~1.7M records); two consecutive runs must produce
  byte-identical SAM bodies. Same for `fgumi duplex` on the resulting
  grouped BAM (must remain deterministic; this command is unchanged but is
  the user-visible symptom).
* Negative: `cargo ci-test` and `cargo ci-lint` clean. The full suite
  (2510 tests at the time of writing) must continue to pass with no
  changes to commands other than `group` and `dedup`.

## Discarded alternatives

* **In-place fix in `serialize_fn`** (closed PR #318): added
  `SerialOrderedArrayQueue` at Q5, `OrderedMiAllocator` (allocate / finalize
  with condvar), and a `serialize_context` thread-local. Correct and
  passing CI, but threaded synchronization through every Q5 consumer and
  required adding a deadlock-avoidance over-capacity branch to the queue
  to keep `test_pipeline_concurrency` happy. Replaced by the cleaner
  staged design above.
* **Mutex around the entire Serialize step.** Would let us keep the
  existing `AtomicU64::fetch_add` pattern, but globally single-threads
  serialization for every command, including `simplex`/`duplex`/`codec`
  whose serialize work is non-trivial.
* **Defer MI assignment to the Write step.** `write_reorder` is already
  serial-ordered, but by the time data reaches it the BAM bytes have
  already been compressed into BGZF blocks. Patching `MI:Z` integer
  strings post-compression would require decompress→patch→recompress for
  every batch.

## Open questions

* Whether `mi_assign_fn` should also receive a `last_in_pipeline: bool` to
  let the hook do per-pipeline cleanup (currently the closure can do that
  itself by capturing state and dropping). Lean: skip until needed.
* Whether the new step warrants its own scheduler priority knob (today's
  scheduler picks the busiest queue; MI Assign is single-threaded but
  cheap, so default round-robin should be fine). Lean: do not add a knob
  until profiling on real data shows starvation.

## History

* Original branch: `nh/consensus-qname-determinism`.
* Supersedes: closed PR fulcrumgenomics/fgumi#318 (in-place fix that
  threaded synchronization through every Q5 consumer; see "Discarded
  alternatives" above).
* Shipped implementation differs from the initial proposal in one
  notable way: the MI Assign stage is realized as a serial-ordered drain
  inside `try_step_serialize` (`try_pop_mi_assigned` +
  `OutputPipelineQueues::mi_assign_reorder`) rather than as a standalone
  `try_step_mi_assign` worker step with its own `Q5_5: ArrayQueue` and
  `PipelineStep::MiAssign` scheduler variant. The shipped form has
  identical semantics, smaller surface area, and avoids touching the
  scheduler / deadlock detector.
