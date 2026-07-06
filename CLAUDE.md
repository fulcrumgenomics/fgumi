# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fgumi is a high-performance Rust CLI tool for UMI (Unique Molecular Identifier) processing in sequencing data. It provides extraction, grouping, correction, and consensus calling functionality.

## Build and Test Commands

```bash
# Build
cargo build --release

# Run all tests (recommended - uses nextest)
cargo ci-test

# Run tests with test-utils feature
cargo t

# Check formatting
cargo ci-fmt

# Run linting
cargo ci-lint

# Run all CI checks
cargo ci-test && cargo ci-fmt && cargo ci-lint

# Run a single test
cargo nextest run test_name

# Run benchmarks
cargo bench
```

## Features

- `consensus` *(default-on)* - All consensus calling: the `simplex`, `duplex`,
  `codec`, `runall`, `simplex-metrics`, and `duplex-metrics` subcommands. This is
  a single umbrella toggle; the granular `simplex`/`duplex`/`codec` features live
  one level down in the `fgumi-consensus` crate for embedders that want a reduced
  library build. Building with `--no-default-features` cleanly drops all consensus
  code (the root crate still compiles); pass `--features consensus` to add it back.
- `simulate` - Enable the simulate command for test data generation (implies `consensus`)
- `compare` - Enable compare subcommand (developer tools)
- `profile-adjacency` - Enable profiling output for adjacency UMI assigner

Build with features: `cargo build --release --features compare,simulate`

The reduced-feature build (`cargo check -p fgumi --no-default-features`) is
exercised in CI by the `fgumi-feature-check` job so always-compiled chain/runall
code can never silently reference a feature-gated consensus module again.

## Code Architecture

### Entry Points
- `src/main.rs` - CLI entry point using clap with enum_dispatch for command routing
- `src/lib/mod.rs` - Library crate (`fgumi_lib`) with all core functionality

### Library Modules (`src/lib/`)

**Core UMI Processing:**
- `umi/` - UMI assignment strategies: identity, edit-distance, adjacency, paired
- `consensus/` - Consensus callers: simplex, duplex, codec
- `grouper.rs` - Read grouping by UMI

**I/O:**
- `bam_io/` - BAM file reading/writing helpers
- `bgzf_reader/`, `bgzf_writer/` - BGZF compression I/O
- `sam/` - SAM/BAM alignment utilities and tag manipulation

**Utilities:**
- `bitenc.rs` - 2-bit DNA encoding (4 bases per byte, 32 per u64)
- `metrics/` - QC metrics for duplex, consensus, grouping
- `validation.rs` - Input file and parameter validation
- `progress.rs` - Progress tracking
- `reference.rs` - Reference genome/dict file handling

**Specialized:**
- `clipper.rs` - Overlapping read pair clipping
- `template.rs` - Template-based read grouping
- `vendored/bam_codec/` - Isolated vendored BAM compression codec

### Commands (`src/commands/`)

Commands implement the `Command` trait dispatched via enum:

1. **Extraction:** `extract` - Extract UMIs from FASTQ to unmapped BAM
2. **Alignment Support:** `fastq`, `zipper`, `sort`
3. **UMI Operations:** `correct`, `group`, `dedup`
4. **Consensus:** `simplex`, `duplex`, `codec`
5. **Post-consensus:** `filter`, `clip`, `duplex-metrics`, `review`
6. **Utilities:** `downsample`, `compare` (feature-gated), `simulate` (feature-gated)

### Key Design Patterns

- **Enum Dispatch:** Commands use `enum_dispatch` for trait-based routing
- **Streaming I/O:** Large files processed without full loading
- **Thread Pooling:** Work-stealing with per-command thread optimization
- **2-bit Encoding:** DNA bases packed efficiently for fast operations
- **Typed-Step Pipeline Framework (`fgumi_pipeline_core`, re-exported as `pipeline::core`):** All multi-threaded commands run through the typed-step framework (`Pipeline::builder().chain(step1).chain(step2).…build().run(...)`). The execution engine lives in the `fgumi-pipeline-core` crate (`crates/fgumi-pipeline-core/`) and is re-exported as `crate::pipeline::core`; the `chains` and `steps` layers remain in the `fgumi` crate. Steps declare `StepKind` (Serial/Parallel/Exclusive), input/output handle types, and an `on_input_drained` callback for end-of-stream cleanup. The `--use-new-pipeline` flag and the legacy `run_bam_pipeline_from_reader{,_with_mi_assign}` drivers were removed in Phase 1 of issue #330; the typed-step framework is the execution mechanism for sort, group, simplex, duplex, codec, correct, zipper, clip, filter, and dedup. New commands should follow the typed-step pattern.
- **Chain-builder façade (`pipeline::chains`):** Phase 2 introduced a single declarative chain-construction entry point — `chains::build_for(spec) -> Result<BuiltPipeline>`. Phase 3 introduced the stage-by-stage `ChainBuilder`: `build_for` validates the spec, then constructs a `ChainBuilder`, calls `add_source()`, walks `spec.stages` calling `chain.add_stage(stage, position)` for each, and finishes with `add_sink()` + `build()`. Each `add_<stage>` method (one per `Stage` variant — `add_dedup`, `add_filter`, `add_clip`, `add_sort`, `add_group`, `add_simplex`, `add_duplex`, `add_codec`, `add_correct`, `add_zipper`, `add_align`, `add_extract`) reads `spec.stage_opts.<stage>` (already validated present), pushes the canonical step sequence via factories in `chains::commands::<command>::build_*_step`, and registers a `<Command>FinalizeHook`. `StagePosition::{Terminal, Intermediate}` gates the serialize step so intermediate stages leave their typed output for the next `add_<stage>`. Shared config wiring (threads, deadlock_timeout, queue_memory, pipeline stats) lives in `chains::build_helpers::build_pipeline_config_for_chain`. The chain-level `StageTimingFinalizeHook` is registered by `ChainBuilder::build()`; `PipelineStatsFinalizeHook` (gated on `--pipeline-stats`) is the next hook in order. `BuiltPipeline::run()` executes the pipeline then drains hooks in registration order. **New commands MUST add: one `Stage` variant, one bag slot, one validator entry, one `add_<stage>` method, and one factory per per-stage step. Do NOT construct chains inline in command `execute` methods.** `runall::execute` routes every chain shape — including chains with an intermediate Sort — through `build_for`; `ChainBuilder::add_sort` performs the intermediate sort in-pipeline (streaming, not file-to-file), so there is no file-to-file fallback dispatcher in `commands::runall`.
- **Multi-Options Macro (`crates/fgumi-cli-macros`):** Per-stage tuning options used by both a standalone command (`fgumi sort`, `fgumi group`, …) and the fused `fgumi runall` command live in a single `<Stage>Options` struct annotated with `#[multi_options("stage", "Help Heading")]`. The standalone command flattens `<Stage>Options` directly so its CLI surface is unchanged (`--max-memory`, `--strategy`, …). The proc-macro generates a sibling `Multi<Stage>Options` struct that runall flattens, exposing the same fields as prefixed `--<stage>::<flag>` flags (`--sort::max-memory`, `--group::strategy`, …) grouped under a `--help` heading. The Multi struct carries a `validate(self) -> Result<<Stage>Options>` method that runall calls before executing a stage; required-without-default fields (e.g. `--group::strategy`) become `Option<T>` on the Multi side and the validator surfaces a clear "required when `<stage>` is selected" error if missing. The convention is established in `commands::{sort,group,duplex,codec}` and is the way to expose new per-stage options on runall going forward.

## Development Practices

### Pre-commit Hooks
Install with `./scripts/install-hooks.sh`. Runs `cargo ci-fmt` and `cargo ci-lint` before commits.

To include tests: `FGUMI_PRECOMMIT_TEST=1 git commit`

### Commit Messages
Use conventional commits format: `<type>[(scope)][!]: <description>`

Types: `feat`, `fix`, `build`, `chore`, `ci`, `config`, `docs`, `example`, `perf`, `refactor`, `style`, `test`

### Branch Naming
`<issue-number>/<user>/<type>-<description>` (e.g., `42/jdidion/fix-fibonacci-calculation`)

### Testing
- Unit tests alongside source in `src/`
- Integration tests in `tests/integration/`
- Uses `rstest` for parameterized tests, `proptest` for property-based testing
- `criterion` for benchmarks in `benches/`

## Rust Version

Minimum: 1.87.0 (uses edition 2024)

## Allocator

Uses `mimalloc` as global allocator for performance.

## Unsafe Code

`#![deny(unsafe_code)]` is set at the crate root of every workspace member. Targeted
`#[allow(unsafe_code)]` blocks are permitted only at the documented sites listed below;
any new `unsafe` block requires updating this section with a written justification.

### Approved non-stdlib FFI exceptions

The following external FFI calls are approved because they back core infrastructure:

- **`libmimalloc_sys`** (`crates/fgumi-sort/src/memory_probe.rs`) — mimalloc is the configured
  global allocator. The FFI wrappers (`force_mi_collect`, `process_rss_bytes`, and the
  `memory-debug`-gated `print_mi_stats` calling `mi_stats_print_out`) are isolated to
  a single `#[allow(unsafe_code)]` sub-module. mimalloc synchronizes these calls
  internally, so they are safe to invoke concurrently.
- **`mach2`** (`crates/fgumi-sort/src/memory_probe.rs`, macOS only) — `task_info(TASK_VM_INFO)`
  is the only way to read `phys_footprint` (the RSS metric mimalloc reports accurately).
  Isolated to the same sub-module as the mimalloc FFI.

### Approved hot-path unsafe (sort engine)

The radix-sort, in-memory record buffer, and natural-order comparator hot paths in
`fgumi-sort` use `unsafe` for performance. Each site is `#[allow(unsafe_code)]` and
covered by a `// SAFETY:` comment that documents the invariant. These are approved
because the hot path runs once per record (BAM workloads are millions to billions
of records) and a safe rewrite measurably regresses sort throughput.

- **`crates/fgumi-sort/src/inline.rs`** — three `#[allow(unsafe_code)]` regions:
  the `radix_sort_record_refs_with_max` (coordinate; `radix_sort_record_refs`
  delegates to it after scanning for the max) and `radix_sort_template_refs` /
  `radix_sort_template_field` (template) LSD radix sorts use `Vec::set_len` to
  skip per-element initialization on the auxiliary scratch buffer, plus
  raw-pointer slice swaps to avoid double-borrow restrictions across the
  source/destination ping-pong. The template path is now generic over the
  `TemplateLaneKey` trait (`TemplateKey24`/`32`/`40`), so the buffers are
  `Vec<TemplateRecordRef<K>>`; the field scatter was extracted into
  `radix_sort_template_field` to keep the per-field getter type fixed across the
  generic recursion. The unsafe invariant is unchanged by going generic. SAFETY
  relies on (a) the buffer being written exactly once per pass before being read,
  and (b) the pointers always referring to disjoint, properly-aligned
  `Vec<RecordRef>` / `Vec<TemplateRecordRef<K>>` storage (`TemplateRecordRef<K>`
  is `Copy`/`Pod` for any `K: TemplateLaneKey`).
- **`crates/fgumi-sort/src/keys.rs`** — `RawQuerynameKey::cmp` calls
  `natural_compare_nul` (defined in `fgumi-raw-bam`) over raw `*const u8`
  pointers. SAFETY: the names that back these pointers are always
  null-terminated by construction (see `extract_queryname_key` and
  `RawQuerynameKey::new`).
- **`crates/fgumi-sort/src/radix.rs`** — internal radix-sort helpers; see file
  comments for the `SAFETY:` invariants.
- **`crates/fgumi-sort/src/ref_sort.rs`** — one `#[allow(unsafe_code)]` site in
  `sort_coordinate_refs`: a pointer cast reinterpreting `&mut [RecordRef]` as
  `&mut [CoordSortRef]` to feed the parallel `voracious_mt_sort` radix on the
  large-input / multi-thread coordinate path. SAFETY: `CoordSortRef` is
  `#[repr(transparent)]` over `RecordRef`, so the two have identical size,
  alignment, and layout; the pointer cast (clippy rejects a ref-to-ref
  `transmute`) is sound. This is approved for the same reason as the other sort
  hot paths — it runs once per coordinate sort over millions of record refs, and
  the parallel radix is measurably (~4×) faster than the safe single-threaded
  fallback; both produce byte-identical output (verified by the ref-sort oracle
  proptests).
- **`crates/fgumi-sort/src/segmented_buf.rs`** — two `#[allow(unsafe_code)]`
  sites backing the arena record buffer (a segmented, append-only `Vec<u8>` the
  sort engine decompresses/frames records into without per-record allocation):
  - `grow_uninit` (≈line 207) — grows a segment's live `len` over
    reserved-but-uninitialized bytes via `Vec::set_len` (after `reserve`), to skip
    zero-filling the slot on the ingest hot path. SAFETY: after `reserve(additional)`,
    `capacity() >= new_len`, so `set_len(new_len)` only extends `len` over
    already-allocated bytes; `u8` has no drop glue and no validity invariant, so
    growing `len` over uninitialized bytes is not itself UB — the documented contract
    requires the caller to fully write the grown `[old_len, new_len)` region (via
    `slice_mut`) before any read, and a read-before-write is the only UB, which the
    contract forbids. Pointer stability is a *separate* invariant: this `reserve`
    MAY reallocate and move the segment, so no `slice_mut` borrow may be live across
    a `grow_uninit` — the borrow checker enforces this in the single-threaded case
    (`&mut self`), and the concurrent-inflate path first calls
    `reserve_full_capacity` once per fresh segment so the per-slot `reserve` is a
    no-op (it does NOT rely on `reserve` alone for stability).
  - `slice_mut` (≈line 267) — synthesizes a `&mut [u8]` slot from a shared `&self`
    (so disjoint slots can be written concurrently by different workers). SAFETY:
    the asserted `seg_offset + len <= seg.len()` keeps the range inside the
    segment's live region (in-bounds pointer + length); the `&mut` is sound only
    under the caller's disjointness contract — each call's range must not overlap
    any other concurrently-borrowed range — so the produced `&mut` aliases no
    other `&`/`&mut`. Approved for the same reason as the other sort hot paths:
    ingest runs once per record over millions-to-billions of records, and a safe
    rewrite (zero-fill on grow, or an owning `Vec` per slot) measurably regresses
    sort throughput.

### Approved natural-order comparator (fgumi-raw-bam)

`natural_compare` and `natural_compare_nul` (the samtools-compatible queryname
comparator) live in `crates/fgumi-raw-bam/src/sort.rs` because they operate on
raw BAM read-name bytes; both are called from `fgumi-sort` via
`RawQuerynameKey::cmp`. They use `unsafe` for the same reason as the sort hot
path: the comparator runs once per sort-key comparison, and the safe form
(re-bounds-checking every byte / re-validating the null terminator) measurably
regresses `samtools sort -n`–style throughput.

- **`crates/fgumi-raw-bam/src/sort.rs`** — four `#[allow(unsafe_code)]` sites:
  - `natural_compare` (line ~80) — `get_unchecked` over `&[u8]` in the digit-run
    hot loop. SAFETY: indices are bounded by the loop invariants `pa < alen` /
    `pb < blen`, asserted in `debug_assert!` for the `at` helper.
  - `natural_compare_nul` (line ~180) — raw `*const u8` walk that mirrors
    samtools' `strnum_cmp`. SAFETY: caller guarantees both pointers are
    null-terminated; `RawQuerynameKey::new` enforces this for every production
    call site.
  - `compare_nul` test helper and the `proptest` agreement test (lines ~273 and
    ~300) — push an explicit NUL into a `Vec<u8>` then take `as_ptr()`.
    SAFETY: the buffers are `to_vec()` + push, so the pointer is valid and
    null-terminated for the call's lifetime.

### Approved typed-step DSL hot path (issue-330 framework)

`TypedStep<S>::resolve_input` / `resolve_outputs` in
`crates/fgumi-pipeline-core/src/erased.rs` cache the typed downcast of
each step's `ctx.input` / `ctx.outputs` after the first dispatch.
`Any::downcast_ref` accounts for ~1.6% of CPU on the dispatch hot path
(profile-measured on a 4-thread CODEC 8M group benchmark), and the
boxes are owned by `ChainContexts` (an `Arc` held alive for the entire
`Pipeline::run` call) which guarantees they outlive every
`TypedStep<S>` instance. Caching is sound because every dispatch passes
the same box reference for a given `step_idx` (see
`run_worker_loop`'s `contexts.inputs[step_idx.0].as_ref()`).

- **`crates/fgumi-pipeline-core/src/erased.rs`** — two
  `#[allow(unsafe_code)]` sites:
  - `TypedStep::resolve_input` (≈line 152) — `mem::transmute` between
    `&'static BranchInputHandle<S::Input>` (cache slot type) and
    `&'a BranchInputHandle<S::Input>` (dispatch lifetime). SAFETY: the
    cached box outlives every `TypedStep<S>` (point 2 in the type-level
    doc); dispatches always pass the same box for the same `step_idx`
    (point 3). The cache is per-`TypedStep<S>` instance — `Parallel`
    clones get fresh caches via `clone_boxed`'s `TypedStep::new` call.
  - `TypedStep::resolve_outputs` (≈line 180) — same pattern, same
    safety argument applied to the typed `OutputHandles<S::Outputs>`
    view.

Any new `unsafe` site must extend this list and explain why the safe
alternative is unacceptable. Do not introduce `unsafe` outside the crates
listed in this section.

## Benchmarking Notes

### samtools sort orders
`samtools sort` supports `--template-coordinate` for template-coordinate sorting. When benchmarking sort orders:
- Coordinate: `samtools sort -@ 4 -m 50M input.bam -o output.bam`
- Queryname: `samtools sort -n -@ 4 -m 50M input.bam -o output.bam`
- Template-coordinate: `samtools sort --template-coordinate -@ 4 -m 50M input.bam -o output.bam`
