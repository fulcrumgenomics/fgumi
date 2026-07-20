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

- `compare` - Enable compare subcommand (developer tools)
- `simulate` - Enable simulate command for test data generation
- `profile-adjacency` - Enable profiling output for adjacency UMI assigner

Build with features: `cargo build --release --features compare,simulate`

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

Minimum supported and pinned build version: 1.93.0 (edition 2024).

fgumi supports exactly the Rust it is developed and tested on: `rust-version`
in `Cargo.toml` (the published MSRV) is kept identical to the `channel` in
`rust-toolchain.toml`, and CI's `msrv-lockstep` job fails if they drift. When
raising the toolchain, bump `rust-version` in the same change and adopt any
idioms the newer compiler unlocks (e.g. let-chains) rather than suppressing the
clippy lints that flag them.

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
- **`libc::umask`** (`crates/fgumi-sort/src/external.rs`, Unix only) — the merge command writes
  its output through an atomic temp+persist (`MergeOutputTarget`) so a rejected/failed merge
  never leaves a partial file. A `NamedTempFile` is created `0600`, so the persisted output must
  be re-stamped with the mode a plain `File::create` would have produced. For a new file that is
  `0o666 & !umask`, and reading the process `umask` requires the libc binding (POSIX `umask(2)`
  only *sets* the mask, returning the previous value, so `process_umask` sets it to `0` and
  restores it in the next call). The call is isolated to a single `#[allow(unsafe_code)]` site
  with a `// SAFETY:` note; it has no preconditions and cannot fail. The read-modify-restore is
  not atomic, so a process-global `Mutex` (`UMASK_LOCK`) serializes concurrent probes — required
  because `fgumi_lib` is a library and callers may run merges concurrently in one process;
  without the lock two interleaved probes could leave the process mask permanently `0`.

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

Any new `unsafe` site must extend this list and explain why the safe
alternative is unacceptable. Do not introduce `unsafe` outside the crates
listed in this section.

## Benchmarking Notes

### samtools sort orders
`samtools sort` supports `--template-coordinate` for template-coordinate sorting. When benchmarking sort orders:
- Coordinate: `samtools sort -@ 4 -m 50M input.bam -o output.bam`
- Queryname: `samtools sort -n -@ 4 -m 50M input.bam -o output.bam`
- Template-coordinate: `samtools sort --template-coordinate -@ 4 -m 50M input.bam -o output.bam`
