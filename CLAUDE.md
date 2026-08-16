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
Install with `./scripts/install-hooks.sh`. Runs `cargo ci-fmt`, `cargo ci-lint`,
`cargo ci-tag-literals`, and `cargo ci-publish-order` before commits.

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

### Approved spare-capacity record read (fgumi-raw-bam)

- **`crates/fgumi-raw-bam/src/raw_bam_record.rs`** — `read_raw_record` reads a BAM record body into the spare capacity of the caller's reused `Vec<u8>` to skip the `resize(_, 0)` zero-fill (gigabytes of wasted zero-stores on large BAMs). SAFETY: `Read::read_exact` only writes to the destination; `spare_capacity_mut()[..block_size]` is valid uninitialized storage; `set_len(block_size)` runs only after `read_exact` returns `Ok`, so every byte covered by the new length is initialized. On error `set_len` is skipped and the Vec length stays 0.

### Approved SIMD intrinsics (fgumi-simd-fastq)

`crates/fgumi-simd-fastq/src/lexer.rs` classifies 64-byte FASTQ blocks through
hand-written NEON (aarch64) and AVX2 (`x86_64`) intrinsics, with a scalar
fallback. All `unsafe` in the crate is confined to this one file. It is approved
because the lexer runs once per 64 bytes of every FASTQ processed (billions of
blocks on production inputs) and the scalar fallback — retained and exercised in
tests as the reference implementation — is materially slower.

Calling a `#[target_feature]` function is `unsafe` because the caller must
guarantee the CPU supports that feature; the intrinsics themselves are then safe
within it. Two obligations therefore apply at every site, and both are recorded
in `// SAFETY:` comments:

- **Feature availability.** `lex_block` / `lex_block_full` dispatch on the
  architecture. On aarch64 NEON is architecturally mandatory, so the precondition
  is unconditional. On `x86_64` every AVX2 call is guarded by
  `is_x86_feature_detected!("avx2")` with a scalar `else` branch.
- **Load bounds.** The only memory accesses are vector loads from
  `block.as_ptr()`, where `block: &[u8; 64]`. NEON issues four 16-byte
  `vld1q_u8` loads at offsets 0/16/32/48; AVX2 issues two 32-byte
  `_mm256_loadu_si256` loads at offsets 0 and 32. Both cover exactly 64 bytes and
  no more. Both use unaligned load forms, so no alignment beyond `u8` is
  required.

Seven `#[allow(unsafe_code)]` sites: the two public dispatch functions
(`lex_block`, `lex_block_full`), the NEON helpers (`neon_movemask`,
`lex_block_newlines_neon`, `lex_block_neon`), and the AVX2 implementations
(`lex_block_newlines_avx2`, `lex_block_avx2`). The nested `neon_two_bits` /
`avx2_two_bits` helpers are covered by the lexically-enclosing allow.

Note that Miri cannot execute these intrinsics, so `miri.yml` does not cover
this crate. Correctness rests instead on the scalar-parity tests in `lexer.rs`,
which assert the SIMD and scalar lexers agree. Fixed inputs pin specific shapes
(`test_simd_matches_scalar_newlines`, `test_simd_full_matches_scalar_all_fields`,
`test_simd_full_matches_scalar_all_acgt`,
`test_simd_full_matches_scalar_every_position`), and `proptest` sweeps cover
arbitrary lane combinations (`simd_newlines_match_scalar_on_arbitrary_blocks`,
`simd_full_matches_scalar_on_arbitrary_blocks`,
`newline_fast_path_matches_full_classifier`,
`simd_matches_scalar_on_uniform_random_blocks`). Any change to the intrinsic
paths must keep these passing on both architectures — note that a dev machine
runs only one of the two, so the other is verified in CI.

The sweeps compare `two_bits` only on lanes where `is_acgt` is set. That is the
contract rather than a weakened assertion: `ENCODE_LUT` maps non-ACGT bytes to a
don't-care value, which the SIMD paths write through unconditionally while the
scalar path leaves those lanes zero.

Any new `unsafe` site must extend this list and explain why the safe
alternative is unacceptable. Do not introduce `unsafe` outside the crates
listed in this section.

## Benchmarking Notes

### Profiling fgumi: what works, in order

Reach for these in order. The first two answer most questions and need no new
code; only drop to sampling when you genuinely need per-function attribution.

**1. `tricorder` — core utilization, memory, and I/O.** The fastest way to learn
whether a command is CPU-bound, I/O-bound, or leaving cores idle:

```sh
tricorder --interval 0.25 --summary --out agg.tsv --trace trace.tsv -- fgumi <cmd> ...
# tricorder: wall=23.68s cpu=23.22s mean_load=98% max_rss=20.3MiB io_in=4.5MiB
```

`mean_load` is the number to read first: **98% means one core busy**, so on a
12-core box the command is single-threaded regardless of what `--threads` was
passed. `--trace` gives cores-over-time, which distinguishes "parallel
throughout" from "parallel phase followed by a serial tail". Note that `io_in`
collapses to near zero once the input is in page cache — a low `io_in` means the
run was CPU-bound *on that invocation*, not that the tool does little I/O.

**2. In-tree phase timers — the primary attribution tool.** `fgumi sort` already
carries `SortPhaseTimer` (`crates/fgumi-sort/src/external.rs`): per-phase `f64`
accumulators, a `time()` helper wrapping each span, and a `log_summary()` that
emits `=== Sort Phase Timing ===` at `info!` (read+decompress / in-memory sort /
spill write / consolidation / k-way merge / write output, each with a percentage).
Prior campaigns got full baseline phase attribution from this with **zero new
code** — just `RUST_LOG=info`.

Mirror that pattern when optimizing any other command: the phases you would name
in a design doc are exactly the phases worth timing, it works identically on
macOS and Linux, and it cannot be defeated by inlining or missing symbols the way
a sampling profiler can. Prefer adding a phase timer over fighting a profiler.

**3. `perf record` / `perf stat` on Linux (EC2)** when per-function or
per-instruction attribution is genuinely required. This is the reliable path for
function-level work; see the EC2 host recipe in the user-level CLAUDE.md.

#### macOS sampling: pitfalls that produce confidently wrong profiles

`samply` works on macOS but has three traps with this codebase. All three
silently yield a plausible-looking, wrong profile:

- **The release binary has no symbols.** `[profile.release]` builds get
  `-C strip=debuginfo`, so everything symbolicates to raw addresses. Use the
  existing `[profile.profiling]` (inherits release, `debug = "line-tables-only"`),
  or `CARGO_PROFILE_PROFILING_DEBUG=2` for full DWARF.
- **Blocked threads are counted as CPU.** samply samples every thread on a timer,
  including ones parked on a channel or condvar. Counting raw samples credited an
  idle BGZF worker pool the same CPU as the working thread and made an 18 s run
  look like 285 s across 21 threads. Weight samples by `threadCPUDelta`
  (microseconds of real CPU between samples) instead of counting them.
- **`atos` mis-resolves addresses from other images.** Feeding every sampled
  address to `atos -o target/profiling/fgumi` resolves `libsystem_kernel` /
  `libsystem_pthread` / `libsystem_malloc` addresses against the fgumi image,
  landing on whatever symbol happens to sit at that offset. This is where
  nonsense like "38% of CPU in `clap_builder::error::Error::print`" on a
  *successful* run comes from. Map each frame to its library first (Gecko
  `funcTable.resource` → `resourceTable` → `libs`) and only symbolicate frames
  belonging to the fgumi image. Under fat LTO, local symbols are gone, so even
  in-image resolution can land on a neighbouring function — treat any
  suspiciously hot dependency symbol as misattribution until corroborated.

Corroborate a sampling result against `tricorder` or a phase timer before acting
on it. In the `compare bams` work, the sampling profile and `tricorder` agreed
that BGZF decode dominates (~61% self time; `mean_load=98%`), which is what made
the conclusion trustworthy — not the profile alone.

### `fgumi compare bams` performance (issue #686, 2026-08-01)

Baseline on M2 Max (12 core), 20.1M records/file, ~2 GB BAMs, page-cached,
`--threads 8`:

| mode | wall | cpu | mean_load | cores busy | peak OS threads |
| --- | --- | --- | --- | --- | --- |
| `--command sort` (coordinate) | 23.7s | 23.2s | 98% | 0.98 / 12 | **1** |
| `--command group` | 24.7s | 24.3s | 98% | 0.98 / 12 | **1** |
| `--command simplex` (content) | 18.3s | 45.1s | 246% | 2.46 / 12 | 21 |

`sort` and `group` ignored `--threads` entirely — the tricorder `--trace`
`n_threads` column stayed at **1** for the whole run, so the flag created no
threads at all. Sampling put ~55% of that single thread in
`libdeflate_deflate_decompress_ex` and ~6% in
`fgumi_bgzf::reader::verify_decompression` (the per-block CRC32): ~61% in BGZF
decode, serialized. `content` reached only 2.46 cores and additionally traversed
each input **twice**, because `verify_records_in_order` ran over both paths before
`positional_compare` reopened and streamed them again.

After the work on #686:

| mode | before | after | note |
| --- | --- | --- | --- |
| `--command sort` `-t 8` | 23.44s | **4.14s** | 5.7x; scales `-t 1` 8.92s → `-t 4` 5.58s → `-t 8` 4.14s |
| `--command group` `-t 8` | 24.36s | ~10s | reader concurrency only; its engine was not otherwise changed |
| `--command simplex` (content) | 47.1 CPU-s | **27.4 CPU-s** | -42% CPU |

`compare --command sort` is now 2.0x *faster* than the `fgumi sort` it validates
(4.14s vs 8.20s), having started 2.7x slower.

Two findings worth keeping:

- **Report content-mode results in CPU-seconds, not wall.** This host's wall time
  is unreliable under interactive load — identical `-t 8` runs measured 3.56s and
  16.56s — while total CPU stayed within 24–27s across the same sweep. Total CPU is
  the metric the content conclusions rest on.
- **Parallel BGZF decode costs ~40% more total CPU** than single-threaded decode
  for the same work (19.3 CPU-s at `-t 1` vs ~27 CPU-s at `-t 4`+), spent in
  `crossbeam_channel::recv` and thread parking. It buys wall-clock on an idle host
  and costs throughput on a busy one. Raising the read-ahead batch from 256 to 4096
  records did **not** reduce it (26.5 vs 27.0 CPU-s) but tripled peak RSS, so that
  time is the consumer waiting on decode, not per-handoff overhead.

### samtools sort orders
`samtools sort` supports `--template-coordinate` for template-coordinate sorting. When benchmarking sort orders:
- Coordinate: `samtools sort -@ 4 -m 50M input.bam -o output.bam`
- Queryname: `samtools sort -n -@ 4 -m 50M input.bam -o output.bam`
- Template-coordinate: `samtools sort --template-coordinate -@ 4 -m 50M input.bam -o output.bam`
