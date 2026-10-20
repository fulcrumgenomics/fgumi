# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [0.8.0] - 2026-10-18

### Bug Fixes

- [**breaking**] Enforce --max-temp-files on the arena spill path ([#993](https://github.com/fulcrumgenomics/fgumi/pull/993))
  - fgumi-pipeline-io's SpillBlockEvent::Block gains a required key_kind field, and SpillWrite now emits SpillReady for its surviving runs at AllAnnounced (with slot_count set to their number) instead of as each run closes. RunMergerDyn::step takes a byte budget.
- Accept bwa's mid-pair -K split; emit BAM from the bwa-mem3 preset; keep mimalloc from purging ([#988](https://github.com/fulcrumgenomics/fgumi/pull/988))
- Bound the mimalloc purge delay instead of never purging ([#1017](https://github.com/fulcrumgenomics/fgumi/pull/1017))
- Sync and check close of output files ([#1047](https://github.com/fulcrumgenomics/fgumi/pull/1047))
- Rename conflicting @RG/@PG IDs per input and rewrite read tags ([#1050](https://github.com/fulcrumgenomics/fgumi/pull/1050))
- [**breaking**] Honour --sort-threads/--merge-threads as phase caps within --threads ([#1048](https://github.com/fulcrumgenomics/fgumi/pull/1048))
  - the sort's memory budget now follows the effective sort-phase count. `--memory-per-thread` multiplies `--max-memory` by `min(--sort-threads, pool)` (by `--threads` when `--sort-threads` is unset) instead of `max(--threads, --sort-threads)`. For standalone `fgumi sort` the pool is `--threads`: a `--sort-threads` below it shrinks the record buffer proportionally and spills sooner, one above it no longer enlarges the buffer (it warns and is clamped), and `--sort-threads 0` gives a 1x budget. In runall the pool is the chain's: an explicit `--sort::sort-threads` up to it (at least 4 with a zipper or subprocess aligner) sizes the budget, and unset keeps it at `--threads`. `--sort-threads` / `--merge-threads` above the pool now warn and are clamped instead of being honoured. A run that relied on `--sort-threads` above `--threads` to enlarge the buffer should raise `--threads` or `--max-memory` instead.

### Documentation

- Remove in-tree design docs and gitignore docs/design/ ([#934](https://github.com/fulcrumgenomics/fgumi/pull/934))

### Features

- Typed-step pipeline foundation ([#870](https://github.com/fulcrumgenomics/fgumi/pull/870))
- Declarative chain-builder layer + group --threads pilot ([#872](https://github.com/fulcrumgenomics/fgumi/pull/872))
- Restore --read-streams on the arena chain via a concurrent scatter reader ([#889](https://github.com/fulcrumgenomics/fgumi/pull/889))
- Give stage option structs a clap::Args surface for runall (PR A) ([#907](https://github.com/fulcrumgenomics/fgumi/pull/907))
- Add the bam-roundtrip subcommand ([#936](https://github.com/fulcrumgenomics/fgumi/pull/936))
- Thread-utilization and per-step bandwidth telemetry for the chain-builder engine ([#951](https://github.com/fulcrumgenomics/fgumi/pull/951))
- Align in process with bwa-mem3 (--aligner::preset bwa-mem3-inproc) ([#990](https://github.com/fulcrumgenomics/fgumi/pull/990))
- [**breaking**] Template key context/batch extraction and reader-excluding pool placement ([#1055](https://github.com/fulcrumgenomics/fgumi/pull/1055))
  - `ErasedStep` gains a required `pool_placement()` method with no default, so an external implementor must add it (the in-tree `TypedStep`, `TypedStep2` and `TypedStepK` adapters forward `Step::pool_placement`, which defaults to `PoolPlacement::AllWorkers`).

### Performance

- Use fixed-seed ahash on hot per-item maps to remove RandomState global-counter contention ([#865](https://github.com/fulcrumgenomics/fgumi/pull/865))
- Extend runs in the arena spill path (t16 parity) ([#886](https://github.com/fulcrumgenomics/fgumi/pull/886))
- Parallel gather for the in-memory sort fast path ([#966](https://github.com/fulcrumgenomics/fgumi/pull/966))
- Stop fdatasync-ing output files on close ([#1054](https://github.com/fulcrumgenomics/fgumi/pull/1054))
- [**breaking**] Wake the merge when its awaited block lands; request a pool worker while it is starved ([#1056](https://github.com/fulcrumgenomics/fgumi/pull/1056))
  - `StepCtx`, `StepCtx2`, `StepCtxK` and `ErasedStepCtx` gain the public field `pool`, so a struct literal of any of them outside this crate no longer compiles. `StepStats` and `StepStatsSnapshot` gain the public fields `pool_requests_woken`, `pool_requests_all_awake` and `pool_requests_pending`.  * perf(sort): wake the merge when its awaited block lands (MergeDemand)  The k-way merge sleeps on its driver's park timer whenever the next record's block is not decompressed yet, and nothing ever woke it early: on the template-coordinate t16 cell it parked ~1.23 M times (~195 s). A per-sort `MergeDemand` (fgumi-sort) now records which spill file the merge waits on and the consumer thread; `MergeDemand::await_slot` is the one consumer sequence (declare the slot awaited, then re-check it under the slot mutex), and `SortSpillDecompress` calls `notify_delivered` after every delivery, EOF or failure, which unparks the merge only when that file is the awaited one.
- [**breaking**] Positional reads and a demand-ordered, budgeted spill supply for the merge ([#1057](https://github.com/fulcrumgenomics/fgumi/pull/1057))
  - the `fgumi_bam_io::scatter_reader` module is now `fgumi_bam_io::pread`; its items keep their names.  * feat(bam-io)!: positional-read primitives, SliceBufferPool, ReadStreamsPolicy ratchet and the PreadSlices step (no caller)  The positional-read primitives the sort's input and spill paths are built on, with no caller yet.  - `fgumi_bam_io::pread`: `PositionalSource` and `read_at_exact`;   `ReadStreamsPolicy`, a monotone ratchet that starts `--read-streams auto` at   one stream and doubles at a window boundary when the device starved the   framer for more than a quarter of the window's wall (`Fixed(n)` pins it);   `SliceBufferPool`/`SliceLease`, a refcounted read-only buffer returned to   its pool on the last drop, reused only for a slice of exactly its size so a   lease's capacity is its length (the gauges charge capacity and count idle   buffers: `resident_bytes`, `idle_bytes`, `peak_held_bytes`);   `RawFrame::{Borrowed, Owned}` with a checked `RawFrame::borrowed`. - `fgumi_bgzf::BgzfSliceFramer` cuts BGZF blocks out of arbitrary slices as   ranges, carrying a block that straddles slices (including a cut inside the   18-byte header) in one allocation of the block's size; `finish()` on a   partial block is `UnexpectedEof`. - `fgumi_pipeline_io::pread::PreadSlices`, one `Parallel` step type for both   read paths, in a module neither `source` nor `sort` owns: it checks for a   request before taking its phase permit, reads each slice into a pooled   buffer and fails closed on a short read. A spill read carries its   read-ahead class on `ReadTarget::Spill { slot, class }`, so an input request   cannot carry a spill class. A queued slice is charged its buffer capacity.  `ScatterReader`, its per-read scoped threads and the stream-count probe are deleted, and with them the openers' `--read-streams` handling (`decide_reader`, `warn_read_streams_unavailable`, `is_scatter_capable`,
- [**breaking**] Merge keeps partial batches across stalls and emits records borrowed in place ([#1059](https://github.com/fulcrumgenomics/fgumi/pull/1059))
  - fgumi-sort removes MergeDemandStats::record_partial_flush and the MergeDemandSnapshot::partial_flushes field; MergeDemandSnapshot:: log_lines returns five lines instead of six.  * perf(sort)!: merge emits records borrowed in place — one copy per record  The slot-backed merge copied every record twice: from the decompressed block into a per-leaf scratch buffer, then from the scratch into the output batch. The shared in-memory residual was copied the same way, out of its arena.  A merge source now reports where its current record lives (RecLoc): a range of its current block, an index into the shared residual chunk, or the leaf's scratch, which holds only a record that spans blocks or an owned residual record (swapped in, as before). MergeStep::Produced borrows the winner's bytes from there. No unsafe: the borrow is sound by the existing deferred-refill protocol. A source's block is replaced only by that source's own read (its refill or priming), which runs only on a later try_step, and Produced borrows &mut self, so the caller has consumed the bytes before any read can replace them.  --sort-stats regains a `Merge output:` line, printed only when the merge ran through MergeDriver: `records borrowed B / assembled A`, which sum to the records merged. The in-memory fast path builds no driver and prints no such line. Production spill blocks hold whole records, so on a real sort every slot record is borrowed. The driver counts only assembled records; borrowed is records_merged minus assembled. The --sort-stats help and the performance guide document the line.  The new merge_driver bench (MergeDriver::from_slots over k slots fed through the slot API, 1 M interleaved 450-byte records, each produced record appended to a reused buffer; three interleaved copy/borrow pairs on a shared laptop): k = 2, 47.6 -> 31.7 ns/record (median; borrow faster in every pair); k = 27, 126.0 -> 102.7 ns/record.  Tests, each in both spill-key layouts (a separate key prefix, and a key embedded in the record as every production key is): produced_bytes_are_borrowed_from_the_block (pointer inside the slot's block, scratch never allocated); incremental_feed_merges_identically and its embedded-key twin (proptest: blocks cut mid-record, fed one at a time on each stall; the output equals a stable merge computed in the test, and the borrowed/assembled split equals the one computed from the cuts); memory_sources_borrow_or_swap (shared residual emitted from the arena's own addresses, owned residual via the scratch). Also empty_body_on_the_borrowed_path_is_the_same_error; fast_path_sort_prints_no_merge_output_line; and the integration check that borrowed + assembled is the merge's total and that the fixture's short records are all borrowed.  Docs updated: InMemoryChunk, MemorySources, MemoryShared, slot_try_next_record, SlotMergeSource::try_next_record, PendingRecord,
- [**breaking**] Phase 1 on a split arena front — pool key extraction, planned spill blocks, borrowed queryname refs ([#1060](https://github.com/fulcrumgenomics/fgumi/pull/1060))
  - `SortBuffer::from_sorter` and `SortBufferConfig::from_sorter` are removed; build a `SortBufferConfig` directly and call `SortBuffer::new`.  * perf(sort)!: arena-borrowed queryname refs — QnKey window + class table, 40 B per record, no allocation  A queryname run sorted 80-byte owned keys: every record copied its name into an inline `SmallVec` (a heap allocation past 44 bytes), and every comparison walked the whole name — for Illumina names the 22-byte
- [**breaking**] Consolidation on the chain, chain-native merge and simulate; delete the v0.7.0 engine ([#1065](https://github.com/fulcrumgenomics/fgumi/pull/1065))
  - fgumi-sort no longer provides RunMergeProgress, RunMergeSpec, RunMergerDyn or new_run_merger: run consolidation is a chain-native tail of the sort (ConsolidateRuns, WriteConsolidatedRun). SpillWrite takes its limit and its consolidation registry in one with_consolidation(limit, registry) call (with_max_temp_files is gone). fgumi-pipeline-core's Step::phase_cap, cap_label and bound_caps are replaced by the single Step::step_cap hook returning a StepCapRef. SpillKeyKind now lives in spill_block.rs and still resolves as fgumi_sort::SpillKeyKind.  * feat(merge)!: chain-native k-way merge over the spill planner  `fgumi merge` now runs on the sort's phase 2. AnnounceMergeInputs opens each input BAM as a body-keyed merge slot: its BGZF stream read from byte 0 through the spill read planner, its decompressed header skipped, its records keyed from their bodies. A regular BGZF file is read positionally; a pipe, FIFO, process substitution, SAM or plain-gzip input is read once, in order, through a SequentialSource (one read in flight per slot), normalized to BGZF as every other input reader does. SortMerge merges the slots with an input-order check (the unsorted-input error names the input and record) and a per-input record rewrite that renames conflicting RG/PG ids before the key is extracted, as the command did before. Each body-keyed record is held to the BAM layout (the 32-byte core and the variable-length fields it declares) before a key is extracted, so a malformed input record is an error naming the input and the record, not a panic. A regular file is classified (stat, then a full format prefix) by the same open_seekable_bgzf that the sort's native input path now uses, and is opened once. The merge builds its runtime through

### Refactor

- [**breaking**] Remove the legacy unified_pipeline engine ([#947](https://github.com/fulcrumgenomics/fgumi/pull/947))
- Stop configuring the chain sort through RawExternalSorter ([#999](https://github.com/fulcrumgenomics/fgumi/pull/999))

<!-- generated by git-cliff -->

## [0.7.0] - 2026-08-24

### Features

- [**breaking**] Put the performance diagnostics behind --sort-stats ([#826](https://github.com/fulcrumgenomics/fgumi/pull/826))
- Split the input reader's time between the disk and framing ([#838](https://github.com/fulcrumgenomics/fgumi/pull/838))

### Performance

- Read ahead deeply on the file the merge is blocked on ([#806](https://github.com/fulcrumgenomics/fgumi/pull/806))
- Measure where the merge consumer's park time actually goes ([#813](https://github.com/fulcrumgenomics/fgumi/pull/813))
- Give read-ahead depth to the file the merge is actually on ([#814](https://github.com/fulcrumgenomics/fgumi/pull/814))
- Name which floor the merge is on, and stop shouting the detail ([#823](https://github.com/fulcrumgenomics/fgumi/pull/823))
- [**breaking**] Cut the Phase 1 ingest thread's per-record cost ([#829](https://github.com/fulcrumgenomics/fgumi/pull/829))
- Overlap key extraction, reuse arena pages, and read with real queue depth ([#846](https://github.com/fulcrumgenomics/fgumi/pull/846))

### Refactor

- Make --sort-stats run-scoped instead of a process-global static ([#855](https://github.com/fulcrumgenomics/fgumi/pull/855))

### Testing

- Bound the Phase 1 input path's cost with two harnesses ([#837](https://github.com/fulcrumgenomics/fgumi/pull/837))

<!-- generated by git-cliff -->

## [0.6.0] - 2026-08-14

### Bug Fixes

- Report the thread counts each phase actually uses ([#691](https://github.com/fulcrumgenomics/fgumi/pull/691))
- Stream `-o -` to stdout instead of creating a file named `-` ([#689](https://github.com/fulcrumgenomics/fgumi/pull/689))

### Features

- [**breaking**] Size the temp-file limit to the process fd budget ([#705](https://github.com/fulcrumgenomics/fgumi/pull/705))

### Performance

- Make `compare bams` faster than the sort it validates ([#696](https://github.com/fulcrumgenomics/fgumi/pull/696))
- Break the k-way merge into its component costs ([#706](https://github.com/fulcrumgenomics/fgumi/pull/706))
- Say why the merge stalled, and trace every stage a spill block passes through ([#708](https://github.com/fulcrumgenomics/fgumi/pull/708))
- Stop the k-way merge idling when the input is already sorted ([#709](https://github.com/fulcrumgenomics/fgumi/pull/709))

<!-- generated by git-cliff -->

## [0.5.0] - 2026-07-28

### Bug Fixes

- Write SS sub-sort tag as <sort-order>:<sub-sort> for fgbio/samtools parity ([#514](https://github.com/fulcrumgenomics/fgumi/pull/514))
- Validate input sort order to prevent silent corruption (MERGE3-01) ([#519](https://github.com/fulcrumgenomics/fgumi/pull/519))
- Place secondary/supplementary reads at exact template coordinate via tc ([#529](https://github.com/fulcrumgenomics/fgumi/pull/529))
- Make the queryname path stable for exact-key ties (C2) ([#540](https://github.com/fulcrumgenomics/fgumi/pull/540))
- Emit spec sub-sort spelling `lexicographical` in @HD SS (SORT3-10) ([#567](https://github.com/fulcrumgenomics/fgumi/pull/567))
- Reserve Phase-1 input serials under the input lock ([#600](https://github.com/fulcrumgenomics/fgumi/pull/600))
- Write the BAI sidecar to the samtools path, not an extension-replaced one ([#607](https://github.com/fulcrumgenomics/fgumi/pull/607))
- Surface worker panics mid-run and fail closed on a dead permit pool ([#609](https://github.com/fulcrumgenomics/fgumi/pull/609))
- Honour --compression-level when the sort fits in memory ([#633](https://github.com/fulcrumgenomics/fgumi/pull/633))
- Honour --compression-level for the tail of a spilling sort ([#653](https://github.com/fulcrumgenomics/fgumi/pull/653))

### Documentation

- Correct --write-index threading note ([#464](https://github.com/fulcrumgenomics/fgumi/pull/464))
- Fix broken intra-doc links and gate rustdoc in CI ([#574](https://github.com/fulcrumgenomics/fgumi/pull/574))

### Features

- Add per-phase --sort-threads / --merge-threads ([#608](https://github.com/fulcrumgenomics/fgumi/pull/608))
- Add --max-temp-files to tune the spill-file consolidation limit ([#643](https://github.com/fulcrumgenomics/fgumi/pull/643))
- Accept uncompressed SAM and stdin anywhere BAM is accepted ([#644](https://github.com/fulcrumgenomics/fgumi/pull/644))

### Miscellaneous Tasks

- [**breaking**] Honest MSRV in lockstep with the toolchain + let-chain adoption ([#575](https://github.com/fulcrumgenomics/fgumi/pull/575))

### Performance

- Stream records into the sort instead of an intermediate BAM ([#593](https://github.com/fulcrumgenomics/fgumi/pull/593))
- Narrow radix passes, reuse the stored queryname NUL, and borrow record bytes on ingest ([#606](https://github.com/fulcrumgenomics/fgumi/pull/606))
- Derive the radix bound inside the first counting pass ([#622](https://github.com/fulcrumgenomics/fgumi/pull/622))
- Restore the unstable queryname sort via a totally ordered key ([#621](https://github.com/fulcrumgenomics/fgumi/pull/621))
- Pool-integrate and optimize the coordinate --write-index merge ([#647](https://github.com/fulcrumgenomics/fgumi/pull/647))

### Refactor

- Pick the compressor from the job, not the pipeline phase ([#658](https://github.com/fulcrumgenomics/fgumi/pull/658))

### Testing

- Wait for step counters to settle in active_worker_limit test ([#646](https://github.com/fulcrumgenomics/fgumi/pull/646))

<!-- generated by git-cliff -->

## [0.4.0] - 2026-06-20

<!-- generated by git-cliff -->

## [0.3.1] - 2026-06-17

<!-- generated by git-cliff -->

## [0.3.0] - 2026-06-10

### Bug Fixes

- Address remaining review feedback from PR #323 ([#325](https://github.com/fulcrumgenomics/fgumi/pull/325))
- Carry library & MI in the completely-unmapped template key ([#386](https://github.com/fulcrumgenomics/fgumi/pull/386))

### Features

- Zstd as the default temp-spill codec ([#341](https://github.com/fulcrumgenomics/fgumi/pull/341))
- Variable-width template-coordinate sort key, zero-copy merge, and slimmer arena ([#375](https://github.com/fulcrumgenomics/fgumi/pull/375))

### Miscellaneous Tasks

- Remove orphaned public-api.txt snapshots ([#376](https://github.com/fulcrumgenomics/fgumi/pull/376))

### Performance

- Follow-ups from the agilent sort profile ([#344](https://github.com/fulcrumgenomics/fgumi/pull/344))
- Share SegmentedBuf across in-memory sort chunks via Arc ([#350](https://github.com/fulcrumgenomics/fgumi/pull/350))

### Refactor

- Extract fgumi-bam-io and fgumi-sort workspace crates ([#323](https://github.com/fulcrumgenomics/fgumi/pull/323))
- Demote algorithm-internal info! logs to debug ([#367](https://github.com/fulcrumgenomics/fgumi/pull/367))

<!-- generated by git-cliff -->

### Refactor

- Extracted the sort engine into the new `fgumi-sort` crate. The main `fgumi` binary now consumes it as a workspace dependency; behavior is unchanged.
