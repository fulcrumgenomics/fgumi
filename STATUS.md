# Code Review Status

Tracking cleanup and simplification work across the codebase.

Last commit: dcc1a78

## How to Resume

To continue this code review, use the following prompt:

> Lets continue reviewing modules in STATUS.md. Go through each pending module one by one looking for opportunities to remove dead code, refactor to use existing crates or library functions, and simplify code while keeping functional equivalence and performance. Prompt me before proceeding to the next module. Keep STATUS.md up to date.

**Next module to review**: `metrics/clip.rs`

**Pending work**: `vanilla_caller.rs` has 8 identified performance optimizations (see consensus/ section below)

### Re-review Queue (modified since last review)

| Module | Status | Reason |
|--------|--------|--------|
| `reorder_buffer.rs` | Done | Simplified `try_pop_next()` - removed 18 lines of redundant checks |
| `bgzf_reader.rs` | Done | Clean. Fixed `base.rs` to use `decompress_block_into()` (avoids per-block allocation) |
| `bgzf_writer.rs` | Done | Clean - buffer pool well-implemented |
| `bam_io.rs` | Done | Clean - stdin/streaming support well-implemented |
| `grouper.rs` | Done | Removed ~188 lines dead code (duplicate MiGroup/MiGrouper - use mi_group.rs instead) |
| `mi_group.rs` | Done | Clean - no changes needed |
| `read_info.rs` | Done | Added docs explaining why LibraryLookup/LibraryIndex and two unclipped position functions coexist |

**Workflow pattern**:
1. Read the module and identify dead code, unused imports, simplification opportunities
2. Present findings and get user approval before making changes
3. Make changes, run `cargo build`, `cargo test`, `cargo ci-lint`
4. Update STATUS.md with what was done
5. Prompt user before moving to next module

## src/lib Modules

### Reviewed

| Module                  | Status | Changes                                                                                                                                                  |
|-------------------------|--------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| `bam_io.rs`             | Done   | Consolidated 5 functions → 3. Renamed `open_bam_reader_auto` → `create_bam_reader`. Removed unused `open_bam_reader`, `open_bam_reader_mt`.             |
| `batched_sam_reader.rs` | Done   | Minor cleanup                                                                                                                                            |
| `bgzf_reader.rs`        | Done   | 486 → 371 lines (~24% reduction). Removed dead code: `BGZF_MAX_UNCOMPRESSED_SIZE`, `decompress_block_verified`, `decompress_blocks`, `RawBgzfBatchReader`. Made `read_raw_block` private. |
| `bgzf_writer.rs`        | Done   | 237 → 210 lines (~11% reduction). Removed: `compression_level` field/getter, `pending_blocks()`, `impl Write for InlineBgzfCompressor`. Made `BGZF_MAX_BLOCK_SIZE` private. |
| `bitenc.rs`             | Done   | Major cleanup (~229 lines removed)                                                                                                                       |
| `dna.rs`                | Done   | Major cleanup (~166 lines removed)                                                                                                                       |
| `errors.rs`             | Done   | Major cleanup (~467 lines removed)                                                                                                                       |
| `grouper.rs`            | Done   | Cleanup (~45 lines removed)                                                                                                                              |
| `logging.rs`            | Done   | Major cleanup (~174 lines removed)                                                                                                                       |
| `mod.rs`                | Done   | Cleanup (~40 lines removed)                                                                                                                              |
| `parallel.rs`           | Done   | File removed (dead code)                                                                                                                                 |
| `progress.rs`           | Done   | Cleanup (~46 lines removed)                                                                                                                              |
| `reference.rs`          | Done   | Cleanup (~82 lines removed)                                                                                                                              |
| `rejection.rs`          | Done   | Major cleanup (~326 lines removed)                                                                                                                       |
| `reorder_buffer.rs`     | Done   | 334 → 227 lines (~32% reduction). Simplified API: `new(capacity)` → `new()`. Removed 5 unused methods. Added `Default` impl.                            |
| `tag_reversal.rs`       | Done   | Major cleanup (~192 lines removed)                                                                                                                       |
| `template_batch.rs`     | Done   | File removed (dead code, ~370 lines)                                                                                                                     |
| `template_reader.rs`    | Done   | File removed (~121 lines). Replaced with `TemplateIterator` from `template.rs`. Moved `TemplateBatch` type alias to `template.rs`.                      |
| `variant_review.rs`     | Done   | Minor cleanup. Removed `BaseCounts::new()` (using `Default`). Moved `extract_mi_base` to `umi/mod.rs`.                                                   |
| `fastq.rs`              | Done   | 755 → 682 lines (~10% reduction). Removed unused `sample_barcode_sequence()`. Refactored `ReadSetIterator::next()` to use `from_record_with_structure()`. |
| `phred.rs`              | Done   | 845 → 695 lines (~18% reduction). Removed unused `PhredLookupTable`. Made `LN_TWO`/`LN_FOUR_THIRDS` private.                                               |
| `mi_group.rs`           | Done   | 1560 → 953 lines (~39% reduction). Removed dead code: `MiGroupReader`, `OwnedMiGroupReader`, `KeyExtractor`, `RecordFilter`, `StripStrandSuffix`, `NoTransform`, `NoFilter`, type aliases. Replaced `strip_strand_suffix` with `extract_mi_base` from `umi/mod.rs`. |
| `position_group.rs`     | Done   | File removed (~729 lines). Entire module was dead code - nothing imported or used it.                                                                       |
| `read_info.rs`          | Done   | 766 → 743 lines. Performance optimizations: (1) `Arc<str>` for library names avoids per-read string cloning, (2) fixed redundant `r1()/r2()` calls, (3) deferred error message allocation, (4) single-pass CIGAR for reverse strand, (5) tuple comparison for position normalization. |
| `validation.rs`         | Done   | No changes needed. Clean validation utilities, all functions in use.                                                                                         |
| `clipper.rs`            | Done   | ~70 lines removed. (1) Added `array_len!`/`slice_array!` macros to reduce repetitive array handling code. (2) `SamRecordClipper::clipped_bases` now delegates to `cigar_utils::clipped_bases`. (3) Added combined `leading_clips`/`trailing_clips` functions using direct slice access (eliminates Vec allocations). (4) Direct slice indexing for sequence/quality clipping. |
| `template.rs`           | Done   | 3560 → 3361 lines (~200 lines removed). Removed unused `BufferedTemplateIterator`. Performance optimizations: (1) `Builder::push` avoids allocation when comparing names - only allocates on first record, (2) `Builder::build` pre-allocates records Vec with exact capacity, (3) `cigar_to_string` pre-allocates String capacity. |

### Pending Review

(No pending top-level modules)

### Subdirectories - Pending Review

#### consensus/
- `base_builder.rs` - Done. Removed unused `_adjusted_error_table` field. Optimized `observations_for_base` to use O(1) lookup table.
- `caller.rs` - Done. Fixed incomplete test for rejection reason codes (added 5 missing variants). No dead code.
- `codec_caller.rs` - Done. 5301→3707 lines (~30% reduction). Removed old reference-based implementation: `consensus_builder` field, `build_codec_duplex_consensus`, `compute_overlap_region`, `get_leading/trailing_soft_clip_*`, `compute_aligned_length_from_cigar`, `get_base_at_ref_pos`, `build_duplex_consensus`, `get_ss_base_at_ref_pos`, `mask_consensus_quals`, `filter_empty_positions`, `build_single_strand_consensus`. Removed 25 tests for deleted functions. Cleaned up unused imports.
- `duplex_caller.rs` - Done. Removed ~70 lines dead code (`partition_reads_by_strand`). Performance optimizations: (1) `are_all_same_strand` takes iterator instead of slice, (2) R1/R2 filtering done once and reused, (3) avoid Vec allocations for strand validation and cell barcode extraction, (4) avoid intermediate x_records/y_records Vecs.
- `filter.rs` - Done. 1917→1870 lines. Removed dead code: `group_by_template` function. Performance optimizations: (1) `mask_bases` and `mask_duplex_bases` collect directly into mutable vectors (avoids 2 clones each), (2) `mean_base_quality` uses iterator `.zip()` instead of allocating vectors, (3) added `extract_consensus_bases` helper to deduplicate ac/bc tag extraction, (4) `mask_bases` uses `extract_per_base_array` helper.
- `mod.rs` - Done. No changes needed. Clean module organization with re-exports.
- `overlapping.rs` - Done. 1232→1196 lines (~3% reduction). Major performance fix: O(n²)→O(n) by cloning seq/qual ONCE at start of `call()` instead of per-position. Other changes: (1) `apply_overlapping_consensus` returns `Result<()>` instead of cloning reads, (2) uses AHashMap, (3) simplified `process_agreement`/`process_disagreement` signatures to take offsets instead of structs, (4) removed unused `ReadMateAndRefPos.ref_pos` field.
- `sequence.rs` - Done. Fixed misleading "# Safety" doc comments → "# Note" (4 occurrences). Optimized `padded()` method to use pre-allocation with `with_capacity` + `resize`/`extend_from_slice` instead of intermediate `vec![]` + `.concat()` allocations.
- `simple_umi.rs` - Done. No changes needed. Clean module with good test coverage.
- `tags.rs` - Done. Converted `tags_to_reverse()`, `tags_to_reverse_complement()`, `all_tags()` from returning `Vec<&str>` to `&'static [&str]` using const arrays (avoids allocation per call).
- `vanilla_caller.rs` - Done. Optimized `VanillaConsensusRead::padded()` to use pre-allocation instead of `vec![]` + `.concat()` allocations. **Pending perf optimizations below.**

##### vanilla_caller.rs - Pending Performance Optimizations

**High Impact:**

1. **Lines 853-855: HashMap `.cloned()` → `.remove()`**
   ```rust
   let fragment_reads = subgroups.get(&ReadType::Fragment).cloned().unwrap_or_default();
   ```
   Each `.cloned()` clones entire `Vec<RecordBuf>`. Use `.remove()` since HashMap is consumed.

2. **Lines 863-864: Redundant clone for rejection tracking**
   ```rust
   let r1_consensus = self.process_subgroup(umi, ReadType::R1, r1_reads.clone())?;
   ```
   Clones only needed if orphan rejection. Restructure to avoid cloning in success case.

3. **Lines 1107-1113: O(positions × reads) error counting**
   ```rust
   let error_count: usize = source_reads.iter().filter(...).count();
   ```
   Iterates ALL source reads again per position. Count during main loop (lines 1089-1099).

4. **Line 27: Use AHashMap instead of std::HashMap**
   Already used elsewhere in codebase.

**Medium Impact:**

5. **Line 713: Clone entire RecordBuf into Arc**
   ```rust
   sam: Some(std::sync::Arc::new(read.clone())),
   ```
   Consider passing ownership or using `Arc<RecordBuf>` from start.

6. **Line 733: CIGAR clone in filter loop**
   ```rust
   .map(|(i, sr)| (i, sr.bases.len(), sr.simplified_cigar.clone()))
   ```
   Could use references instead.

7. **Lines 509-513: Cloning consensus data vectors**
   Could take ownership if `consensus` not used after.

**Lower Impact:**

8. **Lines 1221-1229: String allocation for UMI consensus**
   Creates `Vec<String>` just to call `consensus_umis`. Could use `&str` or `Cow<str>`.

#### metrics/
- `clip.rs`
- `consensus.rs`
- `correct.rs`
- `duplex.rs`
- `group.rs`
- `mod.rs`
- `writer.rs`

#### sam/
- `alignment_tags.rs`
- `builder.rs`
- `mod.rs`
- `record_utils.rs`

#### simulate/
- `family_size.rs`
- `fastq_writer.rs`
- `insert_size.rs`
- `mod.rs`
- `parallel_gzip_writer.rs`
- `quality.rs`
- `rng.rs`
- `strand_bias.rs`

#### umi/
- `assigner.rs`
- `mod.rs`

#### unified_pipeline/
- `bam.rs`
- `base.rs`
- `fastq.rs`
- `mod.rs`
- `scheduler.rs` (new file)

---

## src/commands

### OperationTimer Added

All commands now have consistent timing/logging via `OperationTimer`:

| Command              | Status |
|----------------------|--------|
| `extract.rs`         | Done   |
| `filter.rs`          | Done   |
| `clip.rs`            | Done   |
| `correct.rs`         | Done   |
| `downsample.rs`      | Done   |
| `zipper.rs`          | Done   |
| `duplex_metrics.rs`  | Done   |
| `compare/bams.rs`    | Done   |
| `compare/metrics.rs` | Done   |
| `review.rs`          | Done   |

Commands that already had OperationTimer:
- `group.rs`
- `simplex.rs`
- `duplex.rs`
- `codec.rs`
- `consensus_runner.rs`

---

## Summary

- **Reviewed**: 25 modules (including 4 files removed)
- **Pending**: 2 top-level modules + 6 subdirectories
- **Commands with OperationTimer**: All (15 total)
- **Net reduction**: ~4,956 lines deleted across 46 files
