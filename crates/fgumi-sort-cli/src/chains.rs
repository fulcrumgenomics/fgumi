//! Sort-specific finalize hooks and the `SortBamFile` step factory.
//!
//! These were lifted out of the umbrella's `pipeline::chains::commands::sort`
//! module so that `fgumi sort` can build its pipeline directly via
//! [`fgumi_pipeline_core::Pipeline::builder`] without depending on the
//! umbrella's monolithic `ChainBuilder`.

use std::path::PathBuf;
use std::sync::Arc;

use anyhow::Result;
use bytesize::ByteSize;
use fgumi_cli_common::{MemoryLimit, OperationTimer, format_duration};
use fgumi_pipeline_core::FinalizeHook;
use fgumi_pipeline_io::SortBamFile;
use fgumi_sort::{RawExternalSorter, SortOrder};
use log::info;
use parking_lot::Mutex;

use crate::sort::{SortOptions, TMP_DIRS_ENV, parse_cell_tag, resolve_tmp_dirs};
use crate::version;

// ─────────────────────────────────────────────────────────────────────────────
// SortFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize action for sort. Reads `SortStats` out of the
/// shared slot, logs the Summary block, and calls
/// [`OperationTimer::log_completion`].
///
/// This is the single definition of the sort summary hook, implementing the
/// shared [`FinalizeHook`] contract. The umbrella `fgumi` crate re-exports it
/// (rather than redefining a parallel copy) so its `ChainBuilder::add_sort` can
/// register it alongside the other stages' hooks (X1-005).
pub struct SortFinalizeHook {
    /// Shared slot the `SortBamFile` step fills after `sort()` returns.
    pub stats_slot: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
    /// Path of the finished BAM (logged in the Summary block).
    pub output_path: PathBuf,
    /// Timer started when the sort began (logs wall-time on completion).
    pub timer: OperationTimer,
}

impl FinalizeHook for SortFinalizeHook {
    /// Run the post-pipeline summary logging. Infallible — always returns
    /// `Ok(())` — but typed as `Result` to satisfy the [`FinalizeHook`]
    /// contract so it can share a `Vec<Box<dyn FinalizeHook>>` with the
    /// fallible hooks.
    fn finalize(self: Box<Self>) -> Result<()> {
        let SortFinalizeHook { stats_slot, output_path, timer } = *self;

        let stats = stats_slot.lock().take().unwrap_or_default();
        info!("=== Summary ===");
        info!("Records processed: {}", stats.total_records);
        info!("Records written: {}", stats.output_records);
        if stats.chunks_written > 0 {
            info!("Temporary chunks: {}", stats.chunks_written);
        }
        info!("Output: {}", output_path.display());

        timer.log_completion(stats.total_records);

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// IndexBamFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline action that builds a BAI index for a finished
/// coordinate-sorted BAM, writing the `.bai` sidecar next to the BAM.
///
/// Reads the BAM via `noodles::bam::fs::index` (which walks BGZF block offsets
/// and computes virtual positions from the on-disk layout) and writes the
/// sidecar via `fgumi_bam_io::write_bai_sidecar`, which derives the path by
/// appending `.bai` to the full output path (samtools convention) — so
/// `foo.bam` → `foo.bam.bai` and a non-`.bam` output like `foo.sorted` →
/// `foo.sorted.bai`. This
/// decouples indexing from the sort step itself: the sorter no longer needs
/// single-threaded BGZF compression to track virtual offsets during write, and
/// BAI generation lifts cleanly to a hook that any future BAM-with-index chain
/// can register.
///
/// This is the single definition of the BAI indexer hook, implementing the
/// shared [`FinalizeHook`] contract; the umbrella `fgumi` crate re-exports it
/// rather than redefining a parallel copy (X1-005).
///
/// **Invariant:** the BAM at `output_path` must be writer-closed before
/// `finalize` runs. This is guaranteed by `Pipeline::run` returning (which only
/// happens after every step's writer has dropped); the hook does not call
/// `fsync` because the same process re-reading via the OS page cache will see
/// all flushed bytes on every supported local filesystem. A remote/NFS
/// deployment would need its own close-to-open coherency story; that is out of
/// scope for this hook.
pub struct IndexBamFinalizeHook {
    /// Path of the finished coordinate-sorted BAM to index.
    pub output_path: PathBuf,
}

impl FinalizeHook for IndexBamFinalizeHook {
    /// Build and write the BAI index alongside the BAM.
    ///
    /// # Errors
    ///
    /// Returns an error if indexing or writing the BAI sidecar fails.
    fn finalize(self: Box<Self>) -> Result<()> {
        use std::time::Instant;

        let IndexBamFinalizeHook { output_path } = *self;

        info!("Indexing BAM: {}", output_path.display());
        let start = Instant::now();

        // Index + write the `<output>.bai` sidecar in one call; the path is
        // derived by appending `.bai` to the full BAM path (samtools
        // convention), correct for any path — not just `*.bam`.
        let index_path = fgumi_bam_io::write_bai_sidecar(&output_path)?;
        info!("Wrote BAM index: {} ({})", index_path.display(), format_duration(start.elapsed()));

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factory
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for [`build_sort_step`]. Bundles all captures that the
/// factory needs so the call site stays readable.
pub struct SortStepCaptures {
    /// The resolved per-stage sort options (memory, tmp dirs, order, …).
    pub sort: SortOptions,
    /// Input BAM/SAM path.
    pub input_path: PathBuf,
    /// Output BAM path.
    pub output_path: PathBuf,
    /// Number of threads for the sort engine's worker pool (the overall /
    /// Phase-2 "writing" count).
    pub num_sorter_threads: usize,
    /// Worker count for Phase 1 (accumulate/sort/spill); defaults to
    /// `num_sorter_threads`. Lower to cede CPU to an upstream producer.
    pub num_phase1_threads: usize,
    /// Worker count for Phase 2 (merge/write); defaults to `num_sorter_threads`.
    pub num_phase2_threads: usize,
    /// Resolved memory budget, already computed by the caller via
    /// [`resolve_memory_budget`](fgumi_cli_common::resolve_memory_budget).
    /// Passed in (rather than recomputed here) so
    /// the configured value the caller logged matches the value the sorter
    /// uses, and to avoid a redundant `detect_total_memory()` syscall.
    pub effective_memory: usize,
    /// Output BGZF compression level.
    pub output_compression: u32,
    /// Full command-line string for the `@PG` record.
    pub command_line: String,
    /// Shared slot for post-run stats retrieval. The caller (finalize hook)
    /// holds the other end; the step fills it when `RawExternalSorter::sort`
    /// returns.
    pub stats_slot: Arc<Mutex<Option<fgumi_sort::SortStats>>>,
    /// Wrap the sorter's input in a userspace async prefetch reader (the
    /// standalone `--async-reader` flag). Off for the fused runall sort.
    pub async_reader: bool,
}

/// Build the `SortBamFile` step, constructing and fully configuring the
/// [`RawExternalSorter`].
///
/// # Errors
///
/// Returns an error if the `MemoryLimit::Auto` initial-capacity calculation
/// overflows, or the cell tag is invalid for the sort order.
pub fn build_sort_step(cap: SortStepCaptures) -> Result<SortBamFile> {
    let sort_order: SortOrder = cap.sort.order.into();
    let cell_tag = parse_cell_tag(cap.sort.order)?;

    let effective_memory = cap.effective_memory;

    let mut sorter = RawExternalSorter::new(sort_order)
        .memory_limit(effective_memory)
        .threads(cap.num_sorter_threads)
        .sort_threads(cap.num_phase1_threads)
        .merge_threads(cap.num_phase2_threads)
        .output_compression(cap.output_compression)
        .temp_compression(cap.sort.temp_compression)
        .spill_codec(cap.sort.temp_codec)
        .key_types(cap.sort.key_types.unwrap_or_default())
        .async_reader(cap.async_reader)
        .pg_info(version::version_string(), cap.command_line);

    if matches!(cap.sort.max_memory, MemoryLimit::Auto) {
        let init = 768_usize
            .checked_mul(1024 * 1024)
            .and_then(|b| b.checked_mul(cap.num_sorter_threads))
            .ok_or_else(|| anyhow::anyhow!("initial auto buffer size overflowed"))?;
        sorter = sorter.initial_capacity(effective_memory.min(init));
    }

    if let Some(ct) = cell_tag {
        sorter = sorter.cell_tag(ct);
    }

    let env_value = std::env::var(TMP_DIRS_ENV).ok();
    let resolved_tmp_dirs = resolve_tmp_dirs(&cap.sort.tmp_dirs, env_value.as_deref());
    if !resolved_tmp_dirs.is_empty() {
        sorter = sorter.temp_dirs(resolved_tmp_dirs);
    }

    Ok(SortBamFile::new(sorter, cap.input_path, cap.output_path, cap.stats_slot))
}

// ─────────────────────────────────────────────────────────────────────────────
// Logging helper
// ─────────────────────────────────────────────────────────────────────────────

/// Log the sort startup banner.
pub fn log_sort_start(
    sort: &SortOptions,
    input_path: &std::path::Path,
    output_path: &std::path::Path,
    num_sorter_threads: usize,
    num_phase1_threads: usize,
    num_phase2_threads: usize,
    effective_memory: usize,
) {
    let cell_tag = parse_cell_tag(sort.order).unwrap_or(None);
    info!("Starting Sort");
    info!("Input: {}", input_path.display());
    info!("Output: {}", output_path.display());
    info!("Sort order: {:?}", sort.order);
    if let Some(ct) = cell_tag {
        let ct_bytes = *ct;
        info!("Cell tag: {}{}", ct_bytes[0] as char, ct_bytes[1] as char);
    }
    if let MemoryLimit::Fixed(per_thread) = sort.max_memory {
        if sort.memory_per_thread {
            info!(
                "Max memory: {} ({}/thread x {} threads)",
                ByteSize(effective_memory as u64),
                ByteSize(per_thread as u64),
                num_sorter_threads
            );
        } else {
            info!("Max memory: {} (fixed)", ByteSize(effective_memory as u64));
        }
    }
    if num_phase1_threads == num_sorter_threads && num_phase2_threads == num_sorter_threads {
        info!("Threads: {num_sorter_threads}");
    } else {
        info!(
            "Threads: {num_sorter_threads} (sort phase {num_phase1_threads}, merge phase {num_phase2_threads})"
        );
    }
    info!("Temp compression level: {}", sort.temp_compression);
    let env_value = std::env::var(TMP_DIRS_ENV).ok();
    let resolved_tmp_dirs = resolve_tmp_dirs(&sort.tmp_dirs, env_value.as_deref());
    if !resolved_tmp_dirs.is_empty() {
        let joined = resolved_tmp_dirs
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(", ");
        info!("Temp directories: {joined}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_bam_io::create_raw_bam_writer;
    use noodles::sam::Header;
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
    use std::num::NonZeroUsize;
    use tempfile::tempdir;

    /// Build a minimal two-reference, `SO:coordinate` header for indexer
    /// round-trip tests.
    fn test_header() -> Header {
        let mut builder = Header::builder();
        let chr1 = Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero"));
        let chr2 = Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero"));
        builder = builder.add_reference_sequence(b"chr1", chr1);
        builder = builder.add_reference_sequence(b"chr2", chr2);
        let header = builder.build();
        fgumi_sort::create_output_header(fgumi_sort::SortOrder::Coordinate, &header)
    }

    /// Build the raw BAM bytes for a single 10-base mapped record.
    #[allow(clippy::cast_possible_truncation)]
    fn record_bytes(ref_id: i32, pos: i32, name: &[u8]) -> Vec<u8> {
        let name_with_null = name.len() + 1;
        let padding = (4 - (name_with_null % 4)) % 4;
        let l_read_name = (name_with_null + padding) as u8;

        let mut record = Vec::with_capacity(64);
        record.extend_from_slice(&ref_id.to_le_bytes());
        record.extend_from_slice(&pos.to_le_bytes());
        record.push(l_read_name);
        record.push(60_u8); // mapq
        record.extend_from_slice(&4681_u16.to_le_bytes()); // bin
        record.extend_from_slice(&1_u16.to_le_bytes()); // n_cigar_op
        record.extend_from_slice(&0_u16.to_le_bytes()); // flag
        record.extend_from_slice(&10_u32.to_le_bytes()); // l_seq
        record.extend_from_slice(&(-1_i32).to_le_bytes()); // next_ref_id
        record.extend_from_slice(&(-1_i32).to_le_bytes()); // next_pos
        record.extend_from_slice(&0_i32.to_le_bytes()); // tlen
        record.extend_from_slice(name);
        record.push(0);
        record.extend(std::iter::repeat_n(0_u8, padding));
        let cigar_op: u32 = 10 << 4; // 10M
        record.extend_from_slice(&cigar_op.to_le_bytes());
        record.extend_from_slice(&[0x11_u8; 5]); // packed seq (AAAAAAAAAA)
        record.extend_from_slice(&[30_u8; 10]); // qualities
        record
    }

    /// Smoke test: the hook reads a finished coordinate-sorted BAM and
    /// emits a parseable BAI alongside it.
    #[test]
    fn index_bam_finalize_hook_produces_readable_bai() {
        let dir = tempdir().expect("tempdir");
        let bam_path = dir.path().join("test.bam");

        let header = test_header();
        let mut writer =
            create_raw_bam_writer(&bam_path, &header, 1, 6).expect("create_raw_bam_writer");
        for (ref_id, pos, name) in
            [(0_i32, 100_i32, b"r1".as_ref()), (0, 200, b"r2"), (1, 50, b"r3")]
        {
            writer.write_raw_record(&record_bytes(ref_id, pos, name)).expect("write");
        }
        writer.finish().expect("finish");

        let hook = Box::new(IndexBamFinalizeHook { output_path: bam_path.clone() });
        hook.finalize().expect("hook finalize");

        let bai_path = fgumi_bam_io::bai_sidecar_path(&bam_path);
        assert!(bai_path.exists(), "BAI not created at {}", bai_path.display());
        let index = noodles::bam::bai::fs::read(&bai_path).expect("read bai");
        assert!(!index.reference_sequences().is_empty(), "BAI has no reference sequences");
    }
}
