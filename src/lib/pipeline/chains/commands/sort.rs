//! Chain builder for `Stage::Sort`.
//!
//! Phase 2 (T2.16) held the full ~200-LOC chain construction here.
//! Phase 3 (T3a.6) lifts that logic into
//! [`crate::pipeline::chains::builder::ChainBuilder`]; this module
//! now holds the sort-specific types and step factory that the builder
//! imports: [`SortFinalizeHook`] and [`build_sort_step`].
//!
//! [`build_sort_chain`] shrinks to a ~10-line delegate.
//!
//! ## Sort pipeline topology
//!
//! Sort is architecturally unique among the `Stage` variants: it is
//! implemented as a single `Exclusive` [`SortBamFile`] step that reads
//! from and writes to files directly. It therefore bypasses the normal
//! source-preamble (`ReadBgzfBlocks → … → DecodeRecords`) and
//! `BgzfCompress → WriteBgzfFile` sink that every other stage uses.
//! `build_sort_chain` consequently skips `add_source` and `add_sink`;
//! `ChainBuilder::add_sort` registers `SortBamFile` via
//! `PipelineBuilder::append_source` (legal since `SortBamFile` has
//! `Input = ()`) and sets `ChainBuilder::override_pipeline_threads = 1`
//! so the framework runs a single driver thread while the sorter's own
//! `SortWorkerPool` / rayon pools provide internal concurrency.
//!
//! [`SortBamFile`]: crate::pipeline::steps::sort::SortBamFile

use std::sync::Arc;

use anyhow::{Result, bail};
use log::info;

use crate::logging::OperationTimer;
use crate::pipeline::chains::builder::ChainBuilder;
use crate::pipeline::chains::{BuiltPipeline, ChainSpec, FinalizeHook, SinkSpec, SourceSpec};

// The sort step factory, its captures bundle, and the startup-banner logger now
// live in the `fgumi-sort-cli` crate. Re-export them so the umbrella's
// `ChainBuilder::add_sort` keeps importing them from this module. Their
// signatures are independent of the umbrella's `FinalizeHook` trait, so they
// can be shared verbatim. The two finalize hooks below, by contrast, implement
// the umbrella's `FinalizeHook` trait (which `fgumi-sort-cli`'s plain-method
// versions do not), so they stay umbrella-local.
pub(crate) use fgumi_sort_cli::chains::{SortStepCaptures, build_sort_step, log_sort_start};

// ─────────────────────────────────────────────────────────────────────────────
// SortFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for sort. Reads `SortStats` out of the
/// shared slot, logs the Summary block, and calls
/// `timer.log_completion`.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_sort`.
pub(crate) struct SortFinalizeHook {
    pub(crate) stats_slot: Arc<parking_lot::Mutex<Option<fgumi_sort::SortStats>>>,
    pub(crate) output_path: std::path::PathBuf,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for SortFinalizeHook {
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

/// Post-pipeline finalize hook that builds a BAI index for a finished
/// coordinate-sorted BAM. Registered by [`ChainBuilder::add_sort`] when
/// `spec.sink` is [`SinkSpec::BamWithIndex`].
///
/// Reads the BAM via `noodles::bam::fs::index` (which walks BGZF block
/// offsets and computes virtual positions from the on-disk layout) and
/// writes `<output>.bam.bai` next to the BAM via
/// `fgumi_bam_io::write_bai_index`. This decouples indexing from the
/// sort step itself: the sorter no longer needs single-threaded BGZF
/// compression to track virtual offsets during write, and BAI
/// generation lifts cleanly to a hook that any future
/// `SinkSpec::BamWithIndex`-producing chain can register.
///
/// **Invariant:** the BAM at `output_path` must be writer-closed before
/// `finalize` runs. This is guaranteed by `Pipeline::run` returning
/// (which only happens after every step's writer has dropped); the
/// hook does not call `fsync` because the same process re-reading via
/// the OS page cache will see all flushed bytes on every supported
/// local filesystem. A remote/NFS deployment would need its own
/// close-to-open coherency story; that is out of scope for this hook.
///
/// `pub(crate)` so [`ChainBuilder`] can construct and register it in
/// `add_sort`.
///
/// [`ChainBuilder::add_sort`]: crate::pipeline::chains::builder::ChainBuilder
/// [`SinkSpec::BamWithIndex`]: crate::pipeline::chains::SinkSpec::BamWithIndex
pub(crate) struct IndexBamFinalizeHook {
    pub(crate) output_path: std::path::PathBuf,
}

impl FinalizeHook for IndexBamFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        use crate::logging::format_duration;
        use fgumi_bam_io::write_bai_index;
        use noodles::bam;
        use std::time::Instant;

        let IndexBamFinalizeHook { output_path } = *self;

        info!("Indexing BAM: {}", output_path.display());
        let start = Instant::now();

        let index = bam::fs::index(&output_path)
            .map_err(|e| anyhow::anyhow!("Failed to index {}: {e}", output_path.display()))?;

        // BAI convention: sit next to the BAM with `.bam.bai`, e.g.
        // `foo.bam` → `foo.bam.bai`. Implemented as
        // `Path::with_extension("bam.bai")` to match the legacy
        // in-sort indexer at `crates/fgumi-sort/src/external.rs`
        // (pre-Phase-4) and samtools' default. The contract is
        // "append `.bai` to the BAM path"; `with_extension` happens
        // to produce the right result because the input always
        // carries a single `.bam` extension.
        let index_path = output_path.with_extension("bam.bai");
        write_bai_index(&index_path, &index)
            .map_err(|e| anyhow::anyhow!("Failed to write BAI to {}: {e}", index_path.display()))?;
        info!("Wrote BAM index: {} ({})", index_path.display(), format_duration(start.elapsed()));

        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Step factory
//
// `SortStepCaptures`, `build_sort_step`, and `log_sort_start` are re-exported
// from `fgumi-sort-cli` at the top of this module.
// ─────────────────────────────────────────────────────────────────────────────

// ─────────────────────────────────────────────────────────────────────────────
// build_sort_chain — 10-line delegate
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for a sort-only chain.
///
/// `spec.stages == [Stage::Sort]`. Other layouts are caller errors
/// caught by the dispatch match in [`crate::pipeline::chains::build`].
///
/// Delegates to [`ChainBuilder`]. The full step construction lives in
/// [`ChainBuilder`]'s `add_sort` method. Unlike other stage chains, sort
/// bypasses `add_source` / `add_sink` because `SortBamFile` is a
/// self-contained `Exclusive` step that reads/writes files directly.
///
/// # Errors
///
/// Returns validation errors (missing sort options, invalid source/sink
/// types) or any underlying pipeline construction error.
#[allow(clippy::needless_pass_by_value)]
pub fn build_sort_chain(spec: ChainSpec) -> Result<BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = ChainBuilder::new(&spec)?;
    chain.add_stage(Stage::Sort, StagePosition::Terminal)?;
    chain.build()
}

// ─────────────────────────────────────────────────────────────────────────────
// Source/sink path helpers re-used by add_sort
// ─────────────────────────────────────────────────────────────────────────────

/// Extract the source path from a [`ChainSpec`] for sort (BAM or SAM only).
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_sort`].
pub(crate) fn sort_source_path(spec: &ChainSpec) -> Result<std::path::PathBuf> {
    match &spec.source {
        SourceSpec::Bam(p) | SourceSpec::Sam(p) => Ok(p.clone()),
        other => bail!("sort requires a BAM or SAM source, got {other:?}"),
    }
}

/// Extract the sink path from a [`ChainSpec`] for sort.
///
/// `pub(crate)` — consumed only by [`ChainBuilder::add_sort`].
pub(crate) fn sort_sink_path(spec: &ChainSpec) -> std::path::PathBuf {
    match &spec.sink {
        SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
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
    /// round-trip tests. Uses `fgumi_sort::create_output_header` so the
    /// `@HD SO:coordinate` annotation matches the bytes the real sort
    /// step emits — `noodles::bam::fs::index` (the indexer the hook
    /// invokes) rejects BAMs that don't claim coordinate sort order.
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
    /// Mirrors `fgumi_bam_io::writer::create_test_bam_record` (private to
    /// that test module). 10M CIGAR + 5-byte packed sequence + 10 quality
    /// bytes — enough to exercise the indexer's chunk resolution without
    /// pulling in a fixture file.
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
    ///
    /// Doesn't assert byte-identity with the legacy in-sort indexer
    /// (that would require the same BGZF block layout, which differs
    /// between single- and multi-threaded compression). The
    /// integration tests in `tests/integration/test_sort_write_index.rs`
    /// exercise the end-to-end semantic correctness via samtools-style
    /// region queries; this unit test just pins the hook's mechanical
    /// contract.
    ///
    /// TODO(T4.4): once `--write-index` routes through this hook in
    /// `Sort::execute`, the integration tests in
    /// `tests/integration/test_sort_write_index.rs` provide
    /// content-level cross-checks (samtools region queries return the
    /// same record sets the in-sort indexer would have returned). The
    /// plan's Risk 1 byte-identity check effectively migrates from
    /// T4.2 to T4.4 because of the BGZF-layout divergence between the
    /// legacy single-threaded path and the new multi-threaded
    /// write + post-write index path.
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

        let bai_path = bam_path.with_extension("bam.bai");
        assert!(bai_path.exists(), "BAI not created at {}", bai_path.display());
        let index = noodles::bam::bai::fs::read(&bai_path).expect("read bai");
        assert!(!index.reference_sequences().is_empty(), "BAI has no reference sequences",);
    }
}
