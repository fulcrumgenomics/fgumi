//! The `IndexBamFinalizeHook` post-pipeline BAI indexer for `fgumi sort`.
//!
//! This was lifted out of the umbrella's `pipeline::chains::commands::sort`
//! module so the umbrella's `ChainBuilder::add_sort` can register it from a
//! framework-light crate. The umbrella re-exports it rather than redefining a
//! parallel copy (X1-005).

use std::path::PathBuf;

use anyhow::Result;
use fgumi_cli_common::format_duration;
use fgumi_pipeline_core::FinalizeHook;
use log::info;

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
