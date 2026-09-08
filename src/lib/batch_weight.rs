//! Batch weight trait for template-based batching.
//!
//! Moved out of `unified_pipeline::base` (C5/R6a) so `grouper`/`mi_group` can
//! depend on it without going through the (now-legacy) `unified_pipeline`
//! module. Stays in the root crate rather than `fgumi-bam-io`: its surviving
//! implementations target foreign types (`noodles::sam::alignment::RecordBuf`,
//! `fgumi_raw_bam::RawRecord`, `Vec<u8>`), which are only legal `impl`s while
//! the trait itself is defined in this crate (orphan rule).

/// Trait for groups that can report their "weight" for batching purposes.
///
/// The weight is typically the number of templates in the group, allowing
/// the pipeline to batch groups based on total templates rather than group count.
/// This provides more consistent batch sizes across datasets with varying
/// templates-per-group ratios.
///
/// # Example
///
/// For a position group with 10 templates, `batch_weight()` returns 10.
/// The pipeline accumulates groups until the total weight reaches a threshold
/// (e.g., 500 templates), then flushes the batch.
pub trait BatchWeight {
    /// Returns the weight of this group for batching purposes.
    /// For position groups, this is typically the number of templates.
    fn batch_weight(&self) -> usize;
}
