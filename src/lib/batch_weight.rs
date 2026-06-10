//! The [`BatchWeight`] trait: how the pipeline sizes batches by template count.
//!
//! Lives in the main crate (not `fgumi-bam-io`) because it is implemented for
//! a mix of foreign types (`RawRecord`, `RecordBuf`, `Vec<u8>`) and main-crate
//! types (template/group batches). The orphan rule forbids implementing a
//! foreign trait for foreign types, so the trait must be local to the crate
//! that owns those impls — and conceptually it is a pipeline-batching concern,
//! not a BAM-I/O one.

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
