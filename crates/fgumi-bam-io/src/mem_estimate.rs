//! Memory estimation trait used by pipeline backpressure.
//!
//! Implementors expose a heap-byte estimate so per-stage queues can apply
//! memory-bounded backpressure rather than count-bounded backpressure.

/// Estimates the heap memory used by a value in bytes.
///
/// Implementations should include heap allocations the value owns
/// (`Vec`/`String` capacity, nested boxed/owned data) but should **not**
/// include the size of the struct itself or shared/reference-counted data
/// (which is counted once at the source, not per reference).
///
/// # Examples
///
/// ```
/// use fgumi_bam_io::MemoryEstimate;
///
/// struct Batch { data: Vec<u8> }
///
/// impl MemoryEstimate for Batch {
///     fn estimate_heap_size(&self) -> usize {
///         self.data.capacity()
///     }
/// }
///
/// let batch = Batch { data: vec![0u8; 1024] };
/// assert!(batch.estimate_heap_size() >= 1024);
/// ```
pub trait MemoryEstimate {
    /// Returns an estimate of the heap memory used by this value, in bytes.
    ///
    /// This should include:
    /// - `Vec`/`String` heap allocations (capacity, not just len)
    /// - Nested struct heap allocations
    ///
    /// This should **not** include:
    /// - The size of the struct itself (stack size)
    /// - Shared/reference-counted data (counted once, not per-reference)
    fn estimate_heap_size(&self) -> usize;
}
