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

// ============================================================================
// Built-in impls for foreign/primitive types
// ============================================================================

/// `Vec<T>` heap estimate is the element capacity times element size.
impl<T: MemoryEstimate> MemoryEstimate for Vec<T> {
    fn estimate_heap_size(&self) -> usize {
        self.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.capacity() * std::mem::size_of::<T>()
    }
}

/// `u8` has no heap allocation.
impl MemoryEstimate for u8 {
    fn estimate_heap_size(&self) -> usize {
        0
    }
}

/// Unit type has no heap allocation (used as a zero-weight sentinel in tests).
impl MemoryEstimate for () {
    fn estimate_heap_size(&self) -> usize {
        0
    }
}

/// [`noodles::sam::alignment::RecordBuf`] heap estimate covering name,
/// sequence, quality scores, CIGAR, and auxiliary data fields.
impl MemoryEstimate for noodles::sam::alignment::RecordBuf {
    fn estimate_heap_size(&self) -> usize {
        let name_size = self.name().map_or(0, |n| n.len());
        let seq_len = self.sequence().len();
        let qual_len = self.quality_scores().len();
        let cigar_ops = self.cigar().as_ref().len();
        let cigar_size = cigar_ops * 4;
        let data_fields = self.data().iter().count();
        let entry_capacity = (data_fields * 115) / 100 + 1;
        let entries_size = data_fields * 48;
        let hash_table_size = entry_capacity * 16;
        let value_heap_size = data_fields * 16;
        let data_size = entries_size + hash_table_size + value_heap_size;
        name_size + seq_len + qual_len + cigar_size + data_size
    }
}

/// [`fgumi_raw_bam::RawRecord`] heap estimate: the raw bytes buffer capacity.
impl MemoryEstimate for fgumi_raw_bam::RawRecord {
    fn estimate_heap_size(&self) -> usize {
        self.capacity()
    }
}
