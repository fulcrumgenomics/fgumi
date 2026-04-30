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

/// Heap estimate for a noodles auxiliary-data [`Value`].
///
/// Primitive variants (`Character`, `Int*`, `UInt*`, `Float`) live entirely on
/// the stack and contribute zero. `String`/`Hex` use a `BString` (a `Vec<u8>`
/// under the hood); their heap is the underlying capacity. `Array` variants
/// each wrap a `Vec<T>` whose heap is `capacity() * sizeof::<T>()`.
///
/// [`Value`]: noodles::sam::alignment::record_buf::data::field::Value
impl MemoryEstimate for noodles::sam::alignment::record_buf::data::field::Value {
    fn estimate_heap_size(&self) -> usize {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        match self {
            Value::Character(_)
            | Value::Int8(_)
            | Value::UInt8(_)
            | Value::Int16(_)
            | Value::UInt16(_)
            | Value::Int32(_)
            | Value::UInt32(_)
            | Value::Float(_) => 0,
            Value::String(s) | Value::Hex(s) => s.capacity(),
            Value::Array(array) => match array {
                Array::Int8(v) => v.capacity() * std::mem::size_of::<i8>(),
                Array::UInt8(v) => v.capacity() * std::mem::size_of::<u8>(),
                Array::Int16(v) => v.capacity() * std::mem::size_of::<i16>(),
                Array::UInt16(v) => v.capacity() * std::mem::size_of::<u16>(),
                Array::Int32(v) => v.capacity() * std::mem::size_of::<i32>(),
                Array::UInt32(v) => v.capacity() * std::mem::size_of::<u32>(),
                Array::Float(v) => v.capacity() * std::mem::size_of::<f32>(),
            },
        }
    }
}

/// [`noodles::sam::alignment::RecordBuf`] heap estimate covering name,
/// sequence, quality scores, CIGAR, and auxiliary data fields.
///
/// Auxiliary data is accounted by iterating each `(Tag, Value)` pair and
/// summing the actual `Value` heap usage rather than charging a fixed
/// per-field constant. Variable-length payloads (large `String`s, `Array`
/// values) are non-trivial — under-counting weakens the backpressure signal.
impl MemoryEstimate for noodles::sam::alignment::RecordBuf {
    fn estimate_heap_size(&self) -> usize {
        let name_size = self.name().map_or(0, |n| n.len());
        let seq_len = self.sequence().len();
        let qual_len = self.quality_scores().len();
        let cigar_ops = self.cigar().as_ref().len();
        let cigar_size = cigar_ops * 4;

        // Iterate aux fields once: count + sum of payload heap sizes.
        let mut data_fields = 0usize;
        let mut value_heap_size = 0usize;
        for (_tag, value) in self.data().iter() {
            data_fields += 1;
            value_heap_size += value.estimate_heap_size();
        }
        let entry_capacity = (data_fields * 115) / 100 + 1;
        let entries_size = data_fields * 48;
        let hash_table_size = entry_capacity * 16;
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

#[cfg(test)]
mod tests {
    use super::MemoryEstimate;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::alignment::record_buf::data::field::value::Array;

    #[test]
    fn primitive_value_heap_is_zero() {
        assert_eq!(Value::Character(b'A').estimate_heap_size(), 0);
        assert_eq!(Value::Int8(0).estimate_heap_size(), 0);
        assert_eq!(Value::UInt8(0).estimate_heap_size(), 0);
        assert_eq!(Value::Int16(0).estimate_heap_size(), 0);
        assert_eq!(Value::UInt16(0).estimate_heap_size(), 0);
        assert_eq!(Value::Int32(0).estimate_heap_size(), 0);
        assert_eq!(Value::UInt32(0).estimate_heap_size(), 0);
        assert_eq!(Value::Float(0.0).estimate_heap_size(), 0);
    }

    #[test]
    fn string_value_heap_is_capacity() {
        // BString::from(&[u8]) typically allocates with capacity >= len.
        let s = Value::String(bstr::BString::from(vec![b'x'; 1024]));
        assert!(s.estimate_heap_size() >= 1024);
    }

    #[test]
    fn array_value_heap_scales_with_element_size() {
        // 1024 i32 elements = 4096 bytes minimum.
        let a = Value::Array(Array::Int32(vec![0i32; 1024]));
        assert!(a.estimate_heap_size() >= 4096);

        // 1024 u8 elements = 1024 bytes minimum.
        let b = Value::Array(Array::UInt8(vec![0u8; 1024]));
        assert!(b.estimate_heap_size() >= 1024);
    }

    #[test]
    fn record_buf_heap_includes_aux_payload() {
        use noodles::sam::alignment::RecordBuf;
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::Data;

        // Two RecordBufs identical except one carries a large CO:Z:... aux value.
        let small_data: Data =
            [(Tag::COMMENT, Value::String(bstr::BString::from("x")))].into_iter().collect();
        let large_data: Data =
            [(Tag::COMMENT, Value::String(bstr::BString::from(vec![b'x'; 4096])))]
                .into_iter()
                .collect();

        let small = RecordBuf::builder().set_data(small_data).build();
        let large = RecordBuf::builder().set_data(large_data).build();

        let small_size = small.estimate_heap_size();
        let large_size = large.estimate_heap_size();

        // The large-payload record's estimate must exceed the small one by at
        // least roughly the payload delta (4096 - 1 ≈ 4095 bytes). Without the
        // fix, both estimates collapsed to the same per-field constant.
        assert!(
            large_size >= small_size + 3000,
            "expected large_size >> small_size, got large={large_size} small={small_size}"
        );
    }
}
