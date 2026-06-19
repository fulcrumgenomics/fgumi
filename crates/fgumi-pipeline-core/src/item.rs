//! Item-level orthogonal traits used by the queue layer.
//!
//! These traits are independent of each other:
//!   - [`HeapSize`] is required only for items in a `ByteBoundedQueue<T>`.
//!   - [`Ordered`] is required only for items routed through a `ReorderStage<T>`.
//!
//! An item type may implement both, neither, or just one, depending on which
//! queues and operators it flows through. The framework checks bounds at
//! queue construction time (see `core/queues.rs` and `core/reorder.rs`).

/// Approximate heap footprint, in bytes, of an item.
///
/// Used by `ByteBoundedQueue<T>` to enforce a memory budget rather than an
/// item-count budget. Items that hold variable-size buffers (BAM batches,
/// FASTQ batches, decompressed BGZF blocks) implement this manually.
///
/// **Default impl.** The default method body returns `0`, so an empty
/// `impl HeapSize for MyType {}` opts in but reports zero bytes — every
/// item then fits in any byte-bounded queue regardless of the limit. Step
/// authors must override `heap_size` for any type whose actual footprint
/// matters.
///
/// **No blanket impl.** We do *not* provide `impl<T> HeapSize for T` so that
/// pushing a non-impl type into a `ByteBoundedQueue` is a compile error
/// rather than a silent zero-count.
pub trait HeapSize {
    fn heap_size(&self) -> usize {
        0
    }
}

/// Source-step input type. Carries no heap allocation.
impl HeapSize for () {}

// Primitive impls for the integer/string types that flow through tests
// and pipeline plumbing. Production BAM/FASTQ types impl `HeapSize`
// explicitly with their actual heap footprint; primitives report zero
// (their `Vec` capacity, when applicable, is reported by the wrapping
// type's own impl).
impl HeapSize for u8 {}
impl HeapSize for u16 {}
impl HeapSize for u32 {}
impl HeapSize for u64 {}
impl HeapSize for usize {}
impl HeapSize for i8 {}
impl HeapSize for i16 {}
impl HeapSize for i32 {}
impl HeapSize for i64 {}
impl HeapSize for isize {}
impl HeapSize for bool {}
impl HeapSize for f32 {}
impl HeapSize for f64 {}
impl HeapSize for char {}

impl HeapSize for String {
    fn heap_size(&self) -> usize {
        self.capacity()
    }
}

impl<T: HeapSize> HeapSize for Vec<T> {
    fn heap_size(&self) -> usize {
        self.capacity() * std::mem::size_of::<T>()
            + self.iter().map(HeapSize::heap_size).sum::<usize>()
    }
}

impl<T: HeapSize> HeapSize for Option<T> {
    fn heap_size(&self) -> usize {
        self.as_ref().map_or(0, HeapSize::heap_size)
    }
}

/// Producer-emitted serial ordinal.
///
/// Items routed through a `ReorderStage<T>` carry a monotonically-increasing
/// ordinal assigned by the producer step at push time (the framework
/// allocates from a per-branch `AtomicU64`). The reorder stage uses this
/// ordinal to deliver items to the consumer in producer-emitted order
/// regardless of inter-thread arrival skew.
///
/// Steps don't implement this directly. The framework wraps the user's
/// pushed item in an internal `Sequenced<T>` newtype that carries the
/// ordinal alongside the item; `Sequenced<T>` impls `Ordered`. See the
/// reorder stage for the wrapper type.
pub trait Ordered {
    fn ordinal(&self) -> u64;
}

#[cfg(test)]
mod tests {
    use super::*;

    struct WithHeap(Vec<u8>);
    impl HeapSize for WithHeap {
        fn heap_size(&self) -> usize {
            self.0.len()
        }
    }

    struct WithoutHeap;
    impl HeapSize for WithoutHeap {} // explicit empty impl: opts in, default returns 0

    struct WithOrdinal {
        o: u64,
    }
    impl Ordered for WithOrdinal {
        fn ordinal(&self) -> u64 {
            self.o
        }
    }

    #[test]
    fn heap_size_default_is_zero() {
        assert_eq!(WithoutHeap.heap_size(), 0);
    }

    #[test]
    fn heap_size_override_works() {
        assert_eq!(WithHeap(vec![0; 1024]).heap_size(), 1024);
    }

    #[test]
    fn ordered_returns_ordinal() {
        assert_eq!(WithOrdinal { o: 42 }.ordinal(), 42);
    }
}
