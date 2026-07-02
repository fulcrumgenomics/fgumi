//! `StepOutputs`: type-level description of a step's outputs.
//!
//! Single-output steps declare `type Outputs = Single<T>;`.
//! Multi-output steps declare `type Outputs = (A, B, C);` (positional access).
//! Sinks declare `type Outputs = ();`.
//!
//! **Maximum arity in PR 1: 4.** [`MAX_ARITY`] is the largest tuple shape
//! that has a `StepOutputs` impl. A user step that declares
//! `type Outputs = (A, B, C, D, E);` (arity 5) will fail to compile with a
//! "trait `StepOutputs` is not implemented" error. Larger arities require
//! adding the additional impls in this file plus the matching
//! `build_tupleN_queues` constructor in `handles.rs`.
//!
//! Each `StepOutputs` impl provides:
//!   - `arity()` — number of independent output channels.
//!   - `build_queues(specs, ordering, level) -> (OutputQueueSet, OutputsViewAny)` —
//!     constructs the typed queues + reorder operators (where applicable)
//!     and the type-erased view the framework stores. Implemented in
//!     `handles.rs` (one impl per arity).
//!
//! The `specs` and `ordering` slices both have length `arity()`, sourced
//! from `StepProfile::output_queues` and `StepProfile::branch_ordering`.

use std::marker::PhantomData;

use super::handles::OutputQueueSet;
use super::item::{HeapSize, Ordered};
use super::queues::QueueSpec;
use super::reorder::BranchOrdering;
use super::step::OutputsViewAny;

/// The largest tuple `Outputs` shape that has a `StepOutputs` impl in PR 1.
/// `Single<T>` and `()` are also valid output shapes.
pub const MAX_ARITY: usize = 4;

/// Marker trait for a step's `Outputs` associated type.
pub trait StepOutputs: Send + 'static {
    /// Number of independent output channels.
    fn arity() -> usize;

    /// Construct the typed queues + type-erased view for this Outputs shape.
    ///
    /// `specs[i]` and `ordering[i]` together describe branch `i`. Panics if
    /// `specs.len() != Self::arity()` or `ordering.len() != Self::arity()`
    /// (the builder ensures this invariant before calling).
    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny);

    /// Mark all output branches drained, dispatching through the typed
    /// `OutputHandles<Self>::mark_all_drained` method. Implemented for
    /// each per-arity variant in this file; called by `TypedStep<S>` from
    /// `mark_outputs_drained` in the worker loop's drain-propagation path.
    fn mark_all_drained(handles: &super::step::OutputHandles<Self>)
    where
        Self: Sized;
}

/// Wrapper for single-output steps. `type Outputs = Single<T>;`
pub struct Single<T: Send + HeapSize + 'static>(PhantomData<fn() -> T>);

impl<T: Send + HeapSize + 'static> StepOutputs for Single<T> {
    #[inline]
    fn arity() -> usize {
        1
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_single_queues::<T>(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

/// Wrapper for single-output steps where the item type is heap-aware AND
/// carries its own ordinal. The canonical Phase 3 BAM step output shape:
/// every BAM step's output type impls both `HeapSize` (for byte-bounded
/// queues) and `Ordered` (so a `batch_serial: u64` field carries record-
/// read order through every Parallel transform).
///
/// Steps that need byte-bounded outputs but don't have item-carried serials
/// use `Single<T>` with `BranchOrdering::None` or `ByOrdinal`. Steps that
/// need item-ordinal ordering but not byte-bounded queues are uncommon
/// (and currently not supported as a separate shape) — they can use this
/// shape with `QueueSpec::CountBounded` since `T: HeapSize` is required by
/// the bound but not consulted by count-bounded queues.
pub struct OrderedBytesSingle<T: Send + HeapSize + Ordered + 'static>(PhantomData<fn() -> T>);

impl<T: Send + HeapSize + Ordered + 'static> StepOutputs for OrderedBytesSingle<T> {
    #[inline]
    fn arity() -> usize {
        1
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_single_queues_ordered_bytes::<T>(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

impl<A: Send + HeapSize + 'static, B: Send + HeapSize + 'static> StepOutputs for (A, B) {
    #[inline]
    fn arity() -> usize {
        2
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_tuple2_queues::<A, B>(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

/// Ordered + byte-bounded tuple-2 outputs. Use when a step fans out to
/// two branches that both need `BranchOrdering::ByItemOrdinal` +
/// byte-bounded queues — e.g., `filter` emitting kept and rejected
/// records as parallel ordered streams that downstream
/// `BgzfCompress` / `WriteBgzfFile` sinks can consume.
///
/// The plain `(A, B)` `StepOutputs` impl only bounds `A: HeapSize`
/// (suitable for `Process2`'s `BranchOrdering::None` case). To use
/// `ByItemOrdinal` on either branch, both branches must satisfy
/// `Ordered + HeapSize` — that's what this shape encodes at the type
/// level.
pub struct OrderedBytesTuple2<A, B>(PhantomData<fn() -> (A, B)>)
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static;

impl<A, B> StepOutputs for OrderedBytesTuple2<A, B>
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
{
    #[inline]
    fn arity() -> usize {
        2
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_tuple2_queues_ordered_bytes::<A, B>(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

impl<A, B, C> StepOutputs for (A, B, C)
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
{
    #[inline]
    fn arity() -> usize {
        3
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_tuple3_queues::<A, B, C>(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

impl<A, B, C, D> StepOutputs for (A, B, C, D)
where
    A: Send + HeapSize + 'static,
    B: Send + HeapSize + 'static,
    C: Send + HeapSize + 'static,
    D: Send + HeapSize + 'static,
{
    #[inline]
    fn arity() -> usize {
        4
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_tuple4_queues::<A, B, C, D>(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

impl StepOutputs for () {
    #[inline]
    fn arity() -> usize {
        0
    }

    fn build_queues(
        specs: &[QueueSpec],
        ordering: &[BranchOrdering],
        level: crate::builder::InstrumentationLevel,
    ) -> (OutputQueueSet, OutputsViewAny) {
        super::handles::build_unit_queues(specs, ordering, level)
    }

    fn mark_all_drained(handles: &super::step::OutputHandles<Self>) {
        handles.mark_all_drained();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_arity_is_one() {
        assert_eq!(<Single<u32> as StepOutputs>::arity(), 1);
    }

    #[test]
    fn tuple_2_arity_is_two() {
        assert_eq!(<(u32, u64) as StepOutputs>::arity(), 2);
    }

    #[test]
    fn tuple_3_arity_is_three() {
        assert_eq!(<(u32, u64, String) as StepOutputs>::arity(), 3);
    }

    #[test]
    fn tuple_4_arity_is_four() {
        assert_eq!(<(u32, u64, String, Vec<u8>) as StepOutputs>::arity(), 4);
    }

    #[test]
    fn unit_arity_is_zero() {
        assert_eq!(<() as StepOutputs>::arity(), 0);
    }
}
