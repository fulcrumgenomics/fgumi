//! `StepOutputs`: type-level description of a step's outputs.
//!
//! Single-output steps declare `type Outputs = Single<T>;`.
//! Multi-output steps declare `type Outputs = (A, B, C);` (positional access).
//! Sinks declare `type Outputs = ();`.
//!
//! **Maximum arity: 4.** [`MAX_ARITY`] is the largest tuple shape
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

/// The largest tuple `Outputs` shape that has a `StepOutputs` impl.
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

/// Ordered + byte-bounded tuple-3 outputs. The 3-branch analog of
/// [`OrderedBytesTuple2`]: use when a step fans out to three branches that all
/// need `BranchOrdering::ByItemOrdinal` + byte-bounded queues — e.g. paired
/// FASTQ output splitting one record batch into R1 / R2 / other byte streams,
/// each feeding an ordered `WriteRawFile` sink.
///
/// The plain `(A, B, C)` `StepOutputs` impl only bounds each branch `HeapSize`.
/// To use `ByItemOrdinal` on any branch, all three must satisfy
/// `Ordered + HeapSize` — that's what this shape encodes at the type level.
#[allow(clippy::type_complexity)]
pub struct OrderedBytesTuple3<A, B, C>(PhantomData<fn() -> (A, B, C)>)
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    C: Send + HeapSize + Ordered + 'static;

impl<A, B, C> StepOutputs for OrderedBytesTuple3<A, B, C>
where
    A: Send + HeapSize + Ordered + 'static,
    B: Send + HeapSize + Ordered + 'static,
    C: Send + HeapSize + Ordered + 'static,
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
        super::handles::build_tuple3_queues_ordered_bytes::<A, B, C>(specs, ordering, level)
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
    use crate::step::InputHandle;

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

    #[derive(Clone, Copy)]
    struct OrdU64(u64);
    impl crate::item::HeapSize for OrdU64 {}
    impl crate::item::Ordered for OrdU64 {
        fn ordinal(&self) -> u64 {
            self.0
        }
    }

    #[test]
    fn ordered_bytes_tuple_3_arity_is_three() {
        // OrderedBytesTuple3 is the ordered + byte-bounded 3-way fan-out shape.
        // Its three branches carry Ordered + HeapSize items (here the u64/u32
        // stand-ins just need to satisfy the bounds at the type level).
        fn assert_arity<O: StepOutputs>() -> usize {
            O::arity()
        }
        assert_eq!(assert_arity::<OrderedBytesTuple3<OrdU64, OrdU64, OrdU64>>(), 3);
    }

    #[test]
    fn tuple_4_arity_is_four() {
        assert_eq!(<(u32, u64, String, Vec<u8>) as StepOutputs>::arity(), 4);
    }

    #[test]
    fn unit_arity_is_zero() {
        assert_eq!(<() as StepOutputs>::arity(), 0);
    }

    #[test]
    fn ordered_bytes_tuple_2_arity_is_two() {
        fn assert_arity<O: StepOutputs>() -> usize {
            O::arity()
        }
        assert_eq!(assert_arity::<OrderedBytesTuple2<OrdU64, OrdU64>>(), 2);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // `build_queues` / `mark_all_drained` round-trips, one test per output shape.
    //
    // Every branch a shape declares must come back from `build_queues` as its
    // own live edge: `n_branches()` matches `arity()`, each branch's input handle
    // transports the item type its position declares, each starts OPEN, and
    // `mark_all_drained` closes all of them. The open-before / closed-after pair
    // is what makes these discriminating — asserting only the closed-after state
    // would pass just as well if `build_queues` returned branches already closed
    // and `mark_all_drained` did nothing.
    //
    // These stay separate tests rather than an `#[rstest]` case table: each shape
    // is a distinct type with distinct per-branch item types, so the cases cannot
    // share a function signature. Every shape with a `StepOutputs` impl has one —
    // `Single`, `OrderedBytesSingle`, the 2-/3-/4-tuples, `OrderedBytesTuple2`,
    // `OrderedBytesTuple3`, and `()`. Adding a shape without adding its
    // round-trip leaves its `build_queues` / `mark_all_drained` pair unexercised.
    // ─────────────────────────────────────────────────────────────────────────

    fn specs(n: usize) -> Vec<QueueSpec> {
        vec![QueueSpec::CountBounded { capacity: 4 }; n]
    }

    fn orderings(n: usize) -> Vec<BranchOrdering> {
        vec![BranchOrdering::None; n]
    }

    /// Build a shape's queues and hand back the branch set plus a typed
    /// `OutputHandles` view, the pair the runtime hands to a step.
    fn build<O: StepOutputs>(n: usize) -> (OutputQueueSet, crate::step::OutputHandles<O>) {
        let (queues, view) =
            O::build_queues(&specs(n), &orderings(n), crate::builder::InstrumentationLevel::Off);
        (queues, crate::step::OutputHandles::<O>::new(view))
    }

    #[test]
    fn single_builds_one_branch_and_drains() {
        type Shape = Single<u32>;
        let (mut queues, handles) = build::<Shape>(1);
        assert_eq!(queues.n_branches(), 1, "one branch per declared output");

        let a = queues.take_typed_input::<u32>(0);
        assert!(!a.is_drained(), "branch starts open");

        <Shape as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
    }

    /// The only shape whose `build_queues` routes to
    /// `build_single_queues_ordered_bytes`, so without this its drain path has no
    /// coverage here at all.
    #[test]
    fn ordered_bytes_single_builds_one_branch_and_drains() {
        type Shape = OrderedBytesSingle<OrdU64>;
        let (mut queues, handles) = build::<Shape>(1);
        assert_eq!(queues.n_branches(), 1, "one branch per declared output");

        let a = queues.take_typed_input::<OrdU64>(0);
        assert!(!a.is_drained(), "branch starts open");

        <Shape as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
    }

    #[test]
    fn tuple_2_builds_two_independent_branches() {
        type Shape = (u32, u64);
        let (mut queues, handles) = build::<Shape>(2);
        assert_eq!(queues.n_branches(), 2, "one branch per declared output");

        let a = queues.take_typed_input::<u32>(0);
        let b = queues.take_typed_input::<u64>(1);
        assert!(!a.is_drained(), "branches start open");
        assert!(!b.is_drained());

        <Shape as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
        assert!(b.is_drained(), "mark_all_drained closes branch 1");
    }

    #[test]
    fn tuple_3_builds_three_independent_branches() {
        let (mut queues, handles) = build::<(u32, u64, String)>(3);
        assert_eq!(queues.n_branches(), 3, "one branch per declared output");

        // Each branch's input handle must downcast to that position's type —
        // a mis-wired builder would panic here or hand back the wrong branch.
        let a = queues.take_typed_input::<u32>(0);
        let b = queues.take_typed_input::<u64>(1);
        let c = queues.take_typed_input::<String>(2);
        assert!(!a.is_drained(), "branches start open");
        assert!(!b.is_drained());
        assert!(!c.is_drained());

        <(u32, u64, String) as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
        assert!(b.is_drained(), "mark_all_drained closes branch 1");
        assert!(c.is_drained(), "mark_all_drained closes branch 2");
    }

    #[test]
    fn tuple_4_builds_four_independent_branches() {
        let (mut queues, handles) = build::<(u32, u64, String, Vec<u8>)>(4);
        assert_eq!(queues.n_branches(), 4, "one branch per declared output");

        let a = queues.take_typed_input::<u32>(0);
        let b = queues.take_typed_input::<u64>(1);
        let c = queues.take_typed_input::<String>(2);
        let d = queues.take_typed_input::<Vec<u8>>(3);
        assert!(!a.is_drained(), "branches start open");
        assert!(!b.is_drained());
        assert!(!c.is_drained());
        assert!(!d.is_drained());

        <(u32, u64, String, Vec<u8>) as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
        assert!(b.is_drained(), "mark_all_drained closes branch 1");
        assert!(c.is_drained(), "mark_all_drained closes branch 2");
        assert!(d.is_drained(), "mark_all_drained closes branch 3");
    }

    #[test]
    fn ordered_bytes_tuple_2_builds_two_independent_branches() {
        type Shape = OrderedBytesTuple2<OrdU64, OrdU64>;
        let (mut queues, handles) = build::<Shape>(2);
        assert_eq!(queues.n_branches(), 2, "one branch per declared output");

        let a = queues.take_typed_input::<OrdU64>(0);
        let b = queues.take_typed_input::<OrdU64>(1);
        assert!(!a.is_drained(), "branches start open");
        assert!(!b.is_drained());

        // Identity, not just type. Both branches carry the SAME item type, so the
        // per-position downcast cannot tell them apart, and `mark_all_drained`
        // closes every branch — a builder that aliased branch 0's queue onto both
        // positions would satisfy every other assertion here. Route a distinct
        // value through each and read it back off its own edge.
        let view = handles.view();
        view.a.push(OrdU64(10)).expect("branch 0 accepts one item");
        view.b.push(OrdU64(11)).expect("branch 1 accepts one item");
        assert_eq!(a.pop().map(|v| v.0), Some(10), "branch 0 is its own edge");
        assert_eq!(b.pop().map(|v| v.0), Some(11), "branch 1 is its own edge");

        <Shape as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
        assert!(b.is_drained(), "mark_all_drained closes branch 1");
    }

    #[test]
    fn ordered_bytes_tuple_3_builds_three_independent_branches() {
        type Shape = OrderedBytesTuple3<OrdU64, OrdU64, OrdU64>;
        let (mut queues, handles) = build::<Shape>(3);
        assert_eq!(queues.n_branches(), 3, "one branch per declared output");

        let a = queues.take_typed_input::<OrdU64>(0);
        let b = queues.take_typed_input::<OrdU64>(1);
        let c = queues.take_typed_input::<OrdU64>(2);
        assert!(!a.is_drained(), "branches start open");
        assert!(!b.is_drained());
        assert!(!c.is_drained());

        // Identity, not just type — see the tuple-2 test above for why the
        // same-item-type shapes need this.
        let view = handles.view();
        view.a.push(OrdU64(10)).expect("branch 0 accepts one item");
        view.b.push(OrdU64(11)).expect("branch 1 accepts one item");
        view.c.push(OrdU64(12)).expect("branch 2 accepts one item");
        assert_eq!(a.pop().map(|v| v.0), Some(10), "branch 0 is its own edge");
        assert_eq!(b.pop().map(|v| v.0), Some(11), "branch 1 is its own edge");
        assert_eq!(c.pop().map(|v| v.0), Some(12), "branch 2 is its own edge");

        <Shape as StepOutputs>::mark_all_drained(&handles);
        assert!(a.is_drained(), "mark_all_drained closes branch 0");
        assert!(b.is_drained(), "mark_all_drained closes branch 1");
        assert!(c.is_drained(), "mark_all_drained closes branch 2");
    }

    /// A sink declares no outputs, so its queue set is empty and
    /// `mark_all_drained` is a no-op rather than a panic.
    #[test]
    fn unit_builds_no_branches_and_drains_without_panicking() {
        let (queues, handles) = build::<()>(0);
        assert_eq!(queues.n_branches(), 0, "a sink owns no output branches");
        <() as StepOutputs>::mark_all_drained(&handles);
    }
}
