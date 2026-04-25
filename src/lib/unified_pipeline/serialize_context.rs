//! Per-thread context published by the BAM pipeline's serialize step.
//!
//! When `try_step_serialize` invokes the user-supplied `serialize_fn` for each
//! item in a popped batch, it stores the batch's pipeline serial number in
//! this thread-local cell. Closures that need to allocate IDs in serial order
//! — notably `fgumi group` and `fgumi dedup`, via
//! [`crate::ordered_mi_allocator::OrderedMiAllocator`] — read the cell to
//! know which serial they belong to. Combined with
//! [`crate::unified_pipeline::serial_ordered_array_queue::SerialOrderedArrayQueue`],
//! which makes popping happen in serial order, this gives the closure a
//! deterministic serial to anchor its global counter advance against.
//!
//! The harness clears the cell to `None` after each `serialize_fn` call so
//! reads outside a serialize call get `None` and any closure that requires
//! the context can panic loudly rather than silently ignoring missing data.

use std::cell::Cell;

/// Per-batch serialize context: the batch's pipeline serial plus the item's
/// position within the batch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SerializeContext {
    /// Pipeline serial number of the batch currently being serialized.
    pub serial: u64,
    /// Zero-based index of this item within the batch.
    pub item_idx: usize,
    /// Total number of items in the batch.
    pub batch_len: usize,
}

impl SerializeContext {
    /// True when this is the final item in the batch and any per-batch
    /// finalization (e.g. folding per-batch state into a global counter)
    /// should run before the closure returns.
    #[must_use]
    pub fn is_last(&self) -> bool {
        self.item_idx + 1 == self.batch_len
    }
}

thread_local! {
    static CURRENT: Cell<Option<SerializeContext>> = const { Cell::new(None) };
}

/// Stash the serialize context for the calling thread. The harness must call
/// this immediately before invoking `serialize_fn` for each item in a batch
/// and must call [`clear`] after the for-loop completes (success or error).
pub fn set(ctx: SerializeContext) {
    CURRENT.with(|c| c.set(Some(ctx)));
}

/// Clear the serialize context. Read by tests and by the harness on cleanup.
pub fn clear() {
    CURRENT.with(|c| c.set(None));
}

/// Returns the current serialize context for this thread, if any.
///
/// Closures that require the context to be set should `.expect(...)` with a
/// descriptive message.
#[must_use]
pub fn current() -> Option<SerializeContext> {
    CURRENT.with(Cell::get)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unset_returns_none() {
        clear();
        assert!(current().is_none());
    }

    #[test]
    fn set_and_clear_round_trip() {
        let ctx = SerializeContext { serial: 3, item_idx: 1, batch_len: 4 };
        set(ctx);
        assert_eq!(current(), Some(ctx));
        assert!(!ctx.is_last());
        clear();
        assert!(current().is_none());
    }

    #[test]
    fn is_last_is_true_for_final_item() {
        assert!(SerializeContext { serial: 7, item_idx: 2, batch_len: 3 }.is_last());
        assert!(!SerializeContext { serial: 7, item_idx: 1, batch_len: 3 }.is_last());
    }

    #[test]
    fn thread_local_is_isolated_per_thread() {
        clear();
        let other = std::thread::spawn(|| {
            // Thread starts with cleared cell.
            assert!(current().is_none());
            set(SerializeContext { serial: 99, item_idx: 0, batch_len: 1 });
            current()
        });
        // Main thread is unaffected by the child thread's set().
        assert!(current().is_none());
        let other_ctx = other.join().expect("worker");
        assert_eq!(other_ctx.map(|c| c.serial), Some(99));
    }
}
