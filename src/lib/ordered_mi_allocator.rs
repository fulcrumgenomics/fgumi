//! Serial-ordered allocator for cumulative `MoleculeId` numbering.
//!
//! `fgumi group` and `fgumi dedup` allocate contiguous blocks of integer
//! `MoleculeId`s, one block per item within a pipeline batch, by advancing a
//! shared cumulative counter. Even though
//! [`crate::unified_pipeline::serial_ordered_array_queue::SerialOrderedArrayQueue`]
//! makes the pipeline pop processed batches in pipeline-serial order, the
//! `serialize_fn` invocations can still run concurrently across worker
//! threads — one worker can be writing batch `S+1`'s bytes while another is
//! still inside batch `S`'s `serialize_fn`. A naive `AtomicU64::fetch_add`
//! would therefore observe `S+1` advancing the counter before `S` does,
//! producing identical molecule groupings but different integer `MI:Z`
//! numbering across runs.
//!
//! [`OrderedMiAllocator`] serializes the cumulative-offset advance on the
//! pipeline's per-batch serial number. Callers invoke
//! [`allocate`](OrderedMiAllocator::allocate) once per item in a batch with
//! the batch's `serial` and the item's `count`. Multiple items within the
//! same `serial` get back-to-back base values that share the in-progress
//! per-serial offset. Once the last item in batch `serial` is allocated, the
//! caller invokes [`finalize`](OrderedMiAllocator::finalize) to fold the
//! per-serial offset into the cumulative counter and advance the cursor to
//! `serial + 1`. The pipeline guarantees monotonically increasing serials
//! with no gaps, so the cursor always eventually reaches every requested
//! serial — there is no risk of permanent waiting.

use parking_lot::{Condvar, Mutex};
use std::collections::BTreeMap;

/// Allocator that hands out contiguous `MoleculeId` blocks in serial order.
///
/// See the [module docs][self] for the contract callers must follow.
#[derive(Debug, Default)]
pub struct OrderedMiAllocator {
    inner: Mutex<Inner>,
    cv: Condvar,
}

#[derive(Debug, Default)]
struct Inner {
    /// Next batch serial that may begin allocating.
    cursor: u64,
    /// Total count consumed by serials strictly less than `cursor`.
    cumulative: u64,
    /// Per-serial in-progress offset for batches whose `finalize` has not yet
    /// been called. `BTreeMap` rather than `HashMap` so iteration is
    /// trivially deterministic in tests; in production we only ever look up,
    /// insert, and remove by exact key.
    in_progress: BTreeMap<u64, u64>,
}

impl OrderedMiAllocator {
    /// Create a fresh allocator with cursor and cumulative both at 0.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Allocate the next `count` IDs for an item belonging to batch `serial`.
    ///
    /// Blocks until the cursor reaches `serial`, then returns
    /// `cumulative + per_serial_offset` and advances the per-serial offset
    /// by `count`. Subsequent items within the same `serial` can call
    /// `allocate(serial, _)` again — they will not block because the cursor
    /// is already at `serial`, and they will receive bases that pack onto
    /// the earlier ones in the batch.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if `serial` is strictly less than the
    /// allocator's cursor, which would mean an item from an already-finalized
    /// batch is trying to allocate — a pipeline ordering bug.
    pub fn allocate(&self, serial: u64, count: u64) -> u64 {
        let mut g = self.inner.lock();
        while g.cursor < serial {
            self.cv.wait(&mut g);
        }
        debug_assert!(
            g.cursor == serial,
            "OrderedMiAllocator::allocate called for stale serial {serial} (cursor={})",
            g.cursor,
        );
        let cumulative = g.cumulative;
        let offset = g.in_progress.entry(serial).or_insert(0);
        let base = cumulative + *offset;
        *offset += count;
        base
    }

    /// Mark batch `serial` complete.
    ///
    /// Folds the per-serial offset into `cumulative`, advances the cursor,
    /// and wakes any threads waiting on later serials. Must be called
    /// exactly once per serial, after all items in that batch have called
    /// [`allocate`].
    ///
    /// # Panics
    ///
    /// Panics in debug builds if `serial` does not equal the cursor, which
    /// would indicate a missing or extra finalize.
    pub fn finalize(&self, serial: u64) {
        let mut g = self.inner.lock();
        debug_assert_eq!(
            g.cursor, serial,
            "OrderedMiAllocator::finalize called out of order (cursor={}, serial={serial})",
            g.cursor,
        );
        let folded = g.in_progress.remove(&serial).unwrap_or(0);
        g.cumulative += folded;
        g.cursor += 1;
        drop(g);
        self.cv.notify_all();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::thread;

    #[test]
    fn single_item_per_batch_is_consecutive() {
        let alloc = OrderedMiAllocator::new();
        assert_eq!(alloc.allocate(0, 3), 0);
        alloc.finalize(0);
        assert_eq!(alloc.allocate(1, 5), 3);
        alloc.finalize(1);
        assert_eq!(alloc.allocate(2, 2), 8);
        alloc.finalize(2);
    }

    #[test]
    fn multiple_items_share_serial_and_pack_together() {
        let alloc = OrderedMiAllocator::new();
        // Three items in batch 0, each adding 4 IDs.
        assert_eq!(alloc.allocate(0, 4), 0);
        assert_eq!(alloc.allocate(0, 4), 4);
        assert_eq!(alloc.allocate(0, 4), 8);
        alloc.finalize(0);
        // Next batch starts at 12.
        assert_eq!(alloc.allocate(1, 1), 12);
        alloc.finalize(1);
    }

    #[test]
    fn out_of_order_allocate_blocks_until_finalize() {
        let alloc = Arc::new(OrderedMiAllocator::new());

        let alloc_late = Arc::clone(&alloc);
        let late = thread::spawn(move || {
            // Batch 1 wants to allocate; should block until batch 0 finalizes.
            let base = alloc_late.allocate(1, 7);
            alloc_late.finalize(1);
            base
        });

        // Give the late thread a chance to park on the condvar.
        thread::sleep(std::time::Duration::from_millis(50));

        let early_base = alloc.allocate(0, 3);
        alloc.finalize(0);
        assert_eq!(early_base, 0);

        let late_base = late.join().expect("late thread");
        assert_eq!(late_base, 3, "batch 1 must start after batch 0 finalizes");
    }

    #[test]
    fn many_concurrent_batches_get_contiguous_ranges() {
        const N: usize = 64;
        let alloc = Arc::new(OrderedMiAllocator::new());
        let counts: Vec<u64> =
            (0..N).map(|i| u64::try_from(i + 1).expect("count fits u64")).collect();

        // Spawn one thread per batch. They contend on the condvar; the
        // allocator must hand out [0..1), [1..3), [3..6), ... in order.
        let handles: Vec<_> = (0..N)
            .map(|i| {
                let alloc = Arc::clone(&alloc);
                let count = counts[i];
                let serial = u64::try_from(i).expect("serial fits u64");
                thread::spawn(move || {
                    let base = alloc.allocate(serial, count);
                    alloc.finalize(serial);
                    (serial, base, count)
                })
            })
            .collect();

        let mut results: Vec<(u64, u64, u64)> =
            handles.into_iter().map(|h| h.join().expect("worker")).collect();
        results.sort_by_key(|(s, _, _)| *s);

        let mut expected_base = 0u64;
        for (serial, base, count) in results {
            assert_eq!(
                base, expected_base,
                "serial {serial} got base {base}, expected {expected_base}"
            );
            expected_base += count;
        }
    }
}
