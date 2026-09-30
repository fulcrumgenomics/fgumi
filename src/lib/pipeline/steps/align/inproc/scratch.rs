//! One aligner scratch per pool thread, shared by the seed/extend and
//! pair/emit steps.
//!
//! A bwa-mem3 scratch holds ~60 MB of kernel buffers (seed/chain windows,
//! banded-SW and mate-rescue sequence buffers). When each Parallel step's worker
//! copy owned its own, every pool thread kept two, and the one the other step
//! last used went cold. The pool keys scratches by OS thread instead, so a pool
//! worker uses the same warm buffers whichever of the two steps it is running,
//! and they are all freed with the pool when the align chain drops.
//!
//! A pool worker runs one step's `try_run` at a time and pushes never block, so
//! a thread never asks for its scratch while already holding it; if it did, it
//! would get a second scratch rather than alias the first.

use std::collections::HashMap;
use std::io;
use std::sync::Mutex;
use std::thread::{self, ThreadId};

/// Scratches keyed by the OS thread that last used them.
pub(crate) struct ScratchPool<S> {
    idle: Mutex<HashMap<ThreadId, S>>,
}

impl<S> ScratchPool<S> {
    /// An empty pool; scratches are created on first use per thread.
    pub(crate) fn new() -> Self {
        Self { idle: Mutex::new(HashMap::new()) }
    }

    /// Run `f` on this thread's scratch, creating it with `make` on the
    /// thread's first use. The scratch is taken out of the pool for the
    /// duration of `f` (so the lock is not held across the aligner call) and
    /// put back afterwards, whatever `f` returns. A panic in `f` drops it.
    ///
    /// # Errors
    /// Returns `make`'s error, or `f`'s.
    pub(crate) fn with<R>(
        &self,
        make: impl FnOnce() -> io::Result<S>,
        f: impl FnOnce(&mut S) -> io::Result<R>,
    ) -> io::Result<R> {
        let id = thread::current().id();
        let taken = self.lock().remove(&id);
        let mut scratch = match taken {
            Some(scratch) => scratch,
            None => make()?,
        };
        let result = f(&mut scratch);
        self.lock().insert(id, scratch);
        result
    }

    fn lock(&self) -> std::sync::MutexGuard<'_, HashMap<ThreadId, S>> {
        // A poisoned map only means some `f` panicked; the map itself is intact.
        self.idle.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
    }
}

impl<S> Default for ScratchPool<S> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicU64, Ordering};

    use rstest::rstest;

    use super::ScratchPool;

    /// A scratch that records which creation it came from.
    struct Tagged(u64);

    fn tagged_id(pool: &ScratchPool<Tagged>, created: &AtomicU64) -> u64 {
        pool.with(|| Ok(Tagged(created.fetch_add(1, Ordering::Relaxed))), |scratch| Ok(scratch.0))
            .unwrap()
    }

    #[rstest]
    #[case::same_thread_reuses_its_scratch(true, true)]
    #[case::another_thread_gets_its_own(false, false)]
    fn scratch_is_per_thread(#[case] same_thread: bool, #[case] expect_same: bool) {
        let pool = Arc::new(ScratchPool::new());
        let created = Arc::new(AtomicU64::new(0));
        let first = tagged_id(&pool, &created);
        let second = if same_thread {
            tagged_id(&pool, &created)
        } else {
            let (pool, created) = (Arc::clone(&pool), Arc::clone(&created));
            std::thread::spawn(move || tagged_id(&pool, &created)).join().unwrap()
        };
        assert_eq!(first == second, expect_same, "first={first} second={second}");
    }

    /// A nested use on one thread (which the pool's single-`try_run` workers
    /// never do) gets a scratch of its own rather than aliasing the outer one.
    #[test]
    fn nested_use_gets_its_own_scratch() {
        let pool = ScratchPool::new();
        let created = AtomicU64::new(0);
        let make = || Ok(Tagged(created.fetch_add(1, Ordering::Relaxed)));
        let (outer, inner) = pool
            .with(make, |outer| {
                let inner = pool.with(make, |inner| Ok(inner.0))?;
                Ok((outer.0, inner))
            })
            .unwrap();
        assert_ne!(outer, inner);
    }

    /// A creation failure is reported to the caller and leaves nothing cached.
    #[test]
    fn a_failed_creation_is_an_error() {
        let pool: ScratchPool<Tagged> = ScratchPool::new();
        let err = pool.with(|| Err(std::io::Error::other("no memory")), |_| Ok(())).unwrap_err();
        assert!(err.to_string().contains("no memory"), "{err}");
        let created = AtomicU64::new(7);
        assert_eq!(tagged_id(&pool, &created), 7, "the next use creates afresh");
    }
}
