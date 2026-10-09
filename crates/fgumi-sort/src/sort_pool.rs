//! A lazily built, clone-shared, bounded rayon pool for the per-run sort.

use std::sync::{Arc, OnceLock};

/// A bounded rayon pool for the per-run sort, built on first use and shared by
/// every clone.
///
/// The pool exists so the run sort uses at most the phase-1 thread count
/// instead of the global rayon pool, which would oversubscribe the pipeline's
/// own workers during a spill. When `--sort-threads` caps phase 1, the caller
/// holds the whole cap while a sort runs on this pool. The width is only known
/// at the first seal — after worker copies exist — so the holder is an
/// `Arc<OnceLock<_>>`: copies share one slot and whichever seals first builds
/// the one pool. **First call wins**: a later `bounded(n)` with a different `n`
/// returns the existing pool. The width is a property of the sort, not of a call.
#[derive(Clone)]
pub struct BoundedSortPool {
    pool: Arc<OnceLock<rayon::ThreadPool>>,
    name: &'static str,
}

impl BoundedSortPool {
    /// An unbuilt pool whose threads will be named `"{name}-{i}"`.
    #[must_use]
    pub fn new(name: &'static str) -> Self {
        Self { pool: Arc::new(OnceLock::new()), name }
    }

    /// The pool, built on first use at `sort_threads.max(1)` threads.
    ///
    /// # Panics
    ///
    /// Panics if the OS refuses to create the threads (an infrastructure failure).
    #[must_use]
    pub fn bounded(&self, sort_threads: usize) -> &rayon::ThreadPool {
        let name = self.name;
        self.pool.get_or_init(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(sort_threads.max(1))
                .thread_name(move |i| format!("{name}-{i}"))
                .build()
                .expect("build bounded sort rayon pool")
        })
    }

    /// The built pool's thread count, or `None` before the first
    /// [`bounded`](Self::bounded) call builds it. Never builds the pool.
    #[must_use]
    pub fn built_threads(&self) -> Option<usize> {
        self.pool.get().map(rayon::ThreadPool::current_num_threads)
    }
}

#[cfg(test)]
mod tests {
    use super::BoundedSortPool;

    /// Clones share ONE pool (pointer identity) — the oversubscription guard
    /// `TemplateArenaAccumulator::sort_pool` documented, now a reusable type.
    #[test]
    fn clones_share_one_pool_built_on_first_use() {
        let a = BoundedSortPool::new("tmpl-sort");
        let b = a.clone();
        let pa = std::ptr::from_ref(a.bounded(3));
        let pb = std::ptr::from_ref(b.bounded(3));
        assert!(std::ptr::eq(pa, pb));
        assert_eq!(a.bounded(3).current_num_threads(), 3);
    }

    /// First call wins: a later, different width gets the original pool.
    #[test]
    fn first_width_wins() {
        let p = BoundedSortPool::new("qname-sort");
        assert_eq!(p.bounded(2).current_num_threads(), 2);
        assert_eq!(p.bounded(7).current_num_threads(), 2);
    }

    /// `built_threads` observes the pool without building it.
    #[test]
    fn built_threads_is_none_until_the_first_build() {
        let p = BoundedSortPool::new("x");
        assert_eq!(p.built_threads(), None);
        let _ = p.bounded(3);
        assert_eq!(p.clone().built_threads(), Some(3), "clones observe the shared pool");
    }

    /// Zero is clamped to one thread, never a zero-width pool.
    #[test]
    fn zero_threads_is_one() {
        assert_eq!(BoundedSortPool::new("x").bounded(0).current_num_threads(), 1);
    }

    /// Threads carry the order's name so `top -H` / profiles attribute them.
    #[test]
    fn threads_are_named_after_the_pool() {
        let p = BoundedSortPool::new("tmpl-sort");
        let name = p.bounded(1).install(|| std::thread::current().name().map(str::to_owned));
        assert_eq!(name.as_deref(), Some("tmpl-sort-0"));
    }
}
