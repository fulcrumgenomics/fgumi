//! Allocator for distributing sort spill files across multiple temporary directories.
//!
//! `TmpDirAllocator` hands out base directory paths for spill/merged files in a
//! round-robin fashion across one or more user-supplied temp directories. It
//! also supports:
//!
//! - Dropping a directory from the rotation when it runs out of space
//!   (e.g. `ENOSPC` during chunk creation).
//! - Periodic re-checking of free space; directories whose free space falls
//!   below a minimum threshold are dropped from rotation automatically.

use anyhow::{Result, anyhow, bail};
use std::path::{Path, PathBuf};

/// Minimum free bytes a directory must have to remain in rotation.
pub const DEFAULT_MIN_FREE_BYTES: u64 = 1 << 30; // 1 GiB

/// Number of `next()` calls between periodic free-space re-checks.
pub const DEFAULT_RECHECK_INTERVAL: usize = 8;

/// Function that reports free bytes available on a filesystem for writes.
///
/// The default implementation uses [`fs4::available_space`]. Tests inject a
/// closure so they can simulate a filesystem filling up.
pub type FreeSpaceProbe = Box<dyn FnMut(&Path) -> Result<u64> + Send>;

/// Allocator that distributes spill directory assignments across one or more
/// base directories using free-space-aware round-robin.
pub struct TmpDirAllocator {
    /// Directories currently eligible for allocation. Shrinks as dirs are
    /// dropped via `mark_full` or the periodic free-space recheck.
    active: Vec<PathBuf>,
    cursor: usize,
    min_free_bytes: u64,
    recheck_interval: usize,
    calls_since_recheck: usize,
    probe: FreeSpaceProbe,
}

impl TmpDirAllocator {
    /// Create an allocator over the given base directories.
    ///
    /// At least one directory must be supplied. Directories whose initial free
    /// space is below `DEFAULT_MIN_FREE_BYTES` are dropped immediately.
    ///
    /// # Errors
    ///
    /// Returns an error if `dirs` is empty or if no directory has sufficient
    /// free space.
    pub fn new(dirs: Vec<PathBuf>) -> Result<Self> {
        Self::with_probe(dirs, Box::new(default_free_space_probe), DEFAULT_MIN_FREE_BYTES)
    }

    /// Create an allocator with a custom free-space probe and threshold
    /// (primarily for testing).
    ///
    /// `min_free_bytes` is the single threshold used both for the startup
    /// filter and the periodic recheck.
    ///
    /// # Errors
    ///
    /// Returns an error if `dirs` is empty or if every supplied directory
    /// either fails to probe or has insufficient free space.
    pub fn with_probe(
        dirs: Vec<PathBuf>,
        mut probe: FreeSpaceProbe,
        min_free_bytes: u64,
    ) -> Result<Self> {
        if dirs.is_empty() {
            bail!("TmpDirAllocator requires at least one directory");
        }

        let mut active = Vec::with_capacity(dirs.len());
        for dir in dirs {
            match probe(&dir) {
                Ok(free) if free >= min_free_bytes => active.push(dir),
                Ok(free) => tracing::warn!(
                    "Temp dir {} dropped: only {} free (need {})",
                    dir.display(),
                    free,
                    min_free_bytes
                ),
                Err(e) => tracing::warn!("Temp dir {} dropped: probe failed: {}", dir.display(), e),
            }
        }

        if active.is_empty() {
            bail!("No temp directory has sufficient free space (need {min_free_bytes} bytes)");
        }

        Ok(Self {
            active,
            cursor: 0,
            min_free_bytes,
            recheck_interval: DEFAULT_RECHECK_INTERVAL,
            calls_since_recheck: 0,
            probe,
        })
    }

    /// Override the periodic recheck interval (primarily for testing).
    #[must_use]
    pub fn with_recheck_interval(mut self, interval: usize) -> Self {
        self.recheck_interval = interval.max(1);
        self
    }

    /// Number of directories currently eligible for allocation.
    #[cfg(test)]
    #[must_use]
    pub fn active_count(&self) -> usize {
        self.active.len()
    }

    /// Return the next base directory in round-robin order.
    ///
    /// # Errors
    ///
    /// Returns an error if all directories have been marked full or dropped
    /// by free-space checks.
    #[allow(clippy::should_implement_trait)] // This is not Iterator::next — it is fallible.
    pub fn next(&mut self) -> Result<PathBuf> {
        self.maybe_recheck();

        if self.active.is_empty() {
            return Err(anyhow!(
                "All temp directories have been exhausted (marked full or below min free space)"
            ));
        }

        // `cursor` is maintained in `[0, active.len())` by every mutation of
        // `active`, so indexing is safe without a modulo.
        let dir = self.active[self.cursor].clone();
        self.cursor = (self.cursor + 1) % self.active.len();
        self.calls_since_recheck = self.calls_since_recheck.saturating_add(1);
        Ok(dir)
    }

    /// Drop a directory from rotation (e.g. after `ENOSPC` during a spill write).
    ///
    /// Matches by path equality. A no-op if the path isn't currently active.
    pub fn mark_full(&mut self, dir: &Path) {
        if let Some(pos) = self.active.iter().position(|d| d == dir) {
            self.active.remove(pos);
            self.rebase_cursor_after_remove(pos);
            tracing::warn!("Temp dir {} removed from rotation", dir.display());
        }
    }

    /// Shift `cursor` to compensate for an element removed at `removed_pos`,
    /// preserving round-robin order.
    ///
    /// If the removed index was strictly before the cursor, everything after
    /// it shifted down by one so the cursor must follow. If the cursor now
    /// sits past the end of `active`, wrap it back to the front. When
    /// `active` becomes empty the cursor is reset to zero so the next
    /// `next()` call fails cleanly on the length check rather than
    /// mis-indexing.
    fn rebase_cursor_after_remove(&mut self, removed_pos: usize) {
        if self.active.is_empty() {
            self.cursor = 0;
            return;
        }
        if removed_pos < self.cursor {
            self.cursor -= 1;
        }
        if self.cursor >= self.active.len() {
            self.cursor = 0;
        }
    }

    fn maybe_recheck(&mut self) {
        if self.calls_since_recheck < self.recheck_interval {
            return;
        }
        self.calls_since_recheck = 0;

        let min = self.min_free_bytes;
        let mut idx = 0;
        while idx < self.active.len() {
            let drop_dir = match (self.probe)(&self.active[idx]) {
                Ok(free) if free >= min => false,
                Ok(free) => {
                    tracing::warn!(
                        "Temp dir {} dropped mid-sort: {} free (need {})",
                        self.active[idx].display(),
                        free,
                        min
                    );
                    true
                }
                Err(e) => {
                    tracing::warn!(
                        "Temp dir {} dropped mid-sort: probe failed: {}",
                        self.active[idx].display(),
                        e
                    );
                    true
                }
            };
            if drop_dir {
                self.active.remove(idx);
                self.rebase_cursor_after_remove(idx);
                // Do not advance `idx`: what was at `idx+1` is now at `idx`.
            } else {
                idx += 1;
            }
        }
    }
}

/// Default probe: number of bytes free on the filesystem backing `path`.
fn default_free_space_probe(path: &Path) -> Result<u64> {
    fs4::available_space(path).map_err(anyhow::Error::from)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    /// Probe that always returns a fixed free-space value.
    fn always_free(bytes: u64) -> FreeSpaceProbe {
        Box::new(move |_| Ok(bytes))
    }

    /// Probe backed by a shared map of path -> free-bytes.
    ///
    /// Used to simulate individual directories filling up between calls.
    fn map_probe(state: Arc<Mutex<std::collections::HashMap<PathBuf, u64>>>) -> FreeSpaceProbe {
        Box::new(move |p: &Path| {
            let guard = state.lock().expect("probe map lock");
            Ok(*guard.get(p).unwrap_or(&u64::MAX))
        })
    }

    #[test]
    fn test_new_empty_errors() {
        let result = TmpDirAllocator::new(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_single_dir_always_returns_same() {
        let dir = PathBuf::from("/tmp/a");
        let mut alloc = TmpDirAllocator::with_probe(
            vec![dir.clone()],
            always_free(u64::MAX),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("allocator should be constructable with one dir");
        for _ in 0..5 {
            assert_eq!(alloc.next().expect("next should succeed"), dir);
        }
    }

    #[test]
    fn test_two_dirs_round_robin() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");
        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone()],
            always_free(u64::MAX),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("two-dir alloc should build");

        let seq: Vec<_> = (0..4).map(|_| alloc.next().expect("next")).collect();
        assert_eq!(seq, vec![a.clone(), b.clone(), a.clone(), b]);
    }

    #[test]
    fn test_mark_full_skips_dir() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");
        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone()],
            always_free(u64::MAX),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("two-dir alloc should build");

        assert_eq!(alloc.next().expect("next"), a);
        alloc.mark_full(&b);
        // After dropping b, every subsequent next() must be a.
        for _ in 0..5 {
            assert_eq!(alloc.next().expect("next"), a);
        }
        assert_eq!(alloc.active_count(), 1);
    }

    #[test]
    fn test_mark_full_all_exhausts() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");
        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone()],
            always_free(u64::MAX),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("two-dir alloc should build");

        alloc.mark_full(&a);
        alloc.mark_full(&b);
        assert!(alloc.next().is_err(), "exhausted allocator must error");
    }

    /// Regression: removing a dir at a position strictly before the cursor
    /// must shift the cursor down so round-robin order is preserved.
    #[test]
    fn test_mark_full_preserves_round_robin_order() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");
        let c = PathBuf::from("/tmp/c");
        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone(), c.clone()],
            always_free(u64::MAX),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("three-dir alloc should build");

        // After handing out `a`, the next-in-rotation is `b`.
        assert_eq!(alloc.next().expect("1"), a);
        // Pruning `a` (pos=0, strictly before cursor=1) must leave `b` next,
        // not `c`.
        alloc.mark_full(&a);
        assert_eq!(alloc.next().expect("2"), b);
        assert_eq!(alloc.next().expect("3"), c);
        assert_eq!(alloc.next().expect("4"), b);
    }

    /// Regression: pruning the dir the cursor points at lets the cursor
    /// naturally move to what was the next element, and wrapping still
    /// works cleanly when that was the tail.
    #[test]
    fn test_mark_full_at_cursor_wraps_cleanly() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");
        let c = PathBuf::from("/tmp/c");
        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone(), c.clone()],
            always_free(u64::MAX),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("three-dir alloc should build");

        // Hand out a, then b; cursor now points to c (pos=2).
        assert_eq!(alloc.next().expect("1"), a);
        assert_eq!(alloc.next().expect("2"), b);
        // Prune c at pos=2 == cursor. active=[a,b]; cursor=2 is now past the
        // end and must wrap to 0.
        alloc.mark_full(&c);
        assert_eq!(alloc.next().expect("3"), a);
        assert_eq!(alloc.next().expect("4"), b);
    }

    #[test]
    fn test_initial_free_space_filter() {
        // Dir b has insufficient free space at startup → dropped immediately.
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");

        let mut map = std::collections::HashMap::new();
        map.insert(a.clone(), u64::MAX);
        map.insert(b.clone(), 1024); // tiny
        let state = Arc::new(Mutex::new(map));

        let alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone()],
            map_probe(state),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("should accept at least one valid dir");
        assert_eq!(alloc.active_count(), 1);
    }

    /// The startup filter must honor the `min_free_bytes` argument rather
    /// than the `DEFAULT_MIN_FREE_BYTES` constant.
    #[test]
    fn test_startup_filter_honors_threshold_argument() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");

        // Both dirs report 2048 free, well below the 1 GiB default. A low
        // threshold of 1024 must still admit them.
        let alloc =
            TmpDirAllocator::with_probe(vec![a.clone(), b.clone()], always_free(2048), 1024)
                .expect("low threshold should admit both dirs");
        assert_eq!(alloc.active_count(), 2);
    }

    #[test]
    fn test_periodic_recheck_drops_dir() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");

        let mut map = std::collections::HashMap::new();
        map.insert(a.clone(), u64::MAX);
        map.insert(b.clone(), u64::MAX);
        let state = Arc::new(Mutex::new(map));

        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone()],
            map_probe(Arc::clone(&state)),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("two-dir alloc should build")
        .with_recheck_interval(2);

        // First two calls round-robin between a and b.
        assert_eq!(alloc.next().expect("1"), a);
        assert_eq!(alloc.next().expect("2"), b);

        // Now b runs out of space. Recheck fires on next call.
        state.lock().expect("lock").insert(b.clone(), 1024);

        // After recheck, only a should remain.
        let next = alloc.next().expect("3");
        assert_eq!(next, a);
        // And subsequent calls stay on a.
        for _ in 0..4 {
            assert_eq!(alloc.next().expect("more"), a);
        }
        assert_eq!(alloc.active_count(), 1);
    }

    /// Regression: when `maybe_recheck` drops a dir at a position strictly
    /// before the cursor, the cursor must shift with it.
    #[test]
    fn test_periodic_recheck_preserves_round_robin_order() {
        let a = PathBuf::from("/tmp/a");
        let b = PathBuf::from("/tmp/b");
        let c = PathBuf::from("/tmp/c");

        let mut map = std::collections::HashMap::new();
        map.insert(a.clone(), u64::MAX);
        map.insert(b.clone(), u64::MAX);
        map.insert(c.clone(), u64::MAX);
        let state = Arc::new(Mutex::new(map));

        let mut alloc = TmpDirAllocator::with_probe(
            vec![a.clone(), b.clone(), c.clone()],
            map_probe(Arc::clone(&state)),
            DEFAULT_MIN_FREE_BYTES,
        )
        .expect("three-dir alloc should build")
        .with_recheck_interval(1);

        // Hand out a; cursor is now 1 pointing at b.
        assert_eq!(alloc.next().expect("1"), a);
        // a runs out of space before the next call.
        state.lock().expect("lock").insert(a.clone(), 1024);

        // Recheck drops a at pos=0 (< cursor=1). Cursor must follow so the
        // next allocation stays on b, not c.
        assert_eq!(alloc.next().expect("2"), b);
        assert_eq!(alloc.next().expect("3"), c);
        assert_eq!(alloc.next().expect("4"), b);
    }

    #[test]
    fn test_all_dirs_below_threshold_errors() {
        let a = PathBuf::from("/tmp/a");
        let result =
            TmpDirAllocator::with_probe(vec![a], always_free(1024), DEFAULT_MIN_FREE_BYTES);
        assert!(result.is_err(), "no dirs with enough free space → error");
    }
}
