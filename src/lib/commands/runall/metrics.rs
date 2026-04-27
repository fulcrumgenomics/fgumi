//! Atomic output guard for the runall temp → final rename dance.

use std::path::PathBuf;

/// Drop guard that removes the temporary output file on error or cancellation.
///
/// Call `disarm()` after a successful rename to prevent cleanup.
pub(super) struct AtomicOutputGuard {
    pub(super) tmp_path: PathBuf,
}

impl AtomicOutputGuard {
    /// Disarm the guard so the temp file is not deleted on drop.
    pub(super) fn disarm(self) {
        std::mem::forget(self);
    }
}

impl Drop for AtomicOutputGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.tmp_path);
    }
}
