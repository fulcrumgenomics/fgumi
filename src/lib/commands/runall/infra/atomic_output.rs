//! Atomic output guard for the runall temp → final rename dance.
//!
//! Pipeline runs write to a `.tmp` sibling of the user's `--output` path.
//! On a clean run, [`AtomicOutputGuard::disarm`] is called immediately
//! before the temp is renamed onto the final path. If the run errors out,
//! is cancelled, or panics, the guard's `Drop` removes the temp file so
//! the user is never left with a partial output.

use std::path::{Path, PathBuf};

/// Drop guard that removes the temporary output file on error or cancellation.
///
/// Construct with [`AtomicOutputGuard::new`] before kicking off the
/// pipeline. After a successful rename, call [`AtomicOutputGuard::disarm`]
/// to consume the guard without removing the temp.
///
/// # Cleanup contract
///
/// The guard always attempts to remove `tmp_path` on drop, swallowing any
/// I/O error (the file may legitimately not exist if the pipeline failed
/// before writing anything). It is the caller's job to ensure no other
/// process is concurrently writing to that path; runall is single-pipeline
/// per-process so this is trivially true.
#[derive(Debug)]
pub(crate) struct AtomicOutputGuard {
    /// Path of the temporary file to clean up on drop.
    tmp_path: PathBuf,
}

impl AtomicOutputGuard {
    /// Construct a new guard for `tmp_path`. The path does not need to
    /// exist yet — the guard simply remembers it for cleanup.
    pub(crate) fn new(tmp_path: PathBuf) -> Self {
        Self { tmp_path }
    }

    /// Borrow the temp path for use by chain runners.
    pub(crate) fn tmp_path(&self) -> &Path {
        &self.tmp_path
    }

    /// Disarm the guard so the temp file is not deleted on drop.
    ///
    /// Call this exactly once, after a successful `rename(tmp, final)`.
    pub(crate) fn disarm(self) {
        std::mem::forget(self);
    }
}

impl Drop for AtomicOutputGuard {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.tmp_path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn guard_removes_tmp_on_drop() {
        let dir = TempDir::new().unwrap();
        let tmp = dir.path().join("output.bam.tmp");
        fs::write(&tmp, b"partial output").unwrap();
        assert!(tmp.exists());

        {
            let _guard = AtomicOutputGuard::new(tmp.clone());
        }

        assert!(!tmp.exists(), "guard's Drop should have removed the tmp file");
    }

    #[test]
    fn guard_drops_silently_when_tmp_does_not_exist() {
        let dir = TempDir::new().unwrap();
        let tmp = dir.path().join("never_created.bam.tmp");
        assert!(!tmp.exists());

        // Drop should not panic even though the file does not exist.
        {
            let _guard = AtomicOutputGuard::new(tmp.clone());
        }

        assert!(!tmp.exists());
    }

    #[test]
    fn disarmed_guard_does_not_remove_tmp() {
        let dir = TempDir::new().unwrap();
        let tmp = dir.path().join("output.bam.tmp");
        fs::write(&tmp, b"finalized output").unwrap();
        assert!(tmp.exists());

        let guard = AtomicOutputGuard::new(tmp.clone());
        guard.disarm();

        assert!(tmp.exists(), "disarmed guard must not remove the tmp file");
    }

    #[test]
    fn tmp_path_accessor_returns_input_path() {
        let dir = TempDir::new().unwrap();
        let tmp = dir.path().join("x.bam.tmp");
        let guard = AtomicOutputGuard::new(tmp.clone());
        assert_eq!(guard.tmp_path(), tmp.as_path());
        // Disarm so the test doesn't leave a phantom remove call (path doesn't exist anyway).
        guard.disarm();
    }
}
