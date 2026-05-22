#![deny(unsafe_code)]
//! OS-specific I/O hints.
//!
//! This module provides small, portable wrappers around OS-specific advisory
//! calls. The primary use case is informing the kernel about a file's expected
//! access pattern before streaming through it, which on Linux can enlarge the
//! per-fd read-ahead window well beyond the 128 KiB default set by
//! `/sys/block/*/queue/read_ahead_kb`.
//!
//! All functions are best-effort: they are no-ops on non-supported platforms
//! and never panic. Errors are logged at `debug` level and swallowed so that
//! I/O hints can be wired in at the open site without error-handling
//! boilerplate at every call.

use std::fs::File;

/// Hint to the kernel that `file` will be read from start to finish in a
/// single forward pass. On Linux this enlarges the per-fd read-ahead window;
/// on other platforms this is a no-op.
///
/// Failure is logged at `debug` level and otherwise ignored — the file
/// remains usable, the hint simply wasn't applied.
pub fn advise_sequential(file: &File) {
    #[cfg(target_os = "linux")]
    linux::advise_sequential(file);
    #[cfg(not(target_os = "linux"))]
    let _ = file;
}

/// Hint to the kernel that `len` bytes starting at `offset` on `file` will be
/// needed soon. On Linux this triggers asynchronous page-in via
/// `posix_fadvise(POSIX_FADV_WILLNEED)` — the kernel starts I/O immediately
/// and returns without blocking. On other platforms this is a no-op.
///
/// Used inside the [`PrefetchReader`] producer thread, which owns a clone of
/// the underlying file and uses it solely to issue these hints.
///
/// [`PrefetchReader`]: crate::prefetch_reader::PrefetchReader
pub fn advise_willneed(file: &File, offset: i64, len: i64) {
    #[cfg(target_os = "linux")]
    linux::advise_willneed(file, offset, len);
    #[cfg(not(target_os = "linux"))]
    let _ = (file, offset, len);
}

/// Clone `file` for use as a separate fd that the prefetch producer can keep
/// alive solely to issue `posix_fadvise` hints against. Returns `None` on
/// non-Linux platforms (where the hint is a no-op) and on Linux if the clone
/// fails — the reader still works, the hint is simply skipped.
#[must_use]
#[cfg(target_os = "linux")]
pub fn hint_file(file: &File) -> Option<File> {
    file.try_clone().ok()
}

/// Non-Linux stub — always returns `None`.
#[must_use]
#[cfg(not(target_os = "linux"))]
pub fn hint_file(_file: &File) -> Option<File> {
    None
}

#[cfg(target_os = "linux")]
mod linux {
    use std::fs::File;

    use nix::fcntl::{PosixFadviseAdvice, posix_fadvise};

    pub(super) fn advise_sequential(file: &File) {
        // offset = 0, len = 0 applies the hint to the entire file (current and
        // future contents). Errors are non-fatal — log and move on.
        if let Err(e) = posix_fadvise(file, 0, 0, PosixFadviseAdvice::POSIX_FADV_SEQUENTIAL) {
            log::debug!("posix_fadvise(POSIX_FADV_SEQUENTIAL) failed: {e}");
        }
    }

    pub(super) fn advise_willneed(file: &File, offset: i64, len: i64) {
        if let Err(e) = posix_fadvise(file, offset, len, PosixFadviseAdvice::POSIX_FADV_WILLNEED) {
            log::debug!("posix_fadvise(POSIX_FADV_WILLNEED) failed: {e}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn advise_sequential_does_not_error_on_regular_file() {
        let mut tmp = NamedTempFile::new().expect("create temp file");
        tmp.write_all(b"hello world").expect("write temp file");
        let f = File::open(tmp.path()).expect("reopen temp file for read");
        // The call should succeed or be a no-op. It must not panic.
        advise_sequential(&f);
    }

    #[test]
    fn advise_sequential_is_callable_repeatedly() {
        let tmp = NamedTempFile::new().expect("create temp file");
        let f = File::open(tmp.path()).expect("open temp file");
        for _ in 0..16 {
            advise_sequential(&f);
        }
    }
}
