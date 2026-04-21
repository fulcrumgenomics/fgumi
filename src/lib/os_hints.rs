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

/// Hint to the kernel that `len` bytes starting at `offset` on raw file
/// descriptor `fd` will be needed soon. On Linux this triggers asynchronous
/// page-in via `posix_fadvise(POSIX_FADV_WILLNEED)` — the kernel starts I/O
/// immediately and returns without blocking. On other platforms this is a
/// no-op.
///
/// This is the raw-fd variant for use inside the [`PrefetchReader`] producer
/// thread, which owns the `File` and cannot lend a reference back. The fd
/// remains valid because the producer thread keeps the `File` alive for the
/// duration of its loop.
///
/// [`PrefetchReader`]: crate::prefetch_reader::PrefetchReader
pub fn advise_willneed_raw(fd: i32, offset: i64, len: i64) {
    // Intentionally takes a bare `i32` rather than `RawFd` so callers outside
    // `cfg(unix)` can still compile without importing platform-specific types.
    #[cfg(target_os = "linux")]
    linux::advise_willneed_raw(fd, offset, len);
    #[cfg(not(target_os = "linux"))]
    let _ = (fd, offset, len);
}

/// Extract the raw file descriptor from `file` on Linux so it can be passed to
/// [`advise_willneed_raw`] from a thread that owns the `File`. Returns `None`
/// on non-Linux platforms where `posix_fadvise` is unavailable.
#[must_use]
#[cfg(target_os = "linux")]
pub fn hint_fd(file: &File) -> Option<i32> {
    use std::os::fd::AsRawFd;
    Some(file.as_raw_fd())
}

/// Non-Linux stub — always returns `None`.
#[must_use]
#[cfg(not(target_os = "linux"))]
pub fn hint_fd(_file: &File) -> Option<i32> {
    None
}

#[cfg(target_os = "linux")]
mod linux {
    use std::fs::File;
    use std::os::fd::AsRawFd;

    use nix::fcntl::{PosixFadviseAdvice, posix_fadvise};

    pub(super) fn advise_sequential(file: &File) {
        // offset = 0, len = 0 applies the hint to the entire file (current and
        // future contents). Errors are non-fatal — log and move on.
        let fd = file.as_raw_fd();
        if let Err(e) = posix_fadvise(fd, 0, 0, PosixFadviseAdvice::POSIX_FADV_SEQUENTIAL) {
            tracing::debug!("posix_fadvise(POSIX_FADV_SEQUENTIAL) failed on fd {fd}: {e}");
        }
    }

    pub(super) fn advise_willneed_raw(fd: i32, offset: i64, len: i64) {
        if let Err(e) = posix_fadvise(fd, offset, len, PosixFadviseAdvice::POSIX_FADV_WILLNEED) {
            tracing::debug!("posix_fadvise(POSIX_FADV_WILLNEED) failed on fd {fd}: {e}");
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
