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

/// Hint to the kernel that `len` bytes starting at `offset` on the file behind
/// `fd` will be needed soon. On Linux this triggers asynchronous page-in via
/// `posix_fadvise(POSIX_FADV_WILLNEED)` — the kernel starts I/O immediately and
/// returns without blocking. On other platforms this is a no-op.
///
/// This is the owning-handle variant for use inside the [`PrefetchReader`]
/// producer thread, which owns the source `File` but type-erases it to a
/// generic reader and so cannot lend a `&File` back. The producer instead holds
/// a duplicated descriptor ([`HintFd`]) alive for the lifetime of its loop and
/// passes it here by reference.
///
/// [`PrefetchReader`]: crate::prefetch_reader::PrefetchReader
pub fn advise_willneed(fd: &HintFd, offset: i64, len: i64) {
    #[cfg(target_os = "linux")]
    linux::advise_willneed(fd, offset, len);
    #[cfg(not(target_os = "linux"))]
    let _ = (fd, offset, len);
}

/// Owning descriptor used to issue WILLNEED hints from the prefetch producer
/// thread. On Linux it is a duplicated [`OwnedFd`](std::os::fd::OwnedFd) kept
/// alive for the duration of the producer loop: `posix_fadvise` requires a
/// borrowable descriptor (`AsFd`), and duplicating one is the only way to obtain
/// that safely — without a raw, `unsafe` `BorrowedFd::borrow_raw` — once the
/// original `File` has been type-erased. On non-Linux platforms hinting is
/// unsupported, so the handle is the uninhabited
/// [`Infallible`](std::convert::Infallible) and is always `None`. Using
/// `Infallible` keeps the cross-platform [`PrefetchReader`] signatures free of
/// any `std::os::fd` type (which does not exist on, e.g., Windows).
///
/// [`PrefetchReader`]: crate::prefetch_reader::PrefetchReader
#[cfg(target_os = "linux")]
pub type HintFd = std::os::fd::OwnedFd;

/// Non-Linux variant of the WILLNEED hint handle — see the Linux [`HintFd`] for
/// the rationale.
#[cfg(not(target_os = "linux"))]
pub type HintFd = std::convert::Infallible;

/// Duplicate `file`'s descriptor into an owning [`HintFd`] so it can be passed
/// to [`advise_willneed`] from a producer thread that has moved (and
/// type-erased) the original `File`. The duplicate shares the underlying open
/// file description and is used only for advisory hints, never for reads, so the
/// shared file offset is irrelevant. Returns `None` if duplication fails (e.g.
/// the process is out of descriptors) — hinting then degrades to a no-op.
#[must_use]
#[cfg(target_os = "linux")]
pub fn hint_fd(file: &File) -> Option<HintFd> {
    use std::os::fd::AsFd;
    file.as_fd().try_clone_to_owned().ok()
}

/// Non-Linux stub — always returns `None` because `posix_fadvise` is
/// unavailable.
#[must_use]
#[cfg(not(target_os = "linux"))]
pub fn hint_fd(_file: &File) -> Option<HintFd> {
    None
}

#[cfg(target_os = "linux")]
mod linux {
    use std::fs::File;
    use std::os::fd::{AsRawFd, OwnedFd};

    use nix::fcntl::{PosixFadviseAdvice, posix_fadvise};

    pub(super) fn advise_sequential(file: &File) {
        // offset = 0, len = 0 applies the hint to the entire file (current and
        // future contents). `&File` implements `AsFd`, so pass it straight to
        // `posix_fadvise`. Errors are non-fatal — log and move on.
        if let Err(e) = posix_fadvise(file, 0, 0, PosixFadviseAdvice::POSIX_FADV_SEQUENTIAL) {
            log::debug!(
                "posix_fadvise(POSIX_FADV_SEQUENTIAL) failed on fd {}: {e}",
                file.as_raw_fd()
            );
        }
    }

    pub(super) fn advise_willneed(fd: &OwnedFd, offset: i64, len: i64) {
        // `&OwnedFd` implements `AsFd`; `posix_fadvise` borrows it for the call.
        if let Err(e) = posix_fadvise(fd, offset, len, PosixFadviseAdvice::POSIX_FADV_WILLNEED) {
            log::debug!("posix_fadvise(POSIX_FADV_WILLNEED) failed on fd {}: {e}", fd.as_raw_fd());
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
