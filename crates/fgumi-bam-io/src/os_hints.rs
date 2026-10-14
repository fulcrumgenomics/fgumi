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
/// `fd` will be needed soon, in **one** advisory call (see
/// [`advise_willneed_range`] for what the kernel reads and how a range is
/// clipped). On Linux this is `posix_fadvise(POSIX_FADV_WILLNEED)`, which reads
/// at most the device's read-ahead window (`read_ahead_kb`, 128 KiB by
/// default) from `offset` per call however large `len` is; on macOS it is
/// `fcntl(F_RDADVISE)`, which takes the whole range. Elsewhere, and for a
/// negative `offset` or a `len` that is not positive, this is a no-op.
///
/// This is the owning-handle variant for use inside the [`PrefetchReader`]
/// producer thread, which owns the source `File` but type-erases it to a
/// generic reader and so cannot lend a `&File` back. The producer instead holds
/// a duplicated descriptor ([`HintFd`]) alive for the lifetime of its loop and
/// passes it here by reference.
///
/// [`PrefetchReader`]: crate::prefetch_reader::PrefetchReader
pub fn advise_willneed(fd: &HintFd, offset: i64, len: i64) {
    #[cfg(any(target_os = "linux", target_os = "macos"))]
    {
        let (Ok(offset), Ok(len)) = (u64::try_from(offset), u64::try_from(len)) else {
            return;
        };
        if let Err(e) = willneed_calls(offset, len, u64::MAX, platform::MAX_CALL_BYTES)
            .take(1)
            .try_for_each(|(off, n)| platform::willneed(fd, off, n))
        {
            use std::os::fd::AsRawFd;
            log::debug!("read-ahead advice failed on fd {}: {e}", fd.as_raw_fd());
        }
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    let _ = (fd, offset, len);
}

/// Hint to the kernel that `len` bytes of `file` from `offset` will be read
/// soon, so it starts reading them now. A later positional read of the range
/// copies whatever has arrived from the page cache, waits only for the part
/// still in flight, and reads from the device only what the kernel declined to
/// fetch or has since evicted (the advice is never binding).
///
/// On Linux this is `posix_fadvise(POSIX_FADV_WILLNEED)`, issued in 128 KiB
/// calls over the range (clipped to the file's
/// length): one call reads at most the device's read-ahead window from its own
/// offset (`force_page_cache_ra` clips it to
/// `max(io_pages, ra_pages)`), so a single call over a multi-MiB range would
/// start only its head. On macOS this is one `fcntl(F_RDADVISE)` call over the
/// range (at most `c_int::MAX` bytes), which has no such window. Both queue the
/// reads and return without waiting for them. Elsewhere this is a no-op.
///
/// It allocates nothing: the pages land in the page cache, which the kernel
/// reclaims under pressure. An empty range, or one that is empty once clipped
/// to the file (Linux) or to `off_t`, issues no call; a range straddling the
/// end of the file is clipped to it. macOS refuses a range starting at or past
/// the end of the file (`EFBIG`), which, like any other failure, is logged at
/// `debug` level and otherwise ignored.
pub fn advise_willneed_range(file: &File, offset: u64, len: u64) {
    #[cfg(any(target_os = "linux", target_os = "macos"))]
    if let Err(e) = willneed_range(file, offset, len) {
        use std::os::fd::AsRawFd;
        log::debug!("read-ahead advice failed on fd {}: {e}", file.as_raw_fd());
    }
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    let _ = (file, offset, len);
}

/// The largest range one Linux `POSIX_FADV_WILLNEED` call over a fill covers:
/// 128 KiB, the default `read_ahead_kb` and so the smallest read-ahead window
/// one call is clipped to on a block device in its default configuration.
/// Larger calls would leave part of each chunk unread on such a device; smaller
/// ones only add calls (a 4 MiB fill is 32 of them).
#[cfg(any(target_os = "linux", test))]
const WILLNEED_CHUNK_BYTES: u64 = 128 << 10;

/// The `(offset, len)` advisory calls covering `[offset, offset + len)`, the
/// range first clipped to `end_limit` (the file's length, or `u64::MAX` to
/// leave it to the kernel) and to `off_t`, then split into pieces of at most
/// `max_call` bytes. Every piece is non-empty, so an empty range, or one that
/// is empty once clipped, yields no call: a zero-length `POSIX_FADV_WILLNEED`
/// would otherwise mean "to the end of the file".
#[cfg(any(target_os = "linux", target_os = "macos", test))]
fn willneed_calls(
    offset: u64,
    len: u64,
    end_limit: u64,
    max_call: u64,
) -> impl Iterator<Item = (i64, i64)> {
    const OFF_T_MAX: u64 = i64::MAX.cast_unsigned();
    let end = offset.saturating_add(len).min(end_limit).min(OFF_T_MAX);
    let max_call = max_call.clamp(1, OFF_T_MAX);
    let mut at = offset;
    std::iter::from_fn(move || {
        if at >= end {
            return None;
        }
        let n = (end - at).min(max_call);
        // Both fit `i64`: `at < end <= OFF_T_MAX` and `n <= end - at`.
        let call = (at.cast_signed(), n.cast_signed());
        at += n;
        Some(call)
    })
}

/// [`advise_willneed_range`]'s calls, returning the first failure.
#[cfg(target_os = "linux")]
fn willneed_range(file: &File, offset: u64, len: u64) -> std::io::Result<()> {
    if len == 0 {
        return Ok(());
    }
    // Bound the loop to the file: the kernel would clip each call anyway, but
    // an oversized `len` must not turn into billions of calls.
    let file_len = file.metadata()?.len();
    willneed_calls(offset, len, file_len, WILLNEED_CHUNK_BYTES)
        .try_for_each(|(off, n)| platform::willneed(file, off, n))
}

/// [`advise_willneed_range`]'s one call, returning its failure.
#[cfg(target_os = "macos")]
fn willneed_range(file: &File, offset: u64, len: u64) -> std::io::Result<()> {
    willneed_calls(offset, len, u64::MAX, platform::MAX_CALL_BYTES)
        .take(1)
        .try_for_each(|(off, n)| platform::willneed(file, off, n))
}

#[cfg(target_os = "linux")]
use linux as platform;
#[cfg(target_os = "macos")]
use macos as platform;

/// Owning descriptor used to issue WILLNEED hints from the prefetch producer
/// thread. On Linux and macOS it is a duplicated
/// [`OwnedFd`](std::os::fd::OwnedFd) kept alive for the duration of the
/// producer loop: the advisory calls require a borrowable descriptor (`AsFd`),
/// and duplicating one is the only way to obtain that safely — without a raw,
/// `unsafe` `BorrowedFd::borrow_raw` — once the original `File` has been
/// type-erased. On other platforms hinting is unsupported, so the handle is the
/// uninhabited [`Infallible`](std::convert::Infallible) and is always `None`.
/// Using `Infallible` keeps the cross-platform [`PrefetchReader`] signatures
/// free of any `std::os::fd` type (which does not exist on, e.g., Windows).
///
/// [`PrefetchReader`]: crate::prefetch_reader::PrefetchReader
#[cfg(any(target_os = "linux", target_os = "macos"))]
pub type HintFd = std::os::fd::OwnedFd;

/// Variant of the WILLNEED hint handle for platforms without an advisory call —
/// see the Linux/macOS [`HintFd`] for the rationale.
#[cfg(not(any(target_os = "linux", target_os = "macos")))]
pub type HintFd = std::convert::Infallible;

/// Duplicate `file`'s descriptor into an owning [`HintFd`] so it can be passed
/// to [`advise_willneed`] from a producer thread that has moved (and
/// type-erased) the original `File`. The duplicate shares the underlying open
/// file description and is used only for advisory hints, never for reads, so the
/// shared file offset is irrelevant. Returns `None` if duplication fails (e.g.
/// the process is out of descriptors) — hinting then degrades to a no-op.
#[must_use]
#[cfg(any(target_os = "linux", target_os = "macos"))]
pub fn hint_fd(file: &File) -> Option<HintFd> {
    use std::os::fd::AsFd;
    file.as_fd().try_clone_to_owned().ok()
}

/// Stub for platforms without an advisory call — always returns `None`.
#[must_use]
#[cfg(not(any(target_os = "linux", target_os = "macos")))]
pub fn hint_fd(_file: &File) -> Option<HintFd> {
    None
}

#[cfg(target_os = "linux")]
mod linux {
    use std::fs::File;
    use std::os::fd::{AsFd, AsRawFd};

    use nix::fcntl::{PosixFadviseAdvice, posix_fadvise};

    /// One call takes any `off_t` length (the kernel clips it to its window).
    pub(super) const MAX_CALL_BYTES: u64 = u64::MAX;

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

    /// One `posix_fadvise(POSIX_FADV_WILLNEED)` call; `len` is positive (a
    /// zero length would mean "to the end of the file").
    pub(super) fn willneed(fd: impl AsFd, offset: i64, len: i64) -> std::io::Result<()> {
        debug_assert!(len > 0, "a zero-length WILLNEED advises to the end of the file");
        Ok(posix_fadvise(fd, offset, len, PosixFadviseAdvice::POSIX_FADV_WILLNEED)?)
    }
}

#[cfg(target_os = "macos")]
mod macos {
    use std::os::fd::AsFd;

    use nix::fcntl::{FcntlArg, fcntl};
    use nix::libc;

    /// `ra_count` is a C `int`: at most `c_int::MAX` (2^31 - 1) bytes a call.
    pub(super) const MAX_CALL_BYTES: u64 = (1 << 31) - 1;

    /// One `fcntl(F_RDADVISE)` call; `len` is positive and at most
    /// [`MAX_CALL_BYTES`].
    pub(super) fn willneed(fd: impl AsFd, offset: i64, len: i64) -> std::io::Result<()> {
        let ra_count = libc::c_int::try_from(len).unwrap_or(libc::c_int::MAX);
        fcntl(fd, FcntlArg::F_RDADVISE(libc::radvisory { ra_offset: offset, ra_count }))?;
        Ok(())
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

    /// The advisory calls tile the clipped range in `max_call` pieces, each
    /// non-empty: an empty range, a range past `end_limit`, and a range
    /// starting at `off_t`'s maximum issue nothing (a zero-length WILLNEED
    /// would advise to the end of the file).
    #[rstest::rstest]
    #[case::one_chunk(0, 100, u64::MAX, 128, &[(0, 100)])]
    #[case::exact_chunks(0, 256, u64::MAX, 128, &[(0, 128), (128, 128)])]
    #[case::ragged_tail(10, 300, u64::MAX, 128, &[(10, 128), (138, 128), (266, 44)])]
    #[case::clipped_to_the_file(100, 1 << 20, 300, 128, &[(100, 128), (228, 72)])]
    #[case::empty(5, 0, u64::MAX, 128, &[])]
    #[case::at_the_end_of_the_file(300, 1, 300, 128, &[])]
    #[case::past_the_end_of_the_file(400, 1, 300, 128, &[])]
    #[case::length_overflows_u64(1, u64::MAX, u64::MAX, u64::MAX, &[(1, i64::MAX - 1)])]
    #[case::offset_past_off_t(u64::MAX, u64::MAX, u64::MAX, 128, &[])]
    #[case::zero_after_clipping_to_off_t(i64::MAX.cast_unsigned(), 1, u64::MAX, 128, &[])]
    #[case::last_off_t_byte((i64::MAX - 1).cast_unsigned(), 9, u64::MAX, 128, &[(i64::MAX - 1, 1)])]
    fn willneed_calls_tile_the_clipped_range(
        #[case] offset: u64,
        #[case] len: u64,
        #[case] end_limit: u64,
        #[case] max_call: u64,
        #[case] want: &[(i64, i64)],
    ) {
        let got: Vec<_> = willneed_calls(offset, len, end_limit, max_call).collect();
        assert_eq!(got, want);
    }

    /// A 4 MiB fill is advised on Linux as 32 calls of the chunk size, which
    /// together cover it exactly.
    #[test]
    fn a_hot_fill_is_advised_in_full_chunk_calls() {
        let fill = 4u64 << 20;
        let calls: Vec<_> = willneed_calls(1 << 30, fill, u64::MAX, WILLNEED_CHUNK_BYTES).collect();
        assert_eq!(calls.len(), 32);
        assert!(calls.iter().all(|&(_, n)| n.cast_unsigned() == WILLNEED_CHUNK_BYTES));
        assert_eq!(calls.iter().map(|&(_, n)| n.cast_unsigned()).sum::<u64>(), fill);
    }

    /// The advice never disturbs the file, whatever the range, and the platform
    /// call really accepts the ranges the clipping produces (its errors are not
    /// swallowed here): in the file, straddling its end, and too large for one
    /// call. Ranges at or past the end of the file go only through the
    /// best-effort wrapper, because macOS refuses them (`EFBIG`).
    #[test]
    fn advise_willneed_range_is_safe_for_any_range() {
        let mut tmp = NamedTempFile::new().expect("create temp file");
        let data: Vec<u8> = (0..=255u8).cycle().take(1 << 20).collect();
        tmp.write_all(&data).expect("write temp file");
        let f = File::open(tmp.path()).expect("reopen temp file for read");
        let len = data.len() as u64;
        for (offset, n) in [
            (0, len),
            (4096, 8192),
            (len - 1, 1 << 20),
            (len, 1),
            (len + (1 << 30), 1 << 20),
            (0, 0),
            (0, u64::MAX),
            (u64::MAX, u64::MAX),
            (i64::MAX.cast_unsigned(), 1),
        ] {
            advise_willneed_range(&f, offset, n);
        }
        #[cfg(any(target_os = "linux", target_os = "macos"))]
        for (offset, n) in [(0, len), (4096, 8192), (len - 1, 1 << 20), (0, u64::MAX), (0, 0)] {
            willneed_range(&f, offset, n)
                .unwrap_or_else(|e| panic!("advising ({offset}, {n}) failed: {e}"));
        }
        #[cfg(any(target_os = "linux", target_os = "macos"))]
        {
            let fd = hint_fd(&f).expect("duplicate the descriptor");
            advise_willneed(&fd, 0, i64::MAX);
            advise_willneed(&fd, -1, 1);
            advise_willneed(&fd, 0, 0);
        }
        let mut back = vec![0u8; data.len()];
        std::os::unix::fs::FileExt::read_exact_at(&f, &mut back, 0).expect("read back");
        assert_eq!(back, data, "advising a range leaves its bytes unchanged");
    }
}
