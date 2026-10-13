//! Allocator and RSS helpers: the process resident-set probe the pipeline's
//! thread-telemetry sampler reads, and the mimalloc controls (retaining freed
//! memory, forcing a collection, printing statistics) the sort and align stages
//! use.

#[cfg(feature = "memory-debug")]
#[allow(unused_imports)]
// consumed by main fgumi via the crate-root re-export
pub use platform_ffi::print_mi_stats;
pub use platform_ffi::{
    force_mi_collect, mi_purge_delay_ms, process_rss_bytes, retain_freed_memory,
};

/// mimalloc purge delay, in milliseconds, that [`retain_freed_memory`] sets.
///
/// mimalloc's default (1000 ms) decommits a freed page after a second, so the
/// align stage's per-batch buffers, which are freed and reallocated well within
/// that, fault back in on every batch. Never purging (`-1`) avoids those faults
/// but also never returns any freed page to the OS: in a fused `runall`, each
/// sort spill's freed memory then stays resident, and peak RSS grows by roughly
/// the sort's memory budget. Sixty seconds keeps the batch-to-batch reuse
/// resident (page faults close to never purging, against 25x or more at the
/// default) while memory that sits idle across a spill cycle is still returned.
pub const RETAINED_PURGE_DELAY_MS: i32 = 60_000;

/// Platform-specific FFI helpers isolated from the main module to contain the
/// `#[allow(unsafe_code)]` surface to the smallest possible scope.
#[allow(unsafe_code)]
mod platform_ffi {
    /// Force mimalloc to release retained arena segments back to the OS.
    ///
    /// `mi_collect(true)` is a "force" collect that triggers cross-thread free
    /// list draining and arena purging. Used by the probe to distinguish
    /// "mimalloc is hoarding pages" from "real live allocations" — if
    /// `phys_footprint` drops sharply after a collect, the gap was retained
    /// arena memory.
    pub fn force_mi_collect() {
        // SAFETY: mi_collect(true) is a thread-safe arena maintenance call.
        unsafe {
            libmimalloc_sys::mi_collect(true);
        }
    }

    /// mimalloc's `mi_option_purge_delay`. `libmimalloc-sys` does not export it
    /// as a constant, so this is its position in `mi_option_t`, which is 15 in
    /// both the v2 and v3 `mimalloc.h` enums (deprecated slots keep their
    /// places). `mi_purge_delay_ms_matches_mimalloc_default` pins it against the
    /// linked mimalloc v3's default of 1000 ms.
    const MI_OPTION_PURGE_DELAY: libmimalloc_sys::mi_option_t = 15;

    /// Delay mimalloc returning freed pages to the OS by
    /// [`RETAINED_PURGE_DELAY_MS`](super::RETAINED_PURGE_DELAY_MS), so a page
    /// that is freed and soon reused is not decommitted and faulted back in.
    /// Trades a somewhat higher peak RSS for fewer page faults and less system
    /// time.
    pub fn retain_freed_memory() {
        // SAFETY: mi_option_set only stores the option value in mimalloc's
        // option table; mimalloc reads it on each purge, so setting it after
        // start-up is supported, and the option index is fixed (see above).
        unsafe {
            libmimalloc_sys::mi_option_set(
                MI_OPTION_PURGE_DELAY,
                std::ffi::c_long::from(super::RETAINED_PURGE_DELAY_MS),
            );
        }
    }

    /// mimalloc's current purge delay in milliseconds (`-1`: never purge).
    #[must_use]
    pub fn mi_purge_delay_ms() -> i64 {
        // SAFETY: mi_option_get only reads mimalloc's option table.
        #[allow(clippy::useless_conversion)] // `c_long` is `i32` on Windows
        i64::from(unsafe { libmimalloc_sys::mi_option_get(MI_OPTION_PURGE_DELAY) })
    }

    /// Print mimalloc allocator statistics to stderr via the default output handler.
    ///
    /// `mi_stats_print_out(None, null_mut())` uses mimalloc's internal synchronization,
    /// making it safe to call concurrently with allocation/deallocation on other threads.
    #[cfg(feature = "memory-debug")]
    pub fn print_mi_stats() {
        // SAFETY: mimalloc synchronizes stats collection internally.
        unsafe {
            libmimalloc_sys::mi_stats_print_out(None, std::ptr::null_mut());
        }
    }

    /// Read the current process resident-set size in bytes.
    ///
    /// Returns `None` if the platform sampler is unavailable (e.g. `/proc` is not
    /// mounted on a non-Linux host and the `sysinfo` fallback fails).
    ///
    /// # Implementation
    ///
    /// - **Linux**: reads `VmRSS` from `/proc/self/status`. One short file read,
    ///   no allocations on the hot path beyond the status string itself.
    /// - **macOS**: reads `phys_footprint` from `task_info(TASK_VM_INFO)`. This is
    ///   what Activity Monitor reports as "Memory" and excludes purgeable
    ///   `MADV_FREE` pages, unlike `task_basic_info::resident_size` (which sysinfo
    ///   uses). For sort memory diagnosis we need this metric — Darwin keeps
    ///   freed pages in `resident_size` until kernel pressure forces eviction,
    ///   which produces a misleading ~30 GiB peak when mimalloc has actually
    ///   already purged them.
    #[cfg(target_os = "macos")]
    #[must_use]
    pub fn process_rss_bytes() -> Option<u64> {
        use mach2::message::mach_msg_type_number_t;
        use mach2::task::task_info;
        use mach2::task_info::{TASK_VM_INFO, task_vm_info};
        use mach2::traps::mach_task_self;
        use mach2::vm_types::natural_t;

        let mut info = task_vm_info::default();
        let mut count: mach_msg_type_number_t = mach_msg_type_number_t::try_from(
            std::mem::size_of::<task_vm_info>() / std::mem::size_of::<natural_t>(),
        )
        .ok()?;
        // SAFETY: mach task_info FFI; passes a properly-sized output buffer and count.
        let kr = unsafe {
            task_info(
                mach_task_self(),
                TASK_VM_INFO,
                std::ptr::addr_of_mut!(info).cast(),
                std::ptr::addr_of_mut!(count),
            )
        };
        if kr != 0 {
            return None;
        }
        Some(info.phys_footprint)
    }

    /// Unsupported platform: RSS probing is unavailable here, so this always
    /// returns `None`. (See the `macos` variant above for the metric rationale.)
    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    #[must_use]
    pub fn process_rss_bytes() -> Option<u64> {
        None
    }

    /// Linux: reads `VmRSS` from `/proc/self/status` (kibibytes) and returns it
    /// in bytes. One short file read, no hot-path allocation beyond the status
    /// string. Returns `None` if the field is absent or unparseable.
    #[cfg(target_os = "linux")]
    #[must_use]
    pub fn process_rss_bytes() -> Option<u64> {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        status
            .lines()
            .find(|line| line.starts_with("VmRSS:"))?
            .split_whitespace()
            .nth(1)?
            .parse::<u64>()
            .ok()
            .map(|kb| kb * 1024)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pins [`platform_ffi`]'s purge-delay option index: the linked mimalloc v3
    /// defaults `mi_option_purge_delay` to 1000 ms (v2 used 10 ms), and retaining
    /// freed memory sets it to [`RETAINED_PURGE_DELAY_MS`] (a bounded delay, not
    /// `-1`'s never-purge). Skipped when the environment overrides it.
    #[test]
    fn mi_purge_delay_ms_matches_mimalloc_default() {
        // mimalloc matches its option names case-insensitively.
        if std::env::vars_os().any(|(name, _)| {
            name.eq_ignore_ascii_case("MIMALLOC_PURGE_DELAY")
                || name.eq_ignore_ascii_case("MIMALLOC_RESET_DELAY")
        }) {
            return;
        }
        assert_eq!(mi_purge_delay_ms(), 1000);
        retain_freed_memory();
        assert_eq!(mi_purge_delay_ms(), i64::from(RETAINED_PURGE_DELAY_MS));
    }

    #[test]
    fn test_process_rss_bytes_returns_plausible_value() {
        // Every host we target has either /proc/self/status or a working
        // sysinfo fallback, so this should not return None in CI.
        let rss = process_rss_bytes().expect("RSS sampling should work on supported platforms");
        // 1 MiB is a safe lower bound for any running Rust process.
        assert!(rss >= 1024 * 1024, "RSS {rss} is implausibly small");
        // 1 TiB is a safe upper bound for a test process.
        assert!(rss < 1024u64 * 1024 * 1024 * 1024, "RSS {rss} is implausibly large");
    }
}
