//! Software prefetch for cold-arena forward scans.

/// Bytes ahead of the current frame to software-prefetch in a forward scan of
/// arena frames that extracts keys: the arena boundary scan in
/// `fgumi-pipeline-io` and [`TemplateKeyContext::extract_batch`](crate::TemplateKeyContext::extract_batch).
/// Chosen from a microbench of the boundary+key scan over a cold ~2.6 GiB arena
/// at the production density of ~220 B/record (~9 records of lead): 2 KiB gave
/// the best speedup (~15%), 1 KiB or less gave a negligible one (too little lead
/// time), and 4 KiB performed the same as 2 KiB. The batched template key scan
/// reads the same cold frames, so it uses the same lead.
pub const KEY_PREFETCH_DISTANCE: usize = 2048;

/// Software-prefetch (read, into L1, temporal) the cache line containing `byte`.
///
/// The sort's forward scans walk arena bytes that another thread wrote long
/// before (the inflate workers fill the arena; the boundary scan and the
/// key-extraction workers read it cold), so they are latency-bound on cache
/// misses. Prefetching a fixed distance ahead hides that latency. One copy for
/// the whole sort engine: `fgumi-pipeline-io`'s boundary scan and
/// `TemplateKeyContext::extract_batch` both call this.
///
/// SAFETY note: `prfm pldl1keep` (`aarch64`) and `_mm_prefetch` (`x86_64`) are
/// non-faulting hints — they never read or write observable memory and never
/// trap, even on an unmapped address. `byte` is a live `&u8`, so the pointer is
/// valid to name. A no-op on other architectures. Two cfg-gated
/// `#[allow(unsafe_code)]` sites, listed in CLAUDE.md
/// §"Approved hot-path unsafe (sort engine)".
#[inline]
pub fn prefetch_read_l1(byte: &u8) {
    let ptr: *const u8 = byte;
    #[cfg(target_arch = "aarch64")]
    #[allow(unsafe_code)]
    // SAFETY: `prfm pldl1keep` is a non-faulting prefetch hint over a valid pointer.
    unsafe {
        core::arch::asm!(
            "prfm pldl1keep, [{p}]",
            p = in(reg) ptr,
            options(nostack, readonly, preserves_flags),
        );
    }
    #[cfg(target_arch = "x86_64")]
    #[allow(unsafe_code)]
    // SAFETY: `_mm_prefetch` is a non-faulting prefetch hint over a valid pointer.
    unsafe {
        core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_T0 }>(ptr.cast());
    }
    #[cfg(not(any(target_arch = "aarch64", target_arch = "x86_64")))]
    {
        let _ = ptr; // no portable stable prefetch; the hint is a no-op elsewhere
    }
}

#[cfg(test)]
mod tests {
    use super::prefetch_read_l1;

    /// A prefetch is a hint: it must accept any live byte — first, last, and
    /// one in the middle of a large buffer — and change nothing observable.
    #[test]
    fn prefetch_accepts_any_live_byte_and_changes_nothing() {
        let buf: Vec<u8> = (0..=255u8).cycle().take(1 << 16).collect();
        let before = buf.clone();
        prefetch_read_l1(&buf[0]);
        prefetch_read_l1(&buf[buf.len() / 2]);
        prefetch_read_l1(&buf[buf.len() - 1]);
        assert_eq!(buf, before);
    }
}
