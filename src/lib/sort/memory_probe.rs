//! Per-checkpoint memory instrumentation for the sort pipeline.
//!
//! The probe answers a single question: how much of the process RSS is
//! accounted for by the sort buffer (`RecordBuffer::memory_usage`) versus
//! everything else (pool queues, reorder buffers, mimalloc arenas, noodles
//! internals, etc.)?
//!
//! At each checkpoint it samples process RSS, the caller-supplied tracked
//! byte count, and logs the difference as `residual`. A monotonically
//! increasing residual across spill cycles is the fingerprint for untracked
//! allocation growth (the dominant suspect for issue #238: peak memory usage
//! versus samtools).
//!
//! # Overhead
//!
//! All sampling and formatting is gated on
//! `log::log_enabled!(target: "fgumi_lib::sort::memory_probe", Info)`, so when
//! probe logging is off the hooks are a single atomic load. Enable with
//! `RUST_LOG=fgumi_lib::sort::memory_probe=info` (or at the crate level).
//!
//! # Usage
//!
//! ```ignore
//! use crate::sort::memory_probe::{BufferProbeStats, SpillProbe, log_snapshot};
//!
//! let mut probe = SpillProbe::new("phase1");
//! // ... inside spill loop:
//! probe.pre_spill(buf_stats, Some(pool.phase1_queue_depths()));
//! // ... after drain_pending_spill:
//! probe.post_drain(buf_stats, Some(pool.phase1_queue_depths()));
//! // ... after buffer.clear():
//! probe.post_spill(Some(pool.phase1_queue_depths()));
//! // ... at end of phase 1:
//! probe.phase1_end(buffer.memory_usage() as u64);
//!
//! // ... inside merge:
//! log_snapshot("phase2.start", 0);
//! log_snapshot("phase2.end", 0);
//! ```

use log::{Level, info, log_enabled};

/// Force mimalloc to release retained arena segments back to the OS.
///
/// `mi_collect(true)` is a "force" collect that triggers cross-thread free
/// list draining and arena purging. Used by the probe to distinguish
/// "mimalloc is hoarding pages" from "real live allocations" — if
/// `phys_footprint` drops sharply after a collect, the gap was retained
/// arena memory.
#[allow(unsafe_code)] // mimalloc FFI
pub fn force_mi_collect() {
    // SAFETY: mi_collect(true) is a thread-safe arena maintenance call.
    unsafe {
        libmimalloc_sys::mi_collect(true);
    }
}

/// Log target for the memory probe. Enable with
/// `RUST_LOG=fgumi_lib::sort::memory_probe=info`.
const TARGET: &str = "fgumi_lib::sort::memory_probe";

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
#[allow(unsafe_code)] // mach task_info FFI; required for phys_footprint metric
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

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
#[must_use]
pub fn process_rss_bytes() -> Option<u64> {
    None
}

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

/// Whether the probe is currently enabled. Callers can short-circuit their
/// own work with this when the tracked-bytes computation is non-trivial.
#[must_use]
#[inline]
pub fn enabled() -> bool {
    log_enabled!(target: TARGET, Level::Info)
}

/// Format a byte count as a human-readable value with two-decimal precision.
#[allow(clippy::cast_precision_loss)]
fn format_bytes(bytes: u64) -> String {
    const KIB: f64 = 1024.0;
    const MIB: f64 = 1024.0 * 1024.0;
    const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
    let b = bytes as f64;
    if b >= GIB {
        format!("{:.2}G", b / GIB)
    } else if b >= MIB {
        format!("{:.2}M", b / MIB)
    } else if b >= KIB {
        format!("{:.2}K", b / KIB)
    } else {
        format!("{bytes}B")
    }
}

/// Sample process RSS and log a single memory snapshot under `label`.
///
/// The log line has the form
///
/// ```text
/// MEM[label] rss=8.20G tracked=0.70G residual=7.50G residual_pct=91%
/// ```
///
/// `tracked` is the caller's estimate of accounted-for bytes (typically
/// `RecordBuffer::memory_usage()`). `residual = rss - tracked` is the signal
/// we care about — when it grows monotonically across spill cycles, untracked
/// allocations are piling up.
///
/// When probe logging is off this is a single atomic load via `log_enabled!`;
/// no RSS sampling, no formatting, no allocation.
pub fn log_snapshot(label: &str, tracked_bytes: u64) {
    if !enabled() {
        return;
    }
    let rss = process_rss_bytes();
    let rss_str = rss.map_or_else(|| "?".to_string(), format_bytes);
    let tracked_str = format_bytes(tracked_bytes);
    // Force mimalloc arena collection so we can compare retained vs live.
    force_mi_collect();
    let rss_after = process_rss_bytes();
    let rss_after_str = rss_after.map_or_else(|| "?".to_string(), format_bytes);
    let collected = match (rss, rss_after) {
        (Some(a), Some(b)) => format_bytes(a.saturating_sub(b)),
        _ => "?".to_string(),
    };
    match rss {
        Some(r) => {
            let residual = r.saturating_sub(tracked_bytes);
            let residual_str = format_bytes(residual);
            let pct = if r == 0 { 0 } else { (residual * 100) / r };
            info!(
                target: TARGET,
                "MEM[{label}] rss={rss_str} post_collect={rss_after_str} collected={collected} tracked={tracked_str} residual={residual_str} residual_pct={pct}%"
            );
        }
        None => {
            info!(
                target: TARGET,
                "MEM[{label}] rss=? tracked={tracked_str} residual=? residual_pct=?"
            );
        }
    }
}

/// Per-spill buffer stats for Phase 1 probes.
#[derive(Copy, Clone, Debug)]
pub struct BufferProbeStats {
    /// Logical bytes stored (data + refs).
    pub usage: u64,
    /// Total allocated capacity (segments + refs Vec).
    pub capacity: u64,
    /// Number of records in the buffer.
    pub records: u64,
    /// Number of data segments.
    pub segments: u64,
}

/// Log a Phase 1 snapshot with buffer stats and optional pool queue depths.
///
/// Format: `MEM[label] rss=... buf_use=... buf_cap=... recs=... segs=... [raw_q=... decomp_q=... buf_q=...]`
fn log_phase1_snapshot(
    label: &str,
    buf_stats: Option<BufferProbeStats>,
    pool_depths: Option<(usize, usize, usize)>,
) {
    if !enabled() {
        return;
    }
    let rss = process_rss_bytes();
    let rss_str = rss.map_or_else(|| "?".to_string(), format_bytes);
    force_mi_collect();
    let rss_after = process_rss_bytes();
    let rss_after_str = rss_after.map_or_else(|| "?".to_string(), format_bytes);
    let collected = match (rss, rss_after) {
        (Some(a), Some(b)) => format_bytes(a.saturating_sub(b)),
        _ => "?".to_string(),
    };
    let buf_str = match buf_stats {
        Some(s) => format!(
            " buf_use={} buf_cap={} recs={} segs={}",
            format_bytes(s.usage),
            format_bytes(s.capacity),
            s.records,
            s.segments,
        ),
        None => String::new(),
    };
    let pool_str = match pool_depths {
        Some((raw_q, decomp_q, buf_q)) => {
            format!(" raw_q={raw_q} decomp_q={decomp_q} buf_q={buf_q}")
        }
        None => String::new(),
    };
    let tracked = buf_stats.map_or(0, |s| s.usage);
    let residual_str = match rss {
        Some(r) => {
            let residual = r.saturating_sub(tracked);
            let pct = if r == 0 { 0 } else { (residual * 100) / r };
            format!(" residual={} residual_pct={pct}%", format_bytes(residual))
        }
        None => " residual=? residual_pct=?".to_string(),
    };
    info!(
        target: TARGET,
        "MEM[{label}] rss={rss_str} post_collect={rss_after_str} collected={collected}{buf_str}{pool_str}{residual_str}"
    );
}

/// Spill-cycle probe helper.
///
/// Encapsulates the spill counter so call sites don't format strings at every
/// iteration. Construct at the top of a sort phase, call `pre_spill`/`post_spill`
/// across the spill trigger, and `phase1_end` once before the merge step.
///
/// Each boundary logs a snapshot with a label like `phase1.pre_spill_0`,
/// `phase1.post_spill_0`, `phase1.pre_spill_1`, …
///
/// Also supports periodic mid-read sampling (every `READ_SAMPLE_INTERVAL`
/// records) to track RSS growth between spills.
pub struct SpillProbe {
    phase: &'static str,
    spill_idx: usize,
    read_sample_idx: usize,
    next_read_threshold: u64,
}

impl SpillProbe {
    /// Sample RSS every this-many records during the Phase 1 read loop.
    pub const READ_SAMPLE_INTERVAL: u64 = 1_000_000;

    /// Create a probe for the named phase and log a `"<phase>.start"`
    /// snapshot with tracked=0.
    #[must_use]
    pub fn new(phase: &'static str) -> Self {
        log_snapshot(&format!("{phase}.start"), 0);
        Self {
            phase,
            spill_idx: 0,
            read_sample_idx: 0,
            next_read_threshold: Self::READ_SAMPLE_INTERVAL,
        }
    }

    /// Return true if a mid-read sample should be logged at the current record count.
    #[inline]
    #[must_use]
    pub fn should_sample_read(&self, records_read: u64) -> bool {
        enabled() && records_read >= self.next_read_threshold
    }

    /// Log a mid-read sample with buffer stats and optional pool queue depths.
    ///
    /// Only call this when `should_sample_read()` returned true.
    pub fn log_mid_read(
        &mut self,
        buf_stats: BufferProbeStats,
        pool_depths: Option<(usize, usize, usize)>,
    ) {
        log_phase1_snapshot(
            &format!("{}.mid_read_{}", self.phase, self.read_sample_idx),
            Some(buf_stats),
            pool_depths,
        );
        self.read_sample_idx += 1;
        self.next_read_threshold =
            self.next_read_threshold.saturating_add(Self::READ_SAMPLE_INTERVAL);
    }

    /// Log a pre-spill snapshot with buffer stats and optional pool depths.
    pub fn pre_spill(
        &self,
        buf_stats: BufferProbeStats,
        pool_depths: Option<(usize, usize, usize)>,
    ) {
        log_phase1_snapshot(
            &format!("{}.pre_spill_{}", self.phase, self.spill_idx),
            Some(buf_stats),
            pool_depths,
        );
    }

    /// Log a post-drain snapshot (after previous spill is drained, buffer still full).
    pub fn post_drain(
        &self,
        buf_stats: BufferProbeStats,
        pool_depths: Option<(usize, usize, usize)>,
    ) {
        log_phase1_snapshot(
            &format!("{}.post_drain_{}", self.phase, self.spill_idx),
            Some(buf_stats),
            pool_depths,
        );
    }

    /// Log a post-spill snapshot (tracked=0) and advance the spill index.
    pub fn post_spill(&mut self, pool_depths: Option<(usize, usize, usize)>) {
        log_phase1_snapshot(
            &format!("{}.post_spill_{}", self.phase, self.spill_idx),
            None,
            pool_depths,
        );
        self.spill_idx += 1;
        // Reset mid-read counter for next fill cycle.
        self.read_sample_idx = 0;
        self.next_read_threshold =
            self.next_read_threshold.saturating_add(Self::READ_SAMPLE_INTERVAL);
    }

    /// Log a snapshot marking the end of phase 1 (before merge).
    pub fn phase1_end(&self, tracked: u64) {
        log_snapshot(&format!("{}.end", self.phase), tracked);
    }

    /// Number of spill cycles observed so far.
    #[must_use]
    pub fn spill_count(&self) -> usize {
        self.spill_idx
    }
}

/// Merge-phase probe helper.
///
/// Encapsulates a periodic sampler that logs a snapshot every
/// `SAMPLE_INTERVAL_RECORDS` merged records during phase 2. Like `SpillProbe`,
/// the hot path is gated on `log_enabled!` so probes off = single atomic load.
///
/// Each sample is labelled `phase2.mid_N` where N is the sample index (0..).
/// Optional consumer-side byte accounting sampled at mid-merge boundaries.
#[derive(Copy, Clone, Debug)]
pub struct ConsumerProbeStats {
    pub current_bytes: u64,
    pub current_capacity: u64,
    pub pending_blocks: u64,
    pub pending_bytes: u64,
    pub active_sources: u64,
}

pub struct MergeProbe {
    sample_idx: usize,
    next_threshold: u64,
}

impl MergeProbe {
    /// Sample the RSS every this-many merged records when probe logging is on.
    /// Chosen to give ~20 samples on a 200M-record merge — fine-grained enough
    /// to see the growth-curve shape, sparse enough that the sysinfo/proc
    /// sampling overhead is negligible.
    pub const SAMPLE_INTERVAL_RECORDS: u64 = 1_000_000;

    /// Create a new merge probe. Logs a `"phase2.start"` snapshot with
    /// tracked=0.
    #[must_use]
    pub fn new() -> Self {
        log_snapshot("phase2.start", 0);
        Self { sample_idx: 0, next_threshold: Self::SAMPLE_INTERVAL_RECORDS }
    }

    /// Called once per merged record with the running total. Samples RSS
    /// when `records_merged` crosses the next threshold.
    #[inline]
    pub fn record(&mut self, records_merged: u64) {
        if records_merged < self.next_threshold {
            return;
        }
        if enabled() {
            log_snapshot(&format!("phase2.mid_{}", self.sample_idx), 0);
        }
        self.sample_idx += 1;
        self.next_threshold = self.next_threshold.saturating_add(Self::SAMPLE_INTERVAL_RECORDS);
    }

    /// Number of mid-merge samples logged so far.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.sample_idx
    }

    /// Return true if a mid-merge sample should be logged at the current
    /// record count.
    ///
    /// Call sites use this to gate expensive context gathering (pool queue
    /// depths, buffer pool counts) so those observations only happen on the
    /// ~20-sample boundaries, not every record.
    #[inline]
    #[must_use]
    pub fn should_sample(&self, records_merged: u64) -> bool {
        records_merged >= self.next_threshold
    }

    /// Log a mid-merge sample with per-component queue depths.
    ///
    /// Only call this when `should_sample()` returned true for the same
    /// `records_merged` count. Advances the sample counter and the next
    /// threshold.
    ///
    /// The log line extends the standard `MEM[...]` format with additional
    /// `raw_q=` `decomp_q=` `buf_q=` fields so a single snapshot captures
    /// both the process RSS and the pool plumbing state at the same instant.
    pub fn log_mid_with_depths(
        &mut self,
        raw_chunk_q: usize,
        decompressed_chunk_q: usize,
        buffer_pool_available: usize,
        consumer_stats: Option<ConsumerProbeStats>,
    ) {
        if enabled() {
            let rss = process_rss_bytes();
            let rss_str = rss.map_or_else(|| "?".to_string(), format_bytes);
            // Force mimalloc arena collection and re-sample. The delta
            // (rss - rss_after_collect) is the retained-arena footprint.
            force_mi_collect();
            let rss_after = process_rss_bytes();
            let rss_after_str = rss_after.map_or_else(|| "?".to_string(), format_bytes);
            let collected = match (rss, rss_after) {
                (Some(a), Some(b)) => format_bytes(a.saturating_sub(b)),
                _ => "?".to_string(),
            };
            let consumer_str = match consumer_stats {
                Some(s) => format!(
                    " cur_bytes={} cur_cap={} pend_blocks={} pend_bytes={} active_src={}",
                    format_bytes(s.current_bytes),
                    format_bytes(s.current_capacity),
                    s.pending_blocks,
                    format_bytes(s.pending_bytes),
                    s.active_sources,
                ),
                None => String::new(),
            };
            info!(
                target: TARGET,
                "MEM[phase2.mid_{}] rss={rss_str} post_collect={rss_after_str} collected={collected} raw_q={raw_chunk_q} decomp_q={decompressed_chunk_q} buf_q={buffer_pool_available}{consumer_str}",
                self.sample_idx,
            );
        }
        self.sample_idx += 1;
        self.next_threshold = self.next_threshold.saturating_add(Self::SAMPLE_INTERVAL_RECORDS);
    }
}

impl Default for MergeProbe {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_bytes_units() {
        assert_eq!(format_bytes(0), "0B");
        assert_eq!(format_bytes(512), "512B");
        assert_eq!(format_bytes(1024), "1.00K");
        assert_eq!(format_bytes(1536), "1.50K");
        assert_eq!(format_bytes(1024 * 1024), "1.00M");
        assert_eq!(format_bytes(1024 * 1024 * 1024), "1.00G");
        assert_eq!(format_bytes(2u64 * 1024 * 1024 * 1024 + 512 * 1024 * 1024), "2.50G");
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

    #[test]
    fn test_log_snapshot_does_not_panic_without_logger() {
        // With no logger configured, log_enabled! returns false and the
        // entire function should short-circuit cleanly.
        log_snapshot("test.label", 1024);
        log_snapshot("test.label", 0);
    }

    #[test]
    fn test_spill_probe_increments() {
        let mut probe = SpillProbe::new("test_phase");
        assert_eq!(probe.spill_count(), 0);
        let stats = BufferProbeStats { usage: 1024, capacity: 2048, records: 10, segments: 1 };
        probe.pre_spill(stats, None);
        probe.post_spill(None);
        assert_eq!(probe.spill_count(), 1);
        let stats2 = BufferProbeStats { usage: 2048, capacity: 4096, records: 20, segments: 1 };
        probe.pre_spill(stats2, None);
        probe.post_spill(None);
        assert_eq!(probe.spill_count(), 2);
        probe.phase1_end(0);
    }

    #[test]
    fn test_merge_probe_samples_at_interval() {
        let mut probe = MergeProbe::new();
        assert_eq!(probe.sample_count(), 0);

        // Below threshold: no samples.
        probe.record(MergeProbe::SAMPLE_INTERVAL_RECORDS - 1);
        assert_eq!(probe.sample_count(), 0);

        // At threshold: first sample.
        probe.record(MergeProbe::SAMPLE_INTERVAL_RECORDS);
        assert_eq!(probe.sample_count(), 1);

        // Same interval doesn't re-fire.
        probe.record(MergeProbe::SAMPLE_INTERVAL_RECORDS + 1);
        assert_eq!(probe.sample_count(), 1);

        // Second threshold fires.
        probe.record(MergeProbe::SAMPLE_INTERVAL_RECORDS * 2);
        assert_eq!(probe.sample_count(), 2);
    }
}
