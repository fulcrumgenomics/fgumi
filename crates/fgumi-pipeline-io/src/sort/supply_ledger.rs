//! The merge supply's shared byte ledger.
//!
//! `SpillReadPlanner` adds to it as it issues reads; `SortSpillDecompress`
//! (and the merge, when it serves itself) subtract as slices land and stashed
//! frames are claimed. The planner reads `inflight_slices` and
//! `inflight_cold_bytes` to bound outstanding reads, so the ledger is part of
//! the supply's liveness, not only its diagnostics: the planner and the
//! decompress step must share one, which [`SpillSupply`] guarantees by
//! handing both (and the merge) the same ledger and merge demand. The rest
//! feeds the `--sort-stats` supply lines. All counters are relaxed atomics:
//! each is an independent tally, and the bounds they feed tolerate a
//! momentarily stale value (the planner re-reads on its next pass).

use std::sync::Arc;
use std::sync::atomic::{AtomicI64, AtomicU32, AtomicU64, Ordering};

use fgumi_sort::MergeDemand;

use crate::pread::SpillClass;

/// Shared supply counters (see the module docs).
#[derive(Debug, Default)]
pub struct SupplyLedger {
    inflight_bytes: AtomicU64,
    inflight_cold_bytes: AtomicU64,
    inflight_slices: AtomicU32,
    /// Signed: an ingest publishes its frames to the slot (claimable) before it
    /// adds them here, so a claim's subtraction can precede the addition and
    /// the tally can dip below zero for that interval.
    stash_bytes: AtomicI64,
    stash_peak_bytes: AtomicU64,
    frames_stashed: AtomicU64,
    claims_worker: AtomicU64,
    claims_consumer: AtomicU64,
    read_bytes: AtomicU64,
    read_slices: AtomicU64,
    fills_hot: AtomicU64,
    fills_cold: AtomicU64,
    slices_hot: AtomicU64,
    slices_cold: AtomicU64,
    awaited_starved_reads: AtomicU64,
}

impl SupplyLedger {
    /// The planner issued `slices` read requests of `bytes` in total for
    /// class `c` (one fill).
    pub fn add_inflight(&self, c: SpillClass, bytes: u64, slices: u32) {
        self.inflight_bytes.fetch_add(bytes, Ordering::Relaxed);
        self.inflight_slices.fetch_add(slices, Ordering::Relaxed);
        match c {
            SpillClass::Cold => {
                self.inflight_cold_bytes.fetch_add(bytes, Ordering::Relaxed);
                self.fills_cold.fetch_add(1, Ordering::Relaxed);
                self.slices_cold.fetch_add(u64::from(slices), Ordering::Relaxed);
            }
            SpillClass::Hot => {
                self.fills_hot.fetch_add(1, Ordering::Relaxed);
                self.slices_hot.fetch_add(u64::from(slices), Ordering::Relaxed);
            }
        }
    }

    /// One slice of `bytes` of class `c` landed.
    pub fn land(&self, c: SpillClass, bytes: u64) {
        self.inflight_bytes.fetch_sub(bytes, Ordering::Relaxed);
        self.inflight_slices.fetch_sub(1, Ordering::Relaxed);
        if c == SpillClass::Cold {
            self.inflight_cold_bytes.fetch_sub(bytes, Ordering::Relaxed);
        }
        self.read_bytes.fetch_add(bytes, Ordering::Relaxed);
        self.read_slices.fetch_add(1, Ordering::Relaxed);
    }

    /// Read slices issued and not yet landed.
    #[must_use]
    pub fn inflight_slices(&self) -> u32 {
        self.inflight_slices.load(Ordering::Relaxed)
    }

    /// Cold read bytes issued and not yet landed.
    #[must_use]
    pub fn inflight_cold_bytes(&self) -> u64 {
        self.inflight_cold_bytes.load(Ordering::Relaxed)
    }

    /// `frames` parsed frames entered a slot's stash, changing its charged
    /// bytes by `delta` (`IngestOutcome::stashed_bytes`).
    pub fn add_stash(&self, delta: i64, frames: u64) {
        let now = self.stash_bytes.fetch_add(delta, Ordering::Relaxed) + delta;
        self.stash_peak_bytes.fetch_max(u64::try_from(now).unwrap_or(0), Ordering::Relaxed);
        self.frames_stashed.fetch_add(frames, Ordering::Relaxed);
    }

    /// A stashed frame of `bytes` was claimed, by the merge itself when
    /// `by_consumer`.
    pub fn sub_stash(&self, bytes: u64, by_consumer: bool) {
        self.stash_bytes.fetch_sub(i64::try_from(bytes).unwrap_or(i64::MAX), Ordering::Relaxed);
        if by_consumer {
            self.claims_consumer.fetch_add(1, Ordering::Relaxed);
        } else {
            self.claims_worker.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// The planner found the awaited slot with nothing issued, stashed or
    /// read ahead.
    pub fn note_awaited_starved(&self) {
        self.awaited_starved_reads.fetch_add(1, Ordering::Relaxed);
    }

    /// Times the awaited slot was found with nothing read ahead.
    #[must_use]
    pub fn awaited_starved_reads(&self) -> u64 {
        self.awaited_starved_reads.load(Ordering::Relaxed)
    }

    /// The highest stashed (parsed, unclaimed) byte count seen.
    #[must_use]
    pub fn stash_peak_bytes(&self) -> u64 {
        self.stash_peak_bytes.load(Ordering::Relaxed)
    }

    /// Claims of stashed frames by decompress workers and by the merge.
    #[must_use]
    pub fn claims(&self) -> (u64, u64) {
        (self.claims_worker.load(Ordering::Relaxed), self.claims_consumer.load(Ordering::Relaxed))
    }

    /// The `--sort-stats` supply lines: disk reads, then stash and claims
    /// (`pool_peak_held` is the slice pool's peak of leased + idle capacity).
    #[must_use]
    pub fn supply_lines(&self, pool_peak_held: u64) -> Vec<String> {
        use super::read_ahead_budget::mib;
        let (workers, consumer) = self.claims();
        vec![
            format!(
                "Spill disk read: {} fills [hot {} / cold {}], {} slices, {}",
                self.fills_hot.load(Ordering::Relaxed) + self.fills_cold.load(Ordering::Relaxed),
                self.fills_hot.load(Ordering::Relaxed),
                self.fills_cold.load(Ordering::Relaxed),
                self.read_slices.load(Ordering::Relaxed),
                mib(self.read_bytes.load(Ordering::Relaxed)),
            ),
            format!(
                "Spill supply: frames {} stashed, peak stash {}; claims {} by workers + {} by the \
                 consumer; slice pool peak held {}; awaited-starved reads {}",
                self.frames_stashed.load(Ordering::Relaxed),
                mib(self.stash_peak_bytes.load(Ordering::Relaxed)),
                workers,
                consumer,
                mib(pool_peak_held),
                self.awaited_starved_reads(),
            ),
        ]
    }
}

/// What the merge logs about its supply at the end (`--sort-stats`): the
/// spill `Byte fetch` line, the measured spill block size beside the fixed hot
/// fill, the disk reads, and the stash and claims.
pub struct SupplyDiagnostics {
    /// The shared ledger.
    pub ledger: Arc<SupplyLedger>,
    /// The spill slice pool (peak held: leased + idle).
    pub pool: Arc<fgumi_bam_io::pread::SliceBufferPool>,
    /// `PreadSpillSlices`' request histogram.
    pub hist: Arc<crate::pread::RequestSizeHist>,
    /// The read-stream policy the spill reads adopted.
    pub policy: Arc<fgumi_bam_io::pread::ReadStreamsPolicy>,
}

impl SupplyDiagnostics {
    /// The lines, in order.
    #[must_use]
    pub fn lines(&self) -> Vec<String> {
        let l = &self.ledger;
        let slices_hot = l.slices_hot.load(Ordering::Relaxed);
        let slices_cold = l.slices_cold.load(Ordering::Relaxed);
        let frames = l.frames_stashed.load(Ordering::Relaxed);
        let read = l.read_bytes.load(Ordering::Relaxed);
        let per_block = read.checked_div(frames).unwrap_or(0);
        let hot_fill = super::read_ahead_budget::HOT_FILL_BYTES;
        let mut out = vec![
            format!(
                "Byte fetch: {} slices issued (hot {slices_hot} / cold {slices_cold}), {}",
                slices_hot + slices_cold,
                crate::pread::byte_fetch_summary(&self.hist, &self.policy),
            ),
            format!(
                "Spill blocks: {per_block} B/block measured; hot fill {} MiB (fixed) holds ~{} \
                 blocks",
                hot_fill >> 20,
                hot_fill.checked_div(per_block).unwrap_or(0),
            ),
        ];
        out.extend(l.supply_lines(self.pool.peak_held_bytes()));
        out
    }
}

/// The state the spill supply's three participants — the read planner, the
/// decompress step and the merge — must share: one merge demand (deliveries to
/// the awaited slot wake the merge; the hot set orders reads and decompression)
/// and one ledger (the planner's outstanding-read bound is released as the
/// decompress step lands slices, so two ledgers would stop reads for good).
/// Each takes this at construction, so none can be built with its own copy.
#[derive(Clone)]
pub struct SpillSupply {
    /// The merge-wide demand.
    pub demand: Arc<MergeDemand>,
    /// The shared byte ledger.
    pub ledger: Arc<SupplyLedger>,
}

impl SpillSupply {
    /// A fresh demand and ledger.
    #[must_use]
    pub fn new() -> Self {
        Self { demand: Arc::new(MergeDemand::new()), ledger: Arc::new(SupplyLedger::default()) }
    }
}

impl Default for SpillSupply {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inflight_and_stash_balance() {
        let l = SupplyLedger::default();
        l.add_inflight(SpillClass::Cold, 2 << 20, 2);
        l.add_inflight(SpillClass::Hot, 4 << 20, 4);
        assert_eq!((l.inflight_slices(), l.inflight_cold_bytes()), (6, 2 << 20));
        l.land(SpillClass::Cold, 1 << 20);
        assert_eq!((l.inflight_slices(), l.inflight_cold_bytes()), (5, 1 << 20));
        l.sub_stash(100, false);
        assert_eq!(l.stash_peak_bytes(), 0, "a claim that precedes its ingest's tally");
        l.add_stash(300, 3);
        assert_eq!(l.stash_peak_bytes(), 200);
        l.sub_stash(100, false);
        l.sub_stash(100, true);
        assert_eq!(l.claims(), (2, 1));
        let lines = l.supply_lines(8 << 20);
        assert_eq!(lines[0], "Spill disk read: 2 fills [hot 1 / cold 1], 1 slices, 1 MiB");
        assert!(lines[1].starts_with("Spill supply: frames 3 stashed"), "{}", lines[1]);
        assert!(lines[1].contains("claims 2 by workers + 1 by the consumer"), "{}", lines[1]);
        assert!(lines[1].contains("slice pool peak held 8 MiB"), "{}", lines[1]);
    }
}
