//! The merge's spill read-ahead budget.
//!
//! The read-ahead target `R` is derived from `--max-memory` (`total / 16`,
//! clamped to 64 MiB ..= 512 MiB). It is **not** subtracted from the sort
//! buffer, so it adds to it: the final run is merged from memory (it is never
//! spilled), so the merge's read-ahead coexists with up to one full sort
//! buffer. Charging `R` to the buffer instead would force earlier spills —
//! more runs, more merge work — for every sort, to save at most `R` at the
//! merge.
//!
//! Every live slot may always hold its **cold terms**: a cold allowance, an
//! even share of `R` after the hot reserve (`(R − HOT_SLOTS × hot) / k`,
//! clamped to one 256 KiB storage IOP ..= 2 MiB), read in cold fills.
//! Everything a slot holds above its cold allowance is charged to **one shared
//! pool**, `P = max(R − k × cold, HOT_SLOTS × hot)`: the slots the merge needs
//! next (the hot set, at most [`HOT_SLOTS`]) top up toward the hot allowance
//! only while the pool has room, and a slot that leaves the hot set keeps its
//! charge until the merge drains it — so however often the hot set changes,
//! read-ahead never grows past `k × cold + P`. The awaited slot can always
//! read within its cold terms, so a full pool never starves the merge. Idle
//! slice buffers waiting for reuse are charged to the pool too, and so are
//! decompressed blocks above the cold FIFO cap: a hot slot's FIFO cap rises
//! from [`FIFO_CAP_COLD`] to [`FIFO_CAP_HOT`] only while the pool has room for
//! the difference, and a slot demoted back to the cold cap stays charged for
//! the blocks it still holds above it.
//!
//! A slot is charged what it holds (`SortMergeSlot::stash_bytes` plus its
//! requested bytes): every read slice with a frame still stashed or being
//! decompressed (in full — a borrowed frame pins its whole slice). Its carried
//! partial frame is left out of the charge, because it completes only with the
//! slot's next read and charging it could block that read; it is the one
//! frame per slot the ceiling adds.
//! The peak above the sort buffer is therefore at most
//!
//! ```text
//! k × (cold + one carried frame) + P + decoded
//! ```
//!
//! with `decoded` the cold FIFOs, `k × FIFO_CAP_COLD` blocks
//! ([`ReadAheadBudget::decoded_ceiling`]; every block above that is in the
//! pool).
//! Since `k × cold + P = max(R, k × cold + HOT_SLOTS × hot)`, the read-ahead is
//! `R` unless the IOP floor (very large `k`) or the hot reserve forces more;
//! the `Spill supply` budget line reports every term.

use crate::pread::SpillClass;

/// Slots that can be hot at once: awaited, predicted, frontier.
pub const HOT_SLOTS: u64 = 3;
/// Bytes per hot read (one fill, split into the read-stream count's slices).
pub const HOT_FILL_BYTES: u64 = 4 << 20;
/// Read-ahead (requested + held) a hot slot is topped up to, pool permitting.
pub const HOT_ALLOWANCE_BYTES: u64 = 16 << 20;
/// Smallest cold read: one gp3 IOP (requests ≤ 256 KiB bill as one).
pub const COLD_FILL_MIN_BYTES: u64 = 256 << 10;
/// Largest cold allowance: two 1 MiB fills.
pub const COLD_ALLOWANCE_MAX_BYTES: u64 = 2 << 20;
/// Cap on outstanding cold read bytes, so hot reads never queue behind a wall
/// of cold ones.
pub const COLD_INFLIGHT_BYTES: u64 = 8 << 20;
/// Lower bound of `R`.
pub const R_FLOOR: u64 = 64 << 20;
/// Upper bound of `R`.
pub const R_CEIL: u64 = 512 << 20;
/// `R = total memory / R_DIVISOR`, clamped.
pub const R_DIVISOR: u64 = 16;
/// Decompressed FIFO cap of a hot slot.
pub const FIFO_CAP_HOT: u32 = 32;
/// Decompressed FIFO cap of a cold slot: 8 blocks (≈ 512 KiB of 64 KiB
/// blocks). A cold slot is not next in the merge, so it needs only enough
/// decoded runway to cover the time the merge takes to promote it to the hot
/// set; every block above that is memory held across all `k − 3` cold slots.
pub const FIFO_CAP_COLD: u32 = 8;
/// The carried partial frame a slot may hold above its allowance: a BGZF
/// block is at most 64 KiB, and a zstd frame holds one ~64 KiB sort block
/// (larger only for a single record longer than that).
pub const CARRIED_FRAME_BYTES: u64 = 64 << 10;
/// Bytes of one decompressed block, for the decoded ceiling.
pub const DECOMPRESSED_BLOCK_BYTES: u64 = 64 << 10;

/// The pool charge of a decompressed FIFO of `blocks` (its cap, or its length
/// when that is larger): the blocks above [`FIFO_CAP_COLD`].
#[must_use]
pub fn fifo_charge(blocks: u64) -> u64 {
    blocks.saturating_sub(u64::from(FIFO_CAP_COLD)) * DECOMPRESSED_BLOCK_BYTES
}

/// The resolved split of `R` over `k` merge slots.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadAheadBudget {
    /// The read-ahead target.
    pub r: u64,
    /// Merge slots (spill runs) the budget is split over.
    pub k: usize,
    /// Hot slot allowance.
    pub hot_allowance: u64,
    /// Hot fill size.
    pub hot_fill: u64,
    /// Cold slot allowance.
    pub cold_allowance: u64,
    /// Cold fill size.
    pub cold_fill: u64,
    /// The shared pool every byte held above a cold allowance is charged to.
    pub pool: u64,
}

impl ReadAheadBudget {
    /// Split `R` (from `total_memory`) over `k` slots.
    #[must_use]
    pub fn resolve(total_memory: u64, k: usize) -> Self {
        let r = (total_memory / R_DIVISOR).clamp(R_FLOOR, R_CEIL);
        let cold_pool = r.saturating_sub(HOT_SLOTS * HOT_ALLOWANCE_BYTES);
        let cold_allowance =
            (cold_pool / k.max(1) as u64).clamp(COLD_FILL_MIN_BYTES, COLD_ALLOWANCE_MAX_BYTES);
        let cold_fill = (cold_allowance / 2).max(COLD_FILL_MIN_BYTES).min(cold_allowance);
        Self::with_terms(r, k, HOT_ALLOWANCE_BYTES, HOT_FILL_BYTES, cold_allowance, cold_fill)
    }

    /// A budget of these terms, with its pool derived from them.
    fn with_terms(
        r: u64,
        k: usize,
        hot_allowance: u64,
        hot_fill: u64,
        cold_allowance: u64,
        cold_fill: u64,
    ) -> Self {
        let pool = r.saturating_sub(k as u64 * cold_allowance).max(HOT_SLOTS * hot_allowance);
        Self { r, k, hot_allowance, hot_fill, cold_allowance, cold_fill, pool }
    }

    /// This budget with fills of `hot_fill` / `cold_fill` bytes and the design's
    /// allowance ratios (hot four fills, cold two), for tests that drive the
    /// supply with fills small enough that frames straddle them.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn override_for_test(self, hot_fill: u64, cold_fill: u64) -> Self {
        Self::with_terms(self.r, self.k, 4 * hot_fill, hot_fill, 2 * cold_fill, cold_fill)
    }

    /// This budget with its pool replaced (test support: make the pool bind).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_pool_for_test(self, pool: u64) -> Self {
        Self { pool, ..self }
    }

    /// `(allowance, fill)` for a slot of class `c`.
    #[must_use]
    pub fn for_class(&self, c: SpillClass) -> (u64, u64) {
        match c {
            SpillClass::Hot => (self.hot_allowance, self.hot_fill),
            SpillClass::Cold => (self.cold_allowance, self.cold_fill),
        }
    }

    /// The bound on charged read-ahead: `k × cold + P` (every slot's cold
    /// terms plus the shared pool).
    #[must_use]
    pub fn read_ahead_bound(&self) -> u64 {
        self.k as u64 * self.cold_allowance + self.pool
    }

    /// Worst-case decompressed FIFO bytes outside the pool:
    /// `k × FIFO_CAP_COLD × 64 KiB` (blocks above the cold cap are pool
    /// charges).
    #[must_use]
    pub fn decoded_ceiling(&self) -> u64 {
        self.k as u64 * u64::from(FIFO_CAP_COLD) * DECOMPRESSED_BLOCK_BYTES
    }

    /// The merge's worst-case memory above the sort buffer:
    /// `k × (cold + one carried frame) + P + decoded`.
    #[must_use]
    pub fn ceiling(&self) -> u64 {
        self.read_ahead_bound() + self.k as u64 * CARRIED_FRAME_BYTES + self.decoded_ceiling()
    }

    /// The `--sort-stats` budget line.
    #[must_use]
    pub fn budget_line(&self) -> String {
        format!(
            "Spill supply: read-ahead budget R={} over {} slots (derived from --max-memory, adds \
             to the sort buffer) -> cold {} per slot (fill {}, FIFO {}), hot pool {} (up to {} \
             per hot slot, fill {}, FIFO {}); peak above the sort buffer <= {} = {} x ({} + {} \
             carried frame) + pool {} + decoded {}",
            mib(self.r),
            self.k,
            mib(self.cold_allowance),
            mib(self.cold_fill),
            FIFO_CAP_COLD,
            mib(self.pool),
            mib(self.hot_allowance),
            mib(self.hot_fill),
            FIFO_CAP_HOT,
            mib(self.ceiling()),
            self.k,
            mib(self.cold_allowance),
            kib(CARRIED_FRAME_BYTES),
            mib(self.pool),
            mib(self.decoded_ceiling()),
        )
    }
}

/// `n` bytes in MiB with one decimal (`"1.5 MiB"`), or whole when exact.
#[allow(clippy::cast_precision_loss, reason = "a diagnostic rendering")]
pub(crate) fn mib(n: u64) -> String {
    if n.is_multiple_of(1 << 20) {
        format!("{} MiB", n >> 20)
    } else {
        format!("{:.1} MiB", n as f64 / f64::from(1u32 << 20))
    }
}

/// `n` bytes as whole KiB (`"64 KiB"`).
fn kib(n: u64) -> String {
    format!("{} KiB", n >> 10)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    const MIB: u64 = 1 << 20;
    const KIB: u64 = 1 << 10;

    /// The design table's rows, with the expected integers worked by hand:
    /// `(512 MiB − 48 MiB) / 350 = 1_390_112` (fill = half = `695_056`);
    /// `/ 1024 = 475_136` (fill clamps up to 256 KiB). The pool is
    /// `max(R − k × cold, 48 MiB)`: 512 − 54 = 458 MiB at k = 27; 512 − 178 =
    /// 334 MiB at k = 89; `536_870_912 − 350 × 1_390_112 = 50_331_712` at
    /// k = 350; 512 − 464 = 48 MiB at k = 1024; 256 − 178 = 78 MiB for lowmem;
    /// the 48 MiB reserve for both sweep rows (`k × cold` exceeds `R`).
    #[rstest]
    #[case::t16_k27(768 << 20, 16, 27, 512 << 20, 2 << 20, 1 << 20, 458 * MIB)]
    #[case::t16_k89(768 << 20, 16, 89, 512 << 20, 2 << 20, 1 << 20, 334 * MIB)]
    #[case::t16_k350(768 << 20, 16, 350, 512 << 20, 1_390_112, 695_056, 50_331_712)]
    #[case::t16_k1024(768 << 20, 16, 1024, 512 << 20, 475_136, 256 << 10, 48 * MIB)]
    #[case::lowmem_k89(512 << 20, 8, 89, 256 << 20, 2 << 20, 1 << 20, 78 * MIB)]
    #[case::sweep128_k350(128 << 20, 8, 350, 64 << 20, 256 << 10, 256 << 10, 48 * MIB)]
    #[case::sweep64_k700(64 << 20, 8, 700, 64 << 20, 256 << 10, 256 << 10, 48 * MIB)]
    fn budget_rows_match_the_design_table(
        #[case] per_thread: u64,
        #[case] t: u64,
        #[case] k: usize,
        #[case] r: u64,
        #[case] cold_allowance: u64,
        #[case] cold_fill: u64,
        #[case] pool: u64,
    ) {
        let b = ReadAheadBudget::resolve(per_thread * t, k);
        assert_eq!(
            (b.r, b.cold_allowance, b.cold_fill, b.pool),
            (r, cold_allowance, cold_fill, pool)
        );
        assert!(b.cold_fill >= COLD_FILL_MIN_BYTES);
    }

    /// The ceiling worked by hand: `k × (cold + 64 KiB) + P + k × 8 × 64 KiB`
    /// (hot FIFO blocks above the cold cap are in the pool).
    #[rstest]
    #[case::t16_k27(768 << 20, 16, 27, 27 * (2 * MIB + 64 * KIB) + 458 * MIB + 27 * 512 * KIB)]
    #[case::t16_k1024(
        768 << 20,
        16,
        1024,
        1024 * (475_136 + 64 * KIB) + 48 * MIB + 1024 * 512 * KIB
    )]
    fn ceiling_matches_the_design_formula(
        #[case] per_thread: u64,
        #[case] t: u64,
        #[case] k: usize,
        #[case] ceiling: u64,
    ) {
        assert_eq!(ReadAheadBudget::resolve(per_thread * t, k).ceiling(), ceiling);
    }

    /// `k × cold + P = max(R, k × cold + 48 MiB)`, so the read-ahead bound
    /// exceeds `R` only when the IOP floor forces `k × cold` past
    /// `R − 48 MiB`: beyond `k = (R − 48 MiB) / 256 KiB` = 64, 832, 1856 for
    /// `R` = 64, 256, 512 MiB.
    #[rstest]
    #[case::r64(64 << 20, 64)]
    #[case::r256(256 << 20, 832)]
    #[case::r512(512 << 20, 1856)]
    fn bound_exceeds_r_only_beyond_the_documented_k(#[case] r: u64, #[case] threshold: usize) {
        for k in [1, 3, 27, threshold - 1, threshold, threshold + 1, threshold * 2] {
            let b = ReadAheadBudget::resolve(r * R_DIVISOR, k);
            assert_eq!(b.r, r);
            assert_eq!(b.read_ahead_bound() > r, k > threshold, "R={r} k={k}");
            assert_eq!(
                b.read_ahead_bound(),
                r.max(k as u64 * b.cold_allowance + HOT_SLOTS * HOT_ALLOWANCE_BYTES)
            );
        }
    }

    /// `350 × 8 × 64 KiB = 183_500_800` B (175 MiB) at k = 350.
    #[test]
    fn budget_rows_decoded_ceiling_at_k350_is_175_mib() {
        let b = ReadAheadBudget::resolve(768 << 24, 350);
        assert_eq!(b.decoded_ceiling(), 183_500_800);
    }

    /// At k = 1024 with the defaults (t16 × 768 MiB): `1024 × (475_136 +
    /// 65_536) = 553_648_128` + pool 48 MiB + decoded `1024 × 512 KiB =
    /// 512 MiB` = `1_140_850_688` B (1088 MiB ≈ 1.06 GiB), under 1.1 GiB.
    #[test]
    fn k1024_default_ceiling_is_at_most_1_1_gib() {
        let b = ReadAheadBudget::resolve(768 << 24, 1024);
        assert_eq!(b.ceiling(), 1_140_850_688);
        assert!(b.ceiling() * 10 <= 11 << 30, "{} B > 1.1 GiB", b.ceiling());
    }

    /// The FIFO charge is the blocks above the cold cap.
    #[test]
    fn fifo_charge_counts_blocks_above_the_cold_cap() {
        assert_eq!(fifo_charge(0), 0);
        assert_eq!(fifo_charge(u64::from(FIFO_CAP_COLD)), 0);
        assert_eq!(fifo_charge(u64::from(FIFO_CAP_HOT)), 24 * 64 * KIB);
    }

    #[test]
    fn budget_line_names_target_pool_and_ceiling() {
        let line = ReadAheadBudget::resolve(768 << 24, 27).budget_line();
        assert_eq!(
            line,
            "Spill supply: read-ahead budget R=512 MiB over 27 slots (derived from --max-memory, \
             adds to the sort buffer) -> cold 2 MiB per slot (fill 1 MiB, FIFO 8), hot pool 458 \
             MiB (up to 16 MiB per hot slot, fill 4 MiB, FIFO 32); peak above the sort buffer <= \
             527.2 MiB = 27 x (2 MiB + 64 KiB carried frame) + pool 458 MiB + decoded 13.5 MiB"
        );
    }
}
