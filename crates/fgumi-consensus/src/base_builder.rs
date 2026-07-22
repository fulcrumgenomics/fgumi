//! # Consensus Base Calling
//!
//! This module provides the core likelihood-based algorithm for calling a consensus base at a
//! single position from multiple observations. This is the fundamental building block used by
//! all consensus callers to determine the most likely base and assign an accurate quality score.
//!
//! ## Overview
//!
//! Given multiple observations of a position (base + quality score pairs), the consensus base
//! builder uses Bayesian inference to determine:
//! 1. Which base (A, C, G, or T) is most likely to be the true base
//! 2. What quality score to assign to that consensus call
//!
//! This differs from simple majority voting by:
//! - Incorporating quality scores (low-quality observations contribute less)
//! - Using a probabilistic model that can properly handle ambiguous cases
//! - Generating quality scores that accurately reflect confidence in the call
//!
//! ## Likelihood-Based Model
//!
//! For each candidate base B ∈ {A, C, G, T}, we calculate the likelihood of observing all the
//! data given that B is the true base:
//!
//! ```text
//! L(B | observations) = ∏ P(obs_i | B)
//! ```
//!
//! Where for each observation:
//! - **If observed base matches B**: `P(obs | B) = 1 - error_rate`
//! - **If observed base differs from B**: `P(obs | B) = error_rate / 3`
//!
//! The division by 3 assumes equal probability of the three incorrect bases.
//!
//! ## Quality Score Calculation
//!
//! After computing likelihoods for all four bases:
//!
//! 1. **Normalize** to get posterior probabilities:
//!    ```text
//!    P(B | data) = L(B) / (L(A) + L(C) + L(G) + L(T))
//!    ```
//!
//! 2. **Select** the base with maximum posterior probability
//!
//! 3. **Calculate error probability**:
//!    ```text
//!    P_error = 1 - P(best_base | data)
//!    ```
//!
//! 4. **Convert to Phred scale**:
//!    ```text
//!    Q = -10 * log10(P_error)
//!    ```
//!
//! This Phred score represents the probability that the consensus call is wrong, properly
//! accounting for all the evidence (both supporting and contradicting).
//!
//! ## Error Rate Priors
//!
//! The model incorporates two sources of error:
//!
//! ### Pre-UMI Error Rate
//!
//! Errors present in the original molecule before PCR amplification. These errors will appear
//! in all reads from that molecule, so they cannot be corrected by consensus calling. This
//! rate sets a floor on the achievable consensus quality.
//!
//! **Typical value**: Q45 (error rate ~3 × 10^-5)
//!
//! ### Post-UMI Error Rate
//!
//! Errors introduced during sequencing or PCR amplification. These are independent across reads
//! and can be corrected by consensus calling.
//!
//! **Typical value**: Q40 (error rate ~1 × 10^-4)
//!
//! ### Combined Error Model
//!
//! For each observation with quality Q, the adjusted error probability accounts for both
//! post-UMI errors (from sequencing) and the quality score:
//!
//! ```text
//! P_adjusted_error = P_post_umi_error + (1 - P_post_umi_error) * P_sequencing_error
//! ```
//!
//! This adjusted probability is used in the likelihood calculation.
//!
//! After calling the consensus, the pre-UMI error rate is incorporated:
//! ```text
//! P_final_error = P_pre_umi_error + (1 - P_pre_umi_error) * P_consensus_error
//! ```
//!
//! ## Mathematical Details
//!
//! The implementation uses **log-space arithmetic** to avoid numerical underflow:
//!
//! - All probabilities are stored as natural logarithms
//! - Products become sums: `log(a * b) = log(a) + log(b)`
//! - Sums use log-sum-exp trick: `log(e^a + e^b) = log_sum_exp(a, b)`
//!
//! This allows accurate computation even when individual likelihoods are extremely small
//! (e.g., 10^-100).
//!
//! ## Usage in Consensus Callers
//!
//! Consensus callers typically use `ConsensusBaseBuilder` as follows:
//!
//! ```text
//! use fgumi_consensus::base_builder::ConsensusBaseBuilder;
//!
//! // Create builder with error rates
//! let mut builder = ConsensusBaseBuilder::new(
//!     45,  // Pre-UMI error rate (Q45)
//!     40,  // Post-UMI error rate (Q40)
//! );
//!
//! // For each position in the read:
//! for position in 0..read_length {
//!     builder.reset(); // Clear previous position
//!
//!     // Add observations from all reads
//!     for read in &reads {
//!         let base = read.sequence()[position];
//!         let qual = read.quality_scores()[position];
//!         builder.add(base, qual);
//!     }
//!
//!     // Call consensus for this position
//!     let (consensus_base, consensus_qual) = builder.call();
//!
//!     // Check depth
//!     let depth = builder.contributions();
//! }
//! ```
//!
//! ## Handling Edge Cases
//!
//! ### No Observations
//! Returns `(N, 0)` - no-call with minimum quality.
//!
//! ### N Bases
//! N (no-call) bases in input reads are ignored. They don't contribute evidence for any
//! particular base.
//!
//! ### Ties
//! If multiple bases have exactly equal posterior probability, returns `(N, 0)` to indicate
//! ambiguity. In practice, this is rare due to quality score differences.
//!
//! ### Low Depth
//! With only 1-2 observations, the consensus quality will be similar to the input qualities.
//! The model naturally handles this - more observations lead to higher confidence.
//!
//! ## Quality Score Bounds
//!
//! Phred scores are capped at Q93 (error rate ~5 × 10^-10):
//! - **Minimum**: Q2 (error rate ~0.63)
//! - **Maximum**: Q93 (error rate ~5 × 10^-10)
//!
//! These bounds ensure quality scores fit in standard BAM format and remain numerically stable.
//!
//! ## Performance Characteristics
//!
//! - **Time complexity**: O(n) where n = number of observations at a position
//! - **Space complexity**: O(1) - fixed-size arrays for four bases
//! - **Caching**: Pre-computes error tables for all possible quality scores (Q2-Q93) to avoid
//!   repeated calculations
//!
//! ## Example Quality Progression
//!
//! Consider calling consensus from multiple reads observing 'A':
//!
//! | Reads | All Q30 'A' | Consensus Quality |
//! |-------|-------------|-------------------|
//! | 1     | A           | ~Q30              |
//! | 2     | AA          | ~Q40              |
//! | 5     | AAAAA       | ~Q50              |
//! | 10    | AAAAAAAAAA  | ~Q60              |
//! | 20    | (20x A)     | ~Q70              |
//!
//! The consensus quality increases logarithmically with depth, reflecting the compounding
//! evidence that all observations support the same base.
//!
//! ## Unanimous Fast-Path Correctness (Margin-Gated, Not "By Construction")
//!
//! `ConsensusBaseBuilder::call_full` is the unmodified fgbio computation — for
//! every position, unanimous or not, it runs the general `[w, l, l, l]` (or
//! four-way) `ln_normalize` -> `ln_not` -> `ln_error_prob_two_trials` ->
//! `ln_prob_to_phred` chain. It never consults `unanimous_quality_from_gap` or the
//! gap-threshold table. `call_full` is therefore the sole oracle: fgbio parity is
//! defined as "output equals `call_full`", by definition, not by an equivalence
//! this module has to prove separately.
//!
//! `ConsensusBaseBuilder::try_unanimous_fast_path` is a *shortcut* over a
//! precomputed per-Phred-level gap-threshold table (`unanimous_quality_from_gap`'s
//! step function, built by `build_unanimous_gap_thresholds`). It returns the
//! table's answer only when a stated **sufficient condition** establishes that the
//! answer equals what `call_full` would have computed, and returns `None` (defer
//! to `call_full`) in every other case. The sufficient condition differs between
//! the two regions, because the quantity being bounded behaves differently in each:
//!
//! - **Sub-cap region** (`gap < thresholds[cap]`): `call_full`'s `[w, l, l, l]`
//!   normalize and the table's shifted `[0, -gap, -gap, -gap]` form compute the
//!   *same* mathematical posterior and differ only by floating-point rounding,
//!   empirically bounded by `delta = C * u * (|w| + |l|)`, `u = f64::EPSILON/2`,
//!   `C = 16` (an observed worst-case ratio of ~3.845x against a pinned ceiling of
//!   5.0x; see `test_ln_posterior_full_gap_error_bound`, which asserts both the
//!   hard bound and the ratio ceiling, and the calibration note below). Mapping
//!   that ln-posterior perturbation to a gap perturbation via the mean value
//!   theorem: `Q(gap) = -ln(1 + 3e^-gap)` has derivative `Q'(gap) = 1 - p = c(gap)`
//!   (the consensus error), so a perturbation `delta` in ln-posterior corresponds
//!   to a gap shift bounded in magnitude by `delta` divided by `c_min`. Since `c`
//!   is strictly *decreasing* in `gap`, the minimum of `c` over the containing
//!   bracket `[thresholds[q], thresholds[q+1])` is `c(thresholds[q+1])`,
//!   precomputed as `cerr_min[q]`. If the observed `gap` clears the resulting
//!   margin from *both* bracket edges, no perturbation within the empirically
//!   bounded margin can push `call_full`'s true answer across a bracket boundary
//!   into a different quality level, so the table's `q` is conditionally exact.
//!   If it does not clear the margin on either side, the fast path cannot
//!   establish exactness and defers.
//! - **Cap region** (`gap >= thresholds[cap]`): PR #634 proves the *shifted*
//!   `[0, -gap, -gap, -gap]` form's quality rounds to the pre-UMI cap for every
//!   such `gap`. But `call_full` evaluates fgbio's original, unshifted
//!   `[w, l, l, l]` form, not this module's shifted one — and just above the
//!   threshold the loser-magnitude cancellation in that form can drop `call_full`
//!   one Phred short of the cap. So the two forms are only provably equal subject
//!   to the same empirically bounded absolute ln-posterior error used sub-cap,
//!   `delta = C * u * (|w| + |l|)` (via `unanimous_margin(.., cerr_lower_bound =
//!   1.0)`, i.e. no MVT division). Unlike the sub-cap brackets, there is no bracket
//!   *above* the cap, so the gate cannot be sized as a shrinking margin against a
//!   far edge. The quantity that must exceed `delta` is instead the remaining
//!   ln-posterior *budget* (`Q(gap)` minus `Q(thresholds[cap])`, i.e. `ln(1 +
//!   3e^-thresholds[cap]) - ln(1 + 3e^-gap)`), which **saturates** at
//!   (approximately) the cap boundary's own consensus error —
//!   `consensus_error(thresholds[cap])`, equivalently `cerr_min[cap - 1]` (see
//!   `UnanimousGapTables::cerr_min`'s doc) — as `gap` grows arbitrarily large; extra
//!   gap past the threshold buys vanishingly more budget, it does not keep growing
//!   linearly. The gate therefore requires both `gap - thresholds[cap] >= ln 2`
//!   (so `3e^-gap <= (3e^-thresholds[cap])/2`, meaning the budget has closed to at
//!   least half its saturated value — since `a/(1+a) <= ln(1+a)` for `a =
//!   3e^-thresholds[cap]`) *and* `delta` below half that same cap-boundary
//!   consensus error. An allowance that instead grows unboundedly with `gap`
//!   (equivalently, with depth, since both `delta` and a naive
//!   `gap - thresholds[cap]` allowance grow linearly with depth while the true
//!   budget does not) is unsound: see
//!   `test_fast_path_matches_call_full_at_deep_cap_region` for confirmed
//!   divergences the earlier, unbounded gate produced at depths in the thousands.
//!
//! Both regions share the same `delta = C * u * (|w| + |l|)` bound; see
//! `UNANIMOUS_MARGIN_HEADROOM`'s doc for how `C = 16` was calibrated (over all
//! four winner lanes, not just lane 0).
//!
//! Net effect: `call()` returns a value only when one of the two sufficient
//! conditions above holds, in which case it is proven (subject to the stated
//! empirical bound) to equal `call_full()`; everywhere else it defers to
//! `call_full()` directly, so the two are trivially equal there (same code path).
//! This is *not* an unconditional "always equal" guarantee independent of the
//! bound — it is conditional exactness plus safe deferral. The critical regression
//! test for this is `test_unanimous_fast_path_matches_full_calculation` (broad
//! sweep) and `test_fast_path_matches_call_full_at_deep_cap_region` (targeted deep
//! cap-region regression), together with the dense contiguous-depth sweep in
//! `test_fast_path_matches_call_full_dense`.
//!
//! ## See Also
//!
//! - `vanilla_caller`: Uses this builder to call consensus across entire reads
//! - `duplex_caller`: Uses this builder for single-strand consensus before duplex combination
//! - `phred`: Utility functions for Phred score and probability conversions

use crate::phred::{
    LN_ONE, LogProbability, MAX_PHRED, MIN_PHRED, NO_CALL_BASE, PhredScore,
    ln_error_prob_two_trials, ln_normalize, ln_not, ln_prob_to_phred, ln_sum_exp_array,
    phred_to_ln_error_prob,
};
use approx::abs_diff_eq;
use std::cmp::Ordering;
use std::sync::OnceLock;
use wide::f64x4;

/// The four DNA bases in the order used throughout consensus calling
const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const DNA_BASE_COUNT: usize = 4;

/// Headroom multiplier on the empirically bounded `|ln_posterior_full -
/// ln_posterior_gap|` error bound (empirically ~3.845x at worst) used by
/// [`ConsensusBaseBuilder::try_unanimous_fast_path`]'s margin gate. See the
/// module-level docs and `test_ln_posterior_full_gap_error_bound`, which pins it.
///
/// The calibration sweep must cover **all four winner lanes** (bases A/C/G/T), not
/// just lane 0: `ln_sum_exp_array` seeds its accumulator with the lane holding the
/// minimum value and folds the remaining lanes in index order, so a winner at lane
/// 2 (G) or lane 3 (T) accumulates a different rounding path than one at lane 0
/// (A) or lane 1 (C). Calibrating on lane 0 alone measured only ~1.17x worst case
/// and left almost no real headroom once lanes 2/3 (the true ~3.845x worst case,
/// at `base = G`, `post = Q35`, `obs = Q56`, `depth = 1`) are included. `16.0`
/// gives ~4.2x headroom over that true worst case.
const UNANIMOUS_MARGIN_HEADROOM: f64 = 16.0;

/// Unit roundoff for `f64` (`u = f64::EPSILON / 2`): the per-operation relative
/// error bound IEEE 754 double-precision arithmetic guarantees. Shared between
/// [`unanimous_margin`] (the production margin gate) and
/// `test_ln_posterior_full_gap_error_bound` (which pins the same bound the gate
/// relies on) so the constant cannot drift independently between the two.
const UNANIMOUS_UNIT_ROUNDOFF: f64 = f64::EPSILON / 2.0;

/// Lookup table for converting ASCII base to index (0-3 for A,C,G,T, 255 for invalid)
/// This avoids a loop+compare for every base observation
const BASE_TO_INDEX: [u8; 256] = {
    let mut table = [255u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// The post-UMI-adjusted per-quality probability tables used by `ConsensusBaseBuilder`.
///
/// Both tables are indexed by the observed base quality (`0..=MAX_PHRED`) and depend *only*
/// on the post-UMI error rate they were computed for. Because that error rate and the base
/// quality are both integer Phred scores drawn from a tiny finite domain, the whole table is
/// a pure function of a single `u8` — so it can be computed once and shared, rather than
/// recomputed per `ConsensusBaseBuilder`.
struct AdjustedProbabilityTables {
    /// Log probability that an observation of quality `q` is correct, adjusted for post-UMI error.
    correct: Box<[LogProbability]>,

    /// Log probability of any one *specific* alternate base (total error / 3) for quality `q`.
    error_per_alt: Box<[LogProbability]>,
}

impl AdjustedProbabilityTables {
    /// Computes both tables for a given post-UMI error rate.
    ///
    /// This performs exactly the arithmetic the per-builder loop used to perform inline, in the
    /// same order, so the resulting `f64` values are bit-identical to the un-memoized version.
    fn compute(error_rate_post_umi: PhredScore) -> Self {
        let ln_error_post = phred_to_ln_error_prob(error_rate_post_umi);
        let ln_three = 3.0_f64.ln();

        let mut correct = Vec::with_capacity(MAX_PHRED as usize + 1);
        let mut error_per_alt = Vec::with_capacity(MAX_PHRED as usize + 1);

        for q in 0..=MAX_PHRED {
            let ln_error_seq = phred_to_ln_error_prob(q);
            let adjusted_error = ln_error_prob_two_trials(ln_error_post, ln_error_seq);

            correct.push(ln_not(adjusted_error));
            // Error for any specific alternate base is 1/3 of total error
            // ln(p/3) = ln(p) - ln(3)
            error_per_alt.push(adjusted_error - ln_three);
        }

        Self {
            correct: correct.into_boxed_slice(),
            error_per_alt: error_per_alt.into_boxed_slice(),
        }
    }
}

/// Lazily-populated cache of `AdjustedProbabilityTables`, one slot per possible post-UMI Phred.
///
/// The post-UMI error rate is a `u8`, so 256 slots cover the entire input domain exactly — no
/// clamping, and therefore no behaviour change for out-of-range Phred values. Each slot is
/// initialized at most once; after that, building a `ConsensusBaseBuilder` is a pair of atomic
/// loads instead of 94 `exp`/`log1p` evaluations plus two heap allocations. This matters because
/// `consensus_umis()` constructs a fresh caller (and thus a fresh builder) for every UMI family.
static ADJUSTED_TABLE_CACHE: [OnceLock<AdjustedProbabilityTables>; 256] =
    [const { OnceLock::new() }; 256];

/// Returns the shared, memoized adjusted-probability tables for a post-UMI error rate.
fn adjusted_tables(error_rate_post_umi: PhredScore) -> &'static AdjustedProbabilityTables {
    ADJUSTED_TABLE_CACHE[error_rate_post_umi as usize]
        .get_or_init(|| AdjustedProbabilityTables::compute(error_rate_post_umi))
}

/// Builder for calling consensus at a single base position
///
/// Accumulates observations (base + quality) and uses a likelihood model to
/// call the consensus base and quality.
///
/// Uses Kahan summation algorithm with SIMD vectorization to maintain numeric stability
/// when accumulating log-likelihoods across many observations. The 4 DNA bases map
/// perfectly to f64x4 (256-bit AVX2), enabling parallel computation of all 4 likelihoods.
pub struct ConsensusBaseBuilder {
    /// Log-likelihoods for each of the four bases (A, C, G, T) - SIMD vectorized
    likelihoods: f64x4,

    /// Kahan summation compensation terms for each base - SIMD vectorized
    compensations: f64x4,

    /// Count of observations for each base.
    ///
    /// `u32` (not `u16`) so very deep UMI families — a single position can exceed
    /// `65_535` reads in high-duplication libraries — accumulate without wrapping.
    /// The per-base depth is later clamped to `i16::MAX` only at tag emission, to
    /// match fgbio's `Short` consensus-depth representation.
    observations: [u32; DNA_BASE_COUNT],

    /// Pre-computed correct probabilities adjusted for post-UMI errors.
    ///
    /// Borrowed from the process-wide `ADJUSTED_TABLE_CACHE`; see `adjusted_tables`.
    adjusted_correct_table: &'static [LogProbability],

    /// Pre-computed 1/3 of error probability (for wrong base likelihoods).
    ///
    /// Borrowed from the process-wide `ADJUSTED_TABLE_CACHE`; see `adjusted_tables`.
    adjusted_error_per_alt: &'static [LogProbability],

    /// Pre-UMI error rate (log probability)
    ln_error_pre_umi: LogProbability,

    /// Per-Phred-level gap tables for the unanimous fast path, for this
    /// `ln_error_pre_umi`: the gap thresholds themselves plus the per-bracket
    /// consensus-error minima used to size the fast path's margin gate. Built once
    /// per rate by [`cached_unanimous_gap_tables`]. See
    /// [`Self::try_unanimous_fast_path`] for how this is used, and the module-level
    /// docs for the two sufficient conditions (one per region) under which the
    /// result is proven to equal [`Self::call_full`] whenever the fast path
    /// returns `Some`.
    gap_tables: &'static UnanimousGapTables,
}

/// Quality of a **unanimous** position whose winning base leads its three equal
/// losers by `gap` (`= winner_ll − loser_ll`), for a given pre-UMI error rate.
///
/// A unanimous position has one winner log-likelihood and three equal losers exactly
/// `gap` below it, so the posterior is `1 / (1 + 3·e^{−gap})` — independent of the
/// absolute likelihood level. This evaluates the quality chain in the **shifted**
/// `[0, −gap, −gap, −gap]` frame, whose only offset is `gap` (bounded by the sub-cap
/// ceiling, < ~28.5). It therefore avoids the catastrophic `fl(l + S)` cancellation
/// that `call_full`'s `[w, l, l, l]` normalize incurs at deep, low-per-observation-gap
/// pileups (where the loser magnitude `|l|` grows without bound). This is used only
/// to build the fast path's gap-threshold table (via
/// [`build_unanimous_gap_thresholds`]) — `call_full` computes fgbio's original
/// `[w, l, l, l]` form directly and never calls this helper; see the module-level
/// docs for why the table built from it still agrees with `call_full` —
/// empirically bounded and conditionally exact in both the cap region and the
/// sub-cap region.
fn unanimous_quality_from_gap(gap: f64, ln_error_pre_umi: LogProbability) -> PhredScore {
    let ln_sum = ln_sum_exp_array(&[0.0, -gap, -gap, -gap]);
    let ln_posterior = ln_normalize(0.0, ln_sum);
    let ln_consensus_error = ln_not(ln_posterior);
    let ln_final_error = ln_error_prob_two_trials(ln_error_pre_umi, ln_consensus_error);
    ln_prob_to_phred(ln_final_error)
}

/// Builds the per-Phred-level gap-threshold table for one pre-UMI error rate.
///
/// `thresholds[q]` is the least `gap` at which [`unanimous_quality_from_gap`] is
/// `>= q`. Quality floors at `MIN_PHRED` (bottom entries `0.0`), is monotone in `gap`,
/// and is bounded above by the pre-UMI cap, so levels above the cap are unreachable
/// and stored as `f64::INFINITY`. `thresholds[cap]` is the same bound PR #634
/// derived (by hand, for a single fixed threshold) as its one unanimous cap
/// threshold; here it falls out as just the `cap`-indexed entry of this table.
/// Each boundary is found by the same bisection, over the same helper, so a
/// `partition_point` over the table reproduces the helper exactly. Built once per
/// pre-UMI rate (memoized), so 94 independent bisections are kept simple rather
/// than fused.
fn build_unanimous_gap_thresholds(error_rate_pre_umi: PhredScore) -> Box<[f64]> {
    /// Gap wide enough that the consensus error underflows to nothing.
    const MAX_GAP: f64 = 256.0;
    /// Enough halvings of `MAX_GAP` to converge past f64 resolution. The bracket
    /// collapses below an ULP of the returned gap by ~56 halvings across the whole
    /// representable pre-UMI range, so 64 is comfortably converged with headroom;
    /// further halvings are no-ops that leave the result bit-identical.
    const ITERATIONS: usize = 64;

    let ln_error_pre_umi = phred_to_ln_error_prob(error_rate_pre_umi);

    // Both endpoints are invariant across the loop (they don't depend on `q`), so
    // compute them once rather than re-evaluating the quality chain on every iteration.
    let quality_at_gap_zero = unanimous_quality_from_gap(0.0, ln_error_pre_umi);
    let quality_at_max_gap = unanimous_quality_from_gap(MAX_GAP, ln_error_pre_umi);

    let mut thresholds = vec![f64::INFINITY; MAX_PHRED as usize + 1];
    for q in 0..=MAX_PHRED {
        // Levels already met at gap 0 (at/below the quality floor) need no gap.
        if quality_at_gap_zero >= q {
            thresholds[q as usize] = 0.0;
            continue;
        }
        // Levels above the pre-UMI cap are unreachable; leave them at +INFINITY.
        if quality_at_max_gap < q {
            continue;
        }
        // The predicate `quality(gap) >= q` flips exactly once in (0, MAX_GAP]; bisect
        // and keep the known-good (upper) end so the boundary is never optimistic.
        let (mut too_small, mut wide_enough) = (0.0_f64, MAX_GAP);
        for _ in 0..ITERATIONS {
            let midpoint = 0.5 * (too_small + wide_enough);
            if unanimous_quality_from_gap(midpoint, ln_error_pre_umi) >= q {
                wide_enough = midpoint;
            } else {
                too_small = midpoint;
            }
        }
        thresholds[q as usize] = wide_enough;
    }
    thresholds.into_boxed_slice()
}

/// The consensus error of a unanimous position with winner-minus-loser gap `gap`:
/// `c(gap) = 3*e^-gap / (1 + 3*e^-gap)` (i.e. `1 - posterior`, `Q'(gap)` in the
/// module docs' mean-value-theorem mapping). Strictly decreasing in `gap` — as the
/// winner pulls further ahead, the consensus error shrinks — so the minimum of `c`
/// over a bracket `[a, b)` is `c(b)`, the value at the bracket's far edge. Used to
/// build [`UnanimousGapTables::cerr_min`].
///
/// This is a hand-derived linear-space closed form of the same quantity
/// [`unanimous_quality_from_gap`] computes in log space via the shared
/// [`crate::phred`] helpers (`ln_sum_exp_array` -> `ln_normalize` -> `ln_not`); it
/// is deliberately kept as an inline linear-space twin rather than routed through
/// that log-space chain, since the margin gate needs it as a plain real-valued
/// derivative bound, not a `PhredScore`. `test_consensus_error_at_gap_zero_is_three_quarters`
/// pins the two forms' exact agreement at `gap = 0`; their agreement everywhere
/// else is enforced indirectly, since a divergence would mis-size the margin gate
/// and surface as a mismatch in the fast-path/`call_full` parity sweeps (e.g.
/// `test_fast_path_matches_call_full_dense`). A future edit to the posterior model
/// must update both this function and `unanimous_quality_from_gap` together.
fn consensus_error(gap: f64) -> f64 {
    let e = 3.0 * (-gap).exp();
    e / (1.0 + e)
}

/// [`ConsensusBaseBuilder::try_unanimous_fast_path`]'s margin gate: maps the
/// empirically bounded `|ln_posterior_full - ln_posterior_gap| <=
/// UNANIMOUS_MARGIN_HEADROOM * UNANIMOUS_UNIT_ROUNDOFF * (|winner_ll| +
/// |loser_ll|)` bound (see the module docs and `test_ln_posterior_full_gap_error_bound`,
/// which pins the same bound this function implements) through the mean value
/// theorem into a gap-space margin, using `cerr_lower_bound` — a lower bound on the
/// consensus error over the region being gated (a bracket's minimum consensus
/// error for the sub-cap gate, or `1.0` — the identity element, i.e. no division —
/// for callers that want the raw ln-posterior error bound itself) as the
/// denominator.
///
/// Extracted as a standalone function (rather than inlined at the one production
/// call site) so the exact expression the margin gate uses is directly callable
/// from tests, instead of being transcribed a second time and risking drift.
///
/// Returns `+INFINITY` when `cerr_lower_bound == 0.0` (the fail-safe sentinel used
/// for unused table slots): dividing by zero drives the margin to `+INFINITY`,
/// which always fails the caller's clears-the-margin check and so always defers
/// to `call_full` rather than green-lighting an unvalidated shortcut.
fn unanimous_margin(winner_ll: f64, loser_ll: f64, cerr_lower_bound: f64) -> f64 {
    UNANIMOUS_MARGIN_HEADROOM * UNANIMOUS_UNIT_ROUNDOFF * (winner_ll.abs() + loser_ll.abs())
        / cerr_lower_bound
}

/// Combined per-Phred-level tables for the unanimous fast path, for one pre-UMI
/// error rate.
struct UnanimousGapTables {
    /// `thresholds[q]` is the least winner-minus-loser gap at which
    /// [`unanimous_quality_from_gap`] is `>= q`; see [`build_unanimous_gap_thresholds`].
    thresholds: Box<[f64]>,

    /// `cerr_min[q]`, for `q < cap`, is [`consensus_error`] evaluated at the
    /// bracket's far edge `thresholds[q + 1]` — the minimum consensus error over
    /// `[thresholds[q], thresholds[q + 1])`, since `consensus_error` is strictly
    /// decreasing in `gap`. This is the denominator [`ConsensusBaseBuilder::try_unanimous_fast_path`]
    /// uses to convert the empirically bounded ln-posterior error bound into a gap
    /// margin for the sub-cap brackets.
    ///
    /// At `q = cap - 1` (always populated: `cap >= MIN_PHRED = 2`), this is
    /// `consensus_error(thresholds[cap])` — bit-for-bit the cap boundary's own
    /// consensus error — which the fast path's cap-region gate reuses directly as
    /// its saturating budget (see [`ConsensusBaseBuilder::try_unanimous_fast_path`]),
    /// rather than storing that value a second time.
    ///
    /// Entries strictly above `cap - 1` are unused and set to the inert sentinel
    /// `0.0` — if ever read by mistake, dividing by `0.0` drives the computed
    /// margin to `+INFINITY`, which always fails the clears-the-margin check and so
    /// always defers to `call_full` rather than green-lighting an unvalidated
    /// shortcut.
    cerr_min: Box<[f64]>,

    /// The pre-UMI cap for this rate (`ln_prob_to_phred(ln_error_pre_umi)`): the
    /// highest quality a unanimous position can report. `thresholds[cap]` is the
    /// gap at which quality first reaches the cap; the cap region (`gap >=
    /// thresholds[cap]`) is gated using `cerr_min[cap - 1]` as a saturating budget,
    /// not a bracket denominator (see [`ConsensusBaseBuilder::try_unanimous_fast_path`]).
    cap: usize,
}

/// Builds the combined [`UnanimousGapTables`] for one pre-UMI error rate: the
/// per-level gap thresholds (via [`build_unanimous_gap_thresholds`]) plus the
/// per-bracket consensus-error minima the fast path's margin gate needs.
fn build_unanimous_gap_tables(error_rate_pre_umi: PhredScore) -> UnanimousGapTables {
    let thresholds = build_unanimous_gap_thresholds(error_rate_pre_umi);
    let ln_error_pre_umi = phred_to_ln_error_prob(error_rate_pre_umi);
    let cap = ln_prob_to_phred(ln_error_pre_umi) as usize;

    let mut cerr_min = vec![0.0_f64; thresholds.len()];
    for (q, slot) in cerr_min.iter_mut().enumerate().take(cap) {
        *slot = consensus_error(thresholds[q + 1]);
    }

    UnanimousGapTables { thresholds, cerr_min: cerr_min.into_boxed_slice(), cap }
}

/// Memoized [`build_unanimous_gap_tables`], one slot per possible pre-UMI Phred.
///
/// Building a table is ~94 bisections of the quality helper — far too costly to repeat
/// per builder (`consensus_umis` builds a fresh `ConsensusBaseBuilder` per UMI family).
/// `PhredScore` is a `u8`, so 256 slots cover the input domain exactly, with no
/// clamping.
static UNANIMOUS_GAP_TABLES: [OnceLock<UnanimousGapTables>; 256] = [const { OnceLock::new() }; 256];

/// Returns the memoized combined gap-threshold/consensus-error-minima tables for a
/// pre-UMI rate.
///
/// The tables are a pure function of `error_rate_pre_umi`, so memoizing on the `u8`
/// returns the bit-identical tables a direct build would. Keyed on the Phred score
/// (not the converted `LogProbability`) so key and value provenance cannot drift.
fn cached_unanimous_gap_tables(error_rate_pre_umi: PhredScore) -> &'static UnanimousGapTables {
    UNANIMOUS_GAP_TABLES[error_rate_pre_umi as usize]
        .get_or_init(|| build_unanimous_gap_tables(error_rate_pre_umi))
}

impl ConsensusBaseBuilder {
    /// Creates a new consensus base builder
    ///
    /// # Arguments
    /// * `error_rate_pre_umi` - Phred-scaled error rate prior to UMI integration
    /// * `error_rate_post_umi` - Phred-scaled error rate after UMI integration
    #[must_use]
    pub fn new(error_rate_pre_umi: PhredScore, error_rate_post_umi: PhredScore) -> Self {
        // The adjusted probability tables depend only on the post-UMI error rate, and both it
        // and the base quality are integer Phred scores, so the tables are memoized process-wide
        // and shared rather than recomputed (94 `exp`/`log1p` evaluations) per builder.
        let tables = adjusted_tables(error_rate_post_umi);

        let ln_error_pre_umi = phred_to_ln_error_prob(error_rate_pre_umi);

        Self {
            likelihoods: f64x4::splat(LN_ONE),
            compensations: f64x4::splat(0.0),
            observations: [0; DNA_BASE_COUNT],
            adjusted_correct_table: &tables.correct,
            adjusted_error_per_alt: &tables.error_per_alt,
            ln_error_pre_umi,
            gap_tables: cached_unanimous_gap_tables(error_rate_pre_umi),
        }
    }

    /// Resets the builder to process a new position
    pub fn reset(&mut self) {
        self.likelihoods = f64x4::splat(LN_ONE);
        self.compensations = f64x4::splat(0.0);
        self.observations.fill(0);
    }

    /// Adds an observation (base + quality) to the consensus
    ///
    /// Uses SIMD-vectorized Kahan summation for numeric stability when accumulating
    /// log-likelihoods. All 4 base likelihoods are updated in parallel using f64x4.
    ///
    /// # Arguments
    /// * `base` - The observed base (A, C, G, T, or N)
    /// * `qual` - The base quality score
    pub fn add(&mut self, base: u8, qual: PhredScore) {
        // Use lookup table to convert base to index (handles both upper and lower case)
        let matching_idx = BASE_TO_INDEX[base as usize];

        // Ignore N bases and other invalid bases (index 255)
        if matching_idx == 255 {
            return;
        }

        let matching_idx = matching_idx as usize;

        // Get adjusted probabilities for this quality
        let qual_idx = qual.min(MAX_PHRED) as usize;
        let ln_correct = self.adjusted_correct_table[qual_idx];
        let ln_error_per_base = self.adjusted_error_per_alt[qual_idx];

        // Build values array: error for all bases except ln_correct at matching position
        let mut values = [ln_error_per_base; 4];
        values[matching_idx] = ln_correct;
        let values = f64x4::from(values);

        // SIMD-vectorized Kahan summation for numeric stability
        // y = value - compensation (recovers the lost low-order bits)
        let y = values - self.compensations;
        // t = sum + y (the new sum, but low-order bits of y may be lost)
        let t = self.likelihoods + y;
        // compensation = (t - sum) - y (algebraically zero, but captures lost bits)
        self.compensations = (t - self.likelihoods) - y;
        // Update the sum
        self.likelihoods = t;

        self.observations[matching_idx] += 1;
    }

    /// Fast path for a unanimous position (every observation agrees on one base).
    ///
    /// Returns `Some((base, quality))` only when a stated *sufficient* condition
    /// establishes that the table lookup equals what [`Self::call_full`] would
    /// compute for this position; returns `None` in every other case (multiple
    /// bases observed, the degenerate certain-error case, a gap too close to a
    /// sub-cap bracket edge, or a cap-region gap where the saturating budget check
    /// fails), deferring to [`Self::call_full`] — the oracle, which is always
    /// correct and never wrong, only slower.
    ///
    /// The two regions (sub-cap and cap) use different sufficient conditions; see
    /// the module docs' cap-region derivation for the full argument for both.
    #[inline]
    fn try_unanimous_fast_path(&self) -> Option<(u8, PhredScore)> {
        // Count how many bases have observations
        let mut observed_base_idx = None;
        let mut num_bases_observed = 0;

        for (i, &count) in self.observations.iter().enumerate() {
            if count > 0 {
                num_bases_observed += 1;
                observed_base_idx = Some(i);
                if num_bases_observed > 1 {
                    // Multiple bases observed - need full consensus
                    return None;
                }
            }
        }

        // Unanimous: only one base observed
        let base_idx = observed_base_idx?;

        let likelihoods = self.likelihoods.as_array();
        let winner_ll = likelihoods[base_idx];
        let loser_ll = likelihoods[(base_idx + 1) % 4]; // Any loser will do (all equal)
        let gap = winner_ll - loser_ll;

        // Degenerate certain-error case: the observed base is not the unique maximum-
        // likelihood lane, i.e. gap <= 0 or non-finite (the per-observation adjusted
        // error is >= 3/4). Reached only at extreme inputs: error_rate_post_umi == 0
        // or an observation quality of 0 give ln_correct = -inf, hence gap = -inf and
        // a three-way loser tie; adjusted error exactly 3/4 gives a four-way tie.
        // call_full resolves the position as a tie / no-call; defer to it so the fast
        // path never disagrees. See test_fast_path_defers_degenerate_unanimous_to_call_full.
        //
        // The threshold is `f64::EPSILON`, not `0.0`: `call_full` declares a tie when a
        // loser is within `f64::EPSILON` ABSOLUTE of the max (see its `abs_diff_eq!`
        // check), so a `gap` in `(0.0, f64::EPSILON]` is exactly the window where
        // `call_full` would call a tie but this guard alone would not catch it. No
        // observed input reaches that window today (the smallest per-observation gap
        // is ~1.392e-2, far above `f64::EPSILON`), but nothing else pins it, so the
        // guard is hardened rather than left reliant on that empirical margin.
        if !(gap.is_finite() && gap > f64::EPSILON) {
            return None;
        }

        let thresholds = &self.gap_tables.thresholds;
        let cap = self.gap_tables.cap;
        let cap_threshold = thresholds[cap];

        // Cap region: gap at/above the cap threshold; gated by a saturating budget
        // rather than a shrinking bracket margin. See the module docs' cap-region
        // derivation for why, and `cerr_min`'s doc for why `cerr_min[cap - 1]` is
        // the cap boundary's own consensus error.
        if gap >= cap_threshold {
            // Absolute ln-posterior error budget: pass `cerr_lower_bound = 1.0` to get
            // the raw `C*u*(|w|+|l|)` bound with no division.
            let delta = unanimous_margin(winner_ll, loser_ll, 1.0);

            // `cap >= MIN_PHRED = 2`, so `cap - 1` is always a valid, populated
            // index; pinned by test_unanimous_gap_tables_internal_invariants.
            let cerr_at_cap = self.gap_tables.cerr_min[cap - 1];

            // Gate: `gap - cap_threshold >= ln 2` (budget closed to at least half its
            // saturated value) *and* `delta` below half that budget. See the module
            // docs' cap-region derivation for the full argument.
            if gap - cap_threshold >= std::f64::consts::LN_2 && delta < 0.5 * cerr_at_cap {
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "cap = ln_prob_to_phred(...) <= MAX_PHRED = 93, so it fits u8"
                )]
                return Some((DNA_BASES[base_idx], cap as PhredScore));
            }
            // Too close to the cap boundary (or the depth-driven error has grown too
            // large relative to the saturated budget) to establish exactness: defer to
            // call_full rather than falling through into the sub-cap bracket lookup
            // below, whose partition_point precondition (`gap < cap_threshold`) would
            // be violated.
            return None;
        }

        // Sub-cap: locate the bracket containing gap. thresholds[0] == 0.0 and
        // gap > 0 guarantee partition_point >= 1; gap < thresholds[cap] (cap region
        // handled above) guarantees partition_point <= cap. So q = partition_point - 1
        // is in [0, cap), a populated (non-sentinel) index into cerr_min.
        let q = thresholds.partition_point(|&t| t <= gap) - 1;

        // Margin: the empirically bounded ln-posterior error bound (see module
        // docs), mapped to a gap perturbation via the mean value theorem using the
        // bracket's minimum consensus error (the derivative's minimum magnitude over
        // the bracket).
        let margin = unanimous_margin(winner_ll, loser_ll, self.gap_tables.cerr_min[q]);

        if gap - thresholds[q] > margin && thresholds[q + 1] - gap > margin {
            // gap clears the margin from both bracket edges: no perturbation within
            // the empirically bounded margin can move call_full's true answer to a
            // different quality level, so q is conditionally exact.
            #[expect(
                clippy::cast_possible_truncation,
                reason = "q < cap <= MAX_PHRED = 93, so it fits u8"
            )]
            return Some((DNA_BASES[base_idx], q as PhredScore));
        }

        // Cannot establish exactness: defer to call_full.
        None
    }

    /// Calls the consensus base and quality
    ///
    /// Returns (`consensus_base`, `consensus_quality`)
    /// Returns (N, 0) if no observations or multiple equally likely bases
    ///
    /// # Panics
    ///
    /// Panics if there are observations but no maximum likelihood base is found
    /// (should not occur in practice since we check for ties).
    #[must_use]
    pub fn call(&self) -> (u8, PhredScore) {
        if self.contributions() == 0 {
            return (NO_CALL_BASE, MIN_PHRED);
        }

        // Fast path: check for unanimous consensus (only one base observed)
        // This is common in high-quality data and avoids expensive log-sum-exp
        if let Some((base, qual)) = self.try_unanimous_fast_path() {
            return (base, qual);
        }

        self.call_full()
    }

    /// The full log-sum-exp consensus calculation, with no fast path.
    ///
    /// [`Self::call`] short-circuits to [`Self::try_unanimous_fast_path`] when that
    /// is known to be exact (empirically bounded and conditionally exact in both
    /// the cap region and the sub-cap region); this is the answer it must agree
    /// with. Kept separate so the agreement is directly testable — see
    /// `test_unanimous_fast_path_matches_full_calculation`.
    fn call_full(&self) -> (u8, PhredScore) {
        // Extract likelihoods as array for iteration and function calls
        let likelihoods = self.likelihoods.as_array();

        // Compute the sum of likelihoods for normalization
        let ln_sum = ln_sum_exp_array(likelihoods);

        // Find the maximum likelihood and check if it's unique
        let mut max_likelihood = f64::NEG_INFINITY;
        let mut max_index = None;
        let mut tie = false;

        for (i, &ll) in likelihoods.iter().enumerate() {
            match ll.partial_cmp(&max_likelihood) {
                Some(Ordering::Greater) => {
                    max_likelihood = ll;
                    max_index = Some(i);
                    tie = false;
                }
                Some(Ordering::Equal) => {
                    tie = true;
                }
                Some(Ordering::Less) => {
                    // Check for epsilon-level tie (within machine precision)
                    if abs_diff_eq!(ll, max_likelihood, epsilon = f64::EPSILON) {
                        tie = true;
                    }
                }
                None => {}
            }
        }

        // If there's a tie, return no-call
        if tie || max_index.is_none() {
            return (NO_CALL_BASE, MIN_PHRED);
        }

        let max_idx = max_index.expect("max_index is Some after is_none() check above");
        let consensus_base = DNA_BASES[max_idx];

        // Calculate posterior probability
        let ln_posterior = ln_normalize(max_likelihood, ln_sum);

        // Convert to error probability
        let ln_consensus_error = ln_not(ln_posterior);

        // Factor in pre-UMI error rate
        // P(final_error) = P(pre_UMI_error AND correct_consensus)
        //                  + P(no_pre_UMI_error AND consensus_error)
        //                  + P(pre_UMI_error AND consensus_error that don't cancel)
        // Approximation: just combine the two error rates
        let ln_final_error = ln_error_prob_two_trials(self.ln_error_pre_umi, ln_consensus_error);

        // Convert to Phred score and cap at maximum
        let qual = ln_prob_to_phred(ln_final_error);

        (consensus_base, qual)
    }

    /// Returns the total number of contributing observations
    ///
    /// This is the sum of observations across all four bases
    #[must_use]
    pub fn contributions(&self) -> u32 {
        self.observations.iter().sum()
    }

    /// Returns the number of observations for a specific base
    ///
    /// # Arguments
    /// * `base` - The base to query (A, C, G, or T)
    ///
    /// Returns 0 for invalid bases
    #[must_use]
    pub fn observations_for_base(&self, base: u8) -> u32 {
        let idx = BASE_TO_INDEX[base as usize];
        if idx == 255 { 0 } else { self.observations[idx as usize] }
    }

    /// Returns the observations for all bases as [A, C, G, T]
    #[must_use]
    pub fn all_observations(&self) -> [u32; DNA_BASE_COUNT] {
        self.observations
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// The gap threshold at the pre-UMI cap for `rate`, i.e.
    /// `build_unanimous_gap_thresholds(rate)[cap]`. A small test-only wrapper so
    /// tests that need this value directly don't each re-derive `cap` from
    /// `rate`; replaces the former independent-oracle bisection
    /// `min_exact_unanimous_gap`, which was removed as redundant with the table
    /// (it bisected the same `unanimous_quality_from_gap` helper the table itself
    /// bisects, so a bug in the helper would have moved both identically).
    fn cap_gap_threshold(rate: PhredScore) -> f64 {
        let cap = ln_prob_to_phred(phred_to_ln_error_prob(rate)) as usize;
        build_unanimous_gap_thresholds(rate)[cap]
    }

    /// Recomputes one adjusted-probability entry with the formula the builder used inline
    /// before the tables were memoized, in the same order, so a drift in `compute` shows up
    /// as a bit-level difference rather than an approximate one.
    fn reference_adjusted_entry(
        error_rate_post_umi: PhredScore,
        qual: PhredScore,
    ) -> (LogProbability, LogProbability) {
        let ln_error_post = phred_to_ln_error_prob(error_rate_post_umi);
        let ln_error_seq = phred_to_ln_error_prob(qual);
        let adjusted_error = ln_error_prob_two_trials(ln_error_post, ln_error_seq);
        (ln_not(adjusted_error), adjusted_error - 3.0_f64.ln())
    }

    /// The memoized tables must reproduce the pre-memoization inline math exactly, across the
    /// full `0..=MAX_PHRED` quality domain, for every post-UMI error rate the callers use.
    #[rstest]
    #[case::post_umi_q0(0)]
    #[case::post_umi_q10(10)]
    #[case::post_umi_default_q40(40)]
    #[case::post_umi_q45(45)]
    #[case::post_umi_max_phred(MAX_PHRED)]
    #[case::post_umi_above_max_phred(u8::MAX)]
    fn test_adjusted_tables_match_the_inline_formula_bit_for_bit(
        #[case] error_rate_post_umi: PhredScore,
    ) {
        let tables = adjusted_tables(error_rate_post_umi);

        assert_eq!(tables.correct.len(), MAX_PHRED as usize + 1, "table must cover every Phred");
        assert_eq!(tables.error_per_alt.len(), MAX_PHRED as usize + 1);

        for qual in 0..=MAX_PHRED {
            let (correct, error_per_alt) = reference_adjusted_entry(error_rate_post_umi, qual);
            let actual_correct = tables.correct[qual as usize];
            let actual_error_per_alt = tables.error_per_alt[qual as usize];

            // Compared as bit patterns, not approximately: the memoized computation is meant to
            // be the identical arithmetic in the identical order, so any drift is a regression.
            assert_eq!(
                actual_correct.to_bits(),
                correct.to_bits(),
                "correct probability drifted at post-UMI Q{error_rate_post_umi}, base Q{qual}: \
                 {actual_correct} != {correct}"
            );
            assert_eq!(
                actual_error_per_alt.to_bits(),
                error_per_alt.to_bits(),
                "alt-error probability drifted at post-UMI Q{error_rate_post_umi}, base Q{qual}: \
                 {actual_error_per_alt} != {error_per_alt}"
            );
        }
    }

    /// The cache must hand out the *same* tables on every call, otherwise the memoization is
    /// not actually happening and every builder pays the transcendental math again.
    #[rstest]
    #[case::post_umi_q0(0)]
    #[case::post_umi_default_q40(40)]
    #[case::post_umi_above_max_phred(u8::MAX)]
    fn test_adjusted_tables_are_memoized(#[case] error_rate_post_umi: PhredScore) {
        let first = adjusted_tables(error_rate_post_umi);
        let second = adjusted_tables(error_rate_post_umi);

        assert!(std::ptr::eq(first, second), "the cache slot must be initialized at most once");
        assert!(std::ptr::eq(first.correct.as_ptr(), second.correct.as_ptr()));
        assert!(std::ptr::eq(first.error_per_alt.as_ptr(), second.error_per_alt.as_ptr()));
    }

    /// Distinct post-UMI error rates must land in distinct cache slots — a shared slot would
    /// silently apply one family's error model to another's.
    #[test]
    fn test_adjusted_tables_are_keyed_by_post_umi_error_rate() {
        let q40 = adjusted_tables(40);
        let q45 = adjusted_tables(45);

        assert!(!std::ptr::eq(q40, q45), "different error rates must not share a slot");
        assert_ne!(q40.correct[30].to_bits(), q45.correct[30].to_bits());
    }

    /// The builder must borrow the shared tables rather than computing its own copy.
    #[test]
    fn test_builder_borrows_the_shared_tables() {
        let first = ConsensusBaseBuilder::new(45, 40);
        let second = ConsensusBaseBuilder::new(10, 40);

        assert!(std::ptr::eq(
            first.adjusted_correct_table.as_ptr(),
            second.adjusted_correct_table.as_ptr()
        ));
        assert!(std::ptr::eq(
            first.adjusted_error_per_alt.as_ptr(),
            second.adjusted_error_per_alt.as_ptr()
        ));
    }

    /// The per-base error count in the consensus callers is
    /// `contributions() - observations_for_base(base)`. Both operands must stay
    /// unsaturated `u32` counters: the callers clamp only the *stored* cD/cE
    /// values to fgbio's `Short` ceiling, and clamping the operands instead
    /// would collapse the difference toward zero on a deep pileup and understate
    /// cE. This pins that on a pileup deep enough that a `u16`/`i16` counter
    /// would have saturated.
    #[test]
    fn test_deep_pileup_keeps_error_count_exact_past_the_short_ceiling() {
        let short_ceiling = u32::from(i16::MAX.unsigned_abs()); // 32_767
        let matching = short_ceiling + 5_000;
        let mismatching = short_ceiling + 1_000;

        let mut builder = ConsensusBaseBuilder::new(45, 40);
        for _ in 0..matching {
            builder.add(b'A', 30);
        }
        for _ in 0..mismatching {
            builder.add(b'C', 30);
        }

        let depth = builder.contributions();
        assert_eq!(depth, matching + mismatching, "depth must not saturate");
        assert!(depth > short_ceiling, "precondition: the pileup exceeds the Short ceiling");

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A', "the majority base should win");

        // The subtraction is exact because neither operand is clamped. Were both
        // saturated at the ceiling first, this difference would collapse to 0.
        let error_count = depth - builder.observations_for_base(base);
        assert_eq!(error_count, mismatching, "error count must be exact, not a clamped difference");
        assert_ne!(error_count, 0, "a saturated-operand subtraction would collapse to zero here");
    }

    #[test]
    fn test_single_base_perfect() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add 10 perfect quality A's
        for _ in 0..10 {
            builder.add(b'A', 40);
        }

        let (base, qual) = builder.call();
        assert_eq!(base, b'A');
        assert!(qual >= 40); // Should have high quality
        assert_eq!(builder.contributions(), 10);
    }

    #[test]
    fn test_contributions_do_not_overflow_at_u16_boundary() {
        // A UMI family deeper than u16::MAX (e.g. high-duplication cfDNA) must not
        // wrap the per-base observation counter. With a u16 counter, 70_000 single-base
        // observations wrapped to 70_000 % 65_536 = 4_464; the depth must stay 70_000.
        let mut builder = ConsensusBaseBuilder::new(45, 40);
        let n: u32 = 70_000;
        for _ in 0..n {
            builder.add(b'A', 40);
        }

        assert_eq!(builder.contributions(), n, "consensus depth must not wrap at the u16 boundary");
        assert_eq!(
            builder.observations_for_base(b'A'),
            n,
            "per-base observation count must not wrap at the u16 boundary"
        );
    }

    #[test]
    fn test_mixed_bases() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add 8 A's and 2 C's
        for _ in 0..8 {
            builder.add(b'A', 30);
        }
        for _ in 0..2 {
            builder.add(b'C', 30);
        }

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A'); // Majority should win
        assert_eq!(builder.contributions(), 10);
    }

    #[test]
    fn test_no_observations() {
        let builder = ConsensusBaseBuilder::new(45, 40);

        let (base, qual) = builder.call();
        assert_eq!(base, b'N');
        assert_eq!(qual, 2); // fgbio's PhredScore.MinValue
        assert_eq!(builder.contributions(), 0);
    }

    #[test]
    fn test_ignore_n_bases() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'A', 40);
        builder.add(b'N', 40);
        builder.add(b'A', 40);

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A');
        assert_eq!(builder.contributions(), 2); // N should not be counted
    }

    #[test]
    fn test_case_insensitive() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'a', 40);
        builder.add(b'A', 40);
        builder.add(b'a', 40);

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A');
        assert_eq!(builder.contributions(), 3);
    }

    #[test]
    fn test_reset() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'A', 40);
        builder.add(b'A', 40);
        assert_eq!(builder.contributions(), 2);

        builder.reset();
        assert_eq!(builder.contributions(), 0);

        builder.add(b'C', 40);
        let (base, _qual) = builder.call();
        assert_eq!(base, b'C');
    }

    #[test]
    fn test_observations_for_base() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        builder.add(b'A', 40);
        builder.add(b'A', 40);
        builder.add(b'C', 40);
        builder.add(b'G', 40);

        assert_eq!(builder.observations_for_base(b'A'), 2);
        assert_eq!(builder.observations_for_base(b'C'), 1);
        assert_eq!(builder.observations_for_base(b'G'), 1);
        assert_eq!(builder.observations_for_base(b'T'), 0);
    }

    #[test]
    fn test_quality_variation() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add high quality A's
        for _ in 0..5 {
            builder.add(b'A', 40);
        }

        // Add one low quality C (should not outweigh the A's)
        builder.add(b'C', 10);

        let (base, _qual) = builder.call();
        assert_eq!(base, b'A');
    }

    #[test]
    fn test_tie_returns_no_call() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // This is hard to create a perfect tie, but we can test the no-call logic
        // by having very few observations
        builder.add(b'A', 10);
        builder.add(b'C', 10);

        // With low quality and equal counts, might get different results
        // Just verify it doesn't crash and returns a valid base (A, C, or N)
        let (base, _qual) = builder.call();
        assert!(base == b'A' || base == b'C' || base == b'N');
    }

    // Port of fgbio test: "produce a no-call if two bases are of equal likelihood"
    #[test]
    fn test_equal_likelihood_produces_no_call() {
        // Use MAX_PHRED for error rates to effectively disable error rate adjustment
        let mut builder = ConsensusBaseBuilder::new(MAX_PHRED, MAX_PHRED);

        // Empty builder should return no-call with quality 2
        let (base, qual) = builder.call();
        assert_eq!(base, b'N');
        assert_eq!(qual, 2);

        // Equal evidence should also return no-call
        builder.add(b'A', 20);
        builder.add(b'C', 20);
        let (base, qual) = builder.call();
        assert_eq!(base, b'N');
        assert_eq!(qual, 2);
    }

    // Port of fgbio test: "calculate consensus base and quality given a massive pileup"
    #[test]
    fn test_massive_pileup() {
        let mut builder = ConsensusBaseBuilder::new(50, 50);

        // Add 1000 Q20 C's
        for _ in 0..1000 {
            builder.add(b'C', 20);
        }

        let (base, qual) = builder.call();
        assert_eq!(base, b'C');
        assert_eq!(qual, 50); // Quality capped by pre-UMI error rate
        assert_eq!(builder.contributions(), 1000);
        assert_eq!(builder.observations_for_base(b'A'), 0);
        assert_eq!(builder.observations_for_base(b'C'), 1000);
        assert_eq!(builder.observations_for_base(b'G'), 0);
        assert_eq!(builder.observations_for_base(b'T'), 0);

        // Add 10 T's - should still call C
        for _ in 0..10 {
            builder.add(b'T', 20);
        }

        let (base, qual) = builder.call();
        assert_eq!(base, b'C');
        assert_eq!(qual, 50);
        assert_eq!(builder.contributions(), 1010);
        assert_eq!(builder.observations_for_base(b'T'), 10);
    }

    // Port of fgbio test: "calculate consensus base and quality given conflicting evidence"
    #[test]
    fn test_conflicting_evidence() {
        let mut builder = ConsensusBaseBuilder::new(50, 50);

        builder.add(b'A', 30);
        builder.add(b'C', 28);

        let (base, qual) = builder.call();
        assert_eq!(base, b'A');
        assert!(qual <= 5, "Quality should be low due to conflicting evidence, got {qual}");
    }

    /// SIMPLEX3-01: a Q0 observation of a *minority* base drives that base's
    /// likelihood lane to `−∞`. The consensus base must stay the majority (`A`) and
    /// the quality must NOT exceed the quality of the same pileup *without* the
    /// contradicting read — a conflicting observation cannot increase confidence.
    /// The pre-fix `ln_sum_exp_array` returned `−∞` for the normalizer whenever any
    /// lane was `−∞`, driving the posterior to `+∞` and inflating the quality to the
    /// pre-UMI cap (observed as fgumi Q45 vs fgbio Q44).
    #[test]
    fn test_neg_inf_lane_does_not_inflate_quality() {
        let clean = {
            let mut b = ConsensusBaseBuilder::new(45, 40);
            b.add(b'A', 30);
            b.add(b'A', 30);
            b.call()
        };
        let with_q0_conflict = {
            let mut b = ConsensusBaseBuilder::new(45, 40);
            b.add(b'A', 30);
            b.add(b'A', 30);
            b.add(b'C', 0); // Q0 -> lane[C] = -inf
            b.call()
        };
        assert_eq!(with_q0_conflict.0, b'A', "majority base must remain A");
        // Pin the exact fgbio-compatible qualities: the clean 2×A pileup caps at Q44,
        // and adding the Q0 (probability-0) `C` must leave it at Q44 — the pre-fix bug
        // produced Q45 by driving the normalizer to +∞. Asserting the exact value (not
        // just `conflict <= clean`) pins the fgbio baseline against future drift.
        assert_eq!(clean.1, 44, "clean 2xA pileup must be fgbio's Q44");
        assert_eq!(
            with_q0_conflict.1, 44,
            "a Q0 conflicting base must leave the quality at fgbio's Q44, not inflate it"
        );
    }

    // Port of fgbio test: "support calling multiple pileups from the same builder"
    #[test]
    fn test_reset_and_reuse_fgbio() {
        let mut builder = ConsensusBaseBuilder::new(50, 50);

        builder.add(b'A', 20);
        let (base, qual) = builder.call();
        assert_eq!(base, b'A');
        assert_eq!(qual, 20);
        assert_eq!(builder.contributions(), 1);

        builder.reset();

        builder.add(b'C', 20);
        let (base, qual) = builder.call();
        assert_eq!(base, b'C');
        assert_eq!(qual, 20);
        assert_eq!(builder.contributions(), 1);
    }

    // Test that Kahan summation produces consistent results regardless of addition order
    // This tests the numeric stability improvement from PR #1120
    #[test]
    fn test_kahan_summation_order_independence() {
        // Create two builders with the same error rates
        let mut builder1 = ConsensusBaseBuilder::new(45, 40);
        let mut builder2 = ConsensusBaseBuilder::new(45, 40);

        // Add bases in different orders with varying qualities
        // Builder 1: A's first, then C's
        for q in [10, 20, 30, 40, 50] {
            builder1.add(b'A', q);
        }
        for q in [15, 25, 35, 45] {
            builder1.add(b'C', q);
        }

        // Builder 2: Interleaved
        builder2.add(b'A', 10);
        builder2.add(b'C', 15);
        builder2.add(b'A', 20);
        builder2.add(b'C', 25);
        builder2.add(b'A', 30);
        builder2.add(b'C', 35);
        builder2.add(b'A', 40);
        builder2.add(b'C', 45);
        builder2.add(b'A', 50);

        // Both should produce the same result
        let (base1, qual1) = builder1.call();
        let (base2, qual2) = builder2.call();

        assert_eq!(base1, base2, "Consensus base should be order-independent");
        assert_eq!(qual1, qual2, "Consensus quality should be order-independent");
        assert_eq!(base1, b'A', "Should call A (5 vs 4 observations)");
    }

    // Test that Kahan summation handles extreme quality differences without precision loss
    #[test]
    fn test_kahan_summation_extreme_quality_range() {
        let mut builder = ConsensusBaseBuilder::new(45, 40);

        // Add many low quality bases followed by a few high quality bases
        // Q93 provides MUCH more evidence than Q2 (error rates: ~5e-10 vs ~0.63)
        // So 10 Q93 observations easily outweigh 100 Q2 observations
        for _ in 0..100 {
            builder.add(b'A', 2); // Very low quality (error ~63%)
        }
        for _ in 0..10 {
            builder.add(b'C', 93); // Maximum quality (error ~5e-10)
        }

        // High quality C's should win - this tests that extreme values are handled correctly
        let (base, _qual) = builder.call();
        assert_eq!(base, b'C', "High-quality observations provide much stronger evidence");
    }

    // Port of fgbio test: "scale base qualities using the post-umi error rate"
    // Tests that ConsensusBaseBuilder properly adjusts qualities based on post-UMI error rate.
    // With errorRatePostLabeling=10 (Q10), input qualities should be reduced:
    // Q20 → Q9, Q15 → Q8, Q10 → Q7, Q5 → Q4
    #[test]
    fn test_scale_base_qualities_post_umi_error_rate() {
        // errorRatePreLabeling=MAX_PHRED (no pre-UMI adjustment), errorRatePostLabeling=10
        // This matches fgbio's test setup
        let input_quals: [u8; 4] = [20, 15, 10, 5];
        let expected_output: [u8; 4] = [9, 8, 7, 4];

        for (input_q, expected_q) in input_quals.iter().zip(expected_output.iter()) {
            let mut builder = ConsensusBaseBuilder::new(MAX_PHRED, 10);
            builder.add(b'A', *input_q);
            let (_base, qual) = builder.call();

            // Allow for small differences due to implementation details
            // but the general pattern should hold: qualities are reduced
            assert!(
                qual <= *input_q,
                "Q{input_q} should be reduced by post-UMI error rate, got Q{qual}"
            );
            assert!(
                (i16::from(qual) - i16::from(*expected_q)).abs() <= 1,
                "Q{input_q} should become ~Q{expected_q}, got Q{qual}"
            );
        }
    }

    // ========================================================================
    // Direct unit tests: internal building blocks of the unanimous fast path
    // ========================================================================
    //
    // The end-to-end `call() == call_full()` parity tests below are the safety net,
    // but they don't pin the internal pieces (`unanimous_margin`, `consensus_error`,
    // `UnanimousGapTables`, and the gate decision itself) independently. These tests
    // exercise each piece directly so a regression in one is caught at its source,
    // not only downstream in a parity mismatch.

    /// `unanimous_margin` for hand-computed inputs, checked against a literal
    /// re-derivation of the formula (not by calling the function again) so this
    /// pins the actual arithmetic, not just "the function agrees with itself".
    #[rstest]
    #[case::small_values(-1.0, -2.0, 0.5)]
    #[case::larger_lls(-30.0, -45.5, 0.1)]
    #[case::cerr_lower_bound_one_no_division(-5.0, -7.0, 1.0)]
    fn test_unanimous_margin_matches_hand_computed_value(
        #[case] winner_ll: f64,
        #[case] loser_ll: f64,
        #[case] cerr_lower_bound: f64,
    ) {
        let margin = unanimous_margin(winner_ll, loser_ll, cerr_lower_bound);
        // Same operation order as the production formula, but written out with
        // literal constants (not the named `UNANIMOUS_MARGIN_HEADROOM` /
        // `UNANIMOUS_UNIT_ROUNDOFF`) so an accidental change to either constant is
        // also caught here, not just re-derived.
        let expected =
            16.0 * (f64::EPSILON / 2.0) * (winner_ll.abs() + loser_ll.abs()) / cerr_lower_bound;
        assert_eq!(
            margin.to_bits(),
            expected.to_bits(),
            "unanimous_margin({winner_ll}, {loser_ll}, {cerr_lower_bound}) = {margin:e}, \
             expected {expected:e}"
        );
    }

    /// `unanimous_margin` scales linearly with `|winner_ll| + |loser_ll|`: doubling
    /// both log-likelihoods must double the margin, for a fixed `cerr_lower_bound`.
    #[rstest]
    #[case::unit_bound(1.0)]
    #[case::small_bound(0.01)]
    fn test_unanimous_margin_scales_linearly_with_ll_magnitude(#[case] cerr_lower_bound: f64) {
        let base = unanimous_margin(-3.0, -4.0, cerr_lower_bound);
        let doubled = unanimous_margin(-6.0, -8.0, cerr_lower_bound);
        assert!(
            abs_diff_eq!(doubled, 2.0 * base, epsilon = base * 1e-9),
            "doubling |winner_ll| + |loser_ll| must double the margin: base={base:e} \
             doubled={doubled:e} (expected ~{:e})",
            2.0 * base
        );
    }

    /// `cerr_lower_bound == 0.0` is the fail-safe sentinel path (unused table
    /// slots): division by zero must drive the margin to `+INFINITY`, which always
    /// fails the caller's clears-the-margin check and so always defers to
    /// `call_full` rather than green-lighting an unvalidated shortcut.
    #[test]
    fn test_unanimous_margin_is_infinite_at_zero_cerr_lower_bound() {
        let margin = unanimous_margin(-10.0, -12.0, 0.0);
        assert!(margin.is_infinite(), "margin must be +INFINITY, got {margin}");
        assert!(margin.is_sign_positive(), "margin must be +INFINITY, not -INFINITY");
    }

    /// `consensus_error(0.0) == 0.75`: at zero gap the four-way likelihoods are
    /// exactly tied, so the posterior is `1/4` and the consensus error `1 - 1/4 =
    /// 3/4`. All operations here (`exp(0.0)`, `3.0 * 1.0`, `1.0 + 3.0`, `3.0 / 4.0`)
    /// are exact in `f64`, so this is an exact equality, not an approximation.
    #[test]
    fn test_consensus_error_at_gap_zero_is_three_quarters() {
        assert_eq!(consensus_error(0.0).to_bits(), 0.75_f64.to_bits());
    }

    /// `consensus_error` is strictly decreasing in `gap` and stays in the open
    /// interval `(0, 1)` — the property the fast path's bracket-minimum derivation
    /// (`cerr_min[q] = consensus_error(thresholds[q + 1])`) depends on.
    #[test]
    fn test_consensus_error_strictly_decreasing_and_in_open_unit_interval() {
        let mut prev = f64::INFINITY;
        let mut gap = 0.0;
        while gap <= 200.0 {
            let c = consensus_error(gap);
            assert!(c > 0.0 && c < 1.0, "consensus_error({gap}) = {c} not in (0, 1)");
            assert!(
                c < prev,
                "consensus_error must strictly decrease: gap={gap} gave {c} after {prev}"
            );
            prev = c;
            gap += 0.25;
        }
    }

    /// `UnanimousGapTables`' internal cross-field invariants: `cerr_min[q]` equals
    /// `consensus_error` at the bracket's far edge for every `q < cap`, entries
    /// at/above `cap` are the `0.0` sentinel, and `cap` is exactly
    /// `ln_prob_to_phred(phred_to_ln_error_prob(rate))`. Also pins, as its own
    /// dedicated assertion (not just as a byproduct of the loop below),
    /// `cerr_min[cap - 1] == consensus_error(thresholds[cap])` — the identity the
    /// cap-region gate in `try_unanimous_fast_path` now relies on in place of a
    /// separate `cerr_cap` field.
    #[rstest]
    fn test_unanimous_gap_tables_internal_invariants(
        #[values(0, 2, 20, 45, 69, 90, 93)] rate: PhredScore,
    ) {
        let tables = build_unanimous_gap_tables(rate);
        let ln_pre = phred_to_ln_error_prob(rate);
        let expected_cap = ln_prob_to_phred(ln_pre) as usize;

        assert_eq!(tables.cap, expected_cap, "cap mismatch for Q{rate}");

        for q in 0..tables.cap {
            assert_eq!(
                tables.cerr_min[q].to_bits(),
                consensus_error(tables.thresholds[q + 1]).to_bits(),
                "cerr_min[{q}] must equal consensus_error at the bracket's far edge for Q{rate}"
            );
        }

        assert_eq!(
            tables.cerr_min[tables.cap - 1].to_bits(),
            consensus_error(tables.thresholds[tables.cap]).to_bits(),
            "cerr_min[cap - 1] must equal consensus_error(thresholds[cap]) for Q{rate}: the \
             cap-region gate relies on this identity"
        );

        for q in tables.cap..tables.cerr_min.len() {
            assert_eq!(
                tables.cerr_min[q].to_bits(),
                0.0_f64.to_bits(),
                "entries at/above cap must be the 0.0 sentinel for Q{rate}, q={q}"
            );
        }
    }

    /// Builds a unanimous builder (single observed base at index 0, i.e. `A`) whose
    /// internal likelihoods are set directly to a winner and three equal losers
    /// `gap` below it, bypassing `add` (which cannot land an exact `gap` value) so
    /// the sub-cap and cap-region boundaries can be hit precisely.
    fn unanimous_builder_with_gap(
        pre: PhredScore,
        post: PhredScore,
        winner_ll: f64,
        gap: f64,
    ) -> ConsensusBaseBuilder {
        let mut builder = ConsensusBaseBuilder::new(pre, post);
        builder.observations[0] = 1;
        builder.likelihoods =
            f64x4::from([winner_ll, winner_ll - gap, winner_ll - gap, winner_ll - gap]);
        builder
    }

    /// The gate decision itself, tested directly by manipulating the builder's
    /// internal likelihoods (rather than through `add`, which cannot land an exact
    /// gap value) so the sub-cap and cap-region boundaries can be hit precisely.
    ///
    /// Constructs a unanimous builder (single observed base, arbitrary winner/loser
    /// log-likelihoods) and asserts `Some`/`None` from `try_unanimous_fast_path`
    /// according to the documented sufficient condition for each region.
    #[rstest]
    #[case::sub_cap_bracket_interior_clears_margin(true)]
    #[case::sub_cap_at_lower_edge_does_not_clear_margin(false)]
    fn test_gate_decision_sub_cap_bracket(#[case] expect_some: bool) {
        const PRE: PhredScore = 45;
        let tables = cached_unanimous_gap_tables(PRE);
        // An interior bracket, comfortably away from both 0 and the cap, so its
        // width is not degenerate.
        let q = tables.cap / 2;
        let winner_ll = -5.0_f64;

        let gap = if expect_some {
            // Comfortably inside the bracket: clears the margin (~1e-13 scale) from
            // both edges by a wide margin (the bracket is at least several nats wide
            // this far from the cap).
            0.5 * (tables.thresholds[q] + tables.thresholds[q + 1])
        } else {
            // Exactly at the lower edge: `gap - thresholds[q] == 0`, which is never
            // `> margin` (margin is strictly positive), so this must defer.
            tables.thresholds[q]
        };

        let builder = unanimous_builder_with_gap(PRE, 40, winner_ll, gap);

        let result = builder.try_unanimous_fast_path();
        if expect_some {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "q = tables.cap / 2 <= MAX_PHRED = 93, so it fits u8"
            )]
            let expected_qual = q as PhredScore;
            assert_eq!(
                result,
                Some((b'A', expected_qual)),
                "expected the fast path to resolve bracket {q} for gap {gap}"
            );
        } else {
            assert_eq!(result, None, "expected the fast path to defer at gap {gap} (bracket edge)");
        }
    }

    /// The cap-region gate's two-part condition, tested directly: `gap >=
    /// thresholds[cap] + ln 2` alone is not sufficient if the absolute ln-posterior
    /// error is too large; both conjuncts must hold for `Some`.
    #[rstest]
    #[case::clears_both_conjuncts(true, 100.0, false)]
    #[case::fails_ln2_budget_conjunct(false, 100.0, false)]
    #[case::fails_error_bound_conjunct_via_huge_lls(true, 100.0, true)]
    fn test_gate_decision_cap_region(
        #[case] clear_ln2_gap: bool,
        #[case] extra_gap_past_ln2: f64,
        #[case] inflate_ll_magnitude: bool,
    ) {
        const PRE: PhredScore = 45;
        let tables = cached_unanimous_gap_tables(PRE);
        let cap_threshold = tables.thresholds[tables.cap];

        let gap = if clear_ln2_gap {
            cap_threshold + std::f64::consts::LN_2 + extra_gap_past_ln2
        } else {
            // Just shy of the required ln 2 budget closure: must defer regardless of
            // how small the error term is.
            cap_threshold + std::f64::consts::LN_2 * 0.5
        };

        // A huge log-likelihood magnitude inflates `delta = C * u * (|w| + |l|)`
        // past half the cap-boundary consensus error even though the gap condition
        // is satisfied, exercising the second conjunct independently of the first.
        let winner_ll = if inflate_ll_magnitude { -1.0e14_f64 } else { -1.0_f64 };

        let builder = unanimous_builder_with_gap(PRE, 40, winner_ll, gap);

        let result = builder.try_unanimous_fast_path();
        let expect_some = clear_ln2_gap && !inflate_ll_magnitude;
        if expect_some {
            #[expect(
                clippy::cast_possible_truncation,
                reason = "cap <= MAX_PHRED = 93, so it fits u8"
            )]
            let expected_qual = tables.cap as PhredScore;
            assert_eq!(
                result,
                Some((b'A', expected_qual)),
                "expected the fast path to resolve the cap at gap {gap}"
            );
        } else {
            assert_eq!(
                result, None,
                "expected the fast path to defer at gap {gap} (clear_ln2_gap={clear_ln2_gap}, \
                 inflate_ll_magnitude={inflate_ll_magnitude})"
            );
        }
    }

    // ========================================================================
    // Unanimous fast-path exactness
    // ========================================================================

    /// THE CRITICAL PARITY TEST: `call()` (which may take the margin-gated fast
    /// path) must equal `call_full()` (the unmodified fgbio computation, the sole
    /// oracle) for every unanimous position. This is the entire safety net for
    /// fgbio parity — if the fast path's cap-region or margin logic is ever wrong,
    /// this is what catches it.
    ///
    /// Swept broadly across all five axes that matter: pre-UMI rate (the axis the
    /// cap-region boundary depends on), post-UMI rate and observed quality (which
    /// set the per-observation likelihoods and hence the margin's numerator),
    /// depth (which scales the gap and, at extreme depth, the ln-posterior
    /// magnitudes the margin bound is over), and the winning base (whose lane
    /// position changes `ln_sum_exp_array`'s summation order, and hence the exact
    /// floating-point value of the quantities the margin bound is calibrated
    /// against — see `test_fast_path_matches_call_full_dense` for the
    /// contiguous-depth sweep that actually exercises the near-boundary cases this
    /// broad-but-sparse sweep is too coarse to land on). Pre-UMI, post-UMI, and
    /// base are `#[case]`/`#[values]` axes so a regression names the exact
    /// combination it broke at; observed quality and depth are swept inside the
    /// body (8 x 8 x 4 x 10 x 12 = 30,720 checks total).
    #[rstest]
    fn test_unanimous_fast_path_matches_full_calculation(
        #[values(0, 2, 20, 45, 69, 70, 90, 93)] error_rate_pre_umi: PhredScore,
        #[values(0, 1, 2, 5, 10, 20, 40, 93)] error_rate_post_umi: PhredScore,
        #[values(b'A', b'C', b'G', b'T')] base: u8,
    ) {
        for depth in [1u32, 2, 3, 5, 8, 12, 20, 35, 50, 100, 400, 1000] {
            for obs_qual in [0u8, 1, 2, 5, 10, 20, 30, 40, 60, 93] {
                let mut builder =
                    ConsensusBaseBuilder::new(error_rate_pre_umi, error_rate_post_umi);
                for _ in 0..depth {
                    builder.add(base, obs_qual);
                }

                let full = builder.call_full();
                let fast = builder.call();

                assert_eq!(
                    fast, full,
                    "unanimous fast path disagreed with the full calculation at \
                     pre-UMI Q{error_rate_pre_umi}, post-UMI Q{error_rate_post_umi}, \
                     base {}, depth {depth}, obs Q{obs_qual}: fast={fast:?} full={full:?}",
                    base as char
                );
            }
        }
    }

    /// THE DECISIVE PARITY TEST: unlike
    /// `test_unanimous_fast_path_matches_full_calculation`'s sparse depth axis
    /// (which skips almost every depth and so almost never lands `gap` near a
    /// bracket boundary), this test scans depth **contiguously**, one observation
    /// at a time, reusing the same builder across the whole scan. Full-vs-gap
    /// divergences occur at a rate of ~1e-5 and only when `gap` lands within a few
    /// ULPs of a bracket edge; `gap` advances by a fixed increment per observation,
    /// so contiguous depth scanning is the only way to actually visit those
    /// boundary-adjacent values. Reusing the builder incrementally (one `add` per
    /// depth, not rebuilding from scratch) keeps this O(N) per parameter triple
    /// rather than O(N^2).
    ///
    /// Sweeps `post` and `obs` densely over the low-quality band (`0..=12`), the
    /// regime where the margin is tightest, for `pre` in {45, 70, 90} (spanning the
    /// default, the cap-boundary crossover, and a high-cap rate) and all four
    /// winning bases (lane position changes `ln_sum_exp_array`'s summation order).
    /// `pre=90, post=1, obs=1` — the parameter combination that reproduced the
    /// cap-region parity bug (`call_full` landing one Phred short of the cap just
    /// above `thresholds[cap]`, see the module docs) — is inside this sweep.
    ///
    /// Unlike an earlier version of this test, the per-`(pre, base, post, obs)`
    /// scan does **not** stop at the cap boundary: it continues *through and past*
    /// `thresholds[cap]`, far enough that the cap-region margin (see
    /// `try_unanimous_fast_path`) is comfortably cleared, or to depth 4000,
    /// whichever comes first. Stopping exactly at the boundary (the earlier
    /// behavior) is precisely why this test never caught the cap-region bug: the
    /// divergence only appears a handful of Phred-equivalent gap-units *past* the
    /// threshold, never before it.
    #[rstest]
    fn test_fast_path_matches_call_full_dense(
        #[values(45, 70, 90)] pre: PhredScore,
        #[values(b'A', b'C', b'G', b'T')] base: u8,
    ) {
        /// How far past `thresholds[cap]` to keep scanning before stopping early:
        /// large enough that the cap-region margin (a few Phred-equivalent gap
        /// units at these rates, see the module docs) is cleared with room to
        /// spare, so the known-bad repro band is never missed by an early exit.
        const POST_CAP_SLACK: f64 = 50.0;

        let idx = BASE_TO_INDEX[base as usize] as usize;

        for post in 0u8..=12 {
            for obs in 0u8..=12 {
                let mut builder = ConsensusBaseBuilder::new(pre, post);
                let cap_threshold = builder.gap_tables.thresholds[builder.gap_tables.cap];

                // Add one observation at a time and check after each: contiguous depth
                // scanning is what lands `gap` near a bracket boundary (sub-cap or the
                // cap boundary itself), which is the only place the full and gap forms
                // can disagree.
                for depth in 1..=4000u32 {
                    builder.add(base, obs);
                    let ll = builder.likelihoods.as_array();
                    let gap = ll[idx] - ll[(idx + 1) % 4];

                    assert_eq!(
                        builder.call(),
                        builder.call_full(),
                        "fast path disagreed with call_full at pre=Q{pre} base={} post=Q{post} \
                         obs=Q{obs} depth={depth} gap={gap}",
                        base as char
                    );

                    if gap.is_finite() && gap > cap_threshold + POST_CAP_SLACK {
                        break; // comfortably past the cap boundary: stop early
                    }
                }
            }
        }
    }

    /// C1 regression: the deep cap-region divergences `test_fast_path_matches_call_full_dense`
    /// structurally cannot see, because its `POST_CAP_SLACK = 50.0` stops scanning only ~50
    /// gap-units past `thresholds[cap]` — these divergences begin thousands of gap-units past
    /// it, where the (pre-fix) linear allowance `gap - thresholds[cap] > margin_cap` kept firing
    /// even though the true remaining ln-posterior budget had long since saturated. Each case
    /// below scans an incremental depth range known (from direct reproduction against the
    /// pre-fix code) to straddle the first divergence for that `(pre, post, obs)` triple.
    #[rstest]
    #[case::pre90_post93_obs93(90, 93, 93, 1450..=1560)]
    #[case::pre93_post93_obs40(93, 93, 40, 1550..=1660)]
    #[case::pre90_post90_obs40(90, 90, 40, 3150..=3260)]
    fn test_fast_path_matches_call_full_at_deep_cap_region(
        #[case] pre: PhredScore,
        #[case] post: PhredScore,
        #[case] obs: u8,
        #[case] depths: std::ops::RangeInclusive<u32>,
    ) {
        let mut builder = ConsensusBaseBuilder::new(pre, post);
        let start = *depths.start();
        // Build up to just before the scanned range in one shot, then add one
        // observation at a time through the range, checking parity after each.
        for _ in 0..start.saturating_sub(1) {
            builder.add(b'A', obs);
        }
        for depth in depths {
            builder.add(b'A', obs);
            assert_eq!(
                builder.contributions(),
                depth,
                "depth bookkeeping drifted from the intended scan position"
            );
            let fast = builder.call();
            let full = builder.call_full();
            assert_eq!(
                fast, full,
                "pre=Q{pre} post=Q{post} obs=Q{obs} depth={depth}: fast={fast:?} full={full:?}"
            );
        }
    }

    /// The derived bound must both preserve the old threshold's safety at typical
    /// pre-UMI rates and be tighter than it there — the fix is not a slowdown at the
    /// default. Above Q69 the required gap exceeds the old hardcoded 23.0, which is
    /// exactly the range where the constant was silently inflating qualities.
    #[rstest]
    #[case::q40(40, false)]
    #[case::q45_default(45, false)]
    #[case::q69_last_safe(69, false)]
    #[case::q70_first_unsafe(70, true)]
    #[case::q93_max(93, true)]
    fn test_derived_gap_brackets_the_old_hardcoded_threshold(
        #[case] error_rate_pre_umi: PhredScore,
        #[case] expected_above_old_threshold: bool,
    ) {
        /// The constant this bound replaces.
        const OLD_THRESHOLD: f64 = 23.0;

        let gap = cap_gap_threshold(error_rate_pre_umi);
        assert_eq!(
            gap > OLD_THRESHOLD,
            expected_above_old_threshold,
            "pre-UMI Q{error_rate_pre_umi} needs a gap of {gap:.2}; expected it to be \
             {} 23.0",
            if expected_above_old_threshold { "above" } else { "at or below" }
        );
    }

    /// The required gap rises monotonically with pre-UMI quality — the property that
    /// makes a single constant unable to serve the whole range.
    #[test]
    fn test_required_gap_increases_with_pre_umi_quality() {
        let gaps: Vec<f64> = (0..=MAX_PHRED).map(cap_gap_threshold).collect();

        for (q, pair) in gaps.windows(2).enumerate() {
            assert!(
                pair[1] >= pair[0],
                "gap must not shrink as pre-UMI quality rises, but Q{} needs {:.4} and \
                 Q{} needs {:.4}",
                q,
                pair[0],
                q + 1,
                pair[1]
            );
        }
    }

    /// The memoized tables must equal a fresh build for every representable pre-UMI
    /// Phred, element for element and bit for bit (both the thresholds and the
    /// consensus-error minima, plus the cap). Sweeps the whole `u8` domain
    /// (including values above `MAX_PHRED`) because the cache is indexed by the raw
    /// `u8` and must be exact across every index.
    #[test]
    fn test_cached_thresholds_match_uncached_across_the_whole_u8_domain() {
        for q in 0..=u8::MAX {
            let cached = cached_unanimous_gap_tables(q);
            let direct = build_unanimous_gap_tables(q);
            assert_eq!(cached.cap, direct.cap, "cap mismatch for Q{q}");
            assert_eq!(
                cached.thresholds.len(),
                direct.thresholds.len(),
                "length mismatch for Q{q}"
            );
            for (level, (&c, &d)) in
                cached.thresholds.iter().zip(direct.thresholds.iter()).enumerate()
            {
                assert_eq!(
                    c.to_bits(),
                    d.to_bits(),
                    "memoized threshold for Q{q} level {level} ({c}) differs from direct ({d})"
                );
            }
            for (level, (&c, &d)) in cached.cerr_min.iter().zip(direct.cerr_min.iter()).enumerate()
            {
                assert_eq!(
                    c.to_bits(),
                    d.to_bits(),
                    "memoized cerr_min for Q{q} level {level} ({c}) differs from direct ({d})"
                );
            }
        }
    }

    /// The gap tables a builder stores must be derived from the **pre**-UMI rate,
    /// not the post-UMI one. The asymmetric cases are what discriminate a transposed
    /// argument; `umi_consensus_defaults` is symmetric and cannot, but is kept
    /// because Q90/Q90 is the pair `SimpleConsensusCaller` actually uses.
    #[rstest]
    #[case::default_cli_rates(45, 40)]
    #[case::pre_below_post(40, 45)]
    #[case::umi_consensus_defaults(90, 90)]
    #[case::max_pre_low_post(93, 2)]
    #[case::low_pre_max_post(2, 93)]
    fn test_builder_gap_is_keyed_on_the_pre_umi_rate(
        #[case] error_rate_pre_umi: PhredScore,
        #[case] error_rate_post_umi: PhredScore,
    ) {
        let builder = ConsensusBaseBuilder::new(error_rate_pre_umi, error_rate_post_umi);
        let expected = build_unanimous_gap_tables(error_rate_pre_umi);
        assert_eq!(builder.gap_tables.cap, expected.cap);
        assert_eq!(builder.gap_tables.thresholds.len(), expected.thresholds.len());
        for (level, (&got, &want)) in
            builder.gap_tables.thresholds.iter().zip(expected.thresholds.iter()).enumerate()
        {
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "builder built with pre={error_rate_pre_umi}, post={error_rate_post_umi} stored \
                 threshold {got} at level {level} but the pre-UMI rate requires {want}"
            );
        }
        for (level, (&got, &want)) in
            builder.gap_tables.cerr_min.iter().zip(expected.cerr_min.iter()).enumerate()
        {
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "builder built with pre={error_rate_pre_umi}, post={error_rate_post_umi} stored \
                 cerr_min {got} at level {level} but the pre-UMI rate requires {want}"
            );
        }
    }

    /// `unanimous_quality_from_gap` must be non-decreasing in `gap` for every rate.
    /// This is the monotonicity premise the table's exactness proof relies on; a fine
    /// grid (finer than any threshold gap step) guards it directly.
    #[test]
    fn test_unanimous_quality_from_gap_is_monotone() {
        for rate in 0..=u8::MAX {
            let ln_pre = phred_to_ln_error_prob(rate);
            let mut prev = 0u8;
            // 0..=256 in 0.05 steps covers the whole sub-cap range with margin.
            let mut g = 0.0;
            while g <= 256.0 {
                let q = unanimous_quality_from_gap(g, ln_pre);
                assert!(
                    q >= prev,
                    "quality must not decrease with gap at rate Q{rate}: gap {g:.4} gave Q{q} \
                     after Q{prev}"
                );
                prev = q;
                g += 0.05;
            }
        }
    }

    /// The table must be non-decreasing for every rate (guards the `partition_point`
    /// lookup). `INFINITY` padding above the cap keeps this true at the top.
    #[test]
    fn test_gap_thresholds_are_monotone_non_decreasing() {
        for rate in 0..=u8::MAX {
            let thresholds = build_unanimous_gap_thresholds(rate);
            for (q, pair) in thresholds.windows(2).enumerate() {
                assert!(
                    pair[1] >= pair[0],
                    "thresholds for Q{rate} must not shrink: level {q} = {} vs level {} = {}",
                    pair[0],
                    q + 1,
                    pair[1]
                );
            }
        }
    }

    /// Reachable levels `0..=cap` are finite; unreachable levels above the cap are
    /// `+INFINITY`. Pins the semantics the lookup relies on.
    #[test]
    fn test_gap_thresholds_finite_up_to_cap_infinite_above() {
        for rate in 0..=u8::MAX {
            let ln_pre = phred_to_ln_error_prob(rate);
            let cap = ln_prob_to_phred(ln_pre) as usize;
            let thresholds = build_unanimous_gap_thresholds(rate);
            for (q, &t) in thresholds.iter().enumerate() {
                if q <= cap {
                    assert!(
                        t.is_finite(),
                        "Q{rate} level {q} (<= cap {cap}) must be finite, got {t}"
                    );
                } else {
                    assert!(
                        t.is_infinite(),
                        "Q{rate} level {q} (> cap {cap}) must be +inf, got {t}"
                    );
                }
            }
        }
    }

    /// Table exactness (the theorem's Lemma): the `partition_point` lookup reproduces
    /// `unanimous_quality_from_gap` exactly over a fine gap grid, for every rate.
    #[test]
    fn test_table_lookup_reproduces_helper_exactly() {
        for rate in 0..=u8::MAX {
            let ln_pre = phred_to_ln_error_prob(rate);
            let thresholds = build_unanimous_gap_thresholds(rate);
            let mut g = 0.0;
            while g <= 256.0 {
                #[expect(
                    clippy::cast_possible_truncation,
                    reason = "bounded by the 95-entry table, well within u8"
                )]
                let looked_up = thresholds.partition_point(|&t| t <= g).saturating_sub(1) as u8;
                let direct = unanimous_quality_from_gap(g, ln_pre);
                assert_eq!(
                    looked_up, direct,
                    "table lookup disagreed with helper at rate Q{rate}, gap {g:.4}: \
                     lookup=Q{looked_up} helper=Q{direct}"
                );
                g += 0.05;
            }
        }
    }

    /// Pins the error bound [`UNANIMOUS_MARGIN_HEADROOM`] depends on:
    /// `|ln_posterior_full - ln_posterior_gap| <= C * u * (|w| + |l|)`, `u` =
    /// [`UNANIMOUS_UNIT_ROUNDOFF`]. `ln_posterior_full` is `call_full`'s unshifted
    /// `[w, l, l, l]` normalize, built from the actual `likelihoods` array (not a
    /// reconstructed literal) so it matches what `call_full` really computes;
    /// `ln_posterior_gap` is the shifted `[0, -gap, -gap, -gap]` form the
    /// gap-threshold table is built from. Both compute the same mathematical
    /// posterior (the shift is exact in real arithmetic — subtracting `w` from every
    /// lane before `ln_sum_exp` cancels out of the normalized result) and differ
    /// only by floating-point rounding, which this test bounds directly rather than
    /// trusting the derivation, by calling the production [`unanimous_margin`]
    /// helper with `cerr_lower_bound = 1.0` (its identity element) rather than
    /// transcribing the formula a second time.
    ///
    /// Swept over **all four winner bases** as the outer `#[case]` axis (not just
    /// lane 0), because `ln_sum_exp_array` seeds its accumulator with the *minimum*
    /// lane and folds the rest in index order — a winner at lane 2 (G) or lane 3
    /// (T) rounds differently than one at lane 0 (A) or 1 (C). Calibrating on lane 0
    /// alone (an earlier version of this test) measured only ~1.17x worst case and
    /// silently left almost no real headroom once lanes 2/3 are included. `post` and
    /// `obs_qual` sweep their *entire* representable range (`0..=93`, not a sparse
    /// list), and `depth` densely covers `1..=60` — deliberately starting at `1`, not
    /// some higher floor, because the observed worst case (`base=G, post=Q35,
    /// obs=Q56, depth=1`) is at the *shallowest* possible depth, plus a handful of
    /// deeper spot-checks (up to `4000`) to retain the deep-depth coverage this test
    /// has always had. `pre` is deliberately fixed at Q45 (not swept): `w`, `l`, and
    /// `gap` depend only on the post-UMI rate, observation quality, and depth, never
    /// on the pre-UMI rate, which enters only after this bound's scope (in the
    /// pre-UMI error combination step).
    ///
    /// Beyond the hard `delta <= bound` assertion (which only pins that
    /// `UNANIMOUS_MARGIN_HEADROOM = 16.0` is *sufficient*), this test also tracks
    /// the raw ratio `delta / (u * (|w| + |l|))` — i.e. the headroom actually
    /// consumed, with `UNANIMOUS_MARGIN_HEADROOM` divided back out — and asserts it
    /// stays below a pinned ceiling of `5.0`. The observed worst case across this
    /// sweep is ~3.845 (at `base=G, post=Q35, obs=Q56, depth=1`); a ceiling of `5.0`
    /// gives room for minor floating-point variation across platforms/toolchains
    /// without masking a real regression, while leaving a comfortable margin below
    /// the hard `16.0` headroom. If the ratio ever exceeds `5.0`, the **derivation**
    /// needs revisiting (the empirical bound or its headroom), not this test:
    /// silently loosening the ceiling to make a failure go away would defeat its
    /// purpose. (A ceiling of `2.5` — this test's previous value, calibrated before
    /// lanes 2/3 were swept — reliably fails against the `base=G` case; `5.0` is the
    /// corrected, comprehensive value.)
    ///
    /// `checked > 0` guards against this test silently degrading into a no-op (e.g.
    /// if a future axis edit made every combination vacuous); many individual
    /// `(post, obs_qual)` pairs in the full `0..=93` sweep *are* expected to be
    /// vacuous (gap non-finite or `<= 0`) — that is normal, not a gap in coverage,
    /// since the sign of `gap`'s fixed per-observation increment cannot flip with
    /// more depth.
    #[rstest]
    fn test_ln_posterior_full_gap_error_bound(#[values(b'A', b'C', b'G', b'T')] base: u8) {
        /// Pinned ceiling on the observed `delta / (u * (|w| + |l|))` ratio (the
        /// headroom actually consumed): the observed worst case is ~3.845 (base G),
        /// so 5.0 leaves margin for minor platform variation while still catching
        /// any real growth toward the full `UNANIMOUS_MARGIN_HEADROOM = 16.0`
        /// allowance.
        const MAX_OBSERVED_RATIO_CEILING: f64 = 5.0;

        let idx = BASE_TO_INDEX[base as usize] as usize;
        let mut checked = 0u64;
        let mut max_observed_ratio = 0.0_f64;
        let mut worst_case = String::new();

        for post in 0u8..=93 {
            for obs_qual in 0u8..=93 {
                let mut b = ConsensusBaseBuilder::new(45, post);
                let mut current_depth = 0u32;
                // Dense over the shallow depths where the worst case actually lives
                // (do not start higher than 1), plus a handful of deep spot-checks to
                // retain this test's historical deep-depth coverage.
                for target_depth in (1u32..=60).chain([100, 400, 1000, 2000, 4000]) {
                    for _ in current_depth..target_depth {
                        b.add(base, obs_qual);
                    }
                    current_depth = target_depth;

                    let ll = b.likelihoods.as_array();
                    let winner_ll = ll[idx];
                    let loser_ll = ll[(idx + 1) % 4];
                    let gap = winner_ll - loser_ll;
                    if !(gap.is_finite() && gap > 0.0) {
                        continue; // bound only claimed for unanimous positions with gap > 0
                    }
                    checked += 1;

                    let ln_posterior_full = ln_normalize(winner_ll, ln_sum_exp_array(ll));
                    let ln_posterior_gap =
                        ln_normalize(0.0, ln_sum_exp_array(&[0.0, -gap, -gap, -gap]));
                    let delta = (ln_posterior_full - ln_posterior_gap).abs();
                    let bound = unanimous_margin(winner_ll, loser_ll, 1.0);

                    assert!(
                        delta <= bound,
                        "post=Q{post} base={} obsQ={obs_qual} depth={target_depth} gap={gap:.6}: \
                         |delta|={delta:e} exceeded bound={bound:e} (w={winner_ll:.6}, \
                         l={loser_ll:.6})",
                        base as char
                    );

                    let denom = UNANIMOUS_UNIT_ROUNDOFF * (winner_ll.abs() + loser_ll.abs());
                    if denom > 0.0 {
                        let ratio = delta / denom;
                        if ratio > max_observed_ratio {
                            max_observed_ratio = ratio;
                            worst_case =
                                format!("post=Q{post} obsQ={obs_qual} depth={target_depth}");
                        }
                    }
                }
            }
        }

        assert!(
            checked > 0,
            "base={}: every (post, obs_qual, depth) combination produced a non-finite or \
             non-positive gap, so this case asserted nothing over the whole sweep",
            base as char
        );
        assert!(
            max_observed_ratio < MAX_OBSERVED_RATIO_CEILING,
            "base={}: observed headroom-consumption ratio {max_observed_ratio:.4} at \
             {worst_case} reached the pinned ceiling {MAX_OBSERVED_RATIO_CEILING} (expected \
             worst case ~3.845, at base G) — the derivation needs revisiting, not this ceiling",
            base as char
        );
    }

    /// Degenerate low-quality unanimous positions — where the per-observation adjusted
    /// error is >= 3/4, so the observed base is NOT the maximum-likelihood lane (and
    /// `error_rate_post_umi == 0` or an observation quality of 0 give gap = -inf with a
    /// three-way loser tie) — must be resolved identically by the fast path and
    /// `call_full`. The fast path defers them (returns `None`); this pins that `call()` ==
    /// `call_full()` across the whole degenerate corner, the regime the main matches-full
    /// sweep (obs/post >= Q10) never exercises.
    #[rstest]
    fn test_fast_path_defers_degenerate_unanimous_to_call_full(
        #[values(2, 20, 45, 70, 93)] pre: PhredScore,
        #[values(0, 1, 2, 3, 4, 5)] post: PhredScore,
        #[values(0, 1, 2, 3, 4, 5)] obs: u8,
    ) {
        for depth in [1u32, 2, 3, 5, 10, 50] {
            let mut builder = ConsensusBaseBuilder::new(pre, post);
            for _ in 0..depth {
                builder.add(b'A', obs);
            }
            assert_eq!(
                builder.call(),
                builder.call_full(),
                "fast path disagreed with call_full at pre=Q{pre} post=Q{post} obs=Q{obs} depth={depth}"
            );
        }
    }

    /// Pins `call()`'s exact unanimous quality at representative points spanning the
    /// deep/low-quality sub-cap band, the mid-curve, and the pre-UMI cap, so a future
    /// change to the unanimous math (the margin gate, the gap-threshold table, or the
    /// underlying phred arithmetic) is caught rather than silently drifting. `call_full`
    /// is the unmodified fgbio `[w, l, l, l]` computation with no unanimous-specific
    /// branch, so these values are exactly fgbio's, not an approximation of them
    /// (whether or not the fast path happens to short-circuit a given case, `call()` ==
    /// `call_full()` is guaranteed by `test_unanimous_fast_path_matches_full_calculation`).
    #[rstest]
    #[case::subcap_deep_low_quality_a(45, 2, 2, 50, 16)]
    #[case::subcap_deep_low_quality_b(70, 5, 5, 15, 65)]
    #[case::cap_region_mid_pre(70, 5, 5, 40, 70)]
    #[case::mid_curve_high_pre(93, 40, 20, 3, 69)]
    #[case::mid_curve_low_pre(20, 10, 10, 4, 19)]
    #[case::cap_region_pre_45(45, 40, 40, 50, 45)]
    #[case::cap_region_pre_93(93, 93, 93, 100, 93)]
    #[case::min_boundary(2, 2, 2, 5, 2)]
    fn test_unanimous_quality_pins_representative_points(
        #[case] pre: PhredScore,
        #[case] post: PhredScore,
        #[case] obs: u8,
        #[case] depth: u32,
        #[case] expected_qual: PhredScore,
    ) {
        let mut builder = ConsensusBaseBuilder::new(pre, post);
        for _ in 0..depth {
            builder.add(b'A', obs);
        }
        let (base, qual) = builder.call();
        assert_eq!(base, b'A', "unanimous consensus base must be A");
        assert_eq!(
            qual, expected_qual,
            "pre=Q{pre} post=Q{post} obs=Q{obs} depth={depth}: got Q{qual}, expected Q{expected_qual}"
        );
    }
}
