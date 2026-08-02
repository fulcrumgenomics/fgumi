//! Derivation of a spill-file consolidation limit from the process's file
//! descriptor budget.
//!
//! The external sort consolidates spilled runs once their on-disk count reaches
//! a limit, merging the oldest half into one. That limit exists to bound the
//! number of descriptors held open at once: the final k-way merge opens *every*
//! remaining spill file simultaneously, and consolidation itself opens half the
//! limit. Consolidation is otherwise pure overhead — it rewrites data that is
//! already sorted — so the limit should be as high as the process's descriptor
//! budget genuinely allows, and no higher.
//!
//! This module is a **helper a caller opts into**, not a policy the engine
//! applies. [`RawExternalSorter`](crate::RawExternalSorter) keeps a fixed,
//! portable default so that an embedder who never thinks about the knob gets
//! the same behaviour on every host; deciding to derive the limit from the
//! environment instead belongs to whoever owns the command line. `fgumi sort`
//! makes that decision, so its default is derived and every other sorter
//! construction site is unaffected.
//!
//! Nothing here *raises* the descriptor limit. A process may raise its own soft
//! limit to the hard limit without privileges, but `fgumi-sort` is a library and
//! mutating a process-global resource limit is a side effect its host did not
//! ask for.

/// Limit used when the process descriptor limit cannot be read, or on a
/// non-Unix target. Also the fixed default of a bare
/// [`RawExternalSorter`](crate::RawExternalSorter).
///
/// Matches samtools' implementation, which is where fgumi's original hardcoded
/// default came from.
pub const FALLBACK_MAX_TEMP_FILES: usize = 64;

/// Descriptors reserved for everything the sort holds open besides spill files.
///
/// At peak that is stdio (3), the input BAM, the output BAM and the optional
/// index writer — roughly 7. Nothing scales with `--threads`: the worker pool
/// shares one input reader, spill files are written one at a time, and the
/// temp-directory allocator holds paths rather than descriptors. 32 leaves
/// substantial headroom over the observed peak.
const FD_RESERVE: u64 = 32;

/// Smallest limit worth consolidating at.
///
/// Not the smallest value the consolidation path *accepts* (that is 2) but the
/// smallest at which it still amortizes. At a limit of 2 every spill takes the
/// run count to 2 and immediately merges both, rewriting the entire accumulated
/// output once per spill — quadratic in the number of runs. A host whose
/// descriptor budget cannot cover this floor cannot usefully sort a spilling
/// BAM, and is better served by a loud `EMFILE` than by a silent hundredfold
/// slowdown.
const MIN_TEMP_FILES: usize = 16;

// The floor must stay clear of the degenerate regime described above; a limit of
// 2 or 3 merges on essentially every spill.
const _: () = assert!(MIN_TEMP_FILES > 3, "floor must not sit in the quadratic regime");

/// Ceiling on the derived limit.
///
/// Above this the limit stops binding on any realistic input: the largest
/// workload measured spilled 88 runs at 1.33 billion records, and reaching the
/// ceiling needs either a multi-terabyte input or a pathologically small
/// `--max-memory`. Capping keeps a host with an unusually generous `ulimit -n`
/// inside a merge width that has actually been tested, and bounds the per-file
/// read-ahead buffers the merge allocates (a few MB each, none of which is
/// charged against `--max-memory`).
const MAX_DERIVED_TEMP_FILES: usize = 256;

const _: () = assert!(MIN_TEMP_FILES <= MAX_DERIVED_TEMP_FILES);

/// The process's soft `RLIMIT_NOFILE`, or `None` on a non-Unix target.
///
/// Exposed so a caller that derives its own limit can report which budget the
/// number came from.
#[must_use]
pub fn soft_nofile() -> Option<u64> {
    #[cfg(unix)]
    {
        // `getrlimit` is infallible for a valid resource. A `None` current value
        // means no soft limit is enforced, which we treat as an enormous one and
        // let the ceiling clamp.
        let limit = rustix::process::getrlimit(rustix::process::Resource::Nofile);
        Some(limit.current.unwrap_or(u64::MAX))
    }
    #[cfg(not(unix))]
    {
        None
    }
}

/// Clamp a soft `RLIMIT_NOFILE` value to a usable consolidation limit.
///
/// Split from [`resolve_temp_file_limit`] so the arithmetic is testable without
/// touching the process's real limits. Saturating throughout: a soft limit below
/// the reserve yields the floor rather than underflowing, and `RLIM_INFINITY`
/// (`u64::MAX`) yields the ceiling rather than overflowing.
fn temp_file_limit_from_soft_nofile(soft: u64) -> usize {
    let budget = usize::try_from(soft.saturating_sub(FD_RESERVE)).unwrap_or(usize::MAX);
    budget.clamp(MIN_TEMP_FILES, MAX_DERIVED_TEMP_FILES)
}

/// A consolidation limit derived from the process's descriptor budget.
///
/// Falls back to [`FALLBACK_MAX_TEMP_FILES`] on a non-Unix target.
///
/// This is deliberately silent — see [`fits_fd_budget`] for the check a caller
/// should warn on.
///
/// # Concurrency
///
/// The whole descriptor budget is assumed to be available to one sort. A host
/// running several sorts concurrently in one process should divide the budget
/// itself rather than calling this once per sort.
#[must_use]
pub fn resolve_temp_file_limit() -> usize {
    soft_nofile().map_or(FALLBACK_MAX_TEMP_FILES, temp_file_limit_from_soft_nofile)
}

/// Whether `limit` spill files fit within the process's descriptor budget.
///
/// Deliberately independent of where `limit` came from: an explicitly requested
/// limit can overrun the budget just as a derived one can, and in practice it is
/// more likely to — [`resolve_temp_file_limit`] clamps to the budget, while a
/// hardcoded value cannot. A caller that only checked derived limits would
/// therefore skip the case most likely to fail.
///
/// Always true on a target where the limit cannot be read.
#[must_use]
pub fn fits_fd_budget(limit: usize) -> bool {
    let Some(soft) = soft_nofile() else { return true };
    u64::try_from(limit).is_ok_and(|limit| limit <= soft.saturating_sub(FD_RESERVE))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// The clamp is the whole contract: saturate at both ends, and never return
    /// a value outside `[MIN_TEMP_FILES, MAX_DERIVED_TEMP_FILES]`.
    #[rstest]
    // A budget smaller than the reserve saturates to zero, then clamps up.
    #[case::zero(0, MIN_TEMP_FILES)]
    #[case::exactly_reserve(FD_RESERVE, MIN_TEMP_FILES)]
    // Lower-clamp boundary: the first soft limit that clears the floor on its
    // own merit, and the first that exceeds it.
    #[case::at_floor(FD_RESERVE + MIN_TEMP_FILES as u64, MIN_TEMP_FILES)]
    #[case::just_above_floor(FD_RESERVE + MIN_TEMP_FILES as u64 + 1, MIN_TEMP_FILES + 1)]
    // A soft limit below the legacy hardcoded 64: the derived limit goes *down*,
    // which is the safety half of deriving it at all.
    #[case::below_legacy_default(64, 32)]
    // The historical macOS default.
    #[case::macos_default(256, 224)]
    // A typical Linux container, clamped by the ceiling.
    #[case::linux_container(1024, MAX_DERIVED_TEMP_FILES)]
    // `RLIM_INFINITY` clamps rather than overflowing.
    #[case::rlim_infinity(u64::MAX, MAX_DERIVED_TEMP_FILES)]
    fn test_temp_file_limit_from_soft_nofile(#[case] soft: u64, #[case] expected: usize) {
        assert_eq!(temp_file_limit_from_soft_nofile(soft), expected);
    }

    /// Whatever this host's `ulimit -n` happens to be, the resolved limit must
    /// land in a range the consolidation path can act on.
    ///
    /// The expected range depends on whether the budget is readable at all, so
    /// it is picked up front rather than folded into one disjunction: an
    /// `||` whose right side never runs on a Unix host would leave the fallback
    /// arm silently unasserted.
    #[test]
    fn test_resolve_temp_file_limit_is_in_range() {
        let limit = resolve_temp_file_limit();
        let expected = if soft_nofile().is_some() {
            MIN_TEMP_FILES..=MAX_DERIVED_TEMP_FILES
        } else {
            FALLBACK_MAX_TEMP_FILES..=FALLBACK_MAX_TEMP_FILES
        };
        assert!(expected.contains(&limit), "resolved limit {limit} out of range {expected:?}");
    }

    /// A derived limit fits its own budget by construction, *except* when the
    /// budget was too small to reach the floor — which is exactly the case the
    /// caller must warn about.
    ///
    /// Unix-only: there is no budget to derive from, or fit inside, elsewhere.
    #[test]
    #[cfg(unix)]
    fn test_resolved_limit_fits_unless_budget_is_below_the_floor() {
        let soft = soft_nofile().expect("a Unix host always reports a soft limit");
        let expected = soft.saturating_sub(FD_RESERVE) >= MIN_TEMP_FILES as u64;
        assert_eq!(fits_fd_budget(temp_file_limit_from_soft_nofile(soft)), expected);
    }

    /// An explicit limit is checked against the budget on the same terms as a
    /// derived one. This is the case that matters in practice: callers that pin
    /// a large limit are the ones who can overrun the budget.
    ///
    /// Unix-only: with no readable budget nothing can overrun it, and
    /// [`fits_fd_budget`] admits everything by design.
    #[test]
    #[cfg(unix)]
    fn test_fits_fd_budget_rejects_an_oversized_explicit_limit() {
        assert!(!fits_fd_budget(usize::MAX));
    }
}
