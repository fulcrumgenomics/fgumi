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
pub(crate) const FD_RESERVE: u64 = 32;

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
/// Set from a measured sweep of `--max-temp-files` against spilled-run counts of
/// 86 / 173 / 346 / 692 on a 780M-record WGS BAM (8 threads, coordinate order),
/// obtained by shrinking the total memory budget from 4 GiB to 512 MiB:
///
/// | spilled runs | consolidating (limit 64) | never consolidating | speedup |
/// | --- | --- | --- | --- |
/// | 86 | 789.8s | 596.4s | 1.32x |
/// | 173 | 1540.8s | 607.1s | 2.54x |
/// | 346 | 3026.1s | 619.7s | 4.88x |
///
/// Consolidation cost grows steeply because the number of passes is about
/// `(runs - limit) / (limit / 2)` and each pass rewrites the whole accumulated
/// output: at 346 runs it consumed 2440s of a 3026s sort, 80% of wall clock.
/// A ceiling of 256 would cap the derived limit below the run count in the two
/// deepest cells and re-impose exactly that cost, so it is set above the
/// measured range rather than at it.
///
/// The competing cost is that each spill file open at the final merge carries
/// per-file read-ahead not charged against `--max-memory`. Measured, that is far
/// smaller than a naive per-file estimate suggests, because peak RSS is normally
/// set by the phase-1 sort buffer and the merge's buffers land after that memory
/// is released: widening 346 runs from consolidating to not cost +104 MB. It
/// only bites at the smallest budgets, where phase 1 is too small to dominate --
/// at 692 runs on a 512 MiB total budget, a 1024-wide merge reached 2309 MB
/// against 1257 MB for a 256-wide one. Sizing phase-2 read-ahead as a total
/// budget divided across files, rather than a per-file constant, would remove
/// that and is the reason this is a ceiling rather than no cap at all.
const MAX_DERIVED_TEMP_FILES: usize = 1024;

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
/// Split from [`temp_file_limit_from_nofile`] so the arithmetic is testable
/// without touching the process's real limits. Saturating throughout: a soft
/// limit below the reserve yields the floor rather than underflowing, and
/// `RLIM_INFINITY` (`u64::MAX`) yields the ceiling rather than overflowing.
fn temp_file_limit_from_soft_nofile(soft: u64) -> usize {
    let budget = usize::try_from(soft.saturating_sub(FD_RESERVE)).unwrap_or(usize::MAX);
    budget.clamp(MIN_TEMP_FILES, MAX_DERIVED_TEMP_FILES)
}

/// A consolidation limit derived from an already-taken [`soft_nofile`] reading.
///
/// Takes the snapshot rather than reading the limit itself, because a caller
/// generally needs more than one answer about the same budget — the limit to
/// hand the sorter, the budget to name in a log line, and whether the limit
/// fits. Reading `RLIMIT_NOFILE` once per question lets those answers disagree
/// if the limit changes in between, so the caller reads it once and derives
/// all three from that one value.
#[must_use]
pub fn temp_file_limit_from_nofile(soft: Option<u64>) -> usize {
    soft.map_or(FALLBACK_MAX_TEMP_FILES, temp_file_limit_from_soft_nofile)
}

/// A consolidation limit derived from the process's descriptor budget.
///
/// Falls back to [`FALLBACK_MAX_TEMP_FILES`] on a non-Unix target.
///
/// The convenience form, for a caller that wants only the number. A caller that
/// also reports or checks the budget should take one [`soft_nofile`] reading and
/// use [`temp_file_limit_from_nofile`] and [`fits_nofile_budget`] instead, so
/// all three answers describe the same budget.
///
/// This is deliberately silent — see [`fits_nofile_budget`] for the check a caller
/// should warn on.
///
/// # Concurrency
///
/// The whole descriptor budget is assumed to be available to one sort. A host
/// running several sorts concurrently in one process should divide the budget
/// itself rather than calling this once per sort.
#[must_use]
pub fn resolve_temp_file_limit() -> usize {
    temp_file_limit_from_nofile(soft_nofile())
}

/// Whether `limit` open files fit within the descriptor budget `soft` describes.
///
/// Deliberately independent of where `limit` came from: an explicitly requested
/// limit can overrun the budget just as a derived one can, and in practice it is
/// more likely to — [`temp_file_limit_from_nofile`] clamps to the budget, while
/// a hardcoded value cannot. A caller that only checked derived limits would
/// therefore skip the case most likely to fail.
///
/// Takes a [`soft_nofile`] snapshot for the reason given on
/// [`temp_file_limit_from_nofile`]: the limit being judged here was derived
/// from a budget, and judging it against a second, independently-read budget is
/// how the two come to disagree.
///
/// Always true when `soft` is `None`, i.e. on a target where the limit cannot
/// be read.
#[must_use]
pub fn fits_nofile_budget(soft: Option<u64>, limit: usize) -> bool {
    let Some(soft) = soft else { return true };
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
    // A typical Linux container: below the ceiling, so it passes through.
    #[case::linux_container(1024, 992)]
    // The limit actually observed on the AWS Batch hosts these run on, which is
    // well above the ceiling and therefore clamps.
    #[case::batch_host(65_536, MAX_DERIVED_TEMP_FILES)]
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
        assert_eq!(
            fits_nofile_budget(Some(soft), temp_file_limit_from_soft_nofile(soft)),
            expected
        );
    }

    /// An explicit limit is checked against the budget on the same terms as a
    /// derived one. This is the case that matters in practice: callers that pin
    /// a large limit are the ones who can overrun the budget.
    #[test]
    fn test_fits_nofile_budget_rejects_an_oversized_explicit_limit() {
        assert!(!fits_nofile_budget(Some(1024), usize::MAX));
    }

    /// `fits_nofile_budget` is also the merge path's check, where the count is
    /// an input-file count rather than a spill-file limit. A handful of inputs
    /// must always pass, or every ordinary `fgumi merge` would warn.
    #[test]
    fn test_fits_nofile_budget_accepts_an_ordinary_merge_width() {
        assert!(fits_nofile_budget(Some(1024), 8));
    }

    /// An unreadable budget admits everything, so a non-Unix host never warns.
    #[test]
    fn test_fits_nofile_budget_admits_everything_without_a_budget() {
        assert!(fits_nofile_budget(None, usize::MAX));
    }

    /// An unreadable budget yields the portable fallback rather than a derived
    /// number, which is the whole non-Unix contract.
    #[test]
    fn test_temp_file_limit_from_nofile_falls_back_without_a_budget() {
        assert_eq!(temp_file_limit_from_nofile(None), FALLBACK_MAX_TEMP_FILES);
    }
}
