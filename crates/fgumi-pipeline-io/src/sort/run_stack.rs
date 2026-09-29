//! Which spill runs to consolidate, and when — the `--max-temp-files` policy.
//!
//! [`RunStack`] holds the live spill runs in run order (oldest first) and
//! answers one question after every change: is there a contiguous range of runs
//! that should now be merged into one? It performs no I/O; `SpillWrite` does the
//! merging and reports back through [`RunStack::complete_merge`].
//!
//! # Policy
//!
//! Nothing is merged until the live-run count reaches the limit `L`. A sort that
//! never spills `L` runs is therefore never slowed by consolidation — the limit
//! is a bound on open files, not a target.
//!
//! At the limit, the stack merges the contiguous window of `2..=f` runs (with
//! `f = clamp(L / 2, 2, MAX_FAN_IN)`) that costs the fewest bytes per run
//! eliminated: merging `w` runs holding `B` bytes rewrites `B` bytes and frees
//! `w − 1` slots, so the window minimizing `B / (w − 1)` buys headroom most
//! cheaply. Ties prefer the wider window (more headroom per pass), then the newer
//! one. Merging repeats until the count is back under `L`.
//!
//! Small runs are therefore merged before large ones, the way an optimal merge
//! pattern does, and a large consolidated run is not rewritten again until
//! merging it is the cheapest way to make room. The legacy `RawExternalSorter`
//! engine (which the library and `fgumi simulate` still use) instead merges the
//! oldest `L / 2` runs each time the limit is hit and puts the result back at the
//! front, so its first consolidated run absorbs — and rewrites — the oldest
//! records on every later pass. When the limit is hit once (the common case of a
//! sort just past it) and `L ≤ 2 · MAX_FAN_IN`, both policies do the same single
//! `L / 2`-wide merge; above that this one merges at most `MAX_FAN_IN` runs at a
//! time. When the limit is hit repeatedly this one rewrites several times less.
//!
//! Only contiguous ranges are merged and the merged run takes the range's
//! position, so the final merge's tie-break order (run order) is unchanged and
//! the sorted output is identical to an unconsolidated sort.

use std::ops::Range;
use std::path::PathBuf;

/// Upper bound on a consolidation's fan-in. Each input holds an open descriptor
/// and a read buffer for the duration of the merge; beyond a few dozen inputs
/// the extra width saves little I/O and costs memory the budget does not see.
pub(crate) const MAX_FAN_IN: usize = 32;

/// One live spill run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SpillRun {
    /// Merge-order identity: the final merge orders sources by `file_id`.
    pub(crate) file_id: u32,
    /// The run's file.
    pub(crate) path: PathBuf,
    /// The run's size on disk; the cost of reading it into a consolidation.
    pub(crate) bytes: u64,
    /// Records ingested when the run's last chunk was spilled (progress only).
    pub(crate) records_ingested_so_far: u64,
}

/// The live spill runs, in run order, plus the consolidation policy over them.
#[derive(Debug)]
pub(crate) struct RunStack {
    runs: Vec<SpillRun>,
    /// Maximum live runs; `0` disables consolidation.
    limit: usize,
    /// Widest window a consolidation may merge.
    fan_in: usize,
}

impl RunStack {
    /// A stack enforcing `max_temp_files` live runs. A limit below 2 disables
    /// consolidation: one run cannot be consolidated, and the CLI rejects such
    /// values anyway, so only library callers can reach it.
    pub(crate) fn new(max_temp_files: usize) -> Self {
        let limit = if max_temp_files < 2 { 0 } else { max_temp_files };
        let fan_in = (limit / 2).clamp(2, MAX_FAN_IN);
        Self { runs: Vec::new(), limit, fan_in }
    }

    /// The live runs, oldest first.
    pub(crate) fn runs(&self) -> &[SpillRun] {
        &self.runs
    }

    /// Consume the stack, yielding the surviving runs oldest first.
    pub(crate) fn into_runs(self) -> Vec<SpillRun> {
        self.runs
    }

    /// Add a newly closed run. Runs must arrive in increasing `file_id` order.
    ///
    /// # Panics
    ///
    /// Panics (debug) if `run.file_id` does not exceed the newest live run's.
    pub(crate) fn push(&mut self, run: SpillRun) {
        debug_assert!(
            self.runs.last().is_none_or(|last| last.file_id < run.file_id),
            "spill runs must close in increasing file_id order"
        );
        self.runs.push(run);
    }

    /// The range of runs to merge next, or `None` while the stack is under its
    /// limit. Call repeatedly (merging each answer) until it returns `None`.
    pub(crate) fn next_merge(&self) -> Option<Range<usize>> {
        let n = self.runs.len();
        if self.limit == 0 || n < self.limit {
            return None;
        }
        // Prefix sums of run sizes, so every window's cost is O(1).
        let mut prefix = Vec::with_capacity(n + 1);
        prefix.push(0u128);
        for run in &self.runs {
            prefix.push(prefix.last().copied().unwrap_or(0) + u128::from(run.bytes));
        }
        // Best window so far as (bytes, width, start); compared by bytes per
        // eliminated slot, then width (wider wins), then start (newer wins).
        let mut best: Option<(u128, usize, usize)> = None;
        for width in 2..=self.fan_in.min(n) {
            for start in (0..=n - width).rev() {
                let bytes = prefix[start + width] - prefix[start];
                let better = match best {
                    None => true,
                    Some((best_bytes, best_width, _)) => {
                        // bytes / (width − 1) < best_bytes / (best_width − 1),
                        // cross-multiplied to stay in integers.
                        let lhs = bytes * (best_width as u128 - 1);
                        let rhs = best_bytes * (width as u128 - 1);
                        lhs < rhs || (lhs == rhs && width > best_width)
                    }
                };
                if better {
                    best = Some((bytes, width, start));
                }
            }
        }
        best.map(|(_, width, start)| start..start + width)
    }

    /// Replace the runs in `range` with the run merged from them — written to
    /// `merged_path`, `merged_bytes` long — and return the absorbed runs so their
    /// files can be removed.
    ///
    /// The merged run takes the range's lowest `file_id` — the position of the
    /// range in run order — so the final merge's ordering is unchanged.
    ///
    /// # Panics
    ///
    /// Panics if `range` is empty or out of bounds.
    pub(crate) fn complete_merge(
        &mut self,
        range: Range<usize>,
        merged_path: PathBuf,
        merged_bytes: u64,
    ) -> Vec<SpillRun> {
        assert!(!range.is_empty(), "cannot merge an empty range of runs");
        let start = range.start;
        let absorbed: Vec<SpillRun> = self.runs.drain(range).collect();
        let merged = SpillRun {
            file_id: absorbed[0].file_id,
            path: merged_path,
            bytes: merged_bytes,
            records_ingested_so_far: absorbed
                .iter()
                .map(|r| r.records_ingested_so_far)
                .max()
                .expect("non-empty"),
        };
        self.runs.insert(start, merged);
        absorbed
    }
}

#[cfg(test)]
mod tests;
