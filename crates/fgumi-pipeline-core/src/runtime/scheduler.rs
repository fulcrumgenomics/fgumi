//! Pluggable per-worker dispatch-order policy for the round-robin pool driver.
//!
//! The worker loop ([`run_worker_loop`](crate::runtime::run_worker_loop)) walks
//! each worker's *live* steps once per pass and runs the first that makes
//! progress. A [`Scheduler`] decides the ORDER of that walk — the only thing it
//! controls; it never changes which steps exist, the sticky source/sink
//! fast-path, or the Serial/Exclusive contention rules.
//!
//! Two policies ship:
//!
//! - [`ChainOrderScheduler`] (the default) — walk **upstream-first** (chain
//!   order). A worker attempts the earliest-in-chain step with work, favouring
//!   production. This is the historical behaviour; every command keeps it unless
//!   it opts out, so existing pipelines are byte-for-byte unaffected.
//! - [`DrainFirstScheduler`] — walk **downstream-first** (reverse chain order).
//!   A worker attempts the deepest step with work first, favouring *draining*
//!   buffered work before producing more. Combined with skip-on-Serial-contention
//!   this self-balances: for a Serial step fed by an N-way Parallel producer, one
//!   worker grabs the Serial drain (mutex) while the rest find it contended, skip,
//!   and fall through to the producer — so the drain overlaps production instead
//!   of starving behind it on the shared pool. Sources/sinks are unaffected (they
//!   run on the sticky fast-path, not the round-robin walk).
//!
//! This mirrors main's `Scheduler`-trait design (`BalancedChaseDrainScheduler`
//! among others), but generically over an arbitrary step list rather than a
//! fixed set of named BAM stages: the only lever exposed here is walk direction,
//! which is all the generic driver needs to express drain-first scheduling.

/// The order in which a worker attempts its live steps in one round-robin pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WalkDirection {
    /// Chain order (upstream → downstream): favour production.
    Forward,
    /// Reverse chain order (downstream → upstream): favour draining.
    Reverse,
}

/// A per-worker dispatch-order policy. Selected per pipeline via
/// [`PipelineConfig::with_scheduler`](crate::builder::PipelineConfig::with_scheduler)
/// and shared across all workers (the shipped policies are stateless).
pub trait Scheduler: Send + Sync + std::fmt::Debug {
    /// Direction to walk this worker's live steps this pass. Called once per
    /// round-robin pass, so it must be cheap.
    fn walk(&self) -> WalkDirection;

    /// Human-readable name for diagnostics / `--pipeline-stats`.
    fn name(&self) -> &'static str;
}

/// Default upstream-first (chain-order) walk. Preserves the historical dispatch
/// behaviour for every pipeline that does not opt into a different policy.
#[derive(Debug, Default, Clone, Copy)]
pub struct ChainOrderScheduler;

impl Scheduler for ChainOrderScheduler {
    #[inline]
    fn walk(&self) -> WalkDirection {
        WalkDirection::Forward
    }
    fn name(&self) -> &'static str {
        "chain-order"
    }
}

/// Downstream-first (reverse chain-order) walk — drain buffered work before
/// producing more. Opt-in per pipeline (e.g. the sort chain, to overlap the
/// serial boundary/key scan with the parallel inflate instead of starving it on
/// the shared pool).
#[derive(Debug, Default, Clone, Copy)]
pub struct DrainFirstScheduler;

impl Scheduler for DrainFirstScheduler {
    #[inline]
    fn walk(&self) -> WalkDirection {
        WalkDirection::Reverse
    }
    fn name(&self) -> &'static str {
        "drain-first"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_forward() {
        assert_eq!(ChainOrderScheduler.walk(), WalkDirection::Forward);
        assert_eq!(ChainOrderScheduler.name(), "chain-order");
    }

    #[test]
    fn drain_first_is_reverse() {
        assert_eq!(DrainFirstScheduler.walk(), WalkDirection::Reverse);
        assert_eq!(DrainFirstScheduler.name(), "drain-first");
    }

    #[test]
    fn walk_direction_maps_positions() {
        // The driver maps pass position `i` over `n` live steps to a step index:
        // Forward → i; Reverse → n-1-i. Pin that arithmetic here so the driver
        // and this policy cannot drift.
        let n = 5usize;
        let fwd: Vec<usize> = (0..n).collect();
        let rev: Vec<usize> = (0..n).map(|i| n - 1 - i).collect();
        assert_eq!(fwd, vec![0, 1, 2, 3, 4]);
        assert_eq!(rev, vec![4, 3, 2, 1, 0]);
    }
}
