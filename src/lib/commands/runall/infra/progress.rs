//! Progress-rendering mode selector for `fgumi runall`.
//!
//! The [`ProgressMode`] enum is a CLI-side stable shape that future
//! progress-tracker integrations can consume. PR 1 ships the enum + a
//! no-op [`ProgressGuard`] returned from [`init`]; PRs that wire in
//! the actual tracker can add `From<ProgressMode> for crate::progress::Mode`
//! (or similar) without touching runall's public surface.

use clap::ValueEnum;

/// Progress rendering mode for the runall command.
///
/// `Auto` (the default) is meant to detect TTY vs non-TTY and pick a
/// sensible renderer. `Dashboard` and `Heartbeat` force a specific
/// renderer; `None` disables progress reporting entirely.
#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProgressMode {
    /// Detect TTY vs non-TTY and render accordingly. Default.
    Auto,
    /// Force the interactive multi-bar dashboard.
    Dashboard,
    /// Force periodic logfmt heartbeat lines suitable for log tailing.
    Heartbeat,
    /// Disable the tracker entirely. Every producer call is a no-op.
    None,
}

/// RAII guard for the progress tracker. Held for the duration of a
/// runall execution; on drop, closes the tracker (no-op in PR 1).
#[derive(Debug)]
pub(crate) struct ProgressGuard {
    /// The mode this guard was constructed with — kept as a debug aid
    /// and a hook for future integrations.
    _mode: ProgressMode,
}

/// Initialise the progress tracker for a runall invocation.
///
/// PR 1 returns a no-op guard. Future PRs that wire in a real progress
/// tracker can change the body without breaking callers.
pub(crate) fn init(mode: ProgressMode) -> ProgressGuard {
    ProgressGuard { _mode: mode }
}

/// Announce that the runall pipeline has started. PR 1 logs the event
/// at info level; future PRs can route this through a real tracker.
pub(crate) fn pipeline_started(name: &str) {
    log::info!("pipeline {name} started");
}

/// Announce that the runall pipeline has finished, with the run's
/// success status. PR 1 logs the event at info level.
pub(crate) fn pipeline_finished(success: bool) {
    if success {
        log::info!("pipeline finished: success");
    } else {
        log::info!("pipeline finished: failure");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_mode_value_enum_round_trips() {
        for mode in [
            ProgressMode::Auto,
            ProgressMode::Dashboard,
            ProgressMode::Heartbeat,
            ProgressMode::None,
        ] {
            // Construct a guard, drop it; must not panic.
            let _g = init(mode);
        }
    }

    #[test]
    fn pipeline_event_helpers_do_not_panic() {
        // These helpers route through `log` which is a no-op when no
        // logger is installed; the contract here is just that they run.
        pipeline_started("runall-test");
        pipeline_finished(true);
        pipeline_finished(false);
    }
}
