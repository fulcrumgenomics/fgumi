//! Cross-cutting CLI machinery for `fgumi runall`.
//!
//! Each submodule provides one piece of the runall foundation that
//! exists independent of any specific chain:
//!
//! * [`atomic_output`] — `.tmp` → final-output rename guard.
//! * [`signal`]        — process-wide one-shot SIGINT/SIGTERM handler.
//! * [`progress`]      — `ProgressMode` enum and pipeline lifecycle hooks.
//! * [`explain`]       — `--explain` printer (chain-aware).
//! * [`metrics`]       — `--metrics` TSV writer.

pub(crate) mod atomic_output;
pub(crate) mod explain;
pub(crate) mod metrics;
pub(crate) mod progress;
pub(crate) mod signal;
