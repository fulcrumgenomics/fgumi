//! Integration tests for the `fgumi runall` foundation (PR 1).
//!
//! These exercise the dispatch layer, atomic-output cleanup, and
//! `--explain` short-circuit through the actual `fgumi` binary —
//! the in-process unit tests cover the same code paths but the
//! integration tests guard against CLI-parsing regressions and
//! verify the binary's user-visible behaviour.

mod test_runall_dispatch_errors;
mod test_runall_explain;
