//! Observer for thread names spawned during a test. Used by gate tests
//! to assert that the pipeline doesn't spawn per-stage threads.
//!
//! Delegates to `fgumi_lib::runall::engine::test_support`, which is what the
//! pool (and the pipeline source/sink helpers) feed during their startup.
//! The registry is an append-only log cleared by `clear()`: each pipeline
//! run leaves behind the thread names it spawned, even after the threads
//! have exited — that lets the gate tests assert on the set of names that
//! participated in the run without needing to observe them while alive.

#![allow(dead_code)]

use fgumi_lib::runall::engine::test_support;

/// Register a thread name. Not typically called from tests — the pool does
/// this internally via `pipeline::test_support::register_thread`.
pub fn register(name: &str) {
    test_support::register_thread(name);
}

/// Unregister a thread name. Not typically called from tests.
pub fn unregister(name: &str) {
    test_support::unregister_thread(name);
}

/// Snapshot the current list of registered thread names.
pub fn snapshot() -> Vec<String> {
    test_support::snapshot_threads()
}

/// Clear the registry (call at the start of each gate test so its
/// assertion only observes threads spawned by that test).
pub fn clear() {
    test_support::clear_threads();
}
