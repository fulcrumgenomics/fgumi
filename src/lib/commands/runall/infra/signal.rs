//! Process-wide SIGINT / SIGTERM handler for `fgumi runall`.
//!
//! `ctrlc::set_handler` errors if installed more than once per process,
//! so runall guards installation behind a [`std::sync::Once`]. Each
//! `runall` invocation gets back a clone of the same shared
//! `Arc<AtomicBool>` (the cancel flag); the SIGINT/SIGTERM handler stores
//! `true` into that flag, and the running pipeline polls the flag to shut
//! down cleanly.
//!
//! Re-running runall in the same process (e.g. an integration test that
//! invokes the command in-process several times) is safe: the second
//! call observes the already-installed handler, clears any stale `true`
//! left over from the previous run, and returns the same flag.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result};

/// Once-guard for installing the process-wide signal handler.
static SIGNAL_HANDLER_INIT: std::sync::Once = std::sync::Once::new();
/// The shared cancel flag. Stored in a `OnceLock` so the closure
/// installed by `ctrlc::set_handler` and subsequent runall calls all
/// see the same `Arc<AtomicBool>`.
static GLOBAL_CANCEL_FLAG: std::sync::OnceLock<Arc<AtomicBool>> = std::sync::OnceLock::new();

/// Install (idempotently) a SIGINT/SIGTERM handler that flips the
/// returned cancel flag to `true`.
///
/// On first call, registers the OS-level signal handler via `ctrlc`.
/// On subsequent calls, returns a clone of the same flag without trying
/// to re-register (which would error). Always resets the flag to
/// `false` before returning so a stale `true` from a previous run does
/// not insta-cancel the new pipeline.
///
/// # Errors
///
/// Returns the underlying [`ctrlc::Error`] (wrapped via `anyhow::Context`)
/// only if registration failed on the first call.
pub(crate) fn install_signal_handler_once() -> Result<Arc<AtomicBool>> {
    let flag = GLOBAL_CANCEL_FLAG.get_or_init(|| Arc::new(AtomicBool::new(false))).clone();
    let mut install_err: Option<ctrlc::Error> = None;
    SIGNAL_HANDLER_INIT.call_once(|| {
        let handler_flag = flag.clone();
        if let Err(e) = ctrlc::set_handler(move || {
            log::warn!("received interrupt signal; cancelling pipeline");
            handler_flag.store(true, Ordering::SeqCst);
        }) {
            install_err = Some(e);
        }
    });
    if let Some(e) = install_err {
        return Err(e).context("Failed to set SIGINT/SIGTERM handler");
    }
    // Clear any stale cancel state from a prior invocation in the same
    // process (only relevant in tests that drive runall multiple times).
    flag.store(false, Ordering::SeqCst);
    Ok(flag)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two back-to-back installations must succeed and return flags that
    /// alias the same underlying atomic. Without the `Once` guard, the
    /// second call would error with "handler already set".
    #[test]
    fn install_signal_handler_is_idempotent() {
        let first = install_signal_handler_once().expect("first install must succeed");
        let second = install_signal_handler_once().expect("second install must be a no-op");

        // Mutate via one Arc, observe via the other → same atomic.
        first.store(true, Ordering::SeqCst);
        assert!(second.load(Ordering::SeqCst), "both flags must alias the same atomic");

        // The next call should reset the flag.
        let third = install_signal_handler_once().expect("third install must be a no-op");
        assert!(!third.load(Ordering::SeqCst), "stale cancel state must be cleared");
    }
}
