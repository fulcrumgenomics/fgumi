//! Cancellation token: propagates SIGINT/SIGTERM to pipeline workers.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// Shared cancellation flag.
///
/// The binary crate installs a SIGINT/SIGTERM handler that calls
/// [`CancelToken::cancel`]. Workers check [`CancelToken::is_cancelled`] in
/// their loop and exit promptly if set.
#[derive(Clone, Default)]
pub struct CancelToken {
    flag: Arc<AtomicBool>,
}

impl CancelToken {
    #[must_use]
    pub fn new() -> Self {
        Self { flag: Arc::new(AtomicBool::new(false)) }
    }

    pub fn from_arc(flag: Arc<AtomicBool>) -> Self {
        Self { flag }
    }

    pub fn cancel(&self) {
        self.flag.store(true, Ordering::SeqCst);
    }

    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.flag.load(Ordering::Relaxed)
    }

    /// Shared raw flag (for interop with existing ctrlc handlers).
    #[must_use]
    pub fn inner(&self) -> &Arc<AtomicBool> {
        &self.flag
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initial_not_cancelled() {
        let tok = CancelToken::new();
        assert!(!tok.is_cancelled());
    }

    #[test]
    fn test_cancel_sets_flag() {
        let tok = CancelToken::new();
        tok.cancel();
        assert!(tok.is_cancelled());
    }

    #[test]
    fn test_clone_shares_flag() {
        let tok1 = CancelToken::new();
        let tok2 = tok1.clone();
        tok1.cancel();
        assert!(tok2.is_cancelled());
    }

    #[test]
    fn test_from_arc() {
        let flag = Arc::new(AtomicBool::new(false));
        let tok = CancelToken::from_arc(flag.clone());
        flag.store(true, Ordering::SeqCst);
        assert!(tok.is_cancelled());
    }
}
