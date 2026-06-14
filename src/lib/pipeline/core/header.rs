//! One-shot lazy SAM/BAM header handle.
//!
//! `HeaderHandle` lets a downstream sink open its output file before the
//! full output header is known, and consume the header at first record
//! time once an upstream step has produced it. The primary motivating
//! consumer is the writer downstream of `AlignAndMergeStep`: the
//! aligner's `@PG` (and any `@RG`/`@CO` lines it adds) are runtime
//! contributions that aren't available at `Pipeline::build` time.
//!
//! ### Contract
//!
//! - `set` and `poison` are one-shot — calling either after the handle
//!   has been resolved returns `AlreadySetError` without altering state.
//! - `try_get` is non-blocking and dispatcher-friendly: a step that
//!   needs the header but observes `None` returns `StepOutcome::NoProgress`
//!   and is rescheduled.
//! - Cloning shares state (via `Arc<OnceLock<...>>`), so multiple
//!   readers see the same set-or-poison outcome. The cross-thread
//!   happens-before guarantee comes from `OnceLock`'s internal
//!   synchronization (`set` is Release, `get` is Acquire); the
//!   `Arc<OnceLock<...>>` wrapper exists only for shareability.
//!
//! ### Producer-side invariant
//!
//! Every consumer that polls `try_get` assumes the handle will
//! eventually resolve via `set` or `poison`. A producer step that
//! exits (cleanly or by panic) **without** resolving its handle leaves
//! consumers polling `None` forever — the pipeline deadlocks.
//! Producer steps owning a `HeaderHandle` must therefore poison
//! their handle in their `Drop` impl as a backstop against panics
//! (or any exit path that bypasses normal completion). The
//! framework does not enforce this; it's a step-level convention.
//! A future revision may split the handle into typed setter / reader
//! halves so an orphaned setter is statically detectable.

use std::fmt;
use std::io;
use std::sync::{Arc, OnceLock};

use noodles::sam::Header;

/// Shared, one-shot header slot.
///
/// See module docs for the full contract.
#[derive(Clone, Default, Debug)]
pub struct HeaderHandle {
    inner: Arc<OnceLock<Result<Header, Arc<io::Error>>>>,
}

impl HeaderHandle {
    /// Construct an empty handle. `try_get` returns `None` until `set`
    /// or `poison` is called.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct a handle that already carries `header`. `try_get`
    /// returns `Some(Ok(&header))` immediately. Convenience for call
    /// sites that have an eager header today and want to opt into the
    /// handle-typed sink API without changing behavior.
    ///
    /// # Panics
    /// The internal `set` call is on a freshly-constructed cell and
    /// is therefore unreachable as a failure — the `expect` is a
    /// belt-and-braces assertion documenting the invariant.
    #[must_use]
    #[allow(clippy::missing_panics_doc)] // documented above; clippy
    // doesn't see the # Panics section across this constructor's
    // delegation to `set`.
    pub fn from_header(header: Header) -> Self {
        let handle = Self::new();
        handle.set(header).expect("fresh HeaderHandle accepts its first set");
        handle
    }

    /// Resolve the handle to `header`. Returns `AlreadySetError` if the
    /// handle was previously set or poisoned.
    ///
    /// # Errors
    /// Returns `AlreadySetError` on the second and subsequent calls.
    pub fn set(&self, header: Header) -> Result<(), AlreadySetError> {
        self.inner.set(Ok(header)).map_err(|_| AlreadySetError)
    }

    /// Resolve the handle to a failure. Subsequent `try_get` calls
    /// surface the error.
    ///
    /// # Errors
    /// Returns `AlreadySetError` on the second and subsequent calls.
    pub fn poison(&self, error: io::Error) -> Result<(), AlreadySetError> {
        self.inner.set(Err(Arc::new(error))).map_err(|_| AlreadySetError)
    }

    /// Non-blocking handle probe.
    ///
    /// Returns:
    /// - `None` if neither `set` nor `poison` has been called yet — the
    ///   caller should yield (e.g. return `StepOutcome::NoProgress`).
    /// - `Some(Ok(&header))` if `set(header)` was called.
    /// - `Some(Err(e))` if `poison(e)` was called. A fresh `io::Error`
    ///   is constructed on every call from the stored kind + message
    ///   (because `io::Error` is not `Clone`). **Note:** the
    ///   original error's `source()` chain and any structured
    ///   payload are not preserved — only kind + display string
    ///   round-trip. Producers that need to surface structured
    ///   diagnostic context should encode it into the display
    ///   string before poisoning the handle.
    #[must_use]
    pub fn try_get(&self) -> Option<io::Result<&Header>> {
        self.inner.get().map(|stored| match stored {
            Ok(header) => Ok(header),
            Err(err) => Err(io::Error::new(err.kind(), err.to_string())),
        })
    }

    /// `true` once `set` or `poison` has resolved this handle.
    #[must_use]
    pub fn is_set(&self) -> bool {
        self.inner.get().is_some()
    }
}

/// Returned by `set` / `poison` when the handle was previously
/// resolved.
#[derive(Debug, Clone, Copy)]
pub struct AlreadySetError;

impl fmt::Display for AlreadySetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "HeaderHandle has already been set or poisoned")
    }
}

impl std::error::Error for AlreadySetError {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;

    #[test]
    fn new_handle_has_no_value() {
        let h = HeaderHandle::new();
        assert!(!h.is_set());
        assert!(h.try_get().is_none());
    }

    #[test]
    fn from_header_resolves_immediately() {
        let h = HeaderHandle::from_header(Header::default());
        assert!(h.is_set());
        let got = h.try_get().expect("set").expect("ok");
        assert_eq!(got, &Header::default());
    }

    #[test]
    fn set_then_try_get_returns_ok() {
        let h = HeaderHandle::new();
        h.set(Header::default()).expect("first set");
        assert!(h.is_set());
        let got = h.try_get().expect("set").expect("ok");
        assert_eq!(got, &Header::default());
    }

    #[test]
    fn poison_then_try_get_returns_err() {
        let h = HeaderHandle::new();
        h.poison(io::Error::other("aligner crashed")).expect("first poison");
        assert!(h.is_set());
        let err = h.try_get().expect("set").expect_err("poisoned");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert_eq!(err.to_string(), "aligner crashed");
    }

    #[test]
    fn second_set_returns_already_set() {
        let h = HeaderHandle::new();
        h.set(Header::default()).expect("first");
        let err = h.set(Header::default()).expect_err("second set");
        let _ = err; // type asserted by the binding
    }

    #[test]
    fn set_then_poison_returns_already_set() {
        let h = HeaderHandle::new();
        h.set(Header::default()).expect("first");
        let _err: AlreadySetError =
            h.poison(io::Error::other("late")).expect_err("poison after set");
    }

    #[test]
    fn poison_then_set_returns_already_set() {
        let h = HeaderHandle::new();
        h.poison(io::Error::other("first")).expect("first");
        let _err: AlreadySetError = h.set(Header::default()).expect_err("set after poison");
    }

    #[test]
    fn clones_share_state() {
        let a = HeaderHandle::new();
        let b = a.clone();
        assert!(!a.is_set() && !b.is_set());
        a.set(Header::default()).expect("first");
        assert!(a.is_set() && b.is_set());
        assert!(b.try_get().expect("set").is_ok());
    }

    #[test]
    fn poison_preserves_error_kind_and_message_across_calls() {
        let h = HeaderHandle::new();
        h.poison(io::Error::new(io::ErrorKind::BrokenPipe, "stderr ring: foo")).unwrap();
        for _ in 0..3 {
            let e = h.try_get().unwrap().unwrap_err();
            assert_eq!(e.kind(), io::ErrorKind::BrokenPipe);
            assert_eq!(e.to_string(), "stderr ring: foo");
        }
    }

    #[test]
    fn set_from_another_thread_is_observable() {
        let h = HeaderHandle::new();
        let h2 = h.clone();
        let join = thread::spawn(move || {
            h2.set(Header::default()).expect("first set");
        });
        join.join().expect("thread join");
        assert!(h.is_set());
        assert!(h.try_get().unwrap().is_ok());
    }

    #[test]
    fn from_header_then_set_returns_already_set() {
        let h = HeaderHandle::from_header(Header::default());
        let _err: AlreadySetError =
            h.set(Header::default()).expect_err("from_header consumes the slot");
    }

    #[test]
    fn from_header_then_poison_returns_already_set() {
        let h = HeaderHandle::from_header(Header::default());
        let _err: AlreadySetError =
            h.poison(io::Error::other("late")).expect_err("from_header consumes the slot");
    }

    #[test]
    fn poison_visible_to_all_clones() {
        let a = HeaderHandle::new();
        let b = a.clone();
        let c = a.clone();
        a.poison(io::Error::new(io::ErrorKind::BrokenPipe, "boom")).unwrap();
        for clone in [&a, &b, &c] {
            let err = clone.try_get().expect("set").expect_err("poisoned");
            assert_eq!(err.kind(), io::ErrorKind::BrokenPipe);
            assert_eq!(err.to_string(), "boom");
        }
    }
}
