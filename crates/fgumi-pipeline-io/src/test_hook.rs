//! A process-wide test hook: one optional callback a step fires at a fixed
//! point, which a test installs to observe or perturb it. nextest runs each
//! test in its own process, so one slot per hook point is enough.

use std::sync::OnceLock;

/// The installed callback.
type Callback<A> = Box<dyn Fn(A) + Send + Sync>;

/// A hook point taking an `A` (declare one as a `static`).
pub(crate) struct TestHook<A: 'static>(OnceLock<parking_lot::Mutex<Option<Callback<A>>>>);

impl<A> TestHook<A> {
    /// An empty hook.
    pub(crate) const fn new() -> Self {
        Self(OnceLock::new())
    }

    fn slot(&self) -> &parking_lot::Mutex<Option<Callback<A>>> {
        self.0.get_or_init(|| parking_lot::Mutex::new(None))
    }

    /// Install `f`, replacing any earlier callback.
    pub(crate) fn set(&self, f: impl Fn(A) + Send + Sync + 'static) {
        *self.slot().lock() = Some(Box::new(f));
    }

    /// Remove the callback.
    pub(crate) fn clear(&self) {
        *self.slot().lock() = None;
    }

    /// Run the callback, if one is installed, with `a`.
    pub(crate) fn fire(&self, a: A) {
        if let Some(f) = self.slot().lock().as_ref() {
            f(a);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::TestHook;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static HOOK: TestHook<usize> = TestHook::new();

    #[test]
    fn fires_only_while_installed() {
        let seen = Arc::new(AtomicUsize::new(0));
        HOOK.fire(1);
        let s = Arc::clone(&seen);
        HOOK.set(move |n| {
            s.fetch_add(n, Ordering::Relaxed);
        });
        HOOK.fire(2);
        HOOK.fire(3);
        HOOK.clear();
        HOOK.fire(4);
        assert_eq!(seen.load(Ordering::Relaxed), 5);
    }
}
