//! Generic held-item slot for non-blocking back-pressure.
//!
//! When a step can't push to a downstream queue (it's full), it stashes the
//! item in a `HeldSlot` and returns `StepOutcome::Progress`. The next
//! `try_run` call drains the held item before doing new work.

pub struct HeldSlot<T> {
    inner: Option<T>,
}

impl<T> Default for HeldSlot<T> {
    fn default() -> Self {
        Self { inner: None }
    }
}

impl<T> HeldSlot<T> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn is_held(&self) -> bool {
        self.inner.is_some()
    }

    /// Stash an item.
    ///
    /// # Panics
    ///
    /// Panics if the slot is already occupied (a contract violation: callers
    /// must drain via `take` before calling `put` again). This is a hard
    /// `assert!` rather than `debug_assert!` because silently overwriting would
    /// **drop the previously held item** — a record lost from the pipeline —
    /// which must not pass unnoticed in release builds. The check is a single
    /// predictable branch on the back-pressure path, so the cost is negligible.
    pub fn put(&mut self, item: T) {
        assert!(self.inner.is_none(), "HeldSlot already occupied");
        self.inner = Some(item);
    }

    pub fn take(&mut self) -> Option<T> {
        self.inner.take()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_slot_is_empty() {
        let s: HeldSlot<u32> = HeldSlot::new();
        assert!(!s.is_held());
    }

    #[test]
    fn put_then_take_round_trips() {
        let mut s = HeldSlot::new();
        s.put(42_u32);
        assert!(s.is_held());
        assert_eq!(s.take(), Some(42));
        assert!(!s.is_held());
    }

    #[test]
    fn take_from_empty_returns_none() {
        let mut s: HeldSlot<u32> = HeldSlot::new();
        assert_eq!(s.take(), None);
    }

    #[test]
    fn double_put_panics() {
        let mut s = HeldSlot::new();
        s.put(1_u32);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            s.put(2_u32);
        }));
        assert!(
            result.is_err(),
            "double put must panic (in debug AND release) to avoid silent record loss"
        );
    }
}
