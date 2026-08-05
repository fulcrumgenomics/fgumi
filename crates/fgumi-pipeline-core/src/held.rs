//! Generic held-item slot for non-blocking back-pressure.
//!
//! When a step can't push to a downstream queue (it's full), it stashes the
//! item in a `HeldSlot` and returns `StepOutcome::Progress`. The next
//! `try_run` call drains the held item before doing new work.
//!
//! Draining the held item before new work is a *step-author convention*, not
//! something `HeldSlot` enforces: the type only provides single-slot put/take
//! (with a double-put panic). The step is responsible for checking and draining
//! the slot first on each `try_run`.

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
    #[should_panic(expected = "HeldSlot already occupied")]
    fn double_put_panics() {
        // Pin the message, not just the unwind: `catch_unwind` accepts ANY panic,
        // so a `put` that started failing for an unrelated reason would still pass
        // while the double-put guard itself had silently gone. Must fire in debug
        // AND release — overwriting a held item is silent record loss.
        let mut s = HeldSlot::new();
        s.put(1_u32);
        s.put(2_u32);
    }
}
