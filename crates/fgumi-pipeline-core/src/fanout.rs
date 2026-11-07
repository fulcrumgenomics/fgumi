//! Fan-out primitives: results of work that pool threads do on behalf of a
//! serial step return to it through an [`OrderedInbox`].
//!
//! # Lost-wake protocol
//!
//! Same argument as `fgumi-sort`'s `MergeDemand`: [`OrderedInbox::await_next`]
//! stores `waiting` and THEN re-checks the slot, both under the inbox mutex;
//! [`OrderedInbox::deliver`] inserts under that mutex and THEN reads
//! `waiting`. The mutex orders the two critical sections, so either the
//! owner's re-check sees the item or the deliverer sees `waiting` and wakes
//! the owner (after releasing the lock).

use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use parking_lot::Mutex;

use crate::item::HeapSize;
use crate::runtime::wake::{SideWake, StepWaker};
use crate::step::SideChannel;

struct State<R> {
    /// Delivered results, by ordinal.
    pending: BTreeMap<u64, R>,
    /// Bytes reserved for ordinals not yet delivered.
    reserved: HashMap<u64, u64>,
    /// Bytes held by each delivered result.
    held: HashMap<u64, u64>,
    /// The next ordinal the owner takes in ordered mode.
    next: u64,
    /// The owner declared (via `await_next`) that it is parked.
    waiting: bool,
}

/// A fan-in channel: worker steps [`deliver`](Self::deliver) results
/// `(ordinal, R)` from pool threads; the owning serial step takes them,
/// in ordinal order when `ordered`, and parks when the next one is absent.
///
/// The inbox itself accepts every reservation and delivery. Its bound is the
/// owner's obligation: the owner keeps at most a window of ordinals posted
/// and not yet taken, and that window bounds what the inbox holds.
///
/// Each ordinal is reserved (optionally) and delivered at most once; a
/// duplicate delivery, a reservation of a delivered ordinal, an ordinal below
/// the next one taken (ordered mode), or the take method of the other mode
/// panics, release builds included.
pub struct OrderedInbox<R> {
    ordered: bool,
    state: Mutex<State<R>>,
    owner: Mutex<Option<StepWaker>>,
    /// Reserved plus delivered-and-untaken bytes; lock-free to read.
    bytes: AtomicU64,
    #[cfg(test)]
    wake_attempts: AtomicU64,
}

impl<R: HeapSize + Send + 'static> OrderedInbox<R> {
    /// A new inbox; `ordered` makes [`try_take_next`](Self::try_take_next) yield
    /// ordinals consecutively from 0.
    #[must_use]
    pub fn new(ordered: bool) -> Arc<Self> {
        Arc::new(Self {
            ordered,
            state: Mutex::new(State {
                pending: BTreeMap::new(),
                reserved: HashMap::new(),
                held: HashMap::new(),
                next: 0,
                waiting: false,
            }),
            owner: Mutex::new(None),
            bytes: AtomicU64::new(0),
            #[cfg(test)]
            wake_attempts: AtomicU64::new(0),
        })
    }

    /// Binds the waker of the owning step, used to wake it on delivery.
    pub fn bind_owner(&self, waker: StepWaker) {
        *self.owner.lock() = Some(waker);
    }

    /// Counts `estimate` bytes against the gauge for an ordinal about to be
    /// worked on. Reserving an ordinal again before its delivery replaces the
    /// estimate.
    ///
    /// # Panics
    ///
    /// If `ordinal` was already delivered and not taken, or, in ordered mode,
    /// is below the next ordinal to take.
    pub fn reserve(&self, ordinal: u64, estimate: u64) {
        let mut s = self.state.lock();
        assert!(
            !self.ordered || ordinal >= s.next,
            "OrderedInbox: reserve for ordinal {ordinal} below next {}",
            s.next
        );
        assert!(
            !s.pending.contains_key(&ordinal),
            "OrderedInbox: reserve for ordinal {ordinal}, which is already delivered"
        );
        if let Some(old) = s.reserved.insert(ordinal, estimate) {
            self.bytes.fetch_sub(old, Ordering::Relaxed);
        }
        self.bytes.fetch_add(estimate, Ordering::Relaxed);
    }

    /// Delivers the result for `ordinal`, replacing its reservation by the
    /// result's heap size, and wakes the owner if it is parked on this
    /// delivery. Returns the wake outcome, or [`SideWake::Unavailable`] when no
    /// wake was attempted or no owner is bound.
    ///
    /// # Panics
    ///
    /// If `ordinal` was already delivered and not taken, or, in ordered mode,
    /// is below the next ordinal to take.
    pub fn deliver(&self, ordinal: u64, item: R) -> SideWake {
        let wake = {
            let mut s = self.state.lock();
            assert!(
                !self.ordered || ordinal >= s.next,
                "OrderedInbox: delivery for ordinal {ordinal} below next {}",
                s.next
            );
            if let Some(est) = s.reserved.remove(&ordinal) {
                self.bytes.fetch_sub(est, Ordering::Relaxed);
            }
            let size = item.heap_size() as u64;
            let was_empty = s.pending.is_empty();
            let duplicate = s.pending.insert(ordinal, item).is_some();
            assert!(!duplicate, "OrderedInbox: duplicate delivery for ordinal {ordinal}");
            s.held.insert(ordinal, size);
            self.bytes.fetch_add(size, Ordering::Relaxed);
            let wake = s.waiting && if self.ordered { ordinal == s.next } else { was_empty };
            if wake {
                s.waiting = false;
            }
            wake
        };
        if !wake {
            return SideWake::Unavailable;
        }
        #[cfg(test)]
        self.wake_attempts.fetch_add(1, Ordering::Relaxed);
        let waker = self.owner.lock().clone();
        waker.map_or(SideWake::Unavailable, |w| w.wake())
    }

    fn take(&self, s: &mut State<R>, ordinal: u64) -> Option<(u64, R)> {
        let item = s.pending.remove(&ordinal)?;
        if let Some(b) = s.held.remove(&ordinal) {
            self.bytes.fetch_sub(b, Ordering::Relaxed);
        }
        s.waiting = false;
        Some((ordinal, item))
    }

    /// Takes the result with the next ordinal, if delivered. Ordered inbox
    /// only: an unordered inbox never advances the internal `next` counter this
    /// reads, so use [`Self::try_take_any`] there.
    ///
    /// # Panics
    ///
    /// If the inbox is unordered.
    pub fn try_take_next(&self) -> Option<(u64, R)> {
        assert!(
            self.ordered,
            "OrderedInbox: try_take_next called on an unordered inbox; use try_take_any"
        );
        let mut s = self.state.lock();
        let n = s.next;
        let r = self.take(&mut s, n)?;
        s.next = n + 1;
        Some(r)
    }

    /// Whether a result is ready to take (the next ordinal in ordered mode,
    /// any in unordered mode). When `false`, records that the owner is
    /// waiting, so the delivery it needs wakes it; the check and the record
    /// are atomic.
    pub fn await_next(&self) -> bool {
        let mut s = self.state.lock();
        let ready =
            if self.ordered { s.pending.contains_key(&s.next) } else { !s.pending.is_empty() };
        s.waiting = !ready;
        ready
    }

    /// Takes the lowest-ordinal delivered result regardless of order.
    /// Unordered inbox only: on an ordered inbox it would remove an ordinal
    /// without advancing `next`, leaving the owner waiting on a result that is
    /// gone (and never woken), so use [`Self::try_take_next`] there.
    ///
    /// # Panics
    ///
    /// If the inbox is ordered.
    pub fn try_take_any(&self) -> Option<(u64, R)> {
        assert!(
            !self.ordered,
            "OrderedInbox: try_take_any called on an ordered inbox; use try_take_next"
        );
        let mut s = self.state.lock();
        let k = *s.pending.keys().next()?;
        self.take(&mut s, k)
    }

    /// Delivered results not yet taken.
    #[must_use]
    pub fn pending(&self) -> usize {
        self.state.lock().pending.len()
    }

    /// Bytes reserved or held (lock-free).
    #[must_use]
    pub fn in_flight_bytes(&self) -> u64 {
        self.bytes.load(Ordering::Relaxed)
    }
}

impl<R: HeapSize + Send + 'static> SideChannel for OrderedInbox<R> {
    fn in_flight_bytes(&self) -> u64 {
        Self::in_flight_bytes(self)
    }

    fn describe(&self) -> String {
        let s = self.state.lock();
        format!("inbox(next={}, pending={})", s.next, s.pending.len())
    }
}

#[cfg(test)]
mod tests;
