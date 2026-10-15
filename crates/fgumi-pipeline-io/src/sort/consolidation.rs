//! `ConsolidationRegistry` — the completion side channel from
//! `WriteConsolidatedRun` back to `SpillWrite`.
//!
//! `SpillWrite` starts a consolidation by emitting `ConsolidationStart`; the job
//! runs downstream (`SpillReadPlanner` → `ConsolidateRuns` →
//! `SpillBlockCompressConsolidation` → `WriteConsolidatedRun`) and its end is
//! reported here, off any queue edge. `SpillWrite` takes no input while a job
//! is in flight, so nothing a pool or driver wake would deliver resumes it: it
//! may run on the pinned writer worker (every non-standalone sort) or on the
//! `sort-io` driver, and both idle on a timer park. The registry therefore holds
//! a direct waiter, the pattern of `ArenaPool`'s: `SpillWrite` registers its
//! thread and re-checks before it reports `NoProgress`, and [`complete`] pushes
//! the completion and then unparks the registered thread, so a completion is
//! never missed and never waits on a timer.
//!
//! [`complete`]: ConsolidationRegistry::complete

use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread::Thread;

use parking_lot::Mutex;

/// Hook run between a poll's first empty `take_completed` and its
/// `register_waiter` (tests).
#[cfg(test)]
type RegisterHook = Box<dyn FnOnce(&ConsolidationRegistry) + Send>;

/// Completion side channel with a direct wake (see the module docs).
#[derive(Default)]
pub struct ConsolidationRegistry {
    /// Finished jobs, `(job_id, bytes on disk)`, in completion order.
    done: Mutex<VecDeque<(u32, u64)>>,
    /// The thread to unpark on `complete` (the last one to register).
    waiter: Mutex<Option<Thread>>,
    /// Each started job's output path, until its first block is written.
    outputs: Mutex<HashMap<u32, PathBuf>>,
    #[cfg(test)]
    hook_before_register: Mutex<Option<RegisterHook>>,
}

impl std::fmt::Debug for ConsolidationRegistry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ConsolidationRegistry")
            .field("done", &self.done.lock().len())
            .finish_non_exhaustive()
    }
}

impl ConsolidationRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> Arc<Self> {
        Arc::new(Self::default())
    }

    /// Record a finished job and unpark the registered waiter. The completion
    /// is pushed first, so a waiter that re-checks after registering finds it.
    pub fn complete(&self, job_id: u32, bytes: u64) {
        self.outputs.lock().remove(&job_id);
        self.done.lock().push_back((job_id, bytes));
        if let Some(t) = self.waiter.lock().as_ref() {
            t.unpark();
        }
    }

    /// Pop one completion, if any.
    #[must_use]
    pub fn take_completed(&self) -> Option<(u32, u64)> {
        self.done.lock().pop_front()
    }

    /// Register the calling thread as the one to unpark on `complete`. Cheap
    /// when it is already registered (one `ThreadId` compare).
    pub fn register_waiter(&self) {
        let me = std::thread::current();
        let mut w = self.waiter.lock();
        if w.as_ref().map(Thread::id) != Some(me.id()) {
            *w = Some(me);
        }
    }

    /// Record `job_id`'s output path (`SpillWrite`, at job start).
    pub fn begin(&self, job_id: u32, output: PathBuf) {
        self.outputs.lock().insert(job_id, output);
    }

    /// The output path recorded by [`Self::begin`] for a job not yet
    /// completed (`WriteConsolidatedRun`, on the job's first block).
    #[must_use]
    pub fn output_of(&self, job_id: u32) -> Option<PathBuf> {
        self.outputs.lock().get(&job_id).cloned()
    }

    /// One-shot hook run inside `SpillWrite::poll_registry` between the first
    /// empty `take_completed()` and `register_waiter()` — the lost-wake window.
    #[cfg(test)]
    pub(crate) fn set_hook_before_register_for_test(&self, f: RegisterHook) {
        *self.hook_before_register.lock() = Some(f);
    }

    /// Run (and clear) the hook set by `set_hook_before_register_for_test`.
    #[cfg(test)]
    pub(crate) fn run_hook_before_register_for_test(&self) {
        let hook = self.hook_before_register.lock().take();
        if let Some(f) = hook {
            f(self);
        }
    }
}
