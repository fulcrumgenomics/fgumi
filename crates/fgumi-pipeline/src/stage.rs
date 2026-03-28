//! Core traits and types for the pipeline's work-stealing scheduler.

use anyhow::Result;

/// Trait for Zone 3 pipeline stages that process work items.
///
/// The scheduler pulls items from input queues, calls `process()`, and pushes
/// results to output queues. Any worker thread can execute any stage.
pub trait PipelineStage: Send + Sync {
    /// The input type consumed by this stage.
    type Input: Send;
    /// The output type produced by this stage.
    type Output: Send;

    /// Process a single work item. Called by worker threads.
    ///
    /// # Errors
    ///
    /// Returns an error if processing fails.
    fn process(&self, input: Self::Input) -> Result<Self::Output>;

    /// Estimate memory usage of one output item (bytes).
    /// Used by the scheduler for queue backpressure.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize;
}

/// Execute a closure with a thread-local value, initializing it on first use.
///
/// This is a helper for the common pattern of lazily initializing a per-worker-thread
/// resource stored in `thread_local! { RefCell<Option<T>> }`. On the first call in
/// each thread, `init` creates the value; subsequent calls reuse the stored instance.
///
/// # Panics
///
/// Cannot panic in practice: the `unwrap()` on the inner `Option` is reached only
/// after the `is_none()` check has ensured the value is initialized.
pub fn with_thread_local<T, R>(
    local: &'static std::thread::LocalKey<std::cell::RefCell<Option<T>>>,
    init: impl FnOnce() -> T,
    f: impl FnOnce(&mut T) -> R,
) -> R {
    local.with(|cell| {
        let mut borrow = cell.borrow_mut();
        if borrow.is_none() {
            *borrow = Some(init());
        }
        f(borrow.as_mut().unwrap())
    })
}

/// A work item tagged with a sequence number for ordered reassembly.
#[derive(Debug)]
pub struct SequencedItem<T> {
    /// Sequence number for ordering.
    pub seq: u64,
    /// The wrapped item.
    pub item: T,
    /// Estimated heap memory in bytes for the item, computed at push time.
    /// Used by the scheduler to release the correct amount from the input queue.
    pub memory_estimate: usize,
}

impl<T> SequencedItem<T> {
    /// Create a new sequenced item with a memory estimate.
    pub fn new(seq: u64, item: T, memory_estimate: usize) -> Self {
        Self { seq, item, memory_estimate }
    }
}
