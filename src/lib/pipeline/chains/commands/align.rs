//! Chain builder support for [`crate::pipeline::chains::Stage::Align`].
//!
//! `AlignFinalizeHook` is the post-pipeline hook for the `AlignAndMerge`
//! stage. It logs the records-aligned count and calls
//! `timer.log_completion`.
//!
//! `add_align` (on [`crate::pipeline::chains::builder::ChainBuilder`])
//! registers this hook for every chain that includes `Stage::Align`.

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;

use crate::logging::OperationTimer;
use crate::pipeline::chains::FinalizeHook;

/// Post-pipeline finalize hook for the `AlignAndMerge` stage.
///
/// Logs the number of records that left the aligner and were merged
/// downstream (AAM's `records_emitted` counter), then calls
/// `timer.log_completion` so the per-stage wall-time line appears in
/// the log.
pub(crate) struct AlignFinalizeHook {
    /// Counter atomically incremented by `AlignAndMergeStep` as it
    /// emits `BamTemplateBatch`es downstream.
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for AlignFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let AlignFinalizeHook { records_emitted, timer } = *self;
        let n = records_emitted.load(Ordering::Relaxed);
        info!("AlignAndMerge: pipeline completed successfully ({n} records emitted)");
        timer.log_completion(n);
        Ok(())
    }
}
