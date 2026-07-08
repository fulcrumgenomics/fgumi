//! `GroupBam` mid-step. `Serial + ByItemOrdinal`. Wraps the legacy
//! [`crate::grouper::TemplateGrouper`]: incoming `DecodedRecordBatch`es
//! (records + pre-computed `GroupKey` from the parallel `DecodeRecords`
//! step) are batched by QNAME into `Template`s and emitted as
//! `BamTemplateBatch`es.
//!
//! Templates can span batch boundaries (a template may have records in
//! consecutive `DecodedRecordBatch`es). The Serial mutex is held briefly
//! per call and the per-record `GroupKey` work happens upstream in
//! parallel — matches legacy's "Decode (parallel) → Group (serial)"
//! split (`pipeline/bam.rs:329-403, 2269-2571`).
//!
//! Once the input is drained, `try_run` calls the grouper's `finish()` once
//! (guarded by `finalized`) to flush any final partial template, then reports
//! `Finished`.

use std::io;

use crate::grouper::TemplateGrouper;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BamTemplateBatch, DecodedRecordBatch};
use fgumi_bam_io::Grouper;

/// Max inputs processed per `try_run` invocation. Amortizes the Serial
/// mutex acquisition; matches legacy `bam.rs:2497+` (`MAX_BATCHES_PER_LOCK`).
const MAX_BATCHES_PER_LOCK: usize = 8;

/// `Serial + ByItemOrdinal` template grouper. Holds a [`TemplateGrouper`]
/// and consumes pre-decoded records from upstream `DecodeRecords`.
pub struct GroupBam {
    grouper: TemplateGrouper,
    /// Self-managed output ordinal. Incremented only when a new
    /// `BamTemplateBatch` is emitted. Templates may span several input
    /// batches (a partial template with no `Grouper`-emitted output);
    /// using the input's `batch_serial` would skip ordinals and deadlock
    /// the downstream `ReorderStage`.
    next_output_serial: u64,
    /// Pending output (when push is rejected).
    held: HeldSlot<Unpushed<BamTemplateBatch>>,
    /// Set once the final-flush path has called the (non-idempotent)
    /// `grouper.finish()`; guards against a second call across the
    /// multi-pass completion drain.
    finalized: bool,
    output_byte_limit: u64,
    name: &'static str,
}

impl GroupBam {
    /// Construct a `GroupBam` step.
    ///
    /// * `template_batch_size` — number of templates per emitted batch
    ///   (passed through to `TemplateGrouper::new`).
    /// * `output_byte_limit` — byte budget for the output queue.
    #[must_use]
    pub fn new(template_batch_size: usize, output_byte_limit: u64) -> Self {
        Self {
            grouper: TemplateGrouper::new(template_batch_size),
            next_output_serial: 0,
            held: HeldSlot::new(),
            finalized: false,
            output_byte_limit,
            name: "GroupBam",
        }
    }
}

impl Step for GroupBam {
    type Input = DecodedRecordBatch;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // Process up to `MAX_BATCHES_PER_LOCK` inputs per `try_run` to
        // amortize the Serial mutex acquisition.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(batch) = ctx.input.pop() else { break };
            did_work = true;
            let records = batch.into_records();

            // Feed pre-decoded records into the wrapped grouper. The
            // grouper returns zero or more completed `TemplateBatch`es
            // (= `Vec<Template>`).
            let mut all_templates: Vec<crate::template::Template> = Vec::new();
            for tb in self.grouper.add_records(records)? {
                all_templates.extend(tb);
            }

            if all_templates.is_empty() {
                // No completed templates yet — input absorbed; try next.
                continue;
            }

            let serial = self.next_output_serial;
            self.next_output_serial += 1;
            let out = BamTemplateBatch::new(serial, all_templates);
            match ctx.outputs.push(out) {
                Ok(()) => {}
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // No input this call. If upstream is drained, flush the inner
        // grouper's final partial template once (guarded by `finalized` —
        // `grouper.finish()` is not idempotent), emit it, and report
        // `Finished` once nothing remains. `held` is empty here (step 1
        // returned `Contention` otherwise); a bounced final push is parked
        // in `held` for the held-drain at the top of the next pass.
        if ctx.input.is_drained() {
            if !self.finalized {
                self.finalized = true;
                if let Some(final_batch) = self.grouper.finish()? {
                    if !final_batch.is_empty() {
                        let serial = self.next_output_serial;
                        self.next_output_serial += 1;
                        let out = BamTemplateBatch::new(serial, final_batch);
                        if let Err(unpushed) = ctx.outputs.push(out) {
                            self.held.put(unpushed);
                        }
                        return Ok(StepOutcome::Progress);
                    }
                }
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_serial_byordinal() {
        let s = GroupBam::new(64, 1024);
        let p = s.profile();
        assert_eq!(p.name, "GroupBam");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }
}
