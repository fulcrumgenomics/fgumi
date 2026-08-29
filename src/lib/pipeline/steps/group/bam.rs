//! `GroupBam` mid-step. `Serial + ByItemOrdinal`. Wraps the legacy
//! [`crate::grouper::TemplateGrouper`]: incoming `DecodedRecordBatch`es
//! (records + pre-computed `GroupKey` from the parallel `DecodeRecords`
//! step) are batched by QNAME into `Template`s and emitted as
//! `BamTemplateBatch`es.
//!
//! Templates can span batch boundaries (a template may have records in
//! consecutive `DecodedRecordBatch`es). The `Serial` mutex is held briefly
//! per call and the per-record `GroupKey` work happens upstream in
//! parallel — matches legacy's "Decode (parallel) → Group (serial)"
//! split (`pipeline/bam.rs:329-403, 2269-2571`).
//!
//! Once the input is drained, `try_run` calls the grouper's `finish()` once
//! (guarded by `finalized`) to flush any final partial template, emits whatever
//! batches remain queued, and only then reports `Finished`.
//!
//! ## One emitted batch per formed batch
//!
//! [`TemplateGrouper::add_records`] already returns *correctly sized* batches —
//! its `collect_batches` peels off exactly `batch_size` templates at a time and
//! leaves the remainder pending — so this step must not flatten them back
//! together. An earlier version did, `extend`ing every returned batch into one
//! output, which turned `template_batch_size` from a bound into a mere trigger:
//! the emitted batch was sized by however many templates one input batch
//! happened to complete, measured at 624x the configured size for a 40k-record
//! input. That also let a single item blow the output queue's byte budget wide
//! open, since [`crate::pipeline::core::queues::ByteBoundedQueue`] admits any
//! item once it is under budget (a strict `cur + size <= limit` would deadlock
//! producers whose items legitimately exceed the limit).
//!
//! Batches the output queue has no room for are parked in `pending_batches` and
//! re-offered on later calls, which is why `Finished` is gated on that queue
//! being empty as well as on the input being drained.

use std::collections::VecDeque;
use std::io;

use crate::grouper::TemplateGrouper;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BamTemplateBatch, DecodedRecordBatch};
use crate::template::Template;
use fgumi_bam_io::Grouper;

/// Max inputs processed per `try_run` invocation. Amortizes the `Serial`
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
    /// Batches the grouper has already formed but the output queue has not yet
    /// accepted, in emission order.
    ///
    /// Distinct from `held`, which parks the single batch whose push was
    /// *rejected*: one `add_records` call can complete many batches, and only
    /// one of them can be in flight at a time, so the rest queue here rather
    /// than being merged into an oversized batch.
    pending_batches: VecDeque<Vec<Template>>,
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
            pending_batches: VecDeque::new(),
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

        // Flush batches formed by an earlier call before forming any more, so
        // `pending_batches` drains rather than growing. If the output is still
        // refusing, stop here — taking more input would only deepen the queue.
        if self.emit_pending(ctx) {
            return Ok(StepOutcome::Progress);
        }

        // Process up to `MAX_BATCHES_PER_LOCK` inputs per `try_run` to
        // amortize the Serial mutex acquisition.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(batch) = ctx.input.pop() else { break };
            did_work = true;
            let records = batch.into_records();

            // Feed pre-decoded records into the wrapped grouper. The grouper
            // returns zero or more completed `TemplateBatch`es (= `Vec<Template>`),
            // each already capped at `template_batch_size`. Queue them as-is —
            // merging them here would undo that cap.
            self.pending_batches.extend(self.grouper.add_records(records)?);
        }

        if did_work {
            self.emit_pending(ctx);
            return Ok(StepOutcome::Progress);
        }

        // No input this call. If upstream is drained, flush the inner
        // grouper's final partial template once (guarded by `finalized` —
        // `grouper.finish()` is not idempotent), then drain whatever is still
        // queued. `held` is empty here (step 1 returned `Contention`
        // otherwise); a bounced push is parked in `held` for the held-drain at
        // the top of the next pass.
        if ctx.input.is_drained() {
            if !self.finalized {
                self.finalized = true;
                if let Some(final_batch) = self.grouper.finish()?
                    && !final_batch.is_empty()
                {
                    self.pending_batches.push_back(final_batch);
                }
            }
            // `Finished` must wait for the queue to empty: reporting it while a
            // batch is still queued closes the output with records unsent,
            // which surfaces as short output far from here rather than as an
            // error at the point of loss.
            if !self.pending_batches.is_empty() {
                self.emit_pending(ctx);
                return Ok(StepOutcome::Progress);
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

impl GroupBam {
    /// Push queued batches until the output rejects one or the queue empties.
    ///
    /// Returns `true` if a push was rejected — the batch is parked in `held`
    /// and the caller should yield rather than take on more work. Ordinals are
    /// consumed in emission order, including for the rejected batch, because
    /// `held` retries that same batch with the ordinal it was assigned; the
    /// `ByItemOrdinal` reorder stage downstream releases only contiguously, so
    /// a skipped or reused ordinal would stall it.
    fn emit_pending(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        while let Some(templates) = self.pending_batches.pop_front() {
            let serial = self.next_output_serial;
            self.next_output_serial += 1;
            let out = BamTemplateBatch::new(serial, templates);
            if let Err(unpushed) = ctx.outputs.push(out) {
                self.held.put(unpushed);
                return true;
            }
        }
        false
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
