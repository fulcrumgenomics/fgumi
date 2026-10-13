//! `RecordStreamSource` — the chain source fed by an in-process record
//! generator (`fgumi simulate`) over a bounded channel.
//!
//! A `Serial` step on the reader worker: each `try_run` takes one batch of raw
//! records from the channel and pushes it as a [`RecordBatch`] with the next
//! serial. When the channel is momentarily empty it waits up to a
//! millisecond, so the reader worker does not spin, and a batch that arrives in
//! that wait is processed, not dropped. A disconnected channel (the generator
//! finished and dropped its sender) finishes the step; a batch carrying the
//! generator's error fails the pipeline with that error.

use std::io;
use std::time::Duration;

use crossbeam_channel::{Receiver, RecvTimeoutError, TryRecvError};
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};
use fgumi_raw_bam::RawRecord;

use crate::types::RecordBatch;

/// What one channel message carries: a batch of records, or the generator's
/// terminal error.
pub type RecordStreamItem = anyhow::Result<Vec<RawRecord>>;

/// How long an empty poll waits for the generator before reporting no progress.
const RECV_WAIT: Duration = Duration::from_millis(1);

/// The record-generator source (see the module docs).
pub struct RecordStreamSource {
    rx: Receiver<RecordStreamItem>,
    next_serial: u64,
    held: HeldSlot<Unpushed<RecordBatch>>,
    output_byte_limit: u64,
    /// Run after `try_recv` saw an empty channel, before the timed wait
    /// (tests: place a batch in that window deterministically).
    #[cfg(test)]
    before_recv_timeout: Option<Box<dyn FnMut() + Send>>,
}

impl RecordStreamSource {
    /// A source reading `rx`, bounding its output edge at `output_byte_limit`
    /// bytes.
    #[must_use]
    pub fn new(rx: Receiver<RecordStreamItem>, output_byte_limit: u64) -> Self {
        Self {
            rx,
            next_serial: 0,
            held: HeldSlot::new(),
            output_byte_limit,
            #[cfg(test)]
            before_recv_timeout: None,
        }
    }

    /// Install the hook run between the empty `try_recv` and the timed wait.
    #[cfg(test)]
    pub(crate) fn set_before_recv_timeout_for_test(&mut self, f: Box<dyn FnMut() + Send>) {
        self.before_recv_timeout = Some(f);
    }

    /// Turn one channel message into a pushed batch, or the step's error.
    fn take(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        item: RecordStreamItem,
    ) -> io::Result<StepOutcome> {
        let records = item.map_err(|e| io::Error::other(format!("{e:#}")))?;
        let batch = RecordBatch::new(self.next_serial, &records);
        self.next_serial += 1;
        if let Err(u) = ctx.outputs.push(batch) {
            self.held.put(u);
        }
        Ok(StepOutcome::Progress)
    }
}

impl Step for RecordStreamSource {
    type Input = ();
    type Outputs = OrderedBytesSingle<RecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "RecordStreamSource",
            kind: StepKind::Serial,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        Affinity::Reader
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        match self.rx.try_recv() {
            Ok(item) => self.take(ctx, item),
            Err(TryRecvError::Disconnected) => Ok(StepOutcome::Finished),
            Err(TryRecvError::Empty) => {
                #[cfg(test)]
                if let Some(f) = &mut self.before_recv_timeout {
                    f();
                }
                match self.rx.recv_timeout(RECV_WAIT) {
                    Ok(item) => self.take(ctx, item),
                    Err(RecvTimeoutError::Timeout) => Ok(StepOutcome::NoProgress),
                    Err(RecvTimeoutError::Disconnected) => Ok(StepOutcome::Finished),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_pipeline_core::handles::BranchInputHandle;
    use fgumi_pipeline_core::step::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_raw_bam::SamBuilder;

    fn rec(i: u32) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(format!("r{i}").as_bytes()).sequence(b"ACGT").qualities(&[30; 4]);
        b.build()
    }

    fn drive(
        src: &RecordStreamSource,
    ) -> (StepProbe<RecordStreamSource>, BranchInputHandle<RecordBatch>) {
        let mut p = StepProbe::new(src);
        let out = p.take_output::<RecordBatch>(0);
        (p, out)
    }

    /// Catches: `Finished` before the last batch is pushed, or never on disconnect.
    #[test]
    fn record_stream_source_finishes_on_disconnect() {
        let (tx, rx) = crossbeam_channel::bounded(4);
        let mut src = RecordStreamSource::new(rx, 1 << 20);
        let (p, out) = drive(&src);
        tx.send(Ok(vec![rec(1), rec(2)])).unwrap();
        drop(tx);
        assert_eq!(p.try_run(&mut src).unwrap(), StepOutcome::Progress);
        assert_eq!(out.pop().unwrap().len(), 2);
        assert_eq!(p.try_run(&mut src).unwrap(), StepOutcome::Finished);
    }

    /// Review I13 (r1): a batch returned by `recv_timeout` is processed, not
    /// dropped. The hook runs after `try_recv` saw `Empty`, so the batch
    /// deterministically arrives via `recv_timeout`. Catches: `Ok(_) =>
    /// NoProgress` on that arm (a lost batch).
    #[test]
    fn record_stream_source_keeps_the_batch_recv_timeout_returns() {
        let (tx, rx) = crossbeam_channel::bounded(4);
        let mut src = RecordStreamSource::new(rx, 1 << 20);
        let tx2 = tx.clone();
        src.set_before_recv_timeout_for_test(Box::new(move || tx2.send(Ok(vec![rec(7)])).unwrap()));
        let (p, out) = drive(&src);
        assert_eq!(p.try_run(&mut src).unwrap(), StepOutcome::Progress);
        assert_eq!(out.pop().expect("the recv_timeout batch was pushed").len(), 1);
        drop(tx);
    }

    /// A generator error becomes the step's error, with the message.
    #[test]
    fn record_stream_source_fails_on_a_generator_error() {
        let (tx, rx) = crossbeam_channel::bounded(4);
        let mut src = RecordStreamSource::new(rx, 1 << 20);
        let (p, _out) = drive(&src);
        tx.send(Err(anyhow::anyhow!("generator broke"))).unwrap();
        let err = p.try_run(&mut src).unwrap_err();
        assert!(err.to_string().contains("generator broke"), "{err}");
    }
}
