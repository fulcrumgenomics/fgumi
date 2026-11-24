//! CRAM source steps (feature `cram`): the CRAM counterpart of
//! `ReadBgzfBlocks → BgzfDecompress`.
//!
//! ```text
//! ReadCramSlices    Serial (sticky, reader)  containers/slices off the stream
//!   ↓ CramSlice
//! DecodeCramSlices  Parallel                 slice → framed BAM records
//!   ↓ DecompressedBlock (no BAM header; pair with `FindBamBoundaries::new_no_header`)
//! ```
//!
//! htslib's CRAM reader has this shape internally, with its own thread pool
//! doing the decodes; here the pipeline's workers do them.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use fgumi_htslib::{CramReader, SliceJob};
use fgumi_pipeline_core::{
    BranchOrdering, CounterSpec, HeapSize, HeldRetry, HeldSlot, Ordered, OrderedBytesSingle,
    PhaseCap, QueueSpec, StepCapRef, Unpushed, admit_input, external_wait,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

use crate::types::DecompressedBlock;

pub use fgumi_bam_io::cram_input::CramErrorContext;

/// Estimated decoded-memory allowance shared by all slices from one source.
/// Only the serial source acquires; workers release without blocking.
struct DecodeBudget {
    limit: u64,
    used: AtomicU64,
}

impl DecodeBudget {
    fn try_acquire(self: &Arc<Self>, bytes: u64) -> Option<DecodeLease> {
        let used = self.used.load(Ordering::Acquire);
        // Admit an oversized slice only while idle: otherwise a single long
        // read could never run. Admission in file order avoids reorder deadlock.
        if used != 0 && bytes > self.limit.saturating_sub(used) {
            return None;
        }
        self.used.fetch_add(bytes, Ordering::AcqRel);
        Some(DecodeLease { budget: Arc::clone(self), bytes })
    }
}

/// Holds estimated capacity through decode and any refused output push.
struct DecodeLease {
    budget: Arc<DecodeBudget>,
    bytes: u64,
}

impl Drop for DecodeLease {
    fn drop(&mut self) {
        self.budget.used.fetch_sub(self.bytes, Ordering::Release);
    }
}

/// One CRAM slice awaiting decode.
pub struct CramSlice {
    job: SliceJob,
    lease: DecodeLease,
}

impl HeapSize for CramSlice {
    /// The slice's compressed blocks (its decode output is charged exactly by
    /// `DecompressedBlock`).
    fn heap_size(&self) -> usize {
        usize::try_from(self.job.compressed_bytes()).unwrap_or(usize::MAX)
    }
}

impl Ordered for CramSlice {
    fn ordinal(&self) -> u64 {
        self.job.serial()
    }
}

// ---------------------------------------------------------------------------

/// `Serial + sticky` source: reads CRAM containers and slices in file order.
pub struct ReadCramSlices {
    reader: Option<CramReader>,
    pending: Option<SliceJob>,
    held: HeldSlot<Unpushed<CramSlice>>,
    decode_budget: Arc<DecodeBudget>,
    output_byte_limit: u64,
    context: CramErrorContext,
}

const SLICES: usize = 0;

impl ReadCramSlices {
    /// Use `output_byte_limit` independently for compressed queue bytes and
    /// estimated decoded bytes admitted before dispatch. The decoded allowance
    /// covers queued jobs, running decodes and refused decoded output. Once
    /// output is accepted, the decoded queue's exact byte budget takes over.
    ///
    /// Estimates are conservative metadata-based charges, not a hard RSS bound
    /// on htslib's allocator or reference cache. A slice larger than the whole
    /// allowance is admitted alone, so estimated live decode memory is bounded
    /// by the greater of the allowance and the largest slice estimate.
    #[must_use]
    pub fn new(reader: CramReader, output_byte_limit: u64, context: CramErrorContext) -> Self {
        Self {
            reader: Some(reader),
            pending: None,
            held: HeldSlot::new(),
            decode_budget: Arc::new(DecodeBudget {
                limit: output_byte_limit,
                used: AtomicU64::new(0),
            }),
            output_byte_limit,
            context,
        }
    }
}

impl Step for ReadCramSlices {
    type Input = ();
    type Outputs = OrderedBytesSingle<CramSlice>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReadCramSlices",
            kind: StepKind::Serial,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        Affinity::Reader
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("slices", "slices")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        let job = if let Some(job) = self.pending.take() {
            job
        } else {
            let Some(reader) = self.reader.as_mut() else {
                return Ok(StepOutcome::Finished);
            };
            // A slow pipe writer blocks on the outside world.
            let next =
                external_wait(|| reader.next_slice()).map_err(|e| self.context.wrap_read(&e))?;
            let Some(job) = next else {
                // In-flight slices keep the file alive.
                self.reader = None;
                return Ok(StepOutcome::Finished);
            };
            job
        };
        let estimate = job.decoded_size_estimate().map_err(|e| self.context.wrap_read(&e))?;
        let Some(lease) = self.decode_budget.try_acquire(estimate) else {
            self.pending = Some(job);
            // Capacity returns from decode workers off the queue edge.
            let _ = ctx.pool.request_worker();
            return Ok(StepOutcome::Contention);
        };
        ctx.counters.add(SLICES, 1);
        if let Err(unpushed) = ctx.outputs.push(CramSlice { job, lease }) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

// ---------------------------------------------------------------------------

/// `Parallel` step: decodes one CRAM slice per call into framed BAM records.
pub struct DecodeCramSlices {
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    held_lease: Option<DecodeLease>,
    output_byte_limit: u64,
    context: CramErrorContext,
    cap: Option<Arc<PhaseCap>>,
}

const RECORDS: usize = 0;

impl Clone for DecodeCramSlices {
    fn clone(&self) -> Self {
        Self {
            held: HeldSlot::new(),
            held_lease: None,
            output_byte_limit: self.output_byte_limit,
            context: self.context.clone(),
            cap: self.cap.clone(),
        }
    }
}

impl DecodeCramSlices {
    #[must_use]
    pub fn new(output_byte_limit: u64, context: CramErrorContext) -> Self {
        Self { held: HeldSlot::new(), held_lease: None, output_byte_limit, context, cap: None }
    }

    /// Share a phase admission cap: a sort-first chain's phase-1 cap
    /// (`--sort-threads`), as the BAM and SAM ingest fronts do. `None` (the
    /// default) leaves decoding unrestricted.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        // Accepted output is now charged to its exact byte-bounded queue.
        self.held_lease = None;
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };
        let Some(CramSlice { job, lease }) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let decoded = job.decode().map_err(|e| self.context.wrap_decode(&e))?;
        ctx.counters.add(RECORDS, decoded.records);
        let block = DecompressedBlock { batch_serial: decoded.serial, bytes: decoded.bytes };
        if let Err(unpushed) = ctx.outputs.push(block) {
            self.held.put(unpushed);
            self.held_lease = Some(lease);
        }
        Ok(StepOutcome::Progress)
    }
}

impl Step for DecodeCramSlices {
    type Input = CramSlice;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "DecodeCramSlices",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("records", "records")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn step_cap(&self) -> Option<StepCapRef<'_>> {
        self.cap.as_deref().map(StepCapRef::Fixed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_htslib::{CramWriteOptions, CramWriter};
    use fgumi_pipeline_core::step::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_raw_bam::SamBuilder;

    /// Enough records to cross htslib's default 10,000-record slice boundary.
    fn fixture() -> (tempfile::TempDir, CramReader) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("input.cram");
        let options = CramWriteOptions { no_ref: true, ..CramWriteOptions::default() };
        let mut writer = CramWriter::create(&path, "@HD\tVN:1.6\n", &options).unwrap();
        for i in 0..20_001 {
            let mut builder = SamBuilder::new();
            let record = builder
                .read_name(format!("r{i}").as_bytes())
                .flags(4)
                .sequence(b"ACGT")
                .qualities(&[30; 4])
                .build();
            writer.write_raw(record.as_ref()).unwrap();
            for job in writer.drain_jobs() {
                job.encode().unwrap().write().unwrap();
            }
        }
        let (jobs, finisher) = writer.seal().unwrap();
        for job in jobs {
            job.encode().unwrap().write().unwrap();
        }
        finisher.finish().unwrap();
        let reader = CramReader::open(&path, None).unwrap();
        (dir, reader)
    }

    #[test]
    fn decoded_admission_precedes_dispatch_and_releases_on_slice_drop() {
        let (_dir, reader) = fixture();
        let mut source = ReadCramSlices::new(reader, 1, CramErrorContext::new("fixture", false));
        let mut probe = StepProbe::new(&source);
        let output = probe.take_output::<CramSlice>(0);
        assert_eq!(probe.try_run(&mut source).unwrap(), StepOutcome::Progress);
        let first = output.pop().expect("first oversized slice admitted alone");
        assert_eq!(first.ordinal(), 0);
        assert_eq!(probe.try_run(&mut source).unwrap(), StepOutcome::Contention);
        assert!(output.pop().is_none(), "a second decode must not be admitted");
        drop(first);
        assert_eq!(probe.try_run(&mut source).unwrap(), StepOutcome::Progress);
        assert_eq!(output.pop().expect("next slice admitted after release").ordinal(), 1);
    }

    #[test]
    fn decode_error_releases_admission() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("mapped.cram");
        let reference = dir.path().join("correct.fa");
        let wrong_reference = dir.path().join("wrong.fa");
        std::fs::write(&reference, b">chrT\nACGTACGT\n").unwrap();
        std::fs::write(&wrong_reference, b">chrT\nTTTTTTTT\n").unwrap();
        let options = CramWriteOptions { reference: Some(&reference), ..Default::default() };
        let mut writer =
            CramWriter::create(&path, "@HD\tVN:1.6\n@SQ\tSN:chrT\tLN:8\n", &options).unwrap();
        let record = SamBuilder::new()
            .ref_id(0)
            .pos(0)
            .flags(0)
            .read_name(b"mapped")
            .cigar_ops(&[4 << 4])
            .sequence(b"ACGT")
            .qualities(&[30; 4])
            .build();
        writer.write_raw(record.as_ref()).unwrap();
        let (jobs, finisher) = writer.seal().unwrap();
        for job in jobs {
            job.encode().unwrap().write().unwrap();
        }
        finisher.finish().unwrap();
        let reader = CramReader::open(&path, Some(&wrong_reference)).unwrap();
        let context = CramErrorContext::new("mismatched reference", true);
        let mut source = ReadCramSlices::new(reader, 1, context.clone());
        let mut source_probe = StepProbe::new(&source);
        let slices = source_probe.take_output::<CramSlice>(0);
        let mut decoder = DecodeCramSlices::new(1, context);
        let decode_probe = StepProbe::new(&decoder);
        assert_eq!(source_probe.try_run(&mut source).unwrap(), StepOutcome::Progress);
        decode_probe.push_input(slices.pop().unwrap());
        assert!(source.decode_budget.used.load(Ordering::Acquire) > 0);
        let error = decode_probe.try_run(&mut decoder).unwrap_err();
        assert!(error.to_string().contains("decode"), "{error}");
        assert_eq!(source.decode_budget.used.load(Ordering::Acquire), 0);
    }

    /// Use direct byte backpressure to isolate the worker's refused-output
    /// contract from the ordered queue's separate head-of-line overflow policy.
    struct DirectDecode(DecodeCramSlices);

    impl Step for DirectDecode {
        type Input = CramSlice;
        type Outputs = OrderedBytesSingle<DecompressedBlock>;

        fn profile(&self) -> StepProfile {
            let mut profile = self.0.profile();
            profile.kind = StepKind::Serial;
            profile.branch_ordering = vec![BranchOrdering::None];
            profile
        }

        fn counters(&self) -> &'static [CounterSpec] {
            self.0.counters()
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            self.0.try_run(&mut StepCtx {
                input: ctx.input,
                outputs: ctx.outputs,
                counters: ctx.counters,
                pool: ctx.pool,
                signal: ctx.signal,
            })
        }
    }

    #[rstest::rstest]
    #[case::output_accepted(false)]
    #[case::decoder_abandoned(true)]
    fn decoded_admission_stays_charged_while_output_is_held(#[case] abandon: bool) {
        let (_dir, reader) = fixture();
        let context = CramErrorContext::new("fixture", false);
        let mut source = ReadCramSlices::new(reader, 1, context.clone());
        let mut source_probe = StepProbe::new(&source);
        let slices = source_probe.take_output::<CramSlice>(0);
        let mut decoder = DirectDecode(DecodeCramSlices::new(1, context));
        let mut decode_probe = StepProbe::new(&decoder);
        let blocks = decode_probe.take_output::<DecompressedBlock>(0);
        // The first block fills transport; the second stays with its worker.
        for _ in 0..2 {
            assert_eq!(source_probe.try_run(&mut source).unwrap(), StepOutcome::Progress);
            decode_probe.push_input(slices.pop().unwrap());
            assert_eq!(decode_probe.try_run(&mut decoder).unwrap(), StepOutcome::Progress);
        }
        assert_eq!(source_probe.try_run(&mut source).unwrap(), StepOutcome::Contention);
        if abandon {
            drop(decoder);
        } else {
            assert_eq!(blocks.pop().expect("first decoded block").batch_serial, 0);
            assert_eq!(decode_probe.try_run(&mut decoder).unwrap(), StepOutcome::NoProgress);
            // The retry delivered the held block, exactly once.
            assert_eq!(blocks.pop().expect("held block after retry").batch_serial, 1);
            assert!(blocks.pop().is_none(), "the held block was delivered twice");
        }
        assert_eq!(source_probe.try_run(&mut source).unwrap(), StepOutcome::Progress);
    }
}
