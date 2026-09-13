//! `FindBamBoundaries` mid-step. `Serial + ByItemOrdinal`. Strips the BAM
//! header (on first input) and finds record boundaries within decompressed
//! BGZF block data, emitting `DecompressedBlock`s whose bytes contain only
//! complete records (4-byte `block_size` prefix + record body, repeated).
//!
//! Wraps [`super::state::BoundaryState`] — the boundary-finding state
//! machine handles header skipping, cross-block record carryover, and
//! validation. Reusing it ensures the new framework's boundary semantics
//! match the legacy pipeline's exactly.

use std::io;

use super::state::BoundaryState;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{CounterSpec, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::DecompressedBlock;

/// Max inputs processed per `try_run` invocation. Amortizes the Serial
/// mutex acquisition; matches legacy `bam.rs:2050+` (`MAX_BATCHES_PER_LOCK`).
const MAX_BATCHES_PER_LOCK: usize = 8;

/// Counter slot index: record boundaries found this call.
const RECORDS: usize = 0;

/// Number of records represented by a `BoundaryBatch`'s `offsets`: the vec
/// holds `num_records + 1` entries (a trailing `buffer.len()` sentinel), or is
/// empty when the batch itself carries no offsets at all.
fn record_count(offsets: &[usize]) -> u64 {
    offsets.len().saturating_sub(1) as u64
}

/// `Serial + ByItemOrdinal` boundary finder. Holds `BoundaryState` (which
/// owns the cross-block carryover buffer + header-skip flag).
pub struct FindBamBoundaries {
    state: BoundaryState,
    /// Pending output batch when we found boundaries but the push was
    /// rejected. Held across retries until pushed.
    held: HeldSlot<Unpushed<DecompressedBlock>>,
    /// Self-managed output ordinal. Incremented only when a new
    /// `DecompressedBlock` is emitted. Using the input's `batch_serial`
    /// directly would skip ordinals when the boundary state absorbs an
    /// input without producing output (header bytes, mid-record carryover);
    /// the downstream `ReorderStage` requires consecutive `0, 1, 2, …`.
    next_output_serial: u64,
    /// Set once the final-flush path has called the (non-idempotent)
    /// `state.finish()`; guards against a second call across the multi-pass
    /// completion drain.
    finalized: bool,
    output_byte_limit: u64,
}

impl FindBamBoundaries {
    /// Construct expecting the first input batch to begin with the BAM
    /// header (matches `BAM_MAGIC`). The header bytes are skipped on the
    /// first call; subsequent calls just find record boundaries.
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self {
            state: BoundaryState::new(),
            held: HeldSlot::new(),
            next_output_serial: 0,
            finalized: false,
            output_byte_limit,
        }
    }

    /// Construct expecting the input stream to be already past the BAM
    /// header. Useful for runall-spliced sub-pipelines.
    #[must_use]
    pub fn new_no_header(output_byte_limit: u64) -> Self {
        Self {
            state: BoundaryState::new_no_header(),
            held: HeldSlot::new(),
            next_output_serial: 0,
            finalized: false,
            output_byte_limit,
        }
    }
}

impl Step for FindBamBoundaries {
    type Input = DecompressedBlock;
    type Outputs = OrderedBytesSingle<DecompressedBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "FindBamBoundaries",
            kind: StepKind::Serial,
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
        // Batch total for this call, accumulated across every input processed
        // in the loop (and the final-flush path) below; bumped exactly once —
        // right before returning — regardless of which of `try_run_inner`'s
        // several exit paths fires.
        let mut records_this_call: u64 = 0;
        let outcome = self.try_run_inner(ctx, &mut records_this_call);
        ctx.counters.add(RECORDS, records_this_call);
        outcome
    }
}

impl FindBamBoundaries {
    /// Body of `Step::try_run`, factored out so the counter bump in the trait
    /// method stays a single call regardless of which early-return path below
    /// fires. `records_this_call` accumulates the batch total (record
    /// boundaries found) for this call; the caller bumps `ctx.counters` with
    /// the final tally after this returns.
    fn try_run_inner(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        records_this_call: &mut u64,
    ) -> io::Result<StepOutcome> {
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
        // amortize the Serial mutex acquisition. Each iteration pops one
        // input, runs find_boundaries, and pushes at most one output —
        // header-only or carryover-only inputs are fully absorbed and produce
        // no output (the loop `continue`s, and `did_work` still records the
        // consumed input as Progress). We stop early on push rejection (held
        // the unpushed; subsequent iterations would contend on backpressure)
        // or on input exhaustion.
        let mut did_work = false;
        for _ in 0..MAX_BATCHES_PER_LOCK {
            let Some(block) = ctx.input.pop() else { break };
            did_work = true;

            let boundary_batch = self.state.find_boundaries(&block.bytes)?;
            *records_this_call += record_count(&boundary_batch.offsets);
            if boundary_batch.buffer.is_empty() {
                // Input fully absorbed (header/leftover); try next input.
                continue;
            }

            let serial = self.next_output_serial;
            self.next_output_serial += 1;
            let out = DecompressedBlock { batch_serial: serial, bytes: boundary_batch.buffer };
            match ctx.outputs.push(out) {
                Ok(()) => {}
                Err(unpushed) => {
                    self.held.put(unpushed);
                    // Hold off on more inputs — the held item must clear
                    // before we accept new work.
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if did_work {
            return Ok(StepOutcome::Progress);
        }

        // No input this call. If upstream is drained, flush the boundary
        // state's final partial-record buffer once (guarded by `finalized` —
        // `state.finish()` is not idempotent), emit it, and report `Finished`
        // once nothing remains. `held` is empty here (step 1 returned
        // `Contention` otherwise); a bounced final push is parked in `held`
        // for the held-drain at the top of the next pass.
        if ctx.input.is_drained() {
            if !self.finalized {
                self.finalized = true;
                if let Some(boundary_batch) = self.state.finish()? {
                    *records_this_call += record_count(&boundary_batch.offsets);
                    if !boundary_batch.buffer.is_empty() {
                        let serial = self.next_output_serial;
                        self.next_output_serial += 1;
                        let out = DecompressedBlock {
                            batch_serial: serial,
                            bytes: boundary_batch.buffer,
                        };
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
        let s = FindBamBoundaries::new(1024);
        let p = s.profile();
        assert_eq!(p.name, "FindBamBoundaries");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    // ---------------------------------------------------------------------
    // Domain counter (T-BW2): records, driven through a real pipeline
    // ---------------------------------------------------------------------

    /// Frame `payload` as one record: `[u32 LE block_size][payload]`. Mirrors
    /// `boundaries::state::tests::record` — `BoundaryState` only ever scans
    /// this length-prefix framing, so an opaque payload is sufficient.
    fn record(payload: &[u8]) -> Vec<u8> {
        let mut framed = Vec::with_capacity(4 + payload.len());
        framed.extend_from_slice(&u32::try_from(payload.len()).expect("fits u32").to_le_bytes());
        framed.extend_from_slice(payload);
        framed
    }

    /// `Exclusive` source draining a `Vec<DecompressedBlock>`, one chunk per
    /// `try_run`.
    struct ChunkSource {
        chunks: Vec<DecompressedBlock>,
        held: crate::pipeline::core::held::HeldSlot<
            crate::pipeline::core::Unpushed<DecompressedBlock>,
        >,
    }
    impl Step for ChunkSource {
        type Input = ();
        type Outputs = OrderedBytesSingle<DecompressedBlock>;
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ChunkSource",
                kind: StepKind::Exclusive,
                sticky: true,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                branch_ordering: vec![BranchOrdering::ByItemOrdinal],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(unpushed) = self.held.take()
                && let Err(again) = ctx.outputs.retry(unpushed)
            {
                self.held.put(again);
                return Ok(StepOutcome::Progress);
            }
            let Some(chunk) =
                (if self.chunks.is_empty() { None } else { Some(self.chunks.remove(0)) })
            else {
                return Ok(StepOutcome::Finished);
            };
            if let Err(unpushed) = ctx.outputs.push(chunk) {
                self.held.put(unpushed);
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Serial sink that sleeps briefly per received `DecompressedBlock` before
    /// recording it — mirroring
    /// `fgumi_pipeline_core::tests::CountingSink`'s per-item `thread::sleep`.
    /// `FindBamBoundaries` is the step under test and sits upstream of this
    /// sink, so it can burst through all its `try_run` calls well within the
    /// first sampler tick (see
    /// `fgumi_pipeline_io::source::read_bam::tests::try_run_bumps_blocks_and_bytes_read_counters`
    /// for the same pattern applied to another upstream-of-the-sink counter);
    /// throttling the terminal sink keeps the run open long enough for the
    /// sampler to observe `FindBamBoundaries`' counter at its frozen final
    /// value.
    struct ThrottledChunkSink {
        received: std::sync::Arc<parking_lot::Mutex<Vec<DecompressedBlock>>>,
    }
    impl Step for ThrottledChunkSink {
        type Input = DecompressedBlock;
        type Outputs = ();
        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "ThrottledChunkSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }
        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(block) => {
                    std::thread::sleep(std::time::Duration::from_micros(300));
                    self.received.lock().push(block);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Drives `ChunkSource -> FindBamBoundaries -> ThrottledChunkSink` (via
    /// `new_no_header`, so no BAM header parsing is involved) with telemetry
    /// enabled and asserts the `records` counter lands in the telemetry files
    /// with the exact expected value: `N_RECORDS` framed records, cut into
    /// fixed-size byte chunks that split most records mid-body (exercising the
    /// cross-block carryover path), fed one chunk per `try_run`.
    #[test]
    fn try_run_bumps_records_counter() {
        use fgumi_pipeline_core::builder::{InstrumentationLevel, Pipeline, PipelineConfig};
        use fgumi_pipeline_core::runtime::telemetry::TelemetryConfig;
        use std::time::Duration;

        const N_RECORDS: usize = 500;
        const CHUNK_LEN: usize = 37; // deliberately not record-aligned

        let mut all_bytes = Vec::new();
        for i in 0..N_RECORDS {
            all_bytes.extend_from_slice(&record(format!("payload-{i}").as_bytes()));
        }

        let chunks: Vec<DecompressedBlock> = all_bytes
            .chunks(CHUNK_LEN)
            .enumerate()
            .map(|(i, c)| DecompressedBlock { batch_serial: i as u64, bytes: c.to_vec() })
            .collect();
        assert!(chunks.len() >= 100, "expected many small chunks, got {}", chunks.len());

        let source = ChunkSource { chunks, held: HeldSlot::new() };
        let step = FindBamBoundaries::new_no_header(1024 * 1024);
        let received = std::sync::Arc::new(parking_lot::Mutex::new(Vec::new()));
        let sink = ThrottledChunkSink { received: std::sync::Arc::clone(&received) };

        let dir = std::env::temp_dir()
            .join(format!("fgumi-tbw2-find-bam-boundaries-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let telemetry_stem = dir.join("run");

        let builder = Pipeline::builder();
        builder.chain(source).chain(step).chain(sink).into_sink_marker();
        let pipeline = builder.build().expect("pipeline builds");
        pipeline
            .run(PipelineConfig {
                threads: 1,
                instrumentation: InstrumentationLevel::Summary,
                telemetry: Some(TelemetryConfig {
                    stem: telemetry_stem.clone(),
                    interval: Duration::from_millis(1),
                }),
                ..Default::default()
            })
            .expect("pipeline runs to completion");

        // Ground truth, independent of the sampled telemetry file: every
        // emitted `DecompressedBlock`'s bytes concatenate back to the exact
        // original framed-record stream, so no record was lost, duplicated,
        // or corrupted by the boundary-finding carryover logic.
        let collected = std::mem::take(&mut *received.lock());
        let reassembled: Vec<u8> = collected.iter().flat_map(|b| b.bytes.iter().copied()).collect();
        assert_eq!(reassembled, all_bytes, "reassembled output must match the original input");

        // `FindBamBoundaries` is step index 1 (source=0, step=1, sink=2).
        let names = std::fs::read_to_string(dir.join("run.ticks.counter_names.tsv")).unwrap();
        let name_rows: Vec<&str> = names.lines().skip(1).filter(|l| l.starts_with("1\t")).collect();
        assert_eq!(
            name_rows,
            vec!["1\t0\trecords\trecords"],
            "FindBamBoundaries declares exactly one records counter"
        );

        let counters = std::fs::read_to_string(dir.join("run.ticks.counters.tsv")).unwrap();
        let last_value = counters
            .lines()
            .skip(1)
            .filter(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                f[2] == "1" && f[3] == "0"
            })
            .last()
            .map(|l| l.split('\t').nth(5).unwrap().parse::<u64>().unwrap())
            .expect("records counter recorded at least once");
        assert_eq!(last_value, N_RECORDS as u64, "records counter must equal the true total");

        std::fs::remove_dir_all(&dir).ok();
    }
}
