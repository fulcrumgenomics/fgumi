//! Pool instrumentation: per-stage and per-worker counters for diagnosing
//! performance issues. Modeled on `unified_pipeline::base::PipelineStats`.
//!
//! All counters are atomic u64 so workers can update concurrently without
//! locks. Summarized and logged at pipeline finish.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

use anyhow::{Context, Result};

/// Header written by [`PipelineStats::write_tsv`] and asserted by metrics-
/// output tests. Single source of truth so format drift surfaces as a test
/// failure rather than a silent schema divergence.
pub const METRICS_TSV_HEADER: &str =
    "stage\tkind\twall_time_secs\trecords_in\trecords_out\tbp_in\tbp_out\tseq_busy";

/// Header written by [`PipelineStats::write_queue_tsv`].
pub const QUEUES_TSV_HEADER: &str = "queue_idx\tempty_samples\tmax_depth\tmax_bytes";

/// Shared across all workers and the driver. Atomic counters, no locks.
pub struct PipelineStats {
    num_stages: usize,
    num_workers: usize,

    /// Per-stage names (one per pool stage, index-aligned with counters).
    ///
    /// Populated by the driver via [`PipelineStats::set_stage_names`] before
    /// the pool starts. Used by the pool to tag per-stage progress events
    /// (see `crate::progress::records_out`) by stage name rather than index.
    stage_names: Vec<String>,

    /// Per-stage: successful `Progress` outcomes (stage ran and emitted or suppressed).
    step_count: Vec<AtomicU64>,
    /// Per-stage: wall-time nanoseconds accumulated inside `Stage::process`.
    step_ns: Vec<AtomicU64>,
    /// Per-stage: number of input items consumed from the input queue.
    ///
    /// Incremented once per fresh pop (and once per held-input retry), so the
    /// final value is the total item count flowing INTO the stage.
    step_records_in: Vec<AtomicU64>,
    /// Per-stage: number of output items pushed to the output queue.
    ///
    /// Incremented once per successful `output.push`, whether the push came
    /// from a fresh `process()` emit or from draining the held-output buffer
    /// on a later iteration.
    step_records_out: Vec<AtomicU64>,
    /// Per-stage: total `try_run_stage` attempts (Progress + all backpressure / busy).
    step_attempts: Vec<AtomicU64>,
    /// Per-stage: `BackpressureInput` outcomes.
    step_bp_input: Vec<AtomicU64>,
    /// Per-stage: `BackpressureOutput` outcomes.
    step_bp_output: Vec<AtomicU64>,
    /// Per-stage: `SequentialBusy` outcomes.
    step_seq_busy: Vec<AtomicU64>,

    /// Per-worker: nanoseconds spent in `Backoff::snooze`.
    worker_idle_ns: Vec<AtomicU64>,
    /// Per-worker × per-stage: `Progress` count (which workers did which work).
    worker_step_count: Vec<Vec<AtomicU64>>,
    /// Per-worker: `Progress` events from the eager-exclusive-step path.
    worker_eager_count: Vec<AtomicU64>,

    /// Per-queue: number of iterations the queue was observed empty
    /// (proxy for starvation; sampled by backpressure refresh).
    queue_empty_samples: Vec<AtomicU64>,

    /// Per-queue: the maximum number of items observed in the queue at
    /// any single sample. Tracks contention/backup independent of empty
    /// samples — a queue that is full half the time but never empty
    /// has `queue_empty_samples=0` but a high `queue_max_depth`.
    queue_max_depth: Vec<AtomicU64>,
    /// Per-queue: the maximum number of bytes observed held in the queue
    /// at any single sample.
    queue_max_bytes: Vec<AtomicU64>,

    /// Count of drain-mode activations across all workers.
    drain_mode_activations: AtomicU64,

    /// Per-special-stage wall-time in nanoseconds, keyed by name.
    /// Populated by the driver around `SpecialStage::run_erased`. Special
    /// stages have their own threads, so wall time is a single
    /// total-elapsed measurement (start to finish), not an accumulated
    /// busy-time. Stored with their name so consumers can join with the
    /// pool-stage rows in TSV output.
    special_stage_wall_ns: std::sync::Mutex<Vec<(String, u64)>>,
    /// Per-special-stage records-in observed by the driver. Counted at
    /// the wrapping queue layer (one increment per item entering the
    /// stage's input queue while the special is consuming).
    special_stage_records_in: std::sync::Mutex<Vec<(String, u64)>>,
    /// Per-special-stage records-out observed by the driver. Counted at
    /// the wrapping queue layer (one increment per item the special
    /// stage emits to its output queue).
    special_stage_records_out: std::sync::Mutex<Vec<(String, u64)>>,
}

impl PipelineStats {
    #[must_use]
    pub fn new(num_stages: usize, num_workers: usize, num_queues: usize) -> Self {
        let mk_vec = |n| (0..n).map(|_| AtomicU64::new(0)).collect();
        let worker_step_count: Vec<Vec<AtomicU64>> =
            (0..num_workers).map(|_| mk_vec(num_stages)).collect();
        Self {
            num_stages,
            num_workers,
            step_count: mk_vec(num_stages),
            step_ns: mk_vec(num_stages),
            step_records_in: mk_vec(num_stages),
            step_records_out: mk_vec(num_stages),
            step_attempts: mk_vec(num_stages),
            step_bp_input: mk_vec(num_stages),
            step_bp_output: mk_vec(num_stages),
            step_seq_busy: mk_vec(num_stages),
            worker_idle_ns: mk_vec(num_workers),
            worker_step_count,
            worker_eager_count: mk_vec(num_workers),
            queue_empty_samples: mk_vec(num_queues),
            queue_max_depth: mk_vec(num_queues),
            queue_max_bytes: mk_vec(num_queues),
            drain_mode_activations: AtomicU64::new(0),
            stage_names: (0..num_stages).map(|_| String::new()).collect(),
            special_stage_wall_ns: std::sync::Mutex::new(Vec::new()),
            special_stage_records_in: std::sync::Mutex::new(Vec::new()),
            special_stage_records_out: std::sync::Mutex::new(Vec::new()),
        }
    }

    /// Record the total wall-time spent in a SpecialStage's `run_erased`
    /// call. Called once per Special stage, by the driver, after the
    /// stage's thread has joined.
    pub fn record_special_stage_wall(&self, name: &str, elapsed: Duration) {
        let ns = u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX);
        let mut lock = self.special_stage_wall_ns.lock().unwrap();
        lock.push((name.to_string(), ns));
    }

    /// Record records_in for a Special stage (typically incremented once
    /// per item the wrapping queue consumed during the stage's lifetime).
    pub fn record_special_stage_records_in(&self, name: &str, count: u64) {
        let mut lock = self.special_stage_records_in.lock().unwrap();
        lock.push((name.to_string(), count));
    }

    /// Record records_out for a Special stage.
    pub fn record_special_stage_records_out(&self, name: &str, count: u64) {
        let mut lock = self.special_stage_records_out.lock().unwrap();
        lock.push((name.to_string(), count));
    }

    /// Update the high-water mark for a queue's depth (item count) and
    /// memory bytes. Called by `QueueSummary::summarize`-style code.
    pub fn record_queue_depth_sample(&self, queue_idx: usize, depth: u64, bytes: u64) {
        if queue_idx < self.queue_max_depth.len() {
            // Atomic monotonic max via fetch_update.
            let prev = self.queue_max_depth[queue_idx].load(Ordering::Relaxed);
            if depth > prev {
                let _ = self.queue_max_depth[queue_idx].compare_exchange(
                    prev,
                    depth,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                );
            }
        }
        if queue_idx < self.queue_max_bytes.len() {
            let prev = self.queue_max_bytes[queue_idx].load(Ordering::Relaxed);
            if bytes > prev {
                let _ = self.queue_max_bytes[queue_idx].compare_exchange(
                    prev,
                    bytes,
                    Ordering::Relaxed,
                    Ordering::Relaxed,
                );
            }
        }
    }

    /// Record the per-stage names used when emitting progress events.
    ///
    /// `names.len()` must equal `num_stages`. Called once by the driver
    /// before the pool is started so workers only read the slice after
    /// `Arc<PipelineStats>` has been handed out.
    pub fn set_stage_names(&mut self, names: Vec<String>) {
        debug_assert_eq!(names.len(), self.num_stages);
        self.stage_names = names;
    }

    /// Look up the name registered for `idx`, or `"unknown"` if out of range
    /// or empty.
    #[must_use]
    pub fn stage_name(&self, idx: usize) -> &str {
        match self.stage_names.get(idx).map(String::as_str) {
            Some("") | None => "unknown",
            Some(s) => s,
        }
    }

    pub fn record_step_progress(&self, worker_id: usize, stage_idx: usize, elapsed: Duration) {
        self.step_count[stage_idx].fetch_add(1, Ordering::Relaxed);
        self.step_attempts[stage_idx].fetch_add(1, Ordering::Relaxed);
        // Saturate at u64::MAX rather than truncate; overflow would only occur
        // after ~292 years of accumulated nanoseconds per stage.
        let ns = u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX);
        self.step_ns[stage_idx].fetch_add(ns, Ordering::Relaxed);
        if worker_id < self.num_workers && stage_idx < self.num_stages {
            self.worker_step_count[worker_id][stage_idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Increment the per-stage input counter. Called once per item consumed
    /// from a stage's input queue (fresh pop or held-input retry).
    pub fn record_stage_records_in(&self, stage_idx: usize) {
        if stage_idx < self.step_records_in.len() {
            self.step_records_in[stage_idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Increment the per-stage output counter. Called once per item
    /// successfully pushed to a stage's output queue.
    pub fn record_stage_records_out(&self, stage_idx: usize) {
        if stage_idx < self.step_records_out.len() {
            self.step_records_out[stage_idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn record_eager_progress(&self, worker_id: usize) {
        if worker_id < self.num_workers {
            self.worker_eager_count[worker_id].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn record_bp_input(&self, stage_idx: usize) {
        self.step_attempts[stage_idx].fetch_add(1, Ordering::Relaxed);
        self.step_bp_input[stage_idx].fetch_add(1, Ordering::Relaxed);
    }

    pub fn record_bp_output(&self, stage_idx: usize) {
        self.step_attempts[stage_idx].fetch_add(1, Ordering::Relaxed);
        self.step_bp_output[stage_idx].fetch_add(1, Ordering::Relaxed);
    }

    pub fn record_seq_busy(&self, stage_idx: usize) {
        self.step_attempts[stage_idx].fetch_add(1, Ordering::Relaxed);
        self.step_seq_busy[stage_idx].fetch_add(1, Ordering::Relaxed);
    }

    pub fn record_drained(&self, stage_idx: usize) {
        self.step_attempts[stage_idx].fetch_add(1, Ordering::Relaxed);
    }

    pub fn record_worker_idle(&self, worker_id: usize, elapsed: Duration) {
        if worker_id < self.num_workers {
            let ns = u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX);
            self.worker_idle_ns[worker_id].fetch_add(ns, Ordering::Relaxed);
        }
    }

    pub fn record_queue_empty_sample(&self, queue_idx: usize) {
        if queue_idx < self.queue_empty_samples.len() {
            self.queue_empty_samples[queue_idx].fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn record_drain_activation(&self) {
        self.drain_mode_activations.fetch_add(1, Ordering::Relaxed);
    }

    /// Format a human-readable summary for `tracing::info!`.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn format_summary(&self, stage_names: &[&'static str]) -> String {
        let mut s = String::with_capacity(4096);
        s.push_str("pipeline stats:\n");

        s.push_str("  per-stage:\n");
        s.push_str("    idx  name                 count      attempts    bp_in     bp_out    seq_busy  total_ms   avg_us\n");
        for i in 0..self.num_stages {
            let count = self.step_count[i].load(Ordering::Relaxed);
            let attempts = self.step_attempts[i].load(Ordering::Relaxed);
            let bp_in = self.step_bp_input[i].load(Ordering::Relaxed);
            let bp_out = self.step_bp_output[i].load(Ordering::Relaxed);
            let seq_busy = self.step_seq_busy[i].load(Ordering::Relaxed);
            let ns = self.step_ns[i].load(Ordering::Relaxed);
            let total_ms = ns as f64 / 1_000_000.0;
            let avg_us = if count > 0 { (ns as f64 / count as f64) / 1_000.0 } else { 0.0 };
            let name = stage_names.get(i).copied().unwrap_or("?");
            let _ = writeln!(
                s,
                "    {i:3}  {name:20} {count:10} {attempts:11} {bp_in:9} {bp_out:9} {seq_busy:9} {total_ms:9.1} {avg_us:8.1}"
            );
        }

        s.push_str("  per-worker:\n");
        s.push_str("    id   idle_ms   eager     per-stage progress counts\n");
        for w in 0..self.num_workers {
            let idle = self.worker_idle_ns[w].load(Ordering::Relaxed) as f64 / 1_000_000.0;
            let eager = self.worker_eager_count[w].load(Ordering::Relaxed);
            let counts: Vec<String> = (0..self.num_stages)
                .map(|st| self.worker_step_count[w][st].load(Ordering::Relaxed).to_string())
                .collect();
            let _ = writeln!(s, "    {w:3}  {idle:8.1} {eager:8}   [{}]", counts.join(", "));
        }

        s.push_str("  queue-empty samples (per queue):\n    ");
        let qs: Vec<String> = self
            .queue_empty_samples
            .iter()
            .map(|a| a.load(Ordering::Relaxed).to_string())
            .collect();
        s.push_str(&qs.join(", "));
        s.push('\n');

        let _ = writeln!(
            s,
            "  drain_mode_activations: {}",
            self.drain_mode_activations.load(Ordering::Relaxed)
        );

        s
    }

    /// Write per-stage metrics as a tab-separated file with header:
    /// `stage\twall_time_secs\trecords_in\trecords_out`.
    ///
    /// One row is written per pool stage in the order the stages appear in
    /// the pipeline. `stage_names` is expected to hold one name per pool
    /// stage; if the slice is shorter than `num_stages`, missing entries are
    /// written as `?`.
    ///
    /// Wall time is in seconds with millisecond precision (`{:.3}`). Special
    /// stages (barrier / self-threaded) are not currently instrumented and
    /// are therefore not emitted.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn write_tsv<P: AsRef<Path>>(&self, path: P, stage_names: &[&str]) -> Result<()> {
        let path = path.as_ref();
        let file = File::create(path)
            .with_context(|| format!("failed to create metrics TSV at {}", path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "{METRICS_TSV_HEADER}")
            .with_context(|| format!("failed to write metrics TSV header to {}", path.display()))?;

        // Pool-stage rows: one per stage in the pool stats arrays.
        #[allow(clippy::cast_precision_loss)]
        for i in 0..self.num_stages {
            let name = stage_names.get(i).copied().unwrap_or("?");
            let ns = self.step_ns[i].load(Ordering::Relaxed);
            let wall_time_secs = ns as f64 / 1_000_000_000.0;
            let records_in = self.step_records_in[i].load(Ordering::Relaxed);
            let records_out = self.step_records_out[i].load(Ordering::Relaxed);
            let bp_in = self.step_bp_input[i].load(Ordering::Relaxed);
            let bp_out = self.step_bp_output[i].load(Ordering::Relaxed);
            let seq_busy = self.step_seq_busy[i].load(Ordering::Relaxed);
            writeln!(
                w,
                "{name}\tpool\t{wall_time_secs:.3}\t{records_in}\t{records_out}\t{bp_in}\t{bp_out}\t{seq_busy}"
            )
            .with_context(|| format!("failed to write metrics TSV row to {}", path.display()))?;
        }

        // Special-stage rows: one per recorded special run.
        #[allow(clippy::cast_precision_loss)]
        {
            let walls = self.special_stage_wall_ns.lock().unwrap();
            let ins = self.special_stage_records_in.lock().unwrap();
            let outs = self.special_stage_records_out.lock().unwrap();
            for (i, (name, ns)) in walls.iter().enumerate() {
                let wall_time_secs = (*ns as f64) / 1_000_000_000.0;
                let records_in = ins.get(i).map(|x| x.1).unwrap_or(0);
                let records_out = outs.get(i).map(|x| x.1).unwrap_or(0);
                writeln!(
                    w,
                    "{name}\tspecial\t{wall_time_secs:.3}\t{records_in}\t{records_out}\t0\t0\t0"
                )
                .with_context(|| format!("failed to write metrics TSV row to {}", path.display()))?;
            }
        }

        w.flush().with_context(|| format!("failed to flush metrics TSV to {}", path.display()))?;
        Ok(())
    }

    /// Write per-queue stats as a tab-separated file:
    /// `queue_idx\tempty_samples\tmax_depth\tmax_bytes`.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written.
    pub fn write_queue_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let file = File::create(path)
            .with_context(|| format!("failed to create queues TSV at {}", path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "{QUEUES_TSV_HEADER}")
            .with_context(|| format!("failed to write queues TSV header to {}", path.display()))?;

        for i in 0..self.queue_empty_samples.len() {
            let empty = self.queue_empty_samples[i].load(Ordering::Relaxed);
            let max_depth = self.queue_max_depth[i].load(Ordering::Relaxed);
            let max_bytes = self.queue_max_bytes[i].load(Ordering::Relaxed);
            writeln!(w, "{i}\t{empty}\t{max_depth}\t{max_bytes}").with_context(|| {
                format!("failed to write queues TSV row to {}", path.display())
            })?;
        }

        w.flush().with_context(|| format!("failed to flush queues TSV to {}", path.display()))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufRead;

    #[test]
    fn test_write_tsv_emits_header_and_one_row_per_stage() {
        let stats = PipelineStats::new(3, 1, 4);

        // Stage 0: 10 in, 10 out, 500ms.
        for _ in 0..10 {
            stats.record_stage_records_in(0);
            stats.record_stage_records_out(0);
        }
        stats.record_step_progress(0, 0, Duration::from_millis(500));

        // Stage 1: 10 in, 5 out (suppresses half).
        for _ in 0..10 {
            stats.record_stage_records_in(1);
        }
        for _ in 0..5 {
            stats.record_stage_records_out(1);
        }

        // Stage 2: untouched — 0 in / 0 out.

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let names: Vec<&str> = vec!["StageA", "StageB", "StageC"];
        stats.write_tsv(tmp.path(), &names).unwrap();

        let file = File::open(tmp.path()).unwrap();
        let lines: Vec<String> =
            std::io::BufReader::new(file).lines().map(std::result::Result::unwrap).collect();

        assert_eq!(lines[0], METRICS_TSV_HEADER);
        assert_eq!(lines.len(), 4, "header + 3 stages");

        // Row format: stage, kind, wall_time_secs, records_in, records_out, bp_in, bp_out, seq_busy
        let row_a: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(row_a[0], "StageA");
        assert_eq!(row_a[1], "pool");
        let secs_a: f64 = row_a[2].parse().unwrap();
        assert!(secs_a >= 0.5 - 1e-6, "wall time should be ~0.5s, got {secs_a}");
        assert_eq!(row_a[3], "10");
        assert_eq!(row_a[4], "10");

        let row_b: Vec<&str> = lines[2].split('\t').collect();
        assert_eq!(row_b[0], "StageB");
        assert_eq!(row_b[1], "pool");
        assert_eq!(row_b[3], "10");
        assert_eq!(row_b[4], "5");

        let row_c: Vec<&str> = lines[3].split('\t').collect();
        assert_eq!(row_c[0], "StageC");
        assert_eq!(row_c[1], "pool");
        assert_eq!(row_c[3], "0");
        assert_eq!(row_c[4], "0");
    }

    #[test]
    fn test_write_tsv_errors_on_unwritable_path() {
        let stats = PipelineStats::new(1, 1, 2);
        // Directory does not exist → File::create fails.
        let bad = std::path::Path::new("/nonexistent-dir-a1b2c3d4/metrics.tsv");
        let err = stats.write_tsv(bad, &["X"]).unwrap_err();
        assert!(err.to_string().contains("failed to create metrics TSV"));
    }
}
