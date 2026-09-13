//! Periodic thread/edge/summary telemetry: config, per-worker sample binning,
//! and the six-file TSV writer. Driven by the occupancy sampler (`sampler.rs`).

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::runtime::worker_state::WorkerState;
use crate::step::CounterSpec;
use crate::topology::StepIdx;

/// Where + how often to emit the tick telemetry.
#[derive(Debug, Clone)]
pub struct TelemetryConfig {
    /// File stem; the six TSVs are
    /// `<stem>.ticks.{summary,edges,workers,steps,counters,counter_names}.tsv`
    /// (`summary`/`edges`/`workers`/`counters` carry a row per tick;
    /// `steps`/`counter_names` are written once).
    pub stem: PathBuf,
    /// Emit/sample cadence.
    pub interval: Duration,
}

/// The six per-tick TSV suffixes, in the fixed order the writers create them:
/// `<stem>.ticks.{summary,edges,workers,steps,counters,counter_names}.tsv`.
/// Single source of truth shared by [`TelemetryWriters::open`] (which creates
/// the files) and [`TelemetryConfig::ticks_paths`] (which the chain-build
/// validator uses to reject a telemetry file colliding with the primary output
/// — these paths are *derived* from the stem, so the command-layer
/// output-collision check never sees them).
pub const TICKS_SUFFIXES: [&str; 6] =
    ["summary", "edges", "workers", "steps", "counters", "counter_names"];

/// Derive one per-tick TSV path, `<stem>.ticks.<suffix>.tsv`. This appends to
/// the stem as a filename *suffix*, not as a path component: a bare stem `run`
/// yields `run.ticks.summary.tsv`, and a stem with directories keeps them.
#[must_use]
pub fn ticks_path(stem: &Path, suffix: &str) -> PathBuf {
    let mut p = stem.as_os_str().to_os_string();
    p.push(format!(".ticks.{suffix}.tsv"));
    PathBuf::from(p)
}

impl TelemetryConfig {
    /// The six `<stem>.ticks.*.tsv` files this config's telemetry would write,
    /// derived from [`Self::stem`] via [`ticks_path`] over [`TICKS_SUFFIXES`].
    /// Used by the chain-build validator to reject a telemetry file colliding
    /// with the primary output / rejects / `.bai` sidecar before any writer
    /// truncates a file.
    #[must_use]
    pub fn ticks_paths(&self) -> Vec<PathBuf> {
        TICKS_SUFFIXES.iter().map(|suffix| ticks_path(&self.stem, suffix)).collect()
    }
}

/// Fraction of an emit window a worker spent in each step / non-running state.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowFractions {
    pub per_step: Vec<f32>,
    pub idle: f32,
    pub waiting: f32,
    pub parked: f32,
}

/// Accumulates point-samples for ONE worker over ONE emit window. The sampler
/// calls `record` once per tick per worker, then `fractions` at emit time and
/// `reset` to start the next window.
#[derive(Debug)]
pub struct WorkerBins {
    per_step: Vec<u64>,
    idle: u64,
    waiting: u64,
    parked: u64,
    total: u64,
}

impl WorkerBins {
    #[must_use]
    pub fn new(n_steps: usize) -> Self {
        Self { per_step: vec![0; n_steps], idle: 0, waiting: 0, parked: 0, total: 0 }
    }

    pub fn record(&mut self, state: WorkerState, step: Option<StepIdx>) {
        self.total += 1;
        match (state, step) {
            (WorkerState::Running, Some(idx)) if idx.0 < self.per_step.len() => {
                self.per_step[idx.0] += 1;
            }
            // Running with no/out-of-range step is folded into idle (defensive).
            (WorkerState::Running | WorkerState::Idle, _) => self.idle += 1,
            (WorkerState::Waiting, _) => self.waiting += 1,
            (WorkerState::Parked, _) => self.parked += 1,
        }
    }

    #[must_use]
    pub fn total(&self) -> u64 {
        self.total
    }

    #[allow(clippy::cast_precision_loss)]
    #[must_use]
    pub fn fractions(&self) -> WindowFractions {
        if self.total == 0 {
            return WindowFractions {
                per_step: vec![0.0; self.per_step.len()],
                idle: 0.0,
                waiting: 0.0,
                parked: 0.0,
            };
        }
        let t = self.total as f32;
        WindowFractions {
            per_step: self.per_step.iter().map(|&c| c as f32 / t).collect(),
            idle: self.idle as f32 / t,
            waiting: self.waiting as f32 / t,
            parked: self.parked as f32 / t,
        }
    }

    pub fn reset(&mut self) {
        self.per_step.iter_mut().for_each(|c| *c = 0);
        self.idle = 0;
        self.waiting = 0;
        self.parked = 0;
        self.total = 0;
    }
}

/// One tick's pipeline-wide summary sample (reads processed, memory, queue budget).
pub struct SummarySample {
    pub reads_in: u64,
    pub reads_out: u64,
    pub rss_bytes: Option<u64>,
    pub queue_bytes_used: u64,
    pub queue_bytes_budget: Option<u64>,
}

/// One tick's sample for a single transport edge (queue depth + push/pop deltas).
pub struct EdgeSample {
    pub edge: String,
    pub depth_bytes: Option<u64>,
    pub limit_bytes: Option<u64>,
    pub d_pushed: u64,
    pub d_popped: u64,
    pub d_push_rej: u64,
    pub d_pop_empty: u64,
    /// Bytes pushed this tick. `None` (empty TSV field) for a count/unbounded
    /// edge, which carries no byte accounting; `Some` for a byte-bounded edge.
    pub d_pushed_bytes: Option<u64>,
    /// Bytes popped this tick. `None`/`Some` on the same rule as `d_pushed_bytes`.
    pub d_popped_bytes: Option<u64>,
}

/// One tick's sample for a single worker (current state/step + this-window fractions).
///
/// # v1 worker-state semantics (what the `state`/`f_*` columns actually mean)
///
/// The worker loop (`runtime/driver.rs`) stamps only two of the four
/// [`WorkerState`] variants at runtime:
/// - `Running` — stamped in `dispatch_one_step` immediately before a step's
///   `try_run`, so it covers the whole dispatch **whether or not the step made
///   progress** (a step that returns `NoProgress` is still counted `Running`).
/// - `Parked` — stamped on a no-work iteration when the loop backs off.
/// - `Idle` — only the board's initial pre-run value; at runtime it folds into
///   `Parked`, so it is effectively never sampled once work starts.
/// - `Waiting` — NOT tracked in v1. A worker holding an item it cannot push
///   (producer backpressure) is stamped `Running`, not `Waiting`, because that
///   held item lives inside the step's output handle and is invisible to the
///   loop. Producer backpressure is instead observable in the EDGES TSV via
///   `d_push_rej` (push rejections per tick).
///
/// So `f_waiting` is always 0 and `f_idle` is inert in v1: a reader must not
/// interpret them as meaningful busy/blocked signal. The `waiting`/`idle`
/// columns are retained for schema stability, not because they carry data.
pub struct WorkerSample {
    pub worker: usize,
    /// `"pool"` for a work-stealing pool worker, `"detached"` for a dedicated
    /// detached-driver thread. Decided by the sampler from the slot boundary.
    pub role: &'static str,
    pub state: WorkerState,
    pub step: Option<StepIdx>,
    pub d_serviced: u64,
    pub samples: u64,
    pub fractions: WindowFractions,
}

/// Formats an `Option<u64>` level as its number, or an empty field when `None`
/// (the TSV "NA" encoding for these writers).
fn opt(v: Option<u64>) -> String {
    v.map_or(String::new(), |x| x.to_string())
}

/// Lowercase, stable column value for a [`WorkerState`].
fn state_str(s: WorkerState) -> &'static str {
    match s {
        WorkerState::Running => "running",
        WorkerState::Idle => "idle",
        WorkerState::Waiting => "waiting",
        WorkerState::Parked => "parked",
    }
}

/// Opens and writes rows to the six `<stem>.ticks.{summary,edges,workers,steps,counters,counter_names}.tsv`
/// files. Best-effort: any I/O error is logged once and further rows are silently
/// dropped rather than panicking mid-run (mirrors `sampler.rs::TraceWriter`).
pub struct TelemetryWriters {
    summary: BufWriter<std::fs::File>,
    edges: BufWriter<std::fs::File>,
    workers: BufWriter<std::fs::File>,
    counters: BufWriter<std::fs::File>,
    failed: bool,
}

/// Writes the two static metadata files (`steps`, `counter_names`) that map the
/// numeric ids in the per-tick rows back to names: `steps` = `step_names`
/// indexed `0..N`; `counter_names` = one `step, counter, name, unit` row per
/// declared counter. Any write or flush error is propagated so [`TelemetryWriters::open`]
/// can disable telemetry rather than emit rows whose ids lack a complete mapping.
fn write_static_metadata<S: Write, C: Write>(
    steps: &mut S,
    counter_names: &mut C,
    step_names: &[&'static str],
    step_counters: &[&'static [CounterSpec]],
) -> std::io::Result<()> {
    // Static steps file: one row per step, in index order.
    writeln!(steps, "step\tname")?;
    for (i, name) in step_names.iter().enumerate() {
        writeln!(steps, "{i}\t{name}")?;
    }
    steps.flush()?;
    // Static counter-names file: one row per (step, counter).
    writeln!(counter_names, "step\tcounter\tname\tunit")?;
    for (step, specs) in step_counters.iter().enumerate() {
        for (counter, spec) in specs.iter().enumerate() {
            writeln!(counter_names, "{step}\t{counter}\t{}\t{}", spec.name, spec.unit)?;
        }
    }
    counter_names.flush()
}

// Rows are written with `write!`/`writeln!` and literal tabs rather than through
// the workspace `csv` crate. This is deliberate: the rows are emitted from the
// occupancy sampler on its tick loop (default every 2ms) and the fields are all
// numeric/`&'static str` with no embedded tabs, quotes, or newlines to escape, so
// `csv::Writer`'s per-record quoting-and-buffering machinery would add cost and a
// heap record buffer per tick for no correctness gain. The direct `write!` into
// the already-`BufWriter`-backed file is the cheaper, equally-correct path here.
impl TelemetryWriters {
    /// Creates the six TSVs at
    /// `<stem>.ticks.{summary,edges,workers,steps,counters,counter_names}.tsv`,
    /// writes the two static files once (`steps` = `step_names` indexed 0..N;
    /// `counter_names` = one `step, counter, name, unit` row per declared counter
    /// from `step_counters`), and writes the header row of the four per-tick
    /// files. `step_counters[s]` lists the counters step `s` declared, in slot
    /// order. Returns `None` (after logging a warning) on any create/write error.
    #[must_use]
    pub fn open(
        stem: &Path,
        step_names: &[&'static str],
        step_counters: &[&'static [CounterSpec]],
    ) -> Option<Self> {
        let create =
            |suffix: &str| std::fs::File::create(ticks_path(stem, suffix)).map(BufWriter::new);
        let (
            Ok(mut summary),
            Ok(mut edges),
            Ok(mut workers),
            Ok(mut steps),
            Ok(mut counters),
            Ok(mut counter_names),
        ) = (
            create("summary"),
            create("edges"),
            create("workers"),
            create("steps"),
            create("counters"),
            create("counter_names"),
        )
        else {
            log::warn!(
                "pipeline-telemetry: cannot create one of the .ticks.*.tsv files at stem {}",
                stem.display()
            );
            return None;
        };
        // Static metadata files, written once. A failed static write would leave
        // the sampler emitting rows whose numeric ids have no complete name
        // mapping, so treat any error like a failed per-tick header below and
        // disable telemetry rather than expose incomplete files.
        if write_static_metadata(&mut steps, &mut counter_names, step_names, step_counters).is_err()
        {
            log::warn!(
                "pipeline-telemetry: failed writing a static metadata file at stem {}; disabling telemetry files",
                stem.display()
            );
            return None;
        }
        // Headers.
        let ok_s = writeln!(
            summary,
            "tick\tt_ms\tdt_ms\treads_in\treads_out\trss_bytes\tqueue_bytes_used\tqueue_bytes_budget"
        )
        .is_ok();
        let ok_e = writeln!(
            edges,
            "tick\tt_ms\tedge\tdepth_bytes\tlimit_bytes\td_pushed\td_popped\td_push_rej\td_pop_empty\td_pushed_bytes\td_popped_bytes"
        )
        .is_ok();
        let mut whdr = String::from("tick\tt_ms\tworker\trole\tstate\tstep\td_serviced\tsamples");
        for k in 0..step_names.len() {
            let _ = write!(whdr, "\tf_s{k}");
        }
        whdr.push_str("\tf_idle\tf_waiting\tf_parked");
        let ok_w = writeln!(workers, "{whdr}").is_ok();
        // Per-tick counters file (header-only when no step declared a counter).
        let ok_c = writeln!(counters, "tick\tt_ms\tstep\tcounter\td_value\tvalue").is_ok();
        if !(ok_s && ok_e && ok_w && ok_c) {
            log::warn!("pipeline-telemetry: failed writing a header; disabling telemetry files");
            return None;
        }
        Some(Self { summary, edges, workers, counters, failed: false })
    }

    /// Writes one row to the `summary` TSV.
    pub fn write_summary_row(&mut self, tick: u64, t_ms: f64, dt_ms: f64, s: &SummarySample) {
        if self.failed {
            return;
        }
        let line = format!(
            "{tick}\t{t_ms:.3}\t{dt_ms:.3}\t{}\t{}\t{}\t{}\t{}",
            s.reads_in,
            s.reads_out,
            opt(s.rss_bytes),
            s.queue_bytes_used,
            opt(s.queue_bytes_budget),
        );
        self.write_summary_line(&line);
    }

    /// Writes one row to the `edges` TSV.
    pub fn write_edge_row(&mut self, tick: u64, t_ms: f64, e: &EdgeSample) {
        if self.failed {
            return;
        }
        let line = format!(
            "{tick}\t{t_ms:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.edge,
            opt(e.depth_bytes),
            opt(e.limit_bytes),
            e.d_pushed,
            e.d_popped,
            e.d_push_rej,
            e.d_pop_empty,
            opt(e.d_pushed_bytes),
            opt(e.d_popped_bytes),
        );
        self.write_edges_line(&line);
    }

    /// Writes one row to the `counters` TSV: the per-tick delta (`d_value`) and
    /// cumulative total (`value`) of one step's one domain counter.
    pub fn write_counter_row(
        &mut self,
        tick: u64,
        t_ms: f64,
        step: usize,
        counter: usize,
        d_value: u64,
        value: u64,
    ) {
        if self.failed {
            return;
        }
        let line = format!("{tick}\t{t_ms:.3}\t{step}\t{counter}\t{d_value}\t{value}");
        self.write_counters_line(&line);
    }

    /// Writes one row to the `workers` TSV.
    pub fn write_worker_row(&mut self, tick: u64, t_ms: f64, w: &WorkerSample) {
        if self.failed {
            return;
        }
        let step = w.step.map_or(String::new(), |s| s.0.to_string());
        let mut line = format!(
            "{tick}\t{t_ms:.3}\t{}\t{}\t{}\t{}\t{}\t{}",
            w.worker,
            w.role,
            state_str(w.state),
            step,
            w.d_serviced,
            w.samples,
        );
        for f in &w.fractions.per_step {
            let _ = write!(line, "\t{f:.4}");
        }
        let _ = write!(
            line,
            "\t{:.4}\t{:.4}\t{:.4}",
            w.fractions.idle, w.fractions.waiting, w.fractions.parked
        );
        self.write_workers_line(&line);
    }

    fn write_summary_line(&mut self, line: &str) {
        if writeln!(self.summary, "{line}").is_err() {
            log::warn!("pipeline-telemetry: write failed; dropping further rows");
            self.failed = true;
        }
    }

    fn write_edges_line(&mut self, line: &str) {
        if writeln!(self.edges, "{line}").is_err() {
            log::warn!("pipeline-telemetry: write failed; dropping further rows");
            self.failed = true;
        }
    }

    fn write_workers_line(&mut self, line: &str) {
        if writeln!(self.workers, "{line}").is_err() {
            log::warn!("pipeline-telemetry: write failed; dropping further rows");
            self.failed = true;
        }
    }

    fn write_counters_line(&mut self, line: &str) {
        if writeln!(self.counters, "{line}").is_err() {
            log::warn!("pipeline-telemetry: write failed; dropping further rows");
            self.failed = true;
        }
    }

    /// Flushes all four per-tick writers (the static `steps` / `counter_names`
    /// files are flushed once at `open` and never written again).
    pub fn flush(&mut self) {
        if self.failed {
            return;
        }
        for w in [&mut self.summary, &mut self.edges, &mut self.workers, &mut self.counters] {
            if w.flush().is_err() {
                log::warn!("pipeline-telemetry: flush failed; buffered rows may be lost");
                self.failed = true;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_fractions_sum_to_one_and_split_states() {
        let mut bins = WorkerBins::new(2); // steps 0,1
        bins.record(WorkerState::Running, Some(StepIdx(1)));
        bins.record(WorkerState::Running, Some(StepIdx(1)));
        bins.record(WorkerState::Parked, None);
        bins.record(WorkerState::Idle, None);
        assert_eq!(bins.total(), 4);
        let f = bins.fractions();
        assert_eq!(f.per_step.len(), 2);
        assert!((f.per_step[0] - 0.0).abs() < 1e-6);
        assert!((f.per_step[1] - 0.5).abs() < 1e-6);
        assert!((f.parked - 0.25).abs() < 1e-6);
        assert!((f.idle - 0.25).abs() < 1e-6);
        assert!((f.waiting - 0.0).abs() < 1e-6);
        let sum: f32 = f.per_step.iter().sum::<f32>() + f.idle + f.waiting + f.parked;
        assert!((sum - 1.0).abs() < 1e-5, "fractions sum to 1, got {sum}");
    }

    #[test]
    fn empty_window_is_all_zero() {
        let bins = WorkerBins::new(3);
        assert_eq!(bins.total(), 0);
        let f = bins.fractions();
        assert!(f.per_step.iter().all(|&x| x == 0.0));
        assert_eq!((f.idle, f.waiting, f.parked), (0.0, 0.0, 0.0));
    }

    #[test]
    fn running_with_no_step_counts_as_idle_bucket() {
        // Running but step==None (shouldn't happen, but be defensive) must not
        // index out of bounds; fold it into idle.
        let mut bins = WorkerBins::new(1);
        bins.record(WorkerState::Running, None);
        let f = bins.fractions();
        assert!((f.idle - 1.0).abs() < 1e-6);
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn writers_emit_four_files_with_headers_and_rows() {
        use std::io::Read;
        // Step 1 ("sort") declares one counter; steps 0/2 declare none. A `const`
        // so the inner `&[CounterSpec::new(..)]` slices are `'static` (a const-fn
        // call is not auto-promoted to `'static` in an rvalue context).
        const STEP_COUNTERS: &[&[CounterSpec]] =
            &[&[], &[CounterSpec::new("records", "records")], &[]];
        let dir = std::env::temp_dir().join(format!("fgumi-tele-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let stem = dir.join("run");
        let mut w = TelemetryWriters::open(&stem, &["read", "sort", "write"], STEP_COUNTERS)
            .expect("writers open");
        w.write_summary_row(
            0,
            0.5,
            0.5,
            &SummarySample {
                reads_in: 0,
                reads_out: 0,
                rss_bytes: Some(1024),
                queue_bytes_used: 0,
                queue_bytes_budget: Some(4096),
            },
        );
        w.write_edge_row(
            0,
            0.5,
            &EdgeSample {
                edge: "read__sort#0b0".to_string(),
                depth_bytes: Some(10),
                limit_bytes: Some(100),
                d_pushed: 5,
                d_popped: 3,
                d_push_rej: 0,
                d_pop_empty: 1,
                d_pushed_bytes: Some(500),
                d_popped_bytes: Some(300),
            },
        );
        w.write_edge_row(
            0,
            0.5,
            &EdgeSample {
                edge: "count__edge#1b0".to_string(),
                depth_bytes: None,
                limit_bytes: None,
                d_pushed: 0,
                d_popped: 0,
                d_push_rej: 0,
                d_pop_empty: 0,
                d_pushed_bytes: None,
                d_popped_bytes: None,
            },
        );
        // One counter row for step 1's counter 0 (delta 7, cumulative 7).
        w.write_counter_row(0, 0.5, 1, 0, 7, 7);
        w.write_worker_row(
            0,
            0.5,
            &WorkerSample {
                worker: 0,
                role: "pool",
                state: WorkerState::Running,
                step: Some(StepIdx(1)),
                d_serviced: 3,
                samples: 4,
                fractions: WindowFractions {
                    per_step: vec![0.25, 0.75, 0.0],
                    idle: 0.0,
                    waiting: 0.0,
                    parked: 0.0,
                },
            },
        );
        w.flush();

        let read = |suffix: &str| {
            let mut s = String::new();
            std::fs::File::open(dir.join(format!("run.ticks.{suffix}.tsv")))
                .unwrap()
                .read_to_string(&mut s)
                .unwrap();
            s
        };
        let steps = read("steps");
        assert!(steps.starts_with("step\tname\n"), "steps header");
        assert!(steps.contains("0\tread\n") && steps.contains("2\twrite\n"));
        let summary = read("summary");
        assert!(summary.lines().next().unwrap()
            .starts_with("tick\tt_ms\tdt_ms\treads_in\treads_out\trss_bytes\tqueue_bytes_used\tqueue_bytes_budget"));
        let edges = read("edges");
        let edges_hdr = edges.lines().next().unwrap();
        assert!(
            edges_hdr.ends_with("d_pop_empty\td_pushed_bytes\td_popped_bytes"),
            "byte columns are appended at the end: {edges_hdr}"
        );
        // Byte edge has numeric depth AND numeric byte deltas at the row's tail.
        assert!(edges.contains("read__sort#0b0\t10\t100\t5\t3\t0\t1\t500\t300"));
        // Count edge (None) emits empty fields for depth/limit AND the byte deltas.
        assert!(
            edges.contains("count__edge#1b0\t\t\t0\t0\t0\t0\t\t"),
            "None depth/limit/byte deltas are empty fields"
        );
        // Counter files: static names + per-tick values.
        let counter_names = read("counter_names");
        assert!(counter_names.starts_with("step\tcounter\tname\tunit\n"), "counter_names header");
        assert!(
            counter_names.contains("1\t0\trecords\trecords\n"),
            "step 1's counter 0 is named: {counter_names}"
        );
        let counters = read("counters");
        assert!(
            counters.lines().next().unwrap() == "tick\tt_ms\tstep\tcounter\td_value\tvalue",
            "counters header"
        );
        assert!(counters.contains("0\t0.500\t1\t0\t7\t7"), "counter row: {counters}");
        let workers = read("workers");
        let hdr = workers.lines().next().unwrap();
        assert!(
            hdr.contains("worker\trole\tstate\tstep\td_serviced\tsamples"),
            "role column sits right after worker"
        );
        assert!(
            hdr.contains("f_s0\tf_s1\tf_s2\tf_idle\tf_waiting\tf_parked"),
            "per-step + state cols"
        );
        let row = workers.lines().nth(1).unwrap();
        assert!(row.contains("\tpool\trunning\t1\t3\t4\t0.2500\t0.7500\t0.0000"), "row: {row}");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn write_static_metadata_propagates_write_and_flush_errors() {
        const STEP_COUNTERS: &[&[CounterSpec]] = &[&[], &[CounterSpec::new("records", "records")]];

        // A writer whose every write/flush fails, to exercise the error path a
        // real file only hits on a full disk / I/O error.
        struct FailingWriter;
        impl Write for FailingWriter {
            fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
                Err(std::io::Error::other("boom"))
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Err(std::io::Error::other("boom"))
            }
        }

        // A writer that accepts every write but fails only on `flush`, so the
        // `flush()?` calls (which `FailingWriter` never reaches, since it fails
        // at the first `writeln!`) are exercised on their own.
        struct FlushFailingWriter;
        impl Write for FlushFailingWriter {
            fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
                Ok(buf.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Err(std::io::Error::other("flush boom"))
            }
        }

        // Happy path: both files are in-memory buffers and every write succeeds.
        let mut steps = Vec::<u8>::new();
        let mut counter_names = Vec::<u8>::new();
        assert!(
            write_static_metadata(&mut steps, &mut counter_names, &["read", "sort"], STEP_COUNTERS)
                .is_ok()
        );
        assert!(steps.starts_with(b"step\tname\n"));
        assert!(counter_names.starts_with(b"step\tcounter\tname\tunit\n"));

        // A failing `steps` writer propagates the error (so `open` returns None).
        let mut counter_names = Vec::<u8>::new();
        assert!(
            write_static_metadata(
                &mut FailingWriter,
                &mut counter_names,
                &["read", "sort"],
                STEP_COUNTERS,
            )
            .is_err()
        );

        // A failing `counter_names` writer (steps OK) also propagates.
        let mut steps = Vec::<u8>::new();
        assert!(
            write_static_metadata(
                &mut steps,
                &mut FailingWriter,
                &["read", "sort"],
                STEP_COUNTERS,
            )
            .is_err()
        );

        // A flush-only failure on `steps` propagates from the `steps.flush()?`.
        let mut counter_names = Vec::<u8>::new();
        assert!(
            write_static_metadata(
                &mut FlushFailingWriter,
                &mut counter_names,
                &["read", "sort"],
                STEP_COUNTERS,
            )
            .is_err()
        );

        // A flush-only failure on `counter_names` (steps OK) propagates from the
        // final `counter_names.flush()`.
        let mut steps = Vec::<u8>::new();
        assert!(
            write_static_metadata(
                &mut steps,
                &mut FlushFailingWriter,
                &["read", "sort"],
                STEP_COUNTERS,
            )
            .is_err()
        );
    }
}
