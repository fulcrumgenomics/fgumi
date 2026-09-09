//! Periodic thread/edge/summary telemetry: config, per-worker sample binning,
//! and the four-file TSV writer. Driven by the occupancy sampler (`sampler.rs`).

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::runtime::worker_state::WorkerState;
use crate::topology::StepIdx;

/// Where + how often to emit the tick telemetry.
#[derive(Debug, Clone)]
pub struct TelemetryConfig {
    /// File stem; the four TSVs are `<stem>.ticks.{summary,edges,workers,steps}.tsv`.
    pub stem: PathBuf,
    /// Emit/sample cadence.
    pub interval: Duration,
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
}

/// One tick's sample for a single worker (current state/step + this-window fractions).
pub struct WorkerSample {
    pub worker: usize,
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

/// Opens and writes rows to the four `<stem>.ticks.{summary,edges,workers,steps}.tsv`
/// files. Best-effort: any I/O error is logged once and further rows are silently
/// dropped rather than panicking mid-run (mirrors `sampler.rs::TraceWriter`).
pub struct TelemetryWriters {
    summary: BufWriter<std::fs::File>,
    edges: BufWriter<std::fs::File>,
    workers: BufWriter<std::fs::File>,
    failed: bool,
}

impl TelemetryWriters {
    /// Creates the four TSVs at `<stem>.ticks.{summary,edges,workers,steps}.tsv`,
    /// writes the `steps` file once (`step_names` indexed 0..N), and writes the
    /// header row of the other three. Returns `None` (after logging a warning) on
    /// any create/write error.
    #[must_use]
    pub fn open(stem: &Path, step_names: &[&'static str], n_workers: usize) -> Option<Self> {
        let path = |suffix: &str| {
            let mut p = stem.as_os_str().to_os_string();
            p.push(format!(".ticks.{suffix}.tsv"));
            PathBuf::from(p)
        };
        let create = |suffix: &str| std::fs::File::create(path(suffix)).map(BufWriter::new);
        let (Ok(mut summary), Ok(mut edges), Ok(mut workers), Ok(mut steps)) =
            (create("summary"), create("edges"), create("workers"), create("steps"))
        else {
            log::warn!(
                "pipeline-telemetry: cannot create one of the .ticks.*.tsv files at stem {}",
                stem.display()
            );
            return None;
        };
        // Static steps file, written once.
        let _ = writeln!(steps, "step\tname");
        for (i, name) in step_names.iter().enumerate() {
            let _ = writeln!(steps, "{i}\t{name}");
        }
        let _ = steps.flush();
        // Headers.
        let ok_s = writeln!(
            summary,
            "tick\tt_ms\tdt_ms\treads_in\treads_out\trss_bytes\tqueue_bytes_used\tqueue_bytes_budget"
        )
        .is_ok();
        let ok_e = writeln!(
            edges,
            "tick\tt_ms\tedge\tdepth_bytes\tlimit_bytes\td_pushed\td_popped\td_push_rej\td_pop_empty"
        )
        .is_ok();
        let mut whdr = String::from("tick\tt_ms\tworker\tstate\tstep\td_serviced\tsamples");
        for k in 0..step_names.len() {
            let _ = write!(whdr, "\tf_s{k}");
        }
        whdr.push_str("\tf_idle\tf_waiting\tf_parked");
        let ok_w = writeln!(workers, "{whdr}").is_ok();
        let _ = n_workers; // header shape is per-step; worker count only bounds rows.
        if !(ok_s && ok_e && ok_w) {
            log::warn!("pipeline-telemetry: failed writing a header; disabling telemetry files");
            return None;
        }
        Some(Self { summary, edges, workers, failed: false })
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
            "{tick}\t{t_ms:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.edge,
            opt(e.depth_bytes),
            opt(e.limit_bytes),
            e.d_pushed,
            e.d_popped,
            e.d_push_rej,
            e.d_pop_empty,
        );
        self.write_edges_line(&line);
    }

    /// Writes one row to the `workers` TSV.
    pub fn write_worker_row(&mut self, tick: u64, t_ms: f64, w: &WorkerSample) {
        if self.failed {
            return;
        }
        let step = w.step.map_or(String::new(), |s| s.0.to_string());
        let mut line = format!(
            "{tick}\t{t_ms:.3}\t{}\t{}\t{}\t{}\t{}",
            w.worker,
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

    /// Flushes all three per-tick writers (the `steps` file is flushed once at
    /// `open` and never written again).
    pub fn flush(&mut self) {
        if self.failed {
            return;
        }
        for w in [&mut self.summary, &mut self.edges, &mut self.workers] {
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
    fn writers_emit_four_files_with_headers_and_rows() {
        use std::io::Read;
        let dir = std::env::temp_dir().join(format!("fgumi-tele-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let stem = dir.join("run");
        let mut w =
            TelemetryWriters::open(&stem, &["read", "sort", "write"], 2).expect("writers open");
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
            },
        );
        w.write_worker_row(
            0,
            0.5,
            &WorkerSample {
                worker: 0,
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
        // Byte edge has a numeric depth; count edge (None) emits empty fields → "\t\t".
        assert!(edges.contains("read__sort#0b0\t10\t100\t5\t3\t0\t1"));
        assert!(
            edges.contains("count__edge#1b0\t\t\t0\t0\t0\t0"),
            "None depth/limit are empty fields"
        );
        let workers = read("workers");
        let hdr = workers.lines().next().unwrap();
        assert!(
            hdr.contains("f_s0\tf_s1\tf_s2\tf_idle\tf_waiting\tf_parked"),
            "per-step + state cols"
        );
        assert!(
            workers.lines().nth(1).unwrap().contains("running\t1\t3\t4\t0.2500\t0.7500\t0.0000")
        );

        std::fs::remove_dir_all(&dir).ok();
    }
}
