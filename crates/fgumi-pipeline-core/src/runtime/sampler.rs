//! Background occupancy sampler for `--pipeline-trace`.
//!
//! When instrumentation is on, `Pipeline::run` spawns one
//! [`run_occupancy_sampler`] thread (same lifecycle slot as the deadlock
//! monitor / queue rebalancer) that periodically reads each byte-bounded edge's
//! depth (`current_bytes / limit_bytes`) and feeds it to the edge's
//! [`EdgeMetrics`](super::metrics::EdgeMetrics) occupancy histogram. It is
//! read-only over the live queues (one `Relaxed` load per edge per tick), so it
//! never perturbs the worker hot path; only edges with a `depth_source`
//! (byte-bounded) are sampled — count/unbounded edges still get their push/pop
//! counters, just no occupancy histogram.

use std::fmt::Write as _;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use super::contexts::RegisteredEdge;

/// Default sampling interval. Two milliseconds is cheap (one atomic load per
/// edge) yet fine-grained enough to resolve a chain's phase structure over a
/// multi-second run.
pub const DEFAULT_SAMPLE_INTERVAL: Duration = Duration::from_millis(2);

/// Poll each byte-bounded edge's occupancy into its histogram until `stop` is
/// set. Edges without a `depth_source` (count/unbounded) are skipped. When
/// `trace_path` is `Some` (the `Timeline` level), also append one TSV row per
/// tick — `t_ms` plus, per edge, its depth fraction and cumulative
/// pushed/popped item counts — so the run's phase structure can be plotted.
pub fn run_occupancy_sampler(
    stop: &AtomicBool,
    edges: &[RegisteredEdge],
    interval: Duration,
    trace_path: Option<PathBuf>,
) {
    let mut trace = trace_path.and_then(|p| TraceWriter::open(&p, edges));
    let start = Instant::now();
    while !stop.load(Ordering::Relaxed) {
        sample_once(edges);
        if let Some(t) = trace.as_mut() {
            t.write_row(edges, start.elapsed());
        }
        std::thread::sleep(interval);
    }
    // One final sample so a short run that finishes before the first sleep
    // still records at least one occupancy point per edge.
    sample_once(edges);
    if let Some(mut t) = trace {
        t.write_row(edges, start.elapsed());
        t.flush();
    }
}

/// Per-tick TSV writer for the `Timeline` level. Best-effort: a write error is
/// logged once and further rows are dropped (instrumentation never aborts a run).
struct TraceWriter {
    writer: BufWriter<std::fs::File>,
    failed: bool,
}

impl TraceWriter {
    /// Open `path` and write the header (`t_ms` + three columns per edge).
    /// Returns `None` (with a warning) if the file can't be created.
    fn open(path: &std::path::Path, edges: &[RegisteredEdge]) -> Option<Self> {
        match std::fs::File::create(path) {
            Ok(file) => {
                let mut writer = BufWriter::new(file);
                let mut header = String::from("t_ms");
                for e in edges {
                    let edge =
                        format!("{}__{}", e.producer_name, e.consumer_name.unwrap_or("sink"));
                    let _ = write!(header, "\t{edge}.depth\t{edge}.pushed\t{edge}.popped");
                }
                if writeln!(writer, "{header}").is_err() {
                    log::warn!(
                        "pipeline-trace: failed to write timeline header to {}",
                        path.display()
                    );
                    return None;
                }
                Some(Self { writer, failed: false })
            }
            Err(e) => {
                log::warn!("pipeline-trace: cannot create timeline file {}: {e}", path.display());
                None
            }
        }
    }

    #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
    fn write_row(&mut self, edges: &[RegisteredEdge], elapsed: Duration) {
        if self.failed {
            return;
        }
        let mut row = format!("{}", elapsed.as_millis());
        for e in edges {
            let depth = e.depth_source.as_ref().map_or(0.0, |s| {
                let limit = s.limit_bytes();
                if limit == 0 { 0.0 } else { s.current_bytes() as f32 / limit as f32 }
            });
            let ms = e.metrics.snapshot();
            let _ = write!(row, "\t{depth:.3}\t{}\t{}", ms.pushed_items, ms.popped_items);
        }
        if writeln!(self.writer, "{row}").is_err() {
            log::warn!("pipeline-trace: timeline write failed; dropping further rows");
            self.failed = true;
        }
    }

    fn flush(&mut self) {
        let _ = self.writer.flush();
    }
}

/// One sampling sweep over all edges. Extracted so tests can drive a single
/// deterministic tick without the sleep loop.
#[allow(clippy::cast_precision_loss)]
pub fn sample_once(edges: &[RegisteredEdge]) {
    for edge in edges {
        if let Some(src) = &edge.depth_source {
            let limit = src.limit_bytes();
            if limit > 0 {
                let frac = src.current_bytes() as f32 / limit as f32;
                edge.metrics.record_depth(frac);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::sync::atomic::AtomicBool;

    use crate::item::HeapSize;
    use crate::queues::{BoundedQueueHandle, ByteBoundedQueue, ItemQueue};
    use crate::runtime::metrics::EdgeMetrics;
    use crate::topology::{BranchIdx, StepIdx};

    #[derive(Debug)]
    struct Heavy(Vec<u8>);
    impl HeapSize for Heavy {
        fn heap_size(&self) -> usize {
            self.0.len()
        }
    }

    fn edge_over(
        metrics: Arc<EdgeMetrics>,
        depth_source: Option<Arc<dyn BoundedQueueHandle>>,
    ) -> RegisteredEdge {
        RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "producer",
            consumer_step: Some(StepIdx(1)),
            consumer_name: Some("consumer"),
            branch: BranchIdx(0),
            metrics,
            depth_source,
        }
    }

    #[test]
    fn sample_once_records_occupancy_from_depth_source() {
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 500])).unwrap(); // 50% of the 1000-byte budget
        let edge = edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>));
        for _ in 0..10 {
            sample_once(std::slice::from_ref(&edge));
        }
        let s = m.snapshot();
        assert_eq!(s.depth_samples, 10);
        assert!((s.mean_occupancy - 0.5).abs() < 0.05, "mean ≈ 0.5, got {}", s.mean_occupancy);
    }

    #[test]
    fn count_edge_without_depth_source_is_skipped() {
        // An edge with no depth_source (count/unbounded) records no occupancy.
        let m = EdgeMetrics::new();
        let edge = edge_over(Arc::clone(&m), None);
        for _ in 0..5 {
            sample_once(std::slice::from_ref(&edge));
        }
        assert_eq!(m.snapshot().depth_samples, 0, "no depth source → no occupancy samples");
    }

    #[test]
    fn run_occupancy_sampler_stops_and_records() {
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 800])).unwrap();
        let edges =
            vec![edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>))];
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let handle = std::thread::spawn(move || {
            run_occupancy_sampler(&stop_c, &edges, Duration::from_millis(1), None);
        });
        std::thread::sleep(Duration::from_millis(30));
        stop.store(true, Ordering::Relaxed);
        handle.join().unwrap();
        let s = m.snapshot();
        assert!(s.depth_samples > 0, "sampler recorded at least one tick");
        assert!((s.mean_occupancy - 0.8).abs() < 0.1, "mean ≈ 0.8, got {}", s.mean_occupancy);
    }

    #[test]
    fn timeline_tsv_has_header_and_rows() {
        let m = EdgeMetrics::new();
        let q = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));
        q.try_push(Heavy(vec![0; 400])).unwrap();
        let edges =
            vec![edge_over(Arc::clone(&m), Some(Arc::clone(&q) as Arc<dyn BoundedQueueHandle>))];
        let dir = std::env::temp_dir();
        let path = dir.join(format!("fgumi-trace-test-{}.tsv", std::process::id()));
        let stop = Arc::new(AtomicBool::new(false));
        let stop_c = Arc::clone(&stop);
        let path_c = path.clone();
        let handle = std::thread::spawn(move || {
            run_occupancy_sampler(&stop_c, &edges, Duration::from_millis(2), Some(path_c));
        });
        std::thread::sleep(Duration::from_millis(30));
        stop.store(true, Ordering::Relaxed);
        handle.join().unwrap();

        let content = std::fs::read_to_string(&path).expect("trace file written");
        let _ = std::fs::remove_file(&path);
        let mut lines = content.lines();
        let header = lines.next().expect("header row");
        assert!(header.starts_with("t_ms"), "header begins with t_ms");
        assert!(header.contains("producer__consumer.depth"), "per-edge depth column");
        let rows: Vec<&str> = lines.collect();
        assert!(!rows.is_empty(), "at least one data row");
        // First field of a data row is a monotonic t_ms integer.
        let first_t: u128 = rows[0].split('\t').next().unwrap().parse().expect("t_ms is an int");
        let last_t: u128 = rows.last().unwrap().split('\t').next().unwrap().parse().unwrap();
        assert!(last_t >= first_t, "t_ms is monotonic");
    }
}
