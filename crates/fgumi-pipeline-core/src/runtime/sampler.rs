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
    let mut sampled_in_loop = false;
    while !stop.load(Ordering::Relaxed) {
        sample_once(edges);
        if let Some(t) = trace.as_mut() {
            t.write_row(edges, start.elapsed());
        }
        sampled_in_loop = true;
        std::thread::sleep(interval);
    }
    // Guard the final sample: only take it when the loop never sampled (a run
    // so short `stop` was already set before the first iteration). Sampling
    // unconditionally here would add an extra occupancy point + timeline row
    // taken AFTER the pipeline already drained, biasing the mean toward the
    // empty final state.
    if !sampled_in_loop {
        sample_once(edges);
        if let Some(t) = trace.as_mut() {
            t.write_row(edges, start.elapsed());
        }
    }
    if let Some(mut t) = trace {
        t.flush();
    }
}

/// Bytes buffered in an ordered edge's `ReorderStage` overflow stash (0 for a
/// direct/count/unbounded edge). Added to the transport queue's `current_bytes`
/// when sampling depth so an ordered edge reflects total buffered bytes rather
/// than reading empty while items pile in the reorder buffer awaiting an earlier
/// ordinal.
fn reorder_stash_bytes(edge: &RegisteredEdge) -> u64 {
    edge.reorder_depth.as_ref().map_or(0, |r| r.current_buffer_bytes())
}

/// Column-name prefix for one edge's timeline columns. Includes the producer
/// step index and output branch so fan-out edges (one producer, several
/// branches) and repeated step names produce distinct, collision-free headers —
/// a bare `producer__consumer` prefix duplicates columns whenever two edges
/// share both names. `(producer_step, branch)` uniquely identifies an edge.
fn edge_column_prefix(e: &RegisteredEdge) -> String {
    format!(
        "{}__{}#{}b{}",
        e.producer_name,
        e.consumer_name.unwrap_or("sink"),
        e.producer_step.0,
        e.branch.0,
    )
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
                    let edge = edge_column_prefix(e);
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
            // Count/unbounded edges (no `depth_source`) are unsampled — emit `NA`
            // rather than `0.000`, which would misread as "empty" instead of
            // "not measured". Byte-bounded edges report total buffered depth
            // (transport + reorder stash) as a fraction of the limit.
            let depth = e
                .depth_source
                .as_ref()
                .and_then(|s| {
                    let limit = s.limit_bytes();
                    (limit > 0).then(|| {
                        let occupied = s.current_bytes().saturating_add(reorder_stash_bytes(e));
                        format!("{:.3}", occupied as f32 / limit as f32)
                    })
                })
                .unwrap_or_else(|| "NA".to_string());
            let ms = e.metrics.snapshot();
            let _ = write!(row, "\t{depth}\t{}\t{}", ms.pushed_items, ms.popped_items);
        }
        if writeln!(self.writer, "{row}").is_err() {
            log::warn!("pipeline-trace: timeline write failed; dropping further rows");
            self.failed = true;
        }
    }

    fn flush(&mut self) {
        if self.failed {
            return;
        }
        // A dropped flush error can silently lose buffered rows after every
        // write appeared to succeed — warn and mark the writer failed, matching
        // `write_row`'s best-effort error handling.
        if self.writer.flush().is_err() {
            log::warn!("pipeline-trace: timeline flush failed; buffered rows may be lost");
            self.failed = true;
        }
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
                let occupied = src.current_bytes().saturating_add(reorder_stash_bytes(edge));
                edge.metrics.record_depth(occupied, limit);
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
            reorder_depth: None,
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
    fn edge_columns_are_unique_when_step_names_collide() {
        // Regression: two edges that share producer AND consumer names (fan-out,
        // or duplicate step names) must still produce distinct TSV columns — a
        // bare `producer__consumer` prefix would emit duplicate column headers.
        let m0 = EdgeMetrics::new();
        let m1 = EdgeMetrics::new();
        let e0 = RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "dup",
            consumer_step: Some(StepIdx(1)),
            consumer_name: Some("sink"),
            branch: BranchIdx(0),
            metrics: m0,
            depth_source: None,
            reorder_depth: None,
        };
        // Same names, different (producer_step, branch): a fan-out sibling.
        let e1 = RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "dup",
            consumer_step: Some(StepIdx(2)),
            consumer_name: Some("sink"),
            branch: BranchIdx(1),
            metrics: m1,
            depth_source: None,
            reorder_depth: None,
        };
        let p0 = edge_column_prefix(&e0);
        let p1 = edge_column_prefix(&e1);
        assert_ne!(p0, p1, "colliding names must yield distinct column prefixes");
        assert_eq!(p0, "dup__sink#0b0");
        assert_eq!(p1, "dup__sink#0b1");
    }

    #[test]
    fn ordered_edge_occupancy_includes_reorder_stash() {
        use crate::queues::CountBoundedQueue;
        use crate::reorder::{ReorderCapHandle, ReorderStage, Sequenced};
        // The transport (occupancy depth source) is empty, but the reorder stash
        // holds 400 buffered bytes waiting for an earlier ordinal. Sampled
        // occupancy must reflect the stash (400/1000 = 0.4), not read empty —
        // otherwise a producer-skewed ordered edge looks idle while backed up.
        let m = EdgeMetrics::new();
        let transport = Arc::new(ByteBoundedQueue::<Heavy>::new(1000));

        let reorder_transport: Arc<dyn ItemQueue<Sequenced<Heavy>>> =
            Arc::new(CountBoundedQueue::<Sequenced<Heavy>>::new(8));
        let stage = Arc::new(ReorderStage::new(reorder_transport));
        stage.try_push(1, Heavy(vec![0; 400])).unwrap(); // out-of-order → stashed
        assert!(stage.try_pop_in_order().is_none(), "ordinal 0 absent → nothing pops");
        assert!(stage.current_buffer_bytes() >= 400, "stash holds the buffered bytes");

        let edge = RegisteredEdge {
            producer_step: StepIdx(0),
            producer_name: "p",
            consumer_step: Some(StepIdx(1)),
            consumer_name: Some("c"),
            branch: BranchIdx(0),
            metrics: Arc::clone(&m),
            depth_source: Some(Arc::clone(&transport) as Arc<dyn BoundedQueueHandle>),
            reorder_depth: Some(Arc::clone(&stage) as Arc<dyn ReorderCapHandle>),
        };
        for _ in 0..10 {
            sample_once(std::slice::from_ref(&edge));
        }
        let s = m.snapshot();
        assert_eq!(s.depth_samples, 10);
        assert!(
            (s.mean_occupancy - 0.4).abs() < 0.05,
            "occupancy reflects the reorder stash, got {}",
            s.mean_occupancy
        );
        assert!(
            (s.mean_occupancy_bytes - 400.0).abs() < 1.0,
            "byte mean equals the stashed bytes, got {}",
            s.mean_occupancy_bytes
        );
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
        assert!(header.contains("producer__consumer#0b0.depth"), "per-edge depth column");
        let rows: Vec<&str> = lines.collect();
        assert!(!rows.is_empty(), "at least one data row");
        // First field of a data row is a monotonic t_ms integer.
        let first_t: u128 = rows[0].split('\t').next().unwrap().parse().expect("t_ms is an int");
        let last_t: u128 = rows.last().unwrap().split('\t').next().unwrap().parse().unwrap();
        assert!(last_t >= first_t, "t_ms is monotonic");
    }

    #[test]
    fn timeline_tsv_marks_unsampled_edge_na() {
        // A count/unbounded edge (no depth_source) is unsampled: its depth column
        // must read `NA`, not `0.000` (which would misread as an empty byte edge).
        let m = EdgeMetrics::new();
        let edges = vec![edge_over(Arc::clone(&m), None)];
        let dir = std::env::temp_dir();
        let path = dir.join(format!("fgumi-trace-na-{}.tsv", std::process::id()));
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
        let _header = lines.next().expect("header row");
        let row = lines.next().expect("at least one data row");
        // Columns: t_ms, <edge>.depth, <edge>.pushed, <edge>.popped.
        let depth = row.split('\t').nth(1).expect("depth column");
        assert_eq!(depth, "NA", "unsampled edge's depth column is NA, row: {row}");
    }
}
