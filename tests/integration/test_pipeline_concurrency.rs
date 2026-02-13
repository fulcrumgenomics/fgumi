//! Concurrency tests for the unified pipeline.
//!
//! These tests verify thread-safety, order preservation, data integrity,
//! and error propagation under multi-threaded pipeline execution.
//!
//! # Test categories
//!
//! - **Quick concurrency tests** (Part 3): Always run in CI via `cargo ci-test`.
//! - **Property-based tests** (Part 3.7): Proptest-based queue correctness.
//! - **Stress tests** (Part 4): Behind `#[cfg(feature = "stress-tests")]`, run with
//!   `cargo nextest run --features stress-tests`.
//!
//! # Convention for known bugs
//!
//! Tests that expose known concurrency bugs are annotated with:
//! ```rust,ignore
//! #[ignore = "Known bug: <brief description> — see <file>:<line>"]
//! ```
//! Each such test has a doc comment explaining the root cause and fix strategy.
//! When a bug is fixed, the annotation is removed — turning the test green.

#![allow(clippy::cast_possible_truncation)]

use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs::File;
use std::io::{self, BufReader};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family};
use fgumi_lib::grouper::SingleRecordGrouper;
use fgumi_lib::unified_pipeline::{
    BamPipelineConfig, DecodedRecord, Grouper, MemoryBoundedQueue, OrderedQueue,
    run_bam_pipeline_with_grouper, serialize_bam_record_into,
};

// ============================================================================
// Test Helpers
// ============================================================================

/// Create a test BAM file with the specified UMI families.
///
/// Each family is a `(umi, depth, base_name, sequence)` tuple.
/// Returns the input records for verification.
fn create_test_bam_with_families(
    path: &std::path::Path,
    header: &Header,
    families: &[(&str, usize, &str, &str)],
) -> Vec<RecordBuf> {
    let mut writer = bam::io::Writer::new(File::create(path).expect("Failed to create BAM file"));
    writer.write_header(header).expect("Failed to write header");

    let mut records = Vec::new();
    for &(umi, depth, base_name, sequence) in families {
        let family = create_umi_family(umi, depth, base_name, sequence, 30);
        for record in &family {
            writer.write_alignment_record(header, record).expect("Failed to write record");
            records.push(record.clone());
        }
    }

    writer.finish(header).expect("Failed to finish BAM");
    records
}

/// Create a large test BAM with many records across many UMI families.
///
/// Creates `num_families` families, each with `depth` records.
/// UMI sequences are auto-generated as 8-base strings.
fn create_large_test_bam(
    path: &std::path::Path,
    header: &Header,
    num_families: usize,
    depth: usize,
) -> Vec<RecordBuf> {
    let bases = [b'A', b'C', b'G', b'T'];
    let families: Vec<_> = (0..num_families)
        .map(|i| {
            // Generate unique 8-base UMI from index
            let umi: String = (0..8).map(|j| bases[(i >> (j * 2)) & 3] as char).collect();
            let name = format!("fam{i}");
            (umi, depth, name, "ACGTACGTACGTACGT".to_string())
        })
        .collect();

    let family_refs: Vec<(&str, usize, &str, &str)> = families
        .iter()
        .map(|(umi, depth, name, seq)| (umi.as_str(), *depth, name.as_str(), seq.as_str()))
        .collect();

    create_test_bam_with_families(path, header, &family_refs)
}

/// Read all records from a BAM file, returning them as `RecordBuf`.
fn read_bam_records(path: &std::path::Path) -> Vec<RecordBuf> {
    let file = File::open(path).expect("Failed to open BAM file");
    let mut reader = bam::io::Reader::new(BufReader::new(file));
    let header = reader.read_header().expect("Failed to read header");

    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        records.push(result.expect("Failed to read record"));
    }
    records
}

/// Run the BAM pipeline and return `(records_processed, output_records)`.
///
/// Uses a pass-through `SingleRecordGrouper` — each input record becomes
/// one output record. This tests the full pipeline without domain logic.
fn run_passthrough_pipeline(
    input_path: &std::path::Path,
    output_path: &std::path::Path,
    num_threads: usize,
) -> io::Result<(u64, Vec<RecordBuf>)> {
    let config = BamPipelineConfig::new(num_threads, 6);
    let count = run_bam_pipeline_with_grouper(
        config,
        input_path,
        output_path,
        |_header: &Header| Box::new(SingleRecordGrouper::new()),
        |record: RecordBuf| Ok(record),
        |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
            serialize_bam_record_into(&record, header, output)
        },
    )?;

    let output_records = read_bam_records(output_path);
    Ok((count, output_records))
}

/// Run the BAM pipeline with a custom configuration.
fn run_passthrough_pipeline_with_config(
    input_path: &std::path::Path,
    output_path: &std::path::Path,
    config: BamPipelineConfig,
) -> io::Result<(u64, Vec<RecordBuf>)> {
    let count = run_bam_pipeline_with_grouper(
        config,
        input_path,
        output_path,
        |_header: &Header| Box::new(SingleRecordGrouper::new()),
        |record: RecordBuf| Ok(record),
        |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
            serialize_bam_record_into(&record, header, output)
        },
    )?;

    let output_records = read_bam_records(output_path);
    Ok((count, output_records))
}

/// Extract record names from a list of records, preserving order.
fn record_names(records: &[RecordBuf]) -> Vec<Option<String>> {
    records
        .iter()
        .map(|r| r.name().map(|n| String::from_utf8_lossy(n.as_ref()).to_string()))
        .collect()
}

/// Extract sorted record names from a list of records.
fn sorted_record_names(records: &[RecordBuf]) -> Vec<Option<String>> {
    let mut names = record_names(records);
    names.sort();
    names
}

/// Run a closure with a timeout. Returns Err if the closure doesn't complete in time.
fn run_with_timeout<F, T>(timeout: Duration, f: F) -> Result<T, String>
where
    F: FnOnce() -> T + Send + 'static,
    T: Send + 'static,
{
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let result = f();
        let _ = tx.send(result);
    });
    rx.recv_timeout(timeout).map_err(|_| format!("Operation timed out after {timeout:?}"))
}

// ============================================================================
// Part 3: Quick Concurrency Tests (always run in CI)
// ============================================================================

// ---------- 3.1 Multi-Threaded BAM Pipeline Correctness ----------

/// Verify that the BAM pipeline preserves all records with multiple thread counts.
///
/// Runs the pipeline with 1, 2, 4, and 8 threads and verifies:
/// - Exact same record count in output for all thread counts
/// - Exact same record names (sorted comparison) for all thread counts
///
/// Documents: overall pipeline thread-safety, `is_complete()` correctness (bam.rs:699-723),
/// held-items pattern (base.rs:4229-4231).
#[test]
fn test_bam_pipeline_multithreaded_record_preservation() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let header = create_minimal_header("chr1", 10000);
    let input_records = create_large_test_bam(&input_bam, &header, 50, 20); // 1000 records

    let input_names = sorted_record_names(&input_records);
    let mut results = Vec::new();

    for num_threads in [1, 2, 4, 8] {
        let output_bam = temp_dir.path().join(format!("output_{num_threads}.bam"));
        let (count, output_records) =
            run_passthrough_pipeline(&input_bam, &output_bam, num_threads)
                .unwrap_or_else(|e| panic!("Pipeline with {num_threads} threads failed: {e}"));

        assert_eq!(
            count as usize,
            input_records.len(),
            "Thread count {num_threads}: reported count mismatch"
        );
        assert_eq!(
            output_records.len(),
            input_records.len(),
            "Thread count {num_threads}: output record count mismatch \
             (expected {}, got {})",
            input_records.len(),
            output_records.len()
        );

        let output_names = sorted_record_names(&output_records);
        assert_eq!(input_names, output_names, "Thread count {num_threads}: record names mismatch");

        results.push((num_threads, output_names));
    }

    // Cross-check: all thread counts produce the same output
    for i in 1..results.len() {
        assert_eq!(
            results[0].1, results[i].1,
            "Thread counts {} and {} produced different output",
            results[0].0, results[i].0
        );
    }
}

/// Verify deterministic output: running the same input multiple times with
/// the same thread count produces identical output.
///
/// Documents: reorder buffer correctness, serial number assignment.
#[test]
fn test_bam_pipeline_multithreaded_determinism() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let header = create_minimal_header("chr1", 10000);
    let _input_records = create_large_test_bam(&input_bam, &header, 30, 10); // 300 records

    let num_threads = 4;
    let mut all_names = Vec::new();

    for run in 0..3 {
        let output_bam = temp_dir.path().join(format!("output_run{run}.bam"));
        let (_, output_records) = run_passthrough_pipeline(&input_bam, &output_bam, num_threads)
            .unwrap_or_else(|e| panic!("Run {run} failed: {e}"));

        let names = record_names(&output_records);
        all_names.push(names);
    }

    for i in 1..all_names.len() {
        assert_eq!(all_names[0], all_names[i], "Run 0 and run {i} produced different output order");
    }
}

// ---------- 3.2 Order Preservation ----------

/// Verify that the output record order matches input order.
///
/// Documents: reorder buffer correctness, serial assignment (bam.rs:1968, bam.rs:2220).
#[test]
fn test_bam_pipeline_output_order_matches_input() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let header = create_minimal_header("chr1", 10000);

    // Create records with unique names in a specific order
    let input_records = create_large_test_bam(&input_bam, &header, 30, 5); // 150 records
    let input_names: Vec<_> = input_records
        .iter()
        .map(|r| r.name().map(|n| String::from_utf8_lossy(n.as_ref()).to_string()))
        .collect();

    for num_threads in [1, 2, 4] {
        let output_bam = temp_dir.path().join(format!("output_{num_threads}.bam"));
        let (_, output_records) = run_passthrough_pipeline(&input_bam, &output_bam, num_threads)
            .unwrap_or_else(|e| panic!("Pipeline with {num_threads} threads failed: {e}"));

        let output_names: Vec<_> = output_records
            .iter()
            .map(|r| r.name().map(|n| String::from_utf8_lossy(n.as_ref()).to_string()))
            .collect();

        assert_eq!(
            input_names, output_names,
            "Thread count {num_threads}: output order does not match input order"
        );
    }
}

/// Test that `OrderedQueue` maintains strict serial order under concurrent access.
///
/// Spawns N producer threads that insert items with assigned serials.
/// A consumer thread pops items and verifies they come out in exact serial order.
///
/// Documents: `OrderedQueue::try_pop_next()` ordering guarantee (queue.rs:340-357).
#[test]
fn test_reorder_buffer_concurrent_insert_pop() {
    let num_producers = 4;
    let items_per_producer = 250;
    let total_items = num_producers * items_per_producer;

    let queue = std::sync::Arc::new(OrderedQueue::<u64>::new(u64::MAX));

    // Assign serial ranges to each producer
    let producer_handles: Vec<_> = (0..num_producers)
        .map(|p| {
            let queue = queue.clone();
            std::thread::spawn(move || {
                for i in 0..items_per_producer {
                    let serial = (p * items_per_producer + i) as u64;
                    let value = serial * 10; // distinguishable value
                    // Retry insertion until success (backpressure may reject)
                    loop {
                        match queue.insert(serial, value, 8) {
                            Ok(()) => break,
                            Err(_) => std::thread::yield_now(),
                        }
                    }
                }
            })
        })
        .collect();

    // Consumer thread: pop items in order
    let consumer_queue = queue.clone();
    let consumer = std::thread::spawn(move || {
        let mut received = Vec::with_capacity(total_items);
        while received.len() < total_items {
            if let Some((value, _)) = consumer_queue.try_pop_next() {
                received.push(value);
            } else {
                std::thread::yield_now();
            }
        }
        received
    });

    for handle in producer_handles {
        handle.join().expect("Producer thread panicked");
    }

    let received = consumer.join().expect("Consumer thread panicked");
    assert_eq!(received.len(), total_items);

    // Verify strict serial order
    for (i, &value) in received.iter().enumerate() {
        let expected = (i as u64) * 10;
        assert_eq!(value, expected, "Item at position {i}: expected {expected}, got {value}");
    }
}

// ---------- 3.3 Dropped Items / Data Loss ----------

/// Verify no records are dropped under backpressure (slow process function).
///
/// Uses a process function that simulates slow work to cause backpressure.
/// Verifies all records appear in the output.
///
/// Documents: held-items pattern (base.rs:4229-4231), `is_complete` guard.
#[test]
fn test_bam_pipeline_no_dropped_records_under_backpressure() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let header = create_minimal_header("chr1", 10000);
    let input_records = create_large_test_bam(&input_bam, &header, 50, 10); // 500 records

    // Use a process function that simulates slow work
    let process_count = std::sync::Arc::new(AtomicU64::new(0));
    let process_count_clone = process_count.clone();

    let config = BamPipelineConfig::new(4, 6);
    let result = run_bam_pipeline_with_grouper(
        config,
        &input_bam,
        &output_bam,
        |_header: &Header| Box::new(SingleRecordGrouper::new()),
        move |record: RecordBuf| {
            process_count_clone.fetch_add(1, Ordering::Relaxed);
            // Brief yield to simulate work without excessive delay
            std::thread::yield_now();
            Ok(record)
        },
        |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
            serialize_bam_record_into(&record, header, output)
        },
    );

    assert!(result.is_ok(), "Pipeline failed: {result:?}");
    let count = result.unwrap();
    assert_eq!(count as usize, input_records.len(), "Reported count mismatch");

    let output_records = read_bam_records(&output_bam);
    assert_eq!(
        output_records.len(),
        input_records.len(),
        "Output record count mismatch under backpressure \
         (expected {}, got {})",
        input_records.len(),
        output_records.len()
    );

    let input_names = sorted_record_names(&input_records);
    let output_names = sorted_record_names(&output_records);
    assert_eq!(input_names, output_names, "Record names mismatch under backpressure");
}

/// Verify that held items are properly flushed on pipeline shutdown.
///
/// Configures a pipeline with small queue capacity to ensure backpressure
/// causes held items. Verifies `validate_completion()` passes and output
/// count matches input.
///
/// Documents: held-items pattern, `validate_completion()` (base.rs).
#[test]
fn test_pipeline_held_items_flushed_on_shutdown() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let header = create_minimal_header("chr1", 10000);
    let input_records = create_large_test_bam(&input_bam, &header, 20, 15); // 300 records

    // Small queue capacity to force held items
    let mut config = BamPipelineConfig::new(4, 6);
    config.pipeline.queue_capacity = 8;

    let result = run_passthrough_pipeline_with_config(&input_bam, &output_bam, config);
    assert!(result.is_ok(), "Pipeline failed: {result:?}");

    let (count, output_records) = result.unwrap();
    assert_eq!(count as usize, input_records.len(), "Reported count mismatch with small queues");
    assert_eq!(
        output_records.len(),
        input_records.len(),
        "Output record count mismatch with small queues"
    );
}

// ---------- 3.4 Error Propagation (Multi-Threaded) ----------

/// Verify that a process function error is correctly propagated in multi-threaded mode.
///
/// Documents: error flag (base.rs:1747-1758), `catch_unwind` (bam.rs:3193-3212).
#[test]
fn test_error_propagation_process_fn_multithreaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let header = create_minimal_header("chr1", 10000);
    create_large_test_bam(&input_bam, &header, 20, 10); // 200 records

    let count = std::sync::Arc::new(AtomicU64::new(0));
    let count_clone = count.clone();

    let config = BamPipelineConfig::new(4, 6);
    let result = run_bam_pipeline_with_grouper(
        config,
        &input_bam,
        &output_bam,
        |_header: &Header| Box::new(SingleRecordGrouper::new()),
        move |record: RecordBuf| {
            let n = count_clone.fetch_add(1, Ordering::Relaxed);
            if n == 50 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Injected error at record 50",
                ));
            }
            Ok(record)
        },
        |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
            serialize_bam_record_into(&record, header, output)
        },
    );

    assert!(result.is_err(), "Pipeline should have returned an error");
    let err = result.unwrap_err();
    assert!(
        err.to_string().contains("Injected error") || err.to_string().contains("panicked"),
        "Expected injected error, got: {err}"
    );
}

/// Verify that a serialize function error is correctly propagated in multi-threaded mode.
#[test]
fn test_error_propagation_serialize_fn_multithreaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let header = create_minimal_header("chr1", 10000);
    create_large_test_bam(&input_bam, &header, 20, 10);

    let count = std::sync::Arc::new(AtomicU64::new(0));
    let count_clone = count.clone();

    let config = BamPipelineConfig::new(4, 6);
    let result = run_bam_pipeline_with_grouper(
        config,
        &input_bam,
        &output_bam,
        |_header: &Header| Box::new(SingleRecordGrouper::new()),
        |record: RecordBuf| Ok(record),
        move |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
            let n = count_clone.fetch_add(1, Ordering::Relaxed);
            if n == 50 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Serialize error at record 50",
                ));
            }
            serialize_bam_record_into(&record, header, output)
        },
    );

    assert!(result.is_err(), "Pipeline should have returned an error");
}

/// A grouper that fails after processing a configurable number of records.
struct FailAfterNGrouper {
    count: usize,
    fail_after: usize,
}

impl FailAfterNGrouper {
    fn new(fail_after: usize) -> Self {
        Self { count: 0, fail_after }
    }
}

impl Grouper for FailAfterNGrouper {
    type Group = RecordBuf;

    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>> {
        self.count += records.len();
        if self.count >= self.fail_after {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Grouper error after {} records", self.count),
            ));
        }
        records
            .into_iter()
            .map(|d| {
                d.into_record().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "Expected parsed record")
                })
            })
            .collect()
    }

    fn finish(&mut self) -> io::Result<Option<Self::Group>> {
        Ok(None)
    }

    fn has_pending(&self) -> bool {
        false
    }
}

/// Verify that a grouper error is correctly propagated in multi-threaded mode.
#[test]
fn test_error_propagation_grouper_multithreaded() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let header = create_minimal_header("chr1", 10000);
    create_large_test_bam(&input_bam, &header, 20, 10);

    let config = BamPipelineConfig::new(4, 6);
    let result = run_bam_pipeline_with_grouper(
        config,
        &input_bam,
        &output_bam,
        |_header: &Header| Box::new(FailAfterNGrouper::new(50)),
        |record: RecordBuf| Ok(record),
        |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
            serialize_bam_record_into(&record, header, output)
        },
    );

    assert!(result.is_err(), "Pipeline should have returned an error from grouper");
}

/// Verify that errors at various pipeline stages don't cause deadlocks.
///
/// Documents: error flag propagation (base.rs:1747), deadlock timeout (deadlock.rs).
#[test]
fn test_error_does_not_cause_hang() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let header = create_minimal_header("chr1", 10000);
    create_large_test_bam(&input_bam, &header, 20, 10);

    // Test error in process function doesn't hang
    let output_bam = temp_dir.path().join("output_process.bam");
    let input_clone = input_bam.clone();
    let result = run_with_timeout(Duration::from_secs(30), move || {
        let count = std::sync::Arc::new(AtomicU64::new(0));
        let count_clone = count.clone();
        let config = BamPipelineConfig::new(4, 6);
        run_bam_pipeline_with_grouper(
            config,
            &input_clone,
            &output_bam,
            |_header: &Header| Box::new(SingleRecordGrouper::new()),
            move |record: RecordBuf| {
                let n = count_clone.fetch_add(1, Ordering::Relaxed);
                if n == 30 {
                    return Err(io::Error::other("Injected error"));
                }
                Ok(record)
            },
            |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
                serialize_bam_record_into(&record, header, output)
            },
        )
    });

    assert!(result.is_ok(), "Pipeline hung on process error: {result:?}");
    assert!(result.unwrap().is_err(), "Pipeline should have returned an error");

    // Test error in grouper doesn't hang
    let output_bam2 = temp_dir.path().join("output_grouper.bam");
    let input_clone2 = input_bam.clone();
    let result2 = run_with_timeout(Duration::from_secs(30), move || {
        let config = BamPipelineConfig::new(4, 6);
        run_bam_pipeline_with_grouper(
            config,
            &input_clone2,
            &output_bam2,
            |_header: &Header| Box::new(FailAfterNGrouper::new(30)),
            |record: RecordBuf| Ok(record),
            |record: RecordBuf, header: &Header, output: &mut Vec<u8>| {
                serialize_bam_record_into(&record, header, output)
            },
        )
    });

    assert!(result2.is_ok(), "Pipeline hung on grouper error: {result2:?}");
    assert!(result2.unwrap().is_err(), "Pipeline should have returned a grouper error");
}

// ---------- 3.5 MemoryBoundedQueue Concurrency ----------

/// Verify `MemoryBoundedQueue` concurrent push/pop doesn't lose items.
///
/// Spawns N producer threads pushing items and N consumer threads popping.
/// Verifies: `items_popped` == `items_pushed`, `current_bytes` returns to 0.
///
/// Known bug: `push()` calls `inner.push()` (making the item visible/poppable)
/// before `fetch_add` updates `current_bytes`. A concurrent `pop()` can dequeue
/// the item and `fetch_sub` its size from `current_bytes` before the corresponding
/// `fetch_add` runs, underflowing `current_bytes` to near `u64::MAX`. Then the
/// `fetch_add` return value + `heap_size` overflows u64 in debug mode.
///
/// Root cause: queue.rs:111-113 — item visibility precedes size accounting.
/// Fix: Either (a) `fetch_add` before `inner.push()` (requires rollback on push
/// failure), or (b) use `wrapping_add` and `saturating_sub` for the counter.
#[test]
#[ignore = "Known bug: push/pop accounting race — queue.rs:111-113 item visible before size counted"]
fn test_memory_bounded_queue_concurrent_push_pop() {
    let num_producers = 4;
    let num_consumers = 4;
    let items_per_producer = 500;
    let total_items = num_producers * items_per_producer;

    // Large capacity and limit so we focus on concurrency, not backpressure
    let queue = std::sync::Arc::new(MemoryBoundedQueue::<u64>::new(total_items * 2, u64::MAX));
    let items_pushed = std::sync::Arc::new(AtomicU64::new(0));
    let items_popped = std::sync::Arc::new(AtomicU64::new(0));
    let done = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));

    // Producers
    let producer_handles: Vec<_> = (0..num_producers)
        .map(|p| {
            let queue = queue.clone();
            let items_pushed = items_pushed.clone();
            std::thread::spawn(move || {
                for i in 0..items_per_producer {
                    let value = (p * items_per_producer + i) as u64;
                    loop {
                        match queue.push(value, 8) {
                            Ok(()) => {
                                items_pushed.fetch_add(1, Ordering::Relaxed);
                                break;
                            }
                            Err(_) => std::thread::yield_now(),
                        }
                    }
                }
            })
        })
        .collect();

    // Consumers
    let consumer_handles: Vec<_> = (0..num_consumers)
        .map(|_| {
            let queue = queue.clone();
            let items_popped = items_popped.clone();
            let done = done.clone();
            std::thread::spawn(move || {
                loop {
                    if queue.pop().is_some() {
                        items_popped.fetch_add(1, Ordering::Relaxed);
                    } else if done.load(Ordering::Acquire) {
                        // Drain remaining
                        while queue.pop().is_some() {
                            items_popped.fetch_add(1, Ordering::Relaxed);
                        }
                        break;
                    } else {
                        std::thread::yield_now();
                    }
                }
            })
        })
        .collect();

    for handle in producer_handles {
        handle.join().expect("Producer panicked");
    }
    done.store(true, Ordering::Release);

    for handle in consumer_handles {
        handle.join().expect("Consumer panicked");
    }

    assert_eq!(
        items_pushed.load(Ordering::Relaxed),
        total_items as u64,
        "Not all items were pushed"
    );
    assert_eq!(
        items_popped.load(Ordering::Relaxed),
        total_items as u64,
        "Not all items were popped"
    );
    assert_eq!(queue.current_bytes(), 0, "current_bytes should return to 0 after draining");
}

/// Verify that `MemoryBoundedQueue` TOCTOU overshoot is bounded.
///
/// Multiple threads push simultaneously with a small limit. The TOCTOU race
/// in `push()` (queue.rs:103-108) means `current_bytes` can transiently exceed
/// the limit. This test verifies the overshoot is bounded by `N * max_item_size`.
///
/// Documents: TOCTOU race in queue.rs:103-108 — design tradeoff, not a bug.
#[test]
fn test_memory_bounded_queue_push_overshoot_bounded() {
    let item_size = 100usize;
    let limit = 500u64;
    let num_threads = 8;
    // With 8 threads, max overshoot is (8-1) * 100 = 700, so max is 500 + 700 = 1200
    let max_expected = limit + (num_threads as u64 - 1) * item_size as u64;

    // Large capacity so ArrayQueue doesn't reject
    let queue = std::sync::Arc::new(MemoryBoundedQueue::<Vec<u8>>::new(1000, limit));
    let peak_observed = std::sync::Arc::new(AtomicU64::new(0));

    let barrier = std::sync::Arc::new(std::sync::Barrier::new(num_threads));

    let handles: Vec<_> = (0..num_threads)
        .map(|_| {
            let queue = queue.clone();
            let peak_observed = peak_observed.clone();
            let barrier = barrier.clone();
            std::thread::spawn(move || {
                barrier.wait(); // Synchronize start for maximum contention
                for _ in 0..100 {
                    let item = vec![0u8; item_size];
                    if queue.push(item, item_size).is_ok() {
                        // Record peak
                        let current = queue.current_bytes();
                        peak_observed.fetch_max(current, Ordering::Relaxed);
                    }
                }
            })
        })
        .collect();

    for handle in handles {
        handle.join().expect("Thread panicked");
    }

    let peak = peak_observed.load(Ordering::Relaxed);
    assert!(
        peak <= max_expected,
        "Peak bytes {peak} exceeds expected maximum {max_expected} \
         (limit={limit}, threads={num_threads}, item_size={item_size})"
    );
}

// ---------- 3.6 push_batch Panic Exposure ----------

/// Test that the Group step's push doesn't panic under contention.
///
/// Known bug: `bam.rs:2225` — `is_full()` check is not atomic with `push()`.
/// Multiple threads can pass `is_full()` simultaneously, then one fails the
/// `push()` because the `ArrayQueue` capacity was exhausted by another thread.
///
/// Root cause: bam.rs:2225 — `unwrap_or_else(panic!)` after non-atomic `is_full()` check.
/// Fix: Replace `unwrap_or_else(panic!)` with held-batch recovery or retry loop.
#[test]
fn test_group_step_push_does_not_panic_under_contention() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    let output_bam = temp_dir.path().join("output.bam");
    let header = create_minimal_header("chr1", 10000);
    let input_records = create_large_test_bam(&input_bam, &header, 100, 10); // 1000 records

    // Small queue capacity + many threads to maximize Group step contention
    let mut config = BamPipelineConfig::new(8, 6);
    config.pipeline.queue_capacity = 8;

    let result = run_passthrough_pipeline_with_config(&input_bam, &output_bam, config);

    assert!(result.is_ok(), "Pipeline panicked or failed: {result:?}");
    let (count, output_records) = result.unwrap();
    assert_eq!(count as usize, input_records.len());
    assert_eq!(output_records.len(), input_records.len());
}

// ---------- 3.7 Property-Based Tests ----------

mod proptest_tests {
    use super::*;
    use proptest::prelude::*;

    // Property: Reorder buffer preserves all items in serial order regardless of
    // insertion order.
    //
    // Documents: `OrderedQueue` ordering guarantee (queue.rs:340-357).
    proptest! {
        #[test]
        fn proptest_reorder_buffer_preserves_all_items(
            // Generate a random permutation of 0..n where n is 1..200
            n in 1usize..200,
            seed in any::<u64>(),
        ) {
            use rand::seq::SliceRandom;
            use rand::SeedableRng;

            let queue = OrderedQueue::<u64>::new(u64::MAX);

            // Create shuffled insertion order
            let mut serials: Vec<u64> = (0..n as u64).collect();
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            serials.shuffle(&mut rng);

            // Insert in random order
            for serial in serials {
                queue.insert(serial, serial * 7, 8).expect("insert should succeed");
            }

            // Pop should come out in strict serial order
            let mut popped = Vec::with_capacity(n);
            while let Some((value, _)) = queue.try_pop_next() {
                popped.push(value);
            }

            prop_assert_eq!(popped.len(), n, "Expected all items to be popped");
            for (i, &value) in popped.iter().enumerate() {
                prop_assert_eq!(value, (i as u64) * 7, "Item at position {} out of order", i);
            }
        }
    }

    // Property: OrderedQueue never drops items regardless of concurrent insert/pop patterns.
    //
    // Simulates interleaved insert and pop operations and verifies all items are eventually
    // received in serial order.
    proptest! {
        #[test]
        fn proptest_ordered_queue_never_drops_items(
            n in 1usize..100,
            pop_frequency in 1usize..5,  // pop every N insertions
            seed in any::<u64>(),
        ) {
            use rand::seq::SliceRandom;
            use rand::SeedableRng;

            let queue = OrderedQueue::<u64>::new(u64::MAX);

            let mut serials: Vec<u64> = (0..n as u64).collect();
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            serials.shuffle(&mut rng);

            let mut popped = Vec::new();

            // Interleave inserts and pops
            for (idx, serial) in serials.iter().enumerate() {
                queue.insert(*serial, *serial, 8).expect("insert should succeed");

                // Periodically try to pop
                if idx % pop_frequency == 0 {
                    while let Some((value, _)) = queue.try_pop_next() {
                        popped.push(value);
                    }
                }
            }

            // Drain remaining
            while let Some((value, _)) = queue.try_pop_next() {
                popped.push(value);
            }

            prop_assert_eq!(popped.len(), n, "Expected all items to be popped");
            for (i, &value) in popped.iter().enumerate() {
                prop_assert_eq!(value, i as u64, "Item at position {} out of order", i);
            }
        }
    }
}

// ============================================================================
// Part 4: Stress Tests (behind feature flag)
// ============================================================================

#[cfg(feature = "stress-tests")]
mod stress_tests {
    use super::*;

    // ---------- 4.1 Adversarial Configuration Stress ----------

    /// Stress: tiny queue capacity with many threads and many records.
    ///
    /// Queue capacity: 8, 8 threads, 1000+ records, 10 iterations.
    /// Verifies: all records preserved, no hangs (30s timeout per iteration).
    ///
    /// Documents: backpressure paths, held-items, deadlock avoidance.
    #[test]
    fn test_stress_bam_pipeline_tiny_queues() {
        for iteration in 0..10 {
            let temp_dir = TempDir::new().unwrap();
            let input_bam = temp_dir.path().join("input.bam");
            let output_bam = temp_dir.path().join("output.bam");
            let header = create_minimal_header("chr1", 10000);
            let input_records = create_large_test_bam(&input_bam, &header, 50, 20);

            let input_clone = input_bam.clone();
            let output_clone = output_bam.clone();
            let expected_count = input_records.len();

            let result = run_with_timeout(Duration::from_secs(30), move || {
                let mut config = BamPipelineConfig::new(8, 6);
                config.pipeline.queue_capacity = 8;
                run_passthrough_pipeline_with_config(&input_clone, &output_clone, config)
            });

            let result = result.unwrap_or_else(|e| panic!("Iteration {iteration}: {e}"));
            let (count, output_records) = result.unwrap_or_else(|e| {
                panic!("Iteration {iteration}: pipeline failed: {e}");
            });

            assert_eq!(count as usize, expected_count, "Iteration {iteration}: count mismatch");
            assert_eq!(
                output_records.len(),
                expected_count,
                "Iteration {iteration}: output record count mismatch"
            );
        }
    }

    /// Stress: varied thread counts produce identical output.
    ///
    /// Runs same input with 1-8 threads and verifies identical sorted output.
    #[test]
    fn test_stress_bam_pipeline_varied_thread_counts() {
        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let header = create_minimal_header("chr1", 10000);
        let input_records = create_large_test_bam(&input_bam, &header, 50, 20);

        let input_names = sorted_record_names(&input_records);

        for num_threads in 1..=8 {
            let output_bam = temp_dir.path().join(format!("output_{num_threads}.bam"));

            let result = run_with_timeout(Duration::from_secs(30), {
                let input = input_bam.clone();
                let output = output_bam.clone();
                move || run_passthrough_pipeline(&input, &output, num_threads)
            });

            let result = result.unwrap_or_else(|e| panic!("Thread count {num_threads}: {e}"));
            let (count, output_records) = result.unwrap_or_else(|e| {
                panic!("Thread count {num_threads}: pipeline failed: {e}");
            });

            assert_eq!(
                count as usize,
                input_records.len(),
                "Thread count {num_threads}: count mismatch"
            );

            let output_names = sorted_record_names(&output_records);
            assert_eq!(input_names, output_names, "Thread count {num_threads}: names mismatch");
        }
    }

    // ---------- 4.2 Deadlock Prevention Stress ----------

    /// Stress: tiny memory limits with multiple threads. Should not deadlock.
    ///
    /// Documents: deadlock detection (deadlock.rs), memory limit recovery.
    #[test]
    fn test_stress_pipeline_tiny_memory_no_deadlock() {
        for iteration in 0..10 {
            let temp_dir = TempDir::new().unwrap();
            let input_bam = temp_dir.path().join("input.bam");
            let output_bam = temp_dir.path().join("output.bam");
            let header = create_minimal_header("chr1", 10000);
            create_large_test_bam(&input_bam, &header, 10, 10); // 100 records

            let input_clone = input_bam.clone();
            let output_clone = output_bam.clone();

            let result = run_with_timeout(Duration::from_secs(30), move || {
                let mut config = BamPipelineConfig::new(4, 6);
                config.pipeline.queue_memory_limit = 64 * 1024; // 64KB
                config.pipeline.deadlock_timeout_secs = 5;
                config.pipeline.deadlock_recover_enabled = true;
                run_passthrough_pipeline_with_config(&input_clone, &output_clone, config)
            });

            assert!(result.is_ok(), "Iteration {iteration}: pipeline timed out (deadlock?)");
            // Pipeline may succeed or fail with a recoverable error — but it must not hang
        }
    }

    /// Stress: OrderedQueue with small limit and large items.
    ///
    /// Verifies items drain without deadlock even when individual items exceed the limit.
    ///
    /// Documents: deadlock avoidance in OrderedQueue (queue.rs can_accept).
    #[test]
    fn test_stress_ordered_queue_deadlock_avoidance() {
        let num_items = 100;
        let queue = std::sync::Arc::new(OrderedQueue::<Vec<u8>>::new(100)); // 100 byte limit

        let producer_queue = queue.clone();
        let producer = std::thread::spawn(move || {
            for i in 0..num_items {
                let item = vec![0u8; 1024]; // 1KB items >> 100 byte limit
                loop {
                    match producer_queue.insert(i as u64, item.clone(), 1024) {
                        Ok(()) => break,
                        Err(_) => std::thread::yield_now(),
                    }
                }
            }
        });

        let consumer_queue = queue.clone();
        let consumer = std::thread::spawn(move || {
            let mut count = 0;
            while count < num_items {
                if consumer_queue.try_pop_next().is_some() {
                    count += 1;
                } else {
                    std::thread::yield_now();
                }
            }
            count
        });

        let result = run_with_timeout(Duration::from_secs(10), move || {
            producer.join().expect("Producer panicked");
            consumer.join().expect("Consumer panicked")
        });

        assert!(result.is_ok(), "OrderedQueue deadlocked: {result:?}");
        assert_eq!(result.unwrap(), num_items);
    }

    // ---------- 4.3 Error Recovery Stress ----------

    /// Stress: error at a random record position across 20 iterations.
    ///
    /// Verifies pipeline always returns an error (never hangs), within 10s timeout.
    ///
    /// Documents: error flag propagation (base.rs:1747), catch_unwind (bam.rs:3193).
    #[test]
    fn test_stress_error_at_random_record() {
        use rand::RngExt;

        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let header = create_minimal_header("chr1", 10000);
        let total = create_large_test_bam(&input_bam, &header, 20, 10).len(); // 200

        let mut rng = rand::rng();
        for iteration in 0..20 {
            let fail_at = rng.random_range(1..total as u64);
            let output_bam = temp_dir.path().join(format!("output_{iteration}.bam"));

            let input_clone = input_bam.clone();
            let output_clone = output_bam.clone();

            let result = run_with_timeout(Duration::from_secs(10), move || {
                let count = std::sync::Arc::new(AtomicU64::new(0));
                let count_clone = count.clone();
                let config = BamPipelineConfig::new(4, 6);
                run_bam_pipeline_with_grouper(
                    config,
                    &input_clone,
                    &output_clone,
                    |_| Box::new(SingleRecordGrouper::new()),
                    move |record: RecordBuf| {
                        let n = count_clone.fetch_add(1, Ordering::Relaxed);
                        if n == fail_at {
                            return Err(io::Error::new(io::ErrorKind::Other, "random fail"));
                        }
                        Ok(record)
                    },
                    |r: RecordBuf, h: &Header, o: &mut Vec<u8>| serialize_bam_record_into(&r, h, o),
                )
            });

            assert!(result.is_ok(), "Iteration {iteration} (fail_at={fail_at}): timed out");
            assert!(
                result.unwrap().is_err(),
                "Iteration {iteration} (fail_at={fail_at}): expected error"
            );
        }
    }

    /// Stress: intermittent process failures (10% random failure rate).
    ///
    /// Verifies pipeline always returns an error (not panic), and doesn't hang.
    ///
    /// Documents: error flag propagation (base.rs:1747), catch_unwind (bam.rs:3193).
    #[test]
    fn test_stress_intermittent_process_failure() {
        for iteration in 0..10 {
            let temp_dir = TempDir::new().unwrap();
            let input_bam = temp_dir.path().join("input.bam");
            let output_bam = temp_dir.path().join("output.bam");
            let header = create_minimal_header("chr1", 10000);
            create_large_test_bam(&input_bam, &header, 50, 10); // 500 records

            let input_clone = input_bam.clone();
            let output_clone = output_bam.clone();

            let result = run_with_timeout(Duration::from_secs(10), move || {
                let config = BamPipelineConfig::new(4, 6);
                run_bam_pipeline_with_grouper(
                    config,
                    &input_clone,
                    &output_clone,
                    |_| Box::new(SingleRecordGrouper::new()),
                    |record: RecordBuf| {
                        // ~10% failure rate using a cheap hash of the record name
                        let hash = record
                            .name()
                            .map(|n| {
                                let bytes: &[u8] = n.as_ref();
                                bytes.iter().fold(0u64, |acc, &b| {
                                    acc.wrapping_mul(31).wrapping_add(u64::from(b))
                                })
                            })
                            .unwrap_or(0);
                        if hash % 10 == 0 {
                            return Err(io::Error::new(io::ErrorKind::Other, "intermittent fail"));
                        }
                        Ok(record)
                    },
                    |r: RecordBuf, h: &Header, o: &mut Vec<u8>| serialize_bam_record_into(&r, h, o),
                )
            });

            assert!(result.is_ok(), "Iteration {iteration}: timed out (deadlock?)");
            assert!(
                result.unwrap().is_err(),
                "Iteration {iteration}: expected error from intermittent failures"
            );
        }
    }

    // ---------- 4.4 Edge-Case Pipeline Inputs ----------

    /// Stress: 1 record with many threads.
    ///
    /// Verifies exactly 1 record in output and clean shutdown when most threads
    /// never do real work.
    ///
    /// Documents: shutdown logic, is_complete edge case.
    #[test]
    fn test_stress_single_record_many_threads() {
        for num_threads in [2, 4, 8] {
            let temp_dir = TempDir::new().unwrap();
            let input_bam = temp_dir.path().join("input.bam");
            let output_bam = temp_dir.path().join("output.bam");
            let header = create_minimal_header("chr1", 10000);

            let families = [("AAAAAAAA", 1usize, "single", "ACGTACGT")];
            create_test_bam_with_families(&input_bam, &header, &families);

            let result = run_passthrough_pipeline(&input_bam, &output_bam, num_threads);
            assert!(result.is_ok(), "{num_threads} threads: pipeline failed: {result:?}");
            let (count, output_records) = result.unwrap();
            assert_eq!(count, 1, "{num_threads} threads: count mismatch");
            assert_eq!(
                output_records.len(),
                1,
                "{num_threads} threads: output record count mismatch"
            );
        }
    }

    /// Stress: all records in one UMI family (one massive group).
    ///
    /// Documents: grouper backpressure when Q4 gets one huge batch.
    #[test]
    fn test_stress_all_records_one_umi_family() {
        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let output_bam = temp_dir.path().join("output.bam");
        let header = create_minimal_header("chr1", 10000);

        let families = [("AAAAAAAA", 500usize, "bigfam", "ACGTACGTACGTACGT")];
        let input_records = create_test_bam_with_families(&input_bam, &header, &families);

        let input_clone = input_bam.clone();
        let output_clone = output_bam.clone();
        let result = run_with_timeout(Duration::from_secs(30), move || {
            run_passthrough_pipeline(&input_clone, &output_clone, 4)
        });

        assert!(result.is_ok(), "Pipeline timed out with one big family");
        let (count, output_records) = result.unwrap().expect("Pipeline failed");
        assert_eq!(count as usize, input_records.len());
        assert_eq!(output_records.len(), input_records.len());
    }

    /// Stress: many tiny UMI families (1 record each).
    ///
    /// Tests high-frequency group transitions.
    ///
    /// Documents: pending_groups drain path (bam.rs:2254-2282), grouper.finish() (bam.rs:2285).
    #[test]
    fn test_stress_many_tiny_umi_families() {
        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let output_bam = temp_dir.path().join("output.bam");
        let header = create_minimal_header("chr1", 10000);

        // 1000 families with 1 record each
        let input_records = create_large_test_bam(&input_bam, &header, 1000, 1);
        let expected_len = input_records.len();

        let input_clone = input_bam.clone();
        let output_clone = output_bam.clone();
        let result = run_with_timeout(Duration::from_secs(30), move || {
            run_passthrough_pipeline(&input_clone, &output_clone, 4)
        });
        assert!(result.is_ok(), "Pipeline timed out with many tiny families");
        let (count, output_records) = result.unwrap().expect("Pipeline failed");
        assert_eq!(count as usize, expected_len);
        assert_eq!(output_records.len(), expected_len);
    }

    /// Stress: records with long sequences.
    ///
    /// Documents: boundary finding across blocks, leftover handling.
    #[test]
    fn test_stress_large_records() {
        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let output_bam = temp_dir.path().join("output.bam");
        let header = create_minimal_header("chr1", 100_000);

        // Create records with 10KB sequences
        let long_seq: String = "ACGT".repeat(2500); // 10KB

        let bases = [b'A', b'C', b'G', b'T'];
        let families: Vec<_> = (0..20)
            .map(|i| {
                let umi: String = (0..8).map(|j| bases[(i >> (j * 2)) & 3] as char).collect();
                let name = format!("bigseq{i}");
                (umi, 5usize, name, long_seq.clone())
            })
            .collect();

        let family_refs: Vec<(&str, usize, &str, &str)> = families
            .iter()
            .map(|(umi, depth, name, seq)| (umi.as_str(), *depth, name.as_str(), seq.as_str()))
            .collect();
        let input_records = create_test_bam_with_families(&input_bam, &header, &family_refs);
        let expected_len = input_records.len();
        let input_seqs: Vec<_> =
            input_records.iter().map(|r| r.sequence().as_ref().to_vec()).collect();

        let input_clone = input_bam.clone();
        let output_clone = output_bam.clone();
        let result = run_with_timeout(Duration::from_secs(30), move || {
            run_passthrough_pipeline(&input_clone, &output_clone, 4)
        });
        assert!(result.is_ok(), "Pipeline timed out with large records");
        let (count, output_records) = result.unwrap().expect("Pipeline failed with large records");
        assert_eq!(count as usize, expected_len);
        assert_eq!(output_records.len(), expected_len);

        // Verify sequences are preserved
        let mut sorted_input_seqs = input_seqs;
        sorted_input_seqs.sort();
        let mut output_seqs: Vec<_> =
            output_records.iter().map(|r| r.sequence().as_ref().to_vec()).collect();
        output_seqs.sort();
        assert_eq!(sorted_input_seqs, output_seqs, "Sequences not preserved for large records");
    }

    // ---------- 4.5 Infrastructure Interaction Stress ----------

    /// Stress: deadlock detector does not fire false positives on a slow but healthy pipeline.
    ///
    /// Documents: deadlock timeout logic (deadlock.rs), monitor thread.
    #[test]
    fn test_stress_deadlock_detector_no_false_positive() {
        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let output_bam = temp_dir.path().join("output.bam");
        let header = create_minimal_header("chr1", 10000);
        let input_records = create_large_test_bam(&input_bam, &header, 10, 10); // 100 records

        let input_clone = input_bam.clone();
        let output_clone = output_bam.clone();
        let expected_count = input_records.len();

        let result = run_with_timeout(Duration::from_secs(60), move || {
            let process_count = std::sync::Arc::new(AtomicU64::new(0));
            let process_count_clone = process_count.clone();

            let mut config = BamPipelineConfig::new(4, 6);
            config.pipeline.deadlock_timeout_secs = 10;
            config.pipeline.collect_stats = true;

            run_bam_pipeline_with_grouper(
                config,
                &input_clone,
                &output_clone,
                |_| Box::new(SingleRecordGrouper::new()),
                move |record: RecordBuf| {
                    process_count_clone.fetch_add(1, Ordering::Relaxed);
                    // Brief sleep to simulate slow processing, but not enough to trigger deadlock
                    std::thread::sleep(Duration::from_millis(1));
                    Ok(record)
                },
                |r: RecordBuf, h: &Header, o: &mut Vec<u8>| serialize_bam_record_into(&r, h, o),
            )
        });

        assert!(result.is_ok(), "Pipeline timed out");
        let pipeline_result = result.unwrap();
        assert!(
            pipeline_result.is_ok(),
            "Pipeline failed (possible false deadlock detection): {pipeline_result:?}"
        );
        assert_eq!(pipeline_result.unwrap() as usize, expected_count);
    }

    /// Stress: rapid pipeline create/teardown without leaks.
    ///
    /// Runs 50 short pipelines in sequence and verifies all complete without error.
    ///
    /// Documents: resource cleanup in join_worker_threads, finalize_pipeline.
    #[test]
    fn test_stress_rapid_pipeline_create_teardown() {
        let temp_dir = TempDir::new().unwrap();
        let input_bam = temp_dir.path().join("input.bam");
        let header = create_minimal_header("chr1", 10000);
        let input_records = create_large_test_bam(&input_bam, &header, 5, 2); // 10 records

        for iteration in 0..50 {
            let output_bam = temp_dir.path().join(format!("output_{iteration}.bam"));
            let result = run_passthrough_pipeline(&input_bam, &output_bam, 4);
            assert!(result.is_ok(), "Iteration {iteration}: pipeline failed: {result:?}");
            let (count, _) = result.unwrap();
            assert_eq!(
                count as usize,
                input_records.len(),
                "Iteration {iteration}: count mismatch"
            );
            // Clean up output to avoid accumulating files
            let _ = std::fs::remove_file(&output_bam);
        }
    }
}
