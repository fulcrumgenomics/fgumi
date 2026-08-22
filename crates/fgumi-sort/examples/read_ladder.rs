//! An additive ladder over Phase 1's input path, one stage at a time.
//!
//! Phase 1's in-process instrumentation says what each stage costs *in situ*, but
//! two things it structurally cannot say:
//!
//! 1. **What the disk can actually deliver to this process.** The 605 MB/s
//!    single-stream / 1110 MB/s four-stream figures the campaign quotes were taken
//!    with direct I/O in a different session. `raw` re-measures them buffered, from
//!    the same binary and the same boot, which is also the *actionable* ceiling —
//!    an `O_DIRECT` reader would need aligned buffers and `unsafe`.
//! 2. **What our reader does with nothing downstream.** `record_step` times the
//!    `ReadInputBlocks` step only when it *runs*; a reader sitting ineligible
//!    because `raw_input_blocks` is full (cap = `num_workers * 8` = 128) or because
//!    it is holding blocks is invisible to it. So the measured 358 MB/s is
//!    consistent with both "one buffered stream can do no better" and "the reader
//!    is fine and backpressure is throttling it" — which imply opposite fixes.
//!    `blocks` separates them: it is the same framing code with no consumer at all.
//!
//! # Reading the result — the ladder is serial, the pipeline is not
//!
//! Each rung adds one stage to the one below it, so `rung(n) - rung(n-1)` is that
//! stage's cost **on one thread, with nothing else running**. The real sort does
//! not work that way: reading is one exclusive thread, decompression is spread over
//! all 16 workers, and framing/key/push are the serial main thread. The pipeline's
//! floor is the **max over its serial resources**, not the sum of these rungs.
//!
//! Concretely: decompression is the largest single block of CPU in Phase 1
//! (231.4s of worker busy) and would dominate this ladder, but across 16 workers it
//! is ~14.5s of effective wall and nowhere near binding. A rung being expensive
//! here does **not** mean it binds in production. The `resource` column names where
//! each stage actually lands, and that is the column that decides whether a rung's
//! cost matters.
//!
//! What the ladder is good for: per-stage *ceilings*, and cross-validating the
//! in-process partition by a completely independent method.
//!
//! # Usage
//!
//! One rung per process, so no rung inherits a warm page cache from the one below
//! it — drop caches between invocations (see `scripts/read-ladder.sh`).
//!
//! ```text
//! read_ladder <rung> <bam> [--buf BYTES] [--streams N] [--limit-gb G]
//!
//! rungs, each adding one stage to the previous:
//!   raw         read() into a reused buffer                       [reader thread]
//!   blocks      + BGZF framing (read_raw_blocks)                  [reader thread]
//!   decompress  + inflate each block                              [16 workers]
//!   records     + walk the BAM records in the stream              [main thread]
//!   key         + extract the template-coordinate sort key        [main thread]
//!   push        + copy the record into an arena and keep a ref    [main thread]
//! ```

#![deny(unsafe_code)]

use std::fs::File;
use std::io::BufReader;
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use fgumi_bgzf::reader::{decompress_block_into, read_raw_blocks};
use fgumi_raw_bam::BAM_MAGIC;
use fgumi_sort::{LibraryLookup, cb_hasher, extract_template_key_inline};
use libdeflater::Decompressor;

/// The sort pipeline's own input buffer (`SORT_INPUT_BUFFER_SIZE`).
const DEFAULT_BUF_BYTES: usize = 2 * 1024 * 1024;

/// The sort pipeline's own batch size (`INPUT_READ_BATCH_SIZE`).
const BLOCKS_PER_BATCH: usize = 16;

/// Which resource a stage lands on in the real pipeline. Printed with every rung
/// because a serial ladder over a concurrent pipeline is misread by default.
fn resource_for(rung: &str) -> &'static str {
    match rung {
        "raw" | "blocks" => "reader thread (exclusive, worker 0)",
        "decompress" => "16-worker pool",
        _ => "ingest thread (serial, main)",
    }
}

/// How the `raw` rung's streams divide the file.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Pattern {
    /// Each stream owns one large contiguous span. Every stream is purely
    /// sequential, which is the friendliest case for kernel read-ahead -- and is
    /// *not* what an in-order reader can do, since stream 3 would be delivering
    /// bytes from 30 GB in while the framer still needs byte 0.
    Contiguous,
    /// Streams take turns on fixed-size chunks (stream k takes chunk k, k+N, ...).
    /// This is the access pattern an in-order parallel reader actually generates,
    /// so it is the one worth believing.
    Interleaved,
}

struct Args {
    rung: String,
    path: PathBuf,
    buf_bytes: usize,
    streams: usize,
    limit_bytes: u64,
    pattern: Pattern,
    chunk_bytes: u64,
    writer_mbps: u64,
    writer_dir: Option<PathBuf>,
}

fn parse_args() -> Result<Args> {
    let mut it = std::env::args().skip(1);
    let rung = it
        .next()
        .context("usage: read_ladder <rung> <bam> [--buf N] [--streams N] [--limit-gb G]")?;
    let path = PathBuf::from(it.next().context("missing <bam>")?);
    let mut buf_bytes = DEFAULT_BUF_BYTES;
    let mut streams = 1usize;
    let mut limit_bytes = u64::MAX;
    let mut pattern = Pattern::Contiguous;
    let mut chunk_bytes = 4 * 1024 * 1024u64;
    let mut writer_mbps = 0u64;
    let mut writer_dir = None;
    while let Some(flag) = it.next() {
        let value = it.next().with_context(|| format!("{flag} needs a value"))?;
        match flag.as_str() {
            "--buf" => buf_bytes = value.parse().context("--buf")?,
            "--streams" => streams = value.parse().context("--streams")?,
            "--limit-gb" => {
                let gb: f64 = value.parse().context("--limit-gb")?;
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                {
                    limit_bytes = (gb * 1e9) as u64;
                }
            }
            "--pattern" => {
                pattern = match value.as_str() {
                    "contiguous" => Pattern::Contiguous,
                    "interleaved" => Pattern::Interleaved,
                    other => bail!("--pattern must be contiguous|interleaved, got {other}"),
                };
            }
            "--chunk" => chunk_bytes = value.parse().context("--chunk")?,
            "--writer-mbps" => writer_mbps = value.parse().context("--writer-mbps")?,
            "--writer-dir" => writer_dir = Some(PathBuf::from(value)),
            other => bail!("unknown flag {other}"),
        }
    }
    if streams == 0 {
        bail!("--streams must be at least 1");
    }
    if chunk_bytes == 0 {
        bail!("--chunk must be positive");
    }
    if buf_bytes == 0 {
        bail!("--buf must be positive");
    }
    if writer_mbps > 0 && writer_dir.is_none() {
        bail!("--writer-mbps needs --writer-dir");
    }
    // Only `raw` spawns the concurrent Writer; the pipeline rungs share a single
    // reader loop with no writer, so accepting these flags there would silently
    // ignore them (and, unthrottled, could write a multi-GB temp file for
    // nothing). Reject them instead of pretending to honor them.
    if (writer_mbps > 0 || writer_dir.is_some()) && rung != "raw" {
        bail!("--writer-mbps/--writer-dir only apply to the `raw` rung");
    }
    Ok(Args {
        rung,
        path,
        buf_bytes,
        streams,
        limit_bytes,
        pattern,
        chunk_bytes,
        writer_mbps,
        writer_dir,
    })
}

#[derive(Default)]
struct Tally {
    bytes: u64,
    blocks: u64,
    records: u64,
    wrote_bytes: u64,
}

fn main() -> Result<()> {
    let args = parse_args()?;
    let started = Instant::now();
    let tally = match args.rung.as_str() {
        "raw" => rung_raw(&args)?,
        "blocks" | "decompress" | "records" | "key" | "push" => rung_pipeline(&args)?,
        other => bail!("unknown rung {other}"),
    };
    report(&args, &tally, started.elapsed().as_secs_f64());
    Ok(())
}

#[allow(clippy::cast_precision_loss, reason = "byte and record counts stay far below 2^52")]
fn report(args: &Args, tally: &Tally, secs: f64) {
    let mb_per_sec = if secs > 0.0 { tally.bytes as f64 / secs / 1e6 } else { 0.0 };
    let per = |n: u64| if n == 0 { 0.0 } else { secs * 1e9 / n as f64 };
    let write_mb_per_sec = if secs > 0.0 { tally.wrote_bytes as f64 / secs / 1e6 } else { 0.0 };
    let pattern = match args.pattern {
        Pattern::Contiguous => "contiguous",
        Pattern::Interleaved => "interleaved",
    };
    println!(
        "rung={rung} streams={streams} pattern={pattern} chunk={chunk} buf={buf} \
         secs={secs:.2} MBps={mb:.0} write_MBps={wmb:.0} \
         bytes={bytes} blocks={blocks} ns_per_block={npb:.0} records={records} \
         ns_per_record={npr:.1} resource=\"{res}\"",
        rung = args.rung,
        streams = args.streams,
        chunk = args.chunk_bytes,
        buf = args.buf_bytes,
        wmb = write_mb_per_sec,
        mb = mb_per_sec,
        bytes = tally.bytes,
        blocks = tally.blocks,
        npb = per(tally.blocks),
        records = tally.records,
        npr = per(tally.records),
        res = resource_for(&args.rung),
    );
}

/// `raw`: how fast the file's bytes can be pulled into this process at all.
///
/// `streams > 1` splits the file into contiguous ranges and reads each on its own
/// thread with `read_at`, which is what a parallel reader would do — the volume
/// measured 605 MB/s at one stream and 1110 MB/s at four, and nothing in the
/// current reader can exploit that.
fn rung_raw(args: &Args) -> Result<Tally> {
    let file =
        Arc::new(File::open(&args.path).with_context(|| format!("open {}", args.path.display()))?);
    let total = file.metadata()?.len().min(args.limit_bytes);
    let streams = u64::try_from(args.streams).expect("stream count fits in u64");
    let span = total.div_ceil(streams);

    // Production reads and writes the same volume at once: Phase 1 spills 52.6 GB
    // while it reads 43.1 GB. A read rate measured on an idle device is therefore
    // an upper bound that production can never see, and the gap is the whole
    // question -- the volume is provisioned near 1000 MB/s, so a 1048 MB/s read
    // plus a ~363 MB/s spill stream cannot both be served.
    let writer = args.writer_dir.as_ref().map(|dir| Writer::spawn(dir, args.writer_mbps));

    let mut handles = Vec::with_capacity(args.streams);
    for stream in 0..streams {
        let file = Arc::clone(&file);
        let buf_bytes = args.buf_bytes;
        let pattern = args.pattern;
        let chunk = args.chunk_bytes;
        let start = stream * span;
        let end = ((stream + 1) * span).min(total);
        handles.push(std::thread::spawn(move || -> Result<u64> {
            let mut buf = vec![0u8; buf_bytes];
            let mut read_total = 0u64;
            match pattern {
                Pattern::Contiguous => {
                    let mut offset = start;
                    while offset < end {
                        let want = usize::try_from((end - offset).min(buf_bytes as u64))
                            .expect("clamped to buffer size");
                        let n = file.read_at(&mut buf[..want], offset)?;
                        if n == 0 {
                            break;
                        }
                        offset += n as u64;
                        read_total += n as u64;
                    }
                }
                Pattern::Interleaved => {
                    // Stream k takes chunk k, k + N, k + 2N ... so the file is
                    // covered in order across the pool rather than in N disjoint
                    // spans. Each stream strides by `streams * chunk`, which is
                    // what makes this a different ask of kernel read-ahead.
                    let stride = chunk * streams;
                    let mut chunk_start = stream * chunk;
                    while chunk_start < total {
                        let chunk_end = (chunk_start + chunk).min(total);
                        let mut offset = chunk_start;
                        while offset < chunk_end {
                            let want = usize::try_from((chunk_end - offset).min(buf_bytes as u64))
                                .expect("clamped to buffer size");
                            let n = file.read_at(&mut buf[..want], offset)?;
                            if n == 0 {
                                break;
                            }
                            offset += n as u64;
                            read_total += n as u64;
                        }
                        chunk_start += stride;
                    }
                }
            }
            Ok(read_total)
        }));
    }

    // Join every reader before propagating a failure, so the writer is always
    // stopped and its temp file removed by `finish` -- returning through `??`
    // here would leave the writer thread running and `read_ladder_writer.tmp`
    // (written unthrottled when --writer-mbps is 0) on disk.
    let mut tally = Tally::default();
    let mut read_result = Ok(());
    for handle in handles {
        match handle.join().map_err(|_| anyhow::anyhow!("reader thread panicked"))? {
            Ok(n) => tally.bytes += n,
            Err(e) => read_result = Err(e),
        }
    }
    if let Some(writer) = writer {
        tally.wrote_bytes = writer.finish()?;
    }
    read_result?;
    Ok(tally)
}

/// A rate-limited background writer, so the read rungs can be measured against
/// the spill load production actually runs alongside them.
struct Writer {
    stop: Arc<std::sync::atomic::AtomicBool>,
    handle: std::thread::JoinHandle<Result<u64>>,
    path: PathBuf,
}

impl Writer {
    fn spawn(dir: &Path, target_mbps: u64) -> Self {
        use std::io::Write;
        let path = dir.join("read_ladder_writer.tmp");
        let stop = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let thread_stop = Arc::clone(&stop);
        let thread_path = path.clone();
        let handle = std::thread::spawn(move || -> Result<u64> {
            let mut file = File::create(&thread_path)
                .with_context(|| format!("create {}", thread_path.display()))?;
            let block = vec![0x5Au8; 4 * 1024 * 1024];
            let started = Instant::now();
            let mut written = 0u64;
            while !thread_stop.load(std::sync::atomic::Ordering::Relaxed) {
                file.write_all(&block)?;
                written += block.len() as u64;
                if target_mbps > 0 {
                    // Token bucket: hold the writer at the target rate rather than
                    // letting it race the reader for the device, which would
                    // measure contention at some arbitrary ratio instead of the
                    // one production generates.
                    #[allow(clippy::cast_precision_loss, reason = "byte totals stay below 2^52")]
                    let due = written as f64 / (target_mbps as f64 * 1e6);
                    let elapsed = started.elapsed().as_secs_f64();
                    if due > elapsed {
                        std::thread::sleep(std::time::Duration::from_secs_f64(due - elapsed));
                    }
                }
            }
            file.sync_all()?;
            Ok(written)
        });
        Self { stop, handle, path }
    }

    fn finish(self) -> Result<u64> {
        self.stop.store(true, std::sync::atomic::Ordering::Relaxed);
        let written = self.handle.join().map_err(|_| anyhow::anyhow!("writer thread panicked"))?;
        let _ = std::fs::remove_file(&self.path);
        written
    }
}

/// Everything from `blocks` up, sharing one loop so the rungs differ by exactly the
/// stage each adds and nothing else.
fn rung_pipeline(args: &Args) -> Result<Tally> {
    let rung = args.rung.as_str();
    let want_decompress = matches!(rung, "decompress" | "records" | "key" | "push");
    let want_records = matches!(rung, "records" | "key" | "push");
    let want_key = matches!(rung, "key" | "push");
    let want_push = rung == "push";

    // The header is read through the crate's own opener, so the key rung sees the
    // same library ordinals production does rather than an invented mapping.
    let lookup = if want_key {
        let (_reader, header) = fgumi_sort::open_raw_bam_record_reader_with_header(&args.path)?;
        Some(LibraryLookup::from_header(&header))
    } else {
        None
    };
    let hasher = cb_hasher();

    let file = File::open(&args.path).with_context(|| format!("open {}", args.path.display()))?;
    let mut reader = BufReader::with_capacity(args.buf_bytes, file);

    let mut tally = Tally::default();
    let mut decompressed = Vec::with_capacity(64 * 1024);
    // Reused across every block, as the worker pool does -- constructing one per
    // block would put allocator noise in the decompress rung's delta.
    let mut decompressor = Decompressor::new();
    let mut framer = RecordFramer::default();
    let mut arena: Vec<u8> = Vec::new();
    let mut refs: Vec<(u64, usize)> = Vec::new();

    loop {
        let batch = read_raw_blocks(&mut reader, BLOCKS_PER_BATCH)?;
        if batch.is_empty() {
            break;
        }
        for block in &batch {
            tally.bytes += block.len() as u64;
            tally.blocks += 1;
            if !want_decompress {
                continue;
            }
            // `decompress_block_into` *appends*. Without this the buffer accumulates
            // the whole file, which costs the decompress rung a realloc storm and
            // makes the records rung quadratic -- it re-frames every prior block,
            // which is ~2 TB of copying on a 43 GB input and eats RAM until the
            // OOM killer intervenes. This cost a workstation once; leave it here.
            decompressed.clear();
            decompress_block_into(block, &mut decompressor, &mut decompressed)?;
            if !want_records {
                continue;
            }
            framer.push(&decompressed);
            while let Some(record) = framer.next_record() {
                tally.records += 1;
                if want_key {
                    let key = extract_template_key_inline(
                        record,
                        lookup.as_ref().expect("key rung builds a lookup"),
                        None,
                        &hasher,
                    );
                    if want_push {
                        refs.push((key.primary, arena.len()));
                        arena.extend_from_slice(record);
                        // The real ingest spills at a memory bound; here the arena
                        // is reset instead, so the rung measures the copy rather
                        // than this process's willingness to hold 40 GB.
                        if arena.len() >= 512 * 1024 * 1024 {
                            arena.clear();
                            refs.clear();
                        }
                    } else {
                        std::hint::black_box(&key);
                    }
                }
            }
        }
        if tally.bytes >= args.limit_bytes {
            break;
        }
    }
    Ok(tally)
}

/// Frames BAM records out of the decompressed byte stream.
///
/// Records straddle BGZF block boundaries, so this carries the tail of one block
/// into the next — the same job `PooledInputStream` does in production, minus the
/// reorder buffer, which only exists because production decompresses out of order.
#[derive(Default)]
struct RecordFramer {
    carry: Vec<u8>,
    start: usize,
    header_skipped: bool,
}

/// The carry may hold one partial record plus one decompressed block. A BGZF
/// block is at most 64 KiB and a BAM record in any real file is far below this,
/// so exceeding the bound does not mean "unusual input" -- it means the consumer
/// has stopped draining and the buffer is growing without limit. Asserting it is
/// cheaper than discovering it through the OOM killer.
const CARRY_LIMIT: usize = 16 * 1024 * 1024;

impl RecordFramer {
    fn push(&mut self, bytes: &[u8]) {
        // Reclaim consumed bytes before growing, so the carry stays near one block
        // rather than the whole file.
        if self.start > 0 {
            self.carry.drain(..self.start);
            self.start = 0;
        }
        self.carry.extend_from_slice(bytes);
        assert!(
            self.carry.len() <= CARRY_LIMIT,
            "record framer carry reached {} bytes, past the {CARRY_LIMIT}-byte bound: \
             the consumer is not draining and this is about to exhaust memory",
            self.carry.len(),
        );
    }

    fn next_record(&mut self) -> Option<&[u8]> {
        if !self.header_skipped && !self.try_skip_header() {
            return None;
        }
        let available = &self.carry[self.start..];
        if available.len() < 4 {
            return None;
        }
        let len = u32::from_le_bytes(available[..4].try_into().ok()?) as usize;
        // A real BAM record is far below CARRY_LIMIT (a BGZF block is <= 64 KiB), so
        // an oversized length is a corrupt or misframed record, not one still
        // waiting for more data. Without this guard `available.len() < 4 + len`
        // never clears, `next_record` returns None on every call, and `push` grows
        // the carry until it asserts "the consumer is not draining" -- blaming the
        // wrong cause. Surface the real one here instead.
        assert!(
            len <= CARRY_LIMIT,
            "record framer read a record length of {len} bytes, past the {CARRY_LIMIT}-byte \
             bound: the input is corrupt or misframed",
        );
        if available.len() < 4 + len {
            return None;
        }
        let from = self.start + 4;
        self.start += 4 + len;
        Some(&self.carry[from..from + len])
    }

    /// Skip `magic | l_text | text | n_ref | (l_name, name, l_ref) * n_ref`.
    ///
    /// Returns false while the header is still incomplete, so the caller simply
    /// feeds another block.
    fn try_skip_header(&mut self) -> bool {
        let buf = &self.carry[self.start..];
        if buf.len() < 12 {
            return false;
        }
        assert_eq!(&buf[..4], BAM_MAGIC, "input is not a BAM stream");
        let l_text = u32::from_le_bytes(buf[4..8].try_into().expect("4 bytes")) as usize;
        let Some(after_text) = 8usize.checked_add(l_text) else { return false };
        if buf.len() < after_text + 4 {
            return false;
        }
        let n_ref = u32::from_le_bytes(buf[after_text..after_text + 4].try_into().expect("4 bytes"))
            as usize;
        let mut at = after_text + 4;
        for _ in 0..n_ref {
            if buf.len() < at + 4 {
                return false;
            }
            let l_name = u32::from_le_bytes(buf[at..at + 4].try_into().expect("4 bytes")) as usize;
            // Parity with `next_record`: an oversized name length is corruption, not
            // a header still arriving, so reject it with the real cause.
            assert!(
                l_name <= CARRY_LIMIT,
                "BAM header declared a reference-name length of {l_name} bytes, past the \
                 {CARRY_LIMIT}-byte bound: the input is corrupt or misframed",
            );
            // Checked like `l_text` above: a length that overflows the cursor is
            // treated as an incomplete header (feed another block) rather than
            // wrapping past the buffer.
            let Some(next_at) = at
                .checked_add(4)
                .and_then(|x| x.checked_add(l_name))
                .and_then(|x| x.checked_add(4))
            else {
                return false;
            };
            at = next_at;
            if buf.len() < at {
                return false;
            }
        }
        self.start += at;
        self.header_skipped = true;
        true
    }
}
