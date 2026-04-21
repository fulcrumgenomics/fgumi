//! FASTQ stream coordination — gzip and BGZF variants.
//!
//! This module provides two coordination stages, both
//! [`Parallelism::Sequential`]: one [`FastqPair`] instance sees every
//! input across all streams so the `BTreeMap` / per-stream queues that
//! track pending batches stay consistent.
//!
//! ## [`FastqPair`] — gzip variant
//!
//! ### Input / Output
//!
//! - Input: [`FastqParsedStream`] (records parsed per gzip chunk, tagged
//!   by `stream_idx` and `batch_num`).
//! - Output: [`TemplateBatch`] (up to `MAX_TEMPLATES_PER_EMIT` = 256
//!   `FastqTemplateV2` items per emit, with a fresh monotonic
//!   `ordinal`).
//!
//! ### Ordering guarantees
//!
//! Strictly in-order by `batch_num` — the pending map emits only when
//! the exact `next_batch` is fully populated. Single-end
//! (`num_streams == 1`) emits one template per record immediately.
//! Paired-end waits for both streams to deliver matching `batch_num`.
//! Within a batch, records are zipped position-wise after sorting
//! streams by `stream_idx` for deterministic zipping.
//!
//! ### Memory model
//!
//! Bounded by the upstream backpressure window: `pending` holds at
//! most one `batch_num` per stream in flight. A record-count mismatch
//! between streams is a hard error.
//!
//! ### Determinism
//!
//! Fully deterministic given deterministic inputs: no RNG, no
//! scheduling-dependent branching.
//!
//! ## [`FastqBlockMerge`] — BGZF variant
//!
//! ### Input / Output
//!
//! - Input: [`BlockParsed`] (per-BGZF-block parsed records with
//!   prefix/suffix fragments).
//! - Output: [`TemplateBatch`] (as above).
//!
//! ### Ordering guarantees
//!
//! Per-stream: records are drained in `block_idx` order with
//! cross-block fragment stitching. Zipping across streams emits only
//! when every stream has at least one ready record; if one stream
//! runs ahead, its surplus stays in `ready` until the others catch
//! up.
//!
//! ### Memory model
//!
//! Per-stream buffers (`pending`, `trailing`, `ready`) are bounded by
//! the upstream `ReadAhead`/decompress/parse backlog — no explicit cap
//! here.
//!
//! ### Determinism
//!
//! Fully deterministic: block stitching and position-wise zipping are
//! pure functions of the input sequence.

use std::collections::BTreeMap;

use anyhow::Result;

/// Maximum number of paired templates emitted in a single `TemplateBatch`.
///
/// Caps the fan-out of Sequential coordination stages so downstream Parallel
/// stages (`ExtractStage`, `BgzfCompress`) see enough items to distribute
/// across workers. `unified_pipeline`'s `blocks_per_batch=4` + per-block-pair
/// emit model gives ~100-300 templates/item naturally; 256 matches that
/// granularity for `pipeline`'s aggregation-based merge.
const MAX_TEMPLATES_PER_EMIT: usize = 256;

use crate::runall::engine::fastq_types::{
    FastqParsedStream, FastqTemplateV2, OwnedFastqRecord, TemplateBatch,
};
use crate::runall::engine::stage::{Parallelism, Stage};

pub struct FastqPair {
    num_streams: usize,
    /// `pending[batch_num][stream_idx]` = parsed records for that stream.
    pending: BTreeMap<u64, BTreeMap<usize, FastqParsedStream>>,
    /// Next `batch_num` expected to emit. Batches greater than this must wait
    /// until this exact batch has been populated and emitted (no gap-skipping).
    next_batch: u64,
    /// Monotonic ordinal stamped on each `TemplateBatch` emitted.
    next_ordinal: u64,
}

impl FastqPair {
    #[must_use]
    pub fn new(num_streams: usize) -> Self {
        Self { num_streams, pending: BTreeMap::new(), next_batch: 0, next_ordinal: 0 }
    }
}

impl Clone for FastqPair {
    fn clone(&self) -> Self {
        // Cloning resets per-instance ordering state: the clone is a fresh
        // Sequential worker (pool uses a Mutex so only one instance is ever
        // active at a time, but StageSource::Clone still produces copies).
        Self::new(self.num_streams)
    }
}

impl Stage for FastqPair {
    type Input = FastqParsedStream;
    type Output = TemplateBatch;

    #[tracing::instrument(name = "fastq_pair", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let batch_num = input.batch_num;
        let stream_idx = input.stream_idx;
        self.pending.entry(batch_num).or_default().insert(stream_idx, input);
        let mut templates = self.try_emit()?;
        if templates.is_empty() {
            return Ok(());
        }

        // Emit in chunks of at most MAX_TEMPLATES_PER_EMIT so downstream
        // Parallel stages see enough items to distribute across workers.
        while !templates.is_empty() {
            let take = templates.len().min(MAX_TEMPLATES_PER_EMIT);
            let chunk: Vec<FastqTemplateV2> = templates.drain(..take).collect();
            let ordinal = self.next_ordinal;
            self.next_ordinal += 1;
            out(TemplateBatch { templates: chunk, ordinal });
        }
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Sequential
    }

    fn output_memory_estimate(&self, batch: &Self::Output) -> usize {
        batch
            .templates
            .iter()
            .flat_map(|t| &t.records)
            .map(|r| r.name.len() + r.sequence.len() + r.quality.len())
            .sum()
    }

    fn name(&self) -> &'static str {
        "FastqPair"
    }
}

impl FastqPair {
    fn try_emit(&mut self) -> Result<Vec<FastqTemplateV2>> {
        let mut emitted = Vec::new();
        loop {
            // Only emit the exact next expected batch — no gap-skipping.
            let first_batch = match self.pending.get(&self.next_batch) {
                Some(streams) if streams.len() >= self.num_streams => self.next_batch,
                _ => break,
            };
            let batch = self.pending.remove(&first_batch).unwrap();
            self.next_batch += 1;

            // Sort streams by stream_idx for deterministic zipping.
            let mut streams_vec: Vec<FastqParsedStream> = batch.into_values().collect();
            streams_vec.sort_by_key(|s| s.stream_idx);

            // Validate record counts match.
            let counts: Vec<usize> = streams_vec.iter().map(|s| s.records.len()).collect();
            let min = *counts.iter().min().unwrap_or(&0);
            let max = *counts.iter().max().unwrap_or(&0);
            if min != max {
                anyhow::bail!("stream record count mismatch at batch {first_batch}: {counts:?}");
            }

            // Zip records position-wise.
            for i in 0..min {
                let records: Vec<OwnedFastqRecord> =
                    streams_vec.iter().map(|s| s.records[i].clone()).collect();
                let name = records[0].name.clone();
                emitted.push(FastqTemplateV2 { records, name });
            }
        }
        Ok(emitted)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(n: &[u8], s: &[u8], q: &[u8]) -> OwnedFastqRecord {
        OwnedFastqRecord { name: n.to_vec(), sequence: s.to_vec(), quality: q.to_vec() }
    }

    fn stream(
        batch: u64,
        idx: usize,
        records: Vec<OwnedFastqRecord>,
        is_last: bool,
    ) -> FastqParsedStream {
        FastqParsedStream { batch_num: batch, stream_idx: idx, records, is_last }
    }

    /// Run one `process` call and return the emitted batch (if any).
    /// `FastqPair` suppresses empty batches, so `None` means "waiting for more streams."
    fn process_emit(stage: &mut FastqPair, input: FastqParsedStream) -> Option<TemplateBatch> {
        let mut captured = None;
        stage
            .process(input, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        captured
    }

    #[test]
    fn test_single_end_emits_immediately() {
        let mut stage = FastqPair::new(1);
        let out = process_emit(&mut stage, stream(0, 0, vec![rec(b"r1", b"A", b"I")], false))
            .expect("single-end must emit on first record");
        assert_eq!(out.templates.len(), 1);
        assert_eq!(out.templates[0].records.len(), 1);
        assert_eq!(&out.templates[0].records[0].name, b"r1");
        assert_eq!(out.ordinal, 0, "first emit has ordinal 0");
    }

    #[test]
    fn test_paired_waits_for_both_streams() {
        let mut stage = FastqPair::new(2);

        // Stream 0 arrives first — no emit yet.
        let out0 =
            process_emit(&mut stage, stream(0, 0, vec![rec(b"r1", b"ACGT", b"IIII")], false));
        assert!(out0.is_none(), "should suppress emit while waiting for stream 1");

        // Stream 1 arrives — completes batch 0.
        let out1 =
            process_emit(&mut stage, stream(0, 1, vec![rec(b"r1", b"CGTA", b"JJJJ")], false))
                .expect("batch complete — must emit");
        assert_eq!(out1.templates.len(), 1);
        assert_eq!(out1.templates[0].records[0].sequence, b"ACGT");
        assert_eq!(out1.templates[0].records[1].sequence, b"CGTA");
        assert_eq!(out1.ordinal, 0, "first emit has ordinal 0");
    }

    #[test]
    fn test_multiple_batches_emit_in_order() {
        let mut stage = FastqPair::new(2);

        // Arrive out of order: batch 1 both streams, then batch 0 both streams.
        let out_b1_s0 = process_emit(&mut stage, stream(1, 0, vec![rec(b"r3", b"A", b"I")], false));
        assert!(out_b1_s0.is_none(), "waiting for batch 1 stream 1");

        let out_b1_s1 = process_emit(&mut stage, stream(1, 1, vec![rec(b"r3", b"B", b"J")], false));
        // batch 0 is still missing, so batch 1 must wait — no emit.
        assert!(out_b1_s1.is_none(), "batch 1 must wait for batch 0");

        // Now fulfill batch 0.
        let out_b0_s0 = process_emit(&mut stage, stream(0, 0, vec![rec(b"r1", b"C", b"K")], false));
        assert!(out_b0_s0.is_none(), "batch 0 not yet complete");

        let out_b0_s1 = process_emit(&mut stage, stream(0, 1, vec![rec(b"r1", b"D", b"L")], false))
            .expect("both batches emit now");
        assert_eq!(out_b0_s1.templates.len(), 2);
        assert_eq!(out_b0_s1.templates[0].records[0].sequence, b"C"); // batch 0 first
        assert_eq!(out_b0_s1.templates[1].records[0].sequence, b"A"); // batch 1 second
        assert_eq!(out_b0_s1.ordinal, 0, "first emit has ordinal 0");
    }

    #[test]
    fn test_record_count_mismatch_errors() {
        let mut stage = FastqPair::new(2);
        // Stream 0 arrives with 2 records — suppressed (waiting for stream 1).
        let out0 = process_emit(
            &mut stage,
            stream(0, 0, vec![rec(b"r1", b"A", b"I"), rec(b"r2", b"B", b"J")], false),
        );
        assert!(out0.is_none(), "waiting for stream 1");

        // Stream 1 arrives with 1 record — mismatch detected on emit attempt.
        let result = stage.process(stream(0, 1, vec![rec(b"r1", b"C", b"K")], false), &mut |_| {
            panic!("no output on error")
        });
        assert!(result.is_err(), "mismatched record counts should error");
        let msg = format!("{}", result.unwrap_err());
        assert!(msg.contains("mismatch"), "error message should mention mismatch");
    }

    #[test]
    fn test_multi_record_batch_zips_position_wise() {
        let mut stage = FastqPair::new(2);
        let out0 = process_emit(
            &mut stage,
            stream(0, 0, vec![rec(b"r1", b"A1", b"I1"), rec(b"r2", b"A2", b"I2")], false),
        );
        assert!(out0.is_none(), "waiting for stream 1");

        let out = process_emit(
            &mut stage,
            stream(0, 1, vec![rec(b"r1", b"B1", b"J1"), rec(b"r2", b"B2", b"J2")], false),
        )
        .expect("batch complete — must emit");
        assert_eq!(out.templates.len(), 2);
        assert_eq!(out.templates[0].records[0].sequence, b"A1");
        assert_eq!(out.templates[0].records[1].sequence, b"B1");
        assert_eq!(out.templates[1].records[0].sequence, b"A2");
        assert_eq!(out.templates[1].records[1].sequence, b"B2");
        assert_eq!(out.ordinal, 0, "first emit has ordinal 0");
    }

    #[test]
    fn test_fastq_pair_ordinal_increments_per_emit() {
        let mut stage = FastqPair::new(2);

        // Batch 0: stream 0 suppressed, stream 1 emits with ordinal 0.
        assert!(
            process_emit(&mut stage, stream(0, 0, vec![rec(b"r1", b"A", b"I")], false)).is_none()
        );
        let out1 = process_emit(&mut stage, stream(0, 1, vec![rec(b"r1", b"B", b"J")], false))
            .expect("batch 0 complete");
        assert_eq!(out1.ordinal, 0);

        // Batch 1: stream 0 suppressed, stream 1 emits with ordinal 1 (not 2 — empties don't consume ordinals).
        assert!(
            process_emit(&mut stage, stream(1, 0, vec![rec(b"r2", b"C", b"K")], false)).is_none()
        );
        let out3 = process_emit(&mut stage, stream(1, 1, vec![rec(b"r2", b"D", b"L")], false))
            .expect("batch 1 complete");
        assert_eq!(out3.ordinal, 1, "second emit has ordinal 1 (empties do not consume ordinals)");
    }
}

// ===========================================================================
// FastqBlockMerge (BGZF path)
// ===========================================================================

use std::collections::VecDeque;

use crate::runall::engine::fastq_types::BlockParsed;

struct StreamState {
    pending: BTreeMap<u64, BlockParsed>,
    next_block: u64,
    trailing: Vec<u8>,
    ready: VecDeque<OwnedFastqRecord>,
    eof: bool,
}

impl StreamState {
    fn new() -> Self {
        Self {
            pending: BTreeMap::new(),
            next_block: 0,
            trailing: Vec::new(),
            ready: VecDeque::new(),
            eof: false,
        }
    }
}

/// BGZF variant of FASTQ stream coordination.
///
/// Stitches prefix/suffix fragments across blocks (cross-block records),
/// then zips records across streams position-wise into `FastqTemplateV2`.
pub struct FastqBlockMerge {
    num_streams: usize,
    per_stream: Vec<StreamState>,
    /// Monotonic ordinal stamped on each `TemplateBatch` emitted.
    next_ordinal: u64,
}

impl FastqBlockMerge {
    #[must_use]
    pub fn new(num_streams: usize) -> Self {
        Self {
            num_streams,
            per_stream: (0..num_streams).map(|_| StreamState::new()).collect(),
            next_ordinal: 0,
        }
    }
}

impl Clone for FastqBlockMerge {
    fn clone(&self) -> Self {
        Self::new(self.num_streams)
    }
}

impl Stage for FastqBlockMerge {
    type Input = BlockParsed;
    type Output = TemplateBatch;

    #[tracing::instrument(name = "fastq_block_merge", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let s = input.stream_idx;
        anyhow::ensure!(s < self.num_streams, "stream_idx {s} >= num_streams {}", self.num_streams);
        self.per_stream[s].pending.insert(input.block_idx, input);

        for s in 0..self.num_streams {
            self.drain_stream(s)?;
        }

        let min_ready = self.per_stream.iter().map(|ss| ss.ready.len()).min().unwrap_or(0);
        if min_ready == 0 {
            return Ok(());
        }

        // Emit in chunks of at most MAX_TEMPLATES_PER_EMIT so downstream
        // Parallel stages see enough items to distribute across workers.
        let mut remaining = min_ready;
        while remaining > 0 {
            let take = remaining.min(MAX_TEMPLATES_PER_EMIT);
            let mut templates = Vec::with_capacity(take);
            for _ in 0..take {
                let records: Vec<OwnedFastqRecord> = self
                    .per_stream
                    .iter_mut()
                    .map(|ss| ss.ready.pop_front().expect("min_ready was computed from ready.len"))
                    .collect();
                let name = records[0].name.clone();
                templates.push(FastqTemplateV2 { records, name });
            }
            let ordinal = self.next_ordinal;
            self.next_ordinal += 1;
            out(TemplateBatch { templates, ordinal });
            remaining -= take;
        }
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Sequential
    }

    fn output_memory_estimate(&self, batch: &Self::Output) -> usize {
        batch
            .templates
            .iter()
            .flat_map(|t| &t.records)
            .map(|r| r.name.len() + r.sequence.len() + r.quality.len())
            .sum()
    }

    fn name(&self) -> &'static str {
        "FastqBlockMerge"
    }
}

impl FastqBlockMerge {
    fn drain_stream(&mut self, s: usize) -> Result<()> {
        loop {
            let next = self.per_stream[s].next_block;
            let Some(block) = self.per_stream[s].pending.remove(&next) else { break };

            // Stitch trailing (from previous block) with this block's prefix_bytes.
            let mut stitched = std::mem::take(&mut self.per_stream[s].trailing);
            stitched.extend_from_slice(&block.prefix_bytes);

            if !stitched.is_empty() {
                if has_n_newlines(&stitched, 4) {
                    let rec = parse_cross_block_record(&stitched)?;
                    self.per_stream[s].ready.push_back(rec);
                } else {
                    // Incomplete cross-block record. Stash and stop advancing this stream.
                    // This shouldn't happen normally (prefix_bytes was meant to close the
                    // trailing suffix). If it does, it means we need more prefix bytes from
                    // a later block, but we've already committed to advancing. In practice
                    // FastqBlockParse always includes the full prefix needed.
                    self.per_stream[s].trailing = stitched;
                    break;
                }
            }

            // Extend ready with this block's fully-parsed records.
            for rec in block.records {
                self.per_stream[s].ready.push_back(rec);
            }
            self.per_stream[s].trailing = block.suffix_bytes;
            self.per_stream[s].next_block = next + 1;
            if block.is_last {
                self.per_stream[s].eof = true;
                // If trailing is non-empty AND is_last, check whether it's a complete
                // record (normal — last record of file, no follow-up block to stitch
                // against) or truncated input (incomplete record).
                if !self.per_stream[s].trailing.is_empty() {
                    if !has_n_newlines(&self.per_stream[s].trailing, 4) {
                        anyhow::bail!(
                            "stream {s} reached EOF with incomplete trailing fragment ({} bytes)",
                            self.per_stream[s].trailing.len()
                        );
                    }
                    // The trailing has a complete record — parse it.
                    let rec = parse_cross_block_record(&self.per_stream[s].trailing)?;
                    self.per_stream[s].ready.push_back(rec);
                    self.per_stream[s].trailing.clear();
                }
                break; // no more blocks expected
            }
        }
        Ok(())
    }
}

fn has_n_newlines(data: &[u8], n: usize) -> bool {
    // Short-circuits as soon as the n-th newline is observed.
    memchr::memchr_iter(b'\n', data).take(n).count() >= n
}

fn parse_cross_block_record(bytes: &[u8]) -> Result<OwnedFastqRecord> {
    // Same shape as parse_one_record_owned from fastq_parse.rs but local to this module.
    let mut lines = bytes.split(|&b| b == b'\n');
    let name_line = lines.next().ok_or_else(|| anyhow::anyhow!("missing name line"))?;
    let seq_line = lines.next().ok_or_else(|| anyhow::anyhow!("missing sequence line"))?;
    let _plus = lines.next().ok_or_else(|| anyhow::anyhow!("missing + line"))?;
    let qual_line = lines.next().ok_or_else(|| anyhow::anyhow!("missing quality line"))?;
    let name = name_line
        .strip_prefix(b"@")
        .ok_or_else(|| anyhow::anyhow!("name line missing @ prefix: {name_line:?}"))?;
    Ok(OwnedFastqRecord {
        name: name.to_vec(),
        sequence: seq_line.to_vec(),
        quality: qual_line.to_vec(),
    })
}

#[cfg(test)]
mod merge_tests {
    use super::*;

    fn rec(n: &[u8], s: &[u8], q: &[u8]) -> OwnedFastqRecord {
        OwnedFastqRecord { name: n.to_vec(), sequence: s.to_vec(), quality: q.to_vec() }
    }

    fn block(
        block_idx: u64,
        stream_idx: usize,
        records: Vec<OwnedFastqRecord>,
        prefix: &[u8],
        suffix: &[u8],
        is_last: bool,
    ) -> BlockParsed {
        BlockParsed {
            block_idx,
            stream_idx,
            records,
            prefix_bytes: prefix.to_vec(),
            suffix_bytes: suffix.to_vec(),
            is_last,
        }
    }

    /// Run one `process` call and return the emitted batch (if any).
    /// `FastqBlockMerge` suppresses empty batches, so `None` means "waiting for more streams/blocks."
    fn process_emit(stage: &mut FastqBlockMerge, input: BlockParsed) -> Option<TemplateBatch> {
        let mut captured = None;
        stage
            .process(input, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        captured
    }

    #[test]
    fn test_block_merge_single_stream_no_fragments() {
        let mut stage = FastqBlockMerge::new(1);
        let b = block(0, 0, vec![rec(b"r1", b"A", b"I")], b"", b"", true);
        let out = process_emit(&mut stage, b).expect("single-stream must emit");
        assert_eq!(out.templates.len(), 1);
        assert_eq!(&out.templates[0].records[0].name, b"r1");
        assert_eq!(out.ordinal, 0);
    }

    #[test]
    fn test_block_merge_paired_waits_for_both_streams() {
        let mut stage = FastqBlockMerge::new(2);
        let out0 =
            process_emit(&mut stage, block(0, 0, vec![rec(b"r1", b"A", b"I")], b"", b"", true));
        assert!(out0.is_none(), "should wait for stream 1");
        let out1 =
            process_emit(&mut stage, block(0, 1, vec![rec(b"r1", b"B", b"J")], b"", b"", true))
                .expect("batch complete — must emit");
        assert_eq!(out1.templates.len(), 1);
        assert_eq!(out1.templates[0].records[0].sequence, b"A");
        assert_eq!(out1.templates[0].records[1].sequence, b"B");
        assert_eq!(out1.ordinal, 0);
    }

    #[test]
    fn test_block_merge_cross_block_record_stitching() {
        // Block 0: records=[r1]; suffix="@r2\nACGT\n+\n"  (record 2 starts, missing quality line)
        // Block 1: prefix="IIII\n"; records=[r3]; suffix=""
        // Expected: after processing both blocks, 3 records (r1, r2, r3).
        let mut stage = FastqBlockMerge::new(1);
        let b0 = block(0, 0, vec![rec(b"r1", b"A", b"I")], b"", b"@r2\nACGT\n+\n", false);
        let out0 = process_emit(&mut stage, b0).expect("block 0 must emit r1");
        assert_eq!(out0.templates.len(), 1); // only r1 emitted

        let b1 = block(1, 0, vec![rec(b"r3", b"G", b"K")], b"IIII\n", b"", true);
        let out1 = process_emit(&mut stage, b1).expect("block 1 must emit r2+r3");
        assert_eq!(out1.templates.len(), 2); // r2 (stitched) + r3
        assert_eq!(&out1.templates[0].records[0].name, b"r2");
        assert_eq!(&out1.templates[0].records[0].sequence, b"ACGT");
        assert_eq!(&out1.templates[0].records[0].quality, b"IIII");
        assert_eq!(&out1.templates[1].records[0].name, b"r3");
    }

    #[test]
    fn test_block_merge_out_of_order_blocks_wait() {
        // Block 1 arrives before block 0. Must wait for block 0 — no emit.
        let mut stage = FastqBlockMerge::new(1);
        let b1 = block(1, 0, vec![rec(b"r2", b"B", b"J")], b"", b"", true);
        let out1 = process_emit(&mut stage, b1);
        assert!(out1.is_none(), "should wait for block 0");

        let b0 = block(0, 0, vec![rec(b"r1", b"A", b"I")], b"", b"", false);
        let out0 = process_emit(&mut stage, b0).expect("block 0 arrival drains both blocks");
        assert_eq!(out0.templates.len(), 2);
        assert_eq!(&out0.templates[0].records[0].name, b"r1");
        assert_eq!(&out0.templates[1].records[0].name, b"r2");
    }

    #[test]
    fn test_block_merge_unequal_stream_counts_hold_surplus() {
        // Stream 0 has 3 records, stream 1 has 1. Only 1 template emitted;
        // stream 0's 2 surplus records are held for next zip.
        let mut stage = FastqBlockMerge::new(2);
        let b0_s0 = block(
            0,
            0,
            vec![rec(b"r1", b"A", b"I"), rec(b"r2", b"B", b"J"), rec(b"r3", b"C", b"K")],
            b"",
            b"",
            false,
        );
        assert!(process_emit(&mut stage, b0_s0).is_none(), "waiting for stream 1");

        let b0_s1 = block(0, 1, vec![rec(b"r1", b"D", b"L")], b"", b"", false);
        let out = process_emit(&mut stage, b0_s1).expect("batch complete");
        assert_eq!(out.templates.len(), 1);
        assert_eq!(out.templates[0].records[0].sequence, b"A");
        assert_eq!(out.templates[0].records[1].sequence, b"D");

        // Now more records for stream 1.
        let b1_s1 =
            block(1, 1, vec![rec(b"r2", b"E", b"M"), rec(b"r3", b"F", b"N")], b"", b"", true);
        let out2 = process_emit(&mut stage, b1_s1).expect("surplus+new must emit");
        assert_eq!(out2.templates.len(), 2); // r2 and r3
        assert_eq!(out2.templates[0].records[0].sequence, b"B"); // from stream 0's surplus
        assert_eq!(out2.templates[0].records[1].sequence, b"E");
        assert_eq!(out2.templates[1].records[0].sequence, b"C");
        assert_eq!(out2.templates[1].records[1].sequence, b"F");
    }

    #[test]
    fn test_block_merge_ordinal_increments_per_emit() {
        let mut stage = FastqBlockMerge::new(1);
        let out0 =
            process_emit(&mut stage, block(0, 0, vec![rec(b"r1", b"A", b"I")], b"", b"", false))
                .expect("first emit");
        assert_eq!(out0.templates.len(), 1);
        assert_eq!(out0.ordinal, 0);

        let out1 =
            process_emit(&mut stage, block(1, 0, vec![rec(b"r2", b"C", b"K")], b"", b"", true))
                .expect("second emit");
        assert_eq!(out1.templates.len(), 1);
        assert_eq!(out1.ordinal, 1);
    }
}
