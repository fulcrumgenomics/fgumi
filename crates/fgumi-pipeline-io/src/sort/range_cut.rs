//! Batch cutter for the range-parallel final merge.
//!
//! Merged key ranges come back from the pool in range order. The cutter turns
//! that stream into exactly the batch boundaries the serial merge produces (see
//! `SortMerge::next_batch` and `plan_fast_path_blocks`): per record
//! `count += 1; bytes += len + frame_overhead`, and the batch closes after that
//! record by the serial merge's own rule, [`batch_full`]. The
//! running count, byte total and open segments carry across range boundaries.
//!
//! The cutter emits descriptors only: byte spans into each range's `framed`
//! buffer. It reads `lens` and never touches `framed`; assembling the bytes is
//! done later by an assembly job.
//!
//! A range the merge's own thread merges once every earlier range is cut
//! (a fused range) does not go through `framed`: its records are pushed into
//! the open batch's builder as they are merged ([`BatchCutter::push_fused`]),
//! the serial merge's own loop. The cutter then only counts them; the
//! builder (the *lead* of the open batch, held by the caller) carries over
//! into the next batch plan that closes ([`BatchPlan::lead`]).

use std::collections::VecDeque;
use std::sync::Arc;

use fgumi_sort::{RECORD_LEN_BYTES, RangeOut};

use crate::sort::merge::batch_full;

/// A record-aligned byte span of one range's `framed` buffer.
#[derive(Clone)]
pub struct Segment {
    /// The range the span belongs to.
    pub range: Arc<RangeOut>,
    /// The range's position in range order.
    pub range_idx: usize,
    /// Start of the span in `range.framed` (at a record's 4-byte prefix).
    pub start: usize,
    /// End of the span in `range.framed` (just past a record's body).
    pub end: usize,
    /// Number of records in the span.
    pub records: u32,
}

/// One output batch: its ordinal and the segments that make it up, in order.
pub struct BatchPlan {
    /// Dense batch ordinal, from 0.
    pub ordinal: u64,
    /// The batch starts with the records of the caller's lead builder (a
    /// fused range's open batch), and `segments` follow them.
    pub lead: bool,
    /// The spans, in merge order, that make up the batch.
    pub segments: Vec<Segment>,
    /// Records in the batch.
    pub records: u64,
    /// Bytes of `[u32 len][body]` framing in the batch.
    pub framed_bytes: u64,
}

/// A merged range whose record `i` has body `body(i, lens[i])` (tests).
#[cfg(test)]
pub(crate) fn range_out_for_test(
    lens: &[u32],
    body: impl Fn(usize, u32) -> Vec<u8>,
) -> Arc<RangeOut> {
    let mut framed = Vec::new();
    for (i, &len) in lens.iter().enumerate() {
        framed.extend_from_slice(&len.to_le_bytes());
        framed.extend_from_slice(&body(i, len));
    }
    Arc::new(RangeOut { framed, lens: lens.to_vec(), ..RangeOut::default() })
}

/// Cuts a stream of merged ranges into the serial merge's batches.
pub(crate) struct BatchCutter {
    target_batch_count: usize,
    output_byte_limit: u64,
    frame_overhead: usize,
    next_ordinal: u64,
    /// Records in the open batch.
    count: usize,
    /// Serial-rule byte total of the open batch (`len + frame_overhead` per record).
    bytes: u64,
    /// Actual framed bytes of the open batch.
    framed_bytes: u64,
    /// Segments of the open batch.
    open: Vec<Segment>,
    /// The open batch starts with records held in the caller's lead builder.
    lead_open: bool,
}

impl BatchCutter {
    /// A cutter applying the serial rule with the given caps; `frame_overhead` is
    /// the output builder's `FRAME_OVERHEAD_PER_RECORD`.
    pub fn new(target_batch_count: usize, output_byte_limit: u64, frame_overhead: usize) -> Self {
        Self {
            target_batch_count,
            output_byte_limit,
            frame_overhead,
            next_ordinal: 0,
            count: 0,
            bytes: 0,
            framed_bytes: 0,
            open: Vec::new(),
            lead_open: false,
        }
    }

    /// Feed the next range (ranges must arrive in range order); closed batches are
    /// appended to `out`. Each segment takes its own handle on `range`.
    pub fn cut(&mut self, range_idx: usize, range: &Arc<RangeOut>, out: &mut VecDeque<BatchPlan>) {
        let mut seg_start = 0usize;
        let mut seg_records = 0u32;
        let mut pos = 0usize;
        for &len in &range.lens {
            pos += RECORD_LEN_BYTES + len as usize;
            seg_records += 1;
            self.count += 1;
            self.bytes = self.bytes.saturating_add(u64::from(len) + self.frame_overhead as u64);
            if batch_full(self.count, self.bytes, self.target_batch_count, self.output_byte_limit) {
                self.open.push(Segment {
                    range: Arc::clone(range),
                    range_idx,
                    start: seg_start,
                    end: pos,
                    records: seg_records,
                });
                self.framed_bytes += (pos - seg_start) as u64;
                out.push_back(self.close());
                seg_start = pos;
                seg_records = 0;
            }
        }
        if seg_records > 0 {
            self.open.push(Segment {
                range: Arc::clone(range),
                range_idx,
                start: seg_start,
                end: pos,
                records: seg_records,
            });
            self.framed_bytes += (pos - seg_start) as u64;
        }
    }

    /// The open batch's segments, taken so the caller copies them into its
    /// lead builder before a fused range ([`Self::push_fused`]); the cutter
    /// keeps counting them as the open batch's.
    pub fn take_open_segments(&mut self) -> Vec<Segment> {
        if self.count > 0 {
            self.lead_open = true;
        }
        std::mem::take(&mut self.open)
    }

    /// Count a fused record of body length `len`, which the caller pushed
    /// into its lead builder. Returns the ordinal of the batch it closed,
    /// when it closed one (the caller's lead builder is then that batch).
    ///
    /// # Panics
    ///
    /// If open segments were not taken first ([`Self::take_open_segments`]).
    pub fn push_fused(&mut self, len: u32) -> Option<u64> {
        assert!(self.open.is_empty(), "a fused record follows the open segments' copy");
        self.count += 1;
        self.bytes = self.bytes.saturating_add(u64::from(len) + self.frame_overhead as u64);
        self.framed_bytes += (RECORD_LEN_BYTES + len as usize) as u64;
        self.lead_open = true;
        if !batch_full(self.count, self.bytes, self.target_batch_count, self.output_byte_limit) {
            return None;
        }
        let ordinal = self.next_ordinal;
        drop(self.close());
        Some(ordinal)
    }

    /// The ordinal of the open (next) batch.
    pub fn next_ordinal(&self) -> u64 {
        self.next_ordinal
    }

    /// Emit the non-empty partial batch left after the last range.
    pub fn finish(&mut self) -> Option<BatchPlan> {
        if self.count == 0 { None } else { Some(self.close()) }
    }

    /// The range index of the oldest segment in the open batch, if any.
    pub fn oldest_open_range(&self) -> Option<usize> {
        self.open.first().map(|s| s.range_idx)
    }

    fn close(&mut self) -> BatchPlan {
        let plan = BatchPlan {
            ordinal: self.next_ordinal,
            lead: std::mem::take(&mut self.lead_open),
            segments: std::mem::take(&mut self.open),
            records: self.count as u64,
            framed_bytes: self.framed_bytes,
        };
        self.next_ordinal += 1;
        self.count = 0;
        self.bytes = 0;
        self.framed_bytes = 0;
        plan
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sort::merge::plan_fast_path_blocks;
    use crate::sort::merge::{BlockOutput, RecordBatchOutput};
    use crate::sort::protocol::MemoryChunkErased;
    use fgumi_sort::{InMemoryChunk, RawCoordinateKey};
    use proptest::prelude::*;
    use rstest::rstest;

    const LIMIT: u64 = 4 << 20;

    fn make_range(lens: &[u32]) -> Arc<RangeOut> {
        range_out_for_test(lens, |i, len| vec![i.to_le_bytes()[0]; len as usize])
    }

    /// The serial plan over one chunk holding every record, as global record spans.
    fn serial_spans(
        ranges: &[Vec<u32>],
        target: usize,
        limit: u64,
        overhead: usize,
    ) -> Vec<(usize, usize)> {
        let records: Vec<_> = ranges
            .iter()
            .flatten()
            .enumerate()
            .map(|(i, &len)| (RawCoordinateKey { sort_key: i as u64 }, vec![0u8; len as usize]))
            .collect();
        let chunk = MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(records));
        let blocks = match overhead {
            0 => plan_fast_path_blocks::<RecordBatchOutput>(&chunk, target, limit),
            _ => plan_fast_path_blocks::<BlockOutput>(&chunk, target, limit),
        };
        blocks.iter().map(|b| (b.start, b.end)).collect()
    }

    fn cutter_batches(
        ranges: &[Vec<u32>],
        target: usize,
        limit: u64,
        overhead: usize,
    ) -> Vec<BatchPlan> {
        let mut cutter = BatchCutter::new(target, limit, overhead);
        let mut out = VecDeque::new();
        for (idx, lens) in ranges.iter().enumerate() {
            cutter.cut(idx, &make_range(lens), &mut out);
        }
        out.extend(cutter.finish());
        out.into_iter().collect()
    }

    fn cutter_spans(plans: &[BatchPlan]) -> Vec<(usize, usize)> {
        let mut start = 0usize;
        plans
            .iter()
            .map(|p| {
                let span = (start, start + usize::try_from(p.records).unwrap());
                start = span.1;
                span
            })
            .collect()
    }

    fn check(ranges: &[Vec<u32>], target: usize, limit: u64, overhead: usize) {
        let plans = cutter_batches(ranges, target, limit, overhead);
        assert_eq!(cutter_spans(&plans), serial_spans(ranges, target, limit, overhead));
        for (i, p) in plans.iter().enumerate() {
            assert_eq!(p.ordinal, i as u64);
            assert_eq!(p.records, p.segments.iter().map(|s| u64::from(s.records)).sum::<u64>());
            assert_eq!(
                p.framed_bytes,
                p.segments.iter().map(|s| (s.end - s.start) as u64).sum::<u64>()
            );
        }
    }

    #[rstest]
    #[case::single_range(vec![vec![100; 5000]], 1024, 4 << 20, 4)]
    #[case::carry_across_ranges(vec![vec![100; 700], vec![100; 700]], 1024, 4 << 20, 4)]
    #[case::carry_across_empty_range(vec![vec![100; 700], vec![], vec![100; 700]], 1024, 4 << 20, 4)]
    #[case::range_shorter_than_batch(vec![vec![50; 3]; 1000], 1024, 4 << 20, 4)]
    #[case::byte_cap(vec![vec![5000; 2000]], 1024, 4 << 20, 4)]
    #[case::record_at_least_cap(vec![vec![10, 4 << 20, 10]], 1024, 4 << 20, 4)]
    #[case::no_overhead_record_batch(vec![vec![100; 3000]], 1024, 4 << 20, 0)]
    #[case::all_empty(vec![vec![], vec![]], 1024, 4 << 20, 4)]
    fn cutter_matches_the_serial_rule(
        #[case] ranges: Vec<Vec<u32>>,
        #[case] target: usize,
        #[case] limit: u64,
        #[case] overhead: usize,
    ) {
        check(&ranges, target, limit, overhead);
    }

    #[test]
    fn oldest_open_range_tracks_the_carry() {
        let mut cutter = BatchCutter::new(10, LIMIT, 4);
        let mut out = VecDeque::new();
        assert_eq!(cutter.oldest_open_range(), None);
        cutter.cut(0, &make_range(&[1; 4]), &mut out);
        cutter.cut(1, &make_range(&[1; 4]), &mut out);
        assert_eq!(cutter.oldest_open_range(), Some(0));
        cutter.cut(2, &make_range(&[1; 4]), &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(cutter.oldest_open_range(), Some(2));
        assert!(cutter.finish().is_some());
        assert_eq!(cutter.oldest_open_range(), None);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn cutter_matches_serial_on_random_ranges(
            ranges in prop::collection::vec(prop::collection::vec(1u32..9000, 0..400), 0..8),
            target in 1usize..600,
            overhead in prop::sample::select(vec![0usize, 4]),
        ) {
            check(&ranges, target, 600_000, overhead);
        }

        #[test]
        fn segments_are_record_aligned_and_cover(
            ranges in prop::collection::vec(prop::collection::vec(1u32..9000, 0..300), 0..6),
            target in 1usize..400,
        ) {
            let plans = cutter_batches(&ranges, target, 600_000, 4);
            let mut got = Vec::new();
            for p in &plans {
                for s in &p.segments {
                    let span = &s.range.framed[s.start..s.end];
                    // Walk the span: it must be exactly `records` whole records.
                    let (mut at, mut n) = (0usize, 0u32);
                    while at < span.len() {
                        let len = u32::from_le_bytes(span[at..at + 4].try_into().unwrap());
                        at += 4 + len as usize;
                        n += 1;
                    }
                    prop_assert_eq!(at, span.len());
                    prop_assert_eq!(n, s.records);
                    got.extend_from_slice(span);
                }
            }
            let want: Vec<u8> = ranges
                .iter()
                .flat_map(|r| make_range(r).framed.clone())
                .collect();
            prop_assert_eq!(got, want);
        }
    }
}
