//! Position-batch stage: groups sorted BAM records by template-
//! coordinate position key.
//!
//! [`SpecialStage`] — single-threaded because boundary detection needs
//! in-order per-record scanning. The generic engine drives it
//! alongside pool-parallel stages.
//!
//! ## Input
//!
//! [`SerializedBatch`] — concat-byte BAM records in template-coordinate
//! sort order, as emitted by
//! [`crate::runall::engine::stages::sort::SortStage`].
//!
//! ## Output
//!
//! [`PositionGroupBatch`] — one batch per distinct
//! `(primary, secondary, cb_hash)` sort key encountered; `data` holds
//! all records for that key in arrival order, with a fresh monotonic
//! `ordinal` counter (0..N) stamped on emit. `position_key` is
//! `(primary, secondary)` (cell-barcode hash is consumed internally).
//!
//! ## Ordering guarantees
//!
//! Strictly in-order by upstream sort key: the loop tracks the current
//! key and flushes whenever it changes, guaranteeing one output batch
//! per key and that batches are emitted in sort-key order. Requires
//! that the upstream `SortStage` emits records in canonical sort
//! order (it does).
//!
//! ## Memory model
//!
//! Bounded by the size of the largest position bucket. The in-flight
//! buffer is preallocated to 256 KB and reset after each flush.
//!
//! ## Determinism
//!
//! Fully deterministic given a deterministic (sorted) input stream.
//! Cell-barcode hashing (`cb_hasher()`) is a keyed but deterministic
//! hash.

use anyhow::Result;
use fgumi_raw_bam::RawRecordView;
use fgumi_raw_bam::flags::{SECONDARY, SUPPLEMENTARY};
use noodles::sam::Header;

use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::grouping_types::{PositionGroupBatch, iter_length_prefixed};
use crate::runall::engine::output_types::SerializedBatch;
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::special_stage::SpecialStage;
use crate::runall::engine::stage::SequencedItem;
use crate::sam::SamTag;
use crate::sort::{LibraryLookup, cb_hasher, extract_template_key_inline};

/// Boundary-detecting position-batch stage.
pub struct PositionBatchStage {
    header: Header,
    cell_tag: Option<SamTag>,
    /// When `false` (the default for the consensus path), secondary (`0x100`)
    /// and supplementary (`0x800`) alignments are dropped before they enter
    /// the position grouping logic — matching standalone
    /// `RecordPositionGrouper`'s default (`grouper.rs::process_record`).
    /// Set `true` for any future dedup path, which needs them in templates.
    include_secondary_supplementary: bool,
}

impl PositionBatchStage {
    /// Construct a new `PositionBatchStage` with no cell-barcode hashing.
    /// Secondary/supplementary alignments are skipped (consensus-path default).
    #[must_use]
    pub fn new(header: Header) -> Self {
        Self { header, cell_tag: None, include_secondary_supplementary: false }
    }

    /// Construct with cell-barcode hashing for single-cell workflows.
    /// Secondary/supplementary alignments are skipped (consensus-path default).
    #[must_use]
    pub fn with_cell_tag(header: Header, cell_tag: SamTag) -> Self {
        Self { header, cell_tag: Some(cell_tag), include_secondary_supplementary: false }
    }

    /// Variant that retains secondary/supplementary alignments in the output
    /// position groups. Mirrors
    /// `RecordPositionGrouper::with_secondary_supplementary` for any future
    /// dedup path that needs them in templates. The consensus path uses the
    /// default-skip constructors above.
    #[must_use]
    pub fn with_secondary_supplementary(self) -> Self {
        Self { include_secondary_supplementary: true, ..self }
    }
}

impl SpecialStage for PositionBatchStage {
    type Input = SerializedBatch;
    type Output = PositionGroupBatch;

    #[tracing::instrument(name = "position_batch", skip_all)]
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<SerializedBatch>>,
        output: Box<dyn OutputQueue<PositionGroupBatch>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let lib_lookup = LibraryLookup::from_header(&self.header);
        let hasher = cb_hasher();

        let mut current_key: Option<(u64, u64, u64)> = None;
        let mut current_data: Vec<u8> = Vec::with_capacity(256 * 1024);
        let mut current_count: u64 = 0;
        let mut out_ordinal: u64 = 0;

        let mut backoff = crate::runall::engine::backoff::Backoff::new();

        let flush = |current_key: &mut Option<(u64, u64, u64)>,
                     current_data: &mut Vec<u8>,
                     current_count: &mut u64,
                     out_ordinal: &mut u64|
         -> Result<()> {
            if *current_count == 0 {
                return Ok(());
            }
            let Some((primary, secondary, _cb_hash)) = *current_key else {
                return Ok(());
            };
            let data = std::mem::take(current_data);
            let batch = PositionGroupBatch {
                data,
                record_count: *current_count,
                ordinal: *out_ordinal,
                position_key: (primary, secondary),
            };
            let mem = batch.data.len();
            let item = SequencedItem::new(*out_ordinal, batch, mem);
            // `push_until_cancelled` only returns Err when the shared cancel
            // token has fired. The token is set by whatever stage actually
            // failed (a panicking pool worker, an upstream/downstream error,
            // or the watchdog); surfacing our own "cancelled during push"
            // error here would shadow that real cause when the driver picks
            // the first error in join order. Treat cancellation as a clean
            // early exit and let the originating error propagate.
            if output.push_until_cancelled(item, &cancel).is_err() {
                return Ok(());
            }
            *out_ordinal += 1;
            *current_count = 0;
            *current_data = Vec::with_capacity(256 * 1024);
            *current_key = None;
            Ok(())
        };

        loop {
            if cancel.is_cancelled() {
                break;
            }
            if let Some(item) = input.pop() {
                // Iterate length-prefixed records within the input batch.
                // Upstream `SortStage` emits `SerializedBatch` items with
                // monotonic ordinals in sort-key order, so records within
                // and across batches arrive in sort-key order already.
                for rec_result in iter_length_prefixed(&item.item.primary.data) {
                    let rec = rec_result?;

                    // Drop secondary (0x100) and supplementary (0x800)
                    // alignments before computing position keys, matching
                    // standalone `RecordPositionGrouper`'s default. Without
                    // this, supplementary alignments that share a position
                    // key with a primary (common on chimeric/SA reads in
                    // real WGS/WES data) flow into `GroupAssignStage` and
                    // perturb consensus calls — diverging the runall
                    // pipeline from standalone-equivalent output.
                    if !self.include_secondary_supplementary {
                        let flg = RawRecordView::new(rec).flags();
                        if flg & (SECONDARY | SUPPLEMENTARY) != 0 {
                            continue;
                        }
                    }

                    let tkey =
                        extract_template_key_inline(rec, &lib_lookup, self.cell_tag, &hasher);
                    let pkey = (tkey.primary, tkey.secondary, tkey.cb_hash);

                    match current_key {
                        Some(ck) if ck == pkey => {
                            let bs = u32::try_from(rec.len()).map_err(|_| {
                                anyhow::anyhow!("PositionBatchStage: record exceeds u32::MAX bytes")
                            })?;
                            current_data.extend_from_slice(&bs.to_le_bytes());
                            current_data.extend_from_slice(rec);
                            current_count += 1;
                        }
                        _ => {
                            flush(
                                &mut current_key,
                                &mut current_data,
                                &mut current_count,
                                &mut out_ordinal,
                            )?;
                            current_key = Some(pkey);
                            let bs = u32::try_from(rec.len()).map_err(|_| {
                                anyhow::anyhow!("PositionBatchStage: record exceeds u32::MAX bytes")
                            })?;
                            current_data.extend_from_slice(&bs.to_le_bytes());
                            current_data.extend_from_slice(rec);
                            current_count += 1;
                        }
                    }
                }
                backoff.reset();
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Flush the trailing batch.
        flush(&mut current_key, &mut current_data, &mut current_count, &mut out_ordinal)?;
        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "PositionBatchStage"
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::runall::engine::grouping_types::{iter_length_prefixed, serialize_records};
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::output_types::RawBytes;
    use crate::runall::engine::queue::StageQueue;
    use crate::runall::engine::stage::SequencedItem;
    use fgumi_raw_bam::RawRecordView;

    fn empty_header() -> Header {
        use noodles::sam::header::record::value::Map;
        Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .build()
    }

    /// Build a minimal raw BAM record at `ref_id=0`, `pos=100`,
    /// `mate_pos=200`, `4M` cigar, with the given name, flag, and `RX` UMI
    /// tag. Mirrors the helper in `group_assign.rs::tests`.
    #[allow(clippy::cast_possible_truncation)]
    fn make_raw_record(name: &[u8], flag: u16, umi: &[u8]) -> Vec<u8> {
        let seq_len: usize = 4;
        let l_read_name = (name.len() + 1) as u8;
        let cigar_ops: &[u32] = &[(seq_len as u32) << 4]; // 4M
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_bytes = seq_len.div_ceil(2);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len;
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id
        buf[4..8].copy_from_slice(&100i32.to_le_bytes()); // pos
        buf[8] = l_read_name;
        buf[9] = 60; // mapq
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        buf[20..24].copy_from_slice(&0i32.to_le_bytes()); // mate ref_id
        buf[24..28].copy_from_slice(&200i32.to_le_bytes()); // mate pos
        buf[28..32].copy_from_slice(&150i32.to_le_bytes()); // tlen

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        if !umi.is_empty() {
            buf.extend_from_slice(b"RXZ");
            buf.extend_from_slice(umi);
            buf.push(0);
        }

        buf
    }

    /// Drive a `PositionBatchStage` with a single input batch and collect all
    /// output batches.
    fn run_stage_with(
        stage: Box<PositionBatchStage>,
        records: &[Vec<u8>],
    ) -> Vec<PositionGroupBatch> {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("pb_in", 16, 100_000, tracker.clone()));
        let output: Arc<StageQueue<PositionGroupBatch>> =
            Arc::new(StageQueue::new("pb_out", 16, 100_000, tracker));

        let (data, record_count) = serialize_records(records).unwrap();
        let batch = SerializedBatch {
            primary: RawBytes { data, record_count },
            secondary: None,
            ordinal: 0,
        };
        input.push(SequencedItem::new(0, batch, 0)).unwrap();
        input.close();

        <PositionBatchStage as SpecialStage>::run(
            stage,
            Box::new(input),
            Box::new(output.clone()),
            CancelToken::new(),
        )
        .unwrap();

        let mut out = Vec::new();
        while let Some(item) = output.pop() {
            out.push(item.item);
        }
        out
    }

    #[test]
    fn test_position_batch_empty_input() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("pb_in", 16, 100_000, tracker.clone()));
        let output: Arc<StageQueue<PositionGroupBatch>> =
            Arc::new(StageQueue::new("pb_out", 16, 100_000, tracker));

        input.close();

        let stage = Box::new(PositionBatchStage::new(empty_header()));
        <PositionBatchStage as SpecialStage>::run(
            stage,
            Box::new(input),
            Box::new(output.clone()),
            CancelToken::new(),
        )
        .unwrap();

        assert!(OutputQueue::is_closed(&output));
        assert!(output.pop().is_none());
    }

    /// Regression: secondary (0x100) and supplementary (0x800) alignments
    /// must be dropped by `PositionBatchStage` in its default mode, matching
    /// the behavior of standalone `fgumi group` (which uses
    /// `RecordPositionGrouper` with `include_secondary_supplementary=false`,
    /// see `src/lib/grouper.rs::process_record`).
    ///
    /// Pre-fix: the inner record loop in `run()` extended `current_data` with
    /// every record from `iter_length_prefixed(...)`, with no flag check, so
    /// supplementaries that landed in the same template-coordinate position
    /// group as a primary (common on chimeric/SA reads in real WGS/WES data)
    /// flowed downstream into `GroupAssignStage` and perturbed consensus
    /// calls — diverging the runall pipeline from standalone-equivalent
    /// output by 12 missing-pass records on SRR6109255.
    #[test]
    fn test_position_batch_default_drops_secondary_and_supplementary() {
        const PAIRED_R1: u16 = 0x1 | 0x40;
        const PAIRED_R2: u16 = 0x1 | 0x80;
        const SECONDARY: u16 = 0x100;
        const SUPPLEMENTARY: u16 = 0x800;

        // Three records sharing one qname + one position + one UMI, so they
        // would land in a single output PositionGroupBatch if the stage
        // didn't filter on flags. The two primaries must survive; the
        // supplementary (and a separate secondary) must be dropped.
        let r1_primary = make_raw_record(b"read1", PAIRED_R1, b"ACGT");
        let r2_primary = make_raw_record(b"read1", PAIRED_R2, b"ACGT");
        let r1_supplementary = make_raw_record(b"read1", PAIRED_R1 | SUPPLEMENTARY, b"ACGT");
        let r1_secondary = make_raw_record(b"read1", PAIRED_R1 | SECONDARY, b"ACGT");

        let stage = Box::new(PositionBatchStage::new(empty_header()));
        let batches =
            run_stage_with(stage, &[r1_primary, r2_primary, r1_supplementary, r1_secondary]);

        let total_records: u64 = batches.iter().map(|b| b.record_count).sum();
        assert_eq!(
            total_records,
            2,
            "expected only the 2 primary records to survive PositionBatchStage; \
             secondary/supplementary alignments must be filtered before \
             group_assign so consensus calls match standalone `fgumi group`. \
             Got {total_records} records across {} batch(es).",
            batches.len()
        );

        for batch in &batches {
            for rec_result in iter_length_prefixed(&batch.data) {
                let rec = rec_result.unwrap();
                let flags = RawRecordView::new(rec).flags();
                assert_eq!(
                    flags & (SECONDARY | SUPPLEMENTARY),
                    0,
                    "PositionBatchStage emitted a non-primary record (flags={flags:#06x}) \
                     in default mode; it must filter SECONDARY/SUPPLEMENTARY before emit"
                );
            }
        }
    }
}
