//! Positional (key-lockstep) comparison engine for `fgumi compare bams`.
//!
//! This engine pairs record `i` of `bam1` with record `i` of `bam2` — the two
//! streams are read in strict lockstep by index, exactly as they appear on
//! disk. For each paired index it first checks that the two records carry the
//! same [`RecordKey`](super::super::record_key::RecordKey) *before* ever
//! comparing content. This ordering is the sound core the module exists to
//! protect:
//!
//! - If the keys disagree, the streams have desynchronized (e.g. one file
//!   dropped, added, or reordered a record). [`positional_compare`] records
//!   this as a key mismatch and **stops pairing** — it never attempts to
//!   resync by searching ahead, and it never falls back to comparing the
//!   mismatched records' content. A content match on a key mismatch would be
//!   coincidence, not evidence of parity, so it must never be allowed to mask
//!   the desync.
//! - If the keys agree, the pair is compared under a
//!   [`ContentPredicate`](super::content::ContentPredicate) via
//!   [`content_diffs`](super::content::content_diffs).
//!
//! A record-count mismatch (one file longer than the other) is reported as a
//! presence difference, independent of and in addition to any key mismatch.

use std::path::Path;

use anyhow::{Result, bail};

use fgumi_raw_bam::RawRecord;

use super::super::bams::{RawBatchMessage, start_raw_batch_reader};
use super::super::record_key::record_keys_match;
use super::content::{ContentPredicate, content_diffs};
use super::header::{compare_headers, fold_header_diffs};

/// Outcome of a [`positional_compare`] run.
#[derive(Debug, Default, Clone)]
pub struct PositionalOutcome {
    /// Total number of records read from `bam1`.
    pub bam1_count: u64,
    /// Total number of records read from `bam2`.
    pub bam2_count: u64,
    /// Number of index-paired records whose content differed under the
    /// configured [`ContentPredicate`], among records paired before any key
    /// mismatch was found.
    pub content_diffs: u64,
    /// The index of the first pair whose `RecordKey`s disagreed, if any.
    ///
    /// Once set, no record at or past this index was compared for content:
    /// pairing stopped at the first desync rather than attempting to resync.
    pub key_mismatch_at: Option<u64>,
    /// `true` if the two inputs' `@HD`/`@SQ`/`@RG` headers disagreed on a field
    /// [`compare_headers`](super::header::compare_headers) considers significant (`@PG`/`@CO`
    /// are normalized and never contribute here).
    pub header_mismatch: bool,
    /// Human-readable diff strings (header, key mismatch, content diffs, and/or the
    /// record-count mismatch note), capped at the caller-supplied `max_diffs`.
    pub diff_details: Vec<String>,
}

impl PositionalOutcome {
    /// Returns `true` iff the two streams matched: equal record counts, no
    /// `RecordKey` mismatch, no content diffs among the records paired before any such
    /// mismatch, and no significant header divergence.
    #[must_use]
    pub fn is_match(&self) -> bool {
        self.bam1_count == self.bam2_count
            && self.content_diffs == 0
            && self.key_mismatch_at.is_none()
            && !self.header_mismatch
    }
}

/// Compares two BAM files in strict positional lockstep.
///
/// Record `i` of `bam1` is paired with record `i` of `bam2` purely by index —
/// this engine never reorders, searches ahead, or otherwise resynchronizes.
/// For each pair, `RecordKey` equality
/// is checked first; a mismatch stops pairing immediately (see the module
/// docs), while a match is handed to `content_diffs` under `pred`. Both
/// readers are drained to completion regardless, so `bam1_count`/`bam2_count`
/// always reflect the true totals even after pairing has stopped.
///
/// Also compares the two inputs' headers via
/// [`compare_headers`](super::header::compare_headers) (`@HD`/`@SQ`/`@RG`, normalizing away
/// `@PG`/`@CO`); a significant divergence is folded into
/// `PositionalOutcome::header_mismatch`/`PositionalOutcome::is_match` alongside the
/// record-level findings.
///
/// `threads` and `batch_size` control the underlying double-buffered readers
/// (see `start_raw_batch_reader`); `max_diffs` caps the number of entries
/// collected in `PositionalOutcome::diff_details`.
///
/// # Errors
///
/// Returns an error if either BAM file cannot be opened, or if a read error
/// occurs while streaming either file.
pub fn positional_compare(
    bam1: &Path,
    bam2: &Path,
    threads: usize,
    batch_size: usize,
    max_diffs: usize,
    pred: ContentPredicate,
) -> Result<PositionalOutcome> {
    let mut outcome = PositionalOutcome::default();

    let (rx1, header1) = start_raw_batch_reader(bam1.to_path_buf(), threads, batch_size)?;
    let (rx2, header2) = start_raw_batch_reader(bam2.to_path_buf(), threads, batch_size)?;

    fold_header_diffs(
        compare_headers(&header1, &header2),
        &mut outcome.header_mismatch,
        &mut outcome.diff_details,
        max_diffs,
    );

    let mut bam1_eof = false;
    let mut bam2_eof = false;
    let mut pending_batch1: Option<Vec<RawRecord>> = None;
    let mut pending_batch2: Option<Vec<RawRecord>> = None;
    let mut current_index = 0u64;

    loop {
        if pending_batch1.is_none() && !bam1_eof {
            match rx1.recv() {
                Ok(RawBatchMessage::Batch(batch)) => {
                    outcome.bam1_count += batch.len() as u64;
                    pending_batch1 = Some(batch);
                }
                Ok(RawBatchMessage::Eof) => bam1_eof = true,
                Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                // A clean end-of-stream is always signalled by an explicit `Eof`
                // message (see `start_raw_batch_reader`); a bare channel
                // disconnect therefore means the reader thread died (e.g.
                // panicked) or the stream was truncated without an `Eof`. Treat
                // it as fatal rather than silently accepting a short read as a
                // clean comparison.
                Err(_) => bail!("BAM1 reader disconnected before EOF"),
            }
        }

        if pending_batch2.is_none() && !bam2_eof {
            match rx2.recv() {
                Ok(RawBatchMessage::Batch(batch)) => {
                    outcome.bam2_count += batch.len() as u64;
                    pending_batch2 = Some(batch);
                }
                Ok(RawBatchMessage::Eof) => bam2_eof = true,
                Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
                // See the BAM1 branch above: a disconnect without a prior `Eof`
                // is a dead/truncated reader, not a clean end-of-stream.
                Err(_) => bail!("BAM2 reader disconnected before EOF"),
            }
        }

        match (&pending_batch1, &pending_batch2) {
            (None, None) => break,
            (Some(_), None) | (None, Some(_)) => {
                // One file exhausted before the other: a presence DIFFER,
                // independent of (and in addition to) any key mismatch.
                if outcome.diff_details.len() < max_diffs {
                    outcome
                        .diff_details
                        .push("BAM files have different number of records".to_string());
                }
                // Drain the remaining side to get an accurate final count. A
                // disconnect *before* the explicit `Eof` is a dead/truncated
                // reader (same reasoning as the receive branches above), so bail
                // rather than accept a short count as complete.
                if pending_batch1.is_some() {
                    loop {
                        match rx1.recv() {
                            Ok(RawBatchMessage::Batch(batch)) => {
                                outcome.bam1_count += batch.len() as u64;
                            }
                            Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM1: {e}"),
                            Ok(RawBatchMessage::Eof) => break,
                            Err(_) => bail!("BAM1 reader disconnected before EOF"),
                        }
                    }
                }
                if pending_batch2.is_some() {
                    loop {
                        match rx2.recv() {
                            Ok(RawBatchMessage::Batch(batch)) => {
                                outcome.bam2_count += batch.len() as u64;
                            }
                            Ok(RawBatchMessage::Error(e)) => bail!("Error reading BAM2: {e}"),
                            Ok(RawBatchMessage::Eof) => break,
                            Err(_) => bail!("BAM2 reader disconnected before EOF"),
                        }
                    }
                }
                break;
            }
            (Some(_), Some(_)) => {}
        }

        let batch1 = pending_batch1.take().expect("guarded by (Some, Some) match above");
        let batch2 = pending_batch2.take().expect("guarded by (Some, Some) match above");

        let min_len = batch1.len().min(batch2.len());
        let (cmp_batch1, remainder1) = batch1.split_at(min_len);
        let (cmp_batch2, remainder2) = batch2.split_at(min_len);

        // Once a RecordKey mismatch is found, pairing stops permanently: no
        // further record-key/content comparisons are performed, even though the
        // loop keeps pulling batches (to drain both readers for accurate counts).
        if outcome.key_mismatch_at.is_none() {
            for (offset, (a, b)) in cmp_batch1.iter().zip(cmp_batch2.iter()).enumerate() {
                let index = current_index + offset as u64;

                if !record_keys_match(a, b) {
                    outcome.key_mismatch_at = Some(index);
                    if outcome.diff_details.len() < max_diffs {
                        outcome.diff_details.push(format!(
                            "record {index}: RecordKey mismatch — pairing stopped (no resync)"
                        ));
                    }
                    break;
                }

                if let Some(diffs) = content_diffs(a.as_ref(), b.as_ref(), pred, &header1, &header2)
                {
                    outcome.content_diffs += 1;
                    if outcome.diff_details.len() < max_diffs {
                        outcome
                            .diff_details
                            .push(format!("record {index}: {joined}", joined = diffs.join("; ")));
                    }
                }
            }
        }

        current_index += min_len as u64;

        if !remainder1.is_empty() {
            pending_batch1 = Some(remainder1.to_vec());
        }
        if !remainder2.is_empty() {
            pending_batch2 = Some(remainder2.to_vec());
        }
    }

    Ok(outcome)
}
