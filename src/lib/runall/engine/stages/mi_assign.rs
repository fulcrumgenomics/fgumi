//! `MiAssignStage`: serial-ordered `MoleculeId` offset application.
//!
//! Single-threaded `SpecialStage`. Buffers `MiGroupBatch` arrivals in a
//! `ReorderBuffer` keyed by `ordinal`, pops batches in monotonic pipeline-
//! serial order, applies a cumulative offset via
//! [`fgumi_umi::MoleculeId::with_offset`], and injects the final MI tag
//! bytes into each record in `MiGroup.data`. Output `MiGroupBatch`es carry
//! globally-unique MI integers and MI-tagged record bytes; their `local_mi`
//! fields are cleared.
//!
//! This is runall's analog of the `mi_assign_fn` zone PR #319 added to the
//! unified BAM pipeline for standalone `group`/`dedup`. The unified pipeline
//! has its own scaffolding; runall has its own engine and gets its own copy
//! of the same idea.
//!
//! ## Why a `SpecialStage`
//!
//! The `MoleculeId` cumulative-offset advance is inherently serial in pipeline
//! order. The pool model (`Stage`) doesn't fit because we need a barrier
//! between the parallel `GroupAssignStage` workers and the per-batch offset
//! advance. A single-threaded `SpecialStage` with an internal `ReorderBuffer`
//! is the same shape as `Coalesce`.

use anyhow::Result;
use fgumi_raw_bam::update_string_tag;

use crate::runall::engine::backoff::Backoff;
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::grouping_types::{MiGroupBatch, iter_length_prefixed};
use crate::runall::engine::reorder::ReorderBuffer;
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::special_stage::SpecialStage;
use crate::runall::engine::stage::SequencedItem;

/// Single-threaded stage that converts the local `MoleculeId`s emitted by
/// `GroupAssignStage` into globally-unique IDs in serial pipeline order and
/// writes the final MI tag bytes into every record.
pub struct MiAssignStage;

impl MiAssignStage {
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl Default for MiAssignStage {
    fn default() -> Self {
        Self::new()
    }
}

impl SpecialStage for MiAssignStage {
    type Input = MiGroupBatch;
    type Output = MiGroupBatch;

    #[tracing::instrument(name = "mi_assign", skip_all)]
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        output: Box<dyn OutputQueue<Self::Output>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let mut reorder: ReorderBuffer<MiGroupBatch> = ReorderBuffer::new();
        let mut next_mi_base: u64 = 0;
        let mut backoff = Backoff::new();

        loop {
            if cancel.is_cancelled() {
                break;
            }

            if let Some(item) = input.pop() {
                let ordinal = item.item.ordinal;
                reorder.push(ordinal, item.item);
                backoff.reset();

                while let Some(mut batch) = reorder.pop_ready() {
                    apply_offset_in_place(&mut batch, &mut next_mi_base)?;
                    let mem = estimate_batch_bytes(&batch);
                    let out_ordinal = batch.ordinal;
                    let seq = SequencedItem::new(out_ordinal, batch, mem);
                    if output.push_until_cancelled(seq, &cancel).is_err() {
                        anyhow::bail!("MiAssign: output queue dropped during push");
                    }
                }
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Drain remaining contiguous-prefix batches.
        while let Some(mut batch) = reorder.pop_ready() {
            if cancel.is_cancelled() {
                break;
            }
            apply_offset_in_place(&mut batch, &mut next_mi_base)?;
            let mem = estimate_batch_bytes(&batch);
            let out_ordinal = batch.ordinal;
            let seq = SequencedItem::new(out_ordinal, batch, mem);
            if output.push_until_cancelled(seq, &cancel).is_err() {
                anyhow::bail!("MiAssign: output queue dropped during final drain");
            }
        }

        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "MiAssign"
    }
}

/// Apply the cumulative offset to every `MiGroup` in `batch` and inject the
/// global MI tag bytes into each record. Advances `next_mi_base` by
/// `max(local_id) + 1` so subsequent batches see a fresh, non-overlapping
/// integer range — same accounting as standalone group's `distinct_mi_count`.
fn apply_offset_in_place(batch: &mut MiGroupBatch, next_mi_base: &mut u64) -> Result<()> {
    let mut mi_buf = String::with_capacity(16);
    let mut max_local: Option<u64> = None;

    for group in &mut batch.groups {
        let local = group
            .local_mi
            .take()
            .ok_or_else(|| anyhow::anyhow!("MiAssign: MiGroup missing local_mi"))?;

        if let Some(local_id) = local.id() {
            max_local = Some(max_local.map_or(local_id, |m| m.max(local_id)));
        }

        let global = local.with_offset(*next_mi_base);
        let mi_bytes = global.write_with_offset(0, &mut mi_buf);

        group.mi = global.id().unwrap_or(0);

        let new_data = rewrite_mi_tag(&group.data, batch.assign_tag, mi_bytes)?;
        group.data = new_data;
    }

    if let Some(max_id) = max_local {
        *next_mi_base = next_mi_base.checked_add(max_id + 1).ok_or_else(|| {
            anyhow::anyhow!("MiAssign: cumulative MoleculeId offset overflowed u64")
        })?;
    }

    Ok(())
}

/// Walk `data` as length-prefixed records, write `mi_value` into the
/// `assign_tag` field of each record, and return a fresh concat-byte buffer.
fn rewrite_mi_tag(data: &[u8], assign_tag: [u8; 2], mi_value: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(data.len() + 16);
    for rec in iter_length_prefixed(data) {
        let rec = rec?;
        let mut owned = rec.to_vec();
        update_string_tag(&mut owned, assign_tag, mi_value);
        let block_size = u32::try_from(owned.len()).map_err(|_| {
            anyhow::anyhow!("MiAssign: rewritten record exceeds u32::MAX bytes ({})", owned.len())
        })?;
        out.extend_from_slice(&block_size.to_le_bytes());
        out.extend_from_slice(&owned);
    }
    Ok(out)
}

fn estimate_batch_bytes(batch: &MiGroupBatch) -> usize {
    batch.groups.iter().map(|g| g.data.len()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::grouping_types::{MiGroup, serialize_records};
    use fgumi_raw_bam::SamBuilder;
    use fgumi_umi::MoleculeId;

    /// Build a minimal valid BAM record body (no MI tag yet).
    fn make_record_body(name: &[u8]) -> Vec<u8> {
        let mut b = SamBuilder::new();
        b.read_name(name).ref_id(-1).pos(-1).flags(0).sequence(b"A").qualities(&[30]);
        b.build().into_inner()
    }

    fn make_group(local_mi: MoleculeId, mi_numeric: u64, name: &[u8]) -> MiGroup {
        let body = make_record_body(name);
        let (data, count) = serialize_records(&[body]).unwrap();
        MiGroup { data, record_count: count, mi: mi_numeric, local_mi: Some(local_mi) }
    }

    fn make_two_group_batch(ordinal: u64) -> MiGroupBatch {
        MiGroupBatch {
            groups: vec![
                make_group(MoleculeId::Single(0), 0, b"r0"),
                make_group(MoleculeId::Single(1), 1, b"r1"),
            ],
            ordinal,
            position_key: (0, 0),
            assign_tag: *b"MI",
        }
    }

    /// Read the `MI:Z` value out of every record in a serialized `MiGroup`
    /// data buffer.
    fn extract_mi_values(group: &MiGroup) -> Vec<String> {
        iter_length_prefixed(&group.data)
            .map(|r| {
                let rec = r.expect("valid record");
                let mi =
                    fgumi_raw_bam::find_string_tag_in_record(rec, *b"MI").expect("MI tag present");
                String::from_utf8(mi.to_vec()).unwrap()
            })
            .collect()
    }

    #[test]
    fn apply_offset_advances_global_base() {
        let mut batch = make_two_group_batch(0);
        let mut base: u64 = 0;
        apply_offset_in_place(&mut batch, &mut base).unwrap();
        assert_eq!(base, 2, "max local id 1 -> base advances by 2");
        assert_eq!(batch.groups[0].mi, 0);
        assert_eq!(batch.groups[1].mi, 1);
        assert!(batch.groups[0].local_mi.is_none(), "local_mi cleared after rewrite");
        assert!(batch.groups[1].local_mi.is_none());
    }

    #[test]
    fn apply_offset_global_ids_are_offset_by_base() {
        let mut batch_a = make_two_group_batch(0);
        let mut batch_b = make_two_group_batch(1);
        let mut base: u64 = 0;
        apply_offset_in_place(&mut batch_a, &mut base).unwrap();
        apply_offset_in_place(&mut batch_b, &mut base).unwrap();
        assert_eq!(batch_a.groups[0].mi, 0);
        assert_eq!(batch_a.groups[1].mi, 1);
        assert_eq!(batch_b.groups[0].mi, 2);
        assert_eq!(batch_b.groups[1].mi, 3);
        assert_eq!(base, 4);
    }

    #[test]
    fn missing_local_mi_is_an_error() {
        let mut batch = make_two_group_batch(0);
        batch.groups[0].local_mi = None;
        let mut base: u64 = 0;
        let err = apply_offset_in_place(&mut batch, &mut base);
        assert!(err.is_err(), "missing local_mi must error, not silently succeed");
    }

    #[test]
    fn paired_variants_keep_their_tag_letters_and_share_numeric() {
        let mut batch = MiGroupBatch {
            groups: vec![
                make_group(MoleculeId::PairedA(0), 0, b"r0"),
                make_group(MoleculeId::PairedB(0), 0, b"r1"),
                make_group(MoleculeId::PairedA(1), 1, b"r2"),
            ],
            ordinal: 0,
            position_key: (0, 0),
            assign_tag: *b"MI",
        };
        let mut base: u64 = 10;
        apply_offset_in_place(&mut batch, &mut base).unwrap();
        // numeric ids advance by offset
        assert_eq!(batch.groups[0].mi, 10);
        assert_eq!(batch.groups[1].mi, 10);
        assert_eq!(batch.groups[2].mi, 11);
        // base advances by max_local (1) + 1 = 2
        assert_eq!(base, 12);
        // MI tag string preserves variant suffix
        assert_eq!(extract_mi_values(&batch.groups[0]), vec!["10/A".to_string()]);
        assert_eq!(extract_mi_values(&batch.groups[1]), vec!["10/B".to_string()]);
        assert_eq!(extract_mi_values(&batch.groups[2]), vec!["11/A".to_string()]);
    }

    #[test]
    fn mi_tag_bytes_match_global_id() {
        let mut batch = make_two_group_batch(0);
        let mut base: u64 = 100;
        apply_offset_in_place(&mut batch, &mut base).unwrap();
        assert_eq!(extract_mi_values(&batch.groups[0]), vec!["100".to_string()]);
        assert_eq!(extract_mi_values(&batch.groups[1]), vec!["101".to_string()]);
    }
}
