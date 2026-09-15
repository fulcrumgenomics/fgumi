//! Per-command chain-builder support. Each module holds the pieces
//! [`super::builder::ChainBuilder`] needs for one command's stage: a
//! `FinalizeHook` and any concrete-typed step factories.
//!
//! Stage dispatch is driven by [`super::build::build_for`], which matches on
//! `spec.stages` and drives `ChainBuilder` directly.

pub mod align;
pub mod clip;
#[cfg(feature = "consensus")]
pub mod codec;
pub mod copy_umi;
pub mod correct;
pub mod dedup;
#[cfg(feature = "consensus")]
pub mod duplex;
pub mod extract;
pub mod fastq;
pub mod filter;
pub mod group;
pub mod retag;
#[cfg(feature = "consensus")]
pub mod simplex;
pub mod sort;
pub mod zipper;

use std::io;

use fgumi_raw_bam::RawRecord;

use crate::pipeline::steps::parse::decode::validate_record_for_decode;
use crate::pipeline::steps::types::RecordBatch;

/// Drive a per-record transform over a borrowed [`RecordBatch`] without the
/// per-record heap allocation of the `DecodedRecordBatch` decode path.
///
/// For each record body (a borrowed `&[u8]` slice into the batch's shared
/// backing buffer) this: (1) validates it is a well-formed BAM record with
/// [`validate_record_for_decode`] — the fast path drops `DecodeRecords`, which
/// otherwise runs this guard, so without it a truncated record would panic in
/// the unchecked field accessors the transform closures use; (2) reloads the
/// caller's reusable `scratch` `RawRecord` (cleared, then refilled — never
/// truncated, so a shorter record inherits no tail from a longer predecessor);
/// and (3) hands `&mut scratch` to `f`, which applies the command's per-record
/// mutation and writes the framed output it wants.
///
/// The `scratch` allocation is reused across every record in the batch, so the
/// only per-record cost is the memcpy into it — the `Vec<DecodedRecord>` and the
/// per-record owned `RawRecord` of the decode path are both gone. Today's callers
/// allocate one scratch per batch (inside the `Fn` process closure); a caller
/// that persisted it across batches via `Process*WithWorkerState` would reuse the
/// one allocation for the whole run, but none do so yet.
pub(crate) fn for_each_raw_record(
    batch: &RecordBatch,
    scratch: &mut RawRecord,
    mut f: impl FnMut(&mut RawRecord) -> io::Result<()>,
) -> io::Result<()> {
    for raw in batch.iter_record_bytes() {
        validate_record_for_decode(raw)?;
        scratch.clear();
        scratch.as_mut_vec().extend_from_slice(raw);
        f(scratch)?;
    }
    Ok(())
}

/// Append `rec` to `dst` framed as a BAM record block: a 4-byte little-endian
/// length prefix followed by the record bytes. Shared by the consensus
/// command builders' rejects-serialization paths (codec/simplex/duplex), which
/// emit raw record bytes into a `DecompressedBlock` buffer.
///
/// Thin wrapper over [`fgumi_raw_bam::write_framed_record`], the canonical BAM
/// record framing (also used by `UnmappedSamBuilder::write_with_block_size`), so
/// the length-prefix arithmetic lives in exactly one place and cannot drift.
//
// Only the consensus command builders (codec/simplex/duplex) call this, so it
// is gated under `consensus` to avoid a dead-code warning when consensus is off.
#[cfg(feature = "consensus")]
pub(crate) fn append_framed_bytes(dst: &mut Vec<u8>, rec: &[u8]) -> std::io::Result<()> {
    fgumi_raw_bam::write_framed_record(dst, rec)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_bam_io::{DecodedRecord, GroupKey};
    use fgumi_raw_bam::{RawRecord, SamBuilder};

    use crate::pipeline::steps::types::{DecodedRecordBatch, RecordBatch};

    fn record(name: &[u8], seq: &[u8]) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(name).sequence(seq).qualities(&vec![30u8; seq.len()]);
        b.build()
    }

    /// `for_each_raw_record` loads each record body into the reused scratch and
    /// calls the closure once per record. A shorter record following a longer
    /// one must carry no stale tail — the scratch is cleared, not truncated.
    #[test]
    fn for_each_raw_record_reloads_scratch_without_stale_tail() {
        let long = record(b"longname", b"ACGTACGT");
        let short = record(b"s", b"AC");
        let long_bytes = long.as_ref().to_vec();
        let short_bytes = short.as_ref().to_vec();
        assert!(long_bytes.len() > short_bytes.len(), "long record must exceed short");

        let batch = RecordBatch::new(0, &[long, short]);

        let mut scratch = RawRecord::new();
        let mut seen: Vec<Vec<u8>> = Vec::new();
        for_each_raw_record(&batch, &mut scratch, |rec| {
            seen.push(rec.as_ref().to_vec());
            Ok(())
        })
        .expect("for_each_raw_record");

        assert_eq!(seen.len(), 2, "closure runs once per record");
        assert_eq!(seen[0], long_bytes, "first record bytes match the long record");
        assert_eq!(
            seen[1], short_bytes,
            "second record is exactly the short record, with no stale tail from the long one",
        );
    }

    /// A record body too short to hold the 32-byte BAM fixed header must surface
    /// as a clean `io::Error`, not reach the transform closure (whose unchecked
    /// field accessors would panic-unwind a worker). The decode path this fast
    /// path replaces ran `validate_record_for_decode` for exactly this reason.
    #[test]
    fn for_each_raw_record_errors_on_truncated_record() {
        // 8 bytes: shorter than the 32-byte fixed header.
        let batch = RecordBatch::from_parsed(0, vec![0u8; 8], vec![(0, 8)]);
        let mut scratch = RawRecord::new();

        let mut reached_closure = false;
        let result = for_each_raw_record(&batch, &mut scratch, |_rec| {
            reached_closure = true;
            Ok(())
        });

        assert!(result.is_err(), "a truncated record must surface as an io::Error");
        assert!(!reached_closure, "the transform closure must not see a malformed record");
    }

    /// The borrowed (`RecordBatch`) path frames and routes byte-identically to
    /// the owned (`DecodedRecordBatch`) path under the same per-record transform,
    /// including a length-changing mutation and a keep/reject split. The
    /// per-record transform is the *same function* on both paths, so this pins
    /// the decode-free path's iteration + framing + routing — the part that
    /// actually differs. Full command behavior stays gated by the integration
    /// suites.
    #[test]
    fn raw_path_frames_and_routes_identically_to_decoded_path() {
        // Grow the record by one byte (a stand-in for a tag-adding op — exercises
        // framing under a length change) and keep even-indexed, reject odd.
        fn mutate_and_keep(record: &mut RawRecord, idx: usize) -> bool {
            record.as_mut_vec().push(0xAB);
            idx.is_multiple_of(2)
        }

        let recs =
            vec![record(b"r0", b"ACGT"), record(b"r1", b"ACGTACGTACGT"), record(b"r2", b"A")];

        // Owned path: DecodedRecordBatch → into_records → into_raw_bytes.
        let decoded_batch = DecodedRecordBatch::new(
            0,
            recs.iter()
                .cloned()
                .map(|r| DecodedRecord::from_raw_bytes(r, GroupKey::default()))
                .collect(),
        );
        let (mut d_kept, mut d_rej) = (Vec::new(), Vec::new());
        for (idx, decoded) in decoded_batch.into_records().into_iter().enumerate() {
            let mut record = decoded.into_raw_bytes();
            let keep = mutate_and_keep(&mut record, idx);
            let target = if keep { &mut d_kept } else { &mut d_rej };
            fgumi_raw_bam::write_framed_record(target, record.as_ref()).unwrap();
        }

        // Borrowed path: RecordBatch → for_each_raw_record → reused scratch.
        let raw_batch = RecordBatch::new(0, &recs);
        let mut scratch = RawRecord::new();
        let (mut r_kept, mut r_rej) = (Vec::new(), Vec::new());
        let mut idx = 0usize;
        for_each_raw_record(&raw_batch, &mut scratch, |record| {
            let keep = mutate_and_keep(record, idx);
            idx += 1;
            let target = if keep { &mut r_kept } else { &mut r_rej };
            fgumi_raw_bam::write_framed_record(target, record.as_ref())?;
            Ok(())
        })
        .unwrap();

        assert_eq!(r_kept, d_kept, "kept stream byte-identical across paths");
        assert_eq!(r_rej, d_rej, "rejected stream byte-identical across paths");
    }
}
