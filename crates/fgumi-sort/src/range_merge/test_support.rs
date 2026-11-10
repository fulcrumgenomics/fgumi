//! Test helpers for range merges: indexed spill runs written to disk and the
//! serial-merge oracle.
//!
//! Compiled for this crate's tests and under the `test-utils` feature, so
//! other crates' tests can build indexed runs too.

use std::path::Path;
use std::sync::Arc;

use crate::codec::SpillCodec;
use crate::keys::RawSortKey;
use crate::merge_driver::open_spill_slot;
use crate::merge_slots::SortMergeSlot;
use crate::run_index::{RunIndex, frame_meta, spill_key_width};
use crate::spill_block::{
    SpillBlockCompressor, SpillKeyKind, frame_keyed_record_into, spill_magic, spill_trailer,
};

/// Write `records` (sorted by key) as a zstd spill run `run-<file_id>.zsp`
/// under `dir`, open it as a slot with `file_id`, and attach its run index.
///
/// Frames are cut with the spill gather's rule: a frame takes records while it
/// stays at most `frame_bytes` framed bytes; the record that would push a
/// non-empty frame past that starts the next frame; a lone oversized record
/// forms its own frame; a frame reaching `frame_bytes` exactly closes.
///
/// # Panics
///
/// If `kind` is not indexed, `K` does not match `kind`'s framing, or any I/O
/// fails.
#[must_use]
pub fn write_indexed_run<K: RawSortKey>(
    dir: &Path,
    file_id: u32,
    kind: SpillKeyKind,
    records: &[(K, Vec<u8>)],
    frame_bytes: usize,
) -> Arc<SortMergeSlot> {
    let width = spill_key_width(kind).expect("an indexed spill key kind");
    let key_bytes = kind.serialized_key_bytes();
    let mut compressor = SpillBlockCompressor::new(SpillCodec::Zstd, 1).expect("zstd compressor");
    let mut file = spill_magic(SpillCodec::Zstd).to_vec();
    let mut index = RunIndex::new(width);
    let mut raw = Vec::new();
    let mut unit = Vec::new();
    let mut first_key = Vec::new();

    let mut flush = |raw: &mut Vec<u8>, file: &mut Vec<u8>, index: &mut RunIndex| {
        if raw.is_empty() {
            return;
        }
        unit.clear();
        compressor.compress_block_into(raw, &mut unit).expect("compress frame");
        let count = frame_meta(kind, raw, &mut first_key).expect("frame metadata");
        index.push_frame(
            u32::try_from(unit.len()).expect("frame fits u32"),
            u32::try_from(raw.len()).expect("raw frame fits u32"),
            count,
            &first_key,
        );
        file.extend_from_slice(&unit);
        raw.clear();
    };

    for (key, body) in records {
        let framed = key_bytes + 4 + body.len();
        if !raw.is_empty() && raw.len() + framed > frame_bytes {
            flush(&mut raw, &mut file, &mut index);
        }
        frame_keyed_record_into(&mut raw, key, body).expect("frame record");
        if raw.len() >= frame_bytes {
            flush(&mut raw, &mut file, &mut index);
        }
    }
    flush(&mut raw, &mut file, &mut index);
    file.extend_from_slice(spill_trailer(SpillCodec::Zstd));

    let path = dir.join(format!("run-{file_id}.zsp"));
    std::fs::write(&path, &file).expect("write spill run");
    index.rebase(spill_magic(SpillCodec::Zstd).len() as u64);
    let slot = open_spill_slot(&path, file_id, crate::FramesEnd::AtEof).expect("open spill slot");
    slot.attach_run_index(Arc::new(index));
    slot
}

/// The serial k-way merge of `sources` (each sorted by key, in leaf order):
/// every record's body, ordered by `(key, source index, position in source)`,
/// the loser tree's tie rule.
#[must_use]
pub fn serial_reference<K: Ord + Clone>(sources: &[Vec<(K, Vec<u8>)>]) -> Vec<Vec<u8>> {
    let mut all: Vec<(K, usize, usize, &[u8])> = sources
        .iter()
        .enumerate()
        .flat_map(|(s, recs)| {
            recs.iter().enumerate().map(move |(p, (k, body))| (k.clone(), s, p, body.as_slice()))
        })
        .collect();
    all.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    all.into_iter().map(|(_, _, _, body)| body.to_vec()).collect()
}

/// `bodies` framed back to back as `[u32 LE len][body]`.
///
/// # Panics
///
/// If a body is longer than `u32::MAX`.
#[must_use]
pub fn frame_bodies(bodies: &[Vec<u8>]) -> Vec<u8> {
    let mut out = Vec::new();
    for body in bodies {
        out.extend_from_slice(&u32::try_from(body.len()).expect("body fits u32").to_le_bytes());
        out.extend_from_slice(body);
    }
    out
}
