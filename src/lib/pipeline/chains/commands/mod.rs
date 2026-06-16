//! Per-command chain builders. Each module exports a single
//! `build_<command>_chain(spec: ChainSpec) -> Result<BuiltPipeline>`
//! function. [`super::build`] dispatches to these by matching on
//! `spec.stages`.
//!
//! The bodies are extracted from each command's existing inline
//! chain construction during Phase 2 (T2.12–T2.22 = move). Phase 3
//! will inline these into `super::build`'s match arms directly.

pub mod align;
pub mod clip;
pub mod codec;
pub mod correct;
pub mod dedup;
pub mod duplex;
pub mod extract;
pub mod filter;
pub mod group;
pub mod simplex;
pub mod sort;
pub mod zipper;

/// Append `rec` to `dst` framed as a BAM record block: a 4-byte little-endian
/// length prefix followed by the record bytes. Shared by the consensus
/// command builders' rejects-serialization paths (codec/simplex/duplex), which
/// emit raw record bytes into a `DecompressedBlock` buffer.
///
/// BAM record body size is u32-bounded per the spec; a single record's buffer
/// cannot exceed `u32::MAX` in practice, so the length cast cannot truncate.
pub(crate) fn append_framed_bytes(dst: &mut Vec<u8>, rec: &[u8]) {
    #[allow(clippy::cast_possible_truncation)]
    let block_size = rec.len() as u32;
    dst.extend_from_slice(&block_size.to_le_bytes());
    dst.extend_from_slice(rec);
}
