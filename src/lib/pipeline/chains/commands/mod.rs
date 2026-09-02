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
