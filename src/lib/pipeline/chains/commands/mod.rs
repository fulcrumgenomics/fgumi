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
#[cfg(feature = "consensus")]
pub mod codec;
pub mod correct;
pub mod dedup;
#[cfg(feature = "consensus")]
pub mod duplex;
pub mod extract;
pub mod filter;
pub mod group;
#[cfg(feature = "consensus")]
pub mod simplex;
pub mod sort;
pub mod zipper;

/// Append `rec` to `dst` framed as a BAM record block: a 4-byte little-endian
/// length prefix followed by the record bytes. Shared by the consensus
/// command builders' rejects-serialization paths (codec/simplex/duplex), which
/// emit raw record bytes into a `DecompressedBlock` buffer.
///
/// BAM record body size is u32-bounded per the spec and the canonical builder
/// (`fgumi-raw-bam`) already panics on overflow, so a single record's buffer
/// cannot exceed `u32::MAX` in practice. The length prefix is nonetheless
/// written through a checked `u32::try_from` (mirroring `filter.rs` and the
/// builder) so the length cannot silently truncate if that invariant is ever
/// broken — a loud `InvalidData` error is strictly safer than a corrupt frame.
//
// Only the consensus command builders (codec/simplex/duplex) call this, so it
// is gated under `consensus` to avoid a dead-code warning when consensus is off.
#[cfg(feature = "consensus")]
pub(crate) fn append_framed_bytes(dst: &mut Vec<u8>, rec: &[u8]) -> std::io::Result<()> {
    let block_size = u32::try_from(rec.len()).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("BAM record too large ({} bytes)", rec.len()),
        )
    })?;
    dst.extend_from_slice(&block_size.to_le_bytes());
    dst.extend_from_slice(rec);
    Ok(())
}
