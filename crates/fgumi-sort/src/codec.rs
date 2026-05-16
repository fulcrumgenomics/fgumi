//! Codec selection for spill chunks.
//!
//! Spill files are written by the worker pool in one of two codecs:
//! - **BGZF** (legacy, libdeflate level 1): block-stream with a 1f 8b magic,
//!   compatible with the existing pool-integrated reader.
//! - **Zstd** (default since v0.3): a sequence of length-prefixed zstd frames
//!   preceded by the 4-byte "ZSP1" file magic. Zstd-1 is ~3x faster than
//!   libdeflate-1 to compress and ~6x faster to decompress at a comparable
//!   ratio on BAM-record data.
//!
//! Output BAM files are always BGZF — that's mandated by the BAM spec — and
//! continue to use the separate output writer path.

/// Magic bytes at the start of a zstd-spill file: ASCII "ZSP1".
pub(crate) const ZSPILL_MAGIC: [u8; 4] = *b"ZSP1";

/// First two bytes of a BGZF/gzip stream.
pub(crate) const BGZF_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Codec used to compress spill chunks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SpillCodec {
    /// BGZF blocks with libdeflate. Compatible with all existing tools.
    Bgzf,
    /// Zstd frames with `[u32 LE compressed-len][frame]` records, prefixed
    /// by the file magic `ZSPILL_MAGIC`. Faster than BGZF on every axis.
    #[default]
    Zstd,
}

impl SpillCodec {
    /// Identify the codec implied by a file's first bytes.
    ///
    /// Returns `None` if no known magic matches; callers may treat that as
    /// an uncompressed stream.
    #[must_use]
    pub fn from_magic(bytes: &[u8]) -> Option<Self> {
        if bytes.starts_with(&ZSPILL_MAGIC) {
            Some(Self::Zstd)
        } else if bytes.starts_with(&BGZF_MAGIC) {
            Some(Self::Bgzf)
        } else {
            None
        }
    }
}

impl std::str::FromStr for SpillCodec {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "bgzf" | "gzip" => Ok(Self::Bgzf),
            "zstd" | "zst" => Ok(Self::Zstd),
            other => Err(format!("invalid spill codec '{other}' (expected: bgzf | zstd)")),
        }
    }
}

impl std::fmt::Display for SpillCodec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Bgzf => f.write_str("bgzf"),
            Self::Zstd => f.write_str("zstd"),
        }
    }
}
