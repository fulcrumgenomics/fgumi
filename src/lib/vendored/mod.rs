//! Vendored code that will be removed once upstream PRs are merged.
//!
//! This module contains temporary implementations that are waiting on:
//! - noodles-bgzf: `MultithreadedWriter` with position tracking (<https://github.com/zaeleus/noodles/pull/371>)
//! - noodles-bam: optimized `RecordBuf` encoding (<https://github.com/zaeleus/noodles/pull/367>)
//! - noodles-bam: `Record::as_ref()` raw bytes access (<https://github.com/zaeleus/noodles/pull/373>)
//!
//! ## BAM codec (`bam_codec`)
//!
//! The vendored BAM codec provides `encode_record_buf` used by the raw-byte-mode optimization
//! path that avoids full `RecordBuf` round-tripping (~30% CPU savings).
//!
//! PR [#364](https://github.com/zaeleus/noodles/pull/364) proposed exposing codec encode/decode
//! upstream but was closed — the maintainer recommended `Reader::from(src)` for decoding instead.
//! The decode side now uses the noodles public API. The vendored `encode_record_buf` is still
//! needed for the encode side of raw-byte-mode. The vendored encoder will remain until either:
//! - PR [#367](https://github.com/zaeleus/noodles/pull/367) lands (providing optimized encoding upstream), or
//! - an alternative upstream API for direct `RecordBuf` encoding becomes available.
//!
//! The remaining vendored modules (`bgzf_multithreaded`) should be removed once their
//! respective upstream PRs are merged and released.

pub mod bam_codec;
pub mod bgzf_multithreaded;

// Re-export MultithreadedWriter types
pub use bgzf_multithreaded::{
    BlockInfo, BlockInfoRx, Builder as MultithreadedWriterBuilder, MultithreadedWriter,
};

// Re-export BAM codec functions
pub use bam_codec::{EncodeError, encode, encode_record_buf, encode_with_prealloc};
