//! Vendored code that will be removed once upstream PRs are merged.
//!
//! This module contains temporary implementations that are waiting on:
//! - noodles-bgzf: `MultithreadedWriter` with position tracking (<https://github.com/zaeleus/noodles/pull/371>)
//! - noodles-bam: expose record codec encode/decode (<https://github.com/zaeleus/noodles/pull/364>)
//! - noodles-bam: optimized `RecordBuf` encoding (<https://github.com/zaeleus/noodles/pull/367>)
//! - noodles-bam: `Record::as_ref()` raw bytes access (<https://github.com/zaeleus/noodles/pull/373>)
//!
//! This module should be removed once the upstream PRs are merged and released.

pub mod bam_codec;
pub mod bgzf_multithreaded;

// Re-export MultithreadedWriter types
pub use bgzf_multithreaded::{
    BlockInfo, BlockInfoRx, Builder as MultithreadedWriterBuilder, MultithreadedWriter,
};

// Re-export BAM codec functions
pub use bam_codec::{
    DecodeError, EncodeError, decode, encode, encode_record_buf, encode_with_prealloc,
};
