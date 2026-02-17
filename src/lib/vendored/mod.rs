//! Vendored code that will be removed once upstream PRs are merged.
//!
//! This module contains temporary implementations that are waiting on:
//! - noodles-bgzf PR for `MultithreadedWriter` with position tracking
//! - noodles-bam PR for BAM codec optimizations
//! - noodles-bam PR for `Record::as_ref()` (raw bytes access)
//!
//! This module should be removed once the upstream PRs are merged and released.

pub mod bam_codec;
pub mod bgzf_multithreaded;
pub mod raw_bam_record;

// Re-export MultithreadedWriter types
pub use bgzf_multithreaded::{
    BlockInfo, BlockInfoRx, Builder as MultithreadedWriterBuilder, MultithreadedWriter,
};

// Re-export BAM codec functions
pub use bam_codec::{
    DecodeError, EncodeError, decode, encode, encode_record_buf, encode_with_prealloc,
};

// Re-export raw BAM record types
pub use raw_bam_record::{RawBamReader, RawRecord, read_raw_record};
