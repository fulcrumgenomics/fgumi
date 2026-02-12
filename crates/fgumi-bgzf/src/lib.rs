#![deny(unsafe_code)]

//! BGZF (Blocked GZIP Format) reading and writing utilities.
//!
//! This crate provides low-level BGZF block I/O:
//! - [`reader`] - Raw block reading and decompression using libdeflater
//! - [`writer`] - Inline BGZF compression using the `bgzf` crate

pub mod reader;
pub mod writer;

// Re-export commonly used types
pub use reader::{
    BGZF_EOF, BGZF_FOOTER_SIZE, BGZF_HEADER_SIZE, RawBgzfBlock, decompress_block,
    decompress_block_into, decompress_block_slice_into, read_raw_blocks,
};
pub use writer::{CompressedBlock, InlineBgzfCompressor};
