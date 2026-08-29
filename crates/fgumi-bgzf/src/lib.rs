#![deny(unsafe_code)]

//! BGZF (Blocked GZIP Format) reading and writing utilities.
//!
//! This crate provides low-level BGZF block I/O:
//! - [`header`] - The block-header layout every fgumi reader agrees on
//! - [`reader`] - Raw block reading and decompression using libdeflater
//! - [`writer`] - Inline BGZF compression using the `bgzf` crate

pub mod header;
pub mod reader;
pub mod writer;

// Re-export commonly used types
pub use header::{
    BGZF_HEADER_SIZE, HeaderRejection, MIN_BLOCK_SIZE, block_size, block_size_checked,
    is_bgzf_header,
};
pub use reader::{
    BGZF_EOF, BGZF_FOOTER_SIZE, MAX_UNCOMPRESSED_BLOCK_SIZE, RawBgzfBlock, decompress_block,
    decompress_block_into, decompress_block_into_opts, decompress_block_slice_into,
    decompress_block_slice_into_opts, decompress_into_slice, read_raw_blocks, uncompressed_size,
};
// Re-export the libdeflater decompressor so downstream crates can name the type
// required by `decompress_block_into_opts` without depending on libdeflater
// directly.
pub use libdeflater::Decompressor;
pub use writer::{BGZF_MAX_BLOCK_SIZE, CompressedBlock, InlineBgzfCompressor};
