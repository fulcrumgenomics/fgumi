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
/// The libdeflater decompressor every `decompress_*` entry point here takes.
///
/// Re-exported so a consumer can name the type it must construct and reuse
/// without also declaring a direct `libdeflater` dependency — and so the
/// version it names is necessarily the one this crate decompresses with.
pub use libdeflater::Decompressor;
pub use reader::{
    BGZF_EOF, BGZF_FOOTER_SIZE, RawBgzfBlock, decompress_block, decompress_block_into,
    decompress_block_slice_into, decompress_into_slice, read_raw_blocks, uncompressed_size,
};
pub use writer::{BGZF_MAX_BLOCK_SIZE, CompressedBlock, InlineBgzfCompressor};
