//! Fast BAM record encoding.
//!
//! This module provides optimized functions for BAM record serialization,
//! designed for high-throughput parallel processing pipelines.
//!
//! # Encoding
//!
//! ```rust,ignore
//! use crate::vendored::bam_codec::encode_record_buf;
//! use noodles::sam::{Header, alignment::RecordBuf};
//!
//! let header = Header::default();
//! let record = RecordBuf::default();
//! let mut buf = Vec::new();
//!
//! encode_record_buf(&mut buf, &header, &record).unwrap();
//! ```
//!
//! # Decoding
//!
//! For decoding, use noodles' public `Reader` API directly:
//!
//! ```rust,ignore
//! use noodles::bam::io::Reader;
//! use noodles::sam::{Header, alignment::RecordBuf};
//!
//! let mut reader = Reader::from(&bam_data_with_block_size_prefix[..]);
//! let mut record = RecordBuf::default();
//! reader.read_record_buf(&Header::default(), &mut record)?;
//! ```
//!
//! This is vendored from noodles-bam. The upstream PR to expose codec functions
//! ([#364](https://github.com/zaeleus/noodles/pull/364)) was closed; the maintainer recommended
//! `Reader::from(src)` for decoding (now used). The vendored `encode_record_buf` is still needed
//! for the raw-byte-mode optimization path. See [`super`] module docs for full context.

pub mod encoder;

pub use encoder::{EncodeError, encode, encode_record_buf, encode_with_prealloc};
