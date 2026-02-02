//! Fast BAM record encoding and decoding.
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
//! ```rust,ignore
//! use crate::vendored::bam_codec::decode;
//! use noodles::sam::alignment::RecordBuf;
//!
//! // Raw BAM record bytes (example: unmapped record with default fields)
//! let bam_bytes: &[u8] = &[
//!     0xff, 0xff, 0xff, 0xff, // ref_id = -1
//!     0xff, 0xff, 0xff, 0xff, // pos = -1
//!     0x02, // l_read_name = 2
//!     0xff, // mapq = 255
//!     0x48, 0x12, // bin = 4680
//!     0x00, 0x00, // n_cigar_op = 0
//!     0x04, 0x00, // flag = 4 (UNMAPPED)
//!     0x00, 0x00, 0x00, 0x00, // l_seq = 0
//!     0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
//!     0xff, 0xff, 0xff, 0xff, // next_pos = -1
//!     0x00, 0x00, 0x00, 0x00, // tlen = 0
//!     0x2a, 0x00, // read_name = "*\x00"
//! ];
//!
//! let mut record = RecordBuf::default();
//! decode(&mut &bam_bytes[..], &mut record).unwrap();
//! ```
//!
//! This is vendored code that will be removed once upstream PRs are merged.

pub mod decoder;
pub mod encoder;

pub use decoder::{DecodeError, decode};
pub use encoder::{EncodeError, encode, encode_record_buf, encode_with_prealloc};
