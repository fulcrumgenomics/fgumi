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
//! This is vendored from noodles-bam. The upstream PR to expose codec functions
//! ([#364](https://github.com/zaeleus/noodles/pull/364)) was closed; `encode_record_buf` is
//! retained here for the raw-byte-mode optimization path. This vendored code can be removed
//! once [#367](https://github.com/zaeleus/noodles/pull/367) (optimized `RecordBuf` encoding)
//! lands upstream, or an equivalent API becomes available.

pub mod encoder;

pub use encoder::encode_record_buf;
