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
//! This is vendored code that will be removed once upstream PRs are merged.

pub mod encoder;

pub use encoder::encode_record_buf;
