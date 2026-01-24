//! BAM record codec for encoding and decoding.

pub mod decoder;
pub mod encoder;

// Re-export key items for public use
pub use self::decoder::{decode, DecodeError};
pub use self::encoder::{encode, encode_record_buf, encode_with_prealloc};
