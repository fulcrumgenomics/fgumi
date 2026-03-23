//! SIMD-accelerated FASTQ parsing using Helicase-style bitmask operations.
//!
//! This crate provides high-throughput FASTQ parsing by processing 64 bytes at a time
//! through SIMD registers (NEON on ARM, AVX2 on `x86_64`), classifying newline characters
//! via bitmask operations and finding record boundaries without per-byte branching.
//!
//! # Architecture
//!
//! 1. **Lexer**: Loads 64-byte blocks into SIMD registers, produces a `u64` bitmask where
//!    bit `i` is set if byte `i` is a newline (`\n`).
//! 2. **Parser**: Walks the newline bitmask with `trailing_zeros()` to find record
//!    boundaries. Every 4th newline marks the end of a FASTQ record.
//!
//! # Example
//!
//! ```
//! use fgumi_simd_fastq::{find_record_offsets, parse_records};
//!
//! let fastq = b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
//! let offsets = find_record_offsets(fastq);
//! assert_eq!(offsets, vec![0, 16, 32]);
//!
//! let records: Vec<_> = parse_records(fastq).collect();
//! assert_eq!(records.len(), 2);
//! assert_eq!(records[0].name, b"r1");
//! assert_eq!(records[0].sequence, b"ACGT");
//! ```

mod bitmask;
mod lexer;
mod parser;
mod reader;

pub use bitmask::FastqBitmask;
pub use lexer::lex_block_full;
pub use parser::{FastqRecord, find_record_offsets, parse_records};
pub use reader::SimdFastqReader;
