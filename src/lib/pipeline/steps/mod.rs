//! Typed `Step` library.
//!
//! Layout:
//!
//! - `bgzf/`       — BGZF compress/decompress
//! - `boundaries/` — record-boundary discovery
//! - `parse/`      — record parsing
//! - `sink/`       — write-side step implementations
//! - `source/`     — read-side step implementations (BAM/SAM/FASTQ)
//! - `sort/`       — sort-chain step re-exports
//! - `tuning.rs`   — per-chain byte/queue budgets
//! - `types.rs`    — flowing data types (`HeapSize` + `Ordered`)
//!
//! `group/`, `correct/`, and the closure-driven mid-steps arrive in follow-up
//! ports.

pub mod bgzf;
pub mod boundaries;
#[cfg(test)]
mod chain_tests;
pub mod parse;
pub mod sink;
pub mod sort;
pub mod source;
pub mod tuning;
pub mod types;

pub use tuning::BamPipelineTuning;
pub use types::{
    BamTemplateBatch, BgzfBlock, DecodedRecordBatch, DecompressedBlock, RecordBatch,
    RecordBatchBuilder,
};
