//! Typed `Step` library.
//!
//! Layout:
//!
//! - `bgzf/`                    — BGZF compress/decompress
//! - `boundaries/`              — record-boundary discovery
//! - `correct/`                 — UMI correction against a known-UMI set
//! - `decoded_to_records.rs`    — flatten decoded (SAM-parsed) records back
//!   into records, for the sort ingest
//! - `group/`                   — template grouping (BAM, position, queryname, MI)
//! - `parse/`                   — record parsing
//! - `process.rs`               — closure-driven mid-steps (`Process`,
//!   `ProcessWithWorkerState`, `MiAssign`)
//! - `serialize.rs`             — record serialization
//! - `serialize_processed.rs`   — serialization for processed-template batches
//! - `sink/`                    — write-side step implementations
//! - `sort/`                    — sort-chain step re-exports
//! - `source/`                  — read-side step implementations (BAM/SAM/FASTQ)
//! - `templates_to_records.rs`  — flatten templates back into records
//! - `tuning.rs`                — per-chain byte/queue budgets
//! - `types.rs`                 — flowing data types (`HeapSize` + `Ordered`)
//!
//! `extract.rs` and `align_and_merge.rs` land here with the chain builder
//! (R2), which is their consumer: `extract.rs` needs `ExtractOptions` (ported
//! alongside) and `align_and_merge.rs` needs `crate::aligner` (ported alongside).
//!
//! `coalesce.rs` (a `CoalesceBytes` byte-batching step with no caller anywhere,
//! upstream included) is *not* ported — it would be dead code.

pub mod align_and_merge;
pub mod bgzf;
pub mod boundaries;
#[cfg(test)]
mod chain_tests;
pub mod correct;
pub mod decoded_to_records;
pub mod extract;
pub mod group;
pub mod parse;
pub mod process;
pub mod serialize;
pub mod serialize_processed;
pub mod sink;
pub mod sort;
pub mod source;
pub mod templates_to_records;
pub mod tuning;
pub mod types;

pub use tuning::BamPipelineTuning;
pub use types::{
    BamTemplateBatch, BgzfBlock, DecodedRecordBatch, DecompressedBlock, RecordBatch,
    RecordBatchBuilder,
};

#[cfg(test)]
mod tests;
