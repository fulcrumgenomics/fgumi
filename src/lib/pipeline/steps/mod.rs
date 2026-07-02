//! Typed `Step` library. Phase 3 introduces the BAM step library; Phase 4+
//! adds command-specific compositions; Phase 5 adds FASTQ steps.
//!
//! Layout (per Phase 0 design lines 354–374):
//!
//! - `source/`     — read-side step implementations
//! - `bgzf/`       — BGZF compress/decompress
//! - `boundaries/` — record-boundary discovery
//! - `coalesce.rs` — coalesce small byte-batches into larger ones (BGZF amortization)
//! - `parse/`      — record parsing
//! - `group/`      — template grouping
//! - `process.rs`  — closure-driven mid-steps (`Process`, `ProcessWithWorkerState`, `MiAssign`)
//! - `serialize.rs` — record serialization
//! - `sink/`       — write-side step implementations
//! - `types.rs`    — flowing data types (`HeapSize` + `Ordered`)

pub mod align_and_merge;
pub mod bgzf;
pub mod boundaries;
pub mod coalesce;
pub mod correct;
pub mod extract;
pub mod group;
pub mod parse;
pub mod process;
pub mod roundtrip;
pub mod serialize;
pub mod serialize_processed;
pub mod sink;
pub mod sort;
pub mod source;
pub mod templates_to_records;
pub mod tuning;
pub mod types;

pub use roundtrip::{RoundtripConfig, run_bam_roundtrip};
pub use tuning::BamPipelineTuning;
pub use types::{
    BamTemplateBatch, BgzfBlock, DecodedRecordBatch, DecompressedBlock, RecordBatch,
    RecordBatchBuilder,
};

#[cfg(test)]
mod tests;
