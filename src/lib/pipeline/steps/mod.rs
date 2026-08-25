//! Typed `Step` library.
//!
//! Layout:
//!
//! - `bgzf/`                    — BGZF compress/decompress
//! - `boundaries/`              — record-boundary discovery
//! - `correct/`                 — UMI correction against a known-UMI set
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
//! `extract.rs` and `align_and_merge` arrive in follow-up ports: `extract.rs`
//! still needs `ExtractOptions` from the command-layer options refactor, and
//! `align_and_merge` needs `crate::aligner`, neither of which has landed yet.
//! (`correct/` had the same dependency on `CorrectOptions`; that arrived with
//! the options projection, so it is ported here.)
//!
//! `roundtrip.rs` is a whole-chain convenience whose only consumer is the
//! `compare bam-roundtrip` command; it is ported here alongside that command so
//! it lands caller-complete rather than dormant. `coalesce.rs` (a
//! `CoalesceBytes` byte-batching step) is deliberately *not* ported: it has no
//! caller anywhere, upstream included, so it would be pure dead code.

pub mod bgzf;
pub mod boundaries;
#[cfg(test)]
mod chain_tests;
pub mod correct;
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
