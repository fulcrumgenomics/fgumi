//! BAM I/O factories and pipeline primitives for the fgumi family of crates.
//!
//! This crate provides the layer that bolts [`noodles`]' BAM record types onto fgumi's
//! faster libdeflate-backed BGZF I/O, plus a small set of pipeline primitives
//! (`ProgressTracker`, `ReorderBuffer`, `MemoryEstimate`) used to consume BAM streams
//! efficiently in parallel.
//!
//! Use this crate when you need to read or write BAM files at high throughput in a
//! pipeline that you want to keep independent of the rest of the fgumi binary.

#![deny(missing_docs)]
#![deny(unsafe_code)]

pub mod bam_io;
pub mod progress;
