//! High-performance external merge-sort engine for BAM files.
//!
//! This crate provides the sort engine that powers `fgumi sort`. It sorts BAM files
//! by coordinate, queryname (lexicographic or natural), or template-coordinate, with
//! an external merge-sort pipeline that streams chunks through an N+2 worker pool,
//! parallel radix-sorts in memory, spills to disk under memory pressure, and
//! K-way merges sorted chunks via a loser tree.

#![deny(missing_docs)]
#![deny(unsafe_code)]
