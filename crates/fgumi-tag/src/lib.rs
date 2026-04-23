#![deny(unsafe_code)]

//! Validated SAM/BAM aux-tag identifier — see [`SamTag`].
//!
//! This crate exists to break a dependency cycle between `fgumi-sam` and
//! `fgumi-raw-bam`. Both need access to `SamTag`; putting it here keeps
//! either crate from depending on the other.

mod tag;
pub use tag::SamTag;
