//! Helper utilities for integration tests.
//!
//! Individual helpers are used by different subsets of the test modules
//! (some of which are feature-gated), so some helpers appear unused under
//! certain feature combinations. Silence the resulting warnings wholesale
//! so `--no-default-features` stays clean.

#![allow(dead_code, unused_imports)]

pub mod assertions;
pub mod bam_generator;

pub use assertions::*;
pub use bam_generator::*;
