#![deny(unsafe_code)]

pub mod builder;
pub mod cigar;
pub mod fields;
pub mod overlap;
pub mod raw_bam_record;
pub mod sequence;
pub mod sort;
pub mod tags;

#[cfg(any(test, feature = "test-utils"))]
pub mod testutil;

#[cfg(feature = "noodles")]
pub mod noodles_compat;

// Flat re-exports — callers use fgumi_raw_bam::flags() etc.
pub use builder::*;
pub use cigar::*;
pub use fields::*;
pub use overlap::*;
pub use raw_bam_record::*;
pub use sequence::*;
pub use sort::*;
pub use tags::*;

#[cfg(any(test, feature = "test-utils"))]
pub use testutil::*;

#[cfg(feature = "noodles")]
pub use noodles_compat::*;
