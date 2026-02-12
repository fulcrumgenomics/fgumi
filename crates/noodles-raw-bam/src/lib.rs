#![deny(unsafe_code)]

pub mod fields;
pub mod tags;
pub mod sequence;
pub mod cigar;
pub mod sort;
pub mod overlap;
pub mod builder;

#[cfg(any(test, feature = "test-utils"))]
pub mod testutil;

#[cfg(feature = "noodles")]
pub mod noodles_compat;

// Flat re-exports — callers use noodles_raw_bam::flags() etc.
pub use fields::*;
pub use tags::*;
pub use sequence::*;
pub use cigar::*;
pub use sort::*;
pub use overlap::*;
pub use builder::*;

#[cfg(any(test, feature = "test-utils"))]
pub use testutil::*;

#[cfg(feature = "noodles")]
pub use noodles_compat::*;
