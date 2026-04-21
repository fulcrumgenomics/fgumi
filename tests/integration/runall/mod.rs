//! Integration tests for the `runall` pipeline (`pipeline` engine).

#[cfg(feature = "compare")]
pub mod test_consensus_codec;
#[cfg(feature = "compare")]
pub mod test_consensus_duplex;
#[cfg(feature = "compare")]
pub mod test_consensus_simplex;
#[cfg(feature = "compare")]
pub mod test_correct;
#[cfg(feature = "compare")]
pub mod test_extract;
#[cfg(feature = "compare")]
pub mod test_fastq;
#[cfg(feature = "compare")]
pub mod test_group;
#[cfg(feature = "compare")]
pub mod test_sort;
#[cfg(feature = "compare")]
pub mod test_zipper;
