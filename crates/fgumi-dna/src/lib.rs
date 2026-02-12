#![deny(unsafe_code)]

//! DNA sequence utilities and constants.
//!
//! This crate provides fundamental DNA operations:
//! - Complement and reverse-complement of DNA sequences
//! - 2-bit encoding for efficient UMI comparison (Hamming distance)
//! - Common constants for base quality and no-call handling

pub mod bitenc;
pub mod dna;

// Re-export submodule contents at crate root for convenience
pub use bitenc::BitEnc;
pub use dna::{complement_base, reverse_complement, reverse_complement_str};

/// No-call base character (matches fgbio's `NoCallBase`).
pub const NO_CALL_BASE: u8 = b'N';

/// Lowercase no-call base character.
pub const NO_CALL_BASE_LOWER: u8 = b'n';

/// Minimum Phred score (Q2, matching fgbio's `PhredScore.MinValue`).
pub const MIN_PHRED: u8 = 2;
