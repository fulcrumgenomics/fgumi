//! Simulation utilities for generating synthetic sequencing data.
//!
//! This module provides utilities for generating synthetic FASTQ and BAM files
//! for benchmarking and testing the fgumi pipeline.
//!
//! # Modules
//!
//! - [`rng`] - Seeded random number generator utilities
//! - [`quality`] - Position-dependent quality score models
//! - [`insert_size`] - Insert size distribution models
//! - [`family_size`] - Family size distribution models
//! - [`strand_bias`] - A/B strand ratio models for duplex
//! - [`fastq_writer`] - Gzipped FASTQ file writer

pub mod family_size;
pub mod fastq_writer;
pub mod insert_size;
pub mod parallel_gzip_writer;
pub mod quality;
pub mod rng;
pub mod strand_bias;

pub use family_size::FamilySizeDistribution;
pub use fastq_writer::FastqWriter;
pub use insert_size::InsertSizeModel;
pub use quality::{PositionQualityModel, ReadPairQualityBias};
pub use rng::create_rng;
pub use strand_bias::StrandBiasModel;
