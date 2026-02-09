#![deny(unsafe_code)]
// Clippy lint configuration for CI
// These lints are allowed because:
// - cast_*: Scientific/bioinformatics code intentionally casts between numeric types
// - missing_*_doc: Documentation improvements tracked separately
// - needless_pass_by_value: Some APIs designed for ownership transfer
// - items_after_statements: Some test code uses late item declarations
// - unused_self: Trait implementations may not use self
// - match_same_arms: Sometimes clearer to list arms explicitly
// - unnecessary_wraps: Some Result returns are for API consistency
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::needless_pass_by_value,
    clippy::items_after_statements,
    clippy::unused_self,
    clippy::match_same_arms,
    clippy::unnecessary_wraps,
    clippy::too_many_lines,
    clippy::redundant_closure_for_method_calls,
    clippy::explicit_iter_loop,
    clippy::struct_excessive_bools,
    clippy::map_unwrap_or,
    clippy::uninlined_format_args
)]

//! # fgumi - Fulcrum Genomics UMI Tools Library
//!
//! This library provides core functionality for working with Unique Molecular Identifiers (UMIs)
//! in sequencing data, including grouping, consensus calling, and quality filtering.
//!
//! ## Overview
//!
//! The fgumi library is organized into several key modules:
//!
//! ### Core Functionality
//!
//! - **[`umi`]** - UMI assignment strategies (identity, edit-distance, adjacency, paired)
//! - **[`consensus`]** - Consensus calling algorithms (simplex, duplex, vanilla)
//! - **[`sam`]** - SAM/BAM file utilities and alignment tag manipulation
//!
//! ### Utilities
//!
//! - **[`bam_io`]** - BAM file I/O helpers for reading and writing
//! - **[`validation`]** - Input validation utilities for parameters and files
//! - **[`progress`]** - Progress tracking and logging
//! - **[`logging`]** - Enhanced logging utilities with formatting
//! - **[`metrics`]** - Structured metrics types and file writing utilities
//! - **[`rejection`]** - Rejection reason tracking and statistics
//!
//! ### Specialized Modules
//!
//! - **[`clipper`]** - Read clipping for overlapping pairs
//! - **[`template`]** - Template-based read grouping
//! - **[`reference`][mod@reference]** - Reference genome handling
//!
//! ## Quick Start
//!
//! ### Reading and Writing BAM Files
//!
//! ```no_run
//! use fgumi_lib::bam_io::{create_bam_reader, create_bam_writer};
//!
//! # fn main() -> anyhow::Result<()> {
//! // Open input BAM and get header (path, threads)
//! let (mut reader, header) = create_bam_reader("input.bam", 1)?;
//!
//! // Create output BAM writer (path, header, threads, compression_level)
//! let mut writer = create_bam_writer("output.bam", &header, 1, 6)?;
//! # Ok(())
//! # }
//! ```
//!
//! ### Validating Input Files
//!
//! ```no_run
//! use fgumi_lib::validation::validate_file_exists;
//!
//! # fn main() -> anyhow::Result<()> {
//! // Validate input files exist with clear error messages
//! validate_file_exists("input.bam", "Input BAM")?;
//! validate_file_exists("reference.fa", "Reference FASTA")?;
//! # Ok(())
//! # }
//! ```
//!
//! ### Progress Tracking
//!
//! ```no_run
//! use fgumi_lib::progress::ProgressTracker;
//!
//! # fn main() -> anyhow::Result<()> {
//! let tracker = ProgressTracker::new("Processing records")
//!     .with_interval(100);
//!
//! for _i in 0..1000 {
//!     // Process one record...
//!     tracker.log_if_needed(1);  // Track incremental progress
//! }
//! tracker.log_final();  // Log final count if not exactly on interval
//! # Ok(())
//! # }
//! ```
//!
//! ### UMI Assignment
//!
//! ```
//! use fgumi_lib::umi::{IdentityUmiAssigner, UmiAssigner};
//!
//! let assigner = IdentityUmiAssigner::default();
//! let umis = vec!["ACGTACGT".to_string(), "ACGTACGT".to_string(), "TGCATGCA".to_string()];
//! let assignments = assigner.assign(&umis);
//! // With identity assignment, each unique UMI gets its own molecule ID
//! // So we have 2 unique molecule IDs (ACGTACGT and TGCATGCA)
//! assert_eq!(assignments.iter().collect::<std::collections::HashSet<_>>().len(), 2);
//! ```
//!
//! ## Feature Highlights
//!
//! - **Type-safe BAM I/O** - Headers always paired with readers
//! - **Consistent validation** - Standardized error messages
//! - **Progress tracking** - Uniform logging across tools
//! - **Module organization** - Related functionality grouped logically
//! - **Comprehensive testing** - Extensive test suite ensuring correctness
//!
//! ## Architecture
//!
//! The library follows these design principles:
//!
//! - **Separation of concerns** - Modules have clear, focused responsibilities
//! - **Backward compatibility** - Re-exports maintain existing APIs
//! - **Testability** - Comprehensive unit and integration tests
//! - **Documentation** - All public items documented with examples
//!
//! ## Contributing
//!
//! When adding new functionality:
//!
//! 1. Add to appropriate module group (sam, umi, consensus, etc.)
//! 2. Include comprehensive documentation and examples
//! 3. Add unit tests covering edge cases
//! 4. Maintain backward compatibility via re-exports
//!
//! ## See Also
//!
//! - [fgbio](https://github.com/fulcrumgenomics/fgbio) - Scala implementation
//! - [noodles](https://github.com/zaeleus/noodles) - Rust bioinformatics I/O

pub mod bam_io;
pub mod batched_sam_reader;
pub mod bgzf_reader;
pub mod bgzf_writer;
pub mod bitenc;
pub mod clipper;
pub mod consensus;
pub mod dna;
pub mod errors;
pub mod fastq;
pub mod grouper;
pub mod header;
pub mod logging;
pub mod metrics;
pub mod mi_group;
pub mod phred;
pub mod progress;
pub mod read_info;
pub mod reference;
pub mod rejection;
pub mod reorder_buffer;
pub mod sam;
pub mod sort;
pub mod tag_reversal;
pub mod template;
pub mod umi;
pub mod unified_pipeline;
pub mod validation;
pub mod variant_review;
#[doc(hidden)]
pub mod vendored;

#[cfg(feature = "simulate")]
pub mod simulate;

// Re-export rejection tracking types for convenient access
pub use rejection::RejectionReason;

// Re-export commonly used SAM items for backward compatibility
pub use sam::alignment_tags;

// Re-export UMI items for backward compatibility
pub use umi::assigner;

// Re-export consensus items for backward compatibility
pub use consensus::caller as consensus_caller;
pub use consensus::duplex_caller as duplex_consensus_caller;
pub use consensus::filter as consensus_filter;
pub use consensus::overlapping as overlapping_consensus;
pub use consensus::simple_umi as simple_umi_consensus;
pub use consensus::tags as consensus_tags;
pub use consensus::vanilla_caller as vanilla_consensus_caller;
