//! CLI command implementations for fgumi.
//!
//! This module contains all the command implementations for the fgumi CLI tool.
//! Each submodule implements a specific command (extract, group, filter, etc.).
//!
//! # Command Categories
//!
//! ## Pre-alignment
//! - [`extract`] - Extract UMIs from FASTQ and create unmapped BAM
//!
//! ## Post-alignment
//! - [`group`] - Group reads by UMI into molecular families
//! - [`correct`] - Correct UMI sequences against a fixed list
//! - [`filter`] - Filter reads based on various criteria
//!
//! ## Consensus Calling
//! - [`simplex`] - Call single-strand molecular consensus
//! - [`duplex`] - Call duplex (double-strand) molecular consensus
//! - [`codec`] - CODEC-style consensus calling
//!
//! ## Utilities
//! - [`clip`] - Clip overlapping read pairs
//! - [`duplex_metrics`] - Calculate duplex sequencing metrics

// Blanket clippy pedantic allows for command implementations.
// These will be removed incrementally as commands are refactored.
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::unused_self,
    clippy::unnecessary_wraps,
    clippy::similar_names,
    clippy::needless_pass_by_value,
    clippy::match_same_arms,
    clippy::must_use_candidate,
    clippy::items_after_statements,
    clippy::too_many_lines,
    clippy::fn_params_excessive_bools,
    clippy::redundant_else,
    clippy::manual_let_else,
    clippy::needless_continue,
    clippy::redundant_closure_for_method_calls,
    clippy::explicit_iter_loop,
    clippy::uninlined_format_args,
    clippy::map_unwrap_or
)]

pub mod clip;
pub mod codec;
pub mod command;
pub mod common;
#[cfg(feature = "compare")]
pub mod compare;
pub mod consensus_runner;
pub mod correct;
pub mod dedup;
pub mod downsample;
pub mod duplex;
pub mod duplex_metrics;
pub mod extract;
pub mod fastq;
pub mod filter;
pub mod group;
pub mod review;
pub mod simplex;
#[cfg(feature = "simulate")]
pub mod simulate;
pub mod sort;
pub mod zipper;
