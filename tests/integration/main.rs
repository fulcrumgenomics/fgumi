#![deny(unsafe_code)]

//! Integration tests for fgumi library.
//!
//! These tests validate end-to-end workflows that span multiple modules,
//! ensuring that module interactions work correctly.

mod helpers;
mod test_async_reader;
mod test_bam_pipeline;
mod test_bgzf_eof;
mod test_clip_command;
#[cfg(feature = "codec")]
mod test_codec_command;
#[cfg(feature = "codec")]
mod test_codec_pipeline;
#[cfg(feature = "compare")]
mod test_compare_bams;
mod test_correct_command;
mod test_dedup_command;
mod test_downsample_command;
#[cfg(feature = "duplex")]
mod test_duplex_command;
#[cfg(feature = "duplex")]
mod test_duplex_metrics_command;
#[cfg(all(feature = "compare", feature = "simulate"))]
mod test_e2e_regression;
mod test_error_paths;
mod test_extract_command;
mod test_fastq_command;
mod test_filter_command;
mod test_group_command;
mod test_group_determinism;
mod test_pipeline_concurrency;
mod test_review_command;
#[cfg(feature = "simplex")]
mod test_simplex_command;
#[cfg(feature = "simplex")]
mod test_simplex_metrics_command;
#[cfg(feature = "simplex")]
mod test_simplex_pipeline;
mod test_simulate_sort;
mod test_sort_correctness;
mod test_sort_write_index;
mod test_streaming_input;
mod test_zipper_command;
