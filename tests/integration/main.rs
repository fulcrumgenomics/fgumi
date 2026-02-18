//! Integration tests for fgumi library.
//!
//! These tests validate end-to-end workflows that span multiple modules,
//! ensuring that module interactions work correctly.

mod helpers;
mod test_bam_pipeline;
mod test_clip_command;
mod test_codec_command;
mod test_codec_pipeline;
mod test_correct_command;
mod test_dedup_command;
mod test_downsample_command;
mod test_duplex_command;
mod test_duplex_metrics_command;
mod test_error_paths;
mod test_extract_command;
mod test_fastq_command;
mod test_filter_command;
mod test_group_command;
mod test_pipeline_concurrency;
mod test_review_command;
mod test_simplex_command;
mod test_simplex_pipeline;
mod test_simulate_sort;
mod test_sort_write_index;
mod test_streaming_input;
mod test_zipper_command;
