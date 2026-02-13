//! SAM/BAM file utilities and header manipulation.
//!
//! This module re-exports functionality from the `fgumi-sam` crate, including
//! alignment tag regeneration via the `alignment_tags` submodule.

// Re-export alignment_tags from fgumi-sam
pub use fgumi_sam::alignment_tags;

// Re-export modules from fgumi-sam
pub use fgumi_sam::{builder, record_utils};

// Re-export top-level functions from fgumi-sam
pub use fgumi_sam::{
    buf_value_to_smallest_signed_int, check_sort, is_sorted, is_template_coordinate_sorted,
    revcomp_buf_value, reverse_buf_value, to_smallest_signed_int,
};

// Re-export commonly used items from submodules
pub use alignment_tags::regenerate_alignment_tags;
pub use fgumi_sam::builder::{
    ConsensusTagsBuilder, FragBuilder, MAPPED_PG_ID, PairBuilder, REFERENCE_LENGTH, RecordBuilder,
    SamBuilder, Strand, create_default_test_fasta, create_ref_dict, create_test_fasta,
    degrading_qualities, parse_cigar, repeat_n, uniform_qualities,
};
pub use fgumi_sam::record_utils::{
    PairOrientation, alignment_end, cigar_reference_length, get_pair_orientation, is_fr_pair,
    is_fr_pair_from_tags, leading_clipping, leading_soft_clipping, mate_unclipped_end,
    mate_unclipped_start, num_bases_extending_past_mate, parse_cigar_string, read_pos_at_ref_pos,
    reference_length, trailing_clipping, trailing_soft_clipping, unclipped_end,
    unclipped_five_prime_position, unclipped_start, unclipped_three_prime_position,
};
