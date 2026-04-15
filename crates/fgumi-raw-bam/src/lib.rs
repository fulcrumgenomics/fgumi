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

// Re-exports — callers use fgumi_raw_bam::FooBar etc.
//
// Explicit named list now that fn flags() has been demoted to pub(crate).
// The `flags` entry below re-exports the pub mod (flag-bit constants).
pub use fields::{
    // -- constants --
    BAM_MAGIC,
    MIN_BAM_RECORD_LEN,
    // -- types --
    RawRecordMut,
    RawRecordView,
    TemplateCoordFields,
    TemplateCoordFlags,
    // -- field accessors (free functions) --
    aux_data_offset,
    aux_data_offset_from_record,
    aux_data_slice,
    // -- flag-bit constants module (pub mod, not fn) --
    flags,
    l_seq,
    mapq,
    mate_pos,
    mate_ref_id,
    pos,
    qual_offset,
    read_name,
    ref_id,
    seq_offset,
    // -- field mutators (free functions) --
    set_flags,
    set_mate_pos,
    set_mate_ref_id,
    set_pos,
    set_ref_id,
    set_template_length,
    tag_value_size,
    template_length,
};

// -- builder --
pub use builder::UnmappedBamRecordBuilder;

// -- cigar --
pub use cigar::{
    alignment_end_from_raw, alignment_start_from_raw, cigar_op_kind, cigar_to_string_from_raw,
    clip_cigar_ops_raw, consumes_query, consumes_read, consumes_ref, get_cigar_ops,
    mate_unclipped_5prime, mate_unclipped_5prime_1based, query_length_from_cigar,
    read_pos_at_ref_pos_raw, reference_length_from_cigar, reference_length_from_raw_bam,
    unclipped_5prime_from_raw_bam, unclipped_5prime_raw,
};

// -- overlap --
pub use overlap::{is_fr_pair_raw, num_bases_extending_past_mate_raw};

// -- raw_bam_record --
pub use raw_bam_record::{RawBamReader, RawRecord, read_raw_record};

// -- sequence --
pub use sequence::{
    BAM_BASE_TO_ASCII, extract_sequence, get_base, get_qual, is_base_n, mask_base,
    pack_sequence_into, quality_scores_slice, quality_scores_slice_mut, set_base, set_qual,
};

// -- sort --
pub use sort::{compare_coordinate_raw, compare_names_raw, compare_queryname_raw, natural_compare};

// -- tags --
pub use tags::{
    ArrayTagRef, AuxStringTags, AuxTagsIter, RawTagsEditor, RawTagsMut, RawTagsView, TagEntry,
    TemplateAuxTags, append_i32_array_tag, append_signed_int_tag, array_tag_element_u16,
    array_tag_to_vec_u16, extract_aux_string_tags, extract_template_aux_tags, find_array_tag,
    find_float_tag, find_int_tag, find_mc_tag_in_record, find_string_tag,
    find_string_tag_in_record, find_tag_type, find_uint8_tag, normalize_int_tag_to_smallest_signed,
    remove_tag, reverse_array_tag_in_place, reverse_complement_string_tag_in_place,
    reverse_string_tag_in_place, update_int_tag, update_string_tag,
};

#[cfg(any(test, feature = "test-utils"))]
pub use testutil::*;

#[cfg(feature = "noodles")]
pub use noodles_compat::*;
