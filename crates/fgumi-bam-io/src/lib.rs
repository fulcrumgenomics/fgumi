//! BAM I/O factories and pipeline primitives for the fgumi family of crates.
//!
//! This crate provides the layer that bolts [`noodles`]' BAM record types onto fgumi's
//! faster libdeflate-backed BGZF I/O, plus a small set of pipeline primitives
//! (`ProgressTracker`, `ReorderBuffer`, `MemoryEstimate`) used to consume BAM streams
//! efficiently in parallel.
//!
//! Use this crate when you need to read or write BAM files at high throughput in a
//! pipeline that you want to keep independent of the rest of the fgumi binary.

#![deny(missing_docs)]
#![deny(unsafe_code)]

pub mod grouping;
pub mod header;
pub mod library;
pub mod mem_estimate;
pub mod os_hints;
pub mod paths;
pub mod prefetch_reader;
pub mod progress;
pub mod reader;
pub mod reorder;
pub mod writer;

pub(crate) mod vendored;

pub use grouping::{
    DecodedRecord, GroupKey, GroupKeyConfig, Grouper, compute_group_key_from_raw, name_hash_key,
};
pub use library::{LibraryIndex, LibraryLookup, build_library_lookup};
pub use mem_estimate::MemoryEstimate;
pub use paths::{is_stdin_path, is_stdout_path};
pub use progress::ProgressTracker;
pub use reader::{
    BamReaderAuto, BgzfReaderEnum, ChainedReader, PipelineReaderOpts, RawBamReaderAuto, TeeReader,
    create_bam_reader, create_bam_reader_for_pipeline, create_bam_reader_for_pipeline_with_opts,
    create_bam_reader_with_opts, create_raw_bam_reader, create_raw_bam_reader_with_opts,
};
pub use reorder::{DrainReady, ReorderBuffer};
pub use writer::{
    BamWriter, BgzfWriterEnum, IndexingBamWriter, RawBamWriter, create_bam_writer,
    create_indexing_bam_writer, create_optional_bam_writer, create_raw_bam_writer, write_bai_index,
    write_bam_header,
};
