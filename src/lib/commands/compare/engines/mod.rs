//! Comparison engines for `fgumi compare bams`.
//!
//! An "engine" pairs records from the two input streams (see
//! [`positional`](self::positional), added alongside the positional-alignment work) and
//! decides whether each pair is content-equal via a [`content::ContentPredicate`].
//! Splitting pairing (positional) from equality (content) keeps the pairing logic honest:
//! it can never quietly resync past a mismatch just because two *unrelated* records
//! happen to compare equal only because a lenient predicate ignores an actual divergence.
//!
//! [`header`] is the one engine-adjacent module that isn't about record pairing at all: it
//! compares the two inputs' `noodles::sam::Header` (`@HD`/`@SQ`/`@RG`, normalizing away
//! `@PG`/`@CO`), since a record-level match says nothing about whether the two files declare
//! the same reference dictionary, read groups, or sort order. Every engine below wires its
//! result in as an additional divergence source alongside its own pairing/content checks.

pub mod content;
pub mod header;
pub mod keyjoin;
pub(crate) mod molecule_join;
pub mod positional;
pub mod sort_verify;

use std::fs::File;
use std::path::Path;

use anyhow::Result;
use fgumi_sort::RawBamRecordReader;

/// Append `msg()` to `details` unless it is already at the `max_diffs` cap — lazily, so
/// callers only pay for building the message string when it will actually be kept.
///
/// Shared by every engine's diff-collection loop (`keyjoin`, `metrics`, and, via
/// [`header::fold_header_diffs`], `positional`/`sort_verify`) — the capping/lazy-build
/// discipline is identical everywhere `diff_details` is populated.
pub(crate) fn push_diff(details: &mut Vec<String>, max_diffs: usize, msg: impl FnOnce() -> String) {
    if details.len() < max_diffs {
        details.push(msg());
    }
}

/// Open a raw-byte BAM record reader over `path`, positioned just past the header.
///
/// Shared by every engine that needs to re-open one of its own input paths for
/// sequential raw-record pulling: `keyjoin`'s canonicalized cursors and `sort_verify`'s
/// per-file reader both open + skip-header identically, differing only in what error
/// context they attach to a failure (added by the caller via `.with_context`, if at all).
///
/// # Errors
///
/// Returns an error if `path` cannot be opened, its BAM header cannot be read, or the
/// header cannot be skipped.
pub(crate) fn open_raw_bam_reader<P: AsRef<Path>>(path: P) -> Result<RawBamRecordReader<File>> {
    let mut reader = RawBamRecordReader::new(File::open(path.as_ref())?)?;
    reader.skip_header()?;
    Ok(reader)
}
