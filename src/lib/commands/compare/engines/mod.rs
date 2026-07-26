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
pub(crate) mod molecule_join;
pub mod positional;
pub mod sort_verify;

use std::path::{Path, PathBuf};

use fgumi_sort::OwnedRawBamRecordReader;
use noodles::sam::Header;

/// One input, opened exactly once: its header and its records from the same stream.
///
/// Passed to the single-pass engines instead of a path so the caller's open is the
/// only one. `CompareBams::execute` must read both headers before dispatching — an
/// incompatible pair has to hard-exit ahead of any engine — and if the engine then
/// opened the path again for its records, the command as a whole would require a
/// re-openable input even though each half of it reads the stream once. A FIFO, a
/// process substitution and stdin are all opened once and never again, so that
/// second open is the difference between accepting them and not.
///
/// `path` is carried only to name the input in diagnostics; nothing reopens it.
pub(crate) struct OpenedInput {
    /// Records, positioned just past the header.
    pub(crate) reader: OwnedRawBamRecordReader,
    /// The header parsed off the front of `reader`'s stream.
    pub(crate) header: Header,
    /// What to call this input in messages.
    pub(crate) path: PathBuf,
}

impl OpenedInput {
    /// Open `path` once, taking its header and its records from one stream.
    ///
    /// [`fgumi_sort::open_raw_bam_record_reader_with_header`] parses the header
    /// through a tee and replays the bytes it consumed, which is what lets one
    /// stream serve both halves. Uncompressed SAM is sniffed and normalized there
    /// too, so callers only ever see BAM records.
    ///
    /// This is the only opener the compare engines use — a second, private copy in
    /// each engine is how the two came to disagree about how many opens an input
    /// takes.
    ///
    /// # Errors
    ///
    /// Returns an error if the input cannot be opened, is neither BAM nor SAM, or
    /// its header cannot be parsed.
    pub(crate) fn open(path: &Path) -> anyhow::Result<Self> {
        use anyhow::Context as _;

        let (reader, header) = fgumi_sort::open_raw_bam_record_reader_with_header(path)
            .with_context(|| format!("opening raw BAM reader for {}", path.display()))?;
        Ok(Self { reader, header, path: path.to_path_buf() })
    }
}

/// Append `msg()` to `details` unless it is already at the `max_diffs` cap — lazily, so
/// callers only pay for building the message string when it will actually be kept.
///
/// Shared by every engine's diff-collection loop (`molecule_join`, `metrics`, and, via
/// [`header::fold_header_diffs`], `positional`/`sort_verify`) — the capping/lazy-build
/// discipline is identical everywhere `diff_details` is populated.
pub(crate) fn push_diff(details: &mut Vec<String>, max_diffs: usize, msg: impl FnOnce() -> String) {
    if details.len() < max_diffs {
        details.push(msg());
    }
}
