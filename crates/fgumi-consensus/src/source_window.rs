//! The realign window (`rw`) a consensus read takes from its source reads.
//!
//! `fgumi group --realign-window` stamps each mapped raw read with its molecule's
//! window; the consensus callers carry it onto the consensus read. Every caller
//! validates the `rw` values it reads the same way: a value that does not parse
//! fails the family with an error naming the read, so a malformed or foreign
//! `rw` is never copied onto a consensus read.

use anyhow::{Context, Result};
use fgumi_raw_bam::SamTag;
use fgumi_sam::realign_window::RealignWindow;

/// Parses the `rw` value of `raw`, if it has one.
///
/// # Errors
///
/// Returns an error naming the read if its `rw` does not parse.
pub(crate) fn parse_source_window(raw: &[u8]) -> Result<Option<(RealignWindow<'_>, &[u8])>> {
    let Some(value) =
        fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(raw), SamTag::RW)
    else {
        return Ok(None);
    };
    let window = RealignWindow::parse(value).with_context(|| {
        format!(
            "invalid rw tag on read '{}'",
            String::from_utf8_lossy(fgumi_raw_bam::read_name(raw))
        )
    })?;
    Ok(Some((window, value)))
}

/// The `rw` value of the first of `raws` that carries one, after checking that
/// every `rw` among them parses.
///
/// # Errors
///
/// Returns an error naming the read if any of `raws` carries an `rw` that does
/// not parse.
pub(crate) fn first_source_window<'a>(
    raws: impl IntoIterator<Item = &'a [u8]>,
) -> Result<Option<&'a [u8]>> {
    let mut first = None;
    for raw in raws {
        if let Some((_, value)) = parse_source_window(raw)? {
            first.get_or_insert(value);
        }
    }
    Ok(first)
}
