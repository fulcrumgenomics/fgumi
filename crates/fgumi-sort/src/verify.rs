//! Verification of BAM sort order.
//!
//! Provides [`verify_sort_order`], an iterator-driven check that yields a
//! [`VerifySummary`] reporting total records seen, sort-order violations, and
//! the location of the first violation (if any).

use anyhow::Result;
use std::fs::File;

use crate::reader::RawBamRecordReader;

/// Summary of a sort-order verification pass: `(total_records, violations, first_violation)`.
///
/// `first_violation` is `Some((record_number, read_name))` if any violation was
/// observed and `None` otherwise.
pub type VerifySummary = (u64, u64, Option<(u64, String)>);

/// Iterates the records read by `raw_reader` and verifies that consecutive
/// extracted keys do not violate the sort invariant.
///
/// `extract_key` produces a sort key from each record's raw bytes. `is_violation`
/// returns `true` when the *current* key violates ordering relative to the
/// *previous* key (typically `|cur, prev| cur < prev`).
///
/// # Errors
///
/// Returns any I/O error from the underlying record stream.
pub fn verify_sort_order<K>(
    raw_reader: RawBamRecordReader<File>,
    extract_key: impl Fn(&[u8]) -> K,
    is_violation: impl Fn(&K, &K) -> bool,
) -> Result<VerifySummary> {
    let mut total_records: u64 = 0;
    let mut violations: u64 = 0;
    let mut first_violation: Option<(u64, String)> = None;
    let mut prev_key: Option<K> = None;

    for result in raw_reader {
        let record_bytes = result?;
        total_records += 1;
        let bam: &[u8] = &record_bytes;
        let key = extract_key(bam);

        if let Some(ref prev) = prev_key
            && is_violation(&key, prev)
        {
            violations += 1;
            if first_violation.is_none() {
                let name =
                    String::from_utf8_lossy(fgumi_raw_bam::RawRecordView::new(bam).read_name())
                        .to_string();
                first_violation = Some((total_records, name));
            }
        }
        prev_key = Some(key);
    }

    Ok((total_records, violations, first_violation))
}
