//! Verification of BAM sort order.
//!
//! Provides [`verify_sort_order`], an iterator-driven check that yields a
//! [`VerifySummary`] reporting total records seen, sort-order violations, and
//! the location of the first violation (if any).

use anyhow::Result;

use fgumi_raw_bam::RawRecord;

/// Summary of a sort-order verification pass: `(total_records, violations, first_violation)`.
///
/// `first_violation` is `Some((record_number, read_name))` if any violation was
/// observed and `None` otherwise.
pub type VerifySummary = (u64, u64, Option<(u64, String)>);

/// Iterates the records read by `records` and verifies that consecutive
/// extracted keys do not violate the sort invariant.
///
/// `extract_key` produces a sort key from each record's raw bytes. `is_violation`
/// returns `true` when the *current* key violates ordering relative to the
/// *previous* key (typically `|cur, prev| cur < prev`).
///
/// Takes any stream of `io::Result<RawRecord>` rather than a concrete
/// [`RawBamRecordReader`](crate::reader::RawBamRecordReader) so callers can supply
/// a multithreaded or read-ahead record source. The `Result` item type is
/// deliberate: a read failure must stay distinguishable from end-of-stream, since
/// treating a truncated stream as a clean EOF would report a corrupt file as
/// correctly sorted.
///
/// # Errors
///
/// Returns any I/O error from the underlying record stream.
pub fn verify_sort_order<I, K>(
    records: I,
    extract_key: impl Fn(&[u8]) -> K,
    is_violation: impl Fn(&K, &K) -> bool,
) -> Result<VerifySummary>
where
    I: IntoIterator<Item = std::io::Result<RawRecord>>,
{
    let mut total_records: u64 = 0;
    let mut violations: u64 = 0;
    let mut first_violation: Option<(u64, String)> = None;
    let mut prev_key: Option<K> = None;

    for result in records {
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
