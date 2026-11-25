//! Detects half-mapped templates that have been separated from their molecule.
//!
//! `fgumi group` writes a half-mapped template within its molecule, and `fgumi sort` keeps it there
//! by its `tc` tag. Other tools sort it by its own position, after the molecule, and a reader of
//! consecutive records with one `MI` would then take it for a molecule of its own.

use std::io;

use fgumi_raw_bam::{RawRecord, RawRecordView, TemplateCoordinateTag};
use fgumi_sort::{LibraryLookup, cb_hasher, extract_template_key_inline};
use noodles::sam::Header;

/// Checks a molecule whose `records` are all of half-mapped templates that share it with a fully
/// mapped pair, and come after `previous`, the last record of the molecule before.
///
/// Such a molecule is valid where its template coordinate is in order: the input is then a subset
/// that its pairs are not part of.
///
/// # Errors
///
/// Returns an error if the molecule's template coordinate and `MI` come before those of `previous`.
pub(crate) fn check_copies(
    mi: &str,
    records: &[RawRecord],
    previous: Option<&RawRecord>,
) -> io::Result<()> {
    let is_copy =
        |record: &RawRecord| record.is_unmapped() || RawRecordView::new(record).is_half_mapped();
    if !records.iter().all(is_copy) {
        return Ok(());
    }
    let copy = records.iter().find(|r| TemplateCoordinateTag::from_half_mapped(r).is_some());
    let (Some(copy), Some(previous)) = (copy, previous) else { return Ok(()) };
    if position(copy) >= position(previous) {
        return Ok(());
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "the half-mapped reads of MI {} are separate from its fully mapped pairs, because \
             the input is not in `fgumi sort` template-coordinate order. Run `fgumi sort --order \
             template-coordinate` on it",
            mi.split('\t').next().unwrap_or(mi),
        ),
    ))
}

/// The part of the sort key that does not depend on the library or the cell barcode.
fn position(record: &RawRecord) -> (u64, u64, u64) {
    let libraries = LibraryLookup::from_header(&Header::default());
    let key = extract_template_key_inline(record, &libraries, None, &cb_hasher());
    (key.primary, key.secondary, key.tertiary)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::sam::SamTag;
    use fgumi_raw_bam::SamBuilder;
    use rstest::rstest;

    /// R1 of a pair with ends at `pos` and 1000 bases on. If `half`, only R1 is mapped, and it is in
    /// the molecule of such a pair.
    pub(crate) fn read(mi: &str, pos: i32, half: bool) -> RawRecord {
        let mut record = SamBuilder::new();
        record.read_name(b"read").ref_id(0).pos(pos).mate_ref_id(0);
        record.cigar_ops(&[10 << 4]).sequence(&[b'A'; 10]).qualities(&[30; 10]);
        record.add_string_tag(SamTag::MI, mi.as_bytes());
        if half {
            record.flags(0x49).mate_pos(pos);
            record.add_array_i32(SamTag::TC, &[0, pos + 1, 0, 0, pos + 1001, 0]);
        } else {
            record.flags(0x41).mate_pos(pos + 1000).add_string_tag(SamTag::MC, b"10M");
        }
        record.build()
    }

    #[rstest]
    #[case::first(&[(5, 100, true)], None, true)]
    #[case::after_a_lower_position(&[(5, 100, true)], Some((7, 50, false)), true)]
    #[case::after_a_higher_position(&[(5, 100, true)], Some((4, 150, false)), false)]
    #[case::after_a_lower_molecule(&[(5, 100, true)], Some((4, 100, false)), true)]
    #[case::after_a_higher_molecule(&[(5, 100, true)], Some((6, 100, false)), false)]
    #[case::after_higher_copies(&[(5, 100, true)], Some((6, 150, true)), false)]
    #[case::several(&[(5, 100, true), (5, 100, true)], Some((6, 150, false)), false)]
    #[case::with_a_pair(&[(5, 100, true), (5, 100, false)], Some((6, 150, false)), true)]
    #[case::pair(&[(5, 100, false)], Some((6, 150, false)), true)]
    fn copies_have_to_be_in_order(
        #[case] records: &[(u32, i32, bool)],
        #[case] previous: Option<(u32, i32, bool)>,
        #[case] expected: bool,
    ) {
        let read = |(mi, pos, half): (u32, i32, bool)| read(&mi.to_string(), pos, half);
        let records: Vec<_> = records.iter().copied().map(read).collect();
        let found = check_copies("5\tACGT", &records, previous.map(read).as_ref());
        assert_eq!(found.is_ok(), expected);
        if let Err(e) = found {
            assert!(e.to_string().contains("of MI 5 are separate"), "{e}");
        }
    }
}
