//! Puts the half-mapped templates that share a molecule with a fully mapped pair into that
//! molecule's place in `fgumi group`'s output.
//!
//! A half-mapped template is reached after its molecule has been passed on, so it is kept in a
//! [`CopyStore`] until the end of the input. Molecule ids ascend along the output, so the id alone
//! says which batch of position groups it belongs to; [`insert`] adds it to that batch's records.

use std::collections::BTreeMap;
use std::io::{self, Read, Seek, SeekFrom, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use fgumi_raw_bam::{
    RawRecord, RawRecordView, TemplateCoordinateTag, unclipped_5prime_from_raw_bam,
    write_framed_record,
};
use parking_lot::Mutex;
use tempfile::SpooledTempFile;

use crate::sam::SamTag;

pub(crate) mod hold;
pub(crate) mod splice;

/// Bytes of half-mapped templates kept in memory; the rest go to a temporary file.
const MAX_BYTES_IN_MEMORY: usize = 64 << 20;

/// The half-mapped templates to be put into their molecules.
#[derive(Debug)]
pub(crate) struct CopyStore {
    copies: Mutex<Copies>,
    any: Arc<AtomicBool>,
    size_changes: Mutex<Vec<FamilySizeChange>>,
}

#[derive(Debug)]
struct Copies {
    /// The first molecule id after each batch of position groups.
    batch_ends: Vec<u64>,
    /// Where each template's records are in `records`.
    index: Vec<CopyRef>,
    /// Records, each after its length.
    records: SpooledTempFile,
    records_len: u64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CopyRef {
    molecule: u64,
    offset: u64,
    len: usize,
}

/// A family that went from one number of templates to another; there are none of size zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct FamilySizeChange {
    pub(crate) from: usize,
    pub(crate) to: usize,
}

/// The records of the half-mapped templates that go into one batch of position groups.
pub(crate) type BatchCopies = Vec<RawRecord>;

impl Default for CopyStore {
    fn default() -> Self {
        let copies = Copies {
            batch_ends: Vec::new(),
            index: Vec::new(),
            records: SpooledTempFile::new(MAX_BYTES_IN_MEMORY),
            records_len: 0,
        };
        Self { copies: Mutex::new(copies), any: Arc::default(), size_changes: Mutex::default() }
    }
}

impl CopyStore {
    /// Adds the framed `records` of a half-mapped template, which carry its `molecule`.
    pub(crate) fn add(&self, molecule: u64, records: &[u8]) -> io::Result<()> {
        let mut copies = self.copies.lock();
        copies.records.write_all(records)?;
        let copy = CopyRef { molecule, offset: copies.records_len, len: records.len() };
        copies.records_len += records.len() as u64;
        copies.index.push(copy);
        self.any.store(true, Ordering::Release);
        Ok(())
    }

    /// Ends a batch of position groups, whose molecule ids are all below `next_molecule`.
    pub(crate) fn end_batch(&self, next_molecule: u64) {
        self.copies.lock().batch_ends.push(next_molecule);
    }

    /// How many half-mapped templates there are.
    pub(crate) fn len(&self) -> usize {
        self.copies.lock().index.len()
    }

    /// Set once there is one.
    pub(crate) fn any(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.any)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The copies of each batch that has any, by the batch's index.
    pub(crate) fn by_batch(&self) -> Vec<(usize, Vec<CopyRef>)> {
        let mut copies = self.copies.lock();
        copies.index.sort_by_key(|copy| copy.molecule);
        let mut batches: Vec<(usize, Vec<CopyRef>)> = Vec::new();
        for &copy in &copies.index {
            let batch = copies.batch_ends.partition_point(|&end| end <= copy.molecule);
            match batches.last_mut() {
                Some((last, refs)) if *last == batch => refs.push(copy),
                _ => batches.push((batch, vec![copy])),
            }
        }
        batches
    }

    pub(crate) fn read(&self, refs: &[CopyRef]) -> io::Result<BatchCopies> {
        let mut copies = self.copies.lock();
        let mut bytes = Vec::new();
        for copy in refs {
            copies.records.seek(SeekFrom::Start(copy.offset))?;
            let start = bytes.len();
            bytes.resize(start + copy.len, 0);
            copies.records.read_exact(&mut bytes[start..])?;
        }
        split_records(&bytes)
    }

    pub(crate) fn add_size_changes(&self, changes: &[FamilySizeChange]) {
        self.size_changes.lock().extend_from_slice(changes);
    }

    /// The family sizes that changed as templates were put into their molecules.
    pub(crate) fn take_size_changes(&self) -> Vec<FamilySizeChange> {
        std::mem::take(&mut self.size_changes.lock())
    }
}

/// The bytes of a record's length, which comes before it.
const FRAME_LEN: usize = 4;

/// The framed records in `bytes`.
fn framed_records(mut bytes: &[u8]) -> io::Result<Vec<&[u8]>> {
    let truncated = || io::Error::new(io::ErrorKind::UnexpectedEof, "a record is truncated");
    let mut records = Vec::new();
    while !bytes.is_empty() {
        let (len, rest) = bytes.split_first_chunk::<FRAME_LEN>().ok_or_else(truncated)?;
        let len = u32::from_le_bytes(*len) as usize;
        records.push(rest.get(..len).ok_or_else(truncated)?);
        bytes = &rest[len..];
    }
    Ok(records)
}

/// Splits framed records.
pub(crate) fn split_records(bytes: &[u8]) -> io::Result<Vec<RawRecord>> {
    let records = framed_records(bytes)?;
    Ok(records.into_iter().map(|record| RawRecord::from(record.to_vec())).collect())
}

/// The molecule of a record, and whether it is of its second strand.
fn family(record: &[u8]) -> Option<(u64, bool)> {
    let aux = fgumi_raw_bam::aux_data_slice(record);
    let mi = fgumi_raw_bam::find_string_tag(aux, SamTag::MI)?;
    let (id, strand) = match mi.iter().position(|&b| b == b'/') {
        Some(slash) => (&mi[..slash], &mi[slash + 1..]),
        None => (mi, &b""[..]),
    };
    Some((std::str::from_utf8(id).ok()?.parse().ok()?, strand == b"B"))
}

fn name(record: &[u8]) -> &[u8] {
    RawRecordView::new(record).read_name()
}

/// What the copies of a molecule need from its other templates.
#[derive(Debug, Default)]
struct Molecule {
    coordinate: Option<TemplateCoordinateTag>,
    /// The templates of each strand.
    templates: [usize; 2],
    copies: [usize; 2],
    /// The `rw` tag of each mapped end, if any.
    windows: Vec<(End, Vec<u8>)>,
}

type End = (i32, i32, bool);

fn end(record: &[u8]) -> Option<End> {
    let (view, pos) = (RawRecordView::new(record), unclipped_5prime_from_raw_bam(record));
    (!view.is_unmapped()).then(|| (view.ref_id(), pos, view.is_reverse()))
}

/// Sorts `copies` into the order `fgumi group` writes templates in, and gives each its molecule's
/// template coordinate (`tc`), which is where `fgumi sort` puts it. Returns the index in `records`
/// that each copy comes before, and the families that become larger.
///
/// The `records` are in that order and hold the molecule of each copy.
fn place<R: AsRef<[u8]>>(
    records: &[R],
    copies: &mut BatchCopies,
) -> (Vec<usize>, Vec<FamilySizeChange>) {
    let mut molecules: BTreeMap<u64, Molecule> = BTreeMap::new();
    for copy in copies.chunk_by(|a, b| name(a) == name(b)) {
        if let Some((id, strand)) = family(&copy[0]) {
            molecules.entry(id).or_default().copies[usize::from(strand)] += 1;
        }
    }
    let families: Vec<_> = records.iter().map(|record| family(record.as_ref())).collect();
    let mut start = 0;
    for template in records.chunk_by(|a, b| name(a.as_ref()) == name(b.as_ref())) {
        if let Some((id, strand)) = families[start]
            && let Some(molecule) = molecules.get_mut(&id)
        {
            molecule.add_template(template, strand);
        }
        start += template.len();
    }

    for copy in copies.iter_mut() {
        if let Some(molecule) = family(copy).and_then(|(id, _)| molecules.get(&id)) {
            molecule.stamp(copy);
        }
    }
    copies.sort_by(|a, b| (family(a), name(a)).cmp(&(family(b), name(b))));

    let mut places = Vec::with_capacity(copies.len());
    let mut unplaced = copies.iter().peekable();
    let mut last_family = None;
    for (index, record_family) in families.into_iter().enumerate() {
        let precedes = |copy: &&RawRecord| match record_family {
            Some(known) => {
                (family(copy), name(copy)) < (Some(known), name(records[index].as_ref()))
            }
            None => family(copy) <= last_family,
        };
        while unplaced.next_if(precedes).is_some() {
            places.push(index);
        }
        last_family = record_family.or(last_family);
    }
    places.resize(copies.len(), records.len());

    (places, molecules.values().flat_map(Molecule::size_changes).collect())
}

/// Puts `copies` into `records`, as [`place`] says. Returns the families that become larger.
pub(crate) fn insert(
    records: Vec<RawRecord>,
    mut copies: BatchCopies,
) -> (Vec<RawRecord>, Vec<FamilySizeChange>) {
    let (places, changes) = place(&records, &mut copies);
    let mut merged = Vec::with_capacity(records.len() + copies.len());
    let mut copies = places.into_iter().zip(copies).peekable();
    for (index, record) in records.into_iter().enumerate() {
        while let Some((_, copy)) = copies.next_if(|(place, _)| *place == index) {
            merged.push(copy);
        }
        merged.push(record);
    }
    merged.extend(copies.map(|(_, copy)| copy));
    (merged, changes)
}

/// Puts `copies` into framed records, as [`insert`] does into others.
pub(crate) fn insert_framed(
    bytes: &[u8],
    mut copies: BatchCopies,
) -> io::Result<(Vec<u8>, Vec<FamilySizeChange>)> {
    let records = framed_records(bytes)?;
    let (places, changes) = place(&records, &mut copies);
    let copies_len: usize = copies.iter().map(|copy| copy.len() + FRAME_LEN).sum();
    let mut merged = Vec::with_capacity(bytes.len() + copies_len);
    let (mut record, mut offset) = (0, 0);
    for (place, copy) in places.into_iter().zip(&copies) {
        let skipped: usize = records[record..place].iter().map(|r| r.len() + FRAME_LEN).sum();
        merged.extend_from_slice(&bytes[offset..offset + skipped]);
        (record, offset) = (place, offset + skipped);
        write_framed_record(&mut merged, copy)?;
    }
    merged.extend_from_slice(&bytes[offset..]);
    Ok((merged, changes))
}

impl Molecule {
    fn add_template<R: AsRef<[u8]>>(&mut self, template: &[R], strand: bool) {
        self.templates[usize::from(strand)] += 1;
        let template = template.iter().map(AsRef::as_ref);
        let mut ends = template.clone().filter_map(end);
        if let (None, Some(a), Some(b)) = (&self.coordinate, ends.next(), ends.next()) {
            self.coordinate = Some(TemplateCoordinateTag::from_primaries(a, b));
        }
        for record in template {
            let aux = fgumi_raw_bam::aux_data_slice(record);
            if let (Some(end), Some(window)) =
                (end(record), fgumi_raw_bam::find_string_tag(aux, SamTag::RW))
                && !self.windows.iter().any(|(known, _)| *known == end)
            {
                self.windows.push((end, window.to_vec()));
            }
        }
    }

    fn stamp(&self, copy: &mut RawRecord) {
        if let Some(coordinate) = &self.coordinate {
            let record = copy.as_mut_vec();
            fgumi_raw_bam::remove_tag(record, SamTag::TC);
            fgumi_raw_bam::append_i32_array_tag(record, SamTag::TC, &coordinate.to_values());
        }
        let window = end(copy).and_then(|end| self.windows.iter().find(|(known, _)| *known == end));
        if let Some((_, window)) = window {
            fgumi_raw_bam::update_string_tag(copy.as_mut_vec(), SamTag::RW, window);
        }
    }

    fn size_changes(&self) -> impl Iterator<Item = FamilySizeChange> {
        let changes = self.templates.into_iter().zip(self.copies);
        changes
            .filter(|(_, copies)| *copies > 0)
            .map(|(templates, copies)| FamilySizeChange { from: templates, to: templates + copies })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder;
    use rstest::rstest;

    /// The two reads of a template, with only R1 mapped if `half`.
    pub(super) fn template(name: &str, mi: Option<&str>, half: bool) -> Vec<RawRecord> {
        let unmapped = |bit| if half { bit } else { 0 };
        [(0x41 | unmapped(0x8), 100), (0x81 | unmapped(0x4), 300)]
            .into_iter()
            .map(|(flags, pos)| {
                let mut record = SamBuilder::new();
                record.read_name(name.as_bytes()).flags(flags).ref_id(0).pos(pos);
                record.cigar_ops(&[10 << 4]).sequence(&[b'A'; 10]).qualities(&[30; 10]);
                if let Some(mi) = mi {
                    record.add_string_tag(SamTag::MI, mi.as_bytes());
                }
                record.build()
            })
            .collect()
    }

    fn templates(of: &[(&str, Option<&str>)], half: bool) -> Vec<RawRecord> {
        of.iter().flat_map(|(name, mi)| template(name, *mi, half)).collect()
    }

    fn names(records: &[RawRecord]) -> Vec<String> {
        records
            .chunks(2)
            .map(|template| String::from_utf8(name(&template[0]).to_vec()).unwrap())
            .collect()
    }

    /// Copies go where `fgumi group` would have written them: by molecule, strand, then name.
    #[rstest]
    #[case::first(&[("a", "3/A")], &["a", "b", "d", "f", "x", "h"])]
    #[case::between_names(&[("c", "3/A")], &["b", "c", "d", "f", "x", "h"])]
    #[case::last_of_a_strand(&[("e", "3/A")], &["b", "d", "e", "f", "x", "h"])]
    #[case::first_of_a_strand(&[("a", "3/B")], &["b", "d", "a", "f", "x", "h"])]
    #[case::before_an_unassigned_template(&[("z", "3/B")], &["b", "d", "f", "z", "x", "h"])]
    #[case::new_strand(&[("a", "4/B")], &["b", "d", "f", "x", "h", "a"])]
    #[case::several(
        &[("z", "4/A"), ("a", "4/A"), ("e", "3/A"), ("c", "3/A")],
        &["b", "c", "d", "e", "f", "x", "a", "h", "z"]
    )]
    fn copies_are_put_in_order(#[case] copies: &[(&str, &str)], #[case] expected: &[&str]) {
        let of = [
            ("b", Some("3/A")),
            ("d", Some("3/A")),
            ("f", Some("3/B")),
            ("x", None),
            ("h", Some("4/A")),
        ];
        let copies: Vec<_> = copies.iter().map(|(name, mi)| (*name, Some(*mi))).collect();
        let (merged, changes) = insert(templates(&of, false), templates(&copies, true));
        assert_eq!(names(&merged), expected);

        let mut bytes = Vec::new();
        for record in templates(&of, false) {
            write_framed_record(&mut bytes, &record).unwrap();
        }
        let framed = insert_framed(&bytes, templates(&copies, true)).unwrap();
        assert_eq!(split_records(&framed.0).unwrap(), merged);
        assert_eq!(framed.1, changes);

        for record in &merged {
            let tc = TemplateCoordinateTag::from_aux_data(fgumi_raw_bam::aux_data_slice(record));
            let is_copy = copies.iter().any(|(copy, _)| copy.as_bytes() == name(record));
            assert_eq!(tc.map(|tc| tc.to_values()), is_copy.then_some([0, 101, 0, 0, 301, 0]));
        }
    }

    #[test]
    fn size_changes_are_by_strand() {
        let of = [("a", Some("3/A")), ("b", Some("3/A")), ("c", Some("4/A"))];
        let copies = [("d", Some("3/A")), ("e", Some("3/A")), ("f", Some("4/B"))];
        let (_, mut changes) = insert(templates(&of, false), templates(&copies, true));
        changes.sort_unstable();
        let change = |from, to| FamilySizeChange { from, to };
        assert_eq!(changes, [change(0, 1), change(2, 4)]);
    }

    /// Batches hold the molecules from the end of the batch before up to their own end.
    #[test]
    fn copies_are_found_by_batch() {
        let store = CopyStore::default();
        let any = store.any();
        assert!(store.is_empty() && !any.load(Ordering::Acquire));
        for end in [2, 2, 5, 6] {
            store.end_batch(end);
        }
        for molecule in [5, 1, 4, 2, 0] {
            let mut bytes = Vec::new();
            write_framed_record(&mut bytes, &template(&molecule.to_string(), None, true)[0])
                .unwrap();
            store.add(molecule, &bytes).unwrap();
        }
        assert!(store.len() == 5 && any.load(Ordering::Acquire));
        let found: Vec<_> = store
            .by_batch()
            .into_iter()
            .map(|(batch, refs)| {
                let records = store.read(&refs).unwrap();
                (batch, records.iter().map(|r| name(r).to_vec()).collect::<Vec<_>>())
            })
            .collect();
        let expected: Vec<(usize, Vec<Vec<u8>>)> = vec![
            (0, vec![b"0".to_vec(), b"1".to_vec()]),
            (2, vec![b"2".to_vec(), b"4".to_vec()]),
            (3, vec![b"5".to_vec()]),
        ];
        assert_eq!(found, expected);
    }

    #[rstest]
    #[case::length(2)]
    #[case::record(10)]
    fn truncated_records_are_an_error(#[case] len: usize) {
        let mut bytes = Vec::new();
        write_framed_record(&mut bytes, &template("a", None, true)[0]).unwrap();
        assert_eq!(split_records(&bytes).unwrap().len(), 1);
        assert!(split_records(&bytes[..len]).is_err());
    }
}
