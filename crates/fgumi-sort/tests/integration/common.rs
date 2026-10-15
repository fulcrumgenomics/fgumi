//! Shared fixtures for fgumi-sort's integration tests: a local stable sort
//! (fgumi-sort cannot depend on fgumi-pipeline-io's `testing` oracle) and a
//! helper that sorts a BAM in memory with it.

use std::path::Path;

use anyhow::Result;
use fgumi_sort::RawBamRecordReader;

/// Stable sort by `key`: the local oracle for fgumi-sort's own tests.
pub fn local_reference_sort<K: Ord>(records: &[Vec<u8>], key: impl Fn(&[u8]) -> K) -> Vec<Vec<u8>> {
    let mut idx: Vec<usize> = (0..records.len()).collect();
    idx.sort_by_key(|&i| key(&records[i])); // `sort_by_key` is stable
    idx.into_iter().map(|i| records[i].clone()).collect()
}

/// Sort the BAM at `input` into `output` with [`local_reference_sort`] by
/// `key`, keeping the input's header.
pub fn sort_bam_in_memory<K: Ord>(
    input: &Path,
    output: &Path,
    key: impl Fn(&[u8]) -> K,
) -> Result<()> {
    let (_, header) = fgumi_bam_io::create_bam_reader(input, 1)?;
    let mut reader = RawBamRecordReader::new(std::fs::File::open(input)?)?;
    reader.skip_header()?;
    let mut records = Vec::new();
    while let Some(record) = reader.next_record()? {
        records.push(record.as_ref().to_vec());
    }
    let mut writer = fgumi_bam_io::create_raw_bam_writer(output, &header, 1, 6)?;
    for record in local_reference_sort(&records, key) {
        writer.write_raw_record(&record)?;
    }
    writer.finish()?;
    Ok(())
}

/// Catches: an unstable sort slipping in (ties reordered), which would make the
/// oracle agree with a broken merge.
#[test]
fn local_reference_sort_is_stable() {
    // 200 records, keys 0..10 (20 per key); the payload is the input index.
    let records: Vec<Vec<u8>> =
        (0..200u32).map(|i| [(i % 10).to_le_bytes(), i.to_le_bytes()].concat()).collect();
    let word = |r: &[u8], at: usize| u32::from_le_bytes(r[at..at + 4].try_into().unwrap());
    let out = local_reference_sort(&records, |r| word(r, 0));
    assert_eq!(out.len(), records.len());
    for w in out.windows(2) {
        let (ka, kb) = (word(&w[0], 0), word(&w[1], 0));
        let (ia, ib) = (word(&w[0], 4), word(&w[1], 4));
        assert!(ka < kb || (ka == kb && ia < ib), "unstable at {ka}/{ia} -> {kb}/{ib}");
    }
}
