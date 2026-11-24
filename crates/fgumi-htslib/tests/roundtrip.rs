//! Write CRAM with encodes on worker threads, read it back with slice decodes
//! on worker threads: every record comes back byte-identical (mapped records
//! with exactly MD and NM appended), and the CRAI indexes every slice.
mod common;

use common::{
    appended_md_nm, appended_tags, crai_entries, cram_layout, expected_md_nm, read_cram_slices,
    record, reference, write_cram, write_fasta,
};
use fgumi_htslib::CramWriteOptions;

const REF_LEN: usize = 20_000;

/// A fixture record with its expected regenerated `(MD, NM)` (`None` for an
/// unmapped record), derived from the read and its reference segment.
type Fixture = (Vec<u8>, Option<(String, i64)>);

fn records(seq: &[u8]) -> Vec<Fixture> {
    let mut out = Vec::new();
    for i in 0..25_000usize {
        let pos = (i * 7) % (REF_LEN - 200);
        let mut read = seq[pos..pos + 100].to_vec();
        if i % 5 == 0 {
            read[i % 100] = if read[i % 100] == b'A' { b'C' } else { b'A' };
        }
        let mut tags = b"RXZ".to_vec();
        tags.extend_from_slice(format!("ACGT{:04}", i % 10_000).as_bytes());
        tags.push(0);
        let p = i32::try_from(pos).unwrap();
        let md_nm = expected_md_nm(&read, &seq[pos..pos + 100]);
        out.push((record(&format!("q{i:06}"), 0, 0, p, &[100 << 4], &read, &tags), Some(md_nm)));
    }
    out.sort_by_key(|(r, _)| fgumi_raw_bam::pos(r));
    for i in 0..500 {
        out.push((record(&format!("u{i:04}"), 4, -1, -1, &[], b"ACGTACGTAC", b""), None));
    }
    out
}

#[test]
fn cram_roundtrip_with_external_threads() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0x9E37_79B9_7F4A_7C15);
    let fa = write_fasta(dir.path(), "chr1", &seq);
    let header = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{REF_LEN}\n");
    let out = dir.path().join("out.cram");
    let crai = dir.path().join("out.cram.crai");
    let fixture = records(&seq);
    let recs: Vec<Vec<u8>> = fixture.iter().map(|(r, _)| r.clone()).collect();
    assert!(
        fixture.iter().any(|(_, e)| e.as_ref().is_some_and(|(_, nm)| *nm == 1)),
        "the fixture has reads with a mismatch"
    );

    let opts = CramWriteOptions { reference: Some(&fa), crai: Some(&crai), ..Default::default() };
    let containers = write_cram(&out, &header, &opts, &recs);
    assert!(containers >= 3, "expected several containers, got {containers}");

    let slices = read_cram_slices(&out, Some(&fa)).unwrap();
    let n_slices = slices.len();
    let got: Vec<&Vec<u8>> = slices.iter().flatten().collect();
    assert_eq!(got.len(), recs.len());
    for (i, (g, (e, md_nm))) in got.iter().zip(&fixture).enumerate() {
        match md_nm {
            Some(want) => assert_eq!(&appended_md_nm(g, e), want, "record {i}: MD/NM values"),
            None => assert!(appended_tags(g, e).is_empty(), "record {i}: unmapped, no tags"),
        }
    }

    // The CRAI indexes every slice, in file order, at offsets inside the file.
    let entries = crai_entries(&crai);
    let file_len = std::fs::metadata(&out).unwrap().len();
    assert_eq!(entries.len(), n_slices, "one CRAI entry per slice");
    assert!(entries.windows(2).all(|w| w[0].3 <= w[1].3), "container offsets ascend");
    assert!(entries.iter().all(|e| e.3 > 26 && e.3 < file_len), "offsets inside the file");
    let mapped_starts: Vec<i64> = entries.iter().filter(|e| e.0 == 0).map(|e| e.1).collect();
    assert_eq!(mapped_starts.first(), Some(&1), "first mapped slice starts at position 1");
    assert!(entries.iter().any(|e| e.0 == -1), "the unmapped tail has its own entry");

    // Derive each physical slice range from block types, not from the index
    // or the encoder's landmarks: repeated or misplaced CRAI offsets fail.
    let layout = cram_layout(&std::fs::read(&out).unwrap());
    let mut physical_slices = Vec::new();
    for container in layout {
        let starts: Vec<_> = container
            .blocks
            .iter()
            .filter(|block| block.content_type == 2)
            .map(|block| block.start)
            .collect();
        for (i, &start) in starts.iter().enumerate() {
            let end = starts.get(i + 1).copied().unwrap_or(container.end);
            assert_eq!(
                container
                    .blocks
                    .iter()
                    .rfind(|block| block.start >= start && block.end <= end)
                    .unwrap()
                    .end,
                end,
                "slice ends on a block boundary"
            );
            physical_slices.push((
                container.offset,
                (start - container.payload_offset) as u64,
                (end - start) as u64,
            ));
        }
    }
    let indexed_slices: Vec<_> = entries.iter().map(|entry| (entry.3, entry.4, entry.5)).collect();
    assert_eq!(indexed_slices, physical_slices, "every CRAI entry names its actual slice range");

    // Each row's reference, start and span are those of the records decoded
    // from its slice (single-reference slices: one row per slice, in order),
    // so an indexed query reaches every record of every slice.
    for (i, (entry, slice)) in entries.iter().zip(&slices).enumerate() {
        let want = if fgumi_raw_bam::ref_id(&slice[0]) < 0 {
            (-1, 0, 0)
        } else {
            assert!(slice.iter().all(|r| fgumi_raw_bam::ref_id(r) == 0), "slice {i}: one contig");
            let start = slice.iter().map(|r| i64::from(fgumi_raw_bam::pos(r))).min().unwrap() + 1;
            // Every mapped fixture read spans 100 reference bases.
            let end = slice.iter().map(|r| i64::from(fgumi_raw_bam::pos(r)) + 100).max().unwrap();
            (0, start, end - start + 1)
        };
        assert_eq!((entry.0, entry.1, entry.2), want, "CRAI row {i}: reference, start, span");
    }
}
