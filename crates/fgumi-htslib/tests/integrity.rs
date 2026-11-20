//! A CRAM that is cut short, or whose write was abandoned, must never read as
//! complete: htslib alone would end such a stream cleanly or stamp an EOF
//! container on it.
mod common;

use common::{read_cram, record, reference, write_cram, write_fasta};
use std::io::{self, ErrorKind, Read};

use fgumi_htslib::{CramReader, CramWriteOptions, CramWriter};
use rstest::rstest;

const REF_LEN: usize = 20_000;
const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:20000\n";
/// CRAM 3's EOF container is 38 bytes.
const EOF_LEN: u64 = 38;

fn records(seq: &[u8]) -> Vec<Vec<u8>> {
    let mut v: Vec<Vec<u8>> = (0..30_000usize)
        .map(|i| {
            let pos = (i * 3) % (REF_LEN - 100);
            record(
                &format!("r{i}"),
                0,
                0,
                i32::try_from(pos).unwrap(),
                &[50 << 4],
                &seq[pos..pos + 50],
                b"",
            )
        })
        .collect();
    v.sort_by_key(|r| fgumi_raw_bam::pos(r));
    v
}

/// A cut in the middle of the data: past the header, inside some container.
fn mid_data_offset(path: &std::path::Path) -> u64 {
    (std::fs::metadata(path).unwrap().len() - EOF_LEN) / 2
}

#[rstest]
#[case::missing_eof_container(true, "ends without its EOF container (truncated)")]
#[case::cut_inside_a_container(false, "reading the next CRAM container or slice failed")]
fn truncated_cram_is_an_error(#[case] at_eof_boundary: bool, #[case] why: &str) {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0x2545_F491_4F6C_DD1D);
    let fa = write_fasta(dir.path(), "chrT", &seq);
    let full = dir.path().join("full.cram");
    write_cram(
        &full,
        HEADER,
        &CramWriteOptions { reference: Some(&fa), ..Default::default() },
        &records(&seq),
    );
    let bytes = std::fs::read(&full).unwrap();
    let cut = if at_eof_boundary {
        bytes.len() - usize::try_from(EOF_LEN).unwrap() // every data container, no EOF
    } else {
        usize::try_from(mid_data_offset(&full)).unwrap()
    };
    let truncated = dir.path().join("truncated.cram");
    std::fs::write(&truncated, &bytes[..cut]).unwrap();

    assert!(read_cram(&full, Some(&fa)).is_ok(), "the intact file reads");
    let Err(err) = read_cram(&truncated, Some(&fa)) else {
        panic!("a truncated CRAM must not read as complete");
    };
    assert!(err.to_string().contains(why), "{err}");
}

#[test]
fn abandoned_write_leaves_no_eof_container_and_no_index() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0x2545_F491_4F6C_DD1D);
    let fa = write_fasta(dir.path(), "chrT", &seq);
    let out = dir.path().join("abandoned.cram");
    let crai = dir.path().join("abandoned.cram.crai");
    let opts = CramWriteOptions { reference: Some(&fa), crai: Some(&crai), ..Default::default() };
    let mut w = CramWriter::create(&out, HEADER, &opts).unwrap();
    for r in records(&seq) {
        w.write_raw(&r).unwrap();
    }
    let mut jobs: Vec<_> = w.drain_jobs().collect();
    assert!(jobs.len() >= 2, "need several containers");
    // Write the first container, lose the second (an encode error, say), and
    // abandon everything else: no finisher ever runs.
    jobs.remove(0).encode().unwrap().write().unwrap();
    drop(jobs);
    drop(w);

    assert!(!crai.exists(), "an abandoned write leaves no index");
    let Err(err) = read_cram(&out, Some(&fa)) else {
        panic!("an abandoned CRAM must not read as complete");
    };
    assert!(err.to_string().contains("ends without its EOF container"), "{err}");
}

/// Writes the test records as a reference-based CRAM; returns (dir, fasta, cram).
fn written_cram() -> (tempfile::TempDir, std::path::PathBuf, std::path::PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0x2545_F491_4F6C_DD1D);
    let fa = write_fasta(dir.path(), "chrT", &seq);
    let out = dir.path().join("full.cram");
    write_cram(
        &out,
        HEADER,
        &CramWriteOptions { reference: Some(&fa), ..Default::default() },
        &records(&seq),
    );
    (dir, fa, out)
}

#[rstest]
#[case::one_byte(1)]
#[case::two_bytes(2)]
#[case::three_bytes(3)]
fn stray_bytes_after_the_eof_container_are_an_error(#[case] n_stray: usize) {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(&full).unwrap();
    bytes.extend(std::iter::repeat_n(0x5a, n_stray));
    let padded = dir.path().join("padded.cram");
    std::fs::write(&padded, &bytes).unwrap();

    assert!(read_cram(&full, Some(&fa)).is_ok(), "the intact file reads");
    let Err(err) = read_cram(&padded, Some(&fa)) else {
        panic!("bytes after the EOF container must not read as a clean end");
    };
    assert!(err.to_string().contains("stray bytes after its EOF container"), "{err}");
}

/// A valid reader that fails instead of reporting physical EOF.
struct ErrorAtEnd(io::Cursor<Vec<u8>>);
impl Read for ErrorAtEnd {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.0.position() == self.0.get_ref().len() as u64 {
            Err(io::Error::other("injected read failure after footer"))
        } else {
            self.0.read(buf)
        }
    }
}

#[test]
fn io_error_after_eof_container_is_not_clean_eof() {
    let (_dir, fa, full) = written_cram();
    let bytes = std::fs::read(full).unwrap();
    let mut reader = CramReader::from_reader(
        Box::new(ErrorAtEnd(io::Cursor::new(bytes))),
        "failing-reader",
        Some(&fa),
    )
    .unwrap();
    loop {
        match reader.next_slice() {
            Ok(Some(_)) => {}
            Err(err) => {
                assert!(
                    err.to_string().contains("reading the next CRAM container or slice failed"),
                    "{err}"
                );
                break;
            }
            Ok(None) => panic!("an I/O error after the footer must not read as clean EOF"),
        }
    }
}

#[test]
fn clean_eof_remains_eof_on_repeated_reads() {
    let (_dir, fa, full) = written_cram();
    let mut reader = CramReader::open(&full, Some(&fa)).unwrap();
    while reader.next_slice().unwrap().is_some() {}
    assert!(reader.next_slice().unwrap().is_none());
}

#[test]
fn duplicate_eof_container_is_an_error() {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(full).unwrap();
    bytes.extend_from_within(bytes.len() - usize::try_from(EOF_LEN).unwrap()..);
    let padded = dir.path().join("duplicate-eof.cram");
    std::fs::write(&padded, bytes).unwrap();
    let Err(err) = read_cram(&padded, Some(&fa)) else {
        panic!("a duplicate EOF container must not read as complete");
    };
    assert!(err.to_string().contains("stray bytes after its EOF container"), "{err}");
}

/// CRC32 for a synthetic CRAM 3 container header, independent of htslib.
fn crc32(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0xedb8_8320 & 0u32.wrapping_sub(crc & 1));
        }
    }
    !crc
}

#[rstest]
#[case::zero_length(0, 0)]
#[case::wrong_span(13, 1)]
#[case::wrong_record_counter(15, 1)]
#[case::wrong_num_bases(16, 1)]
#[case::wrong_num_blocks(17, 2)]
fn malformed_eof_container_is_an_error(#[case] field: usize, #[case] value: u8) {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(full).unwrap();
    let start = bytes.len() - usize::try_from(EOF_LEN).unwrap();
    // CRAM 3's footer has a 19-byte header followed by its CRC32.
    bytes[start + field] = value;
    let checksum = crc32(&bytes[start..start + 19]);
    bytes[start + 19..start + 23].copy_from_slice(&checksum.to_le_bytes());
    if field == 0 {
        bytes.truncate(start + 23);
    }
    let malformed = dir.path().join("malformed-eof.cram");
    std::fs::write(&malformed, bytes).unwrap();
    assert!(read_cram(&malformed, Some(&fa)).is_err(), "invalid EOF footer accepted");
}

#[rstest]
#[case::v1_0(1, 0)]
#[case::v2_0(2, 0)]
fn cram_before_2_1_is_refused_at_open(#[case] major: u8, #[case] minor: u8) {
    let dir = tempfile::tempdir().unwrap();
    // The file definition alone (magic, version, 20-byte file ID): the version
    // is checked before htslib parses anything else.
    let mut bytes = b"CRAM".to_vec();
    bytes.extend_from_slice(&[major, minor]);
    bytes.extend_from_slice(&[0; 20]);
    bytes.extend_from_slice(&[0; 64]);
    let old = dir.path().join("old.cram");
    std::fs::write(&old, &bytes).unwrap();

    let Err(err) = CramReader::open(&old, None) else {
        panic!("a CRAM {major}.{minor} input must be refused");
    };
    assert_eq!(err.kind(), ErrorKind::InvalidData);
    let msg = err.to_string();
    assert!(msg.contains(&format!("CRAM version {major}.{minor} is not supported")), "{msg}");
    assert!(msg.contains("samtools view -C"), "{msg}");
}

#[test]
fn out_of_order_container_write_is_refused_and_leaves_the_file_incomplete() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0x2545_F491_4F6C_DD1D);
    let fa = write_fasta(dir.path(), "chrT", &seq);
    let out = dir.path().join("out_of_order.cram");
    let opts = CramWriteOptions { reference: Some(&fa), ..Default::default() };
    let mut w = CramWriter::create(&out, HEADER, &opts).unwrap();
    for r in records(&seq) {
        w.write_raw(&r).unwrap();
    }
    let mut jobs: Vec<_> = w.drain_jobs().collect();
    assert!(jobs.len() >= 2, "need several containers");
    let second = jobs.remove(1).encode().unwrap();
    assert!(second.write().is_err(), "container 1 written before container 0");
    drop(jobs);
    let (rest, fin) = w.seal().unwrap();
    drop(rest);
    assert!(fin.finish().is_err(), "finish with unwritten containers");

    let err = read_cram(&out, Some(&fa)).unwrap_err();
    assert!(err.to_string().contains("ends without its EOF container"), "{err}");
}

#[test]
fn finish_with_unwritten_containers_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference(REF_LEN, 0x2545_F491_4F6C_DD1D);
    let fa = write_fasta(dir.path(), "chrT", &seq);
    let out = dir.path().join("unfinished.cram");
    let opts = CramWriteOptions { reference: Some(&fa), ..Default::default() };
    let mut w = CramWriter::create(&out, HEADER, &opts).unwrap();
    for r in records(&seq) {
        w.write_raw(&r).unwrap();
    }
    let mut jobs: Vec<_> = w.drain_jobs().collect();
    assert!(!jobs.is_empty(), "need a container");
    // Write every container but the final one handed over by `seal`.
    let (tail, fin) = w.seal().unwrap();
    assert!(!tail.is_empty(), "seal hands over the final container");
    jobs.drain(..).for_each(|j| j.encode().unwrap().write().unwrap());
    drop(tail);

    let err = fin.finish().unwrap_err();
    assert!(err.to_string().contains("containers written"), "{err}");
    let mut reader = CramReader::open(&out, Some(&fa)).unwrap();
    let first = reader.next_slice().unwrap().expect("written slices survive abandonment");
    assert!(first.decode().unwrap().records > 0);
    let err = read_cram(&out, Some(&fa)).unwrap_err();
    assert!(err.to_string().contains("ends without its EOF container"), "{err}");
}

#[test]
fn ordinary_mapped_slices_leave_room_for_parallel_decode() {
    let (_dir, fa, full) = written_cram();
    let mut reader = CramReader::open(&full, Some(&fa)).unwrap();
    let job = reader.next_slice().unwrap().unwrap();
    let estimate = job.decoded_size_estimate().unwrap();
    // Two ordinary short-read slices fit a 512 MiB decoded allowance.
    assert!(estimate < 256 * 1024 * 1024, "ordinary slice estimate: {estimate}");
    assert!(estimate >= job.decode().unwrap().bytes.len() as u64 * 3);
}

/// Generate a header-only CRAM 2.1 independently of the writer under test.
fn empty_cram_2_1() -> Vec<u8> {
    let header = b"@HD\tVN:1.6\n";
    let payload_len = u8::try_from(4 + header.len()).unwrap();
    let block_len = 5 + u32::from(payload_len);
    let mut bytes = b"CRAM\x02\x01".to_vec();
    bytes.extend_from_slice(&[0; 20]);
    bytes.extend_from_slice(&block_len.to_le_bytes());
    bytes.extend_from_slice(&[0xff, 0xff, 0xff, 0xff, 0x0f]); // reference -1
    bytes.extend_from_slice(&[0, 0, 0, 0, 0, 1, 0]);
    bytes.extend_from_slice(&[0, 0, 0, payload_len, payload_len]);
    bytes.extend_from_slice(&u32::try_from(header.len()).unwrap().to_le_bytes());
    bytes.extend_from_slice(header);
    bytes.extend_from_slice(&[
        11, 0, 0, 0, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0, 0, 0, 0, 1, 0, 0, 1,
        0, 6, 6, 1, 0, 1, 0, 1, 0,
    ]);
    bytes
}

#[rstest]
#[case::complete(0, true)]
#[case::missing_footer(30, false)]
#[case::missing_payload(11, false)]
#[case::cut_payload(1, false)]
fn cram_2_1_footer_integrity(#[case] cut: usize, #[case] complete: bool) {
    let mut bytes = empty_cram_2_1();
    bytes.truncate(bytes.len() - cut);
    let mut reader =
        CramReader::from_reader(Box::new(io::Cursor::new(bytes)), "generated-2.1", None).unwrap();
    let result = reader.next_slice();
    if complete {
        assert!(result.unwrap().is_none());
    } else {
        assert!(result.is_err(), "truncated CRAM 2.1 must not read as complete");
    }
}

#[test]
fn decoded_size_estimate_covers_long_reads_before_decoding() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("long-reads.cram");
    let seq = reference(20_000, 0x00ab_cdef);
    let recs: Vec<_> =
        (0..32).map(|i| record(&format!("long{i}"), 4, -1, -1, &[], &seq, b"")).collect();
    write_cram(&out, HEADER, &CramWriteOptions { no_ref: true, ..Default::default() }, &recs);
    let mut reader = CramReader::open(&out, None).unwrap();
    let job = reader.next_slice().unwrap().unwrap();
    let estimate = job.decoded_size_estimate().unwrap();
    let decoded = job.decode().unwrap();
    assert!(
        estimate >= decoded.bytes.len() as u64 * 3,
        "native, copied and decoded buffers need admission"
    );
}

/// A raw record with its fixed fields; `l_seq` and `n_cigar` are as given, so
/// the record can disagree with itself.
fn raw(tid: i32, cigar: &[u32], seq_len: usize) -> Vec<u8> {
    let seq = vec![b'A'; seq_len];
    record("bad", 0, tid, 10, cigar, &seq, b"")
}

#[rstest]
#[case::too_short_for_a_bam_record(vec![0u8; 10], "a malformed BAM record")]
#[case::reference_id_outside_the_header(raw(5, &[50 << 4], 50), "a reference id outside the header")]
#[case::cigar_and_sequence_lengths_differ(raw(0, &[50 << 4], 40), "a CIGAR whose query length")]
fn write_raw_rejects_malformed_records(#[case] rec: Vec<u8>, #[case] why: &str) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("rejects.cram");
    let opts = CramWriteOptions { no_ref: true, ..Default::default() };
    let mut w = CramWriter::create(&out, HEADER, &opts).unwrap();

    let err = w.write_raw(&rec).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidData);
    assert!(err.to_string().contains(why), "{err}");
}

#[test]
fn decoding_against_a_different_reference_is_an_error() {
    let (dir, _fa, full) = written_cram();
    // Same contig name and length, different bases: the slices' reference MD5s
    // do not match it.
    let other_dir = dir.path().join("other");
    std::fs::create_dir(&other_dir).unwrap();
    let other = write_fasta(&other_dir, "chrT", &reference(REF_LEN, 0x9E37_79B9_7F4A_7C15));

    assert!(read_cram(&full, Some(&other)).is_err(), "a reference MD5 mismatch must fail");
}

#[rstest]
#[case::missing_reference(Some("missing.fa"), None, HEADER)]
#[case::nul_in_file_id(None, Some("id\0x"), HEADER)]
#[case::nul_in_header_text(None, None, "@HD\tVN:1.6\0\n")]
fn bad_writer_inputs_leave_an_existing_output_untouched(
    #[case] reference: Option<&str>,
    #[case] file_id: Option<&str>,
    #[case] header: &str,
) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("existing.cram");
    std::fs::write(&out, b"precious").unwrap();
    let reference = reference.map(|r| dir.path().join(r));
    let opts = CramWriteOptions {
        reference: reference.as_deref(),
        no_ref: reference.is_none(),
        file_id,
        ..Default::default()
    };

    assert!(CramWriter::create(&out, header, &opts).is_err());
    assert_eq!(std::fs::read(&out).unwrap(), b"precious", "the output was truncated");
}
