//! A CRAM that is cut short, or whose write was abandoned, must never read as
//! complete: htslib alone would end such a stream cleanly or stamp an EOF
//! container on it.
mod common;

use common::{read_cram, record, reference, write_cram, write_fasta};
use std::io::{self, ErrorKind, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

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

/// A valid reader that fails instead of reporting physical EOF, and records
/// that it did.
struct ErrorAtEnd(io::Cursor<Vec<u8>>, Arc<AtomicBool>);
impl Read for ErrorAtEnd {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.0.position() == self.0.get_ref().len() as u64 {
            self.1.store(true, Ordering::SeqCst);
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
    let injected = Arc::new(AtomicBool::new(false));
    let mut reader = CramReader::from_reader(
        Box::new(ErrorAtEnd(io::Cursor::new(bytes), Arc::clone(&injected))),
        "failing-reader",
        Some(&fa),
    )
    .unwrap();
    loop {
        match reader.next_slice() {
            Ok(Some(_)) => {}
            Err(err) => {
                assert!(
                    injected.load(Ordering::SeqCst),
                    "the error must come from the read after the footer: {err}"
                );
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

#[test]
fn framing_error_is_not_cleared_by_a_retry() {
    // Two valid EOF containers: the first read past the footer fails, and a
    // retry must not read on to the second footer and report a clean end.
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(full).unwrap();
    bytes.extend_from_within(bytes.len() - usize::try_from(EOF_LEN).unwrap()..);
    let padded = dir.path().join("duplicate-eof.cram");
    std::fs::write(&padded, bytes).unwrap();
    let mut reader = CramReader::open(&padded, Some(&fa)).unwrap();
    let first = loop {
        match reader.next_slice() {
            Ok(Some(_)) => {}
            Ok(None) => panic!("a duplicate EOF container must not read as complete"),
            Err(e) => break e,
        }
    };
    assert!(first.to_string().contains("stray bytes after its EOF container"), "{first}");
    for _ in 0..2 {
        let Err(again) = reader.next_slice() else {
            panic!("a retry after a framing error must fail again");
        };
        assert_eq!(again.kind(), ErrorKind::InvalidData);
        assert!(again.to_string().contains("stray bytes after its EOF container"), "{again}");
        assert!(again.to_string().contains("an earlier read already failed"), "{again}");
    }
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

#[test]
fn data_container_with_unfilled_length_is_an_error() {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(full).unwrap();
    // The first data container (the header container comes first): declare
    // 8 more payload bytes than its blocks fill, and pad it to match, so
    // every block and record still reads.
    let c = &common::cram_layout(&bytes)[1];
    let (offset, header_end, end) =
        (usize::try_from(c.offset).unwrap(), c.payload_offset - 4, c.end);
    let length = u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap()) + 8;
    bytes[offset..offset + 4].copy_from_slice(&length.to_le_bytes());
    let checksum = crc32(&bytes[offset..header_end]);
    bytes[header_end..header_end + 4].copy_from_slice(&checksum.to_le_bytes());
    bytes.splice(end..end, [0u8; 8]);
    let padded = dir.path().join("padded-container.cram");
    std::fs::write(&padded, bytes).unwrap();
    let Err(err) = read_cram(&padded, Some(&fa)) else {
        panic!("a container whose blocks do not fill its length must not read as complete");
    };
    assert_eq!(err.kind(), ErrorKind::InvalidData);
    assert!(err.to_string().contains("do not fill its declared length"), "{err}");
}

/// A data container whose declared record count disagrees with its slices is
/// rejected: with zero records it is skipped only if its blocks still fill its
/// length (they do not: its slices follow), and otherwise the slices' counts
/// must sum to the container's.
#[rstest]
#[case::zero_records(|_| 0, "do not fill its declared length")]
#[case::one_more_record(|n| n + 1, "do not hold the number of records it declares")]
#[case::one_fewer_record(|n| n - 1, "do not hold the number of records it declares")]
fn data_container_with_a_wrong_record_count_is_an_error(
    #[case] count: fn(u32) -> u32,
    #[case] why: &str,
) {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(full).unwrap();
    let c = &common::cram_layout(&bytes)[1];
    let offset = usize::try_from(c.offset).unwrap();
    let (start, end, n) = common::container_record_count(&bytes, offset);
    assert!(n > 1, "the first data container holds records");
    let encoded = common::itf8_bytes(count(n));
    let header_end = c.payload_offset - 4 + encoded.len() - (end - start);
    bytes.splice(start..end, encoded);
    let checksum = crc32(&bytes[offset..header_end]);
    bytes[header_end..header_end + 4].copy_from_slice(&checksum.to_le_bytes());
    let bad = dir.path().join("bad-record-count.cram");
    std::fs::write(&bad, bytes).unwrap();

    let Err(err) = read_cram(&bad, Some(&fa)) else {
        panic!("a container whose record count disagrees with its slices must not read");
    };
    assert_eq!(err.kind(), ErrorKind::InvalidData);
    assert!(err.to_string().contains(why), "{err}");
}

/// A data container's declared block count is not checked: htsjdk (Picard,
/// GATK) writes it without the compression-header and slice-header blocks, so
/// a valid htsjdk CRAM declares two fewer blocks than it holds and must read.
#[test]
fn a_block_count_without_header_blocks_still_reads() {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(&full).unwrap();
    let c = &common::cram_layout(&bytes)[1];
    let offset = usize::try_from(c.offset).unwrap();
    let (start, end, n) = common::container_block_count(&bytes, offset);
    assert_eq!(n as usize, c.blocks.len(), "the writer declares every block");
    let encoded = common::itf8_bytes(n - 2);
    let header_end = c.payload_offset - 4 + encoded.len() - (end - start);
    bytes.splice(start..end, encoded);
    let checksum = crc32(&bytes[offset..header_end]);
    bytes[header_end..header_end + 4].copy_from_slice(&checksum.to_le_bytes());
    let htsjdk_style = dir.path().join("htsjdk-block-count.cram");
    std::fs::write(&htsjdk_style, bytes).unwrap();

    assert_eq!(
        read_cram(&htsjdk_style, Some(&fa)).unwrap().0,
        read_cram(&full, Some(&fa)).unwrap().0
    );
}

/// A container header with no payload, inserted before the first data
/// container: skipped when it declares no records, rejected when it declares
/// records it cannot hold.
#[rstest]
#[case::declares_no_records(0, None)]
#[case::declares_records(5, Some("do not hold the number of records it declares"))]
fn zero_length_container(#[case] n_records: u32, #[case] why: Option<&str>) {
    let (dir, fa, full) = written_cram();
    let mut bytes = std::fs::read(&full).unwrap();
    let offset = usize::try_from(common::cram_layout(&bytes)[1].offset).unwrap();
    let mut header = 0u32.to_le_bytes().to_vec(); // length
    header.extend_from_slice(&[0xff, 0xff, 0xff, 0xff, 0x0f]); // reference id -1
    header.extend_from_slice(&[0, 0]); // start, span
    header.extend(common::itf8_bytes(n_records));
    header.extend_from_slice(&[0, 0]); // record counter, bases (LTF8)
    header.extend_from_slice(&[0, 0]); // no blocks, no landmarks
    let checksum = crc32(&header);
    header.extend_from_slice(&checksum.to_le_bytes());
    bytes.splice(offset..offset, header);
    let empty = dir.path().join("empty-container.cram");
    std::fs::write(&empty, bytes).unwrap();

    let result = read_cram(&empty, Some(&fa));
    match why {
        None => assert_eq!(result.unwrap().0, read_cram(&full, Some(&fa)).unwrap().0),
        Some(why) => {
            let Err(err) = result else {
                panic!("an empty container that declares contents must not read");
            };
            assert_eq!(err.kind(), ErrorKind::InvalidData);
            assert!(err.to_string().contains(why), "{err}");
        }
    }
}

/// A reader that fails once, then serves `bytes`.
struct FailsFirst(bool, io::Cursor<Vec<u8>>);
impl Read for FailsFirst {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if std::mem::replace(&mut self.0, false) {
            Err(io::Error::other("injected read failure before the file definition"))
        } else {
            self.1.read(buf)
        }
    }
}

/// A read failure while probing the version must fail the open, not let
/// htslib probe again and open a CRAM version the probe would have refused.
#[test]
fn read_failure_while_probing_the_version_fails_the_open() {
    let mut bytes = b"CRAM\x02\x00".to_vec();
    bytes.extend_from_slice(&[0; 20]);
    bytes.extend(empty_cram_2_1().split_off(26));
    let Err(err) = CramReader::from_reader(
        Box::new(FailsFirst(true, io::Cursor::new(bytes))),
        "fails-first",
        None,
    ) else {
        panic!("an open whose version probe failed must fail");
    };
    assert!(err.to_string().contains("reading the CRAM file definition failed"), "{err}");
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
    let err = second.write().unwrap_err();
    assert!(
        err.to_string().contains("CRAM container 1 written out of order (expected 0)"),
        "{err}"
    );
    drop(jobs);
    let (rest, fin) = w.seal().unwrap();
    drop(rest);
    let err = fin.finish().unwrap_err();
    assert!(err.to_string().contains("containers written"), "{err}");

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

#[cfg(unix)]
#[test]
fn failed_container_write_refuses_later_writes_and_the_finish() {
    const CHILD: &str = "FGUMI_TEST_FAILED_CONTAINER_WRITE";
    if std::env::var_os(CHILD).is_some() {
        let dir = tempfile::tempdir().unwrap();
        let seq = reference(REF_LEN, 0x2545_F491_4F6C_DD1D);
        let fa = write_fasta(dir.path(), "chrT", &seq);
        let out = dir.path().join("failed_write.cram");
        let opts = CramWriteOptions { reference: Some(&fa), ..Default::default() };
        let mut w = CramWriter::create(&out, HEADER, &opts).unwrap();
        for r in records(&seq) {
            w.write_raw(&r).unwrap();
        }
        let jobs: Vec<_> = w.drain_jobs().map(|j| j.encode().unwrap()).collect();
        let (tail, fin) = w.seal().unwrap();
        let tail: Vec<_> = tail.into_iter().map(|j| j.encode().unwrap()).collect();
        assert!(jobs.len() + tail.len() >= 2, "need several containers");
        // No byte can reach the file from here on.
        rustix::process::setrlimit(
            rustix::process::Resource::Fsize,
            rustix::process::Rlimit { current: Some(0), maximum: Some(0) },
        )
        .unwrap();
        let mut all = jobs.into_iter().chain(tail);
        // htslib buffers output, so the first failing write may be a later one.
        let failed = loop {
            let c = all.next().expect("some container write must fail");
            let serial = c.serial();
            if c.write().is_err() {
                break serial;
            }
        };
        let rest: Vec<_> = all.collect();
        assert!(!rest.is_empty(), "container {failed} failed last: no later write to refuse");
        for c in rest {
            let err = c.write().unwrap_err();
            assert!(err.to_string().contains("an earlier container write failed"), "{err}");
        }
        let err = fin.finish().unwrap_err();
        assert!(err.to_string().contains("after a failed container write"), "{err}");
        eprintln!("container {failed} failed");
        return;
    }

    let output = std::process::Command::new("sh")
        .args(["-c", "trap '' XFSZ; exec \"$1\" --exact failed_container_write_refuses_later_writes_and_the_finish --nocapture", "sh"])
        .arg(std::env::current_exe().unwrap())
        .env(CHILD, "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "failed-write child failed: {}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
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

    // A rejected record leaves the writer usable: a valid one still writes
    // and the output finishes holding only it.
    let valid = raw(0, &[50 << 4], 50);
    w.write_raw(&valid).unwrap();
    let mut jobs: Vec<_> = w.drain_jobs().collect();
    let (tail, fin) = w.seal().unwrap();
    jobs.extend(tail);
    for job in jobs {
        job.encode().unwrap().write().unwrap();
    }
    fin.finish().unwrap();
    let (decoded, _) = read_cram(&out, None).unwrap();
    assert_eq!(decoded, vec![valid], "the output holds exactly the valid record");
}

#[test]
fn decoding_against_a_different_reference_is_an_error() {
    let (dir, fa, full) = written_cram();
    // Same contig name and length, different bases: the slices' reference MD5s
    // do not match it.
    let other_dir = dir.path().join("other");
    std::fs::create_dir(&other_dir).unwrap();
    let other = write_fasta(&other_dir, "chrT", &reference(REF_LEN, 0x9E37_79B9_7F4A_7C15));

    // The file and its first (mapped) slice read with either reference; only
    // decoding the slice checks the reference, and only the right one passes.
    let first_slice = |reference: &std::path::Path| {
        let mut reader = CramReader::open(&full, Some(reference)).unwrap();
        reader.next_slice().unwrap().expect("the file holds a slice")
    };
    assert!(first_slice(&fa).decode().unwrap().records > 0);
    let Err(err) = first_slice(&other).decode() else {
        panic!("a reference MD5 mismatch must fail the slice decode");
    };
    assert!(err.to_string().contains("CRAM slice decode failed"), "{err}");
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
