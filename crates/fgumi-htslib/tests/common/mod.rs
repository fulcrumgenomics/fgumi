//! Shared fixtures for the CRAM tests: a deterministic reference, raw BAM
//! record encoding (incl. BAM's long-CIGAR form, which `fgumi_raw_bam`'s
//! `SamBuilder` cannot build: it caps a CIGAR at 65535 ops), and write/read helpers that
//! run the codec the way fgumi's pipeline does (encodes on worker threads,
//! in-order writes, parallel slice decodes).
#![allow(dead_code)] // each test binary uses a different subset

use std::io::Write;
use std::path::{Path, PathBuf};

use fgumi_htslib::{CramReader, CramWriteOptions, CramWriter, EncodeJob};
use fgumi_raw_bam::{RawBamReader, RawRecord, RawTagsView, pack_sequence_into, reg2bin};
use fgumi_tag::SamTag;

/// Deterministic pseudo-random ACGT of length `len`.
pub fn reference(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed | 1;
    (0..len)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[usize::try_from(x % 4).unwrap()]
        })
        .collect()
}

/// Write `seq` as a one-contig FASTA named `contig` in `dir`.
pub fn write_fasta(dir: &Path, contig: &str, seq: &[u8]) -> PathBuf {
    let fa = dir.join("ref.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">{contig}").unwrap();
    for line in seq.chunks(60) {
        f.write_all(line).unwrap();
        writeln!(f).unwrap();
    }
    fa
}

/// One raw BAM record (no `block_size`). Mapped when `tid >= 0`; more than
/// 65535 CIGAR ops are stored in BAM's long-CIGAR form (`<qlen>S<rlen>N`
/// placeholder plus a trailing `CG:B,I` tag).
pub fn record(
    name: &str,
    flag: u16,
    tid: i32,
    pos: i32,
    cigar: &[u32],
    seq: &[u8],
    tags: &[u8],
) -> Vec<u8> {
    let long = cigar.len() > 0xffff;
    let rlen: u32 =
        cigar.iter().filter(|op| matches!(*op & 0xf, 0 | 2 | 3 | 7 | 8)).map(|op| op >> 4).sum();
    let mapped = flag & 4 == 0;
    let bin = if mapped {
        u32::from(reg2bin(pos, pos + i32::try_from(rlen.max(1)).unwrap()))
    } else {
        4680
    };
    let mapq: u32 = if mapped { 60 } else { 0 };
    let n_cigar = if long { 2 } else { u32::try_from(cigar.len()).unwrap() };
    let mut r = Vec::new();
    r.extend_from_slice(&tid.to_le_bytes());
    r.extend_from_slice(&pos.to_le_bytes());
    r.extend_from_slice(
        &(bin << 16 | mapq << 8 | u32::try_from(name.len() + 1).unwrap()).to_le_bytes(),
    );
    r.extend_from_slice(&(u32::from(flag) << 16 | n_cigar).to_le_bytes());
    r.extend_from_slice(&i32::try_from(seq.len()).unwrap().to_le_bytes());
    r.extend_from_slice(&(-1i32).to_le_bytes());
    r.extend_from_slice(&(-1i32).to_le_bytes());
    r.extend_from_slice(&0i32.to_le_bytes());
    r.extend_from_slice(name.as_bytes());
    r.push(0);
    if long {
        r.extend_from_slice(&(u32::try_from(seq.len()).unwrap() << 4 | 4).to_le_bytes());
        r.extend_from_slice(&(rlen << 4 | 3).to_le_bytes());
    } else {
        for op in cigar {
            r.extend_from_slice(&op.to_le_bytes());
        }
    }
    pack_sequence_into(&mut r, seq);
    r.extend(seq.iter().enumerate().map(|(i, _)| 20 + u8::try_from(i % 20).unwrap()));
    r.extend_from_slice(tags);
    if long {
        r.extend_from_slice(b"CGBI");
        r.extend_from_slice(&u32::try_from(cigar.len()).unwrap().to_le_bytes());
        for op in cigar {
            r.extend_from_slice(&op.to_le_bytes());
        }
    }
    r
}

/// Encode jobs on worker threads, then write them in order.
pub fn encode_and_write(jobs: Vec<EncodeJob>) {
    let encoded: Vec<_> = std::thread::scope(|s| {
        let handles: Vec<_> =
            jobs.into_iter().map(|j| s.spawn(move || j.encode().unwrap())).collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for e in encoded {
        e.write().unwrap();
    }
}

/// Write `records` as CRAM at `path`; returns the number of containers.
pub fn write_cram(
    path: &Path,
    header: &str,
    opts: &CramWriteOptions<'_>,
    records: &[Vec<u8>],
) -> u64 {
    let mut w = CramWriter::create(path, header, opts).unwrap();
    for chunk in records.chunks(4096) {
        for r in chunk {
            w.write_raw(r).unwrap();
        }
        encode_and_write(w.drain_jobs().collect());
    }
    let (jobs, fin) = w.seal().unwrap();
    encode_and_write(jobs);
    let n = fin.total_containers();
    fin.finish().unwrap();
    n
}

/// Read every record of a CRAM (slices decoded in parallel), grouped by the
/// slice it was decoded from, in file order.
pub fn read_cram_slices(
    path: &Path,
    reference: Option<&Path>,
) -> std::io::Result<Vec<Vec<Vec<u8>>>> {
    let mut r = CramReader::open(path, reference)?;
    let mut slices = Vec::new();
    while let Some(s) = r.next_slice()? {
        slices.push(s);
    }
    let decoded: Vec<_> = std::thread::scope(|s| {
        let hs: Vec<_> = slices.into_iter().map(|j| s.spawn(move || j.decode())).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect::<std::io::Result<Vec<_>>>()
    })?;
    let mut out = Vec::new();
    for d in &decoded {
        let mut reader = RawBamReader::new(&d.bytes[..]);
        let mut rec = RawRecord::new();
        let mut slice = Vec::new();
        while reader.read_record(&mut rec)? > 0 {
            slice.push(rec.to_vec());
        }
        out.push(slice);
    }
    Ok(out)
}

/// Read every record of a CRAM (slices decoded in parallel); returns the raw
/// records (no `block_size`) and the number of slices.
pub fn read_cram(path: &Path, reference: Option<&Path>) -> std::io::Result<(Vec<Vec<u8>>, usize)> {
    let slices = read_cram_slices(path, reference)?;
    let n_slices = slices.len();
    Ok((slices.into_iter().flatten().collect(), n_slices))
}

/// The `MD` string and `NM` count of a read aligned without indels (CIGAR
/// all `M`) at the start of `ref_segment`, computed independently of htslib.
pub fn expected_md_nm(read: &[u8], ref_segment: &[u8]) -> (String, i64) {
    let mut md = String::new();
    let mut run = 0usize;
    let mut nm = 0i64;
    for (&q, &r) in read.iter().zip(ref_segment) {
        if q == r {
            run += 1;
        } else {
            md.push_str(&run.to_string());
            md.push(char::from(r));
            run = 0;
            nm += 1;
        }
    }
    md.push_str(&run.to_string());
    (md, nm)
}

/// The values of the regenerated `MD:Z` and `NM` tags `decoded` carries
/// beyond `original`; panics unless exactly those two tags were appended.
pub fn appended_md_nm(decoded: &[u8], original: &[u8]) -> (String, i64) {
    let appended: Vec<[u8; 2]> = appended_tags(decoded, original).iter().map(|(t, _)| *t).collect();
    assert_eq!(
        appended,
        vec![<[u8; 2]>::from(SamTag::MD), SamTag::NM.into()],
        "exactly MD and NM are appended"
    );
    let extra = RawTagsView::new(&decoded[original.len()..]);
    let md = String::from_utf8(extra.find_string(SamTag::MD).unwrap().to_vec()).unwrap();
    (md, extra.find_int(SamTag::NM).unwrap())
}

/// The aux tags `decoded` carries beyond `original` (which it must start with),
/// as `(tag, type)` pairs. CRAM 3.x decoding (like `samtools view`) appends
/// regenerated `MD:Z` and `NM:i` to mapped records.
pub fn appended_tags(decoded: &[u8], original: &[u8]) -> Vec<([u8; 2], u8)> {
    assert!(decoded.starts_with(original), "decoded record does not extend the original");
    let extra = &decoded[original.len()..];
    let entries: Vec<_> = RawTagsView::new(extra).iter().collect();
    let parsed: usize = entries.iter().map(|e| 3 + e.value_bytes.len()).sum();
    assert_eq!(parsed, extra.len(), "the appended tags do not parse to the record's end");
    entries.iter().map(|e| (e.tag, e.type_byte)).collect()
}

/// A CRAM 3 container's physical block boundaries, read independently of the codec.
#[derive(Debug)]
pub struct ContainerLayout {
    pub offset: u64,
    pub payload_offset: usize,
    pub end: usize,
    pub blocks: Vec<BlockLayout>,
}

/// One on-disk block, including its serialized header, payload and checksum.
#[derive(Debug)]
pub struct BlockLayout {
    pub start: usize,
    pub end: usize,
    pub method: u8,
    pub content_type: u8,
}

/// Consume a CRAM 3 ITF8 integer, retaining its unsigned bit pattern.
fn itf8(bytes: &[u8], cursor: &mut usize) -> u32 {
    let first = bytes[*cursor];
    *cursor += 1;
    let trailing = first.leading_ones().min(4);
    let mask = if trailing == 4 { 0xf } else { 0x7f >> trailing };
    let mut value = u32::from(first & mask);
    for i in 0..trailing {
        let byte = bytes[*cursor];
        *cursor += 1;
        if trailing == 4 && i == 3 {
            value = (value << 4) | u32::from(byte & 0xf);
        } else {
            value = (value << 8) | u32::from(byte);
        }
    }
    value
}

/// The byte range and value of the record count in the container header at
/// `offset` (after its length, reference id, start and span).
pub fn container_record_count(bytes: &[u8], offset: usize) -> (usize, usize, u32) {
    let mut cursor = offset + 4;
    for _ in 0..3 {
        itf8(bytes, &mut cursor);
    }
    let start = cursor;
    let value = itf8(bytes, &mut cursor);
    (start, cursor, value)
}

/// The byte range and value of the block count in the container header at
/// `offset` (after its length, reference id, start, span, record count,
/// record counter and base count).
pub fn container_block_count(bytes: &[u8], offset: usize) -> (usize, usize, u32) {
    let mut cursor = offset + 4;
    for _ in 0..4 {
        itf8(bytes, &mut cursor);
    }
    skip_ltf8(bytes, &mut cursor);
    skip_ltf8(bytes, &mut cursor);
    let start = cursor;
    let value = itf8(bytes, &mut cursor);
    (start, cursor, value)
}

/// `value` as a CRAM 3 ITF8 integer (values below 2^21 suffice here).
pub fn itf8_bytes(value: u32) -> Vec<u8> {
    let [_, b1, b2, b3] = value.to_be_bytes();
    match value {
        0..0x80 => vec![b3],
        0x80..0x4000 => vec![0x80 | b2, b3],
        0x4000..0x20_0000 => vec![0xc0 | b1, b2, b3],
        _ => panic!("ITF8 value {value} too large for this helper"),
    }
}

/// Skip an LTF8 integer (the container's record/base counts).
fn skip_ltf8(bytes: &[u8], cursor: &mut usize) {
    let trailing = bytes[*cursor].leading_ones();
    *cursor += 1 + usize::try_from(trailing).unwrap();
}

/// Walk serialized containers and blocks without consulting the emitted CRAI.
/// Only CRAM 3 is needed: it is the writer's output format under test.
pub fn cram_layout(bytes: &[u8]) -> Vec<ContainerLayout> {
    assert_eq!(&bytes[..4], b"CRAM");
    assert_eq!(bytes[4], 3);
    let mut containers = Vec::new();
    let mut cursor = 26;
    while cursor < bytes.len() {
        let offset = cursor as u64;
        let length =
            usize::try_from(u32::from_le_bytes(bytes[cursor..cursor + 4].try_into().unwrap()))
                .unwrap();
        cursor += 4;
        // Reference id, start, span, record count.
        for _ in 0..4 {
            itf8(bytes, &mut cursor);
        }
        skip_ltf8(bytes, &mut cursor);
        skip_ltf8(bytes, &mut cursor);
        let n_blocks = itf8(bytes, &mut cursor);
        let n_landmarks = itf8(bytes, &mut cursor);
        for _ in 0..n_landmarks {
            itf8(bytes, &mut cursor);
        }
        cursor += 4; // Container header CRC32.
        let payload_offset = cursor;
        let end = payload_offset + length;
        let mut blocks = Vec::new();
        for _ in 0..n_blocks {
            let start = cursor;
            let method = bytes[cursor];
            let content_type = bytes[cursor + 1];
            cursor += 2;
            itf8(bytes, &mut cursor); // Content id.
            let compressed = itf8(bytes, &mut cursor);
            itf8(bytes, &mut cursor); // Uncompressed size.
            cursor += usize::try_from(compressed).unwrap() + 4; // Payload and CRC32.
            blocks.push(BlockLayout { start, end: cursor, method, content_type });
        }
        assert_eq!(cursor, end, "blocks must cover the container payload");
        containers.push(ContainerLayout { offset, payload_offset, end, blocks });
    }
    assert_eq!(cursor, bytes.len());
    containers
}

/// CRAI entries: (ref id, alignment start, span, container offset, slice offset, slice size).
pub fn crai_entries(crai: &Path) -> Vec<(i64, i64, i64, u64, u64, u64)> {
    use std::io::Read;
    let mut text = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(crai).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text.lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 6, "CRAI line has 6 fields: {l}");
            (
                f[0].parse().unwrap(),
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
                f[5].parse().unwrap(),
            )
        })
        .collect()
}
