//! FASTQ-writing helpers shared by the `runall` integration tests.
//!
//! `runall`'s extract-fusion tests (`test_runall_command.rs`,
//! `test_runall_chain_transitions.rs`) both need to stage small gzip-compressed
//! FASTQ fixtures; this used to be duplicated verbatim in both files.

use std::fs;
use std::io::Write as _;
use std::path::Path;

/// Write a gzip-compressed FASTQ from `(name, seq, qual)` records.
pub fn write_gzip_fastq(path: &Path, records: &[(&str, &str, &str)]) {
    let file = fs::File::create(path).expect("create gzip fastq");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for (name, seq, qual) in records {
        writeln!(encoder, "@{name}\n{seq}\n+\n{qual}").expect("write fastq record");
    }
    encoder.finish().expect("finish gzip fastq");
}

/// Write `fastq` (raw FASTQ text) BGZF-compressed to `path`, then flip one
/// byte of the CRC32 footer of the last *data* block. Callers pass enough
/// records to span several blocks, so the corruption can sit past any reader's
/// sampling window.
pub fn write_bgzf_fastq_with_corrupt_last_crc(path: &Path, fastq: &[u8]) {
    let mut bytes = Vec::new();
    {
        let mut writer = noodles::bgzf::io::Writer::new(&mut bytes);
        writer.write_all(fastq).expect("write bgzf fastq");
        writer.finish().expect("finish bgzf fastq");
    }
    let blocks = {
        let mut cursor: &[u8] = &bytes;
        fgumi_bgzf::read_raw_blocks(&mut cursor, 1_000_000).expect("read bgzf blocks")
    };
    // `read_raw_blocks` skips the empty EOF marker, so the returned blocks are
    // all data blocks and their total length ends at the last one's footer.
    assert!(blocks.len() >= 2, "input must span >= 2 data blocks; got {}", blocks.len());
    let data_end: usize = blocks.iter().map(fgumi_bgzf::RawBgzfBlock::len).sum();
    bytes[data_end - fgumi_bgzf::BGZF_FOOTER_SIZE] ^= 0x01;
    fs::write(path, bytes).expect("write corrupted bgzf fastq");
}
