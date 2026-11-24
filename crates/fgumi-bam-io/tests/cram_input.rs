//! CRAM input through the normalization boundary (feature `cram`): a CRAM is
//! read as the same BGZF BAM stream every fgumi reader consumes, from a file
//! and from a non-seekable stream whose leading bytes were already sniffed.
#![cfg(feature = "cram")]

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use fgumi_htslib::{CramWriteOptions, CramWriter};
use fgumi_raw_bam::{SamBuilder, reg2bin};
use rstest::rstest;

const REF_LEN: usize = 5_000;
const N_READS: usize = 12_000;

/// One shared fixture: the CRAM reference is process-global, so every case
/// must decode against the same FASTA.
struct Fixture {
    _dir: tempfile::TempDir,
    cram: PathBuf,
    records: Vec<Vec<u8>>,
}

fn fixture() -> &'static Fixture {
    static FIXTURE: OnceLock<Fixture> = OnceLock::new();
    FIXTURE.get_or_init(|| {
        let dir = tempfile::tempdir().unwrap();
        let seq = reference();
        let fa = write_fasta(dir.path(), &seq);
        let cram = dir.path().join("in.cram");
        let header = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:{REF_LEN}\n");
        let mut w = CramWriter::create(
            &cram,
            &header,
            &CramWriteOptions { reference: Some(&fa), ..Default::default() },
        )
        .unwrap();
        let mut records = Vec::new();
        for i in 0..N_READS {
            let pos = (i * 3) % (REF_LEN - 100);
            let r = record(i, i32::try_from(pos).unwrap(), &seq[pos..pos + 50]);
            w.write_raw(&r).unwrap();
            records.push(r);
            if i % 1000 == 999 {
                for j in w.drain_jobs() {
                    j.encode().unwrap().write().unwrap();
                }
            }
        }
        let (jobs, fin) = w.seal().unwrap();
        for j in jobs {
            j.encode().unwrap().write().unwrap();
        }
        fin.finish().unwrap();
        fgumi_bam_io::cram_input::set_cram_reference(&fa);
        Fixture { _dir: dir, cram, records }
    })
}

fn reference() -> Vec<u8> {
    let mut x: u64 = 0x2545_F491_4F6C_DD1D;
    (0..REF_LEN)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[usize::try_from(x % 4).unwrap()]
        })
        .collect()
}

fn write_fasta(dir: &Path, seq: &[u8]) -> PathBuf {
    let fa = dir.join("ref.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">chrT").unwrap();
    for line in seq.chunks(60) {
        f.write_all(line).unwrap();
        writeln!(f).unwrap();
    }
    fa
}

/// One mapped 50M record (no `block_size`).
fn record(i: usize, pos: i32, seq: &[u8]) -> Vec<u8> {
    let len = i32::try_from(seq.len()).unwrap();
    SamBuilder::new()
        .ref_id(0)
        .pos(pos)
        .mapq(30)
        .bin(reg2bin(pos, pos + len))
        .read_name(format!("r{i:05}").as_bytes())
        .cigar_ops(&[u32::try_from(seq.len()).unwrap() << 4]) // 50M
        .sequence(seq)
        .qualities(&vec![35u8; seq.len()])
        .build()
        .into_inner()
}

/// A reader that hands out at most 7 bytes per call (a slow pipe).
struct Trickle<R>(R);
impl<R: Read> Read for Trickle<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = buf.len().min(7);
        self.0.read(&mut buf[..n])
    }
}

/// The reference dictionary (name, length) and the raw records (no
/// `block_size`) of an uncompressed BAM byte stream.
fn bam_records(bam: &[u8]) -> (Vec<(String, usize)>, Vec<Vec<u8>>) {
    assert_eq!(&bam[..4], b"BAM\x01");
    let u32_at = |b: &[u8], i: usize| u32::from_le_bytes(b[i..i + 4].try_into().unwrap()) as usize;
    let mut i = 8 + u32_at(bam, 4);
    let n_ref = u32_at(bam, i);
    i += 4;
    let mut refs = Vec::with_capacity(n_ref);
    for _ in 0..n_ref {
        let l_name = u32_at(bam, i);
        // The name is NUL-terminated within `l_name` bytes.
        let name = String::from_utf8(bam[i + 4..i + 4 + l_name - 1].to_vec()).unwrap();
        refs.push((name, u32_at(bam, i + 4 + l_name)));
        i += 4 + l_name + 4;
    }
    let mut out = Vec::new();
    while i < bam.len() {
        let n = u32_at(bam, i);
        out.push(bam[i + 4..i + 4 + n].to_vec());
        i += 4 + n;
    }
    (refs, out)
}

#[rstest]
#[case::file(false)]
#[case::trickling_stream(true)]
fn cram_reads_as_bam(#[case] stream: bool) {
    let fx = fixture();
    let input: Box<dyn Read + Send> = if stream {
        Box::new(Trickle(std::fs::File::open(&fx.cram).unwrap()))
    } else {
        Box::new(std::fs::File::open(&fx.cram).unwrap())
    };
    let normalized = fgumi_bam_io::sam_input::normalize_to_bgzf(input, &fx.cram).unwrap();
    let mut bam = Vec::new();
    noodles::bgzf::io::Reader::new(normalized).read_to_end(&mut bam).unwrap();

    let (refs, got) = bam_records(&bam);
    // Reference ID 0 in every record must still name the sequence it was
    // written against.
    assert_eq!(refs, vec![("chrT".to_string(), REF_LEN)]);
    assert_eq!(got.len(), N_READS);
    for (i, (g, e)) in got.iter().zip(&fx.records).enumerate() {
        // Each record is the one written, with exactly MD:Z and NM appended
        // (CRAM 3 decoding regenerates them, as `samtools view` does). Every
        // read matches the reference over all 50 bases: MD:Z:50 and NM:0.
        assert!(g.starts_with(e), "record {i} changed");
        let extra = &g[e.len()..];
        assert_eq!(&extra[..6], b"MDZ50\0", "record {i}: MD:Z:50 not appended");
        assert_eq!(
            &extra[6..8],
            AsRef::<[u8; 2]>::as_ref(&fgumi_tag::SamTag::NM),
            "record {i}: NM not appended"
        );
        let nm = &extra[9..];
        let nm_len = match extra[8] {
            b'c' | b'C' => 1,
            b's' | b'S' => 2,
            b'i' | b'I' => 4,
            t => panic!("record {i}: NM has non-integer type {:?}", t as char),
        };
        assert_eq!(nm.len(), nm_len, "record {i}: unexpected appended tags");
        assert!(nm.iter().all(|&b| b == 0), "record {i}: NM is not 0");
    }
}

/// The CRAM 3.x EOF container's length in bytes.
const CRAM_EOF_LEN: usize = 38;

/// How a truncated copy of the fixture is cut.
#[derive(Debug, Clone, Copy)]
enum Truncation {
    /// Every container intact, but the trailing EOF container is gone.
    MissingEof,
    /// The stream stops halfway through, inside a data container.
    MidContainer,
}

/// A truncated CRAM must read as an error, never as a clean (short) BAM: the
/// reader thread's hang-up must not be mistaken for end of file.
#[rstest]
fn truncated_cram_is_an_error(
    #[values(Truncation::MissingEof, Truncation::MidContainer)] cut: Truncation,
    #[values(false, true)] stream: bool,
) {
    let fx = fixture();
    let full = std::fs::read(&fx.cram).unwrap();
    let len = match cut {
        Truncation::MissingEof => full.len() - CRAM_EOF_LEN,
        Truncation::MidContainer => full.len() / 2,
    };
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("truncated.cram");
    std::fs::write(&path, &full[..len]).unwrap();

    let input: Box<dyn Read + Send> = if stream {
        Box::new(Trickle(std::fs::File::open(&path).unwrap()))
    } else {
        Box::new(std::fs::File::open(&path).unwrap())
    };
    let result = fgumi_bam_io::sam_input::normalize_to_bgzf(input, &path).and_then(|normalized| {
        let mut bam = Vec::new();
        noodles::bgzf::io::Reader::new(normalized).read_to_end(&mut bam)?;
        Ok(bam)
    });
    let Err(err) = result else {
        panic!("{cut:?} (stream: {stream}): a truncated CRAM read as a complete BAM");
    };
    let msg = format!("{err:#}");
    assert!(msg.contains("truncated.cram"), "error does not name the input: {msg}");
    // A read failure is not a missing reference (and one is configured here).
    assert!(!msg.contains("no reference configured"), "misleading reference hint: {msg}");
}
