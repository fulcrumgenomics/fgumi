//! Byte-comparison: runall `--start-from extract --stop-after filter` with
//! `--consensus duplex` matches the equivalent standalone chain:
//! `fgumi extract --read-structures 8M+T 8M+T` | `fgumi fastq` | `bwa mem` |
//! `fgumi zipper` | `fgumi sort --order template-coordinate` |
//! `fgumi group --strategy paired` | `fgumi duplex --min-reads 1` |
//! `fgumi filter --min-reads 1 --min-base-quality 2`.
//!
//! Requires `bwa` and `samtools` on PATH. Skips gracefully if either is
//! missing. Compares via `fgumi compare bams --command filter` so header `@PG`
//! differences between the two command chains are ignored.
//!
//! Paired UMIs are produced by extracting an 8-base UMI from the start of both
//! R1 and R2, yielding an `RX` tag of the form `"AAAAAAAA-CCCCCCCC"`. With
//! `--group::strategy paired`, reads from opposite strands are assigned MI
//! tags with `/A` and `/B` suffixes so the duplex caller can combine them.

use std::io::Write;
use std::path::Path;
use std::process::Command;

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

const DNA: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn ref_sequence() -> Vec<u8> {
    (0..1000).map(|i| DNA[i % 4]).collect()
}

fn write_reference(path: &Path) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, ">chr_test").unwrap();
    let seq = ref_sequence();
    for chunk in seq.chunks(80) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
}

fn bwa_index(ref_path: &Path) {
    let output = Command::new("bwa")
        .args(["index", ref_path.to_str().unwrap()])
        .output()
        .expect("spawn bwa index");
    assert!(
        output.status.success(),
        "bwa index failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn samtools_dict(ref_path: &Path) {
    let dict_path = ref_path.with_extension("dict");
    let output = Command::new("samtools")
        .args(["dict", ref_path.to_str().unwrap(), "-o", dict_path.to_str().unwrap()])
        .output()
        .expect("spawn samtools dict");
    assert!(
        output.status.success(),
        "samtools dict failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Generate paired FASTQs with an 8-base UMI prefix on both R1 and R2.
///
/// The fixture emits pairs in groups of four per starting position: two "A"
/// pairs (UMI `a-b`, R1 forward from `start`) and two "B" pairs (UMI `b-a`,
/// i.e. the same two UMI halves swapped, R2 forward from `start`). After
/// extraction, `--group::strategy paired` canonicalizes `a-b` and `b-a` into
/// the same duplex molecule with `/A` and `/B` strand suffixes on the MI tag,
/// giving the duplex consensus caller two reads per strand to work with.
fn write_paired_fastq(r1: &Path, r2: &Path, num_positions: usize) {
    use flate2::{Compression, write::GzEncoder};

    let ref_seq = ref_sequence();
    let f1 = std::fs::File::create(r1).unwrap();
    let mut e1 = GzEncoder::new(f1, Compression::default());
    let f2 = std::fs::File::create(r2).unwrap();
    let mut e2 = GzEncoder::new(f2, Compression::default());

    let revcomp = |seq: &[u8]| -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b,
            })
            .collect()
    };

    let mut read_id: usize = 0;
    for pos_idx in 0..num_positions {
        // Two distinct 8-mers used as the paired UMI halves.
        let umi_a: Vec<u8> = (0..8).map(|k| DNA[(pos_idx + k) % 4]).collect();
        let umi_b: Vec<u8> = (0..8).map(|k| DNA[(pos_idx + k + 2) % 4]).collect();

        let start = (pos_idx * 7) % 900;
        let fwd: Vec<u8> = ref_seq[start..start + 50].to_vec();
        let rev: Vec<u8> = revcomp(&ref_seq[start + 20..start + 70]);

        // Strand-A fragments: R1 forward with UMI `a`, R2 reverse with UMI `b`.
        // Two pairs with the same UMI so the /A single-strand consensus has
        // depth > 1.
        for _ in 0..2 {
            let r1_seq: Vec<u8> = [umi_a.as_slice(), fwd.as_slice()].concat();
            let r2_seq: Vec<u8> = [umi_b.as_slice(), rev.as_slice()].concat();
            emit_pair(&mut e1, &mut e2, read_id, &r1_seq, &r2_seq);
            read_id += 1;
        }

        // Strand-B fragments: swapped halves, R1/R2 sequences flipped so
        // R1 reads the opposite strand. UMI pair `b-a` pairs with `a-b`
        // under the paired assigner.
        for _ in 0..2 {
            let r1_seq: Vec<u8> = [umi_b.as_slice(), rev.as_slice()].concat();
            let r2_seq: Vec<u8> = [umi_a.as_slice(), fwd.as_slice()].concat();
            emit_pair(&mut e1, &mut e2, read_id, &r1_seq, &r2_seq);
            read_id += 1;
        }
    }
    e1.finish().unwrap();
    e2.finish().unwrap();
}

fn emit_pair<W: Write>(e1: &mut W, e2: &mut W, id: usize, r1_seq: &[u8], r2_seq: &[u8]) {
    writeln!(e1, "@read{id}/1").unwrap();
    e1.write_all(r1_seq).unwrap();
    writeln!(e1).unwrap();
    writeln!(e1, "+").unwrap();
    writeln!(e1, "{}", "I".repeat(r1_seq.len())).unwrap();

    writeln!(e2, "@read{id}/2").unwrap();
    e2.write_all(r2_seq).unwrap();
    writeln!(e2).unwrap();
    writeln!(e2, "+").unwrap();
    writeln!(e2, "{}", "I".repeat(r2_seq.len())).unwrap();
}

/// Run the standalone chain: extract -> fastq -> bwa mem -> zipper -> sort ->
/// group (paired) -> duplex -> filter.
#[allow(
    clippy::too_many_lines,
    reason = "linear chain of subprocess invocations; readability > splitting"
)]
fn run_v1_pipeline(r1: &Path, r2: &Path, reference: &Path, out: &Path, tmp_dir: &Path) {
    let unmapped = tmp_dir.join("unmapped.bam");
    let fastq = tmp_dir.join("interleaved.fastq");
    let mapped_sam = tmp_dir.join("mapped.sam");
    let zippered = tmp_dir.join("zippered.bam");
    let sorted = tmp_dir.join("sorted.bam");
    let grouped = tmp_dir.join("grouped.bam");
    let consensus = tmp_dir.join("duplex.bam");

    let output = Command::new(fgumi_bin())
        .args([
            "extract",
            "--inputs",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            unmapped.to_str().unwrap(),
            "--sample",
            "S",
            "--library",
            "L",
            "--read-structures",
            "8M+T",
            "8M+T",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi extract");
    assert!(
        output.status.success(),
        "v1 fgumi extract failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let fq_output = Command::new(fgumi_bin())
        .args(["fastq", "--input", unmapped.to_str().unwrap()])
        .output()
        .expect("spawn fgumi fastq");
    assert!(
        fq_output.status.success(),
        "v1 fgumi fastq failed:\n{}",
        String::from_utf8_lossy(&fq_output.stderr)
    );
    std::fs::write(&fastq, &fq_output.stdout).unwrap();

    let bwa_output = Command::new("bwa")
        .args([
            "mem",
            "-p",
            "-K",
            "150000000",
            "-t",
            "1",
            reference.to_str().unwrap(),
            fastq.to_str().unwrap(),
        ])
        .output()
        .expect("spawn bwa mem");
    assert!(
        bwa_output.status.success(),
        "bwa mem failed:\n{}",
        String::from_utf8_lossy(&bwa_output.stderr)
    );
    std::fs::write(&mapped_sam, &bwa_output.stdout).unwrap();

    let zip_output = Command::new(fgumi_bin())
        .args([
            "zipper",
            "--unmapped",
            unmapped.to_str().unwrap(),
            "--input",
            mapped_sam.to_str().unwrap(),
            "--output",
            zippered.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi zipper");
    assert!(
        zip_output.status.success(),
        "v1 fgumi zipper failed:\n{}",
        String::from_utf8_lossy(&zip_output.stderr)
    );

    let sort_output = Command::new(fgumi_bin())
        .args([
            "sort",
            "--input",
            zippered.to_str().unwrap(),
            "--output",
            sorted.to_str().unwrap(),
            "--order",
            "template-coordinate",
            "--threads",
            "1",
        ])
        .output()
        .expect("spawn fgumi sort");
    assert!(
        sort_output.status.success(),
        "v1 fgumi sort failed:\n{}",
        String::from_utf8_lossy(&sort_output.stderr)
    );

    let group_output = Command::new(fgumi_bin())
        .args([
            "group",
            "--input",
            sorted.to_str().unwrap(),
            "--output",
            grouped.to_str().unwrap(),
            "--strategy",
            "paired",
        ])
        .output()
        .expect("spawn fgumi group");
    assert!(
        group_output.status.success(),
        "v1 fgumi group failed:\n{}",
        String::from_utf8_lossy(&group_output.stderr)
    );

    let duplex_output = Command::new(fgumi_bin())
        .args([
            "duplex",
            "--input",
            grouped.to_str().unwrap(),
            "--output",
            consensus.to_str().unwrap(),
            "--min-reads",
            "1",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi duplex");
    assert!(
        duplex_output.status.success(),
        "v1 fgumi duplex failed:\n{}",
        String::from_utf8_lossy(&duplex_output.stderr)
    );

    let filter_output = Command::new(fgumi_bin())
        .args([
            "filter",
            "--input",
            consensus.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--min-reads",
            "1",
            "--min-base-quality",
            "2",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi filter");
    assert!(
        filter_output.status.success(),
        "v1 fgumi filter failed:\n{}",
        String::from_utf8_lossy(&filter_output.stderr)
    );
}

fn run_v2_runall(r1: &Path, r2: &Path, reference: &Path, out: &Path) {
    let output = Command::new(fgumi_bin())
        .args([
            "runall",
            "--start-from",
            "extract",
            "--input",
            r1.to_str().unwrap(),
            r2.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--extract::sample",
            "S",
            "--extract::library",
            "L",
            "--extract::read-structures",
            "8M+T",
            "8M+T",
            "--reference",
            reference.to_str().unwrap(),
            "--aligner::preset",
            "bwa-mem",
            "--aligner::threads",
            "1",
            "--consensus",
            "duplex",
            "--group::strategy",
            "paired",
            "--duplex::min-reads",
            "1",
            "--filter::min-reads",
            "1",
            "--filter::min-base-quality",
            "2",
            "--threads",
            "1",
            "--compression-level",
            "1",
        ])
        .output()
        .expect("spawn fgumi runall");
    assert!(
        output.status.success(),
        "v2 fgumi runall failed (exit {}):\nstdout:\n{}\nstderr:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

fn compare_bams(a: &Path, b: &Path) -> bool {
    let output = Command::new(fgumi_bin())
        .args(["compare", "bams", a.to_str().unwrap(), b.to_str().unwrap(), "--command", "filter"])
        .output()
        .expect("spawn fgumi compare bams");
    if !output.status.success() {
        eprintln!("fgumi compare bams stdout:\n{}", String::from_utf8_lossy(&output.stdout));
        eprintln!("fgumi compare bams stderr:\n{}", String::from_utf8_lossy(&output.stderr));
    }
    output.status.success()
}

#[test]
fn test_extract_to_filter_duplex_v2_matches_v1() {
    if !crate::helpers::references::tool_on_path("bwa")
        || !crate::helpers::references::tool_on_path("samtools")
    {
        crate::helpers::references::skip_or_fail(
            "bwa/samtools not on PATH; skipping consensus-duplex byte-identity gate",
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let reference = tmp.path().join("ref.fa");
    let r1 = tmp.path().join("R1.fastq.gz");
    let r2 = tmp.path().join("R2.fastq.gz");
    let v1_out = tmp.path().join("v1.bam");
    let v2_out = tmp.path().join("v2.bam");
    let v1_tmp = tmp.path().join("v1_tmp");
    std::fs::create_dir_all(&v1_tmp).unwrap();

    write_reference(&reference);
    bwa_index(&reference);
    samtools_dict(&reference);
    // 20 starting positions x 4 pairs each (2 per strand) = 80 pairs.
    write_paired_fastq(&r1, &r2, 20);

    run_v1_pipeline(&r1, &r2, &reference, &v1_out, &v1_tmp);
    run_v2_runall(&r1, &r2, &reference, &v2_out);

    assert!(
        compare_bams(&v1_out, &v2_out),
        "runall --start-from extract --consensus duplex output must match \
         v1 extract+fastq+bwa+zipper+sort+group(paired)+duplex+filter record-for-record"
    );
}
