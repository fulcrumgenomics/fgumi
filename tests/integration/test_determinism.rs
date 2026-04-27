//! `runall` produces byte-identical output across independent runs with the
//! same input. Regression guard for non-deterministic iteration order (e.g.
//! `HashMap`) leaking into output.
//!
//! All runs use `-t 1` (and `--aligner::threads 1` where applicable) so that
//! only the pipeline's own ordering is being tested — not `bwa`'s thread
//! non-determinism.

use std::io::Write;
use std::path::Path;
use std::process::Command;

use crate::helpers::grouped_bam_fixture::build_small_grouped_bam;
use crate::helpers::references::{skip_or_fail, tool_on_path};

fn fgumi_bin() -> &'static str {
    env!("CARGO_BIN_EXE_fgumi")
}

/// Run runall twice with identical args and return the two output payloads.
///
/// Both runs use the same output filename so that the `CL:` field of the
/// `@PG` header line (which includes the command line verbatim) is identical
/// across runs. The first run's output is moved aside before the second run.
fn run_twice_same_args(
    args: &[&std::ffi::OsStr],
    work_dir: &Path,
    output_name: &str,
) -> (Vec<u8>, Vec<u8>) {
    let out_path = work_dir.join(output_name);

    let run_once = || {
        let status = Command::new(fgumi_bin())
            .args(args)
            .current_dir(work_dir)
            .status()
            .expect("spawn fgumi runall");
        assert!(status.success(), "fgumi runall failed with exit {status}");
        std::fs::read(&out_path).unwrap()
    };

    let a = run_once();
    std::fs::remove_file(&out_path).unwrap();
    let b = run_once();
    (a, b)
}

#[test]
fn runall_group_to_filter_is_deterministic() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("input.bam");
    build_small_grouped_bam(&input);

    // Use relative paths via `current_dir` so the @PG CL field is identical
    // across runs (absolute tmpdir paths are injected there too).
    let args: Vec<std::ffi::OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "group".into(),
        "--threads".into(),
        "1".into(),
        "--input".into(),
        "input.bam".into(),
        "--output".into(),
        "out.bam".into(),
        // --filter::min-reads is required when the plan runs filter (default
        // --stop-after); match standalone fgumi filter semantics.
        "--filter::min-reads".into(),
        "1".into(),
    ];
    let arg_refs: Vec<&std::ffi::OsStr> = args.iter().map(AsRef::as_ref).collect();

    let (a, b) = run_twice_same_args(&arg_refs, tmp.path(), "out.bam");
    assert_eq!(
        a,
        b,
        "runall --start-from group produced different bytes on re-run \
         (len_a={}, len_b={})",
        a.len(),
        b.len()
    );
}

// ----------------------------------------------------------------------------
// extract-to-sort determinism (requires bwa + samtools)
// ----------------------------------------------------------------------------
//
// The helpers below mirror those in `tests/integration/runall/test_sort.rs`
// (file-local there). We duplicate rather than lift out to keep this commit
// small and focused — the duplication is confined to this single test.

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

fn write_paired_fastq(r1: &Path, r2: &Path, num_pairs: usize) {
    use flate2::{Compression, write::GzEncoder};

    let ref_seq = ref_sequence();
    let f1 = std::fs::File::create(r1).unwrap();
    let mut e1 = GzEncoder::new(f1, Compression::default());
    let f2 = std::fs::File::create(r2).unwrap();
    let mut e2 = GzEncoder::new(f2, Compression::default());

    for i in 0..num_pairs {
        let umi: Vec<u8> = (0..8).map(|k| DNA[(i + k) % 4]).collect();
        let start = (i * 7) % 900;
        let r1_template: Vec<u8> = ref_seq[start..start + 50].to_vec();
        let r2_template: Vec<u8> = ref_seq[start + 20..start + 70]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b,
            })
            .collect();

        let r1_seq: Vec<u8> = [umi.as_slice(), r1_template.as_slice()].concat();

        writeln!(e1, "@read{i}/1").unwrap();
        e1.write_all(&r1_seq).unwrap();
        writeln!(e1).unwrap();
        writeln!(e1, "+").unwrap();
        writeln!(e1, "{}", "I".repeat(r1_seq.len())).unwrap();

        writeln!(e2, "@read{i}/2").unwrap();
        e2.write_all(&r2_template).unwrap();
        writeln!(e2).unwrap();
        writeln!(e2, "+").unwrap();
        writeln!(e2, "{}", "I".repeat(r2_template.len())).unwrap();
    }
    e1.finish().unwrap();
    e2.finish().unwrap();
}

#[test]
fn runall_extract_to_sort_is_deterministic() {
    if !tool_on_path("bwa") || !tool_on_path("samtools") {
        skip_or_fail("bwa or samtools not on PATH; skipping extract-to-sort determinism test");
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let reference = tmp.path().join("ref.fa");
    let r1 = tmp.path().join("R1.fastq.gz");
    let r2 = tmp.path().join("R2.fastq.gz");

    write_reference(&reference);
    bwa_index(&reference);
    samtools_dict(&reference);
    write_paired_fastq(&r1, &r2, 50);

    // Use relative paths via `current_dir` so the @PG CL field is identical
    // across runs.
    let args: Vec<std::ffi::OsString> = vec![
        "runall".into(),
        "--start-from".into(),
        "extract".into(),
        "--stop-after".into(),
        "sort".into(),
        "--input".into(),
        "R1.fastq.gz".into(),
        "R2.fastq.gz".into(),
        "--output".into(),
        "out.bam".into(),
        "--extract::sample".into(),
        "S".into(),
        "--extract::library".into(),
        "L".into(),
        "--extract::read-structures".into(),
        "8M+T".into(),
        "+T".into(),
        "--reference".into(),
        "ref.fa".into(),
        "--aligner::preset".into(),
        "bwa-mem".into(),
        "--aligner::threads".into(),
        "1".into(),
        "--threads".into(),
        "1".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    let arg_refs: Vec<&std::ffi::OsStr> = args.iter().map(AsRef::as_ref).collect();

    let (a, b) = run_twice_same_args(&arg_refs, tmp.path(), "out.bam");
    assert_eq!(
        a,
        b,
        "runall --start-from extract --stop-after sort produced different bytes on re-run \
         (len_a={}, len_b={})",
        a.len(),
        b.len()
    );
}
