//! Integration tests for the `--async-reader` flag across all supported input paths.
//!
//! Each test runs a command twice — once without `--async-reader` and once with it —
//! then asserts that the outputs are identical record-by-record. This ensures the
//! prefetch reader is transparent across:
//!
//! - **Extract:** BGZF, gzip, and plain FASTQ inputs (single- and multi-threaded)
//! - **BAM commands:** file input with/without `--threads`, and piped stdin

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use flate2::Compression;
use flate2::write::GzEncoder;
use noodles::bam;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles_bgzf::io::Writer as BgzfWriter;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, create_umi_family, to_record_buf};

// ============================================================================
// Helpers
// ============================================================================

/// Type alias for FASTQ test records (name, sequence, quality).
type FastqRecords = Vec<(&'static str, &'static str, &'static str)>;

/// Standard paired-end test records with a 5-base UMI prefix on each read.
fn paired_end_records() -> (FastqRecords, FastqRecords) {
    let r1 = vec![
        ("read1", "AAAAACGTACGTAAAA", "IIIIIIIIIIIIIIII"),
        ("read2", "TTTTTCCCCCCCCCCC", "IIIIIIIIIIIIIIII"),
        ("read3", "CCCCCTTTTTAAAAAG", "IIIIIIIIIIIIIIII"),
    ];
    let r2 = vec![
        ("read1", "GGGGGCGTACGTCCCC", "IIIIIIIIIIIIIIII"),
        ("read2", "AAAAATTTTTTTTGGG", "IIIIIIIIIIIIIIII"),
        ("read3", "TTTTTCCCCGGGGAAA", "IIIIIIIIIIIIIIII"),
    ];
    (r1, r2)
}

fn create_plain_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.path().join(name);
    let mut file = File::create(&path).unwrap();
    for (name, seq, qual) in records {
        writeln!(file, "@{name}").unwrap();
        writeln!(file, "{seq}").unwrap();
        writeln!(file, "+").unwrap();
        writeln!(file, "{qual}").unwrap();
    }
    path
}

fn create_gzip_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.path().join(name);
    let file = File::create(&path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    for (name, seq, qual) in records {
        writeln!(encoder, "@{name}").unwrap();
        writeln!(encoder, "{seq}").unwrap();
        writeln!(encoder, "+").unwrap();
        writeln!(encoder, "{qual}").unwrap();
    }
    encoder.finish().unwrap();
    path
}

fn create_bgzf_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
    let path = dir.path().join(name);
    let file = File::create(&path).unwrap();
    let mut writer = BgzfWriter::new(file);
    for (name, seq, qual) in records {
        writeln!(writer, "@{name}").unwrap();
        writeln!(writer, "{seq}").unwrap();
        writeln!(writer, "+").unwrap();
        writeln!(writer, "{qual}").unwrap();
    }
    writer.finish().unwrap();
    path
}

/// Read BAM records from a file.
fn read_bam_records(path: &Path) -> Vec<RecordBuf> {
    let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
    let header = reader.read_header().unwrap();
    reader.record_bufs(&header).map(|r| r.expect("Failed to read BAM record")).collect()
}

/// Read the BAM header from a file.
fn read_bam_header(path: &Path) -> noodles::sam::Header {
    let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
    reader.read_header().unwrap()
}

/// Assert that two BAM files contain identical records and identical headers
/// (ignoring `@PG` entries, which differ by command line).
fn assert_bam_records_equal(path1: &Path, path2: &Path) {
    // Compare headers with @PG entries stripped: those differ by command line
    // (--async-reader is in one and not the other) but everything else —
    // @HD, @SQ, @RG, @CO, sort order — must match.
    let mut header1 = read_bam_header(path1);
    let mut header2 = read_bam_header(path2);
    header1.programs_mut().as_mut().clear();
    header2.programs_mut().as_mut().clear();
    assert_eq!(header1, header2, "BAM headers differ (excluding @PG)");

    let records1 = read_bam_records(path1);
    let records2 = read_bam_records(path2);

    assert_eq!(
        records1.len(),
        records2.len(),
        "Record count differs: {} (without --async-reader) vs {} (with --async-reader)",
        records1.len(),
        records2.len(),
    );

    for (i, (r1, r2)) in records1.iter().zip(records2.iter()).enumerate() {
        assert_eq!(r1.name(), r2.name(), "Record {i} name mismatch");
        assert_eq!(r1.flags(), r2.flags(), "Record {i} flags mismatch");
        assert_eq!(
            r1.alignment_start(),
            r2.alignment_start(),
            "Record {i} alignment_start mismatch"
        );
        assert_eq!(
            r1.mapping_quality(),
            r2.mapping_quality(),
            "Record {i} mapping_quality mismatch"
        );
        assert_eq!(r1.cigar(), r2.cigar(), "Record {i} cigar mismatch");
        assert_eq!(r1.sequence(), r2.sequence(), "Record {i} sequence mismatch");
        assert_eq!(r1.quality_scores(), r2.quality_scores(), "Record {i} quality mismatch");
        assert_eq!(r1.data(), r2.data(), "Record {i} aux tags mismatch");
    }
}

/// Run `fgumi extract` with the given args, return the output BAM path.
fn run_extract(
    tmp: &TempDir,
    r1: &Path,
    r2: &Path,
    output_name: &str,
    threads: Option<usize>,
    async_reader: bool,
) -> PathBuf {
    let output = tmp.path().join(output_name);
    let mut args = vec![
        "extract".to_string(),
        "--inputs".into(),
        r1.to_str().unwrap().into(),
        r2.to_str().unwrap().into(),
        "--output".into(),
        output.to_str().unwrap().into(),
        "--read-structures".into(),
        "5M+T".into(),
        "5M+T".into(),
        "--sample".into(),
        "test_sample".into(),
        "--library".into(),
        "test_library".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
    }
    if async_reader {
        args.push("--async-reader".into());
    }

    let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(&args)
        .status()
        .expect("Failed to execute extract command");
    assert!(status.success(), "Extract command failed (async_reader={async_reader})");
    output
}

/// Create a test BAM with UMI families for grouping.
fn create_test_bam(path: &Path) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    let family1 = create_umi_family("AAAAAAAA", 5, "family1", "ACGTACGT", 30);
    let family2 = create_umi_family("CCCCCCCC", 3, "family2", "TGCATGCA", 30);
    let family3 = create_umi_family("GGGGGGGG", 4, "family3", "TTTTAAAA", 30);

    for record in family1.iter().chain(family2.iter()).chain(family3.iter()) {
        writer
            .write_alignment_record(&header, &to_record_buf(record))
            .expect("Failed to write record");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Create a test BAM with grouped reads (MI tag set) for simplex.
fn create_grouped_bam(path: &Path) {
    let header = create_minimal_header("chr1", 10000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");

    let mi_tag = Tag::from(fgumi_lib::sam::SamTag::MI);

    for (mi, umi, depth, base) in
        [(1, "AAAAAAAA", 5, "mol1"), (2, "CCCCCCCC", 3, "mol2"), (3, "GGGGGGGG", 4, "mol3")]
    {
        let records = create_umi_family(umi, depth, base, "ACGTACGT", 30);
        for record in records {
            let mut rec = to_record_buf(&record);
            rec.data_mut().insert(mi_tag, Value::from(mi));
            writer.write_alignment_record(&header, &rec).expect("Failed to write record");
        }
    }

    writer.try_finish().expect("Failed to finish BAM");
}

// ============================================================================
// Extract: BGZF FASTQ inputs
// ============================================================================

#[test]
fn test_extract_bgzf_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_recs);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", None, false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", None, true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_extract_bgzf_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_bgzf_fastq(&tmp, "r1.fq.bgz", &r1_recs);
    let r2 = create_bgzf_fastq(&tmp, "r2.fq.bgz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", Some(4), false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", Some(4), true);
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// Extract: gzip FASTQ inputs
// ============================================================================

#[test]
fn test_extract_gzip_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_recs);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", None, false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", None, true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_extract_gzip_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_gzip_fastq(&tmp, "r1.fq.gz", &r1_recs);
    let r2 = create_gzip_fastq(&tmp, "r2.fq.gz", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", Some(4), false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", Some(4), true);
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// Extract: plain (uncompressed) FASTQ inputs
// ============================================================================

#[test]
fn test_extract_plain_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_plain_fastq(&tmp, "r1.fq", &r1_recs);
    let r2 = create_plain_fastq(&tmp, "r2.fq", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", None, false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", None, true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_extract_plain_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let (r1_recs, r2_recs) = paired_end_records();
    let r1 = create_plain_fastq(&tmp, "r1.fq", &r1_recs);
    let r2 = create_plain_fastq(&tmp, "r2.fq", &r2_recs);

    let baseline = run_extract(&tmp, &r1, &r2, "baseline.bam", Some(4), false);
    let async_out = run_extract(&tmp, &r1, &r2, "async.bam", Some(4), true);
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// BAM commands: group (file input, with and without --threads)
// ============================================================================

/// Run `fgumi group` and return the output BAM path.
fn run_group(
    tmp: &TempDir,
    input: &str,
    output_name: &str,
    threads: Option<usize>,
    async_reader: bool,
    stdin: Option<Stdio>,
) -> PathBuf {
    let output = tmp.path().join(output_name);
    let mut args = vec![
        "group".to_string(),
        "--input".into(),
        input.into(),
        "--output".into(),
        output.to_str().unwrap().into(),
        "--strategy".into(),
        "identity".into(),
        "--edits".into(),
        "0".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
    }
    if async_reader {
        args.push("--async-reader".into());
    }

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.args(&args);
    if let Some(stdin) = stdin {
        cmd.stdin(stdin);
    }
    let status = cmd.status().expect("Failed to execute group command");
    assert!(status.success(), "Group command failed (async_reader={async_reader})");
    output
}

#[test]
fn test_group_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    create_test_bam(&input);

    let baseline = run_group(&tmp, input.to_str().unwrap(), "baseline.bam", Some(2), false, None);
    let async_out = run_group(&tmp, input.to_str().unwrap(), "async.bam", Some(2), true, None);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_group_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    create_test_bam(&input);

    let baseline = run_group(&tmp, input.to_str().unwrap(), "baseline.bam", None, false, None);
    let async_out = run_group(&tmp, input.to_str().unwrap(), "async.bam", None, true, None);
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// BAM commands: group with stdin input
// ============================================================================

#[test]
fn test_group_async_reader_stdin() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    create_test_bam(&input);

    // Baseline: stdin without --async-reader
    let baseline = run_group(
        &tmp,
        "-",
        "baseline.bam",
        Some(2),
        false,
        Some(Stdio::from(File::open(&input).unwrap())),
    );

    // Async: stdin with --async-reader
    let async_out = run_group(
        &tmp,
        "-",
        "async_stdin.bam",
        Some(2),
        true,
        Some(Stdio::from(File::open(&input).unwrap())),
    );
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// BAM commands: simplex (file input, with and without --threads)
// ============================================================================

/// Run `fgumi simplex` and return the output BAM path.
fn run_simplex(
    tmp: &TempDir,
    input: &str,
    output_name: &str,
    threads: Option<usize>,
    async_reader: bool,
    stdin: Option<Stdio>,
) -> PathBuf {
    let output = tmp.path().join(output_name);
    let mut args = vec![
        "simplex".to_string(),
        "--input".into(),
        input.into(),
        "--output".into(),
        output.to_str().unwrap().into(),
        "--min-reads".into(),
        "1".into(),
        "--min-consensus-base-quality".into(),
        "0".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
    }
    if async_reader {
        args.push("--async-reader".into());
    }

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.args(&args);
    if let Some(stdin) = stdin {
        cmd.stdin(stdin);
    }
    let status = cmd.status().expect("Failed to execute simplex command");
    assert!(status.success(), "Simplex command failed (async_reader={async_reader})");
    output
}

#[test]
fn test_simplex_async_reader_multithreaded() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    create_grouped_bam(&input);

    let baseline = run_simplex(&tmp, input.to_str().unwrap(), "baseline.bam", Some(2), false, None);
    let async_out = run_simplex(&tmp, input.to_str().unwrap(), "async.bam", Some(2), true, None);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_simplex_async_reader_single_threaded() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    create_grouped_bam(&input);

    let baseline = run_simplex(&tmp, input.to_str().unwrap(), "baseline.bam", None, false, None);
    let async_out = run_simplex(&tmp, input.to_str().unwrap(), "async.bam", None, true, None);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
fn test_simplex_async_reader_stdin() {
    let tmp = TempDir::new().unwrap();
    let input = tmp.path().join("input.bam");
    create_grouped_bam(&input);

    // Baseline: stdin without --async-reader
    let baseline = run_simplex(
        &tmp,
        "-",
        "baseline.bam",
        Some(2),
        false,
        Some(Stdio::from(File::open(&input).unwrap())),
    );

    // Async: stdin with --async-reader
    let async_out = run_simplex(
        &tmp,
        "-",
        "async_stdin.bam",
        Some(2),
        true,
        Some(Stdio::from(File::open(&input).unwrap())),
    );
    assert_bam_records_equal(&baseline, &async_out);
}

// ============================================================================
// Sort: async-reader regression tests
// ============================================================================

/// Check if samtools is available for sort tests.
fn samtools_available() -> bool {
    Command::new("samtools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

/// Create an unsorted BAM with many records spread across chromosomes using samtools.
fn create_unsorted_bam(dir: &Path, num_reads: usize) -> PathBuf {
    use std::fmt::Write as _;

    let bam_path = dir.join("unsorted.bam");
    let mut sam_content = String::new();
    sam_content.push_str("@HD\tVN:1.6\tSO:unsorted\n");
    sam_content.push_str("@SQ\tSN:chr1\tLN:100000\n");
    sam_content.push_str("@SQ\tSN:chr2\tLN:100000\n");
    sam_content.push_str("@SQ\tSN:chr3\tLN:100000\n");
    sam_content.push_str("@RG\tID:rg1\tSM:sample1\tLB:lib1\n");

    let seq = "ACGTACGTAC";
    let qual = "IIIIIIIIII";

    for i in 0..num_reads {
        let chrom = match i % 3 {
            0 => "chr2",
            1 => "chr1",
            _ => "chr3",
        };
        let pos = 50000 - (i % 1000) * 50 + 1;
        let name = format!("read_{i}");
        let mate_pos = pos + 200;
        let cell = format!("cell{}", i % 4);

        // R1
        writeln!(
            sam_content,
            "{name}\t99\t{chrom}\t{pos}\t60\t10M\t=\t{mate_pos}\t210\t{seq}\t{qual}\tRG:Z:rg1\tMI:Z:umi_{i}\tCB:Z:{cell}"
        )
        .unwrap();
        // R2
        writeln!(
            sam_content,
            "{name}\t147\t{chrom}\t{mate_pos}\t60\t10M\t=\t{pos}\t-210\t{seq}\t{qual}\tRG:Z:rg1\tMI:Z:umi_{i}\tCB:Z:{cell}"
        )
        .unwrap();
    }

    // Add unmapped reads
    for i in 0..10 {
        let name = format!("unmapped_{i:04}");
        writeln!(sam_content, "{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1").unwrap();
        writeln!(sam_content, "{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:rg1").unwrap();
    }

    let sam_path = dir.join("unsorted.sam");
    std::fs::write(&sam_path, sam_content).expect("write SAM");

    let status = Command::new("samtools")
        .args(["view", "-b", "-o", bam_path.to_str().unwrap(), sam_path.to_str().unwrap()])
        .status()
        .expect("samtools view");
    assert!(status.success(), "samtools view failed");
    bam_path
}

/// Run `fgumi sort` and return the output BAM path.
fn run_sort(
    tmp: &TempDir,
    input: &Path,
    output_name: &str,
    order: &str,
    threads: usize,
    max_memory: &str,
    async_reader: bool,
) -> PathBuf {
    let output = tmp.path().join(output_name);
    let mut args: Vec<std::ffi::OsString> = vec![
        "sort".into(),
        "-i".into(),
        input.as_os_str().to_owned(),
        "-o".into(),
        output.as_os_str().to_owned(),
        "--order".into(),
        order.into(),
        "--threads".into(),
        threads.to_string().into(),
        "-m".into(),
        max_memory.into(),
        "--temp-compression".into(),
        "1".into(),
        "--compression-level".into(),
        "1".into(),
    ];
    if async_reader {
        args.push("--async-reader".into());
    }

    let result =
        Command::new(env!("CARGO_BIN_EXE_fgumi")).args(&args).output().expect("fgumi sort");
    assert!(
        result.status.success(),
        "fgumi sort failed (order={order}, threads={threads}, async={async_reader}):\n{}",
        String::from_utf8_lossy(&result.stderr),
    );
    output
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_async_reader_single_threaded() {
    if !samtools_available() {
        return;
    }
    let tmp = TempDir::new().unwrap();
    let input = create_unsorted_bam(tmp.path(), 500);

    let baseline = run_sort(&tmp, &input, "baseline.bam", "coordinate", 1, "100M", false);
    let async_out = run_sort(&tmp, &input, "async.bam", "coordinate", 1, "100M", true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_coordinate_async_reader_multithreaded() {
    if !samtools_available() {
        return;
    }
    let tmp = TempDir::new().unwrap();
    let input = create_unsorted_bam(tmp.path(), 2000);

    let baseline = run_sort(&tmp, &input, "baseline.bam", "coordinate", 4, "50K", false);
    let async_out = run_sort(&tmp, &input, "async.bam", "coordinate", 4, "50K", true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_queryname_async_reader_multithreaded() {
    if !samtools_available() {
        return;
    }
    let tmp = TempDir::new().unwrap();
    let input = create_unsorted_bam(tmp.path(), 2000);

    let baseline = run_sort(&tmp, &input, "baseline.bam", "queryname", 4, "50K", false);
    let async_out = run_sort(&tmp, &input, "async.bam", "queryname", 4, "50K", true);
    assert_bam_records_equal(&baseline, &async_out);
}

#[test]
#[ignore = "requires samtools"]
fn test_sort_template_coordinate_async_reader_multithreaded() {
    if !samtools_available() {
        return;
    }
    let tmp = TempDir::new().unwrap();
    let input = create_unsorted_bam(tmp.path(), 2000);

    let baseline = run_sort(&tmp, &input, "baseline.bam", "template-coordinate", 4, "50K", false);
    let async_out = run_sort(&tmp, &input, "async.bam", "template-coordinate", 4, "50K", true);
    assert_bam_records_equal(&baseline, &async_out);
}
