//! Every command that accepts BAM must accept uncompressed SAM.
//!
//! fgumi opens input through several different readers — the single-threaded
//! raw-byte reader, the multi-threaded pipeline reader, and the sort engine's
//! pool-integrated reader — and a command picks between them based on flags
//! like `--threads` that have nothing to do with the input format. That makes
//! "does this command accept SAM?" easy to answer differently per command, and
//! even per flag combination within one command.
//!
//! These tests pin the invariant instead: for each command, the same data as
//! `.bam` and as `.sam` must both succeed and produce the same output — header
//! and records alike.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_minimal_header, create_umi_family, transcode_bam_to_sam, write_bam,
};
use crate::helpers::read_bam_output as read_output;

/// Writes a small grouped-input BAM and the byte-for-byte equivalent SAM.
///
/// Both files hold the same records, so any difference in a command's output
/// between the two is a difference in how the input was read.
fn write_input_pair(dir: &Path) -> (PathBuf, PathBuf) {
    let bam_path = dir.join("input.bam");
    let sam_path = dir.join("input.sam");

    // `create_minimal_header` already advertises SO:unsorted GO:query
    // SS:template-coordinate, which is what the grouping commands require.
    let header = create_minimal_header("chr1", 10_000);

    let mut records = create_umi_family("ACGTACGT", 3, "fam_a", "ACGTACGTAC", 35);
    records.extend(create_umi_family("TTTTGGGG", 2, "fam_b", "ACGTACGTAC", 35));
    write_bam(&bam_path, &header, &records);

    transcode_bam_to_sam(&bam_path, &sam_path);

    (bam_path, sam_path)
}

/// Runs a command with `input` substituted for the `{input}` placeholder.
fn run(args: &[&str], input: &Path, output: &Path) -> std::process::Output {
    let input = input.to_str().expect("input path is UTF-8");
    let output = output.to_str().expect("output path is UTF-8");
    let args: Vec<&str> = args
        .iter()
        .map(|arg| match *arg {
            "{input}" => input,
            "{output}" => output,
            other => other,
        })
        .collect();

    Command::new(env!("CARGO_BIN_EXE_fgumi")).args(&args).output().expect("failed to run fgumi")
}

/// Runs a command with `input` fed in on stdin rather than named as a path.
fn run_with_stdin(args: &[&str], input: &Path, output: &Path) -> std::process::Output {
    let output = output.to_str().expect("output path is UTF-8");
    let args: Vec<&str> =
        args.iter().map(|arg| if *arg == "{output}" { output } else { *arg }).collect();

    let stdin = File::open(input).expect("open input to pipe");
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(&args)
        .stdin(std::process::Stdio::from(stdin))
        .output()
        .expect("failed to run fgumi")
}

/// A command reading BAM must read the same data as SAM, and must produce the
/// same records from both.
///
/// `--threads` is varied per case because it selects between entirely separate
/// orchestrations (single-threaded fast path vs. pipeline); input format
/// support must not depend on which one runs.
#[rstest]
#[case::sort_single_threaded(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"])]
#[case::sort_multi_threaded(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname", "--threads", "2"])]
#[case::group_single_threaded(&["group", "-i", "{input}", "-o", "{output}", "--strategy", "identity", "--edits", "0"])]
#[case::group_multi_threaded(&["group", "-i", "{input}", "-o", "{output}", "--strategy", "identity", "--edits", "0", "--threads", "2"])]
#[case::dedup(&["dedup", "-i", "{input}", "-o", "{output}"])]
#[case::sort_coordinate(&["sort", "-i", "{input}", "-o", "{output}", "--order", "coordinate"])]
#[case::sort_template_coordinate(&["sort", "-i", "{input}", "-o", "{output}", "--order", "template-coordinate"])]
fn command_accepts_sam_input_identically_to_bam(#[case] args: &[&str]) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());

    let bam_output = dir.path().join("bam_output.bam");
    let bam_run = run(args, &bam_input, &bam_output);
    assert!(
        bam_run.status.success(),
        "BAM input failed: {}",
        String::from_utf8_lossy(&bam_run.stderr)
    );

    let sam_output = dir.path().join("sam_output.bam");
    let sam_run = run(args, &sam_input, &sam_output);
    assert!(
        sam_run.status.success(),
        "SAM input failed where BAM succeeded: {}",
        String::from_utf8_lossy(&sam_run.stderr)
    );

    assert_eq!(
        read_output(&sam_output),
        read_output(&bam_output),
        "SAM and BAM input produced different output",
    );
}

/// Stdin is the one input that cannot be rewound or reopened, which is why the
/// header is parsed through a tee and replayed. Piping SAM in is therefore the
/// case that a regression to seek-based header parsing would break first.
#[rstest]
#[case::sam_on_stdin("input.sam")]
#[case::bam_on_stdin("input.bam")]
fn sort_accepts_piped_input_identically_to_a_file(#[case] piped: &str) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, _) = write_input_pair(dir.path());

    let piped_output = dir.path().join("piped_output.bam");
    let piped_run = run_with_stdin(
        &["sort", "-i", "-", "-o", "{output}", "--order", "queryname"],
        &dir.path().join(piped),
        &piped_output,
    );
    assert!(
        piped_run.status.success(),
        "piped input failed: {}",
        String::from_utf8_lossy(&piped_run.stderr)
    );

    let file_output = dir.path().join("file_output.bam");
    let file_run = run(
        &["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"],
        &bam_input,
        &file_output,
    );
    assert!(
        file_run.status.success(),
        "BAM file input failed: {}",
        String::from_utf8_lossy(&file_run.stderr)
    );

    assert_eq!(
        read_output(&piped_output),
        read_output(&file_output),
        "piped input produced different output than the BAM file run",
    );
}

/// Format detection reads the input's magic bytes, so a file's extension
/// cannot make a command misread it. A `.bam` holding SAM text and a `.sam`
/// holding BGZF both have to work.
#[rstest]
#[case::sam_text_named_bam("misnamed.bam", false)]
#[case::bgzf_named_sam("misnamed.sam", true)]
fn input_format_follows_content_not_extension(#[case] name: &str, #[case] write_bgzf: bool) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());
    let args = &["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"];

    let misnamed = dir.path().join(name);
    let source = if write_bgzf { &bam_input } else { &sam_input };
    std::fs::copy(source, &misnamed).expect("copy to misnamed path");

    let output = dir.path().join("out.bam");
    let result = run(args, &misnamed, &output);
    assert!(
        result.status.success(),
        "misnamed input rejected: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    // Succeeding is not enough: the run must produce what the correctly named
    // input produces, so a misread that silently drops records fails here.
    let canonical_output = dir.path().join("canonical.bam");
    let canonical = run(args, &bam_input, &canonical_output);
    assert!(
        canonical.status.success(),
        "canonical BAM input rejected: {}",
        String::from_utf8_lossy(&canonical.stderr)
    );

    assert_eq!(
        read_output(&output),
        read_output(&canonical_output),
        "misnamed input produced different output than the canonical BAM run",
    );
}

/// `fastq` reads through the single-threaded raw-byte reader and writes FASTQ,
/// so it exercises a reader surface the BAM-output cases do not and its output
/// can be compared byte for byte.
#[test]
fn fastq_command_accepts_sam_input_identically_to_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());

    let bam_output = dir.path().join("bam_output.fq");
    let bam_run = run(&["fastq", "-i", "{input}", "-o", "{output}"], &bam_input, &bam_output);
    assert!(
        bam_run.status.success(),
        "BAM input failed: {}",
        String::from_utf8_lossy(&bam_run.stderr)
    );

    let sam_output = dir.path().join("sam_output.fq");
    let sam_run = run(&["fastq", "-i", "{input}", "-o", "{output}"], &sam_input, &sam_output);
    assert!(
        sam_run.status.success(),
        "SAM input failed where BAM succeeded: {}",
        String::from_utf8_lossy(&sam_run.stderr)
    );

    assert_eq!(
        std::fs::read(&sam_output).expect("read SAM-derived FASTQ"),
        std::fs::read(&bam_output).expect("read BAM-derived FASTQ"),
        "SAM and BAM input produced different FASTQ",
    );
}

/// `sort --verify` reads records without writing a BAM, through a reader that
/// used to be opened separately from the one that reads the header. That second
/// open bypassed format detection, so `--verify` rejected the very SAM files
/// `sort` itself accepted.
#[rstest]
#[case::coordinate("coordinate")]
#[case::queryname("queryname")]
fn sort_verify_accepts_sam_input(#[case] order: &str) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, _) = write_input_pair(dir.path());

    // `--verify` checks an already-sorted file, so produce one first.
    let sorted = dir.path().join("sorted.bam");
    let sort_run =
        run(&["sort", "-i", "{input}", "-o", "{output}", "--order", order], &bam_input, &sorted);
    assert!(
        sort_run.status.success(),
        "sort failed: {}",
        String::from_utf8_lossy(&sort_run.stderr)
    );

    let sorted_sam = dir.path().join("sorted.sam");
    transcode_bam_to_sam(&sorted, &sorted_sam);

    let summaries: Vec<(usize, usize)> = [&sorted, &sorted_sam]
        .iter()
        .map(|input| {
            let result = verify_sort_order(input, order);
            assert!(
                result.status.success(),
                "sort --verify rejected {}: {}",
                input.display(),
                String::from_utf8_lossy(&result.stderr)
            );
            verify_summary(&result)
        })
        .collect();

    // Exiting zero is not enough on its own: a `--verify` that read no records
    // at all would also exit zero. Pinning the reported record count — and
    // requiring it to be non-zero and equal across the two formats — is what
    // makes "it accepted SAM" mean "it read the SAM records".
    assert!(summaries[0].0 > 0, "verify checked no records even for BAM input");
    assert_eq!(
        summaries[1], summaries[0],
        "sort --verify reported (checked, violations) {:?} for SAM but {:?} for the equivalent BAM",
        summaries[1], summaries[0]
    );
}

/// Runs `sort --verify` over `input` in `order`.
fn verify_sort_order(input: &Path, order: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["sort", "-i", input.to_str().expect("path is UTF-8"), "--verify"])
        .args(["--order", order])
        .output()
        .expect("failed to run fgumi")
}

/// Extracts `(records checked, sort-order violations)` from a `--verify` run.
fn verify_summary(result: &std::process::Output) -> (usize, usize) {
    let stderr = String::from_utf8_lossy(&result.stderr);
    (crate::helpers::records_checked(&stderr), crate::helpers::violations(&stderr))
}

/// `sort --verify` must detect the *same violations* in SAM and BAM, not merely
/// read the same number of records.
///
/// A SAM-specific record-boundary bug could preserve the record total while
/// corrupting order comparison, and a parity test that checks only the count
/// would pass anyway. Feeding an input that is deliberately out of queryname
/// order makes the violation count the thing under test: both formats must
/// report it, and both must fail.
#[test]
fn sort_verify_reports_the_same_violations_for_sam_and_bam() {
    let dir = TempDir::new().expect("create temp dir");
    let unsorted_bam = dir.path().join("unsorted.bam");
    let equivalent_sam = dir.path().join("unsorted.sam");

    // `fam_b` before `fam_a` is descending queryname order, so every boundary
    // between the two families is a violation.
    let header = create_minimal_header("chr1", 10_000);
    let mut records = create_umi_family("TTTTGGGG", 2, "fam_b", "ACGTACGTAC", 35);
    records.extend(create_umi_family("ACGTACGT", 3, "fam_a", "ACGTACGTAC", 35));
    write_bam(&unsorted_bam, &header, &records);
    transcode_bam_to_sam(&unsorted_bam, &equivalent_sam);

    let summaries: Vec<(usize, usize)> = [&unsorted_bam, &equivalent_sam]
        .iter()
        .map(|input| {
            let result = verify_sort_order(input, "queryname");
            assert!(
                !result.status.success(),
                "sort --verify accepted out-of-order {}: {}",
                input.display(),
                String::from_utf8_lossy(&result.stderr)
            );
            verify_summary(&result)
        })
        .collect();

    assert!(summaries[0].1 > 0, "no violations reported for out-of-order BAM input");
    assert_eq!(
        summaries[1], summaries[0],
        "sort --verify reported (checked, violations) {:?} for SAM but {:?} for the equivalent BAM",
        summaries[1], summaries[0]
    );
}

/// A FIFO cannot be rewound or re-opened, so any byte read to classify the
/// input is gone for good.
///
/// `zipper` picks its mapped reader by sniffing magic bytes. Sniffing by path
/// and then letting the chosen reader re-open it works for a regular file and
/// silently corrupts a pipe — which is what `-i <(samtools view -h ...)`, a
/// routine invocation, hands it. Both formats are covered because the two
/// sniff outcomes take different branches.
#[cfg(unix)]
#[rstest]
#[case::sam_over_fifo(false)]
#[case::bam_over_fifo(true)]
fn zipper_reads_a_fifo_without_losing_the_leading_bytes(#[case] feed_bgzf: bool) {
    use crate::helpers::bam_generator::create_test_reference;

    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());
    let reference = create_test_reference(dir.path());

    let source = if feed_bgzf { &bam_input } else { &sam_input };

    // A regular-file run of the same data is the oracle: the FIFO run has to
    // match it record for record, not merely exit zero.
    let expected = dir.path().join("expected.bam");
    let baseline = run_zipper(source, &bam_input, &reference, &expected);
    assert!(
        baseline.status.success(),
        "regular-file zipper run failed: {}",
        String::from_utf8_lossy(&baseline.stderr)
    );

    let fifo = dir.path().join("mapped.fifo");
    make_fifo(&fifo);

    // The writer must run concurrently: opening a FIFO for writing blocks until
    // a reader attaches, and vice versa. That blocking open is also the hand-off
    // — opening `O_RDWR` instead would return immediately, and this thread could
    // then copy the whole (small) file into the pipe buffer and close both ends
    // before zipper ever opened it, leaving the bytes discarded and zipper
    // waiting on a writer that no longer exists.
    let source = source.clone();
    let fifo_path = fifo.clone();
    let feeder = std::thread::spawn(move || {
        let mut src = File::open(&source).expect("open FIFO source");
        let mut sink = File::create(&fifo_path).expect("open FIFO for writing");
        std::io::copy(&mut src, &mut sink).expect("feed FIFO");
    });

    let actual = dir.path().join("actual.bam");
    let result = run_zipper(&fifo, &bam_input, &reference, &actual);

    // Assert *before* joining. A zipper that dies before opening the FIFO leaves
    // the feeder blocked in its `open` forever, so joining first would turn a
    // zipper regression into a CI hang with no failing assertion to read. Failing
    // here instead leaves the blocked thread to be reaped at process exit. On the
    // success path zipper has drained the FIFO, so the join below returns at once.
    assert!(
        result.status.success(),
        "zipper rejected a FIFO it accepts as a regular file: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    feeder.join().expect("FIFO feeder thread panicked");

    assert_eq!(
        read_output(&actual),
        read_output(&expected),
        "FIFO input produced different output than the same data as a regular file",
    );
}

/// Creates a FIFO at `path`.
#[cfg(unix)]
fn make_fifo(path: &Path) {
    let status = Command::new("mkfifo").arg(path).status().expect("failed to run mkfifo");
    assert!(status.success(), "mkfifo failed for {}", path.display());
}

/// Runs `zipper` with `mapped` as the mapped input and `unmapped` as the unmapped input.
#[cfg(unix)]
fn run_zipper(
    mapped: &Path,
    unmapped: &Path,
    reference: &Path,
    output: &Path,
) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["zipper", "-i"])
        .arg(mapped)
        .arg("-u")
        .arg(unmapped)
        .arg("-r")
        .arg(reference)
        .arg("-o")
        .arg(output)
        .output()
        .expect("failed to run fgumi")
}

/// Input that is neither BGZF nor SAM must be rejected outright.
///
/// A SAM header parse consumes only `@` lines, so arbitrary text yields an
/// *empty* header rather than an error — which would transcode into a
/// valid-looking BAM carrying no reference sequences at all. Before SAM was
/// accepted these inputs failed the BGZF magic check, and they must keep
/// failing.
#[rstest]
#[case::plain_text(b"this is not a SAM file\nnor is this\n" as &[u8], "neither BAM nor SAM")]
#[case::empty(b"" as &[u8], "empty")]
fn non_sam_non_bam_input_is_rejected(#[case] contents: &[u8], #[case] expected: &str) {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("input.sam");
    std::fs::write(&input, contents).expect("write input");

    let output = dir.path().join("out.bam");
    let result =
        run(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"], &input, &output);

    assert!(!result.status.success(), "malformed input was accepted");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains(expected), "error should mention {expected:?}, got: {stderr}");
}

/// A malformed *record* must be reported as a record fault.
///
/// The transcode is lazy, so records are pulled while the consumer is still
/// reading the header; without care a bad record on line 4 surfaces as "failed
/// to read header", sending whoever reads the error to the wrong end of the
/// file.
#[test]
fn a_malformed_record_is_not_reported_as_a_header_error() {
    let dir = TempDir::new().expect("create temp dir");
    let input = dir.path().join("badrec.sam");
    std::fs::write(
        &input,
        "@HD\tVN:1.6\tSO:unsorted\n\
         @SQ\tSN:chr1\tLN:100\n\
         r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
         NOT_A_RECORD\n",
    )
    .expect("write input");

    let output = dir.path().join("out.bam");
    let result =
        run(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"], &input, &output);

    assert!(!result.status.success(), "malformed record was accepted");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("SAM record"), "error should name the record, got: {stderr}");
    assert!(
        !stderr.contains("Failed to read header"),
        "record fault misreported as a header failure: {stderr}"
    );
}

/// A header fault on zipper's stream path must name the mapped input.
///
/// `zipper` reads two inputs, and only the mapped one can arrive as a stream.
/// `normalize_to_bgzf` names the input in its own errors, but the BAM header read
/// that follows is noodles', which does not — so a fault there used to surface as a
/// bare "invalid BAM header" that read as if it came from `--unmapped`.
///
/// A BGZF stream whose contents are not BAM gets past format detection (it *is*
/// BGZF) and fails at exactly that header read, which is the fault this covers.
#[test]
fn a_stream_header_fault_names_the_mapped_input() {
    use std::io::Write as _;

    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, _) = write_input_pair(dir.path());
    let reference = crate::helpers::bam_generator::create_test_reference(dir.path());

    // Valid BGZF framing, contents that are not a BAM header.
    let bgzf = {
        let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"NOTBAM\x00\x00").expect("write BGZF payload");
        writer.finish().expect("finish BGZF stream")
    };
    let piped = dir.path().join("notbam.bgzf");
    std::fs::write(&piped, &bgzf).expect("write BGZF input");

    let output = dir.path().join("out.bam");
    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["zipper", "-i", "-"])
        .arg("-u")
        .arg(&bam_input)
        .arg("-r")
        .arg(&reference)
        .arg("-o")
        .arg(&output)
        .stdin(std::process::Stdio::from(File::open(&piped).expect("open BGZF input")))
        .output()
        .expect("failed to run fgumi");

    assert!(!result.status.success(), "a non-BAM BGZF stream was accepted");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("mapped BAM header") && stderr.contains("stdin"),
        "a header fault on the mapped stream must name it, got: {stderr}"
    );
}

// ============================================================================
// Plain gzip, accepted at the normalization boundary
// ============================================================================

/// gzip-compresses `source` into `dest`.
fn gzip_file(source: &Path, dest: &Path) {
    use std::io::Write as _;

    let bytes = std::fs::read(source).expect("read gzip source");
    let mut encoder = flate2::write::GzEncoder::new(
        File::create(dest).expect("create gzip output"),
        flate2::Compression::fast(),
    );
    encoder.write_all(&bytes).expect("gzip the input");
    encoder.finish().expect("finish gzip stream");
}

/// Rewrites every BGZF block in `source` to carry a second extra subfield.
///
/// `BC` stays first and `BSIZE` stays correct, so the result is a spec-legal BGZF
/// file per RFC 1952 — and one that no fgumi decoder can frame, because they all
/// require `XLEN` to be exactly 6 (as do htslib and noodles). It is still a valid
/// gzip member, which is the whole point: the boundary decompresses it with a
/// spec-complete gzip reader and re-frames what comes out.
fn add_trailing_extra_subfield(source: &Path, dest: &Path) {
    let data = std::fs::read(source).expect("read BGZF source");
    let mut out = Vec::with_capacity(data.len());
    let mut offset = 0;

    while offset < data.len() {
        let bsize = u16::from_le_bytes([data[offset + 16], data[offset + 17]]) as usize;
        let total = bsize + 1;
        let new_total = total + 4; // the foreign subfield's four bytes

        out.extend_from_slice(&data[offset..offset + 10]); // magic, CM, FLG, MTIME, XFL, OS
        out.extend_from_slice(&10u16.to_le_bytes()); // XLEN: 6 -> 10
        out.extend_from_slice(&data[offset + 12..offset + 16]); // BC, SLEN
        out.extend_from_slice(&u16::try_from(new_total - 1).expect("BSIZE fits").to_le_bytes());
        out.extend_from_slice(b"XY"); // the foreign subfield, behind `BC`
        out.extend_from_slice(&0u16.to_le_bytes()); // ...of zero length
        out.extend_from_slice(&data[offset + 18..offset + total]); // payload and footer

        offset += total;
    }

    std::fs::write(dest, out).expect("write reframed BGZF");
}

/// gzip is accepted at the boundary, whatever it turns out to be wrapping.
///
/// fgumi reads BGZF so it can decompress blocks in parallel, and plain gzip gives
/// that up — but rejecting it outright made `fgumi` the only tool in the pipeline
/// that could not read a file `gzip` produced. It is now decompressed at the
/// boundary and re-classified, so the rest of the pipeline still sees only BGZF.
///
/// Each case must produce exactly what the equivalent BAM input produces; matching
/// exit status alone would not catch a transcode that dropped records.
#[rstest]
// `gzip -c reads.sam`: text inside, so it takes the SAM transcode.
#[case::gzipped_sam(false, false)]
// `gzip -c reads.bam`: BGZF inside, which needs no transcode at all.
#[case::gzipped_bam(true, false)]
// A BGZF file carrying a second extra subfield. Spec-legal, unframeable by every
// decoder in the pipeline, and a plain gzip member — so it arrives here.
#[case::bgzf_with_a_trailing_subfield(true, true)]
fn gzip_input_is_accepted_and_matches_bam(#[case] from_bam: bool, #[case] reframe: bool) {
    let dir = TempDir::new().expect("create temp dir");
    let (bam_input, sam_input) = write_input_pair(dir.path());
    let args = &["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"];

    let gzipped = dir.path().join("input.gz");
    if reframe {
        add_trailing_extra_subfield(&bam_input, &gzipped);
    } else {
        gzip_file(if from_bam { &bam_input } else { &sam_input }, &gzipped);
    }

    let gzip_output = dir.path().join("from_gzip.bam");
    let result = run(args, &gzipped, &gzip_output);
    assert!(
        result.status.success(),
        "gzip input rejected: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    // Succeeding is not enough: the records have to survive the round trip.
    let bam_output = dir.path().join("from_bam.bam");
    let baseline = run(args, &bam_input, &bam_output);
    assert!(baseline.status.success(), "BAM input failed, so the gzip result proves nothing");
    assert_eq!(
        read_output(&gzip_output),
        read_output(&bam_output),
        "gzip input produced different records than the same data as BGZF",
    );
}

/// Two gzip layers are named rather than peeled.
///
/// Unwrapping until something recognisable appears would let a small crafted input
/// cost unbounded work for one open, so the boundary decompresses exactly once and
/// says what it found.
#[test]
fn doubly_gzipped_input_is_rejected() {
    let dir = TempDir::new().expect("create temp dir");
    let (_, sam_input) = write_input_pair(dir.path());

    let once = dir.path().join("once.gz");
    gzip_file(&sam_input, &once);
    let twice = dir.path().join("twice.gz");
    gzip_file(&once, &twice);

    let output = dir.path().join("out.bam");
    let result =
        run(&["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"], &twice, &output);

    assert!(!result.status.success(), "doubly-gzipped input was accepted");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("gzip-compressed twice"),
        "the error should name the double compression, got: {stderr}"
    );
}

/// A gzip member holding nothing is named as such, not as a SAM parse failure.
///
/// The empty-input check at the boundary runs on the bytes *as they arrive*, so a
/// gzip member that decompresses to nothing slips past it and is caught on the
/// second classification instead. Without that arm the empty result would take the
/// SAM transcode and fail with a message about an absent `@` header line, which
/// says nothing about the file being empty.
#[test]
fn empty_gzip_member_is_rejected() {
    let dir = TempDir::new().expect("create temp dir");

    let empty = dir.path().join("empty");
    File::create(&empty).expect("create empty source");
    let gzipped = dir.path().join("empty.gz");
    gzip_file(&empty, &gzipped);

    let output = dir.path().join("out.bam");
    let result = run(
        &["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"],
        &gzipped,
        &output,
    );

    assert!(!result.status.success(), "an empty gzip member was accepted");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("empty gzip member"),
        "the error should name the empty member, got: {stderr}"
    );
}
