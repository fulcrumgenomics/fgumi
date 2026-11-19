//! CRAM through the CLI. With feature `cram` (module `with_cram`), fgumi
//! commands write CRAM when the output ends in `.cram`, and read CRAM back
//! by content. Builds without the feature reject CRAM output names.

use std::io::Write;
use std::path::Path;
use std::process::Command;

const REF_LEN: usize = 4_000;

fn reference() -> Vec<u8> {
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    (0..REF_LEN)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[usize::try_from(x % 4).unwrap()]
        })
        .collect()
}

/// An unsorted SAM of 400 mapped pairs plus one unmapped pair, and its FASTA.
fn write_inputs(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf, usize) {
    write_inputs_with(dir, 400)
}

/// An unsorted SAM of `pairs` mapped pairs plus one unmapped pair, and its
/// FASTA.
fn write_inputs_with(dir: &Path, pairs: usize) -> (std::path::PathBuf, std::path::PathBuf, usize) {
    let seq = reference();
    let fa = dir.join("ref.fa");
    let mut f = std::fs::File::create(&fa).unwrap();
    writeln!(f, ">chrT").unwrap();
    for line in seq.chunks(60) {
        f.write_all(line).unwrap();
        writeln!(f).unwrap();
    }
    let sam = dir.join("in.sam");
    let mut s = std::fs::File::create(&sam).unwrap();
    writeln!(s, "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:{REF_LEN}").unwrap();
    let mut n = 0;
    for i in 0..pairs {
        // Scrambled positions so the sort has work to do.
        let p1 = (i * 37) % (REF_LEN - 400);
        let p2 = p1 + 150;
        let r1 = std::str::from_utf8(&seq[p1..p1 + 60]).unwrap();
        let r2 = std::str::from_utf8(&seq[p2..p2 + 60]).unwrap();
        let q = "I".repeat(60);
        writeln!(
            s,
            "q{i}\t99\tchrT\t{}\t60\t60M\t=\t{}\t210\t{r1}\t{q}\tRX:Z:ACGT",
            p1 + 1,
            p2 + 1
        )
        .unwrap();
        writeln!(
            s,
            "q{i}\t147\tchrT\t{}\t60\t60M\t=\t{}\t-210\t{r2}\t{q}\tRX:Z:ACGT",
            p2 + 1,
            p1 + 1
        )
        .unwrap();
        n += 2;
    }
    writeln!(s, "u0\t77\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII").unwrap();
    writeln!(s, "u0\t141\t*\t0\t0\t*\t*\t0\t0\tTTTTACGTAC\tIIIIIIIIII").unwrap();
    n += 2;
    (sam, fa, n)
}

/// Run fgumi with the top-level `--reference` (before the subcommand).
fn fgumi(args: &[&str], reference: &Path) {
    let out = fgumi_output(args, reference);
    assert!(
        out.status.success(),
        "fgumi {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// [`fgumi`]'s run, returned rather than asserted.
fn fgumi_output(args: &[&str], reference: &Path) -> std::process::Output {
    fgumi_command().arg("--reference").arg(reference).args(args).output().expect("run fgumi")
}

/// Ignore host reference lookup configuration so fixtures are self-contained.
fn fgumi_command() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    command.env_remove("FGUMI_REFERENCE").env_remove("REF_PATH").env_remove("REF_CACHE");
    command
}

/// Tests that need fgumi built with the `cram` feature.
#[cfg(feature = "cram")]
mod with_cram {
    use super::*;
    use crate::helpers::chain_spec::logged_count;
    use crate::helpers::cli::utf8_path;
    #[cfg(feature = "simulate")]
    use fgumi_lib::sam::SamTag;

    use noodles::bam;
    use noodles::sam;

    /// Every record of a BAM, in file order, as `RecordBuf`s, with `MD`/`NM`
    /// removed when `strip_md_nm` (CRAM decoding regenerates them, as `samtools
    /// view` does; that is checked separately) — everything else must match.
    fn bam_records(path: &Path, strip_md_nm: bool) -> Vec<sam::alignment::RecordBuf> {
        let mut r = bam::io::reader::Builder.build_from_path(path).unwrap();
        let header: sam::Header = r.read_header().unwrap();
        r.records()
            .map(|rec| {
                normalize_record(
                    sam::alignment::RecordBuf::try_from_alignment_record(&header, &rec.unwrap())
                        .unwrap(),
                    strip_md_nm,
                )
            })
            .collect()
    }

    /// Decode slices sequentially so the returned records retain emitted order.
    /// Pass the reference explicitly, without changing process-global lookup settings.
    fn cram_records(path: &Path, reference: &Path) -> Vec<sam::alignment::RecordBuf> {
        let mut reader = fgumi_htslib::CramReader::open(path, Some(reference)).unwrap();
        let header: sam::Header = reader.header_text().parse().unwrap();
        let mut records = Vec::new();
        while let Some(slice) = reader.next_slice().unwrap() {
            let decoded = slice.decode().unwrap();
            let mut framed = bam::io::Reader::from(&decoded.bytes[..]);
            for record in framed.records() {
                records.push(normalize_record(
                    sam::alignment::RecordBuf::try_from_alignment_record(&header, &record.unwrap())
                        .unwrap(),
                    true,
                ));
            }
        }
        records
    }

    /// Normalize only the documented record-level differences of a CRAM round trip.
    fn normalize_record(
        mut record: sam::alignment::RecordBuf,
        strip_md_nm: bool,
    ) -> sam::alignment::RecordBuf {
        use sam::alignment::record::data::field::Tag;
        if strip_md_nm {
            record.data_mut().remove(&Tag::MISMATCHED_POSITIONS);
            record.data_mut().remove(&Tag::EDIT_DISTANCE);
        }
        // CRAM does not preserve MAPQ/CIGAR for unmapped records.
        if record.flags().is_unmapped() {
            *record.mapping_quality_mut() = None;
            record.cigar_mut().as_mut().clear();
        }
        // Aux tag order is not content; htslib moves RG/MD/NM to the end.
        let mut fields: Vec<_> = record.data().iter().map(|(t, v)| (t, v.clone())).collect();
        fields.sort_by_key(|(t, _)| <[u8; 2]>::from(*t));
        *record.data_mut() = fields.into_iter().collect();
        record
    }

    /// Whether every mapped record carries `MD` and `NM`.
    fn mapped_have_md_nm(records: &[sam::alignment::RecordBuf]) -> bool {
        use sam::alignment::record::data::field::Tag;
        records.iter().filter(|r| !r.flags().is_unmapped()).all(|r| {
            r.data().get(&Tag::MISMATCHED_POSITIONS).is_some()
                && r.data().get(&Tag::EDIT_DISTANCE).is_some()
        })
    }

    /// CRAI entries' container offsets (field 4 of each gunzipped line).
    fn crai_offsets(path: &Path) -> Vec<u64> {
        use std::io::Read;
        let mut text = String::new();
        flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text.lines().map(|l| l.split('\t').nth(3).unwrap().parse().unwrap()).collect()
    }

    /// The local reference is consumed by CRAM I/O, even without a methylation
    /// or alignment stage. A renamed CRAM must be recognized by content.
    #[rstest::rstest]
    #[case::cram_output(false, true, false, None)]
    #[case::renamed_cram_input(true, false, false, None)]
    #[case::cram_input_and_output(true, true, false, None)]
    #[case::unused_bam_reference(false, false, false, Some("--ref requires --methylation-mode"))]
    #[case::unused_methylation(false, true, true, Some("dead on a runall chain"))]
    fn runall_sort_reference_consumers(
        #[case] cram_input: bool,
        #[case] cram_output: bool,
        #[case] methylation: bool,
        #[case] error: Option<&str>,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, _) = write_inputs(dir.path());
        let expected = dir.path().join("expected.bam");
        sort(&sam, &expected, "template-coordinate", &fa);
        let input = if cram_input {
            let cram = dir.path().join("input.cram");
            sort(&sam, &cram, "template-coordinate", &fa);
            let renamed = dir.path().join("input.data");
            std::fs::rename(cram, &renamed).unwrap();
            renamed
        } else {
            expected.clone()
        };
        let output = dir.path().join(if cram_output { "runall.cram" } else { "runall.bam" });
        let mut command = fgumi_command();
        command.args([
            "runall",
            "--start-from",
            "sort",
            "--stop-after",
            "sort",
            "-i",
            utf8_path(&input),
            "-o",
            utf8_path(&output),
            "--ref",
            utf8_path(&fa),
        ]);
        if methylation {
            command.args(["--methylation-mode", "em-seq"]);
        }
        let result = command.output().unwrap();
        if let Some(error) = error {
            assert!(!result.status.success());
            assert!(
                String::from_utf8_lossy(&result.stderr).contains(error),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            assert!(!output.exists());
        } else {
            assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
            let actual =
                if cram_output { cram_records(&output, &fa) } else { bam_records(&output, true) };
            assert_eq!(actual, bam_records(&expected, true));
        }
    }

    /// Streaming input is classified once, and an unused local reference must
    /// fail before an existing destination is touched.
    #[rstest::rstest]
    #[case::cram(true)]
    #[case::bam(false)]
    fn runall_stream_reference_consumers(#[case] cram: bool) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, _) = write_inputs(dir.path());
        let expected = dir.path().join("expected.bam");
        sort(&sam, &expected, "template-coordinate", &fa);
        let input = dir.path().join(if cram { "stream.cram" } else { "stream.bam" });
        sort(&sam, &input, "template-coordinate", &fa);
        let output = dir.path().join("out.bam");
        std::fs::write(&output, b"preserve existing output").unwrap();
        let mut child = fgumi_command()
            .args([
                "runall",
                "--start-from",
                "sort",
                "--stop-after",
                "sort",
                "-i",
                "-",
                "-o",
                utf8_path(&output),
                "--ref",
                utf8_path(&fa),
            ])
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        let mut stdin = child.stdin.take().unwrap();
        let bytes = std::fs::read(input).unwrap();
        let writer = std::thread::spawn(move || {
            let _ = stdin.write_all(&bytes);
        });
        let result = child.wait_with_output().unwrap();
        writer.join().unwrap();
        if cram {
            assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
            assert_eq!(bam_records(&output, true), bam_records(&expected, true));
        } else {
            assert!(!result.status.success());
            assert!(
                String::from_utf8_lossy(&result.stderr)
                    .contains("--ref requires --methylation-mode")
            );
            assert_eq!(std::fs::read(&output).unwrap(), b"preserve existing output");
        }
    }

    /// A `.bam` FIFO still uses content detection, without reopening and losing
    /// its prefix. A nonblocking producer also exits when validation rejects input.
    #[cfg(unix)]
    #[rstest::rstest]
    #[case::cram(true)]
    #[case::bam(false)]
    fn runall_fifo_reference_consumers(#[case] cram: bool) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, _) = write_inputs_with(dir.path(), 3);
        let expected = dir.path().join("expected.bam");
        sort(&sam, &expected, "template-coordinate", &fa);
        let input = dir.path().join(if cram { "stream.cram" } else { "stream.bam" });
        sort(&sam, &input, "template-coordinate", &fa);
        let fifo = dir.path().join("input.bam");
        assert!(Command::new("mkfifo").arg(&fifo).status().unwrap().success());
        let output = dir.path().join("out.bam");
        std::fs::write(&output, b"preserve existing output").unwrap();
        let mut child = fgumi_command()
            .args([
                "runall",
                "--start-from",
                "sort",
                "--stop-after",
                "sort",
                "-i",
                utf8_path(&fifo),
                "-o",
                utf8_path(&output),
                "--ref",
                utf8_path(&fa),
            ])
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        let bytes = std::fs::read(&input).unwrap();
        assert!(bytes.len() < 4096, "fixture fits a pipe even when input is rejected");
        // WRONLY opens only once the real consumer has opened its read end.
        // Nonblocking retries also detect a consumer rejecting input before open.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        loop {
            match rustix::fs::open(
                &fifo,
                rustix::fs::OFlags::WRONLY | rustix::fs::OFlags::NONBLOCK,
                rustix::fs::Mode::empty(),
            ) {
                Ok(fd) => {
                    std::fs::File::from(fd).write_all(&bytes).unwrap();
                    break;
                }
                Err(rustix::io::Errno::NXIO) => {
                    if child.try_wait().unwrap().is_some() {
                        break;
                    }
                    if std::time::Instant::now() >= deadline {
                        child.kill().unwrap();
                        child.wait().unwrap();
                        panic!("runall never opened its FIFO input");
                    }
                    std::thread::yield_now();
                }
                Err(error) => panic!("opening FIFO producer: {error}"),
            }
        }
        let result = child.wait_with_output().unwrap();
        if cram {
            assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
            assert_eq!(bam_records(&output, true), bam_records(&expected, true));
        } else {
            assert!(!result.status.success());
            assert!(
                String::from_utf8_lossy(&result.stderr)
                    .contains("--ref requires --methylation-mode")
            );
            assert_eq!(std::fs::read(&output).unwrap(), b"preserve existing output");
        }
    }

    #[test]
    fn sort_writes_cram_with_crai_and_reads_it_back() {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, n) = write_inputs(dir.path());
        let p = |name: &str| dir.path().join(name).to_str().unwrap().to_string();

        // SAM -> coordinate-sorted CRAM (+ CRAI), and the same sort to BAM.
        let sam_path = sam.to_str().unwrap();
        fgumi(
            &[
                "sort",
                "-i",
                sam_path,
                "-o",
                &p("out.cram"),
                "--order",
                "coordinate",
                "--write-index",
            ],
            &fa,
        );
        fgumi(&["sort", "-i", sam_path, "-o", &p("out.bam"), "--order", "coordinate"], &fa);
        let cram = std::fs::read(p("out.cram")).unwrap();
        assert_eq!(&cram[..4], b"CRAM", "a .cram output is CRAM");
        assert_eq!(
            &cram[6..13],
            b"out.cra",
            "the file ID is the destination name, not the staged temp"
        );
        assert!(!dir.path().join("out.cram.bai").exists());
        let offsets = crai_offsets(&dir.path().join("out.cram.crai"));
        assert!(!offsets.is_empty(), "the CRAI indexes the file");
        assert!(
            offsets.iter().all(|&o| o > 26 && o < cram.len() as u64),
            "CRAI offsets inside the file"
        );

        assert_eq!(
            assert_cram_records_match(Path::new(&p("out.cram")), Path::new(&p("out.bam")), &fa)
                .len(),
            n
        );

        // CRAM input is detected by content: re-sort it to BAM and compare whole
        // records (MD/NM, which CRAM decoding adds, are pinned separately).
        fgumi(&["sort", "-i", &p("out.cram"), "-o", &p("back.bam"), "--order", "coordinate"], &fa);
        let via_cram = bam_records(Path::new(&p("back.bam")), false);
        assert!(mapped_have_md_nm(&via_cram), "CRAM decoding regenerates MD/NM on mapped records");
        let via_cram = bam_records(Path::new(&p("back.bam")), true);
        let direct = bam_records(Path::new(&p("out.bam")), true);
        assert_eq!(via_cram.len(), n);
        assert_eq!(via_cram, direct, "CRAM round trip changed records or their order");

        // The decode-free fast path (a per-record first stage) reads CRAM too.
        fgumi(&["retag", "-i", &p("out.cram"), "-o", &p("retag_cram.bam"), "RX::delete"], &fa);
        fgumi(&["retag", "-i", &p("out.bam"), "-o", &p("retag_bam.bam"), "RX::delete"], &fa);
        assert_eq!(
            bam_records(Path::new(&p("retag_cram.bam")), true),
            bam_records(Path::new(&p("retag_bam.bam")), true),
            "retag over CRAM differs from retag over BAM"
        );
    }

    /// An in-file CRAI offset is not enough: indexed reads must return exactly
    /// the fixture alignments overlapping the requested region, including none.
    #[rstest::rstest]
    #[case::populated_region("chrT:1-500", 1, 500, 113)]
    #[case::empty_region("chrT:3901-4000", 3901, 4000, 0)]
    #[ignore = "requires samtools; gated by the cram CI job"]
    fn cram_index_queries_return_expected_records(
        #[case] region: &str,
        #[case] start: usize,
        #[case] end: usize,
        #[case] count: usize,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, _) = write_inputs(dir.path());
        let cram = dir.path().join("indexed.cram");
        fgumi(
            &[
                "sort",
                "-i",
                utf8_path(&sam),
                "-o",
                utf8_path(&cram),
                "--order",
                "coordinate",
                "--write-index",
            ],
            &fa,
        );
        let offsets = crai_offsets(&dir.path().join("indexed.cram.crai"));
        let length = std::fs::metadata(&cram).unwrap().len();
        assert!(!offsets.is_empty());
        assert!(offsets.iter().all(|&offset| offset > 26 && offset < length));

        let queried = dir.path().join("queried.bam");
        let result = Command::new("samtools")
            .args([
                "view",
                "-b",
                "-T",
                utf8_path(&fa),
                "-o",
                utf8_path(&queried),
                utf8_path(&cram),
                region,
            ])
            .env_remove("REF_PATH")
            .env_remove("REF_CACHE")
            .output()
            .expect("samtools must be installed for indexed CRAM tests");
        assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
        let mut input =
            sam::io::Reader::new(std::io::BufReader::new(std::fs::File::open(&sam).unwrap()));
        let header = input.read_header().unwrap();
        let mut expected: Vec<_> = input
            .record_bufs(&header)
            .map(Result::unwrap)
            .filter(|record| {
                !record.flags().is_unmapped()
                    && record.alignment_start().is_some_and(|position| usize::from(position) <= end)
                    && record.alignment_end().is_some_and(|position| usize::from(position) >= start)
            })
            .map(|record| normalize_record(record, true))
            .collect();
        let actual = bam_records(&queried, true);
        assert_eq!(expected.len(), count, "independent generated fixture oracle");
        assert_eq!(actual.len(), count);
        for record in actual {
            let index = expected
                .iter()
                .position(|expected| expected == &record)
                .expect("indexed query returned an extra or changed record");
            expected.swap_remove(index);
        }
        assert!(expected.is_empty(), "indexed query missed expected alignments");
    }

    /// `extract -o *.cram` writes unmapped CRAM through the pipeline's CRAM sink
    /// (no reference: reference-free) with the same records as `-o *.bam`.
    #[test]
    fn extract_writes_unmapped_cram() {
        let dir = tempfile::tempdir().unwrap();
        let fq = dir.path().join("r.fq");
        let mut f = std::fs::File::create(&fq).unwrap();
        for i in 0..300 {
            writeln!(f, "@q{i}\nACGTACGTTTGGCCAAACGTACGTACGTACGTACGTACGT\n+\n{}", "I".repeat(40))
                .unwrap();
        }
        drop(f);
        let run = |out: &Path| {
            let o = fgumi_command()
                .args(["extract", "-i"])
                .arg(&fq)
                .args(["-r", "8M+T", "--sample", "s", "--library", "l", "-o"])
                .arg(out)
                .env_remove("FGUMI_REFERENCE")
                .env_remove("REF_PATH")
                .output()
                .expect("run fgumi");
            assert!(o.status.success(), "extract failed: {}", String::from_utf8_lossy(&o.stderr));
        };
        let (cram, bam) = (dir.path().join("u.cram"), dir.path().join("u.bam"));
        run(&cram);
        run(&bam);
        assert_eq!(&std::fs::read(&cram).unwrap()[..4], b"CRAM");
        let (_, fa, _) = write_inputs(dir.path());
        assert_eq!(assert_cram_records_match(&cram, &bam, &fa).len(), 300);
    }

    /// `fgumi sort -i <input> -o <output> --order <order>`, asserted.
    fn sort(input: &Path, output: &Path, order: &str, reference: &Path) {
        fgumi(
            &["sort", "-i", utf8_path(input), "-o", utf8_path(output), "--order", order],
            reference,
        );
    }

    /// A sort that spills runs and merges them writes the same records to
    /// CRAM as to BAM.
    #[test]
    fn a_spilled_sort_writes_cram() {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, n) = write_inputs_with(dir.path(), 5_000);
        let (cram, bam) = (dir.path().join("out.cram"), dir.path().join("out.bam"));
        for out in [&cram, &bam] {
            let o = fgumi_output(
                &[
                    "sort",
                    "-i",
                    utf8_path(&sam),
                    "-o",
                    utf8_path(out),
                    "--order",
                    "coordinate",
                    "--threads",
                    "2",
                    "--max-memory",
                    "64K",
                    "--sort-stats",
                    "--tmp-dir",
                    utf8_path(dir.path()),
                ],
                &fa,
            );
            let stderr = String::from_utf8_lossy(&o.stderr);
            assert!(o.status.success(), "{stderr}");
            let runs = logged_count(&stderr, "Spill runs:")
                .unwrap_or_else(|| panic!("no `Spill runs: N` line:\n{stderr}"));
            assert!(runs > 1, "the sort must spill and merge:\n{stderr}");
        }
        assert_eq!(&std::fs::read(&cram).unwrap()[..4], b"CRAM");
        assert_eq!(assert_cram_records_match(&cram, &bam, &fa).len(), n);
    }

    /// A CRAM on stdin (`sort -i -`) is detected by content and decoded.
    #[test]
    fn sort_reads_cram_from_stdin() {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, n) = write_inputs(dir.path());
        let (cram, bam) = (dir.path().join("in.cram"), dir.path().join("direct.bam"));
        sort(&sam, &cram, "coordinate", &fa);
        sort(&sam, &bam, "coordinate", &fa);
        let back = dir.path().join("back.bam");
        let o = fgumi_command()
            .arg("--reference")
            .arg(&fa)
            .args(["sort", "-i", "-", "-o", utf8_path(&back), "--order", "coordinate"])
            .stdin(std::fs::File::open(&cram).unwrap())
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        let via_stdin = bam_records(&back, true);
        assert_eq!(via_stdin.len(), n);
        assert_eq!(via_stdin, bam_records(&bam, true));
    }

    /// A truncated CRAM, read from a file or from stdin, fails the run and
    /// publishes nothing at the destination: one cut before its EOF container
    /// (the last 38 bytes of a CRAM 3), one inside its only data container.
    #[rstest::rstest]
    fn a_truncated_cram_fails_and_writes_nothing(
        #[values("missing_eof", "mid_container")] cut: &str,
        #[values(false, true)] stdin: bool,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, _) = write_inputs(dir.path());
        let cram = dir.path().join("whole.cram");
        sort(&sam, &cram, "coordinate", &fa);
        let bytes = std::fs::read(&cram).unwrap();
        let keep = if cut == "missing_eof" { bytes.len() - 38 } else { bytes.len() / 2 };
        let truncated = dir.path().join("truncated.cram");
        std::fs::write(&truncated, &bytes[..keep]).unwrap();

        let dest = dir.path().join("dest.bam");
        let mut cmd = fgumi_command();
        cmd.arg("--reference").arg(&fa).args(["sort", "-i"]);
        if stdin {
            cmd.arg("-").stdin(std::fs::File::open(&truncated).unwrap());
        } else {
            cmd.arg(&truncated);
        }
        let o = cmd.args(["-o", utf8_path(&dest), "--order", "coordinate"]).output().unwrap();
        assert!(
            !o.status.success(),
            "a truncated CRAM must fail the run: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        assert!(!dest.exists(), "nothing is published at the destination");
    }

    /// The unmapped reads of [`write_inputs`]'s SAM: each record with its
    /// alignment removed and its read restored to sequencing orientation.
    fn write_unmapped(sam: &Path, out: &Path) {
        let text = std::fs::read_to_string(sam).unwrap();
        let mut f = std::fs::File::create(out).unwrap();
        writeln!(f, "@HD\tVN:1.6\tSO:unsorted").unwrap();
        for line in text.lines().filter(|l| !l.starts_with('@')) {
            let fields: Vec<&str> = line.split('\t').collect();
            let flag: u16 = fields[1].parse().unwrap();
            let (mut seq, mut qual) = (fields[9].to_string(), fields[10].to_string());
            if flag & 0x10 != 0 {
                seq = fgumi_dna::reverse_complement_str(&seq);
                qual = qual.chars().rev().collect();
            }
            let unmapped_flag = if flag & 0x40 != 0 { 77 } else { 141 };
            let tags = fields[11..].join("\t");
            let sep = if tags.is_empty() { "" } else { "\t" };
            writeln!(
                f,
                "{}\t{unmapped_flag}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}{sep}{tags}",
                fields[0]
            )
            .unwrap();
        }
    }

    /// Zipper reads and writes CRAM on the serial and chain paths, including
    /// mapped input streamed on stdin.
    #[rstest::rstest]
    fn zipper_reads_cram_inputs(
        #[values("1", "2", "4")] threads: &str,
        #[values(false, true)] unmapped_cram: bool,
        #[values(false, true)] stdin: bool,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, n) = write_inputs(dir.path());
        std::fs::write(
            dir.path().join("ref.dict"),
            format!("@HD\tVN:1.6\n@SQ\tSN:chrT\tLN:{REF_LEN}\n"),
        )
        .unwrap();
        let p = |name: &str| dir.path().join(name);
        sort(&sam, &p("mapped.cram"), "queryname", &fa);
        sort(&sam, &p("mapped.bam"), "queryname", &fa);
        write_unmapped(&sam, &p("u.sam"));
        let unmapped = if unmapped_cram { p("u.cram") } else { p("u.bam") };
        sort(&p("u.sam"), &unmapped, "queryname", &fa);

        let zip = |mapped: &Path, out: &Path| {
            let mut command = fgumi_command();
            command.arg("--reference").arg(&fa);
            let mapped_input = if stdin { "-" } else { utf8_path(mapped) };
            if stdin {
                command.stdin(std::fs::File::open(mapped).unwrap());
            }
            let result = command
                .args([
                    "zipper",
                    "-i",
                    mapped_input,
                    "-u",
                    utf8_path(&unmapped),
                    "-r",
                    utf8_path(&fa),
                    "-o",
                    utf8_path(out),
                    "--threads",
                    threads,
                ])
                .output()
                .unwrap();
            assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
        };
        zip(&p("mapped.cram"), &p("z_from_cram.bam"));
        zip(&p("mapped.bam"), &p("z_from_bam.bam"));
        let from_cram = bam_records(&p("z_from_cram.bam"), true);
        assert_eq!(from_cram.len(), n);
        assert_eq!(from_cram, bam_records(&p("z_from_bam.bam"), true));

        zip(&p("mapped.bam"), &p("z.cram"));
        assert_eq!(&std::fs::read(p("z.cram")).unwrap()[..4], b"CRAM");
        assert_eq!(assert_cram_records_match(&p("z.cram"), &p("z_from_bam.bam"), &fa).len(), n);
    }

    /// Compare normalized records in emitted order, without a sort that could
    /// conceal an ordering regression.
    fn assert_cram_records_match(
        cram: &Path,
        bam: &Path,
        reference: &Path,
    ) -> Vec<sam::alignment::RecordBuf> {
        assert_eq!(&std::fs::read(cram).unwrap()[..4], b"CRAM");
        let records = cram_records(cram, reference);
        assert_eq!(records, bam_records(bam, true), "CRAM changed records or emitted order");
        records
    }

    /// Merge stages its destination under a temporary name; format selection and
    /// the CRAM file ID must still use the requested output name.
    #[rstest::rstest]
    #[case::coordinate("coordinate", [2, 3, 0, 1])]
    #[case::queryname("queryname", [0, 1, 2, 3])]
    fn merge_writes_cram(
        #[values("1", "4")] threads: &str,
        #[case] order: &str,
        #[case] expected_indices: [usize; 4],
    ) {
        use std::fmt::Write as _;

        let dir = tempfile::tempdir().unwrap();
        let (_, fa, _) = write_inputs(dir.path());
        let bases = reference();
        let mut expected = Vec::new();
        let inputs = [dir.path().join("source_a.bam"), dir.path().join("source_b.bam")];
        // Distinct source names, positions, sequences and tags catch a source
        // being lost and the other duplicated, even at the right total count.
        for (source, indices) in [[0, 2], [1, 3]].into_iter().enumerate() {
            let mut text = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:{REF_LEN}\n");
            for index in indices {
                let pos = [301, 401, 101, 201][index];
                let seq = std::str::from_utf8(&bases[pos - 1..pos + 19]).unwrap();
                writeln!(
                    text,
                    "a{index}\t0\tchrT\t{pos}\t60\t20M\t*\t0\t0\t{seq}\t{}\tZS:Z:source_{source}",
                    "I".repeat(20)
                )
                .unwrap();
            }
            let source_sam = dir.path().join(format!("source_{source}.sam"));
            std::fs::write(&source_sam, &text).unwrap();
            let mut reader = sam::io::Reader::new(std::io::Cursor::new(text));
            let header = reader.read_header().unwrap();
            expected
                .extend(reader.record_bufs(&header).map(|r| normalize_record(r.unwrap(), true)));
            sort(&source_sam, &inputs[source], order, &fa);
        }
        // Expected record order is hand-derived from the fixture, not a second
        // invocation of the merge/sort implementation under test.
        let by_name: std::collections::HashMap<_, _> =
            expected.into_iter().map(|r| (r.name().unwrap().to_string(), r)).collect();
        let expected: Vec<_> =
            expected_indices.iter().map(|i| by_name[&format!("a{i}")].clone()).collect();
        let bam = dir.path().join("merged.bam");
        let cram = dir.path().join("merged.cram");
        for output in [&bam, &cram] {
            fgumi(
                &[
                    "merge",
                    "--order",
                    order,
                    "--threads",
                    threads,
                    "-o",
                    utf8_path(output),
                    utf8_path(&inputs[0]),
                    utf8_path(&inputs[1]),
                ],
                &fa,
            );
        }
        let bytes = std::fs::read(&cram).unwrap();
        assert_eq!(&bytes[6..18], b"merged.cram\0", "file ID must not be the staged name");
        assert_eq!(assert_cram_records_match(&cram, &bam, &fa), expected);
    }

    /// Generate template-coordinate families without relying on fixture files.
    #[cfg(feature = "simulate")]
    fn simulate_grouped(output: &Path, truth: &Path, reference: &Path, duplex: bool) {
        fgumi(
            &[
                "simulate",
                "grouped-reads",
                "-o",
                utf8_path(output),
                "--truth",
                utf8_path(truth),
                "--reference",
                utf8_path(reference),
                "--num-molecules",
                "40",
                "--seed",
                "73",
                "--threads",
                "1",
                "--duplex",
                if duplex { "true" } else { "false" },
            ],
            reference,
        );
    }

    #[cfg(feature = "simulate")]
    #[rstest::rstest]
    fn simulate_grouped_reads_writes_cram(#[values(false, true)] duplex: bool) {
        let dir = tempfile::tempdir().unwrap();
        let (_, fa, _) = write_inputs(dir.path());
        let bam = dir.path().join("grouped.bam");
        let cram = dir.path().join("grouped.cram");
        let bam_truth = dir.path().join("bam.truth.tsv");
        let cram_truth = dir.path().join("cram.truth.tsv");
        simulate_grouped(&bam, &bam_truth, &fa, duplex);
        simulate_grouped(&cram, &cram_truth, &fa, duplex);
        let records = assert_cram_records_match(&cram, &bam, &fa);
        assert!(!records.is_empty());
        assert!(records.iter().all(|r| r.data().get(&SamTag::MI.to_noodles_tag()).is_some()));
        assert_eq!(std::fs::read(&cram_truth).unwrap(), std::fs::read(&bam_truth).unwrap());
    }

    #[cfg(feature = "simulate")]
    #[rstest::rstest]
    #[case::consensus("consensus-reads", SamTag::CD)]
    #[case::correct("correct-reads", SamTag::RX)]
    // These generators emit channel-arrival order, not an ordered parallel
    // sequence. One generation thread makes cross-run emitted-order parity valid.
    fn simulate_reads_writes_cram(#[case] command: &str, #[case] tag: SamTag) {
        let dir = tempfile::tempdir().unwrap();
        let (_, fa, _) = write_inputs(dir.path());
        let bam = dir.path().join("reads.bam");
        let cram = dir.path().join("reads.cram");
        for (format, output) in [("bam", &bam), ("cram", &cram)] {
            let truth = dir.path().join(format!("{format}.truth.tsv"));
            let includelist = dir.path().join(format!("{format}.umis.txt"));
            let mut args = vec![
                "simulate",
                command,
                "-o",
                utf8_path(output),
                "--truth",
                utf8_path(&truth),
                "--num-reads",
                "20",
                "--seed",
                "73",
                "--threads",
                "1",
            ];
            if command == "correct-reads" {
                // correct-reads generates the includelist as well as the alignments.
                args.extend(["--includelist", utf8_path(&includelist), "--num-umis", "10"]);
            } else {
                args.extend(["--reference", utf8_path(&fa), "--duplex"]);
            }
            fgumi(&args, &fa);
        }
        let records = assert_cram_records_match(&cram, &bam, &fa);
        assert_eq!(records.len(), 40);
        assert!(records.iter().all(|r| r.data().get(&tag.to_noodles_tag()).is_some()));
        assert_eq!(
            std::fs::read(dir.path().join("cram.truth.tsv")).unwrap(),
            std::fs::read(dir.path().join("bam.truth.tsv")).unwrap()
        );
        if command == "correct-reads" {
            assert_eq!(
                std::fs::read(dir.path().join("cram.umis.txt")).unwrap(),
                std::fs::read(dir.path().join("bam.umis.txt")).unwrap()
            );
        }
    }

    #[cfg(feature = "simulate")]
    #[test]
    fn downsample_writes_cram_main_and_rejects() {
        let dir = tempfile::tempdir().unwrap();
        let (_, fa, _) = write_inputs(dir.path());
        let input = dir.path().join("grouped.bam");
        simulate_grouped(&input, &dir.path().join("truth.tsv"), &fa, true);
        for extension in ["bam", "cram"] {
            let kept = dir.path().join(format!("kept.{extension}"));
            let rejected = dir.path().join(format!("rejected.{extension}"));
            fgumi(
                &[
                    "downsample",
                    "-i",
                    utf8_path(&input),
                    "-o",
                    utf8_path(&kept),
                    "--rejects",
                    utf8_path(&rejected),
                    "--fraction",
                    "0.5",
                    "--seed",
                    "73",
                ],
                &fa,
            );
        }
        let mut remaining = bam_records(&input, true);
        for name in ["kept", "rejected"] {
            let records = assert_cram_records_match(
                &dir.path().join(format!("{name}.cram")),
                &dir.path().join(format!("{name}.bam")),
                &fa,
            );
            assert!(!records.is_empty(), "exercise both output writers");
            assert!(records.iter().all(|r| r.data().get(&SamTag::MI.to_noodles_tag()).is_some()));
            // Remove one occurrence at a time so duplicates cannot hide a missing record.
            for record in records {
                let index = remaining
                    .iter()
                    .position(|expected| expected == &record)
                    .expect("output partition contains an extra or changed record");
                remaining.swap_remove(index);
            }
        }
        assert!(remaining.is_empty(), "input records missing from the combined partitions");
    }

    #[cfg(feature = "compare")]
    #[rstest::rstest]
    fn compare_bam_roundtrip_writes_cram(#[values("1", "4")] threads: &str) {
        let dir = tempfile::tempdir().unwrap();
        let (sam, fa, n) = write_inputs(dir.path());
        let input = dir.path().join("in.bam");
        sort(&sam, &input, "queryname", &fa);
        let bam = dir.path().join("roundtrip.bam");
        let cram = dir.path().join("roundtrip.cram");
        for output in [&bam, &cram] {
            fgumi(
                &[
                    "compare",
                    "bam-roundtrip",
                    utf8_path(&input),
                    "--output",
                    utf8_path(output),
                    "--threads",
                    threads,
                ],
                &fa,
            );
        }
        let records = assert_cram_records_match(&cram, &bam, &fa);
        assert_eq!(records.len(), n);
        assert_eq!(
            records,
            bam_records(&input, true),
            "roundtrip changed records or emitted order"
        );
    }

    /// An align-like stage defers the output header until the aligner has
    /// emitted its own, so a `realign -o *.cram` writes through the CRAM
    /// sink's deferred-header path. A stand-in aligner (reports every read
    /// unmapped) keeps the run free of a real aligner.
    #[test]
    fn realign_writes_cram_with_a_deferred_header() {
        use crate::helpers::bam_generator::create_test_reference;
        use crate::test_realign_command::{UNMAPPED_ALIGNER, write_consensus};

        let dir = tempfile::tempdir().unwrap();
        let input =
            write_consensus(dir.path(), 3, |_| (Some("chr1:1-100:+"), Some("chr1:1-100:-")));
        let fa = create_test_reference(dir.path());
        let p = |name: &str| dir.path().join(name);
        for out in [p("out.cram"), p("out.bam")] {
            let o = fgumi_command()
                .args([
                    "realign",
                    "-i",
                    utf8_path(&input),
                    "-o",
                    utf8_path(&out),
                    "-r",
                    utf8_path(&fa),
                ])
                .args(["--aligner::command", UNMAPPED_ALIGNER])
                .env_remove("FGUMI_REFERENCE")
                .output()
                .unwrap();
            assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        }
        assert_eq!(&std::fs::read(p("out.cram")).unwrap()[..4], b"CRAM");
        assert_eq!(assert_cram_records_match(&p("out.cram"), &p("out.bam"), &fa).len(), 6);
    }
}

/// A build without CRAM support fails helpfully rather than publishing BAM
/// under a CRAM output name.
#[cfg(not(feature = "cram"))]
#[test]
fn merge_refuses_a_cram_output_without_the_feature() {
    let dir = tempfile::tempdir().unwrap();
    let (sam, fa, _) = write_inputs(dir.path());
    let sorted = dir.path().join("in.bam");
    fgumi(
        &[
            "sort",
            "-i",
            sam.to_str().unwrap(),
            "-o",
            sorted.to_str().unwrap(),
            "--order",
            "coordinate",
        ],
        &fa,
    );
    let out = fgumi_command()
        .args(["merge", "--order", "coordinate", "-o"])
        .arg(dir.path().join("m.cram"))
        .arg(&sorted)
        .output()
        .unwrap();
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--features cram"), "unexpected error: {stderr}");
    assert!(!dir.path().join("m.cram").exists(), "no misnamed BAM is left behind");
}

/// Refusing an unsupported alignment format must precede all output side effects.
#[cfg(all(feature = "simulate", not(feature = "cram")))]
#[rstest::rstest]
#[case::lowercase("out.cram")]
#[case::uppercase("out.CRAM")]
fn correct_reads_cram_rejection_preserves_existing_includelist(#[case] output_name: &str) {
    let dir = tempfile::tempdir().unwrap();
    let includelist = dir.path().join("umis.txt");
    let original = b"existing includelist\nACGT\n";
    std::fs::write(&includelist, original).unwrap();
    let output = dir.path().join(output_name);
    let truth = dir.path().join("truth.tsv");
    let result = fgumi_command()
        .args(["simulate", "correct-reads", "-o"])
        .arg(&output)
        .arg("--includelist")
        .arg(&includelist)
        .arg("--truth")
        .arg(&truth)
        .args(["--num-reads", "2", "--num-umis", "2", "--threads", "1"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(!result.status.success(), "{stderr}");
    assert!(stderr.contains("--features cram"), "unexpected error: {stderr}");
    assert_eq!(std::fs::read(&includelist).unwrap(), original, "includelist was overwritten");
    assert!(!output.exists(), "no alignment output is published");
    assert!(!truth.exists(), "no truth output is published");
}

#[cfg(not(feature = "cram"))]
#[rstest::rstest]
#[case::lowercase("rejects.cram")]
#[case::uppercase("rejects.CRAM")]
fn downsample_cram_rejection_preserves_existing_main_output(#[case] rejects_name: &str) {
    let dir = tempfile::tempdir().unwrap();
    let (input, _, _) = write_inputs_with(dir.path(), 1);
    let text = std::fs::read_to_string(&input).unwrap();
    let text = text
        .lines()
        .map(|line| {
            if line.starts_with("@HD") {
                "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate".to_string()
            } else if line.starts_with('@') {
                line.to_string()
            } else {
                format!("{line}\tMI:Z:1")
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(&input, format!("{text}\n")).unwrap();
    let output = dir.path().join("kept.bam");
    let original = b"existing main output\n";
    std::fs::write(&output, original).unwrap();
    let rejects = dir.path().join(rejects_name);
    let result = fgumi_command()
        .args(["downsample", "-i"])
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .arg("--rejects")
        .arg(&rejects)
        .args(["--fraction", "0.5"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(!result.status.success(), "{stderr}");
    assert!(stderr.contains("--features cram"), "unexpected error: {stderr}");
    assert_eq!(std::fs::read(&output).unwrap(), original, "main output was overwritten");
    assert!(!rejects.exists(), "no rejected alignment output is published");
}

/// One reference per run: `fgumi --reference` and a subcommand's own reference
/// (here `filter --ref`) must agree; a different pair is rejected (clap's
/// usage error, exit 2) naming the subcommand's flag, rather than one silently
/// winning. Checked in every build, with or without `cram`.
#[test]
fn conflicting_references_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("missing.bam");
    let input = input.to_str().unwrap();
    let run = |args: &[&str]| fgumi_command().args(args).output().expect("run fgumi");
    let out =
        run(&["--reference", "a.fa", "filter", "-i", input, "-o", "out.bam", "--ref", "b.fa"]);
    assert_eq!(out.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr
            .contains("two different references: --reference a.fa and the subcommand's --ref b.fa"),
        "unexpected error: {stderr}"
    );

    // The same reference twice is fine: the run gets past the arguments and
    // the references, and fails on the missing input.
    let same =
        run(&["--reference", "a.fa", "filter", "-i", input, "-o", "out.bam", "--ref", "a.fa"]);
    let stderr = String::from_utf8_lossy(&same.stderr);
    assert_eq!(same.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("File does not exist") && stderr.contains(input), "{stderr}");
    assert!(!stderr.contains("two different references"), "{stderr}");

    // `--reference` belongs before the subcommand: after it, it is not a
    // subcommand option.
    let after = run(&["filter", "-i", input, "-o", "out.bam", "--reference", "a.fa"]);
    assert_eq!(after.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&after.stderr).contains("unexpected argument '--reference'"));
}

/// Without the `cram` feature a `.cram` output is refused, naming the user's
/// destination (not the staged temp the sort writes), and nothing is created.
#[cfg(not(feature = "cram"))]
#[test]
fn sort_refuses_a_cram_output_without_the_feature() {
    let dir = tempfile::tempdir().unwrap();
    let (sam, _, _) = write_inputs(dir.path());
    let out = dir.path().join("x.cram");
    let o = fgumi_command()
        .args(["sort", "-i"])
        .arg(&sam)
        .arg("-o")
        .arg(&out)
        .args(["--order", "coordinate"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(
        stderr.contains(&format!("CRAM output ({})", out.display()))
            && stderr.contains("--features cram"),
        "the error must name x.cram: {stderr}"
    );
    let left: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n != "in.sam" && n != "ref.fa")
        .collect();
    assert!(left.is_empty(), "nothing is created: {left:?}");
}

/// Simulation must not publish CRAM or truth output without CRAM support.
#[cfg(all(feature = "simulate", not(feature = "cram")))]
#[test]
fn simulate_grouped_reads_refuses_a_cram_output_without_the_feature() {
    let dir = tempfile::tempdir().unwrap();
    let (_, fa, _) = write_inputs(dir.path());
    let out = dir.path().join("x.cram");
    let o = fgumi_command()
        .args(["simulate", "grouped-reads", "-o"])
        .arg(&out)
        .arg("--truth")
        .arg(dir.path().join("truth.tsv"))
        .arg("--reference")
        .arg(&fa)
        .args(["--num-molecules", "10"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(!o.status.success());
    assert!(stderr.contains("--features cram"), "unexpected error: {stderr}");
    assert!(!out.exists(), "no misnamed BAM is left behind");
    assert!(!dir.path().join("truth.tsv").exists(), "no truth output is published");
}
