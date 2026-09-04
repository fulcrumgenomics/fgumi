//! Helper utilities for integration tests.
//!
//! Individual helpers are used by different subsets of the test modules
//! (some of which are feature-gated), so some helpers appear unused under
//! certain feature combinations. Silence the resulting warnings wholesale
//! so `--no-default-features` stays clean.

#![allow(dead_code, unused_imports)]

pub mod aligner;
pub mod assertions;
pub mod bam_generator;
pub mod fastq;

pub use aligner::*;
pub use assertions::*;
pub use bam_generator::*;
pub use fastq::*;

/// Writes `contents` to `dir/name` and returns the path.
///
/// Shared by the `compare` integration tests (`test_compare_metrics_command.rs`,
/// `test_compare_mutation.rs`) for building ad hoc TSV fixtures.
pub fn write_tsv(dir: &std::path::Path, name: &str, contents: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, contents).expect("failed to write temp TSV");
    path
}

/// Extracts the record count from `sort --verify`'s "Records checked: N" line.
///
/// `--verify` writes no output file, so its log is the only place the number of
/// records it actually read is observable. That couples this to the log wording
/// in `Sort::execute_verify`; the panic below says so explicitly, so a reworded
/// log line reads as "update this marker" rather than a mysterious parse failure.
pub const RECORDS_CHECKED_MARKER: &str = "Records checked:";

pub fn records_checked(stderr: &str) -> usize {
    let Some(count) = stderr.lines().find_map(|line| line.split(RECORDS_CHECKED_MARKER).nth(1))
    else {
        panic!(
            "no {RECORDS_CHECKED_MARKER:?} line in sort --verify output — if the log wording \
             changed, update RECORDS_CHECKED_MARKER. Full output:\n{stderr}"
        );
    };

    count.trim().parse().unwrap_or_else(|e| {
        panic!("{RECORDS_CHECKED_MARKER:?} was not followed by a number ({e}): {count:?}")
    })
}

/// Extracts the violation count from `sort --verify`'s "Violations: N" line.
///
/// Paired with [`records_checked`]: a format-specific record-boundary bug could
/// preserve the record total while corrupting order-violation detection, so a
/// parity test that compares only the count would still pass. Same coupling to
/// the log wording, and the same explicit panic when it changes.
pub const VIOLATIONS_MARKER: &str = "Sort order violations:";

pub fn violations(stderr: &str) -> usize {
    let Some(count) = stderr.lines().find_map(|line| line.split(VIOLATIONS_MARKER).nth(1)) else {
        panic!(
            "no {VIOLATIONS_MARKER:?} line in sort --verify output — if the log wording \
             changed, update VIOLATIONS_MARKER. Full output:\n{stderr}"
        );
    };

    count.trim().parse().unwrap_or_else(|e| {
        panic!("{VIOLATIONS_MARKER:?} was not followed by a number ({e}): {count:?}")
    })
}

/// Reads a BAM's entire header and records, for comparing two runs' outputs.
///
/// Comparing the header as well as the records is what makes a dropped `@SQ`,
/// `@RG` or sort-order line a test failure rather than something only the
/// records would have to expose. The one part that legitimately differs between
/// runs is normalized away: fgumi stamps a `@PG` line whose `CL` field is the
/// command line, which names the input path.
pub fn read_bam_output(
    path: &std::path::Path,
) -> (noodles::sam::Header, Vec<noodles::sam::alignment::RecordBuf>) {
    use noodles::sam::header::record::value::map::program::tag;

    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(path).expect("open output BAM"),
    ));
    let mut header = reader.read_header().expect("read output header");
    let records = reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .expect("read output records");

    for (_, program) in header.programs_mut().as_mut() {
        if let Some(command_line) = program.other_fields_mut().get_mut(&tag::COMMAND_LINE) {
            *command_line = bstr::BString::from("<normalized>");
        }
    }

    (header, records)
}
