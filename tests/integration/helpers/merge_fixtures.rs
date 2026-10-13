//! Fixtures for the `fgumi merge` tests: sorted inputs partitioned from the
//! sort oracle's tie-heavy fixture, the merge oracle, and a CLI runner.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use fgumi_bam_io::create_raw_bam_reader;
use fgumi_pipeline_io::sort::testing::{
    reference_sort, reference_sort_full_template, tie_heavy_records, xi_of,
};
use fgumi_raw_bam::RawRecord;
use fgumi_sort::{RawReadAheadReader, SortOrder, create_output_header};
use noodles::sam::Header;

/// Write `records` verbatim (raw bytes, no re-encoding) under `header`.
pub fn write_raw_bam<R: AsRef<[u8]>>(path: &Path, header: &Header, records: &[R]) {
    let mut w = fgumi_bam_io::create_raw_bam_writer(path, header, 1, 1).expect("create BAM");
    for r in records {
        w.write_raw_record(r.as_ref()).expect("write record");
    }
    w.finish().expect("finish BAM");
}

/// Every record's bytes and the header, as written.
pub fn read_raw(path: &Path) -> (Header, Vec<Vec<u8>>) {
    let (reader, header) = create_raw_bam_reader(path, 1).expect("open");
    (header, RawReadAheadReader::new(reader).map(|r| r.as_ref().to_vec()).collect())
}

/// `n_inputs` BAMs partitioning `tie_heavy_records(n, seed)` round-robin, each
/// stably sorted in `order` (the reference sort) and written under
/// `create_output_header(order, &header)`. Returns the paths, the common header
/// and the records concatenated in input order (input 0's sorted records, then
/// input 1's, …): the sequence whose stable sort is the merge's oracle.
pub fn sorted_inputs(
    dir: &Path,
    order: SortOrder,
    n: usize,
    n_inputs: usize,
    seed: u64,
) -> (Vec<PathBuf>, Header, Vec<RawRecord>) {
    let (header, records) = tie_heavy_records(n, seed);
    let header = create_output_header(order, &header);
    let by_xi: HashMap<u32, &RawRecord> = records.iter().map(|r| (xi_of(r.as_ref()), r)).collect();
    let mut paths = Vec::new();
    let mut concat = Vec::new();
    for i in 0..n_inputs {
        let part: Vec<RawRecord> = records.iter().skip(i).step_by(n_inputs).cloned().collect();
        let sorted: Vec<RawRecord> =
            sort_for(order, &header, &part).iter().map(|b| by_xi[&xi_of(b)].clone()).collect();
        let p = dir.join(format!("in{i}.bam"));
        write_raw_bam(&p, &header, &sorted);
        concat.extend(sorted);
        paths.push(p);
    }
    (paths, header, concat)
}

/// The oracle: stable sort in `order`; template uses the full key (review B4).
pub fn sort_for(order: SortOrder, header: &Header, records: &[RawRecord]) -> Vec<Vec<u8>> {
    match order {
        SortOrder::TemplateCoordinate => reference_sort_full_template(header, records, None),
        _ => reference_sort(order, header, records, None),
    }
}

/// `fgumi merge` as a subprocess; returns stderr (asserting success).
pub fn merge_cli(inputs: &[PathBuf], output: &Path, order: &str, threads: usize) -> String {
    let out = std::process::Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .env("RUST_LOG", "info")
        .args(["merge", "--order", order, "--threads", &threads.to_string(), "-o"])
        .arg(output)
        .args(inputs)
        .output()
        .expect("spawn fgumi merge");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "fgumi merge failed:\n{stderr}");
    stderr
}

/// A header rendered as SAM text.
pub fn header_text(header: &Header) -> String {
    let mut writer = noodles::sam::io::Writer::new(Vec::new());
    writer.write_header(header).expect("render header");
    String::from_utf8(writer.into_inner()).expect("SAM header text is UTF-8")
}

/// [`header_text`] without its `@PG` lines.
pub fn header_text_without_pg(header: &Header) -> String {
    let mut out = String::new();
    for line in header_text(header).lines().filter(|l| !l.starts_with("@PG")) {
        out.push_str(line);
        out.push('\n');
    }
    out
}
