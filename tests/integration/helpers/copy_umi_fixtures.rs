//! Shared fixtures for the `copy-umi` command's integration tests
//! (`test_copy_umi_command.rs`, `test_copy_umi_cutover_parity.rs`).

use std::path::Path;

use fgumi_lib::commands::copy_umi::CopyUmiMetric;
use fgumi_raw_bam::{RawRecord, SamBuilder};

use super::bam_generator::{create_minimal_header, transcode_bam_to_sam, write_bam};

/// A minimal mapped record carrying the given read name (sequence/quals are fixed
/// filler — `copy-umi` never touches them).
pub fn record_named(name: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[4 << 4]);
    b.build()
}

/// As [`record_named`], but pre-populated with an `RX` tag to exercise the
/// overwrite / `--fail-if-tag-present` policy.
pub fn record_named_with_rx(name: &str, rx: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[4 << 4])
        .add_string_tag(fgumi_lib::sam::SamTag::RX, rx.as_bytes());
    b.build()
}

/// Write `records` to a BAM in `dir` and return its path.
pub fn write_input(dir: &Path, records: &[RawRecord]) -> std::path::PathBuf {
    let path = dir.join("in.bam");
    let header = create_minimal_header("chr1", 10_000);
    write_bam(&path, &header, records);
    path
}

/// Read (`name`, `RX`) for every record in a BAM by transcoding to SAM text.
pub fn read_name_and_rx(dir: &Path, bam: &Path) -> Vec<(String, Option<String>)> {
    let sam = dir.join("read_name_and_rx.sam");
    transcode_bam_to_sam(bam, &sam);
    std::fs::read_to_string(&sam)
        .expect("read sam")
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let name = fields[0].to_string();
            let rx = fields.iter().find_map(|f| f.strip_prefix("RX:Z:").map(str::to_string));
            (name, rx)
        })
        .collect()
}

/// Reads the single-row `--metrics` TSV `copy-umi` writes, as a typed
/// [`CopyUmiMetric`] (via `fgumi_metrics::read_metrics_auto`), rather than a
/// stringly `HashMap`: a typed reader can't silently accept a column rename or
/// a reordered header row the way a hand-parsed `header[i] -> value[i]` zip can.
///
/// # Panics
/// Panics if the file cannot be read/parsed as a `CopyUmiMetric` TSV, or does
/// not contain exactly one data row (`copy-umi --metrics` always writes one).
pub fn read_copy_umi_metrics(path: &Path) -> CopyUmiMetric {
    let mut rows = fgumi_metrics::read_metrics_auto::<_, CopyUmiMetric>(path)
        .unwrap_or_else(|e| panic!("read copy-umi metrics {}: {e}", path.display()));
    assert_eq!(
        rows.len(),
        1,
        "copy-umi --metrics TSV at {} must contain exactly one data row, got {}",
        path.display(),
        rows.len()
    );
    rows.remove(0)
}
