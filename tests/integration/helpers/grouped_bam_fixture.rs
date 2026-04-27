//! Small grouped-BAM fixture for tests that start from `runall --start-from group`.
//!
//! Builds three UMI families at a single reference position, tags each record
//! with `RX` (UMI) and `MI` (molecule identifier), and writes them to a BAM
//! file at the given path. Header is a minimal template-coordinate-sorted
//! header with a single `chr_test` reference sequence so `runall
//! --start-from group` accepts the input.
//!
//! The fixture is intentionally small and deterministic: the same inputs and
//! same tag values every call, so it can be re-used as input for byte-compare
//! tests (e.g. determinism guards).

use std::path::Path;

use fgumi_lib::sam::builder::RecordBuilder;
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;

use super::bam_generator::{create_minimal_header, create_umi_family};

/// Build a small grouped BAM at `path`.
///
/// The resulting BAM contains three UMI families (each with 4 single-end
/// reads) tagged with `RX` and `MI`. The header declares `chr_test` (length
/// 1000) and a template-coordinate sort order.
///
/// # Panics
///
/// Panics if the file cannot be created or any record cannot be written.
pub fn build_small_grouped_bam(path: &Path) {
    let header = create_minimal_header("chr_test", 1000);

    let file = std::fs::File::create(path).expect("Failed to create fixture BAM file");
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(&header).expect("Failed to write BAM header");

    // Three UMI families, each identified by a distinct MI value. All records
    // share reference `chr_test` at position 100 (set by `create_umi_family`).
    let families: &[(u32, &str)] = &[(1, "AAAAAAAA"), (2, "CCCCCCCC"), (3, "GGGGGGGG")];

    for &(mi, umi) in families {
        let base_name = format!("fam{mi}");
        let family = create_umi_family(umi, 4, &base_name, "ACGTACGTACGTACGT", 30);
        for record in &family {
            // `create_umi_family` already sets RX. Re-build to also add MI so
            // the fixture is valid as a grouped-BAM input (RX + MI tags).
            let mi_string = mi.to_string();
            let raw_bytes: &[u8] = record.as_ref();
            let name =
                std::str::from_utf8(fgumi_raw_bam::read_name(raw_bytes)).expect("name is ASCII");
            let sequence_bytes = record.sequence_vec();
            let sequence = std::str::from_utf8(&sequence_bytes).expect("sequence is ASCII");
            let qualities: Vec<u8> = fgumi_raw_bam::quality_scores_slice(raw_bytes).to_vec();

            let tagged = RecordBuilder::new()
                .name(name)
                .sequence(sequence)
                .qualities(&qualities)
                .reference_sequence_id(0)
                .alignment_start(100)
                .mapping_quality(60)
                .tag("RX", umi)
                .tag("MI", mi_string.as_str())
                .build();

            writer.write_alignment_record(&header, &tagged).expect("write record");
        }
    }

    writer.try_finish().expect("Failed to finish BAM");

    assert!(path.exists(), "fixture BAM was not written to {}", path.display());
}
