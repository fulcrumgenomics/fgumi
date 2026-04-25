//! Shared fixture builders for MI-determinism regression tests.
//!
//! Provides the synthetic paired-end BAM fixture used by both
//! `test_group_determinism` and `test_runall_mi_determinism` to ensure the two
//! test suites exercise the same input data and that any future changes to the
//! fixture automatically apply to all callers.

use std::path::Path;

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Build a paired-end template for `--strategy paired` with explicit orientation.
///
/// `is_rf = false` produces FR orientation: R1 forward at `r1_pos`, R2 reverse
/// at `r2_pos`. `is_rf = true` flips the strand bits for both reads.
pub fn build_paired_template(
    name: &str,
    umi: &str,
    sequence: &str,
    quality: u8,
    r1_pos: i32,
    r2_pos: i32,
    is_rf: bool,
) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let cigar_op = u32::try_from(seq.len()).expect("seq len fits u32") << 4;
    let cigar_str = format!("{}M", seq.len());
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "small synthetic positions fit in i32"
    )]
    let tlen = (r2_pos - r1_pos + seq.len() as i32).abs();

    let r1_flags = flags::PAIRED
        | flags::FIRST_SEGMENT
        | if is_rf { flags::REVERSE } else { 0 }
        | if is_rf { 0 } else { flags::MATE_REVERSE };
    let r2_flags = flags::PAIRED
        | flags::LAST_SEGMENT
        | if is_rf { 0 } else { flags::REVERSE }
        | if is_rf { flags::MATE_REVERSE } else { 0 };

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .ref_id(0)
            .pos(r1_pos - 1)
            .mapq(60)
            .flags(r1_flags)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r2_pos - 1)
            .template_length(if is_rf { -tlen } else { tlen })
            .sequence(seq)
            .qualities(&vec![quality; seq.len()])
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::MC, cigar_str.as_bytes());
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .ref_id(0)
            .pos(r2_pos - 1)
            .mapq(60)
            .flags(r2_flags)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r1_pos - 1)
            .template_length(if is_rf { tlen } else { -tlen })
            .sequence(seq)
            .qualities(&vec![quality; seq.len()])
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::MC, cigar_str.as_bytes());
        b.build()
    };
    (r1, r2)
}

/// Write a BAM containing the given templates against a minimal `chr1` header.
pub fn write_bam(path: &Path, templates: &[(RawRecord, RawRecord)]) {
    use std::fs;
    let header = create_minimal_header("chr1", 100_000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for (r1, r2) in templates {
        writer.write_alignment_record(&header, &to_record_buf(r1)).expect("Failed to write R1");
        writer.write_alignment_record(&header, &to_record_buf(r2)).expect("Failed to write R2");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Build a BAM with paired-end UMI families spread across multiple
/// orientations, multiple UMIs per orientation, and multiple positions.
///
/// The combination of distinct orientations and distinct paired UMIs is what
/// triggers the orientation-subgroup iteration-order bug: the assigner sees
/// more than one orientation subgroup, and within each subgroup more than one
/// UMI family, so `MoleculeId` numbering depends on subgroup iteration order.
///
/// We materialize a large number of distinct molecules per subgroup so that
/// even small variances in `AHashMap::values()` iteration order become
/// observable as different `MoleculeId` numbering across runs.
///
/// Produces: 40 position groups × 256 paired UMIs × 2 orientations
/// = 20480 templates = 40960 reads.
pub fn build_mixed_orientation_bam(path: &Path) {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut umis: Vec<String> = Vec::new();
    for &a in &bases {
        for &b in &bases {
            for &c in &bases {
                for &d in &bases {
                    let prefix = format!("{}{}{}{}", a as char, b as char, c as char, d as char,);
                    let suffix = format!("{}{}{}{}", d as char, c as char, b as char, a as char,);
                    umis.push(format!("{prefix}{prefix}-{suffix}{suffix}"));
                }
            }
        }
    }

    let mut templates = Vec::new();
    let mut counter = 0;
    let mut emit = |name: String, umi: &str, pos: i32, is_rf: bool| {
        templates.push(build_paired_template(
            &name,
            umi,
            "ACGTACGTACGTACGT",
            30,
            pos,
            pos + 100,
            is_rf,
        ));
    };

    // Many position groups, each carrying templates from BOTH orientations.
    // Position-grouping is per-coordinate, so distinct positions form
    // independent subgroup AHashMaps in the assigner — each one a separate
    // chance for iteration order to vary.
    for group_idx in 0..40 {
        let base_pos = 1_000 + group_idx * 1_000;
        for (u_idx, umi) in umis.iter().enumerate() {
            for is_rf in [false, true] {
                counter += 1;
                let name = format!("read_g{group_idx}_u{u_idx}_{is_rf}_{counter}");
                emit(name, umi, base_pos, is_rf);
            }
        }
    }
    write_bam(path, &templates);
}
