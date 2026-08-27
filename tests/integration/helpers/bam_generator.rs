//! Utilities for generating test BAM data programmatically.

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::sam::Header;
use noodles::sam::alignment::record_buf::RecordBuf;

/// Rewrites a BAM as uncompressed SAM text.
///
/// Done in-process with noodles rather than by shelling out to `samtools`, so
/// tests do not depend on an external tool being installed. Shared by the
/// SAM-input tests, which need the same records in both encodings to prove a
/// command reads them identically.
pub fn transcode_bam_to_sam(bam_path: &std::path::Path, sam_path: &std::path::Path) {
    use noodles::sam::alignment::io::Write as _;

    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(bam_path).expect("open BAM to transcode"),
    ));
    let header = reader.read_header().expect("read BAM header");

    let mut writer =
        noodles::sam::io::Writer::new(std::fs::File::create(sam_path).expect("create SAM"));
    writer.write_header(&header).expect("write SAM header");
    for result in reader.record_bufs(&header) {
        let record = result.expect("read BAM record");
        writer.write_alignment_record(&header, &record).expect("write SAM record");
    }
}

/// Creates a UMI family that is already grouped: every read carries both `RX`
/// and the same `MI` molecule id.
///
/// `downsample` and the other post-grouping commands reject reads without `MI`.
/// All reads share one position and one molecule, so the result is trivially in
/// template-coordinate order and can also be fed to `merge`, which refuses an
/// input that is not sorted.
pub fn create_grouped_family(
    umi: &str,
    mi: &str,
    depth: usize,
    base_name: &str,
    sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    create_family_with_tags(
        &[(SamTag::RX, umi.as_bytes()), (SamTag::MI, mi.as_bytes())],
        depth,
        base_name,
        sequence,
        quality,
    )
}

/// Creates one duplex-grouped molecule: `depth` read pairs on each strand.
///
/// The duplex consensus caller and both metrics commands key off an `MI` tag
/// carrying a `/A` or `/B` strand suffix, and a molecule seen on only one strand
/// yields no duplex output at all — so the minimum shape that makes those
/// commands emit anything is both strands of one molecule. Reads also carry
/// `RX`, which is what the `umi_counts` table is keyed on; with `MI` alone the
/// metrics commands write an empty UMI table.
///
/// Orientation mirrors a real duplex pair: `/A` is FR (R1 forward at
/// `ref_start`, R2 reverse 100bp downstream) and `/B` is RF. The two strands see
/// the paired UMI in opposite order, so `RX` is swapped between them.
///
/// `test_duplex_command` has its own private builder of the same shape that
/// tags `MI` only. This one lives here because several test modules need it and
/// because it also tags `RX`; folding that one into this would add `RX` to the
/// records those duplex tests assert on, which is a change to them, not to this.
pub fn create_duplex_grouped_family(
    mi: &str,
    depth: usize,
    sequence: &str,
    quality: u8,
    ref_start: i32,
) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(depth * 4);
    for (strand, umi, is_b_strand) in [("A", "AAA-TTT", false), ("B", "TTT-AAA", true)] {
        records.extend(strand_reads(
            mi,
            strand,
            umi,
            depth,
            sequence,
            quality,
            ref_start,
            is_b_strand,
        ));
    }
    records
}

/// Creates a single-strand (`/A`-only) grouped molecule of `depth` read pairs.
///
/// `simplex-metrics` *rejects* input whose base UMI appears on both strands —
/// "received duplex-UMI data … run duplex-metrics" — so it cannot share
/// [`create_duplex_grouped_family`]'s fixture. This is the same shape with the
/// `/B` half left off.
pub fn create_single_strand_family(
    mi: &str,
    depth: usize,
    sequence: &str,
    quality: u8,
    ref_start: i32,
) -> Vec<RawRecord> {
    strand_reads(mi, "A", "AAA-TTT", depth, sequence, quality, ref_start, false)
}

/// `depth` read pairs on one strand of molecule `mi`.
#[expect(clippy::too_many_arguments, reason = "internal builder; every field varies per strand")]
fn strand_reads(
    mi: &str,
    strand: &str,
    umi: &str,
    depth: usize,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    is_b_strand: bool,
) -> Vec<RawRecord> {
    let mut records = Vec::with_capacity(depth * 2);
    for i in 0..depth {
        let (r1, r2) = duplex_read_pair(
            &format!("{mi}{strand}{i}"),
            &format!("{mi}/{strand}"),
            umi,
            sequence,
            quality,
            ref_start,
            is_b_strand,
        );
        records.push(r1);
        records.push(r2);
    }
    records
}

/// One `RX`+`MI`-tagged read pair on the `/A` (FR) or `/B` (RF) strand.
///
/// The strand's orientation is the whole point: the duplex caller pairs an
/// FR template with the RF template at the same coordinates, so getting the
/// flags or the mate positions wrong silently produces two single-strand
/// molecules rather than one duplex.
fn duplex_read_pair(
    name: &str,
    mi: &str,
    umi: &str,
    sequence: &str,
    quality: u8,
    ref_start: i32,
    is_b_strand: bool,
) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let read_len = seq.len();
    let cigar_op = u32::try_from(read_len).expect("read_len fits u32") << 4;
    let span = i32::try_from(read_len + 100).expect("template span fits i32");

    // B strand is the same molecule read the other way round: R1 reverse at the
    // far position, R2 forward at the near one.
    let (r1_start, r2_start, r1_rev, r2_rev) = if is_b_strand {
        (ref_start + 100, ref_start, true, false)
    } else {
        (ref_start, ref_start + 100, false, true)
    };
    let tlen = if is_b_strand { -span } else { span };

    let build = |first: bool| {
        let (pos, mate_pos, rev, mate_rev) = if first {
            (r1_start, r2_start, r1_rev, r2_rev)
        } else {
            (r2_start, r1_start, r2_rev, r1_rev)
        };
        let segment = if first { flags::FIRST_SEGMENT } else { flags::LAST_SEGMENT };
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .sequence(seq)
            .qualities(&vec![quality; read_len])
            .flags(
                flags::PAIRED
                    | segment
                    | if rev { flags::REVERSE } else { 0 }
                    | if mate_rev { flags::MATE_REVERSE } else { 0 },
            )
            .ref_id(0)
            .pos(pos - 1)
            .mapq(60)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(mate_pos - 1)
            .template_length(if first { tlen } else { -tlen })
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::MI, mi.as_bytes());
        b.build()
    };

    (build(true), build(false))
}

/// Creates mapped consensus records carrying the per-read and per-base consensus
/// tags (`cD`/`cE` and `cd`/`ce`).
///
/// `filter` rejects any read without consensus-calling tags, so this is the
/// minimum shape that command accepts. Records are placed at increasing
/// positions on reference 0 so they are also coordinate-ordered.
pub fn create_consensus_family(depth: usize, base_name: &str, sequence: &str) -> Vec<RawRecord> {
    let seq = sequence.as_bytes();
    let cigar_op = u32::try_from(seq.len()).expect("seq.len() fits u32") << 4; // M
    (0..depth)
        .map(|i| {
            let name = format!("{base_name}_{i}");
            let pos = i32::try_from(99 + i * 100).expect("position fits i32");
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[cigar_op])
                .sequence(seq)
                .qualities(&vec![35; seq.len()])
                .add_int_tag(SamTag::CD, 10)
                .add_float_tag(SamTag::CE, 0.0_f32)
                .add_array_u16(SamTag::CD_BASES, &vec![10; seq.len()])
                .add_array_u16(SamTag::CE_BASES, &vec![0; seq.len()]);
            b.build()
        })
        .collect()
}

/// Convert a `RawRecord` to a noodles `RecordBuf` using the default (empty) header.
///
/// Used to bridge the raw-byte builder API with the noodles BAM writer for test
/// file creation.
pub fn to_record_buf(raw: &RawRecord) -> RecordBuf {
    fgumi_raw_bam::raw_record_to_record_buf(raw, &noodles::sam::Header::default())
        .expect("raw_record_to_record_buf should succeed in test")
}

/// Records whose UMI lives in the last `:`-delimited field of the read name and
/// which carry *no* `RX` tag — the input shape `copy-umi` consumes.
///
/// Two mapped reads at one position with Illumina-shaped names ending in a UMI
/// field, so `copy-umi` writes a distinct `RX` per read (and an empty input
/// yields a header-only file, keeping the output-depends-on-input axis honest).
pub fn create_read_name_umi_records() -> Vec<RawRecord> {
    ["inst:1:FC:1:1101:5:1:ACGT", "inst:1:FC:1:1101:5:2:TTGG"]
        .iter()
        .map(|name| {
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8])
                .ref_id(0)
                .pos(100)
                .mapq(60)
                .cigar_ops(&[8 << 4]);
            b.build()
        })
        .collect()
}

/// Creates a UMI family with specified parameters.
///
/// All reads in the family are mapped to reference sequence 0 at position 100
/// with a simple match CIGAR. This ensures they pass the group command's
/// unmapped filter.
///
/// # Arguments
///
/// * `umi` - The UMI sequence to assign
/// * `depth` - Number of reads in the family
/// * `base_name` - Base name for reads (will be suffixed with index)
/// * `sequence` - The read sequence (all reads will have this sequence)
/// * `quality` - Base quality score for all bases
///
/// # Returns
///
/// Vector of `RawRecord` representing the UMI family
pub fn create_umi_family(
    umi: &str,
    depth: usize,
    base_name: &str,
    sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    create_family_with_tags(&[(SamTag::RX, umi.as_bytes())], depth, base_name, sequence, quality)
}

/// Like [`create_umi_family`] but tags each read with an arbitrary `SamTag`
/// (e.g. `SamTag::RX` for UMI mode or `SamTag::BC` for barcode mode).
pub fn create_family_with_tag(
    tag: SamTag,
    umi: &str,
    depth: usize,
    base_name: &str,
    sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    create_family_with_tags(&[(tag, umi.as_bytes())], depth, base_name, sequence, quality)
}

/// Creates a mapped family where every read carries each of `tags`.
///
/// The shared builder behind [`create_umi_family`], [`create_family_with_tag`]
/// and [`create_grouped_family`] — all reads sit at one position on reference 0
/// with a full-length match CIGAR, differing only in which aux tags they carry.
pub fn create_family_with_tags(
    tags: &[(SamTag, &[u8])],
    depth: usize,
    base_name: &str,
    sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    let seq = sequence.as_bytes();
    let cigar_op = u32::try_from(seq.len()).expect("seq.len() fits u32") << 4; // M
    (0..depth)
        .map(|i| {
            let name = format!("{base_name}_{i}");
            let mut b = SamBuilder::new();
            b.read_name(name.as_bytes())
                .ref_id(0)
                .pos(99)
                .mapq(60)
                .flags(0)
                .cigar_ops(&[cigar_op])
                .sequence(seq)
                .qualities(&vec![quality; seq.len()]);
            for (tag, value) in tags {
                b.add_string_tag(*tag, value);
            }
            b.build()
        })
        .collect()
}

/// Creates a paired-end UMI family.
///
/// All reads are mapped to reference sequence 0. R1 is mapped at position 100,
/// R2 at position 200 (to simulate insert size).
///
/// # Arguments
///
/// * `umi` - The UMI sequence (for paired UMIs, use "AAAA-CCCC" format)
/// * `depth` - Number of read pairs in the family
/// * `base_name` - Base name for reads
/// * `r1_sequence` - R1 sequence
/// * `r2_sequence` - R2 sequence
/// * `quality` - Base quality score
///
/// # Returns
///
/// Vector of `RawRecord` with R1 and R2 reads properly flagged
pub fn create_paired_umi_family(
    umi: &str,
    depth: usize,
    base_name: &str,
    r1_sequence: &str,
    r2_sequence: &str,
    quality: u8,
) -> Vec<RawRecord> {
    let r1_seq = r1_sequence.as_bytes();
    let r2_seq = r2_sequence.as_bytes();
    let r1_cigar = u32::try_from(r1_seq.len()).expect("r1_seq.len() fits u32") << 4;
    let r2_cigar = u32::try_from(r2_seq.len()).expect("r2_seq.len() fits u32") << 4;

    // R1 at pos 99 (0-based); R2 at pos 199. Template spans from R1 start through end
    // of R2 — 100bp gap + R2 length — and R2 carries the negated length.
    let template_len = i32::try_from(100 + r2_seq.len()).expect("template length fits i32");

    let mut records = Vec::new();

    for i in 0..depth {
        let read_name = format!("{base_name}_{i}");

        // R1: paired + first segment
        let mut b1 = SamBuilder::new();
        b1.read_name(read_name.as_bytes())
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(199)
            .template_length(template_len)
            .cigar_ops(&[r1_cigar])
            .sequence(r1_seq)
            .qualities(&vec![quality; r1_seq.len()])
            .add_string_tag(SamTag::RX, umi.as_bytes());
        records.push(b1.build());

        // R2: paired + last segment
        let mut b2 = SamBuilder::new();
        b2.read_name(read_name.as_bytes())
            .ref_id(0)
            .pos(199)
            .mapq(60)
            .flags(flags::PAIRED | flags::LAST_SEGMENT)
            .mate_ref_id(0)
            .mate_pos(99)
            .template_length(-template_len)
            .cigar_ops(&[r2_cigar])
            .sequence(r2_seq)
            .qualities(&vec![quality; r2_seq.len()])
            .add_string_tag(SamTag::RX, umi.as_bytes());
        records.push(b2.build());
    }

    records
}

/// Creates a consensus read with specified metrics.
///
/// # Arguments
///
/// * `name` - Read name
/// * `sequence` - Consensus sequence
/// * `depth_max` - Maximum depth (cD tag)
/// * `depth_min` - Minimum depth (cM tag)
/// * `error_rate` - Error rate (cE tag)
/// * `mean_quality` - Mean base quality score
///
/// # Returns
///
/// `RawRecord` representing a consensus read
pub fn create_consensus_read(
    name: &str,
    sequence: &str,
    depth_max: i32,
    depth_min: i32,
    error_rate: f32,
    mean_quality: u8,
) -> RawRecord {
    let seq = sequence.as_bytes();
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(seq)
        .qualities(&vec![mean_quality; seq.len()])
        .add_int_tag(SamTag::CD, depth_max)
        .add_int_tag(SamTag::CM, depth_min)
        .add_float_tag(SamTag::CE, error_rate);
    b.build()
}

/// Builds a SAM header with the given header tags and one reference sequence.
///
/// This is a shared helper used by the public header constructors to avoid
/// duplicating the noodles header building boilerplate.
fn build_header_with_tags(ref_name: &str, ref_len: usize, tags: &[([u8; 2], &str)]) -> Header {
    use bstr::BString;
    use noodles::sam::header::record::value::map::Map as HeaderRecordMap;
    use noodles::sam::header::record::value::map::header::tag::Tag as HeaderTag;
    use noodles::sam::header::record::value::{
        Map, map::Header as HeaderRecord, map::ReferenceSequence,
    };
    use std::num::NonZeroUsize;

    let mut builder = HeaderRecordMap::<HeaderRecord>::builder();
    for &(tag_bytes, value) in tags {
        let HeaderTag::Other(tag) = HeaderTag::from(tag_bytes) else { unreachable!() };
        builder = builder.insert(tag, value);
    }
    let header_map = builder.build().expect("valid header map");

    let reference_sequence = Map::<ReferenceSequence>::new(
        NonZeroUsize::new(ref_len).expect("reference length must be non-zero"),
    );

    Header::builder()
        .set_header(header_map)
        .add_reference_sequence(BString::from(ref_name), reference_sequence)
        .build()
}

/// Creates a minimal SAM header with one reference sequence.
///
/// The header is configured with template-coordinate sort order tags
/// (SO:unsorted, GO:query, SS:template-coordinate) to be compatible with
/// the group command.
///
/// # Arguments
///
/// * `ref_name` - Reference sequence name (e.g., "chr1")
/// * `ref_len` - Reference sequence length
///
/// # Returns
///
/// Configured SAM Header with template-coordinate sort order
pub fn create_minimal_header(ref_name: &str, ref_len: usize) -> Header {
    build_header_with_tags(
        ref_name,
        ref_len,
        &[(*b"SO", "unsorted"), (*b"GO", "query"), (*b"SS", "template-coordinate")],
    )
}

/// Creates a coordinate-sorted SAM header with one reference sequence.
///
/// # Arguments
///
/// * `ref_name` - Reference sequence name (e.g., "chr1")
/// * `ref_len` - Reference sequence length
///
/// # Returns
///
/// Configured SAM Header with SO:coordinate sort order
pub fn create_coordinate_sorted_header(ref_name: &str, ref_len: usize) -> Header {
    build_header_with_tags(ref_name, ref_len, &[(*b"SO", "coordinate")])
}

/// Creates a UMI family with intentional sequencing errors.
///
/// # Arguments
///
/// * `base_umi` - The "true" UMI sequence
/// * `error_umi` - The error variant of the UMI
/// * `base_depth` - Number of reads with `base_umi`
/// * `error_depth` - Number of reads with `error_umi`
/// * `sequence` - Read sequence
///
/// # Returns
///
/// Combined vector of reads with both UMI variants
pub fn create_umi_family_with_errors(
    base_umi: &str,
    error_umi: &str,
    base_depth: usize,
    error_depth: usize,
    sequence: &str,
) -> Vec<RawRecord> {
    let mut records = Vec::new();

    // Add base UMI reads
    records.extend(create_umi_family(
        base_umi,
        base_depth,
        &format!("base_{base_umi}"),
        sequence,
        30,
    ));

    // Add error UMI reads
    records.extend(create_umi_family(
        error_umi,
        error_depth,
        &format!("error_{error_umi}"),
        sequence,
        30,
    ));

    records
}

/// Create a reference FASTA + FAI + sequence dictionary that matches the test header
/// (chr1, 10000bp).
///
/// Writes `ref.fa`, `ref.fa.fai`, and `ref.dict` into the given directory.
/// Returns the path to `ref.fa`.
pub fn create_test_reference(dir: &std::path::Path) -> std::path::PathBuf {
    use std::io::Write;

    let ref_path = dir.join("ref.fa");
    let fai_path = dir.join("ref.fa.fai");
    let dict_path = dir.join("ref.dict");

    let ref_seq = "ACGTACGT".repeat(1250); // 10000 bases
    let mut fasta = std::fs::File::create(&ref_path).unwrap();
    writeln!(fasta, ">chr1").unwrap();
    writeln!(fasta, "{ref_seq}").unwrap();
    fasta.flush().unwrap();

    let fai_content = "chr1\t10000\t6\t10000\t10001\n";
    std::fs::write(&fai_path, fai_content).unwrap();

    let mut dict = std::fs::File::create(&dict_path).unwrap();
    writeln!(dict, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(dict, "@SQ\tSN:chr1\tLN:10000").unwrap();
    dict.flush().unwrap();

    ref_path
}

// ---------------------------------------------------------------------------
// `compare`-feature-only helpers.
//
// Originally shared by `test_compare_bams.rs` and `test_compare_mutation.rs`
// (previously duplicated in both files) since both build small BAMs from
// `RawRecord` fixtures and run them through the `compare` engines.
// `write_bam` is now used outside those too, so only the compare-specific
// record builders below stay feature-gated.
// ---------------------------------------------------------------------------

use noodles::sam::alignment::io::Write as AlignmentWrite;

/// Writes a BAM file from the given header and records.
///
/// Not gated on `compare`: `test_sam_input.rs` builds its fixtures with this
/// too, and that test runs in every feature configuration.
pub fn write_bam(path: &std::path::Path, header: &Header, records: &[RawRecord]) {
    let mut writer =
        noodles::bam::io::Writer::new(std::fs::File::create(path).expect("create test BAM"));
    writer.write_header(header).expect("write test header");
    for record in records {
        writer.write_alignment_record(header, &to_record_buf(record)).expect("write test record");
    }
    writer.try_finish().expect("finish test BAM");
}

/// Builds a simple mapped record with a given name, position, and MI tag.
#[cfg(feature = "compare")]
pub fn mi_record(name: &[u8], pos: i32, mi: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name)
        .sequence(b"ACGTACGT")
        .qualities(&[30; 8])
        .ref_id(0)
        .pos(pos - 1) // pos is 1-based in tests, BAM uses 0-based
        .mapq(60)
        .add_string_tag(SamTag::MI, mi.as_bytes());
    b.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_umi_family() {
        let family = create_umi_family("ACGTACGT", 5, "test", "AAAA", 30);
        assert_eq!(family.len(), 5);

        for (i, record) in family.iter().enumerate() {
            let buf = to_record_buf(record);
            assert_eq!(
                buf.name().map(std::convert::AsRef::as_ref),
                Some(format!("test_{i}").as_bytes())
            );
            assert_eq!(buf.sequence().as_ref(), b"AAAA");
        }
    }

    #[test]
    fn test_create_paired_umi_family() {
        let family = create_paired_umi_family("ACGT-TGCA", 3, "pair", "AAAA", "TTTT", 30);

        // Should have 6 records (3 pairs)
        assert_eq!(family.len(), 6);

        // Check R1/R2 flags
        let r1_flags = family[0].flags();
        let r2_flags = family[1].flags();
        assert_ne!(r1_flags & flags::FIRST_SEGMENT, 0, "R1 should have FIRST_SEGMENT flag");
        assert_eq!(r2_flags & flags::FIRST_SEGMENT, 0, "R2 should not have FIRST_SEGMENT flag");
        assert_ne!(r2_flags & flags::LAST_SEGMENT, 0, "R2 should have LAST_SEGMENT flag");
    }

    #[test]
    fn test_create_consensus_read() {
        let consensus = create_consensus_read("cons1", "ACGT", 10, 5, 0.01, 35);
        let buf = to_record_buf(&consensus);

        assert_eq!(buf.name().map(std::convert::AsRef::as_ref), Some(b"cons1".as_ref()));
        assert_eq!(buf.sequence().as_ref(), b"ACGT");
        assert_eq!(buf.quality_scores().as_ref(), &[35, 35, 35, 35]);
    }

    #[test]
    fn test_create_minimal_header() {
        let header = create_minimal_header("chr1", 1000);
        assert_eq!(header.reference_sequences().len(), 1);
    }
}
