//! Phase 3 cross-step integration tests. The block-level chain:
//!
//! ```text
//! ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries → DecodeRecords
//!     → GroupBam → SerializeBamRecords → BgzfCompress → WriteBgzfFile
//! ```
//!
//! Drives the chain via `Pipeline::run` against a synthetic BAM, then
//! reads back the output to verify record-count invariance and
//! ordering. Threads-{1, 4} variants exercise the framework's drain
//! protocol under both Serial-only and Parallel + Serial mixes.

use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::steps::bgzf::compress::BgzfCompress;
use crate::pipeline::steps::bgzf::decompress::BgzfDecompress;
use crate::pipeline::steps::boundaries::bam::FindBamBoundaries;
use crate::pipeline::steps::group::bam::GroupBam;
use crate::pipeline::steps::parse::decode::DecodeRecords;
use crate::pipeline::steps::serialize::SerializeBamRecords;
use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;
use crate::pipeline::steps::source::read_bam::{DEFAULT_BLOCKS_PER_BATCH, read_bam};
use fgumi_bam_io::GroupKeyConfig;

/// Build an in-memory BAM file with `n_records` synthetic records using
/// fgumi-bam-io's writer. Returns the path (a `tempfile::TempPath`) so
/// the file is cleaned up when dropped.
fn build_test_bam(n_records: usize) -> tempfile::TempPath {
    use fgumi_raw_bam::testutil::make_bam_bytes;
    use noodles::sam::Header;

    let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    let header = Header::default();
    let mut writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
    for i in 0..n_records {
        let name_str = format!("read_{i}");
        let record_bytes = make_bam_bytes(
            -1,  // tid (unmapped)
            -1,  // pos
            0x4, // flag (unmapped)
            name_str.as_bytes(),
            &[], // cigar_ops
            0,   // seq_len
            -1,  // mate_tid
            -1,  // mate_pos
            &[], // aux_data
        );
        writer.write_raw_record(&record_bytes).unwrap();
    }
    writer.finish().unwrap();
    path
}

/// Read back a BAM file and return the count of records.
fn count_records(path: &std::path::Path) -> usize {
    use fgumi_bam_io::PipelineReaderOpts;
    use fgumi_raw_bam::RawRecord;

    let (mut reader, _hdr) =
        fgumi_bam_io::create_raw_bam_reader_with_opts(path, 1, PipelineReaderOpts::default())
            .unwrap();
    let mut record = RawRecord::default();
    let mut n = 0;
    loop {
        match reader.read_record(&mut record).unwrap() {
            0 => break,
            _ => n += 1,
        }
    }
    n
}

/// Read every record's bytes into a `Vec<Vec<u8>>` in input order. Used
/// for byte-level + ordering equivalence checks.
fn collect_record_bytes(path: &std::path::Path) -> Vec<Vec<u8>> {
    use fgumi_bam_io::PipelineReaderOpts;
    use fgumi_raw_bam::RawRecord;

    let (mut reader, _hdr) =
        fgumi_bam_io::create_raw_bam_reader_with_opts(path, 1, PipelineReaderOpts::default())
            .unwrap();
    let mut record = RawRecord::default();
    let mut all = Vec::new();
    while reader.read_record(&mut record).unwrap() > 0 {
        all.push(record.as_ref().to_vec());
    }
    all
}

/// Assert input and output have byte-identical records in input order.
/// Catches reorderings, dropped records, and per-record mutations that a
/// pure count check would miss.
fn assert_records_round_trip_in_order(input: &std::path::Path, output: &std::path::Path) {
    let in_recs = collect_record_bytes(input);
    let out_recs = collect_record_bytes(output);
    assert_eq!(
        in_recs.len(),
        out_recs.len(),
        "record count mismatch: input={}, output={}",
        in_recs.len(),
        out_recs.len()
    );
    for (i, (a, b)) in in_recs.iter().zip(out_recs.iter()).enumerate() {
        assert_eq!(a, b, "record {i}: input bytes != output bytes");
    }
}

/// Wire the full block-level chain end-to-end and run it.
fn run_full_chain(input: &std::path::Path, output: &std::path::Path, threads: usize) {
    let opts = fgumi_bam_io::PipelineReaderOpts::default();
    let (read_step, header) =
        read_bam(input, opts, DEFAULT_BLOCKS_PER_BATCH, 4 * 1024 * 1024).unwrap();
    let bytes_4mb: u64 = 4 * 1024 * 1024;

    let builder = Pipeline::builder();
    builder
        .chain(read_step)
        .chain(BgzfDecompress::new(bytes_4mb))
        .chain(FindBamBoundaries::new(bytes_4mb))
        // DecodeRecords now does parse + group-key extraction in one pass
        // (matches legacy bam.rs Decode step); the previous Parse → Decode
        // split was unused outside tests.
        .chain(DecodeRecords::new(GroupKeyConfig::default(), bytes_4mb))
        .chain(GroupBam::new(64, bytes_4mb))
        .chain(SerializeBamRecords::new(bytes_4mb))
        .chain(BgzfCompress::new(1, bytes_4mb))
        .chain(WriteBgzfFile::new(output, &header, 1).unwrap())
        .into_sink_marker();
    let pipeline = builder.build().unwrap();
    pipeline.run(PipelineConfig { threads, ..Default::default() }).unwrap();
}

#[test]
fn full_chain_threads_2() {
    // Source/sink are Serial+Affinity (Reader/Writer): worker 0 reads,
    // worker 1 writes, all workers can do interior Parallel/Serial steps.
    // The chain works at threads >= 1; this test exercises a small case.
    let input = build_test_bam(50);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    run_full_chain(&input, &output, 2);
    let bytes = std::fs::read(&output).unwrap();
    assert!(bytes.len() >= 28, "output BAM must contain at least the BGZF EOF: {}", bytes.len());
    assert_eq!(&bytes[0..2], &[0x1f, 0x8b], "BGZF/gzip magic");
    assert_records_round_trip_in_order(&input, &output);
}

#[test]
fn full_chain_threads_4() {
    // Source/sink are Serial+Affinity (Reader on T0, Writer on TN-1);
    // Boundaries/Group are Serial (any worker eligible). Decode/parse/
    // serialize/compress are Parallel. With threads=4 every worker
    // rotates through whatever has work.
    let input = build_test_bam(200);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    run_full_chain(&input, &output, 4);
    let bytes = std::fs::read(&output).unwrap();
    assert!(bytes.len() >= 28);
    assert_eq!(&bytes[0..2], &[0x1f, 0x8b]);
    assert_records_round_trip_in_order(&input, &output);
}

#[test]
fn full_chain_many_records_spans_many_bgzf_blocks() {
    // 50_000 unmapped records → tens of BGZF blocks. Exercises the
    // boundary state's cross-block carryover and the reorder stages
    // under realistic block counts.
    let input = build_test_bam(50_000);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    run_full_chain(&input, &output, 4);
    assert_records_round_trip_in_order(&input, &output);
}

#[test]
#[ignore = "big — run with --ignored or in CI"]
fn full_chain_500k_records() {
    let input = build_test_bam(500_000);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    run_full_chain(&input, &output, 4);
    assert_eq!(count_records(&output), 500_000);
}

/// Build a BAM with mapped paired-end records (more realistic record
/// sizes — sequence + cigar + tags). 100bp reads, 4bp cigar, RG tag.
fn build_paired_test_bam(n_pairs: usize) -> tempfile::TempPath {
    use fgumi_raw_bam::testutil::make_bam_bytes;
    use noodles::sam::Header;

    let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    let header = Header::default();
    let mut writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
    // Flags: PAIRED | UNMAPPED | MATE_UNMAPPED + (FIRST_SEGMENT | LAST_SEGMENT).
    let r1_flags: u16 = 1 | 4 | 8 | 64;
    let r2_flags: u16 = 1 | 4 | 8 | 128;
    for i in 0..n_pairs {
        let name_str = format!("read_{i:08}");
        let r1 = make_bam_bytes(-1, -1, r1_flags, name_str.as_bytes(), &[], 100, -1, -1, &[]);
        writer.write_raw_record(&r1).unwrap();
        let r2 = make_bam_bytes(-1, -1, r2_flags, name_str.as_bytes(), &[], 100, -1, -1, &[]);
        writer.write_raw_record(&r2).unwrap();
    }
    writer.finish().unwrap();
    path
}

#[test]
fn full_chain_paired_records() {
    // 10K pairs = 20K records, paired by name (GroupBam should produce
    // one Template per pair).
    let input = build_paired_test_bam(10_000);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    run_full_chain(&input, &output, 4);
    assert_records_round_trip_in_order(&input, &output);
}

#[test]
fn full_chain_empty_bam() {
    // Header-only input; the chain should still produce a valid BAM with
    // a header + EOF and zero records.
    let input = build_test_bam(0);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    run_full_chain(&input, &output, 2);
    let bytes = std::fs::read(&output).unwrap();
    assert!(bytes.len() >= 28);
    assert_eq!(count_records(&output), 0);
}

#[test]
fn full_chain_with_process_mutates_mq() {
    // Plan Task 12 Step 4: insert a `Process` mid-step that bumps every
    // record's MAPQ to a known value, run the chain end-to-end, and
    // verify the output reflects the change. This is the regression test
    // for the closure-driven `ProcessOrdered` family + the held-slot /
    // drain contract under a non-trivial Parallel transform.
    use crate::pipeline::steps::process::process_ordered;
    use crate::pipeline::steps::types::DecodedRecordBatch;

    let input = build_test_bam(2_000);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();

    let opts = fgumi_bam_io::PipelineReaderOpts::default();
    let (read_step, header) =
        read_bam(&input, opts, DEFAULT_BLOCKS_PER_BATCH, 4 * 1024 * 1024).unwrap();
    let bytes_4mb: u64 = 4 * 1024 * 1024;

    let bump_mq = process_ordered::<DecodedRecordBatch, DecodedRecordBatch, _>(
        "BumpMq",
        bytes_4mb,
        |batch: DecodedRecordBatch| {
            // BAM record body byte 9 is `mapq` (raw layout: refID(4) +
            // pos(4) + l_read_name(1) + mapq(1) + ...). Mutating mapq
            // doesn't change the record's identity (qname / library /
            // cell barcode), so the pre-computed GroupKey stays valid.
            // Take ownership, edit in place, and rebuild via `new` so the
            // cached `total_bytes` is recomputed for the edited records.
            let serial = batch.batch_serial();
            let mut records = batch.into_records();
            for record in &mut records {
                let buf = record.raw_bytes_mut().as_mut_vec();
                if buf.len() > 9 {
                    buf[9] = 42;
                }
            }
            Ok(DecodedRecordBatch::new(serial, records))
        },
    );

    let builder = Pipeline::builder();
    builder
        .chain(read_step)
        .chain(BgzfDecompress::new(bytes_4mb))
        .chain(FindBamBoundaries::new(bytes_4mb))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), bytes_4mb))
        .chain(bump_mq)
        .chain(GroupBam::new(64, bytes_4mb))
        .chain(SerializeBamRecords::new(bytes_4mb))
        .chain(BgzfCompress::new(1, bytes_4mb))
        .chain(WriteBgzfFile::new(&output, &header, 1).unwrap())
        .into_sink_marker();
    let pipeline = builder.build().unwrap();
    pipeline.run(PipelineConfig { threads: 4, ..Default::default() }).unwrap();

    // Every output record's mapq byte must be 42.
    let recs = collect_record_bytes(&output);
    assert_eq!(recs.len(), 2_000);
    for (i, r) in recs.iter().enumerate() {
        assert!(r.len() > 9, "record {i} too short");
        assert_eq!(r[9], 42, "record {i} mapq not mutated by Process step");
    }
}

#[test]
fn pipeline_dag_renders_for_full_chain() {
    let input = build_test_bam(0);
    let output = tempfile::NamedTempFile::new().unwrap().into_temp_path();
    let opts = fgumi_bam_io::PipelineReaderOpts::default();
    let (read_step, header) =
        read_bam(&input, opts, DEFAULT_BLOCKS_PER_BATCH, 4 * 1024 * 1024).unwrap();
    let bytes_4mb: u64 = 4 * 1024 * 1024;

    let builder = Pipeline::builder();
    builder
        .chain(read_step)
        .chain(BgzfDecompress::new(bytes_4mb))
        .chain(FindBamBoundaries::new(bytes_4mb))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), bytes_4mb))
        .chain(GroupBam::new(64, bytes_4mb))
        .chain(SerializeBamRecords::new(bytes_4mb))
        .chain(BgzfCompress::new(1, bytes_4mb))
        .chain(WriteBgzfFile::new(&output, &header, 1).unwrap())
        .into_sink_marker();
    let pipeline = builder.build().unwrap();

    let dag = pipeline.dag();
    for name in [
        "ReadBgzfBlocks",
        "BgzfDecompress",
        "FindBamBoundaries",
        "DecodeRecords",
        "GroupBam",
        "SerializeBamRecords",
        "BgzfCompress",
        "WriteBgzfFile",
    ] {
        assert!(dag.contains(name), "DAG missing {name}: {dag}");
    }
    assert!(dag.contains("(sink)"), "DAG missing sink marker: {dag}");
}
