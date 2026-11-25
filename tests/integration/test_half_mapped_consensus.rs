//! Consensus calling of molecules that have both fully mapped and half-mapped pairs.

use std::ffi::OsStr;
use std::fmt::Write as _;
use std::fs;

use fgumi_lib::sam::SamTag;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::cli::{run_args_and_capture_logs, utf8_path as p};
use crate::helpers::read_bam_output;

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate\n\
                      @SQ\tSN:chr1\tLN:1000000\n\
                      @RG\tID:A\tSM:s\tLB:l1\n";

/// The reads of a pair: both mapped, or only the lower or upper one.
#[derive(Clone, Copy, PartialEq)]
enum Ends {
    Both,
    Lower,
}

/// A pair in molecule 0, whose R1 is the lower read if it is of strand A.
fn pair(name: &str, strand: &str, umi: &str, ends: Ends) -> String {
    let mut sam = String::new();
    for (segment, is_lower) in [(0x40, strand == "A"), (0x80, strand == "B")] {
        let is_mapped = |lower: bool| match ends {
            Ends::Both => true,
            Ends::Lower => lower,
        };
        let (own, mate) = (is_mapped(is_lower), is_mapped(!is_lower));
        let pos = |lower: bool| if lower { 501 } else { 701 };
        let mut flags = 0x1 | segment;
        for (is_set, bit) in
            [(!own, 0x4), (!mate, 0x8), (own && !is_lower, 0x10), (mate && is_lower, 0x20)]
        {
            flags |= if is_set { bit } else { 0 };
        }
        let place = pos(if own { is_lower } else { !is_lower });
        let mate_place = pos(if mate { !is_lower } else { is_lower });
        let (mapq, cigar) = if own { (60, "10M") } else { (0, "*") };
        let (bases, quals) = ("ACGTACGTAC", "IIIIIIIIII");
        write!(
            sam,
            "{name}\t{flags}\tchr1\t{place}\t{mapq}\t{cigar}\t=\t{mate_place}\t0\t{bases}\t{quals}\
             \tRG:Z:A\tRX:Z:{umi}\tMI:Z:0/{strand}"
        )
        .unwrap();
        sam.push_str(if mate { "\tMC:Z:10M\tMQ:i:60\n" } else { "\n" });
    }
    sam
}

fn call(command: &str, pairs: &[String], args: &[&str]) -> Vec<RecordBuf> {
    let dir = TempDir::new().unwrap();
    let (input, output) = (dir.path().join("in.sam"), dir.path().join("out.bam"));
    fs::write(&input, HEADER.to_string() + &pairs.concat()).unwrap();
    let args = [&[command, "-i", p(&input), "-o", p(&output)], args].concat();
    run_args_and_capture_logs(&args.iter().map(OsStr::new).collect::<Vec<_>>());
    read_bam_output(&output).1
}

/// The consensus UMI does not use the UMI segment of a read that did not map.
#[rstest]
fn the_consensus_umi_is_that_of_the_pairs(#[values("simplex", "duplex")] command: &str) {
    let pairs = [
        pair("full", "A", "AAAA-CCCC", Ends::Both),
        pair("half1", "A", "AAAA-GTGT", Ends::Lower),
        pair("half2", "A", "AAAA-TATA", Ends::Lower),
    ];
    let min_reads = if command == "duplex" { "1,1,0" } else { "1" };
    for record in call(command, &pairs, &["--min-reads", min_reads]) {
        match record.data().get(&SamTag::RX.to_noodles_tag()) {
            Some(Value::String(umi)) => assert_eq!(umi.to_string(), "AAAA-CCCC"),
            other => panic!("no UMI: {other:?}"),
        }
    }
}
