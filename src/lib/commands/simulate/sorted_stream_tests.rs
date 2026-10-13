//! Tests of `run_simulated_sort`, the chain sort the simulate commands stream
//! their records through: failure handling (a generator error, a sort failure
//! while the generator is blocked, an aborted generator) and, through the
//! test-only tee, each subcommand's output against the independent reference
//! sort of the records it generated.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::Duration;

use clap::Parser;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use tempfile::TempDir;

use super::common::{TEE_FOR_TEST, run_simulated_sort};
use super::{GroupedReads, MappedReads};
use crate::commands::command::Command;
use crate::commands::common::MemoryLimit;
use crate::commands::sort::{SortOptions, SortOrderArg};
use crate::sam::SamTag;

/// A one-contig header (`chr1`, 1 Mb).
fn header() -> noodles::sam::Header {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    let len = std::num::NonZeroUsize::new(1_000_000).expect("nonzero");
    noodles::sam::Header::builder()
        .add_reference_sequence("chr1", Map::<ReferenceSequence>::new(len))
        .build()
}

/// `n` unsorted single-end records with unique names on `chr1` (every 17th
/// unmapped, strands alternating, a unique `RX` each): the shape of the
/// integration suite's `uniquely_named_records`.
fn records(n: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let mut b = SamBuilder::new();
            b.read_name(format!("q{i}").as_bytes()).sequence(b"ACGTACGTAC").qualities(&[30u8; 10]);
            if i % 17 == 0 {
                b.ref_id(-1).pos(-1).mapq(0).flags(fgumi_raw_bam::flags::UNMAPPED);
            } else {
                let pos = i32::try_from(1 + (i * 7919 % 500) * 10).expect("pos fits i32");
                let strand = if i % 2 == 0 { 0 } else { fgumi_raw_bam::flags::REVERSE };
                b.ref_id(0).pos(pos).mapq(60).flags(strand).cigar_ops(&[10u32 << 4]);
            }
            b.add_string_tag(SamTag::RX, format!("u{i}").as_bytes());
            b.build()
        })
        .collect()
}

/// A template-coordinate sort whose 64 KiB budget spills every few hundred
/// records, with its spills under `tmp`.
fn tiny_sort(tmp: &Path) -> SortOptions {
    SortOptions {
        order: SortOrderArg::TemplateCoordinate,
        max_memory: MemoryLimit::Fixed(64 * 1024),
        memory_per_thread: false,
        tmp_dirs: vec![tmp.to_path_buf()],
        ..SortOptions::default()
    }
}

/// Catches: the generator's `Err` swallowed (a short, "successful" output).
#[test]
fn generator_error_fails_the_pipeline() {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("o.bam");
    let err = run_simulated_sort(
        |tx| {
            tx.send(Ok(records(1_000))).unwrap();
            tx.send(Err(anyhow::anyhow!("generator broke after 1000 records"))).unwrap();
            Ok(())
        },
        &header(),
        &out,
        &tiny_sort(dir.path()),
        2,
        1,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("generator broke after 1000 records"), "{err:#}");
    assert!(!out.exists(), "no output on failure");
}

/// The sort fails (its tmp dir vanishes) while the generator is blocked on the
/// bounded channel. Catches: joining the generator before dropping the
/// pipeline (a deadlock: the receiver stays alive, `send` never fails).
#[test]
fn sort_failure_mid_stream_does_not_hang() {
    let tmp = TempDir::new().unwrap();
    let dir = TempDir::new().unwrap();
    let (out, tmp_path) = (dir.path().join("o.bam"), tmp.path().to_path_buf());
    let sort = tiny_sort(tmp.path());
    let (done_tx, done_rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let r = run_simulated_sort(
            move |tx| {
                tx.send(Ok(records(2_000))).unwrap();
                // Wait for the sort's own temp dir under `tmp`: removing `tmp`
                // before the chain creates it would let the chain re-create it
                // (`create_dir_all`), and the sort would never fail.
                let deadline = std::time::Instant::now() + Duration::from_secs(20);
                while std::fs::read_dir(&tmp_path).unwrap().next().is_none() {
                    assert!(
                        std::time::Instant::now() < deadline,
                        "the sort never made its temp dir"
                    );
                    std::thread::sleep(Duration::from_millis(1));
                }
                // The next spill (or a read of an earlier one) fails.
                std::fs::remove_dir_all(&tmp_path).unwrap();
                for _ in 0..100_000 {
                    if tx.send(Ok(records(2_000))).is_err() {
                        return Ok(()); // the sort hung up
                    }
                }
                panic!("the sort never hung up on the generator");
            },
            &header(),
            &out,
            &sort,
            2,
            1,
        );
        done_tx.send((r.is_err(), out.exists())).unwrap();
    });
    let (failed, out_exists) = match done_rx.recv_timeout(Duration::from_secs(30)) {
        Ok(done) => done,
        Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {
            panic!("simulate hung after a sort failure")
        }
        Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
            panic!("the simulate thread panicked; see its output above")
        }
    };
    assert!(failed, "the sort failure was not reported");
    assert!(!out_exists, "an output was left behind");
}

/// Twin of the retired engine's empty-stream test: a generator that sends
/// nothing yields a valid, header-only output.
#[test]
fn an_empty_generator_writes_a_header_only_output() {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("o.bam");
    let ((), stats) =
        run_simulated_sort(|_tx| Ok(()), &header(), &out, &tiny_sort(dir.path()), 2, 1).unwrap();
    assert_eq!(stats.total_records, 0);
    assert!(read_bam(&out).1.is_empty());
}

/// Twin of the retired engine's first-record-error test: an error before any
/// record fails the sort and leaves no output.
#[test]
fn an_error_before_any_record_fails_the_pipeline() {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("o.bam");
    let err = run_simulated_sort(
        |tx| {
            tx.send(Err(anyhow::anyhow!("generator broke at once"))).unwrap();
            Ok(())
        },
        &header(),
        &out,
        &tiny_sort(dir.path()),
        2,
        1,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("generator broke at once"), "{err:#}");
    assert!(!out.exists(), "no output on failure");
}

/// Twin of the retired engine's producer-error test: the generator itself
/// returns `Err` (not through the channel). Catches: an output left behind.
#[test]
fn aborted_simulate_leaves_no_output() {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("o.bam");
    let err = run_simulated_sort(
        |tx| -> anyhow::Result<()> {
            tx.send(Ok(records(500))).unwrap();
            anyhow::bail!("generator aborted")
        },
        &header(),
        &out,
        &tiny_sort(dir.path()),
        2,
        1,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("generator aborted"), "{err:#}");
    assert!(!out.exists(), "no output on failure");
}

/// A failed simulate publishes nothing over an existing destination: it keeps
/// its bytes. Catches: the sink writing the destination in place.
#[test]
fn a_failed_simulate_leaves_an_existing_destination_as_it_was() {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("o.bam");
    std::fs::write(&out, b"previous contents").unwrap();
    let err = run_simulated_sort(
        |tx| -> anyhow::Result<()> {
            tx.send(Ok(records(500))).unwrap();
            anyhow::bail!("generator aborted")
        },
        &header(),
        &out,
        &tiny_sort(dir.path()),
        2,
        1,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("generator aborted"), "{err:#}");
    assert_eq!(std::fs::read(&out).unwrap(), b"previous contents");
    let staged: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.starts_with(".fgumi-"))
        .collect();
    assert!(staged.is_empty(), "the staged output is removed: {staged:?}");
}

/// A deterministic two-contig reference for the subcommands.
fn write_reference(dir: &Path) -> PathBuf {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut fasta = Vec::new();
    for (name, salt) in [("chr1", 1u64), ("chr2", 7u64)] {
        fasta.extend_from_slice(format!(">{name}\n").as_bytes());
        let mut state = salt;
        for _ in 0..100_000 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            fasta.push(BASES[usize::try_from((state >> 33) & 3).unwrap()]);
        }
        fasta.push(b'\n');
    }
    let path = dir.join("ref.fa");
    std::fs::write(&path, fasta).unwrap();
    path
}

/// The records of a BAM, with its header.
fn read_bam(path: &Path) -> (noodles::sam::Header, Vec<RecordBuf>) {
    let mut reader = noodles::bam::io::reader::Builder.build_from_path(path).unwrap();
    let header = reader.read_header().unwrap();
    let records = reader.record_bufs(&header).map(Result::unwrap).collect();
    (header, records)
}

/// MI bases renumbered from 1 by first appearance in file order (as
/// `grouped-reads` mints them); `/A`, `/B` suffixes kept.
fn renumber_mi_in_order(records: &mut [RecordBuf]) {
    let mut map: HashMap<String, usize> = HashMap::new();
    for r in records {
        let Some(Value::String(mi)) = r.data().get(&Tag::from(SamTag::MI)).cloned() else {
            continue;
        };
        let mi = mi.to_string();
        let (base, suffix) = mi.split_once('/').map_or((mi.as_str(), None), |(b, s)| (b, Some(s)));
        let next = map.len() + 1;
        let new = *map.entry(base.to_string()).or_insert(next);
        let v = suffix.map_or_else(|| new.to_string(), |s| format!("{new}/{s}"));
        r.data_mut().insert(Tag::from(SamTag::MI), Value::String(v.into()));
    }
}

/// Each subcommand's output is the independent reference stable sort
/// (`fgumi_pipeline_io::sort::testing::reference_sort`, not the chain under
/// test) of the records it generated — captured unsorted by the test-only tee
/// — plus, for `grouped-reads`, its MI renumbering in sorted order, with a
/// spilling budget. Catches: wrong key types wired into the chain, dropped or
/// reordered records, MI renumbering applied before the sort.
#[rstest::rstest]
#[case::grouped_reads("grouped-reads", true)]
#[case::mapped_reads("mapped-reads", false)]
fn simulate_output_is_the_reference_sort_of_its_records(#[case] sub: &str, #[case] renumber: bool) {
    let dir = TempDir::new().unwrap();
    let (tee, out, truth) =
        (dir.path().join("tee.bam"), dir.path().join("out.bam"), dir.path().join("t.tsv"));
    let reference = write_reference(dir.path());
    let args = [
        sub,
        "-o",
        out.to_str().unwrap(),
        "--truth",
        truth.to_str().unwrap(),
        "--reference",
        reference.to_str().unwrap(),
        "--num-molecules",
        "1500",
        "--seed",
        "42",
        "--max-memory",
        "1M",
        "--threads",
        "2",
        "--tmp-dir",
        dir.path().to_str().unwrap(),
    ];
    TEE_FOR_TEST.with(|t| *t.borrow_mut() = Some(tee.clone()));
    let run = if sub == "grouped-reads" {
        GroupedReads::try_parse_from(args).unwrap().execute("test")
    } else {
        MappedReads::try_parse_from(args).unwrap().execute("test")
    };
    TEE_FOR_TEST.with(|t| *t.borrow_mut() = None);
    run.unwrap();

    let (mut reader, tee_header) = fgumi_bam_io::create_raw_bam_reader(&tee, 1).unwrap();
    let mut tee_raw = Vec::new();
    let mut record = RawRecord::default();
    while reader.read_record(&mut record).unwrap() != 0 {
        tee_raw.push(record.clone());
    }
    assert!(!tee_raw.is_empty(), "the tee is empty");
    let sorted = fgumi_pipeline_io::sort::testing::reference_sort(
        fgumi_sort::SortOrder::TemplateCoordinate,
        &tee_header,
        &tee_raw,
        Some(SamTag::CB),
    );
    let mut want: Vec<RecordBuf> = sorted
        .iter()
        .map(|r| {
            fgumi_raw_bam::raw_record_to_record_buf(&RawRecord::from(r.clone()), &tee_header)
                .unwrap()
        })
        .collect();
    if renumber {
        renumber_mi_in_order(&mut want);
    }
    assert_eq!(read_bam(&out).1, want);
}
