//! Integration tests: `fgumi simulate` emits template-coordinate sorted BAMs.
//!
//! These run the **prebuilt** `fgumi` binary (via `CARGO_BIN_EXE_fgumi`) — no
//! recursive `cargo run`, so they are fast and hermetic — and verify the output
//! with `fgumi sort --verify --order template-coordinate`, fgumi's own canonical
//! sort-order checker (the same engine `fgumi sort`/`fgumi group` use, validated
//! against `samtools sort --template-coordinate`). No external `samtools` is
//! required. Gated on the `simulate` feature so the binary actually has the
//! subcommand.
//!
//! This is what catches ordering regressions in `simulate` — e.g. the historical
//! bug where reverse-R1 (F2R1) pairs were mis-positioned because `simulate`
//! computed its own template-coordinate key instead of delegating to `fgumi-sort`.

#![cfg(feature = "simulate")]

use std::path::{Path, PathBuf};
use std::process::Command;

use rstest::rstest;
use tempfile::TempDir;

/// The `fgumi` binary built for this test run (has the `simulate` feature).
const FGUMI: &str = env!("CARGO_BIN_EXE_fgumi");

/// Write a small, deterministic two-contig reference FASTA and return its path.
///
/// The sequence is generated programmatically (a fixed LCG over `ACGT`) so the
/// test needs no committed data file yet has ample varied mapped positions.
fn write_reference(dir: &Path) -> PathBuf {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    const CONTIG_LEN: usize = 100_000;

    let mut fasta = Vec::with_capacity(2 * CONTIG_LEN + 32);
    for (name, salt) in [("chr1", 1u64), ("chr2", 7u64)] {
        fasta.extend_from_slice(format!(">{name}\n").as_bytes());
        let mut state = salt;
        for _ in 0..CONTIG_LEN {
            // SplitMix-ish LCG; deterministic and well-distributed over 2 bits.
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            fasta.push(BASES[((state >> 33) & 3) as usize]);
        }
        fasta.push(b'\n');
    }

    let path = dir.join("ref.fa");
    std::fs::write(&path, fasta).expect("write reference FASTA");
    path
}

/// Run the prebuilt `fgumi` binary and assert it exits successfully.
fn fgumi(args: &[&str]) {
    let status = Command::new(FGUMI).args(args).status().expect("failed to spawn fgumi");
    assert!(status.success(), "`fgumi {}` failed: {status}", args.join(" "));
}

/// Run `fgumi simulate <subcommand>` into `dir`, returning the output BAM path.
///
/// `max_memory` is passed straight through to `--max-memory`, which bounds the
/// streaming sort's in-memory buffer: a large value keeps the whole run in
/// memory, a tiny one forces spill-to-disk plus a k-way merge.
fn simulate(
    dir: &Path,
    subcommand: &str,
    num_molecules: usize,
    duplex: bool,
    max_memory: &str,
    name: &str,
) -> PathBuf {
    let reference = write_reference(dir);
    let bam = dir.join(format!("{name}.bam"));
    let truth = dir.join(format!("{name}.truth.tsv"));

    let num_molecules = num_molecules.to_string();
    let bam_arg = bam.to_str().unwrap().to_string();
    let reference = reference.to_str().unwrap().to_string();
    let truth = truth.to_str().unwrap().to_string();
    // Keep spill files inside the test's temp dir rather than the system temp.
    let tmp_dir = dir.to_str().unwrap().to_string();

    let mut args = vec![
        "simulate",
        subcommand,
        "-o",
        &bam_arg,
        "--truth",
        &truth,
        "--reference",
        &reference,
        "--num-molecules",
        &num_molecules,
        "--seed",
        "42",
        "--max-memory",
        max_memory,
        "--tmp-dir",
        &tmp_dir,
    ];
    if duplex {
        args.push("--duplex");
    }
    fgumi(&args);

    bam
}

/// `simulate` generates in molecule order and streams into the canonical
/// `fgumi-sort` engine, so its output must be template-coordinate sorted for every
/// subcommand and read orientation. `fgumi sort --verify` asserts exactly that
/// (0 sort-order violations), and it is the check that fails if the ordering
/// regresses.
///
/// Each case runs twice — once with a memory budget that holds the whole run
/// (the sort's in-memory fast path) and once with a tiny budget that forces
/// spill-to-disk and a k-way merge — because streaming generation feeds both
/// paths and only the second exercises the spill/merge code.
#[rstest]
#[case::mapped_reads("mapped-reads", 1_000, false)]
#[case::grouped_reads_simplex("grouped-reads", 1_000, false)]
#[case::grouped_reads_duplex("grouped-reads", 1_000, true)]
#[case::grouped_reads_duplex_large("grouped-reads", 10_000, true)]
fn simulate_output_is_template_coordinate_sorted(
    #[case] subcommand: &str,
    #[case] num_molecules: usize,
    #[case] duplex: bool,
    #[values("768M", "1M")] max_memory: &str,
) {
    let dir = TempDir::new().expect("create temp dir");
    let bam = simulate(dir.path(), subcommand, num_molecules, duplex, max_memory, "out");

    // The canonical verifier: the output must already be in template-coordinate
    // order. Exit code is non-zero (and this asserts) if any records are out of
    // order — which is how the previous reverse-R1 mis-ordering was caught.
    fgumi(&["sort", "-i", bam.to_str().unwrap(), "--verify", "--order", "template-coordinate"]);
}

/// Read every alignment record of a BAM, excluding the header.
///
/// The header cannot be part of an output comparison: its `@PG` record embeds
/// the command line, which differs between any two runs invoked with different
/// flags.
fn read_records(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    let mut reader = noodles::bam::io::Reader::new(
        std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display())),
    );
    let header = reader.read_header().expect("read BAM header");
    reader.record_bufs(&header).map(|r| r.expect("read BAM record")).collect()
}

/// The output must not depend on how much memory the sort was given: spilling
/// to disk and merging has to reproduce the in-memory result exactly.
///
/// This is what pins the streaming design — records enter the sort in
/// generation order regardless of budget, so tie-breaking between equal sort
/// keys cannot drift with the spill boundaries.
#[rstest]
#[case::mapped_reads("mapped-reads", 2_000, false)]
#[case::grouped_reads_duplex("grouped-reads", 2_000, true)]
fn simulate_output_is_identical_regardless_of_spilling(
    #[case] subcommand: &str,
    #[case] num_molecules: usize,
    #[case] duplex: bool,
) {
    let dir = TempDir::new().expect("create temp dir");
    let in_memory = simulate(dir.path(), subcommand, num_molecules, duplex, "768M", "in_memory");
    let spilled = simulate(dir.path(), subcommand, num_molecules, duplex, "1M", "spilled");

    assert_eq!(
        read_records(&in_memory),
        read_records(&spilled),
        "spilling changed the simulated output"
    );
}
