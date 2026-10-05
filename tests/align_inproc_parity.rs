#![cfg(feature = "aligner-bwa-mem3")]
#![deny(unsafe_code)]

//! End-to-end byte-parity + determinism acceptance gate for the in-process
//! bwa-mem3 aligner backend.
//!
//! For a fixed host, bwa-mem3 version, reference index,
//! `--aligner::chunk-size`, and unmapped input, the in-process backend's merged
//! BAM must be byte-identical to the subprocess `bwa-mem3` preset's — for every
//! `--threads`, every `--aligner::sub-batch-templates`, and every
//! `--pool-scheduler` (the default `auto`, `drain-first`, and `chain-order`) —
//! modulo exactly the normalizations listed on [`align_common::normalize`]: the
//! `@PG` `CL` field, a git dev suffix on the `@PG` `VN` field, and the order of
//! each record's optional tags. The comparison is on decoded BAM records, so
//! every core field and every tag's type (including integer width) and value is
//! compared byte for byte. The parity matrix ([`inproc_matches_subprocess`])
//! drives both presets over the same simulated unmapped BAM; the determinism
//! cases assert two in-process runs agree and that the in-process backend is
//! thread-invariant (`--threads 16` == `1`).
//!
//! The fixture (see `align_common`) is built to exercise real alignment
//! behaviour — substitutions, indels, unmapped reads, chimeric reads with
//! supplementary alignments, discordant pairs, and multi-mapping reads in a
//! duplicated segment — and every case asserts the output actually contains
//! those shapes, so the gate cannot degrade into comparing trivial output. Every
//! input record carries an `RX` tag (some are QC-fail), and every case also
//! checks the merge's contract independently of parity: each input read yields
//! one primary record carrying its `RX` and QC-fail state.
//!
//! ## Version-parity is only defined against the vendored commit
//!
//! Byte-parity is defined ONLY against a `bwa-mem3` CLI built from the SAME
//! commit the in-process backend links (the vendored v0.14.0 `5c1d5e39…`, via
//! `bwa-mem3-rs` 0.4.2). A *different* CLI version — e.g. the pixi
//! `bwa-mem3` a dev box happens to carry — will show spurious divergence, so
//! this test must never be pointed at an arbitrary PATH binary and trusted. The
//! reference binary is supplied by the CI e2e-parity job, which builds it from
//! the vendored commit; see `.github/workflows/check.yml`.
//!
//! ## Environment gating (skip-as-pass, unless tools are required)
//!
//! Every case needs one thing: `BWA_MEM3_BIN`, a `bwa-mem3` binary built from the
//! vendored commit. It indexes the fixture reference (so the in-process backend,
//! which links the same commit, can load the index) and runs the subprocess leg.
//! There is no PATH fallback, because a PATH `bwa-mem3` on a dev box is very
//! likely a different version and would show spurious divergence instead of a
//! clear skip. When it is absent the test prints a skip line and returns (a
//! pass); setting `FGUMI_BWA_MEM3_REQUIRE_TOOLS` flips that to a hard failure, so
//! the dedicated e2e-parity job cannot silently degrade into a no-op.

mod align_common;

use std::path::Path;
use std::process::Command;

use align_common::{
    Bam, Input, NormBam, TagOrder, Template, assert_fixture_coverage, assert_reads_preserved,
    assert_same_bam, assert_tags_transferred, bwa_mem3_bin_from_env, mid_pair_chunk, normalize,
    normalized_pg_version, read_bam, simulate_templates, split_names, tools_required,
    without_git_suffix, write_fixture_reference, write_unmapped_bam,
};
use rstest::rstest;
use tempfile::TempDir;

/// Templates simulated per input. Small enough that the full matrix runs fast,
/// large enough that every read kind appears several times per template shape
/// and the small `--aligner::chunk-size` values below split the input into
/// several / many cohorts.
const N_TEMPLATES: usize = 64;

/// Resolve `BWA_MEM3_BIN`, or skip (panic under `FGUMI_BWA_MEM3_REQUIRE_TOOLS`)
/// with a line naming the case.
macro_rules! bin_or_skip {
    ($label:expr) => {
        match bwa_mem3_bin_from_env() {
            Some(bin) => bin,
            None => {
                let reason = "no reference bwa-mem3 binary (set BWA_MEM3_BIN)";
                assert!(!tools_required(), "{}: required tool missing: {reason}", $label);
                eprintln!("skip {}: {reason}", $label);
                return;
            }
        }
    };
}

/// A prepared case: temp dir, indexed fixture reference, simulated templates and
/// the unmapped input BAM written from them.
struct Case {
    dir: TempDir,
    reference: std::path::PathBuf,
    unmapped: std::path::PathBuf,
    templates: Vec<Template>,
    /// `runall --methylation-mode` for both legs, or `None`.
    methylation: Option<&'static str>,
}

impl Case {
    fn new(bin: &str, input: Input) -> Self {
        let dir = TempDir::new().expect("create temp dir");
        let fixture = write_fixture_reference(dir.path(), bin);
        let templates = simulate_templates(&fixture.seq, input, N_TEMPLATES);
        let unmapped = dir.path().join("unmapped.bam");
        write_unmapped_bam(&unmapped, &templates);
        Self { reference: fixture.fasta, unmapped, templates, dir, methylation: None }
    }

    fn out(&self, file: &str) -> std::path::PathBuf {
        self.dir.path().join(file)
    }
}

// ---------------------------------------------------------------------------
// Running runall
// ---------------------------------------------------------------------------

/// The `--pool-scheduler` value under test.
#[derive(Debug, Clone, Copy)]
enum Scheduler {
    /// Omit `--pool-scheduler`: the production default (`auto`), which is the
    /// only setting that installs `RefillDrainScheduler` for the in-process
    /// backend.
    Auto,
    DrainFirst,
    ChainOrder,
}
impl Scheduler {
    /// The `--pool-scheduler` value to pass, or `None` to leave it at its
    /// default.
    fn flag(self) -> Option<&'static str> {
        match self {
            Self::Auto => None,
            Self::DrainFirst => Some("drain-first"),
            Self::ChainOrder => Some("chain-order"),
        }
    }

    /// A label for test messages.
    fn label(self) -> &'static str {
        self.flag().unwrap_or("auto")
    }
}

/// The `--aligner::chunk-size` under test, chosen so a tiny input still spans
/// one, several, and many cohorts.
///
/// Production-scale values (`200k`, `5M`) would need a multi-megabase input to
/// split at all, so these are small base counts that give the `N_TEMPLATES`
/// input the same one/several/many cohort structure. Both legs use the identical
/// value, so parity is preserved regardless of the number.
#[derive(Debug, Clone, Copy)]
enum Chunk {
    /// The pipeline default (`-K` unset): a single cohort.
    One,
    /// Several cohorts.
    Several,
    /// Many cohorts.
    Many,
    /// A value computed from the input ([`mid_pair_chunk`]) so the first `-K`
    /// cut lands between a pair's two reads, forcing bwa-mem3's mid-pair split.
    /// Only meaningful for [`Input::Mixed`]: paired-only input never splits.
    MidPair,
}
impl Chunk {
    /// `None` leaves `--aligner::chunk-size` at its default; `Some(v)` passes it.
    fn value(self, templates: &[Template]) -> Option<u64> {
        match self {
            Self::One => None,
            Self::Several => Some(5000),
            Self::Many => Some(1000),
            Self::MidPair => Some(
                mid_pair_chunk(templates)
                    .expect("Chunk::MidPair needs an R1 at an even running read count")
                    .chunk_size,
            ),
        }
    }
}

/// Which alignment backend a leg runs.
#[derive(Debug, Clone, Copy)]
enum Leg<'a> {
    /// Subprocess `bwa-mem3` preset, using the reference binary via
    /// `--aligner-bin`.
    Subprocess { reference_bin: &'a str },
    /// In-process `bwa-mem3-inproc` preset with the given sub-batch size, and
    /// `--aligner::dedup-reads <dedup>` when set (else its default, `on`).
    InProcess { sub_batch_templates: usize, dedup: Option<&'a str> },
}

/// The `runall --start-from align --stop-after zipper` command for one leg of
/// `case`, writing `out`.
fn leg_command(
    case: &Case,
    out: &Path,
    leg: Leg<'_>,
    threads: usize,
    chunk: Chunk,
    scheduler: Scheduler,
) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.arg("runall")
        .args(["--start-from", "align"])
        .args(["--stop-after", "zipper"])
        .args(["--threads", &threads.to_string()])
        .arg("--input")
        .arg(&case.unmapped)
        .arg("--ref")
        .arg(&case.reference)
        .arg("--output")
        .arg(out);
    if let Some(flag) = scheduler.flag() {
        cmd.args(["--pool-scheduler", flag]);
    }
    if let Some(k) = chunk.value(&case.templates) {
        cmd.args(["--aligner::chunk-size", &k.to_string()]);
    }
    if let Some(mode) = case.methylation {
        cmd.args(["--methylation-mode", mode]);
    }
    match leg {
        Leg::Subprocess { reference_bin } => {
            cmd.args(["--aligner::preset", "bwa-mem3"]);
            cmd.arg("--aligner-bin").arg(reference_bin);
        }
        Leg::InProcess { sub_batch_templates, dedup } => {
            cmd.args(["--aligner::preset", "bwa-mem3-inproc"]);
            cmd.args(["--aligner::sub-batch-templates", &sub_batch_templates.to_string()]);
            if let Some(dedup) = dedup {
                cmd.args(["--aligner::dedup-reads", dedup]);
            }
        }
    }
    cmd
}

/// Run one leg (see [`leg_command`]) and return the decoded merged BAM. Panics
/// on a non-zero exit so a broken run never masquerades as a parity result.
fn run_leg(
    case: &Case,
    out: &Path,
    leg: Leg<'_>,
    threads: usize,
    chunk: Chunk,
    scheduler: Scheduler,
) -> Bam {
    run_leg_logged(case, out, leg, threads, chunk, scheduler).0
}

/// [`run_leg`], also returning the run's stderr (its log).
fn run_leg_logged(
    case: &Case,
    out: &Path,
    leg: Leg<'_>,
    threads: usize,
    chunk: Chunk,
    scheduler: Scheduler,
) -> (Bam, String) {
    let output = leg_command(case, out, leg, threads, chunk, scheduler)
        .output()
        .expect("run `fgumi runall`");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "`fgumi runall` ({leg:?}) failed with status {}:\n{stderr}",
        output.status
    );
    (read_bam(out), stderr)
}

/// Normalize for a cross-backend comparison (tag order ignored).
fn for_parity(bam: &Bam) -> NormBam {
    normalize(bam, TagOrder::Sort)
}

#[rstest]
#[case::release("0.12.0", "0.12.0")]
#[case::off_tag("0.12.0-ea58288", "0.12.0")]
#[case::off_tag_dirty("0.12.0-ea58288-dirty", "0.12.0")]
#[case::at_tag_dirty("0.12.0-dirty", "0.12.0")]
#[case::full_sha("0.7.0-e4b1f7c8e176cd27799a2732a1a4029dd29bc406", "0.7.0")]
#[case::prerelease_kept("1.0.0-rc1", "1.0.0-rc1")]
fn without_git_suffix_strips_only_the_dev_suffix(#[case] version: &str, #[case] expected: &str) {
    assert_eq!(without_git_suffix(version), expected);
}

/// The `--meth` `@PG` version (`<version>-meth`) loses the git dev suffix from
/// before its `-meth`, so a CLI built from a git checkout compares equal to
/// the vendored build.
#[rstest]
#[case::release("0.14.0", "0.14.0")]
#[case::dev("0.14.0-5c1d5e3", "0.14.0")]
#[case::meth_release("0.14.0-meth", "0.14.0-meth")]
#[case::meth_dev("0.14.0-5c1d5e3-meth", "0.14.0-meth")]
#[case::meth_dev_dirty("0.14.0-5c1d5e3-dirty-meth", "0.14.0-meth")]
fn normalized_pg_version_strips_the_dev_suffix_before_meth(
    #[case] version: &str,
    #[case] expected: &str,
) {
    assert_eq!(normalized_pg_version(version), expected);
}

// ---------------------------------------------------------------------------
// Parity matrix
// ---------------------------------------------------------------------------

/// Byte-parity of the in-process backend against the subprocess `bwa-mem3`
/// preset over a curated matrix.
///
/// A full cartesian (`threads {1,2,4,16}` × `sub-batch {1,7,64,256,1024}` ×
/// `chunk {one,several,many}` × `scheduler {auto,drain-first,chain-order}` ×
/// `input {PE,SE,mixed}` = 540 cells, each two `runall` runs) is impractical for
/// CI wall time, so this is a representative subset. It keeps every
/// parity-critical path: threads 1 vs 16; the smallest sub-batches (1, 7), which
/// give many sub-batches per cohort and exercise both id formulas (sub-batching
/// never splits a pair); a sub-batch larger than the cohort (1024);
/// one/several/many cohorts; a forced mid-pair `-K` cut (`Chunk::MidPair`, on
/// mixed input); all three schedulers (including the default `auto`, which runs
/// `RefillDrainScheduler`); RC reads; and `Layout::Mixed`. It still covers every
/// value of each individual dimension across the cases.
#[rstest]
// Base cell, then one-factor-at-a-time variation off it (PE unless noted).
#[case::pe_base(Input::Paired, 1, 256, Chunk::One, Scheduler::DrainFirst)]
#[case::pe_threads16(Input::Paired, 16, 256, Chunk::One, Scheduler::DrainFirst)]
#[case::pe_threads2_sb64(Input::Paired, 2, 64, Chunk::One, Scheduler::DrainFirst)]
#[case::pe_threads4_sb64_chainorder(Input::Paired, 4, 64, Chunk::One, Scheduler::ChainOrder)]
#[case::pe_sb1(Input::Paired, 1, 1, Chunk::One, Scheduler::DrainFirst)]
#[case::pe_sb7_odd(Input::Paired, 1, 7, Chunk::One, Scheduler::DrainFirst)]
#[case::pe_sb1024_over(Input::Paired, 1, 1024, Chunk::One, Scheduler::DrainFirst)]
#[case::pe_chunk_several(Input::Paired, 1, 256, Chunk::Several, Scheduler::DrainFirst)]
#[case::pe_chunk_many(Input::Paired, 1, 256, Chunk::Many, Scheduler::DrainFirst)]
#[case::pe_chainorder(Input::Paired, 1, 256, Chunk::One, Scheduler::ChainOrder)]
#[case::pe_threads16_auto(Input::Paired, 16, 256, Chunk::Several, Scheduler::Auto)]
// Single-end: SE id formula + RC in the SE path.
#[case::se_base(Input::Single, 1, 256, Chunk::One, Scheduler::DrainFirst)]
#[case::se_threads16_sb1_many(Input::Single, 16, 1, Chunk::Many, Scheduler::DrainFirst)]
#[case::se_chainorder(Input::Single, 1, 256, Chunk::One, Scheduler::ChainOrder)]
// Mixed SE/PE: Layout::Mixed and both id formulas together.
#[case::mixed_base(Input::Mixed, 1, 256, Chunk::One, Scheduler::DrainFirst)]
#[case::mixed_stress(Input::Mixed, 1, 1, Chunk::Many, Scheduler::DrainFirst)]
#[case::mixed_threads16_sb7_chainorder(Input::Mixed, 16, 7, Chunk::Many, Scheduler::ChainOrder)]
#[case::mixed_threads16_sb7_auto(Input::Mixed, 16, 7, Chunk::Many, Scheduler::Auto)]
#[case::mixed_threads16_sb1024_several(
    Input::Mixed,
    16,
    1024,
    Chunk::Several,
    Scheduler::DrainFirst
)]
// Forced mid-pair -K cut: the split pair's halves must match and keep their tags.
#[case::mixed_midpair(Input::Mixed, 1, 256, Chunk::MidPair, Scheduler::DrainFirst)]
#[case::mixed_midpair_threads16_auto(Input::Mixed, 16, 7, Chunk::MidPair, Scheduler::Auto)]
fn inproc_matches_subprocess(
    #[case] input: Input,
    #[case] threads: usize,
    #[case] sub_batch_templates: usize,
    #[case] chunk: Chunk,
    #[case] scheduler: Scheduler,
) {
    let label = format!(
        "inproc_matches_subprocess[{input:?} t{threads} sb{sub_batch_templates} {chunk:?} {}]",
        scheduler.label()
    );
    let reference_bin = bin_or_skip!(label);
    let case = Case::new(&reference_bin, input);

    let subprocess = run_leg(
        &case,
        &case.out("subprocess.bam"),
        Leg::Subprocess { reference_bin: &reference_bin },
        threads,
        chunk,
        scheduler,
    );
    let inproc = run_leg(
        &case,
        &case.out("inproc.bam"),
        Leg::InProcess { sub_batch_templates, dedup: None },
        threads,
        chunk,
        scheduler,
    );

    assert_fixture_coverage(&subprocess, &label);
    if matches!(chunk, Chunk::MidPair) {
        let cut = mid_pair_chunk(&case.templates).expect("Chunk::MidPair has a cut point");
        assert!(
            split_names(&subprocess).contains(cut.template.as_bytes()),
            "{label}: the forced mid-pair -K cut did not split its target pair '{}' — the \
             case no longer exercises the mid-pair split",
            cut.template
        );
    }
    // Checked on both legs: parity alone cannot catch a defect both share.
    assert_tags_transferred(&subprocess, &case.templates, &format!("{label} subprocess"));
    assert_tags_transferred(&inproc, &case.templates, &format!("{label} in-process"));
    assert_reads_preserved(&subprocess, &case.templates, &format!("{label} subprocess"));
    assert_reads_preserved(&inproc, &case.templates, &format!("{label} in-process"));
    assert_same_bam(
        &for_parity(&subprocess),
        &for_parity(&inproc),
        &format!("{label}: in-process diverged from the subprocess bwa-mem3 preset"),
    );
}

// ---------------------------------------------------------------------------
// Duplicate read pairs (--aligner::dedup-reads)
// ---------------------------------------------------------------------------

impl Case {
    /// [`Case::new`] over paired input where about four in nine templates
    /// recur (every third, plus every ninth): each copy has its own name, lands both near its original and in later
    /// cohorts, and some recur twice, the shape of a PCR-duplicate-rich UMI
    /// library.
    fn with_duplicates(bin: &str) -> Self {
        let dir = TempDir::new().expect("create temp dir");
        let fixture = write_fixture_reference(dir.path(), bin);
        let base = simulate_templates(&fixture.seq, Input::Paired, N_TEMPLATES);
        let mut templates = Vec::with_capacity(base.len() * 3 / 2);
        // `t`'s reads, RX and flags under `name`.
        let named = |t: &Template, name: String| Template {
            name,
            reads: t.reads.clone(),
            paired: t.paired,
            rx: t.rx.clone(),
            qc_fail: t.qc_fail,
        };
        for (i, t) in base.iter().enumerate() {
            templates.push(named(t, t.name.clone()));
            if i % 3 == 1 {
                let dup = &base[i - 1];
                templates.push(named(dup, format!("{}_dup1", dup.name)));
            }
            if i % 9 == 8 {
                let dup = &base[i / 2];
                templates.push(named(dup, format!("{}_dup2", dup.name)));
            }
        }
        let unmapped = dir.path().join("unmapped.bam");
        write_unmapped_bam(&unmapped, &templates);
        Self { reference: fixture.fasta, unmapped, templates, dir, methylation: None }
    }
}

/// On a duplicate-rich input, the in-process backend with
/// `--aligner::dedup-reads` on or off is byte-identical to the subprocess
/// preset, and with it on the memo actually copied duplicates (its end-of-run
/// log line reports them).
#[rstest]
#[case::base(1, 256, Chunk::One, Scheduler::DrainFirst)]
#[case::threads16_sb7_many(16, 7, Chunk::Many, Scheduler::Auto)]
#[case::threads4_sb64_several(4, 64, Chunk::Several, Scheduler::ChainOrder)]
fn inproc_dedup_reads_matches_subprocess_on_duplicates(
    #[case] threads: usize,
    #[case] sub_batch_templates: usize,
    #[case] chunk: Chunk,
    #[case] scheduler: Scheduler,
    #[values("on", "off")] dedup: &str,
) {
    let label = format!(
        "inproc_dedup_reads_matches_subprocess_on_duplicates[t{threads} sb{sub_batch_templates} \
         {chunk:?} {} dedup-{dedup}]",
        scheduler.label()
    );
    let reference_bin = bin_or_skip!(label);
    let case = Case::with_duplicates(&reference_bin);
    let subprocess = run_leg(
        &case,
        &case.out("subprocess.bam"),
        Leg::Subprocess { reference_bin: &reference_bin },
        threads,
        chunk,
        scheduler,
    );
    let (inproc, log) = run_leg_logged(
        &case,
        &case.out("inproc.bam"),
        Leg::InProcess { sub_batch_templates, dedup: Some(dedup) },
        threads,
        chunk,
        scheduler,
    );
    let memo_line = log.lines().find(|l| l.contains("read-pair memo"));
    match dedup {
        "on" => {
            let line = memo_line.unwrap_or_else(|| panic!("{label}: no memo line in:\n{log}"));
            let dup_pairs: u64 = line
                .split_once("memo: ")
                .and_then(|(_, rest)| rest.split_once(" duplicate pairs"))
                .and_then(|(n, _)| n.parse().ok())
                .unwrap_or_else(|| panic!("{label}: unparseable memo line: {line}"));
            assert!(dup_pairs > 0, "{label}: the memo found no duplicates: {line}");
        }
        _ => assert!(memo_line.is_none(), "{label}: memo ran with dedup-reads off"),
    }
    assert_reads_preserved(&inproc, &case.templates, &format!("{label} in-process"));
    assert_tags_transferred(&inproc, &case.templates, &format!("{label} in-process"));
    assert_same_bam(
        &for_parity(&subprocess),
        &for_parity(&inproc),
        &format!("{label}: in-process diverged from the subprocess bwa-mem3 preset"),
    );
}

// ---------------------------------------------------------------------------
// Determinism
// ---------------------------------------------------------------------------

/// Two in-process runs of the same input at the same settings produce identical
/// output (tags in emitted order), under both the default scheduler (`auto`,
/// i.e. `RefillDrainScheduler`) and `drain-first`. Both legs are in-process;
/// `BWA_MEM3_BIN` is needed only to index the fixture reference.
#[rstest]
#[case::paired(Input::Paired)]
#[case::single(Input::Single)]
#[case::mixed(Input::Mixed)]
fn inproc_is_deterministic(
    #[case] input: Input,
    #[values(Scheduler::Auto, Scheduler::DrainFirst)] scheduler: Scheduler,
) {
    let label = format!("inproc_is_deterministic[{input:?} {}]", scheduler.label());
    let bin = bin_or_skip!(label);
    let case = Case::new(&bin, input);

    let leg = Leg::InProcess { sub_batch_templates: 64, dedup: None };
    let first = run_leg(&case, &case.out("first.bam"), leg, 4, Chunk::Several, scheduler);
    let second = run_leg(&case, &case.out("second.bam"), leg, 4, Chunk::Several, scheduler);
    assert_fixture_coverage(&first, &label);
    assert_same_bam(
        &normalize(&first, TagOrder::Keep),
        &normalize(&second, TagOrder::Keep),
        &format!("{label}: two in-process runs diverged"),
    );
}

/// The in-process backend is thread-invariant: `--threads 16` produces the same
/// bytes as `--threads 1`, under both the default scheduler (`auto`) and
/// `drain-first`. `--threads 1` runs the fused single-thread runtime; at
/// `--threads 16` the pooled runtime's parallel merge restores input order
/// through its ordinal-ordered output.
#[rstest]
#[case::paired(Input::Paired)]
#[case::single(Input::Single)]
#[case::mixed(Input::Mixed)]
fn inproc_is_thread_invariant(
    #[case] input: Input,
    #[values(Scheduler::Auto, Scheduler::DrainFirst)] scheduler: Scheduler,
) {
    let label = format!("inproc_is_thread_invariant[{input:?} {}]", scheduler.label());
    let bin = bin_or_skip!(label);
    let case = Case::new(&bin, input);

    let leg = Leg::InProcess { sub_batch_templates: 7, dedup: None };
    let one = run_leg(&case, &case.out("t1.bam"), leg, 1, Chunk::Many, scheduler);
    let sixteen = run_leg(&case, &case.out("t16.bam"), leg, 16, Chunk::Many, scheduler);
    assert_fixture_coverage(&one, &label);
    assert_same_bam(
        &normalize(&one, TagOrder::Keep),
        &normalize(&sixteen, TagOrder::Keep),
        &format!("{label}: --threads 16 diverged from --threads 1"),
    );
}

// ---------------------------------------------------------------------------
// Bisulfite-aware alignment (--methylation-mode)
// ---------------------------------------------------------------------------

impl Case {
    /// [`Case::new`] over directional EM-seq reads, with the
    /// `bwa-mem3 index --meth` dual index built beside the fixture's plain one.
    ///
    /// Each read is converted in its own orientation, as a directional library
    /// is: R1 (and a single-end read) C→T, R2 G→A, at every cytosine outside a
    /// `CpG` (EM-seq leaves methylated `CpG` cytosines unconverted).
    fn with_emseq(bin: &str) -> Self {
        let dir = TempDir::new().expect("create temp dir");
        let fixture = write_fixture_reference(dir.path(), bin);
        let status = Command::new(bin)
            .args(["index", "--meth"])
            .arg(&fixture.fasta)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .expect("run `bwa-mem3 index --meth`");
        assert!(status.success(), "`bwa-mem3 index --meth` failed with status {status}");

        let mut templates = simulate_templates(&fixture.seq, Input::Mixed, N_TEMPLATES);
        for t in &mut templates {
            for (i, read) in t.reads.iter_mut().enumerate() {
                // R2 reads the complementary strand, so its conversions are G→A.
                let r2 = t.paired && i == 1;
                *read = emseq_convert(read, r2);
            }
        }
        let unmapped = dir.path().join("unmapped.bam");
        write_unmapped_bam(&unmapped, &templates);
        Self { reference: fixture.fasta, unmapped, templates, dir, methylation: Some("em-seq") }
    }
}

/// Convert one read as sequenced. For an R1 (`r2 == false`) a C converts to T;
/// for an R2 the read is the reverse complement of the converted strand, so a
/// G converts to A. A `CpG` stays unconverted, tested in the read's own
/// orientation: C followed by G for an R1, G preceded by C for an R2.
fn emseq_convert(read: &[u8], r2: bool) -> Vec<u8> {
    let (from, to) = if r2 { (b'G', b'A') } else { (b'C', b'T') };
    (0..read.len())
        .map(|i| {
            let cpg =
                if r2 { i > 0 && read[i - 1] == b'C' } else { read.get(i + 1) == Some(&b'G') };
            if read[i] == from && !cpg { to } else { read[i] }
        })
        .collect()
}

/// How many primary records carry bwa-mem3's `--meth` strand tag `XG:Z`.
fn meth_tagged(bam: &Bam) -> usize {
    let xg = fgumi_raw_bam::SamTag::new(b'X', b'G');
    bam.records
        .iter()
        .filter(|b| align_common::is_primary(b) && align_common::string_tag(b, xg).is_some())
        .count()
}

/// Under `runall --methylation-mode em-seq`, the in-process backend aligns
/// bisulfite-aware exactly as the subprocess preset's `bwa-mem3 mem --meth`
/// does, record for record. (TAPS aligns plain, which the parity tests above
/// already cover; on an align-only chain runall rejects it as dead.)
#[rstest]
#[case::base(1, 256, Chunk::One, Scheduler::DrainFirst)]
#[case::threads4_sb7_many(4, 7, Chunk::Many, Scheduler::Auto)]
#[case::threads4_sb64_several(4, 64, Chunk::Several, Scheduler::ChainOrder)]
fn inproc_matches_subprocess_under_emseq(
    #[case] threads: usize,
    #[case] sub_batch_templates: usize,
    #[case] chunk: Chunk,
    #[case] scheduler: Scheduler,
) {
    let label = format!(
        "inproc_matches_subprocess_under_emseq[t{threads} sb{sub_batch_templates} {chunk:?} {}]",
        scheduler.label()
    );
    let reference_bin = bin_or_skip!(label);
    let case = Case::with_emseq(&reference_bin);

    let subprocess = run_leg(
        &case,
        &case.out("subprocess.bam"),
        Leg::Subprocess { reference_bin: &reference_bin },
        threads,
        chunk,
        scheduler,
    );
    let inproc = run_leg(
        &case,
        &case.out("inproc.bam"),
        Leg::InProcess { sub_batch_templates, dedup: None },
        threads,
        chunk,
        scheduler,
    );

    // Both legs really aligned with --meth: it alone writes XG:Z.
    for (leg, bam) in [("subprocess", &subprocess), ("in-process", &inproc)] {
        assert!(meth_tagged(bam) > 0, "{label}: the {leg} leg wrote no XG:Z, so it ignored --meth");
    }
    assert_tags_transferred(&subprocess, &case.templates, &format!("{label} subprocess"));
    assert_tags_transferred(&inproc, &case.templates, &format!("{label} in-process"));
    // --meth reports the original read, not the converted one it seeded with.
    assert_reads_preserved(&subprocess, &case.templates, &format!("{label} subprocess"));
    assert_reads_preserved(&inproc, &case.templates, &format!("{label} in-process"));
    assert_same_bam(
        &for_parity(&subprocess),
        &for_parity(&inproc),
        &format!("{label}: in-process diverged from the subprocess bwa-mem3 --meth preset"),
    );
}
