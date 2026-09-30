#![deny(unsafe_code)]

//! Guard for the subprocess align+merge path.
//!
//! The subprocess align path is a `SubprocessAlignStep` (emits `ZipperBatch`)
//! plus a shared, `Parallel` `MergeAlignedStep`, behind an `AlignBackend` trait
//! (formerly a single `AlignAndMergeStep`). This test runs `runall --start-from
//! align --stop-after zipper --aligner::preset bwa-mem3` over the shared
//! align-parity fixture (see `align_common`: substitutions, indels, unmapped,
//! chimeric, discordant and multi-mapping reads, each input record carrying an
//! `RX` tag, some QC-fail) and asserts:
//!
//! - the output holds every alignment shape the fixture is built to produce, so
//!   the comparisons below never run on trivial output;
//! - an independent expectation against the input, not just self-consistency:
//!   every input read yields exactly one primary record carrying its template's
//!   `RX` and QC-fail state;
//! - the decompressed records (with the `@PG` `CL` field stripped) are identical
//!   across `--threads` values and across a re-run. The aligner output is
//!   `-t`-invariant when `-K` is fixed, and the merge runs on a `Parallel` step
//!   whose `ByItemOrdinal` output restores input order, so a lost ordinal,
//!   mis-paired half or dropped batch shows up as a divergence;
//! - with mixed single/paired input and a `-K` chosen to cut between a pair's
//!   reads, the targeted pair is the one that splits, both halves keep their `RX`
//!   and QC-fail state and their own reads, and the output is byte-identical
//!   across `--threads` at that `-K`.
//!
//! ## Environment gating (skip-as-pass, unless tools are required)
//!
//! The test needs a real `bwa-mem3`: `BWA_MEM3_BIN` if set, else `bwa-mem3` on
//! `PATH`. It indexes the fixture reference and is the aligner the preset runs.
//! When neither is available the test prints a skip line and returns (a pass);
//! setting `FGUMI_BWA_MEM3_REQUIRE_TOOLS` makes that a hard failure, so the CI
//! e2e-parity job cannot silently skip it.

mod align_common;

use std::path::PathBuf;
use std::process::Command;

use align_common::{
    Bam, Input, TagOrder, Template, assert_fixture_coverage, assert_reads_preserved,
    assert_same_bam, assert_tags_transferred, bwa_mem3_bin_from_env, mid_pair_chunk, normalize,
    read_bam, simulate_templates, split_names, tools_required, write_fixture_reference,
    write_unmapped_bam,
};
use rstest::rstest;
use tempfile::TempDir;

/// Templates simulated per input.
const N_TEMPLATES: usize = 64;

/// The `bwa-mem3` to use: `BWA_MEM3_BIN`, else `bwa-mem3` on `PATH`.
fn bwa_mem3_bin() -> Option<String> {
    bwa_mem3_bin_from_env()
        .or_else(|| which::which("bwa-mem3").ok().map(|p| p.to_string_lossy().into_owned()))
}

/// Resolve the `bwa-mem3` binary, or skip (panic under
/// `FGUMI_BWA_MEM3_REQUIRE_TOOLS`) with a line naming the case.
macro_rules! bin_or_skip {
    ($label:expr) => {
        match bwa_mem3_bin() {
            Some(bin) => bin,
            None => {
                let reason = "no bwa-mem3 (set BWA_MEM3_BIN or put bwa-mem3 on PATH)";
                assert!(!tools_required(), "{}: required tool missing: {reason}", $label);
                eprintln!("skip {}: {reason}", $label);
                return;
            }
        }
    };
}

/// A prepared reference + unmapped-BAM fixture for the subprocess runs.
struct TestEnv {
    dir: TempDir,
    bin: String,
    reference: PathBuf,
    unmapped: PathBuf,
    templates: Vec<Template>,
}

impl TestEnv {
    /// Build and index the fixture reference with `bin`, and write the unmapped
    /// input BAM of `input`-shaped templates drawn from it.
    fn setup(bin: String, input: Input) -> Self {
        let dir = TempDir::new().expect("create temp dir");
        let fixture = write_fixture_reference(dir.path(), &bin);
        let templates = simulate_templates(&fixture.seq, input, N_TEMPLATES);
        let unmapped = dir.path().join("unmapped.bam");
        write_unmapped_bam(&unmapped, &templates);
        Self { bin, reference: fixture.fasta, unmapped, templates, dir }
    }

    /// Runs `runall --start-from align --stop-after zipper --aligner::preset
    /// bwa-mem3` at `threads` (and `--aligner::chunk-size` when given) and returns
    /// the decoded output.
    fn run(&self, threads: usize, chunk_size: Option<u64>) -> Bam {
        let out = self.dir.path().join(format!("aligned.t{threads}.bam"));
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
        cmd.arg("runall")
            .args(["--start-from", "align"])
            .args(["--stop-after", "zipper"])
            .args(["--aligner::preset", "bwa-mem3"])
            .arg("--aligner-bin")
            .arg(&self.bin)
            .args(["--threads", &threads.to_string()])
            .arg("--input")
            .arg(&self.unmapped)
            .arg("--ref")
            .arg(&self.reference)
            .arg("--output")
            .arg(&out);
        if let Some(k) = chunk_size {
            cmd.args(["--aligner::chunk-size", &k.to_string()]);
        }
        let status = cmd.status().expect("run `fgumi runall`");
        assert!(status.success(), "`fgumi runall` (threads={threads}) failed with status {status}");
        read_bam(&out)
    }
}

/// The subprocess `bwa-mem3` align path must merge every read with its own tags
/// and produce byte-identical output for every `--threads`, and the same bytes
/// on a re-run.
#[rstest]
#[case::t1(1)]
#[case::t2(2)]
#[case::t8(8)]
fn subprocess_align_output_is_thread_invariant(#[case] threads: usize) {
    let label = format!("subprocess_align_output_is_thread_invariant[threads={threads}]");
    let bin = bin_or_skip!(label);
    let env = TestEnv::setup(bin, Input::Paired);
    let out = env.run(threads, None);

    assert_fixture_coverage(&out, &label);
    assert_tags_transferred(&out, &env.templates, &label);
    assert_reads_preserved(&out, &env.templates, &label);
    let out = normalize(&out, TagOrder::Keep);
    assert_same_bam(
        &normalize(&env.run(1, None), TagOrder::Keep),
        &out,
        &format!("{label}: diverged from threads=1"),
    );
    assert_same_bam(
        &normalize(&env.run(threads, None), TagOrder::Keep),
        &out,
        &format!("{label}: not stable across a re-run"),
    );
}

/// With mixed single/paired input and a `-K` that cuts between a pair's two
/// reads, bwa-mem3 aligns that pair as two unpaired reads. The pair the cut
/// targets must be the one that splits, both halves must keep the pair's `RX`
/// and QC-fail state and their own reads after the merge, and the output must be
/// byte-identical across `--threads` at the same `-K`, since the split path is
/// where the merge's reordering changes.
#[test]
fn subprocess_mid_pair_split_keeps_unmapped_tags() {
    let label = "subprocess_mid_pair_split_keeps_unmapped_tags";
    let bin = bin_or_skip!(label);
    let env = TestEnv::setup(bin, Input::Mixed);
    let cut = mid_pair_chunk(&env.templates).expect("mixed input has a mid-pair cut point");
    let out = env.run(1, Some(cut.chunk_size));

    assert!(
        split_names(&out).contains(cut.template.as_bytes()),
        "{label}: the forced mid-pair -K cut did not split its target pair '{}'",
        cut.template
    );
    assert_tags_transferred(&out, &env.templates, label);
    assert_reads_preserved(&out, &env.templates, label);
    assert_same_bam(
        &normalize(&out, TagOrder::Keep),
        &normalize(&env.run(8, Some(cut.chunk_size)), TagOrder::Keep),
        &format!("{label}: threads=8 diverged from threads=1"),
    );
}
