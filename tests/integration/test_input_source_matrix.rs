//! Every command declares — and is held to — its I/O contract.
//!
//! Three properties are easy to get wrong per-command and impossible to eyeball
//! across the whole CLI:
//!
//! 1. **Uncompressed SAM** must be accepted anywhere BAM is, because input
//!    format is a property of the data, not of the command.
//! 2. **stdin** (`-i -`) must be accepted by anything that streams, so commands
//!    compose in a pipeline.
//! 3. **stdout** (`-o -`) must be written by anything that streams, for the same
//!    reason — a pipeline needs both ends.
//!
//! All three were previously true of some commands and not others, and nothing
//! failed when a new command picked the wrong reader or writer. The table below
//! is the contract: every command the binary advertises must appear in it with an
//! explicit stance, and [`every_command_declares_an_io_contract`] fails when
//! one does not — so adding a command forces a decision rather than inheriting
//! whichever default its reader or writer happened to have.
//!
//! The stdout axis was added after `sort`, `fastq`, and `extract --threads N`
//! were found to create a regular file *named* `-` instead of streaming: they
//! reached for `File::create` directly rather than the stdout-aware writer the
//! other commands share. Every one of them exited zero while doing it, which is
//! why the check compares what arrived on the pipe rather than trusting a status.
//!
//! Checking any of these means running a command two ways and comparing what it
//! wrote, which is only evidence if what it writes depends on what it read. A
//! command handed input it cannot form a molecule from exits zero and writes a
//! header-only file, so the comparisons pass however badly the reader behaved.
//! That is why the table carries a fourth axis, `output_depends_on_input`, and
//! why [`every_command_output_depends_on_its_input`] drains each fixture and
//! demands the output change: it is the assertion that keeps the other three
//! honest.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use tempfile::TempDir;

use crate::helpers::bam_generator::{
    create_consensus_family, create_duplex_grouped_family, create_grouped_family,
    create_minimal_header, create_read_name_umi_records, create_single_strand_family,
    create_test_reference, create_umi_family, transcode_bam_to_sam, write_bam,
};

/// Whether a command is held to one of the properties [`CONTRACTS`] declares.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Support {
    /// The command must satisfy it.
    Required,
    /// The command legitimately cannot; the string is why.
    NotApplicable(&'static str),
}

use Support::{NotApplicable, Required};

/// One command's declared I/O contract.
struct IoContract {
    /// Subcommand name as it appears in `fgumi --help`.
    command: &'static str,
    /// Whether `-i in.sam` must work.
    sam: Support,
    /// Whether `-i -` must work.
    stdin: Support,
    /// Whether `-o -` must stream to stdout rather than write a file named `-`.
    ///
    /// `NotApplicable` here is for commands whose `-o` is not a single stream:
    /// the metrics commands and `review` take a *prefix* and write several files
    /// from it, which no pipe can carry.
    stdout: Support,
    /// Whether emptying the input must change the output.
    ///
    /// This is what keeps the other three axes honest. A command handed input it
    /// cannot form a molecule from exits zero and writes a header-only file, so
    /// its SAM, stdin and stdout runs match the oracle no matter how badly the
    /// reader behaved. Declaring `Required` here puts that failure mode under
    /// test — see [`every_command_output_depends_on_its_input`].
    output_depends_on_input: Support,
}

/// The contract for every command the binary advertises.
///
/// A `NotApplicable` reason is a claim about the command's design, not a
/// to-do — if the reason stops being true, the entry should change.
const CONTRACTS: &[IoContract] = &[
    IoContract {
        command: "extract",
        sam: NotApplicable("reads FASTQ, not alignment records"),
        // Only with a single input: one stdin cannot supply two FASTQs.
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "correct",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "fastq",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "zipper",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: NotApplicable(
            "emits one output template per uBAM template, and the uBAM arrives on `-u`; an \
             empty `-i` yields the same templates, unmapped",
        ),
    },
    IoContract {
        command: "sort",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "merge",
        sam: Required,
        stdin: NotApplicable("merges N pre-sorted inputs named in a list file"),
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "group",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "dedup",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    #[cfg(feature = "simplex")]
    IoContract {
        command: "simplex",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    #[cfg(feature = "duplex")]
    IoContract {
        command: "duplex",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    #[cfg(feature = "codec")]
    IoContract {
        command: "codec",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "filter",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "clip",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "retag",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    IoContract {
        command: "copy-umi",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    #[cfg(feature = "duplex")]
    IoContract {
        command: "duplex-metrics",
        sam: Required,
        stdin: Required,
        stdout: NotApplicable("`-o` names a prefix for several metrics files, not one stream"),
        output_depends_on_input: Required,
    },
    #[cfg(feature = "simplex")]
    IoContract {
        command: "simplex-metrics",
        sam: Required,
        stdin: Required,
        stdout: NotApplicable("`-o` names a prefix for several metrics files, not one stream"),
        output_depends_on_input: Required,
    },
    IoContract {
        command: "review",
        sam: NotApplicable("requires BAI indexes and does random access, which needs BGZF"),
        stdin: NotApplicable("does random access against a BAI index"),
        stdout: NotApplicable("`-o` names a prefix for a directory of BAMs and their indexes"),
        output_depends_on_input: NotApplicable("declares no axis this harness invokes"),
    },
    IoContract {
        command: "downsample",
        sam: Required,
        stdin: Required,
        stdout: Required,
        output_depends_on_input: Required,
    },
    #[cfg(feature = "compare")]
    IoContract {
        command: "compare",
        sam: Required,
        stdin: NotApplicable("compares two named files against each other"),
        stdout: NotApplicable("writes its report to stdout already and takes no `--output`"),
        output_depends_on_input: Required,
    },
    #[cfg(feature = "simulate")]
    IoContract {
        command: "simulate",
        sam: NotApplicable("generates data rather than reading it"),
        stdin: NotApplicable("generates data rather than reading it"),
        stdout: NotApplicable("takes several named outputs, not one stream"),
        output_depends_on_input: NotApplicable("generates data rather than reading it"),
    },
];

/// Subcommands that `fgumi --help` lists but that are not fgumi commands.
const NOT_A_COMMAND: &[&str] = &["help"];

/// Parses the subcommand names out of `fgumi --help`.
///
/// Reading the binary's own help keeps this test honest against the real CLI
/// rather than a hand-maintained copy of it: a command added to the clap enum
/// shows up here immediately.
fn advertised_commands() -> Vec<String> {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .arg("--help")
        .output()
        .expect("failed to run fgumi --help");
    let help = String::from_utf8_lossy(&output.stdout);

    let mut commands = Vec::new();
    let mut in_commands = false;
    for line in help.lines() {
        if line.starts_with("Commands:") {
            in_commands = true;
            continue;
        }
        if in_commands {
            // The section ends at the first blank line before `Options:`.
            if line.trim().is_empty() {
                break;
            }
            // Entries are indented; continuation lines are indented further and
            // carry no new name in the first column.
            let Some(name) = line.split_whitespace().next() else { continue };
            if line.starts_with("  ") && !name.starts_with('-') {
                commands.push(name.to_string());
            }
        }
    }

    assert!(!commands.is_empty(), "parsed no commands out of `fgumi --help`:\n{help}");
    commands
}

/// The contract table must cover exactly the commands the binary advertises.
///
/// This is the test that makes the rest of the file exhaustive: a new command
/// fails here until someone states its position on every axis [`IoContract`]
/// carries — SAM, stdin, stdout, and whether its output depends on its input.
#[test]
fn every_command_declares_an_io_contract() {
    let advertised: Vec<String> = advertised_commands()
        .into_iter()
        .filter(|name| !NOT_A_COMMAND.contains(&name.as_str()))
        .collect();

    let declared: Vec<&str> = CONTRACTS.iter().map(|c| c.command).collect();

    let undeclared: Vec<&String> =
        advertised.iter().filter(|name| !declared.contains(&name.as_str())).collect();
    assert!(
        undeclared.is_empty(),
        "these commands are advertised by `fgumi --help` but declare no I/O contract in \
         CONTRACTS: {undeclared:?}. Add an entry stating each command's stance on all four \
         axes: sam, stdin, stdout, and output_depends_on_input."
    );

    let stale: Vec<&&str> =
        declared.iter().filter(|name| !advertised.contains(&(**name).to_string())).collect();
    assert!(
        stale.is_empty(),
        "these commands declare an I/O contract but are no longer advertised by \
         `fgumi --help`: {stale:?}. Remove their CONTRACTS entries."
    );
}

/// Test data shaped for whichever command is under test.
///
/// Commands reject inputs that don't carry what they consume (MI tags, a
/// query-grouped header), so "does it accept SAM/stdin" can only be answered
/// with an input the command would otherwise be happy with.
struct Fixture {
    bam: PathBuf,
    sam: PathBuf,
    /// One-line list files naming `bam` / `sam`, for commands that take `-b`.
    bam_list: PathBuf,
    sam_list: PathBuf,
    /// A small uncompressed FASTQ, for `extract`, which reads FASTQ not BAM.
    fastq: PathBuf,
}

impl Fixture {
    /// The path to hand a command as its input, in the requested format.
    ///
    /// `merge` names its inputs in a list file rather than via `-i`, so its
    /// format coverage comes from the list pointing at a SAM path.
    fn input_for(&self, command: &str, sam: bool) -> &Path {
        match (command, sam) {
            // `extract` reads FASTQ; it has no BAM/SAM axis at all.
            ("extract", _) => &self.fastq,
            ("merge", true) => &self.sam_list,
            ("merge", false) => &self.bam_list,
            (_, true) => &self.sam,
            (_, false) => &self.bam,
        }
    }
}

/// The record shape a command accepts.
///
/// A command rejects input that doesn't carry what it consumes, so "does it
/// accept SAM/stdin" can only be answered with an input it would otherwise be
/// happy with — otherwise a fixture mismatch masquerades as a format failure.
#[derive(Debug, Clone, Copy)]
enum Shape {
    /// Two UMI families tagged `RX`, un-grouped.
    Ungrouped,
    /// One already-grouped family tagged `RX` + `MI`, trivially sorted.
    Grouped,
    /// Consensus-called reads tagged `cD`/`cE`.
    Consensus,
    /// Both strands of one molecule, tagged `RX` + strand-suffixed `MI`.
    DuplexGrouped,
    /// Two single-strand (`/A`-only) molecules of differing depth.
    SimplexGrouped,
    /// Overlapping read pairs tagged `MI` + `MC`, the CODEC shape.
    Codec,
    /// Reads with the UMI in the last `:`-field of the name and no `RX` tag.
    ReadNameUmi,
}

/// The record shape `command` requires.
///
/// "Requires" here means more than "is accepted by": a command handed input it
/// cannot form a molecule from exits zero and writes a header-only file, which
/// passes every comparison in this file vacuously. So the shape has to be one
/// the command produces *records* from — see
/// [`every_command_output_depends_on_its_input`], which enforces exactly that.
fn shape_for(command: &str) -> Shape {
    match command {
        // Rejects any read lacking cD/cE consensus-calling tags.
        "filter" => Shape::Consensus,
        // `downsample` requires MI; `merge` refuses an input that is not already
        // in template-coordinate order, which a single one-position family is.
        "downsample" | "merge" => Shape::Grouped,
        // Given ungrouped reads these emit nothing at all, so they need a real
        // molecule: `duplex` and `duplex-metrics` need both strands present.
        "duplex" | "duplex-metrics" => Shape::DuplexGrouped,
        // `simplex-metrics` rejects a base UMI seen on both strands, so it and
        // the simplex caller get single-strand molecules instead.
        "simplex" | "simplex-metrics" => Shape::SimplexGrouped,
        // CODEC gets duplex from a single overlapping pair rather than from two
        // strand-tagged templates, so it needs its own shape.
        "codec" => Shape::Codec,
        // `copy-umi` reads the UMI from the read name and rejects names without
        // one, so it needs names carrying a UMI field rather than tagged families.
        "copy-umi" => Shape::ReadNameUmi,
        _ => Shape::Ungrouped,
    }
}

/// Writes the input `command` accepts, as both BAM and the equivalent SAM.
fn write_fixture(dir: &Path, command: &str) -> Fixture {
    let bam = dir.join("input.bam");
    let sam = dir.join("input.sam");

    // `create_minimal_header` advertises SO:unsorted GO:query
    // SS:template-coordinate, which the grouping and consensus commands require.
    let header = create_minimal_header("chr1", 10_000);
    let records = match shape_for(command) {
        Shape::Consensus => create_consensus_family(2, "cons", "ACGTACGT"),
        Shape::Grouped => create_grouped_family("ACGTACGT", "1", 3, "fam_a", "ACGTACGTAC", 35),
        Shape::Ungrouped => {
            let mut records = create_umi_family("ACGTACGT", 3, "fam_a", "ACGTACGTAC", 35);
            records.extend(create_umi_family("TTTTGGGG", 2, "fam_b", "ACGTACGTAC", 35));
            records
        }
        Shape::DuplexGrouped => create_duplex_grouped_family("1", 2, "ACGTACGT", 30, 100),
        // Two families of different depth at different coordinates, so the
        // family-size table has more than one row to get wrong.
        Shape::SimplexGrouped => {
            let mut records = create_single_strand_family("1", 3, "ACGTACGT", 30, 100);
            records.extend(create_single_strand_family("2", 1, "ACGTACGT", 30, 400));
            records
        }
        // Reuses the CODEC command's own pair builder rather than a second copy
        // of the MI/MC/overlap invariants it encodes.
        Shape::Codec => (0..2)
            .flat_map(|i| {
                let (r1, r2) = crate::test_codec_command::create_codec_read_pair(
                    &format!("codec_{i}"),
                    b"ACGTACGT",
                    b"ACGTACGT",
                    &[30; 8],
                    &[30; 8],
                    100,
                    "UMI001",
                    None,
                );
                [r1, r2]
            })
            .collect(),
        Shape::ReadNameUmi => create_read_name_umi_records(),
    };
    write_bam(&bam, &header, &records);

    // `merge` refuses an input that is not already in template-coordinate order,
    // and that order tie-breaks the name lane by hash — so hand-ordering the
    // records is not reliably sorted. Sort with the tool rather than guessing.
    if command == "merge" {
        let sorted = dir.join("sorted.bam");
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args(["sort", "-i"])
            .arg(&bam)
            .arg("-o")
            .arg(&sorted)
            .args(["--order", "template-coordinate"])
            .status()
            .expect("failed to run fgumi sort");
        assert!(status.success(), "failed to template-coordinate sort the merge fixture");
        std::fs::rename(&sorted, &bam).expect("replace fixture with its sorted form");
    }

    transcode_bam_to_sam(&bam, &sam);

    // Four reads: enough for `extract` to emit records, small enough to stay cheap.
    let fastq = dir.join("input.fq");
    let fastq_text = (0..4).fold(String::new(), |mut text, i| {
        use std::fmt::Write as _;
        writeln!(text, "@read_{i}\nACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIII")
            .expect("writing to a String cannot fail");
        text
    });
    std::fs::write(&fastq, fastq_text).expect("write FASTQ");

    let bam_list = dir.join("bam.list");
    std::fs::write(&bam_list, format!("{}\n", bam.display())).expect("write BAM list");
    let sam_list = dir.join("sam.list");
    std::fs::write(&sam_list, format!("{}\n", sam.display())).expect("write SAM list");

    Fixture { bam, sam, bam_list, sam_list, fastq }
}

/// The arguments that exercise `command`'s record-reading path.
///
/// `{input}` and `{output}` are substituted. Returns `None` for commands whose
/// contract is `NotApplicable` on both axes, which are never invoked here.
fn invocation(command: &str) -> Option<Vec<&'static str>> {
    let args: Vec<&'static str> = match command {
        "correct" => vec!["correct", "-i", "{input}", "-o", "{output}", "-U", "{umis}", "-d", "1"],
        "fastq" => vec!["fastq", "-i", "{input}", "-o", "{output}"],
        "zipper" => {
            vec!["zipper", "-i", "{input}", "-u", "{bam}", "-r", "{reference}", "-o", "{output}"]
        }
        "sort" => vec!["sort", "-i", "{input}", "-o", "{output}", "--order", "queryname"],
        "merge" => {
            vec!["merge", "-b", "{input}", "-o", "{output}", "--order", "template-coordinate"]
        }
        "group" => vec!["group", "-i", "{input}", "-o", "{output}", "--strategy", "adjacency"],
        "dedup" => vec!["dedup", "-i", "{input}", "-o", "{output}"],
        "simplex" => vec!["simplex", "-i", "{input}", "-o", "{output}", "--min-reads", "1"],
        "duplex" => vec!["duplex", "-i", "{input}", "-o", "{output}", "--min-reads", "1"],
        "codec" => vec![
            "codec",
            "-i",
            "{input}",
            "-o",
            "{output}",
            "--min-reads",
            "1",
            "--min-duplex-length",
            "1",
        ],
        "filter" => {
            vec![
                "filter",
                "-i",
                "{input}",
                "-o",
                "{output}",
                "-r",
                "{reference}",
                "--min-reads",
                "1",
            ]
        }
        "clip" => {
            vec![
                "clip",
                "-i",
                "{input}",
                "-o",
                "{output}",
                "-r",
                "{reference}",
                "--clip-overlapping-reads",
            ]
        }
        "retag" => vec!["retag", "-i", "{input}", "-o", "{output}", "RX::copy::BX"],
        "copy-umi" => vec!["copy-umi", "-i", "{input}", "-o", "{output}"],
        "duplex-metrics" => vec!["duplex-metrics", "-i", "{input}", "-o", "{output}"],
        "simplex-metrics" => vec!["simplex-metrics", "-i", "{input}", "-o", "{output}"],
        "downsample" => {
            vec!["downsample", "-i", "{input}", "-o", "{output}", "-f", "1.0", "--seed", "42"]
        }
        "compare" => vec!["compare", "bams", "{input}", "{bam}"],
        "extract" => vec![
            "extract",
            "-i",
            "{input}",
            "-o",
            "{output}",
            "-r",
            "8M+T",
            "--sample",
            "S1",
            "--library",
            "L1",
        ],
        _ => return None,
    };
    Some(args)
}

/// Runs `command`'s invocation with `input` substituted in.
///
/// `run_tag` names a per-run subdirectory to write `{output}` into. Two runs of
/// the same command in one fixture directory would otherwise write to the same
/// path, so the second silently overwrites the first — which is fatal when one
/// run is meant to be the other's oracle. It also keeps the multi-file outputs
/// (the metrics commands write several files off one prefix) separated without
/// having to know each command's file naming.
fn run_command(
    command: &str,
    input: &str,
    dir: &Path,
    fixture: &Fixture,
    pipe: Option<&Path>,
    run_tag: &str,
) -> std::process::Output {
    run_command_writing(command, input, dir, fixture, pipe, run_tag, OutputSpec::NamedPath)
}

/// Where a run is told to put its output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputSpec {
    /// `{output}` is a path inside the run directory — the normal case.
    NamedPath,
    /// `{output}` is the carried stdout spelling, so the run is expected to
    /// stream to stdout.
    ///
    /// The spelling travels in the variant because `open_output_writer` honours
    /// two of them, `-` and `/dev/stdout`, and a writer that special-cased only
    /// `-` would pass a matrix that never asked for the other one.
    ///
    /// The child runs in a scratch working directory (see [`stray_file_dir`]), so
    /// a command that mistakenly treats `-` as a filename creates its file there
    /// rather than in the crate root: the repository stays clean and the stray
    /// file is somewhere the test can look for it by name. That directory is
    /// deliberately *not* the digested run directory, because a build with the
    /// heap profiler enabled also drops a `dhat-heap.json` in the working
    /// directory, which would otherwise read as an output difference.
    Stdout(&'static str),
}

/// The working directory a [`OutputSpec::Stdout`] run is given.
///
/// Anything the command leaves here it created by writing to a path it should
/// have streamed instead — see [`declared_stdout_support_holds`].
fn stray_file_dir(dir: &Path, run_tag: &str) -> PathBuf {
    dir.join(format!("{run_tag}-cwd"))
}

/// [`run_command`], with control over what `{output}` is substituted with.
fn run_command_writing(
    command: &str,
    input: &str,
    dir: &Path,
    fixture: &Fixture,
    pipe: Option<&Path>,
    run_tag: &str,
    output_spec: OutputSpec,
) -> std::process::Output {
    let args = invocation(command).expect("command has an invocation");

    let run_dir = dir.join(run_tag);
    std::fs::create_dir_all(&run_dir).expect("create per-run output dir");
    let output = match output_spec {
        OutputSpec::NamedPath => run_dir.join(format!("{command}.out")),
        OutputSpec::Stdout(spelling) => PathBuf::from(spelling),
    };
    let reference = create_test_reference(dir);
    let umis = dir.join("umis.txt");
    std::fs::write(&umis, "ACGTACGT\nTTTTGGGG\n").expect("write UMI list");

    let substituted: Vec<String> = args
        .iter()
        .map(|arg| match *arg {
            "{input}" => input.to_string(),
            "{output}" => output.display().to_string(),
            "{bam}" => fixture.bam.display().to_string(),
            "{reference}" => reference.display().to_string(),
            "{umis}" => umis.display().to_string(),
            other => other.to_string(),
        })
        .collect();

    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fgumi"));
    cmd.args(&substituted);
    if matches!(output_spec, OutputSpec::Stdout(_)) {
        let cwd = stray_file_dir(dir, run_tag);
        std::fs::create_dir_all(&cwd).expect("create the stdout run's working dir");
        cmd.current_dir(&cwd);
    }
    if let Some(path) = pipe {
        cmd.stdin(Stdio::from(std::fs::File::open(path).expect("open input to pipe")));
    }
    cmd.output().expect("failed to run fgumi")
}

/// Runs `check` over every command whose contract selects it, reporting *all*
/// failures rather than dying on the first.
///
/// Driving the command list from [`CONTRACTS`] instead of a hand-written
/// `#[values(...)]` list is what makes these matrices exhaustive: a command that
/// declares `Required` is exercised because it declared it, so the declaration
/// and the coverage cannot drift apart. (Feature-gated commands compile out of
/// `CONTRACTS` along with themselves, so they are simply absent here.)
fn for_each_command_requiring(
    axis: impl Fn(&IoContract) -> Support,
    check: impl Fn(&str, &Fixture, &Path) -> Result<(), String>,
) {
    let mut failures = Vec::new();
    let mut exercised = 0;

    for contract in CONTRACTS.iter().filter(|c| axis(c) == Required) {
        let dir = TempDir::new().expect("create temp dir");
        let fixture = write_fixture(dir.path(), contract.command);
        exercised += 1;
        if let Err(failure) = check(contract.command, &fixture, dir.path()) {
            failures.push(failure);
        }
    }

    assert!(exercised > 0, "no command declared this axis Required — the matrix ran nothing");
    assert!(
        failures.is_empty(),
        "{} command(s) failed:\n\n{}",
        failures.len(),
        failures.join("\n\n")
    );
}

/// A run's output, reduced to something comparable across input formats.
///
/// BAM outputs are rendered back to SAM text, so a dropped, reordered or
/// altered record fails while two things that legitimately differ are
/// normalized away:
///
/// - the `@PG` `CL` line, which names the input path;
/// - the **width** of integer aux values. SAM's `i` type carries no width, so a
///   `cD:i:10` read from SAM text is re-encoded as `Int32` while the same tag in
///   the original BAM may have been stored as `Int8`. Every tool picks its own
///   width (htslib narrows to the smallest that fits; noodles does not), so
///   requiring byte-identical aux encoding across a SAM round-trip would be
///   requiring something no BAM writer guarantees. Rendering to SAM compares the
///   values, which is the property that actually has to hold.
///
/// Everything else — FASTQ, metrics tables — is compared as text with the parts
/// that name the run scrubbed.
#[derive(Debug, PartialEq)]
enum OutputDigest {
    Bam(String),
    Text(String),
}

/// Renders a BAM output back to SAM text for comparison. See [`OutputDigest`].
fn render_bam_as_sam(path: &Path) -> String {
    use noodles::sam::alignment::io::Write as _;

    let (header, records) = crate::helpers::read_bam_output(path);
    let mut text = Vec::new();
    {
        let mut writer = noodles::sam::io::Writer::new(&mut text);
        writer.write_header(&header).expect("render SAM header");
        for record in &records {
            writer.write_alignment_record(&header, record).expect("render SAM record");
        }
    }
    String::from_utf8(text).expect("rendered SAM is UTF-8")
}

/// The digest key holding a run's stdout. See [`digest_run`].
///
/// Angle-bracketed so it reads as "not a file" in a failure message. No command's
/// `invocation` writes an output by this name, so it cannot collide with one; a key
/// present twice would surface as a difference rather than silently overwrite,
/// since the digest is a sorted `Vec` compared element-wise.
const STDOUT_KEY: &str = "<stdout>";

/// Reduces a run to a name-keyed digest: every file it wrote, plus its stdout.
///
/// Keyed by name rather than compared as a blob because the metrics commands write
/// several files from one `--output` prefix; a missing file must read as a
/// difference, not as a shorter list.
///
/// `stdout` is folded in under [`STDOUT_KEY`] because a command's output is not
/// necessarily a file: `compare bams` writes its report to stdout and takes no
/// `--output` at all, so a files-only digest was empty for both runs and the SAM
/// comparison silently degenerated to comparing exit statuses — the very weakness
/// this comparison exists to close. Empty stdout contributes no entry, so commands
/// that only write files are unaffected.
fn digest_run(run_dir: &Path, stdout: &[u8]) -> Vec<(String, OutputDigest)> {
    let mut digests = digest_files(run_dir);
    if !stdout.is_empty() {
        digests.push((
            STDOUT_KEY.to_string(),
            OutputDigest::Text(scrub_run_identity(&String::from_utf8_lossy(stdout))),
        ));
    }
    digests.sort_by(|a, b| a.0.cmp(&b.0));
    digests
}

/// Reduces every file `run_dir` holds to a name-keyed digest, in no particular
/// order — [`digest_run`] owns the sort, since it also contributes an entry. See
/// [`digest_run`].
fn digest_files(run_dir: &Path) -> Vec<(String, OutputDigest)> {
    let Ok(entries) = std::fs::read_dir(run_dir) else { return Vec::new() };

    let digests: Vec<(String, OutputDigest)> = entries
        .map(|entry| entry.expect("read output dir entry").path())
        .filter(|path| path.is_file())
        .map(|path| {
            // The file name is the same in both runs — only the parent
            // directory differs — so it is a stable key.
            let name =
                path.file_name().expect("output file has a name").to_string_lossy().to_string();
            let bytes = std::fs::read(&path).expect("read output file");
            let digest = if bytes.starts_with(&[0x1f, 0x8b]) {
                OutputDigest::Bam(render_bam_as_sam(&path))
            } else {
                OutputDigest::Text(scrub_run_identity(&String::from_utf8_lossy(&bytes)))
            };
            (name, digest)
        })
        .collect();
    digests
}

/// Removes the parts of a text output that name *which* run produced it.
///
/// The runs differ only in their input source and their output directory, and
/// commands that echo either into their output would otherwise never compare
/// equal. Anything else that differs is a real difference.
fn scrub_run_identity(text: &str) -> String {
    // Longest tag first: `STDIN_SAM_RUN` contains `SAM_RUN`, so scrubbing the
    // short one first would leave `stdin-<RUN>` behind and every stdin-SAM
    // comparison would fail on a difference that is not one.
    text.replace("input.sam", "<INPUT>")
        .replace("input.bam", "<INPUT>")
        .replace("input.fq", "<INPUT>")
        // The drained fixture is a different file, so its name must scrub to the
        // same token: a difference that is only the input's name would let
        // `every_command_output_depends_on_its_input` pass without the records
        // having mattered.
        .replace("empty.bam", "<INPUT>")
        .replace("empty.fq", "<INPUT>")
        .replace("sam.list", "<LIST>")
        .replace("bam.list", "<LIST>")
        .replace("empty.list", "<LIST>")
        .replace(STDIN_SAM_RUN, "<RUN>")
        .replace(STDIN_RUN, "<RUN>")
        .replace(STDOUT_DEV_RUN, "<RUN>")
        .replace(STDOUT_RUN, "<RUN>")
        .replace(SAM_RUN, "<RUN>")
        .replace(BAM_RUN, "<RUN>")
        .replace(EMPTY_RUN, "<RUN>")
}

/// Per-run output subdirectory names (see `run_command`'s `run_tag`).
const BAM_RUN: &str = "bam-run";
const SAM_RUN: &str = "sam-run";
const STDIN_RUN: &str = "stdin-run";
const STDIN_SAM_RUN: &str = "stdin-sam-run";
const STDOUT_RUN: &str = "stdout-run";
const STDOUT_DEV_RUN: &str = "stdout-dev-run";
const EMPTY_RUN: &str = "empty-run";

/// Every spelling of stdout `open_output_writer` honours, paired with the
/// per-run output subdirectory that spelling's run is given.
///
/// Kept in lockstep with `is_stdout_path`: a spelling the writer accepts but
/// this table omits is a spelling no command in the matrix is ever asked for.
const STDOUT_SPELLINGS: [(&str, &str); 2] = [("-", STDOUT_RUN), ("/dev/stdout", STDOUT_DEV_RUN)];

/// Describes the first way two runs' outputs differ, or `None` if they match.
///
/// A bare `assert_eq!` on the digests would dump two multi-megabyte record
/// dumps into the failure, which is unreadable. This names the file and the
/// offset that first diverges and quotes a window around it.
///
/// The oracle is always the run that read a named path (see
/// [`run_named_path_oracle`]); `candidate_label` names the input source being
/// held to it, so a failure says which axis broke rather than always blaming
/// SAM.
fn describe_difference(
    oracle: &[(String, OutputDigest)],
    candidate: &[(String, OutputDigest)],
    candidate_label: &str,
) -> Option<String> {
    let oracle_files: Vec<&String> = oracle.iter().map(|(name, _)| name).collect();
    let candidate_files: Vec<&String> = candidate.iter().map(|(name, _)| name).collect();
    if oracle_files != candidate_files {
        return Some(format!(
            "  files from a named path: {oracle_files:?}\n  \
             files from {candidate_label}: {candidate_files:?}"
        ));
    }

    for ((name, bam), (_, sam)) in oracle.iter().zip(candidate) {
        if bam == sam {
            continue;
        }
        let (kind, oracle_text, candidate_text) = match (bam, sam) {
            (OutputDigest::Bam(b), OutputDigest::Bam(s)) => ("BAM records", b, s),
            (OutputDigest::Text(b), OutputDigest::Text(s)) => ("text", b, s),
            _ => return Some(format!("  {name}: one run wrote BAM, the other wrote text")),
        };
        let at = oracle_text
            .char_indices()
            .zip(candidate_text.chars())
            .find(|((_, b), s)| b != s)
            .map_or_else(|| oracle_text.len().min(candidate_text.len()), |((i, _), _)| i);
        let window =
            |text: &str| -> String { text.chars().skip(at.saturating_sub(60)).take(180).collect() };
        return Some(format!(
            "  {name} ({kind}) diverges at byte {at}\n    \
             from a named path: …{}…\n    from {candidate_label}: …{}…",
            window(oracle_text),
            window(candidate_text)
        ));
    }
    None
}

/// Runs `command` over its fixture named on the command line — the input source
/// every other axis is held to.
///
/// That fixture is a BAM for every command but `extract`, which reads FASTQ, so
/// the wording here stays "a named path" rather than "BAM".
///
/// This is the control as much as the oracle: if the command cannot process the
/// fixture at all, whatever the SAM or stdin run did proves nothing, so that is
/// reported as the failure rather than a difference nobody can interpret.
fn run_named_path_oracle(
    command: &str,
    fixture: &Fixture,
    dir: &Path,
    candidate_label: &str,
) -> Result<std::process::Output, String> {
    let named_input = fixture.input_for(command, false).display().to_string();
    let oracle_run = run_command(command, &named_input, dir, fixture, None, BAM_RUN);
    if oracle_run.status.success() {
        Ok(oracle_run)
    } else {
        Err(format!(
            "{command}: failed on a named path, so the {candidate_label} comparison would \
             prove nothing: {}",
            String::from_utf8_lossy(&oracle_run.stderr)
        ))
    }
}

/// Holds a run to the named-path oracle's output, file by file.
///
/// Exit status alone is too weak for any of these contracts: a command that read
/// the input a different way and dropped, reordered or altered records still
/// exits zero, and on the stdin axis a lost prefix or a double-read of fd 0
/// leaves a header-only output that also exits zero. Comparing against the
/// named-path run is what turns "it ran" into "it produced the same thing".
///
/// # What the `oracle.is_empty()` guard does and does not prove
///
/// It proves the harness can *see* the run's artifacts, not that those artifacts
/// are sensitive to the records that went in — a command that writes a
/// header-only BAM from any input compares equal here no matter what the reader
/// did. That is the vacuity this comparison exists to avoid, so it is checked
/// separately rather than assumed: every command declares
/// `output_depends_on_input`, and [`every_command_output_depends_on_its_input`]
/// holds it to that by draining the fixture and requiring the output to change.
fn compare_to_oracle(
    command: &str,
    dir: &Path,
    oracle_run: &std::process::Output,
    candidate_run: &std::process::Output,
    candidate_tag: &str,
    candidate_label: &str,
) -> Result<(), String> {
    let oracle = digest_run(&dir.join(BAM_RUN), &oracle_run.stdout);
    let candidate = digest_run(&dir.join(candidate_tag), &candidate_run.stdout);

    // A command whose artifacts this harness cannot see would compare two empty
    // digests and pass, which is indistinguishable from a genuine match. Fail
    // loudly instead, so a new command that writes somewhere unexpected has to
    // be taught about rather than silently exempted.
    if oracle.is_empty() {
        return Err(format!(
            "{command}: the named-path run produced no comparable output, so the \
             {candidate_label} comparison would be vacuous — give this command an \
             `{{output}}` in `invocation`, or have it write something this digest can see"
        ));
    }

    if let Some(difference) = describe_difference(&oracle, &candidate, candidate_label) {
        return Err(format!(
            "{command}: accepted {candidate_label} but produced different output than the \
             same records from a named path\n{difference}"
        ));
    }
    Ok(())
}

/// Every command declaring `sam: Required` must accept a `.sam` input **and
/// produce the same output from it**.
///
/// SAM must be accepted *anywhere BAM is*, so the BAM run is the oracle and the
/// SAM run is compared against it file by file. See [`compare_to_oracle`].
#[test]
fn declared_sam_support_holds() {
    for_each_command_requiring(
        |c| c.sam,
        |command, fixture, dir| {
            let oracle_run = run_named_path_oracle(command, fixture, dir, "SAM")?;

            let sam_input = fixture.input_for(command, true).display().to_string();
            let sam_run = run_command(command, &sam_input, dir, fixture, None, SAM_RUN);
            if !sam_run.status.success() {
                return Err(format!(
                    "{command}: declares SAM support but rejected SAM input: {}",
                    String::from_utf8_lossy(&sam_run.stderr)
                ));
            }

            compare_to_oracle(command, dir, &oracle_run, &sam_run, SAM_RUN, "SAM")
        },
    );
}

/// Every command declaring `stdin: Required` must accept `-i -` **and produce
/// the same output from it as from the same records in a named file**.
///
/// stdin is the non-rewindable path this harness exists to cover: the format is
/// sniffed by tee-and-replay rather than by seeking, so a prefix consumed and
/// not replayed, or a second reader that re-reads an already-drained fd 0, is
/// the live failure mode — and every one of those still exits zero over a
/// header-only output. See [`compare_to_oracle`].
#[test]
fn declared_stdin_support_holds() {
    for_each_command_requiring(
        |c| c.stdin,
        |command, fixture, dir| {
            let oracle_run = run_named_path_oracle(command, fixture, dir, "stdin")?;

            let piped_source = fixture.input_for(command, false).to_path_buf();
            let piped = run_command(command, "-", dir, fixture, Some(&piped_source), STDIN_RUN);
            if !piped.status.success() {
                return Err(format!(
                    "{command}: declares stdin support but rejected `-i -`: {}",
                    String::from_utf8_lossy(&piped.stderr)
                ));
            }

            compare_to_oracle(command, dir, &oracle_run, &piped, STDIN_RUN, "stdin")
        },
    );
}

/// Every command declaring `stdout: Required` must write **every** stdout
/// spelling in [`STDOUT_SPELLINGS`] to stdout, and put the same bytes there as
/// it would have put in a named file.
///
/// The failure this exists to catch is silent: a writer that reaches for
/// `File::create` instead of the stdout-aware one creates a regular file *named*
/// `-` in the working directory, writes the output into it, and exits zero. The
/// pipe is simply empty, so a status check — or anything that only looks at the
/// declared output path — sees a clean run. Two things are therefore asserted:
/// that no file named after the spelling was created, and that what arrived on
/// the pipe matches the named-path oracle byte for byte. See
/// [`compare_to_oracle`].
///
/// Both spellings run for every command rather than only `-`, because they part
/// company inside `is_stdout_path`: `/dev/stdout` reaching `File::create` still
/// lands on the pipe on Linux and macOS, so the stray-file check alone cannot
/// catch it — only the byte-for-byte comparison can.
#[test]
fn declared_stdout_support_holds() {
    for_each_command_requiring(
        |c| c.stdout,
        |command, fixture, dir| {
            let oracle_run = run_named_path_oracle(command, fixture, dir, "a stdout spelling")?;
            let named_input = fixture.input_for(command, false).display().to_string();

            // Both spellings, every command: the two differ only inside
            // `is_stdout_path`, so a writer that special-cased `-` alone would
            // pass a matrix that only ever asked for `-`.
            for (spelling, run_tag) in STDOUT_SPELLINGS {
                let label = format!("`-o {spelling}`");
                let streamed = run_command_writing(
                    command,
                    &named_input,
                    dir,
                    fixture,
                    None,
                    run_tag,
                    OutputSpec::Stdout(spelling),
                );
                if !streamed.status.success() {
                    return Err(format!(
                        "{command}: declares stdout support but rejected {label}: {}",
                        String::from_utf8_lossy(&streamed.stderr)
                    ));
                }

                // Checked before the digests, because this is the actual bug and
                // deserves to say so rather than surfacing as an empty output.
                // Only a relative spelling can land in the run's working
                // directory; `/dev/stdout` names a path outside it either way.
                let run_dir = dir.join(run_tag);
                if Path::new(spelling).is_relative()
                    && stray_file_dir(dir, run_tag).join(spelling).is_file()
                {
                    return Err(format!(
                        "{command}: {label} created a regular file named `{spelling}` instead of \
                         writing to stdout — its writer needs to go through `open_output_writer`"
                    ));
                }

                // The oracle wrote `<command>.out`; this run wrote the same bytes
                // to the pipe. Materializing them under that name is what lets the
                // two runs digest under the same key, so the existing file-by-file
                // comparison applies unchanged. Skipped when the pipe is empty, so
                // that case reads as a missing file rather than an empty one.
                if !streamed.stdout.is_empty() {
                    std::fs::write(run_dir.join(format!("{command}.out")), &streamed.stdout)
                        .expect("materialize the stdout run's stream");
                }
                let materialized = std::process::Output {
                    status: streamed.status,
                    stdout: Vec::new(),
                    stderr: streamed.stderr.clone(),
                };

                compare_to_oracle(command, dir, &oracle_run, &materialized, run_tag, &label)?;
            }
            Ok(())
        },
    );
}

/// SAM on stdin is the compound case: the format is sniffed from a stream that
/// cannot be rewound, so it exercises the tee-and-replay header parse rather
/// than the seek-based one — and the replayed prefix is SAM text rather than a
/// BGZF block, so both normalizations run over a pipe at once.
///
/// Only commands that require *both* axes qualify — `extract` streams stdin but
/// reads FASTQ, so there is no SAM form of its input.
#[test]
fn declared_stdin_support_accepts_piped_sam() {
    for_each_command_requiring(
        |c| {
            if c.sam == Required && c.stdin == Required {
                Required
            } else {
                NotApplicable("needs both a SAM form and stdin support")
            }
        },
        |command, fixture, dir| {
            let oracle_run = run_named_path_oracle(command, fixture, dir, "piped SAM")?;

            let piped = run_command(command, "-", dir, fixture, Some(&fixture.sam), STDIN_SAM_RUN);
            if !piped.status.success() {
                return Err(format!(
                    "{command}: accepts piped BAM but rejected piped SAM: {}",
                    String::from_utf8_lossy(&piped.stderr)
                ));
            }

            compare_to_oracle(command, dir, &oracle_run, &piped, STDIN_SAM_RUN, "piped SAM")
        },
    );
}

/// An input with no records, in whatever form `command` reads.
///
/// Same header and same shape as the real fixture, minus the records — so the
/// only thing that changed between the two runs is the data itself.
///
/// The header is read back from the fixture rather than re-minted with
/// `create_minimal_header`, because for `merge` those two are *not* the same
/// header: `write_fixture` re-writes the merge fixture through `fgumi sort` to
/// get template-coordinate order, and that stamps an `@PG ID:fgumi` line a fresh
/// header does not have. `merge` propagates input `@PG` records into its output
/// and `read_bam_output` normalizes only each `@PG`'s `CL` value, not its
/// presence — so a re-minted header left the populated and drained outputs
/// differing by that one line, and `merge` satisfied this whole test on a
/// provenance difference with zero records on either side. That is precisely the
/// vacuity the test exists to reject.
fn drained_input(command: &str, dir: &Path) -> PathBuf {
    /// Write a record-free copy of the fixture's own header to `path`.
    fn empty_bam_matching_fixture(dir: &Path, path: &Path) {
        let (header, records) = crate::helpers::read_bam_output(&dir.join("input.bam"));
        assert!(
            !records.is_empty(),
            "the fixture must have records for draining it to mean anything"
        );
        write_bam(path, &header, &[]);
    }

    match command {
        // `extract` reads FASTQ, where "no records" is an empty file.
        "extract" => {
            let path = dir.join("empty.fq");
            std::fs::write(&path, "").expect("write empty FASTQ");
            path
        }
        // `merge` names its inputs in a list file, so the list must point at an
        // empty BAM rather than being empty itself — an empty list is a
        // different error, not a drained input.
        "merge" => {
            let bam = dir.join("empty.bam");
            empty_bam_matching_fixture(dir, &bam);
            let list = dir.join("empty.list");
            std::fs::write(&list, format!("{}\n", bam.display())).expect("write empty list");
            list
        }
        _ => {
            let path = dir.join("empty.bam");
            empty_bam_matching_fixture(dir, &path);
            path
        }
    }
}

/// Commands that answer a drained input with a non-zero exit, and the phrase
/// whose presence in stderr identifies that answer.
///
/// A non-zero exit is a legitimate — and stronger — form of "the output depends
/// on the input" than writing different bytes, but only when the command exits
/// for the declared reason. The reason is pinned, not just the command name,
/// because the two entries here fail in genuinely different ways and neither
/// would notice being replaced by an unrelated failure:
///
/// - `extract` *refuses* the input: it cannot infer a quality encoding with no
///   records to look at.
/// - `compare` does not refuse anything — it runs to completion and exits 1
///   because the drained file and the populated oracle **differ**, which is
///   itself the dependence being tested.
///
/// A command that starts failing for a new reason — a bad flag combination, a
/// rejected header, a panic — no longer matches its phrase and surfaces here
/// instead of being absorbed.
/// Each phrase must be one the command emits regardless of `RUST_LOG`. `extract`
/// fails with an `anyhow` error, which reaches stderr unconditionally. `compare`
/// is pinned to its **stdout** `RESULT:` line rather than the matching
/// `info!("BAM files differ")`: log records are filtered by `RUST_LOG`
/// (`main.rs` builds the logger with `Env::default().default_filter_or("info")`),
/// and `run_command` inherits the parent environment, so a developer running
/// under `RUST_LOG=warn` would lose the log line and see this test fail on a
/// correctly-behaving binary. The `RESULT:` line is a `println!` gated only on
/// `--quiet`, which this harness never passes.
const DRAINED_INPUT_NONZERO_EXITS: &[(&str, &str)] = &[
    ("extract", "Cannot detect quality encoding: no records provided"),
    #[cfg(feature = "compare")]
    ("compare", "RESULT: BAM files DIFFER"),
];

/// Hold a drained run's non-zero exit to what [`DRAINED_INPUT_NONZERO_EXITS`]
/// declares for `command`.
///
/// Both streams are searched: which one carries the phrase is an implementation
/// detail of how the command reports (an `anyhow` bail on stderr, a `println!`
/// report on stdout), and pinning the wrong stream is a test failure that says
/// nothing about the binary.
fn check_declared_nonzero_exit(command: &str, stdout: &[u8], stderr: &[u8]) -> Result<(), String> {
    let stdout = String::from_utf8_lossy(stdout);
    let stderr = String::from_utf8_lossy(stderr);

    let Some((_, expected)) =
        DRAINED_INPUT_NONZERO_EXITS.iter().find(|(declared, _)| *declared == command)
    else {
        return Err(format!(
            "{command}: exited non-zero on a drained input, so the output comparison never ran \
             and this command's matrices could still be passing vacuously. Either fix the \
             invocation so it succeeds, or — if refusing a drained input is the contract — add \
             it to `DRAINED_INPUT_NONZERO_EXITS` with the phrase that identifies the refusal. \
             stdout:\n{stdout}\nstderr:\n{stderr}"
        ));
    };

    if stdout.contains(expected) || stderr.contains(expected) {
        return Ok(());
    }
    Err(format!(
        "{command}: exited non-zero on a drained input, but not for the reason \
         `DRAINED_INPUT_NONZERO_EXITS` declares — expected stdout or stderr to contain \
         {expected:?}. A different failure means the declared exemption is now hiding something \
         else; update the declaration or fix the failure. stdout:\n{stdout}\nstderr:\n{stderr}"
    ))
}

/// Every command declaring `output_depends_on_input: Required` must produce
/// different output when its input has no records.
///
/// This is the test that stops the other two from passing vacuously. Both of
/// them compare a SAM or stdin run against a named-path run, which proves
/// nothing about a command whose output is the same either way: `simplex` on
/// ungrouped reads, or `duplex-metrics` with no consensus tags to count, emits
/// an identical header-only file whether the reader delivered every record or
/// none, so a reader that silently dropped the whole stream would still pass.
///
/// Rather than curating a list of which commands are safe, this drives the same
/// fixture through each command twice — once populated, once drained — and
/// fails if the two agree. A future fixture that stops producing records for a
/// command breaks here, at the cause, instead of quietly hollowing out the
/// matrices.
///
/// The one thing that *is* declared rather than derived is a non-zero exit on
/// the drained run: it can only be compared as bytes if the command produced
/// bytes, so [`DRAINED_INPUT_NONZERO_EXITS`] names the commands that legitimately
/// exit non-zero and the reason each must give.
#[test]
fn every_command_output_depends_on_its_input() {
    // Which declared exemptions actually fired, so a stale one cannot outlive the
    // contract it documents — the same guard `CONTRACTS` gets from its `stale`
    // check. If a command grows a fallback and starts exiting zero on a drained
    // input, it takes the byte-comparison branch below, its entry is never
    // consulted, and the declaration silently becomes a lie.
    let consulted = std::cell::RefCell::new(Vec::<&str>::new());

    for_each_command_requiring(
        |c| c.output_depends_on_input,
        |command, fixture, dir| {
            let populated = run_named_path_oracle(command, fixture, dir, "drained input")?;

            let drained_path = drained_input(command, dir);
            let drained = run_command(
                command,
                &drained_path.display().to_string(),
                dir,
                fixture,
                None,
                EMPTY_RUN,
            );

            // A non-zero exit can itself be the dependence on the input — but only
            // when it is the exit this command is declared to give, for the declared
            // reason. Accepting *any* failure would re-open the vacuity hole this
            // test exists to close: a bad flag combination, a rejected header or a
            // panic would all be absorbed, and the output comparison below would
            // silently never run.
            if !drained.status.success() {
                let verdict =
                    check_declared_nonzero_exit(command, &drained.stdout, &drained.stderr);
                if verdict.is_ok()
                    && let Some((declared, _)) = DRAINED_INPUT_NONZERO_EXITS
                        .iter()
                        .find(|(declared, _)| *declared == command)
                {
                    consulted.borrow_mut().push(declared);
                }
                return verdict;
            }

            let with_records = digest_run(&dir.join(BAM_RUN), &populated.stdout);
            let without_records = digest_run(&dir.join(EMPTY_RUN), &drained.stdout);
            if describe_difference(&with_records, &without_records, "drained input").is_none() {
                return Err(format!(
                    "{command}: produced identical output from its fixture and from an input \
                     with no records, so every comparison in this file passes vacuously for \
                     it — give it a fixture it emits records from (see `shape_for`), or \
                     declare `output_depends_on_input: NotApplicable` with the reason"
                ));
            }
            Ok(())
        },
    );

    let consulted = consulted.into_inner();
    let stale: Vec<&str> = DRAINED_INPUT_NONZERO_EXITS
        .iter()
        .map(|(command, _)| *command)
        .filter(|command| !consulted.contains(command))
        .collect();
    assert!(
        stale.is_empty(),
        "these commands are declared in `DRAINED_INPUT_NONZERO_EXITS` but exited zero on a \
         drained input, so their declared exemption was never used: {stale:?}. They now go \
         through the output comparison like every other command — remove their entries."
    );
}

/// One stdin cannot supply two FASTQs, so `extract` must reject `-` alongside a
/// second input rather than reading the same stream twice and silently emitting
/// records for only one of them.
#[test]
fn extract_rejects_stdin_when_more_than_one_input_is_given() {
    let dir = TempDir::new().expect("create temp dir");
    let fixture = write_fixture(dir.path(), "extract");
    let output = dir.path().join("extract.bam");

    let result = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["extract", "-i", "-", "-i"])
        .arg(&fixture.fastq)
        .args(["-r", "8M+T", "-r", "8M+T", "-o"])
        .arg(&output)
        .args(["--sample", "S1", "--library", "L1"])
        .stdin(Stdio::from(std::fs::File::open(&fixture.fastq).expect("open FASTQ to pipe")))
        .output()
        .expect("failed to run fgumi");

    assert!(!result.status.success(), "extract accepted stdin alongside a second input");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("exactly one --input"),
        "error should explain the single-input restriction, got: {stderr}"
    );
}
