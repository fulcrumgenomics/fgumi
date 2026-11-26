//! The `star` aligner preset: STAR, the splice-aware RNA-seq aligner, run as
//! the align stage's subprocess.
//!
//! STAR has no interleaved-FASTQ input, so the preset feeds it unaligned SAM
//! (`--readFilesType SAM PE`) on stdin and reads SAM on stdout in input order
//! (`--outSAMorder PairedKeepInputOrder`, `--outSAMunmapped Within`), which is
//! what the fused zipper merge pairs reads by. Everything STAR writes besides
//! its alignments (`Log.*`, `SJ.out.tab`, the per-chunk files of
//! `PairedKeepInputOrder`) goes to a scratch directory fgumi creates and
//! removes.
//!
//! Before anything is spawned, the preset checks that the STAR index's
//! `chrName.txt` and `chrLength.txt` name the same sequences as the `--ref`
//! dictionary, so a mismatched index fails before the genome load.

use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, bail};

use super::{
    AlignerOptions, AlignerPreset, AlignerScratchDir, ResolvedAligner, ResolvedAlignerMode,
    ResolvedBackend, SubprocessRuntime,
};
use crate::commands::fastq::UnalignedReadFormat;
use crate::commands::sort::{TMP_DIRS_ENV, resolve_tmp_dirs};

/// STAR's default `--limitIObufferSize` first value: the bytes of input each
/// mapping thread reads per chunk.
const STAR_DEFAULT_INPUT_CHUNK_BYTES: u64 = 30_000_000;

/// fgumi heap bytes held per byte of SAM written to STAR. An in-flight
/// unmapped read is held as its raw BAM record (by capacity) inside a template;
/// measured at about 2.4 times the stdin bytes for 100 bp reads, and more for
/// short reads with long names, where per-record overhead dominates.
const STAR_HEAP_BYTES_PER_INPUT_BYTE: u64 = 3;

/// Margin over STAR's estimated read-ahead. STAR counts its own rewritten
/// per-chunk buffer, not the bytes fgumi wrote, and short reads cost more heap
/// per input byte than [`STAR_HEAP_BYTES_PER_INPUT_BYTE`].
const STAR_READ_AHEAD_MARGIN: u64 = 2;

/// The in-flight budget, in fgumi heap bytes, that covers STAR's read-ahead
/// before its first output.
///
/// With `--outSAMorder PairedKeepInputOrder` each of the `threads` mapping
/// threads reads one chunk of `input_chunk_bytes` (`--limitIObufferSize`'s
/// first value) before thread 0, the only thread that writes output, has mapped
/// its own chunk; so STAR may read `threads × input_chunk_bytes` of input
/// before it emits a byte. That is converted to heap with
/// [`STAR_HEAP_BYTES_PER_INPUT_BYTE`] and doubled ([`STAR_READ_AHEAD_MARGIN`]):
/// at the default 30 MB chunk, 180 MB of budget per thread (5.6 GiB at 32
/// threads). The budget bounds memory only by what STAR reads ahead, and the
/// in-flight gate still grows past it, up to `--aligner::max-read-ahead`, if
/// this estimate falls short.
#[must_use]
fn star_in_flight_budget(threads: usize, input_chunk_bytes: u64) -> u64 {
    (threads as u64)
        .saturating_mul(input_chunk_bytes)
        .saturating_mul(STAR_HEAP_BYTES_PER_INPUT_BYTE)
        .saturating_mul(STAR_READ_AHEAD_MARGIN)
}

/// STAR parameters the preset sets itself; passing one in
/// `--aligner::extra-args` would make STAR exit on a duplicate parameter, or
/// break the input/output contract the merge relies on.
const PRESET_OWNED_PARAMETERS: &[&str] = &[
    "--runMode",
    "--runThreadN",
    "--genomeDir",
    "--readFilesIn",
    "--readFilesType",
    "--readFilesCommand",
    "--readFilesManifest",
    "--readFilesPrefix",
    "--readNameSeparator",
    "--outStd",
    "--outSAMtype",
    "--outSAMorder",
    "--outSAMunmapped",
    "--outTmpDir",
    "--outFileNamePrefix",
];

/// The extra STAR parameters (`--aligner::extra-args`), split into words.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
struct StarExtraArgs {
    parameters: Vec<StarParameter>,
}

/// A normalized parameter and its values, including an inline first value.
#[derive(Debug, Clone, PartialEq, Eq)]
struct StarParameter {
    name: String,
    values: Vec<String>,
}

impl StarExtraArgs {
    /// Split `--aligner::extra-args` on whitespace and reject what the streaming
    /// preset cannot run with: a parameter the preset sets itself, two-pass
    /// mapping (`--twopassMode`, which re-reads the input), `--outFilterType
    /// BySJout` (incompatible with `PairedKeepInputOrder`), `--outSAMprimaryFlag
    /// AllBestScore` (several primaries per read), `--outSAMfilter` (drops reads
    /// from the output), `--readMapNumber` (stops early), and the parameters whose
    /// output goes to a side file in the scratch directory, which is removed
    /// (`--quantMode`, `--outReadsUnmapped Fastx`, chimeric detection via
    /// `--chimSegmentMin` > 0 or `--chimOutType`, `--soloType`, `--sjdbInsertSave
    /// All`). Also rejected: `--parametersFiles`, which could set any of these out
    /// of sight of this check; `--outSAMreadID` other than `Standard` and
    /// `--outSAMmode` other than `Full`, which change the read names or records
    /// the merge pairs; `--outSAMmultNmax` other than -1 or a positive limit;
    /// flag masks that change pairing, mapping, strand or alignment-classification
    /// bits; shared-memory genome loads (`--genomeLoad` other than `NoSharedMemory`);
    /// and quote characters, since words are split on
    /// whitespace (a value cannot contain a space). Both `--name value` and
    /// `--name=value` are normalized before validation, sizing and quoting.
    ///
    /// # Errors
    ///
    /// Returns an error naming the offending parameter, or if the first word is
    /// not a `--parameter`.
    fn parse(extra_args: &str) -> Result<Self> {
        let words: Vec<&str> = extra_args.split_whitespace().collect();
        if let Some(word) = words.iter().find(|w| w.contains(['\'', '"'])) {
            bail!(
                "--aligner::extra-args: {word:?} contains a quote character; the parameters are \
                 split on whitespace and passed to STAR verbatim, so quotes are not interpreted \
                 and a value cannot contain a space"
            );
        }
        if let Some(first) = words.first()
            && !is_parameter(first)
        {
            bail!(
                "--aligner::extra-args must start with a STAR parameter (e.g. \
                 `--outSAMmapqUnique 60`), got {first:?}"
            );
        }
        let mut parsed = Self::default();
        for word in words {
            if is_parameter(word) {
                let (name, inline) = word.split_once('=').map_or((word, None), |(name, value)| {
                    (name, (!value.is_empty()).then_some(value))
                });
                parsed.parameters.push(StarParameter {
                    name: name.to_string(),
                    values: inline.into_iter().map(str::to_string).collect(),
                });
            } else {
                parsed
                    .parameters
                    .last_mut()
                    .expect("first word is a parameter")
                    .values
                    .push(word.to_string());
            }
        }
        for (name, values) in parsed.parameters() {
            check_parameter(name, values)?;
        }
        Ok(parsed)
    }

    /// Each normalized parameter and its values.
    fn parameters(&self) -> impl Iterator<Item = (&str, &[String])> {
        self.parameters.iter().map(|p| (p.name.as_str(), p.values.as_slice()))
    }

    /// The bytes of input each STAR mapping thread reads per chunk: the first
    /// value of `--limitIObufferSize` when it is passed, else STAR's default.
    ///
    /// # Errors
    ///
    /// Returns an error if the value is missing or not a positive integer.
    fn input_chunk_bytes(&self) -> Result<u64> {
        let Some((_, values)) = self.parameters().find(|(name, _)| *name == "--limitIObufferSize")
        else {
            return Ok(STAR_DEFAULT_INPUT_CHUNK_BYTES);
        };
        let value = values.first().context(
            "--aligner::extra-args: --limitIObufferSize needs a value (input and output \
             buffer sizes in bytes)",
        )?;
        match value.parse::<u64>() {
            Ok(bytes) if bytes > 0 => Ok(bytes),
            _ => bail!(
                "--aligner::extra-args: --limitIObufferSize's first value must be a positive \
                 number of bytes, got {value:?}"
            ),
        }
    }

    /// The words, each single-quoted for `/bin/bash -c`.
    fn quoted(&self) -> impl Iterator<Item = String> + '_ {
        self.parameters.iter().flat_map(|parameter| {
            let mut name = parameter.name.clone();
            let mut values = parameter.values.as_slice();
            // STAR does not treat an inline value beginning with `--` as an option.
            if let Some(first) = values.first().filter(|v| v.starts_with("--")) {
                name.push('=');
                name.push_str(first);
                values = &values[1..];
            }
            std::iter::once(shell_quote(&name))
                .chain(values.iter().map(|v| shell_quote(v)))
                .collect::<Vec<_>>()
        })
    }
}

/// Whether `word` names a STAR parameter (`--name`), as opposed to a value
/// such as `-1`.
fn is_parameter(word: &str) -> bool {
    word.len() > 2 && word.starts_with("--")
}

/// Why a parameter whose output STAR writes beside `Log.out` is rejected: those
/// files land in the preset's scratch directory, which is removed after the run,
/// so the output would be lost without an error.
const SIDE_OUTPUT_REASON: &str = "it writes its output to a file in the preset's scratch \
     directory, which is removed after the run; run STAR through --aligner::command with \
     --aligner::input-format sam to keep such files";

/// SAM bits used to pair, orient and classify alignments. QC-fail and duplicate
/// bits may be changed without breaking the streaming merge.
const ALIGNMENT_FLAGS: u16 = {
    use fgumi_raw_bam::flags;
    flags::PAIRED
        | flags::PROPER_PAIR
        | flags::UNMAPPED
        | flags::MATE_UNMAPPED
        | flags::REVERSE
        | flags::MATE_REVERSE
        | flags::FIRST_SEGMENT
        | flags::LAST_SEGMENT
        | flags::SECONDARY
        | flags::SUPPLEMENTARY
};

/// Reject one extra STAR parameter the preset cannot run with (see
/// [`StarExtraArgs::parse`]).
fn check_parameter(name: &str, values: &[String]) -> Result<()> {
    let first = values.first().map(String::as_str);
    if PRESET_OWNED_PARAMETERS.contains(&name) {
        bail!(
            "--aligner::extra-args must not set {name}: --aligner::preset star sets it to \
             stream reads through STAR in input order"
        );
    }
    let reason = match (name, first) {
        ("--twopassMode", Some(v)) if v != "None" => Some(
            "two-pass mapping re-reads the input, which is a stream; build the genome with \
             annotated junctions, or pass a junction file with --sjdbFileChrStartEnd",
        ),
        ("--outFilterType", Some("BySJout")) => {
            Some("BySJout holds reads back and cannot keep them in input order")
        }
        ("--outSAMprimaryFlag", Some("AllBestScore")) => {
            Some("AllBestScore marks several primary alignments per read")
        }
        ("--outSAMfilter", Some(v)) if v != "None" => {
            Some("it drops reads from the output, and every read must come back")
        }
        ("--readMapNumber", Some(v)) if v != "-1" => {
            Some("it stops mapping before the end of the input")
        }
        ("--quantMode", Some(v)) if v != "-" => Some(SIDE_OUTPUT_REASON),
        ("--outReadsUnmapped", Some(v)) if v != "None" => Some(SIDE_OUTPUT_REASON),
        ("--chimSegmentMin", Some(v)) if v.parse::<i64>().map_or(true, |n| n > 0) => {
            Some(SIDE_OUTPUT_REASON)
        }
        ("--soloType", Some(v)) if v != "None" => Some(SIDE_OUTPUT_REASON),
        ("--chimOutType", _) | ("--sjdbInsertSave", Some("All")) => Some(SIDE_OUTPUT_REASON),
        ("--parametersFiles", _) => Some(
            "a parameters file can set parameters this check cannot see; pass them in \
             --aligner::extra-args instead",
        ),
        ("--outSAMreadID", Some(v)) if v != "Standard" => {
            Some("the merge pairs STAR's records with the input by read name")
        }
        ("--outSAMmode", Some(v)) if v != "Full" => {
            Some("the merge needs STAR's full record for every read")
        }
        ("--genomeLoad", value) if value != Some("NoSharedMemory") => Some(
            "shared-memory genome loads can leave a partial segment after termination; use NoSharedMemory",
        ),
        ("--outSAMmultNmax", value)
            if !value.is_some_and(|v| v.parse::<i64>().is_ok_and(|n| n == -1 || n > 0)) =>
        {
            Some("the merge needs at least one alignment per read; use -1 or a positive limit")
        }
        ("--outSAMflagAND", value)
            if !value.is_some_and(|v| {
                v.parse::<u16>().is_ok_and(|mask| mask & ALIGNMENT_FLAGS == ALIGNMENT_FLAGS)
            }) =>
        {
            Some(
                "the flag mask must preserve pairing, mapping, strand and alignment-classification bits",
            )
        }
        ("--outSAMflagOR", value)
            if !value.is_some_and(|v| {
                v.parse::<u16>().is_ok_and(|mask| mask & ALIGNMENT_FLAGS == 0)
            }) =>
        {
            Some(
                "the flag mask must not set pairing, mapping, strand or alignment-classification bits",
            )
        }
        _ => None,
    };
    if let Some(reason) = reason {
        let spelled = if values.is_empty() {
            name.to_string()
        } else {
            format!("{name} {}", values.join(" "))
        };
        bail!(
            "--aligner::extra-args: {spelled} is not supported by --aligner::preset star: {reason}"
        );
    }
    Ok(())
}

/// Single-quote `word` for `/bin/bash -c`, so a path or value is passed as one
/// argument whatever characters it holds.
fn shell_quote(word: &str) -> String {
    format!("'{}'", word.replace('\'', r"'\''"))
}

/// Create the scratch directory a STAR run writes into: `--outFileNamePrefix`
/// (its `Log.*` and `SJ.out.tab` files) and, beneath it, `--outTmpDir` (which
/// STAR must create itself). It is placed in the first `FGUMI_TMP_DIRS` entry,
/// else the system temporary directory (`TMPDIR`); runall's `--sort::tmp-dir`
/// does not apply to it.
///
/// # Errors
///
/// Returns an error if the directory cannot be created.
fn create_work_dir() -> Result<AlignerScratchDir> {
    let roots = resolve_tmp_dirs(&[], std::env::var(TMP_DIRS_ENV).ok().as_deref());
    let builder = {
        let mut b = tempfile::Builder::new();
        b.prefix("fgumi-star.");
        b
    };
    let dir = match roots.first() {
        Some(root) => {
            std::fs::create_dir_all(root)
                .with_context(|| format!("creating temp directory root {}", root.display()))?;
            builder.tempdir_in(root)
        }
        None => builder.tempdir(),
    }
    .context("creating the STAR scratch directory")?;
    Ok(Arc::new(dir))
}

/// What [`build_star_command`] needs.
struct StarCommand<'a> {
    /// The STAR binary: `STAR` on `PATH`, or `--aligner-bin`.
    binary: &'a str,
    /// `--aligner::index`, the STAR genome directory.
    index: &'a Path,
    /// `--runThreadN`.
    threads: usize,
    /// The scratch directory (see [`AlignerScratchDir`]).
    work_dir: &'a Path,
    /// `--aligner::extra-args`.
    extra_args: &'a StarExtraArgs,
}

/// The STAR scratch subdirectory passed as `--outTmpDir`; it must not exist
/// when STAR starts.
const STAR_TMP_SUBDIR: &str = "tmp";

/// Build the STAR shell command: unaligned SAM pairs on stdin, SAM in input
/// order on stdout, every unmapped read reported, and all side files in
/// `work_dir`. Paths and extra arguments are single-quoted.
///
/// # Errors
///
/// Returns an error if an index or scratch path is not UTF-8.
fn build_star_command(cmd: &StarCommand<'_>) -> Result<String> {
    let tmp = cmd.work_dir.join(STAR_TMP_SUBDIR);
    // `--outFileNamePrefix` is a prefix, not a directory: the trailing
    // separator puts `Log.out` etc. inside `work_dir`.
    let prefix = format!("{}/", utf8_path(cmd.work_dir)?);
    let mut words = vec![
        shell_quote(cmd.binary),
        "--runMode alignReads".to_string(),
        format!("--runThreadN {}", cmd.threads),
        format!("--genomeDir {}", shell_quote(utf8_path(cmd.index)?)),
        "--readFilesType SAM PE".to_string(),
        "--readFilesIn /dev/stdin".to_string(),
        "--readNameSeparator none".to_string(),
        "--outStd SAM".to_string(),
        "--outSAMtype SAM".to_string(),
        "--outSAMorder PairedKeepInputOrder".to_string(),
        "--outSAMunmapped Within".to_string(),
        format!("--outTmpDir {}", shell_quote(utf8_path(&tmp)?)),
        format!("--outFileNamePrefix {}", shell_quote(&prefix)),
    ];
    words.extend(cmd.extra_args.quoted());
    Ok(words.join(" "))
}

/// Shell commands are UTF-8 strings; never substitute replacement path bytes.
fn utf8_path(path: &Path) -> Result<&str> {
    path.to_str().with_context(|| format!("STAR paths must be valid UTF-8: {}", path.display()))
}

/// The file in a STAR genome directory listing its sequences and lengths.
const CHR_NAME_LENGTH_FILE: &str = "chrNameLength.txt";

/// Files every STAR genome directory holds.
const STAR_INDEX_FILES: &[&str] = &[
    "Genome",
    "SA",
    "SAindex",
    "genomeParameters.txt",
    "chrName.txt",
    "chrLength.txt",
    "chrStart.txt",
    CHR_NAME_LENGTH_FILE,
];

/// Check that `index` is a STAR genome directory and names the same sequences,
/// with the same lengths and in the same order, as the sequence dictionary
/// beside `reference`. STAR's genomic `@SQ` lines come from `chrName.txt` and
/// `chrLength.txt`; also check the combined file for consistency. A mismatch
/// would otherwise surface only after STAR loaded the genome.
///
/// # Errors
///
/// Returns an error if `index` is not a directory, lacks one of the index
/// files, the dictionary is missing or unreadable, or the two disagree.
fn check_index(index: &Path, reference: &Path) -> Result<()> {
    utf8_path(index)?;
    if !index.is_dir() {
        bail!(
            "--aligner::index is not a directory: {} (--aligner::preset star needs the STAR \
             genome directory built by `STAR --runMode genomeGenerate`)",
            index.display()
        );
    }
    for file in STAR_INDEX_FILES {
        let path = index.join(file);
        if !path.is_file() {
            bail!(
                "STAR index file not found: {} (is --aligner::index a STAR genome directory?)",
                path.display()
            );
        }
    }
    let dict = crate::reference::find_dict_path(reference).with_context(|| {
        format!(
            "no sequence dictionary found next to --ref {} (expected `<ref>.dict` or \
             `<ref-without-extension>.dict`); generate one with `samtools dict`",
            reference.display()
        )
    })?;
    let combined_header = read_chr_name_length(&index.join(CHR_NAME_LENGTH_FILE))?;
    let index_header = read_star_chromosomes(index)?;
    let dict_header = read_dict(&dict)?;
    [&index_header, &combined_header].into_iter().try_for_each(|header| {
        crate::pipeline::steps::align::validate_sq_consistency(&dict_header, header).with_context(
            || {
                format!(
                    "the STAR index {} does not match the sequence dictionary {}: build the index \
                 from the same FASTA as --ref",
                    index.display(),
                    dict.display()
                )
            },
        )
    })
}

/// Read STAR's `chrNameLength.txt` (one `name<TAB>length` line per sequence)
/// as a header holding those `@SQ` lines, in order.
fn read_chr_name_length(path: &Path) -> Result<noodles::sam::Header> {
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading STAR index file {}", path.display()))?;
    let mut builder = noodles::sam::Header::builder();
    for (i, line) in text.lines().enumerate().filter(|(_, line)| !line.is_empty()) {
        let parsed = line.split_once('\t').and_then(|(name, length)| {
            let length = length.trim().parse::<usize>().ok()?;
            Some((name, std::num::NonZeroUsize::new(length)?))
        });
        let (name, length) = parsed.with_context(|| {
            format!(
                "{} line {}: expected `name<TAB>length` with a positive length, got {line:?}",
                path.display(),
                i + 1
            )
        })?;
        builder = builder.add_reference_sequence(name, Map::<ReferenceSequence>::new(length));
    }
    Ok(builder.build())
}

/// Read the chromosome files STAR uses to build its genomic `@SQ` lines.
fn read_star_chromosomes(index: &Path) -> Result<noodles::sam::Header> {
    use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

    let names_path = index.join("chrName.txt");
    let lengths_path = index.join("chrLength.txt");
    let names = std::fs::read_to_string(&names_path)
        .with_context(|| format!("reading STAR index file {}", names_path.display()))?;
    let lengths = std::fs::read_to_string(&lengths_path)
        .with_context(|| format!("reading STAR index file {}", lengths_path.display()))?;
    // STAR reads names with getline, not formatted whitespace extraction.
    // Reject malformed lines instead of normalizing a name STAR would preserve
    // or skipping a blank line where STAR would stop reading chromosomes.
    let names: Vec<_> = names.split_terminator('\n').collect();
    if names.iter().any(|name| name.is_empty() || name.chars().any(char::is_whitespace)) {
        bail!(
            "{}: expected one nonempty chromosome name without whitespace per line",
            names_path.display()
        );
    }
    let lengths: Vec<_> = lengths.split_ascii_whitespace().collect();
    if names.len() != lengths.len() {
        bail!("STAR chrName.txt and chrLength.txt entry counts disagree in {}", index.display());
    }
    let mut header = noodles::sam::Header::default();
    for (name, length) in names.into_iter().zip(lengths) {
        let length = length.parse::<std::num::NonZeroUsize>().with_context(|| {
            format!(
                "{}: expected a positive length for {name:?}, got {length:?}",
                lengths_path.display()
            )
        })?;
        if header
            .reference_sequences_mut()
            .insert(name.into(), Map::<ReferenceSequence>::new(length))
            .is_some()
        {
            bail!("{}: duplicate sequence name {name:?}", names_path.display());
        }
    }
    Ok(header)
}

/// Read a sequence dictionary's header.
fn read_dict(path: &Path) -> Result<noodles::sam::Header> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("opening sequence dictionary {}", path.display()))?;
    noodles::sam::io::Reader::new(std::io::BufReader::new(file))
        .read_header()
        .with_context(|| format!("reading sequence dictionary {}", path.display()))
}

/// Resolve `--aligner::preset star` into a [`ResolvedAligner`]: validate the
/// reference, the `STAR` binary, the index against the `.dict`, and the extra
/// parameters; create the scratch directory; and build the command.
///
/// Called by [`AlignerOptions::resolve`] after its mode-independent checks.
///
/// # Errors
///
/// Returns an error if `--methylation-mode em-seq` is set, `--ref` or the
/// binary is missing, `--aligner::threads` is zero, `--aligner::index` is
/// missing or disagrees with the dictionary, an extra parameter is rejected,
/// or the scratch directory cannot be created.
pub(super) fn resolve_star(
    opts: &AlignerOptions,
    reference: &Path,
    top_threads: usize,
    aligner_bin: Option<&Path>,
) -> Result<ResolvedAligner> {
    let preset = AlignerPreset::Star;
    preset.validate_with_methylation(reference, aligner_bin, opts.methylation_mode)?;
    if opts.threads == Some(0) {
        bail!("--aligner::threads must be greater than 0 (got 0)");
    }
    let threads = opts.threads.unwrap_or(top_threads);
    let index = opts.index.as_deref().context(
        "--aligner::preset star requires --aligner::index <DIR>, the STAR genome directory \
         built by `STAR --runMode genomeGenerate` from the same FASTA as --ref",
    )?;
    let startup_timeout = opts.startup_timeout.unwrap_or(900);
    if startup_timeout == 0 {
        bail!("--aligner::startup-timeout must be at least 1 second");
    }
    check_index(index, reference)?;
    let extra_args = StarExtraArgs::parse(opts.extra_args.as_deref().unwrap_or(""))?;
    let input_chunk_bytes = extra_args.input_chunk_bytes()?;
    let work_dir = create_work_dir()?;
    let command = build_star_command(&StarCommand {
        binary: aligner_bin.map(utf8_path).transpose()?.unwrap_or(preset.binary_name()),
        index,
        threads,
        work_dir: work_dir.path(),
        extra_args: &extra_args,
    })?;
    Ok(ResolvedAligner {
        backend: ResolvedBackend::Subprocess {
            command,
            // STAR pairs mates itself; it never splits a pair as bwa's `-K`
            // chunking can.
            accept_mid_pair_split: false,
            input_format: UnalignedReadFormat::Sam,
            max_read_ahead: opts.max_read_ahead,
            runtime: SubprocessRuntime {
                work_dir: Some(work_dir),
                require_complete_pairs: true,
                min_in_flight_budget: star_in_flight_budget(threads, input_chunk_bytes),
                startup_timeout: Some(std::time::Duration::from_secs(startup_timeout)),
            },
        },
        chunk_size: opts.chunk_size,
        threads: Some(threads),
        mode: ResolvedAlignerMode::Preset(preset),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;

    #[rstest]
    #[case::default_chunk(32, STAR_DEFAULT_INPUT_CHUNK_BYTES, 32 * 30_000_000 * 3 * 2)]
    #[case::one_thread(1, STAR_DEFAULT_INPUT_CHUNK_BYTES, 180_000_000)]
    #[case::small_chunk(8, 4_000_000, 8 * 4_000_000 * 6)]
    #[case::saturates(usize::MAX, u64::MAX, u64::MAX)]
    fn star_budget_scales_with_threads_and_chunk(
        #[case] threads: usize,
        #[case] chunk: u64,
        #[case] expected: u64,
    ) {
        assert_eq!(star_in_flight_budget(threads, chunk), expected);
    }

    /// The 32-thread probe that stalled at the default 1.8 GiB budget read
    /// about 780 MB of input; the STAR budget covers that with room to spare.
    #[test]
    fn star_budget_covers_the_stalled_probe() {
        let budget = star_in_flight_budget(32, STAR_DEFAULT_INPUT_CHUNK_BYTES);
        assert!(budget > 2 * 1_800_000_000, "got {budget}");
    }

    #[rstest]
    #[case::empty("", STAR_DEFAULT_INPUT_CHUNK_BYTES)]
    #[case::unrelated("--outSAMmapqUnique 60", STAR_DEFAULT_INPUT_CHUNK_BYTES)]
    #[case::set("--limitIObufferSize 4000000 50000000", 4_000_000)]
    #[case::equals("--limitIObufferSize=100000000 50000000", 100_000_000)]
    #[case::empty_inline("--limitIObufferSize= 1000000 50000000", 1_000_000)]
    #[case::set_among_others(
        "--outSAMattributes NH HI AS nM --limitIObufferSize 1000 2000 --alignIntronMax 1000",
        1000
    )]
    fn input_chunk_bytes_reads_limit_io_buffer_size(#[case] args: &str, #[case] expected: u64) {
        let parsed = StarExtraArgs::parse(args).unwrap();
        assert_eq!(parsed.input_chunk_bytes().unwrap(), expected);
    }

    #[rstest]
    #[case::missing("--limitIObufferSize")]
    #[case::zero("--limitIObufferSize 0 50000000")]
    #[case::not_a_number("--limitIObufferSize big 50000000")]
    fn input_chunk_bytes_rejects_bad_values(#[case] args: &str) {
        let err = StarExtraArgs::parse(args).unwrap().input_chunk_bytes().unwrap_err();
        assert!(err.to_string().contains("--limitIObufferSize"), "got: {err}");
    }

    #[rstest]
    #[case::two_pass("--twopassMode Basic", "--twopassMode Basic is not supported")]
    #[case::by_sj_out("--outFilterType BySJout", "--outFilterType BySJout is not supported")]
    #[case::all_best("--outSAMprimaryFlag AllBestScore", "AllBestScore is not supported")]
    #[case::filter("--outSAMfilter KeepOnlyAddedReferences", "--outSAMfilter")]
    #[case::read_map_number("--readMapNumber 1000", "--readMapNumber 1000 is not supported")]
    #[case::owned_threads("--runThreadN 4", "must not set --runThreadN")]
    #[case::owned_order("--outSAMorder Paired", "must not set --outSAMorder")]
    #[case::owned_prefix("--outFileNamePrefix out/", "must not set --outFileNamePrefix")]
    #[case::owned_read_command("--readFilesCommand zcat", "must not set --readFilesCommand")]
    #[case::later_parameter(
        "--outSAMmapqUnique 60 --twopassMode Basic",
        "--twopassMode Basic is not supported"
    )]
    #[case::leading_value("60 --outSAMmapqUnique", "must start with a STAR parameter")]
    #[case::gene_counts("--quantMode GeneCounts", "--quantMode GeneCounts is not supported")]
    #[case::transcriptome_bam(
        "--quantMode TranscriptomeSAM GeneCounts",
        "--quantMode TranscriptomeSAM GeneCounts is not supported"
    )]
    #[case::unmapped_fastx("--outReadsUnmapped Fastx", "--outReadsUnmapped Fastx is not supported")]
    #[case::chimeric_segments("--chimSegmentMin 12", "--chimSegmentMin 12 is not supported")]
    #[case::chimeric_output("--chimOutType Junctions", "--chimOutType Junctions is not supported")]
    #[case::chimeric_output_bare("--chimOutType", "--chimOutType is not supported")]
    #[case::parameters_file("--parametersFiles p.txt", "--parametersFiles p.txt is not supported")]
    #[case::owned_manifest("--readFilesManifest m.tsv", "must not set --readFilesManifest")]
    #[case::owned_prefix_dir("--readFilesPrefix /data/", "must not set --readFilesPrefix")]
    #[case::read_id_number("--outSAMreadID Number", "--outSAMreadID Number is not supported")]
    #[case::sam_mode_none("--outSAMmode None", "--outSAMmode None is not supported")]
    #[case::sam_mode_equals("--outSAMmode=None", "--outSAMmode None is not supported")]
    #[case::owned_equals("--runThreadN=4", "must not set --runThreadN")]
    #[case::shared_keep("--genomeLoad LoadAndKeep", "--genomeLoad")]
    #[case::shared_remove("--genomeLoad LoadAndRemove", "--genomeLoad")]
    #[case::shared_exit("--genomeLoad=LoadAndExit", "--genomeLoad")]
    #[case::shared_unload("--genomeLoad Remove", "--genomeLoad")]
    #[case::shared_missing("--genomeLoad", "--genomeLoad")]
    #[case::solo("--soloType CB_UMI_Simple", "--soloType CB_UMI_Simple is not supported")]
    #[case::sjdb_insert_save("--sjdbInsertSave All", "--sjdbInsertSave All is not supported")]
    #[case::single_quote("--outSAMattrRGline 'ID:a SM:b'", "contains a quote character")]
    #[case::double_quote("--outSAMattrRGline \"ID:a\"", "contains a quote character")]
    #[case::read_name_separator("--readNameSeparator /", "must not set --readNameSeparator")]
    #[case::no_alignments("--outSAMmultNmax 0", "--outSAMmultNmax 0 is not supported")]
    #[case::negative_limit("--outSAMmultNmax -2", "--outSAMmultNmax -2 is not supported")]
    #[case::invalid_limit("--outSAMmultNmax all", "--outSAMmultNmax all is not supported")]
    #[case::cleared_flags("--outSAMflagAND 0", "--outSAMflagAND 0 is not supported")]
    fn extra_args_rejected(#[case] args: &str, #[case] needle: &str) {
        let err = StarExtraArgs::parse(args).unwrap_err().to_string();
        assert!(err.contains(needle), "got: {err}");
    }

    #[rstest]
    #[case::empty("")]
    #[case::mapq("--outSAMmapqUnique 60")]
    #[case::two_pass_none("--twopassMode None")]
    #[case::filter_type_normal("--outFilterType Normal")]
    #[case::read_map_all("--readMapNumber -1")]
    #[case::attributes("--outSAMattributes NH HI AS nM NM MD --alignIntronMax 500000")]
    #[case::quant_mode_off("--quantMode -")]
    #[case::unmapped_none("--outReadsUnmapped None")]
    #[case::chimeric_off("--chimSegmentMin 0")]
    #[case::read_id_standard("--outSAMreadID Standard")]
    #[case::sam_mode_full("--outSAMmode Full")]
    #[case::no_shared_memory("--genomeLoad NoSharedMemory")]
    #[case::no_shared_memory_equals("--genomeLoad=NoSharedMemory")]
    #[case::solo_off("--soloType None")]
    #[case::sjdb_insert_basic("--sjdbInsertSave Basic")]
    #[case::all_alignments("--outSAMmultNmax -1")]
    #[case::one_alignment("--outSAMmultNmax 1")]
    #[case::several_alignments("--outSAMmultNmax 10")]
    #[case::unchanged_flags("--outSAMflagAND 65535 --outSAMflagOR 0")]
    #[case::qc_and_duplicate_flags("--outSAMflagAND 63999 --outSAMflagOR 1536")]
    fn extra_args_accepted(#[case] args: &str) {
        StarExtraArgs::parse(args).unwrap();
    }

    #[rstest]
    fn extra_args_preserve_alignment_flags(
        #[values(1_u16, 2, 4, 8, 16, 32, 64, 128, 256, 2048)] bit: u16,
        #[values("--outSAMflagAND", "--outSAMflagOR")] parameter: &str,
    ) {
        let mask = if parameter == "--outSAMflagAND" { !bit } else { bit };
        let args = format!("{parameter} {mask}");
        assert!(StarExtraArgs::parse(&args).is_err(), "must reject {args}");
    }

    #[rstest]
    #[case::spaces("--outSAMmapqUnique 60 --outSAMattributes NH HI")]
    #[case::equals("--outSAMmapqUnique=60 --outSAMattributes=NH HI")]
    fn command_streams_sam_in_input_order_into_the_work_dir(#[case] args: &str) {
        let extra = StarExtraArgs::parse(args).unwrap();
        let command = build_star_command(&StarCommand {
            binary: "STAR",
            index: Path::new("/idx/star with spaces"),
            threads: 12,
            work_dir: Path::new("/scratch/fgumi star.abc"),
            extra_args: &extra,
        })
        .unwrap();
        assert_eq!(
            command,
            "'STAR' --runMode alignReads --runThreadN 12 --genomeDir '/idx/star with spaces' \
             --readFilesType SAM PE --readFilesIn /dev/stdin --readNameSeparator none \
             --outStd SAM --outSAMtype SAM \
             --outSAMorder PairedKeepInputOrder --outSAMunmapped Within \
             --outTmpDir '/scratch/fgumi star.abc/tmp' \
             --outFileNamePrefix '/scratch/fgumi star.abc/' \
             '--outSAMmapqUnique' '60' '--outSAMattributes' 'NH' 'HI'"
        );
    }

    #[rstest]
    #[case::embedded_equals("--outSAMattrRGline=ID:a=b", "'--outSAMattrRGline' 'ID:a=b'")]
    #[case::value_equals("--outSAMattrRGline ID:a=b", "'--outSAMattrRGline' 'ID:a=b'")]
    #[case::double_equals("--outSAMattrRGline==", "'--outSAMattrRGline' '='")]
    #[case::inline_option_like("--outSAMattrRGline=--literal", "'--outSAMattrRGline=--literal'")]
    fn extra_args_keep_inline_value_semantics(#[case] args: &str, #[case] expected: &str) {
        let parsed = StarExtraArgs::parse(args).unwrap();
        assert_eq!(parsed.quoted().collect::<Vec<_>>().join(" "), expected);
    }

    #[rstest]
    #[case::plain("abc", "'abc'")]
    #[case::space("my dir", "'my dir'")]
    #[case::quote("it's", r"'it'\''s'")]
    #[case::meta("$(rm -rf x);", "'$(rm -rf x);'")]
    fn shell_quote_makes_one_word(#[case] word: &str, #[case] expected: &str) {
        assert_eq!(shell_quote(word), expected);
    }

    /// Write a STAR-like index directory naming `sequences`.
    fn write_index(dir: &Path, sequences: &[(&str, u64)]) -> PathBuf {
        let index = dir.join("star");
        std::fs::create_dir_all(&index).unwrap();
        for file in STAR_INDEX_FILES {
            std::fs::write(index.join(file), b"").unwrap();
        }
        let text: String =
            sequences.iter().map(|(n, l)| [n, "\t", &l.to_string(), "\n"].concat()).collect();
        std::fs::write(index.join(CHR_NAME_LENGTH_FILE), text).unwrap();
        let mut names = String::new();
        let mut lengths = String::new();
        for (name, length) in sequences {
            use std::fmt::Write as _;
            writeln!(names, "{name}").unwrap();
            writeln!(lengths, "{length}").unwrap();
        }
        std::fs::write(index.join("chrName.txt"), names).unwrap();
        std::fs::write(index.join("chrLength.txt"), lengths).unwrap();
        index
    }

    /// Write `ref.fa` and its `ref.dict` naming `sequences`.
    fn write_reference(dir: &Path, sequences: &[(&str, u64)]) -> PathBuf {
        use std::fmt::Write as _;
        let reference = dir.join("ref.fa");
        std::fs::write(&reference, b">chr1\nACGT\n").unwrap();
        let mut dict = String::from("@HD\tVN:1.6\n");
        for (name, length) in sequences {
            writeln!(dict, "@SQ\tSN:{name}\tLN:{length}").unwrap();
        }
        std::fs::write(dir.join("ref.dict"), dict).unwrap();
        reference
    }

    const DICT: &[(&str, u64)] = &[("chr1", 1000), ("chr2", 500)];

    #[rstest]
    #[case::matching(DICT, None)]
    #[case::name(
        &[("chr1", 1000), ("chrX", 500)],
        Some("name mismatch at position 1: dict='chr2' aligner='chrX'")
    )]
    #[case::length(
        &[("chr1", 1001), ("chr2", 500)],
        Some("length mismatch for 'chr1': dict=1000 aligner=1001")
    )]
    #[case::order(
        &[("chr2", 500), ("chr1", 1000)],
        Some("name mismatch at position 0: dict='chr1' aligner='chr2'")
    )]
    #[case::fewer(&[("chr1", 1000)], Some("@SQ count (1) does not match reference dict @SQ count (2)"))]
    #[case::more(
        &[("chr1", 1000), ("chr2", 500), ("ERCC-00002", 1061)],
        Some("@SQ count (3) does not match reference dict @SQ count (2)")
    )]
    fn check_index_compares_chr_name_length_with_the_dict(
        #[case] index_sequences: &[(&str, u64)],
        #[case] needle: Option<&str>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let index = write_index(tmp.path(), index_sequences);
        let reference = write_reference(tmp.path(), DICT);
        match (check_index(&index, &reference), needle) {
            (Ok(()), None) => {}
            (Err(e), Some(needle)) => {
                let msg = format!("{e:#}");
                assert!(msg.contains(needle), "got: {msg}");
                assert!(msg.contains("does not match the sequence dictionary"), "got: {msg}");
            }
            (result, needle) => panic!("expected {needle:?}, got {result:?}"),
        }
    }

    #[rstest]
    #[case::sa("SA")]
    #[case::parameters("genomeParameters.txt")]
    #[case::names("chrName.txt")]
    #[case::lengths("chrLength.txt")]
    #[case::starts("chrStart.txt")]
    fn check_index_requires_the_index_files(#[case] file: &str) {
        let tmp = tempfile::tempdir().unwrap();
        let index = write_index(tmp.path(), DICT);
        std::fs::remove_file(index.join(file)).unwrap();
        let reference = write_reference(tmp.path(), DICT);
        let msg = check_index(&index, &reference).unwrap_err().to_string();
        assert!(msg.contains("STAR index file not found"), "got: {msg}");
    }

    #[rstest]
    #[case::name("chrName.txt", "chr1\nchrX\n")]
    #[case::length("chrLength.txt", "1001\n500\n")]
    #[case::too_few_lengths("chrLength.txt", "1000\n")]
    #[case::too_many_names("chrName.txt", "chr1\nchr2\nchr3\n")]
    #[case::zero_length("chrLength.txt", "0\n500\n")]
    #[case::invalid_length("chrLength.txt", "big\n500\n")]
    #[case::duplicate_name("chrName.txt", "chr1\nchr1\n")]
    #[case::blank_name("chrName.txt", "chr1\n\nchr2\n")]
    #[case::leading_space("chrName.txt", " chr1\nchr2\n")]
    #[case::trailing_space("chrName.txt", "chr1 \nchr2\n")]
    #[case::tab_separated_names("chrName.txt", "chr1\tchr2\n")]
    #[case::crlf_names("chrName.txt", "chr1\r\nchr2\r\n")]
    #[case::unicode_length_separator("chrLength.txt", "1000\u{2003}500\n")]
    fn check_index_checks_star_header_metadata(#[case] file: &str, #[case] text: &str) {
        let tmp = tempfile::tempdir().unwrap();
        let index = write_index(tmp.path(), DICT);
        std::fs::write(index.join("chrName.txt"), "chr1\nchr2\n").unwrap();
        std::fs::write(index.join("chrLength.txt"), "1000\n500\n").unwrap();
        std::fs::write(index.join(file), text).unwrap();
        let reference = write_reference(tmp.path(), DICT);
        assert!(check_index(&index, &reference).is_err(), "must reject {file}: {text:?}");
    }

    #[rstest]
    #[case::index(true)]
    #[case::scratch(false)]
    fn command_rejects_non_utf8_paths(#[case] bad_index: bool) {
        use std::os::unix::ffi::OsStrExt;
        let invalid = Path::new(std::ffi::OsStr::from_bytes(b"/tmp/non-utf8-\xff"));
        let extra = StarExtraArgs::default();
        let result = build_star_command(&StarCommand {
            binary: "STAR",
            index: if bad_index { invalid } else { Path::new("/idx/star") },
            threads: 1,
            work_dir: if bad_index { Path::new("/tmp/star") } else { invalid },
            extra_args: &extra,
        });
        let error = result.unwrap_err();
        assert!(error.to_string().contains("UTF-8"), "got {error:#}");
    }

    #[test]
    fn check_index_requires_a_directory() {
        let tmp = tempfile::tempdir().unwrap();
        let reference = write_reference(tmp.path(), DICT);
        let msg = check_index(&tmp.path().join("missing"), &reference).unwrap_err().to_string();
        assert!(msg.contains("--aligner::index is not a directory"), "got: {msg}");
    }

    #[test]
    fn check_index_requires_a_dict() {
        let tmp = tempfile::tempdir().unwrap();
        let index = write_index(tmp.path(), DICT);
        let reference = write_reference(tmp.path(), DICT);
        std::fs::remove_file(tmp.path().join("ref.dict")).unwrap();
        let msg = check_index(&index, &reference).unwrap_err().to_string();
        assert!(msg.contains("no sequence dictionary found"), "got: {msg}");
    }

    #[rstest]
    #[case::no_tab("chr1 1000\n")]
    #[case::zero_length("chr1\t0\n")]
    #[case::not_a_number("chr1\tbig\n")]
    fn malformed_chr_name_length_is_reported(#[case] text: &str) {
        let tmp = tempfile::tempdir().unwrap();
        let index = write_index(tmp.path(), DICT);
        std::fs::write(index.join(CHR_NAME_LENGTH_FILE), text).unwrap();
        let reference = write_reference(tmp.path(), DICT);
        let msg = format!("{:#}", check_index(&index, &reference).unwrap_err());
        assert!(msg.contains("expected `name<TAB>length`"), "got: {msg}");
    }

    /// An executable stand-in for the `STAR` binary, passed as `--aligner-bin`.
    fn fake_binary(dir: &Path) -> PathBuf {
        use std::os::unix::fs::PermissionsExt;
        // Its own directory: on a case-insensitive filesystem `STAR` would
        // collide with the `star` index directory.
        let bin_dir = dir.join("bin");
        std::fs::create_dir_all(&bin_dir).unwrap();
        let bin = bin_dir.join("STAR");
        std::fs::write(&bin, "#!/bin/sh\n").unwrap();
        std::fs::set_permissions(&bin, std::fs::Permissions::from_mode(0o755)).unwrap();
        bin
    }

    /// A star-preset `AlignerOptions` with a matching index under `dir`, plus
    /// the reference and a stand-in binary.
    fn star_setup(dir: &Path) -> (AlignerOptions, PathBuf, PathBuf) {
        let index = write_index(dir, DICT);
        let reference = write_reference(dir, DICT);
        let opts = AlignerOptions {
            preset: Some(AlignerPreset::Star),
            index: Some(index),
            ..AlignerOptions::default()
        };
        (opts, reference, fake_binary(dir))
    }

    #[rstest]
    #[case::default_chunk(None, Some(6), STAR_DEFAULT_INPUT_CHUNK_BYTES, 6)]
    #[case::top_threads(Some("--limitIObufferSize 1000000 50000000"), None, 1_000_000, 16)]
    #[case::inline_chunk(Some("--limitIObufferSize=100000000 50000000"), Some(2), 100_000_000, 2)]
    fn resolve_builds_a_sam_streaming_star_backend(
        #[case] extra_args: Option<&str>,
        #[case] threads: Option<usize>,
        #[case] chunk: u64,
        #[case] expected_threads: usize,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let (mut opts, reference, bin) = star_setup(tmp.path());
        opts.extra_args = extra_args.map(str::to_string);
        opts.threads = threads;
        let resolved = opts.resolve(&reference, 16, Some(&bin)).unwrap();
        assert_eq!(resolved.threads, Some(expected_threads));
        assert!(matches!(resolved.mode, ResolvedAlignerMode::Preset(AlignerPreset::Star)));
        let ResolvedBackend::Subprocess {
            command,
            accept_mid_pair_split,
            input_format,
            runtime,
            ..
        } = resolved.backend
        else {
            panic!("star resolves to the subprocess backend");
        };
        assert_eq!(input_format, UnalignedReadFormat::Sam);
        assert!(!accept_mid_pair_split);
        assert_eq!(runtime.startup_timeout, Some(std::time::Duration::from_mins(15)));
        assert_eq!(runtime.min_in_flight_budget, star_in_flight_budget(expected_threads, chunk));
        let work_dir = runtime.work_dir.expect("star has a scratch directory");
        assert!(work_dir.path().is_dir());
        assert!(command.starts_with(&format!("'{}' --runMode alignReads", bin.display())));
        assert!(command.contains(&format!("--runThreadN {expected_threads} ")), "{command}");
        assert!(command.contains(&format!("--outFileNamePrefix '{}/'", work_dir.path().display())));
    }

    #[rstest]
    #[case::no_index(|o: &mut AlignerOptions| o.index = None, "requires --aligner::index")]
    #[case::zero_threads(
        |o: &mut AlignerOptions| o.threads = Some(0),
        "--aligner::threads must be greater than 0"
    )]
    #[case::fastq_input(
        |o: &mut AlignerOptions| o.input_format = Some(UnalignedReadFormat::Fastq),
        "--aligner::input-format fastq is not valid with --aligner::preset star"
    )]
    #[case::em_seq(
        |o: &mut AlignerOptions| o.methylation_mode = fgumi_consensus::MethylationMode::EmSeq,
        "and `star` has none"
    )]
    #[case::two_pass(
        |o: &mut AlignerOptions| o.extra_args = Some("--twopassMode Basic".into()),
        "--twopassMode Basic is not supported"
    )]
    #[case::zero_startup_timeout(
        |o: &mut AlignerOptions| o.startup_timeout = Some(0),
        "--aligner::startup-timeout must be at least 1 second"
    )]
    #[case::command_too(
        |o: &mut AlignerOptions| o.command = Some("STAR {ref}".into()),
        "mutually exclusive"
    )]
    fn resolve_rejects(#[case] change: fn(&mut AlignerOptions), #[case] needle: &str) {
        let tmp = tempfile::tempdir().unwrap();
        let (mut opts, reference, bin) = star_setup(tmp.path());
        change(&mut opts);
        let msg = format!("{:#}", opts.resolve(&reference, 4, Some(&bin)).unwrap_err());
        assert!(msg.contains(needle), "got: {msg}");
    }

    /// `--aligner::chunk-size` never reaches STAR (it only sizes the budget),
    /// so bwa's 32-bit `-K` limit does not apply, and SAM input is accepted.
    #[test]
    fn resolve_accepts_a_large_chunk_size_and_sam_input() {
        let tmp = tempfile::tempdir().unwrap();
        let (mut opts, reference, bin) = star_setup(tmp.path());
        opts.chunk_size = 1 << 32;
        opts.input_format = Some(UnalignedReadFormat::Sam);
        opts.resolve(&reference, 4, Some(&bin)).unwrap();
    }

    #[test]
    fn resolve_fails_fast_on_a_mismatched_index() {
        let tmp = tempfile::tempdir().unwrap();
        let (opts, reference, bin) = star_setup(tmp.path());
        std::fs::write(tmp.path().join("ref.dict"), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:999\n").unwrap();
        let msg = format!("{:#}", opts.resolve(&reference, 4, Some(&bin)).unwrap_err());
        assert!(msg.contains("does not match the sequence dictionary"), "got: {msg}");
    }

    #[test]
    fn resolve_looks_for_star_on_path() {
        let tmp = tempfile::tempdir().unwrap();
        let (opts, reference, _) = star_setup(tmp.path());
        // Only meaningful where no STAR is installed; the message names the
        // binary either way when it is missing.
        if which::which("STAR").is_err() {
            let msg = opts.resolve(&reference, 4, None).unwrap_err().to_string();
            assert!(msg.contains("aligner binary `STAR` not found on PATH"), "got: {msg}");
        }
    }

    #[rstest]
    #[case::index_with_bwa_mem3(Some(AlignerPreset::BwaMem3), true, false)]
    #[case::extra_args_with_bwa(Some(AlignerPreset::Bwa), false, true)]
    #[case::index_with_command(None, true, false)]
    #[case::extra_args_with_inproc(Some(AlignerPreset::BwaMem3InProc), false, true)]
    fn star_options_rejected_without_star(
        #[case] preset: Option<AlignerPreset>,
        #[case] index: bool,
        #[case] extra_args: bool,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let reference = write_reference(tmp.path(), DICT);
        let opts = AlignerOptions {
            preset,
            command: preset.is_none().then(|| "cat {ref}".to_string()),
            index: index.then(|| tmp.path().join("star")),
            extra_args: extra_args.then(|| "--outSAMmapqUnique 60".to_string()),
            ..AlignerOptions::default()
        };
        let msg = opts.resolve(&reference, 4, None).unwrap_err().to_string();
        let flag = if index { "--aligner::index" } else { "--aligner::extra-args" };
        assert!(
            msg.contains(&format!("{flag} is only valid with --aligner::preset star")),
            "{msg}"
        );
    }

    #[rstest]
    fn startup_timeout_rejected_without_star(
        #[values(
            None,
            Some(AlignerPreset::Bwa),
            Some(AlignerPreset::BwaMem3),
            Some(AlignerPreset::BwaMem3InProc)
        )]
        preset: Option<AlignerPreset>,
    ) {
        let opts = AlignerOptions { preset, startup_timeout: Some(1), ..AlignerOptions::default() };
        let error = opts.check_backend_specific_flags().unwrap_err();
        assert!(
            error
                .to_string()
                .contains("--aligner::startup-timeout is only valid with --aligner::preset star")
        );
    }

    #[rstest]
    #[case::star(Some(AlignerPreset::Star), true)]
    #[case::bwa_mem3(Some(AlignerPreset::BwaMem3), false)]
    #[case::command(None, false)]
    fn realign_rejects_only_star(#[case] preset: Option<AlignerPreset>, #[case] rejected: bool) {
        let opts = AlignerOptions { preset, ..AlignerOptions::default() };
        assert_eq!(opts.check_realign_support().is_err(), rejected);
    }

    #[test]
    fn preset_names() {
        use clap::ValueEnum;
        assert_eq!(AlignerPreset::from_str("star", false).unwrap(), AlignerPreset::Star);
        assert_eq!(AlignerPreset::Star.to_string(), "star");
        assert_eq!(AlignerPreset::Star.binary_name(), "STAR");
    }

    /// `--aligner::max-read-ahead` and `--aligner::chunk-size` reach the
    /// resolved star backend.
    #[rstest]
    fn resolve_carries_max_read_ahead_and_chunk_size_for_star(
        #[values(None, Some(8 << 30))] max_read_ahead: Option<u64>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let (mut opts, reference, bin) = star_setup(tmp.path());
        opts.max_read_ahead = max_read_ahead;
        opts.chunk_size = 12_345_678;
        let resolved = opts.resolve(&reference, 4, Some(&bin)).unwrap();
        assert_eq!(resolved.chunk_size, 12_345_678);
        let ResolvedBackend::Subprocess { max_read_ahead: resolved_max, .. } = resolved.backend
        else {
            panic!("star resolves to the subprocess backend");
        };
        assert_eq!(resolved_max, max_read_ahead);
    }

    /// `AlignerPreset::validate` cannot check a STAR index (a separate
    /// directory) and says so instead of passing.
    #[test]
    fn preset_validate_rejects_star() {
        let tmp = tempfile::tempdir().unwrap();
        let (_, reference, bin) = star_setup(tmp.path());
        let msg = AlignerPreset::Star.validate(&reference, Some(&bin)).unwrap_err().to_string();
        assert!(msg.contains("cannot check a star index"), "got: {msg}");
    }

    #[rstest]
    #[case::default(None, 900)]
    #[case::short(Some(1), 1)]
    #[case::long(Some(3600), 3600)]
    #[case::largest(Some(u64::MAX), u64::MAX)]
    fn resolve_carries_star_startup_timeout(#[case] seconds: Option<u64>, #[case] expected: u64) {
        let tmp = tempfile::tempdir().unwrap();
        let (mut opts, reference, bin) = star_setup(tmp.path());
        opts.startup_timeout = seconds;
        let resolved = opts.resolve(&reference, 4, Some(&bin)).unwrap();
        let ResolvedBackend::Subprocess { runtime, .. } = resolved.backend else {
            panic!("star resolves to a subprocess");
        };
        assert_eq!(runtime.startup_timeout, Some(std::time::Duration::from_secs(expected)));
    }

    #[test]
    fn work_dir_is_created_and_removed() {
        let dir = create_work_dir().unwrap();
        let path = dir.path().to_path_buf();
        assert!(path.is_dir());
        assert!(!path.join(STAR_TMP_SUBDIR).exists(), "STAR must create --outTmpDir itself");
        drop(dir);
        assert!(!path.exists(), "the work dir is removed when dropped");
    }
}
