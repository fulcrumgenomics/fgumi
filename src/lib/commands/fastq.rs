//! Convert BAM to FASTQ format.
//!
//! Reads a BAM file and writes FASTQ, either interleaved to stdout (the default,
//! for piping to `bwa mem -p`) or split into per-read files (`--out1`/`--out2`,
//! plus optional `--out0`). Input should be queryname-sorted or
//! template-coordinate sorted. The conversion runs on the typed-step pipeline.

use crate::commands::common::{
    CompressionOptions, MemoryLimit, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
    reject_output_collisions,
};
use crate::sam::SamTag;
use crate::validation::validate_input_exists;
use anyhow::Result;
use clap::Parser;
use fgumi_bam_io::{ReadStreams, is_stdin_path};
use fgumi_raw_bam::{
    RawRecord, aux_data_slice, extract_sequence_into, find_string_tag, quality_scores_slice,
    read_name as raw_read_name,
};
use log::info;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::commands::command::Command;

/// Lookup table for Phred to Phred+33 ASCII conversion (clamped to 126)
static QUAL_TO_ASCII: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let val = (i as u8).saturating_add(33);
        table[i] = if val > 126 { 126 } else { val };
        i += 1;
    }
    table
};

/// Lookup table for base complement (A<->T, C<->G, others->N)
static COMPLEMENT: [u8; 256] = {
    // Values are shared from `fgumi_dna::COMPLEMENT` (IUPAC-aware, case-preserving),
    // but this table keeps the FASTQ-validity policy: any byte outside the valid
    // IUPAC alphabet folds to 'N', since an emitted FASTQ base must be a valid
    // nucleotide code (unlike the consensus/tag paths, which pass unknowns through).
    let mut table = [b'N'; 256];
    const VALID: &[u8] = b"ACGTURYSWKMBVDHNacgturyswkmbvdhn";
    let mut i = 0;
    while i < VALID.len() {
        let base = VALID[i];
        table[base as usize] = fgumi_dna::COMPLEMENT[base as usize];
        i += 1;
    }
    table
};

/// Convert BAM to FASTQ format.
#[derive(Debug, Parser)]
#[command(
    name = "fastq",
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mConvert BAM to FASTQ format\x1b[0m",
    long_about = r#"
Convert a BAM file to FASTQ, either interleaved (the default, for `bwa mem -p`)
or split into per-read files (--out1/--out2, plus optional --out0).

Reads BAM records and, by default, writes interleaved FASTQ to stdout for piping
to aligners. Input should be queryname-sorted or template-coordinate sorted. Any
`.gz`/`.bgz` output path is written as BGZF (gzip-compatible, block-indexable);
stdout is always plain text. The conversion runs on the typed-step pipeline, so
`--threads` parallelizes both BAM decompression and BGZF output compression.

EXAMPLES:

  # Pipe interleaved FASTQ to bwa mem for alignment
  fgumi fastq -i unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y ref.fa -

  # Paired split output to gzipped files (mirrors `samtools fastq -1/-2/-0`)
  fgumi fastq -i unmapped.bam -@ 4 \
    --out1 R1.fastq.gz --out2 R2.fastq.gz --out0 other.fastq.gz

  # Embed the UMI (from the RX tag) in the read name, DRAGEN-ready
  fgumi fastq -i extracted.bam --annotate-read-names \
    --out1 R1.fastq.gz --out2 R2.fastq.gz --out0 /dev/null

  # Exclude secondary and supplementary alignments (default)
  fgumi fastq -i aligned.bam -F 0x900 | bwa mem ...

NOTES:

  Read-name suffixes (/1, /2): appended from the FLAG unless --no-read-suffix,
  matching `samtools fastq -N`. A QNAME that already carries a mate suffix
  (a re-imported BAM; the SAM spec forbids it) is not stripped, so the output
  would double it (`name/1/1`) -- same as samtools -N. Pass --no-read-suffix
  to leave names untouched. Paired split (--out1/--out2) mode always omits the
  suffix, since the R1/R2 file already identifies the mate.

  Missing base qualities: a read with no stored quality (all 0xFF, the SAM
  no-quality sentinel) is emitted with a fixed Q33 ('B') per base, matching
  `samtools fastq` -- not Q93, which would falsely claim near-perfect quality.

  Interleaved pairing: for `bwa mem -p` the input must be queryname or
  template-coordinate sorted with both mates present; a lone mate (its pair
  absent or dropped by --exclude-flags / --require-flags) desyncs the R1/R2
  stream. Prefer paired split (--out1/--out2) when mates may be missing.
"#
)]
pub struct Fastq {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output FASTQ file, or `-` / `/dev/stdout` for stdout. If omitted, the
    /// FASTQ stream is written to stdout (the default, intended for piping
    /// straight to an aligner). A `.gz`/`.bgz` path is BGZF-compressed; stdout
    /// is always plain text.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Don't append /1 and /2 to read names.
    #[arg(short = 'n', long = "no-read-suffix", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub no_suffix: bool,

    /// Exclude reads with any of these flags present [0x900 = secondary|supplementary].
    #[arg(short = 'F', long = "exclude-flags", default_value_t = 0x900, value_parser = parse_flags)]
    pub exclude_flags: u16,

    /// Only include reads with all of these flags present.
    #[arg(short = 'f', long = "require-flags", default_value_t = 0, value_parser = parse_flags)]
    pub require_flags: u16,

    /// Number of threads for BAM decompression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Deprecated and ignored: output batching is managed by the typed-step
    /// pipeline (tune with `--max-memory`). Retained so existing `-K …`
    /// invocations keep parsing; a non-default value logs a deprecation notice.
    #[arg(short = 'K', long = "bwa-chunk-size", default_value_t = DEFAULT_BWA_CHUNK_SIZE, hide = true)]
    pub bwa_chunk_size: u64,

    /// Append the record's UMI to the read name, before any /1 or /2 suffix.
    ///
    /// With the default delimiters this matches `samtools fastq -U`
    /// (`readname:AAAA+CCCC`), the layout DRAGEN expects.
    #[arg(short = 'a', short_alias = 'U', long = "annotate-read-names", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub annotate_read_names: bool,

    /// Tags to read the UMI from, in priority order; the first present wins.
    #[arg(long = "umi-tag", default_value = "RX,OX", value_delimiter = ',')]
    pub umi_tag: Vec<String>,

    /// Delimiter between the read name and the UMI.
    #[arg(long = "umi-name-delim", default_value = ":")]
    pub umi_name_delim: String,

    /// Separator between the two halves of a duplex UMI in the read name.
    ///
    /// fgumi stores duplex UMIs as `AAAA-CCCC`; `samtools fastq -U` and DRAGEN
    /// expect `AAAA+CCCC`, so the stored `-` is rewritten to this value.
    #[arg(long = "umi-sep", default_value = "+")]
    pub umi_sep: String,

    /// Write read 1 (R1) to this file instead of the interleaved stream; requires
    /// `--out2`. A `.gz`/`.bgz` path is written as BGZF. Mirrors `samtools fastq -1`.
    #[arg(short = '1', long = "out1", requires = "out2", conflicts_with = "output")]
    pub out1: Option<PathBuf>,

    /// Write read 2 (R2) to this file; requires `--out1`. A `.gz`/`.bgz` path is
    /// written as BGZF. Mirrors `samtools fastq -2`.
    #[arg(short = '2', long = "out2", requires = "out1", conflicts_with = "output")]
    pub out2: Option<PathBuf>,

    /// Write reads that are neither cleanly R1 nor R2 (single-end / ambiguous)
    /// here; requires `--out1`. If omitted, such reads go to stdout, matching
    /// `samtools fastq` without `-0`. A `.gz`/`.bgz` path is written as BGZF.
    #[arg(short = '0', long = "out0", requires = "out1")]
    pub out0: Option<PathBuf>,

    /// Pipeline scheduler diagnostics (`--pipeline-stats`, deadlock detection).
    #[command(flatten)]
    pub scheduler: SchedulerOptions,

    /// Pipeline queue-memory limits (`--max-memory`, …). The conversion runs on
    /// the typed-step pipeline; these cap the in-flight queue memory. Defaults to
    /// a lean fastq budget (see `FASTQ_DEFAULT_QUEUE_MEMORY_MB`); pass
    /// `--max-memory` to override.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

/// Parse flag values supporting both decimal and hex (0x) notation.
fn parse_flags(s: &str) -> Result<u16, String> {
    if s.starts_with("0x") || s.starts_with("0X") {
        u16::from_str_radix(&s[2..], 16).map_err(|e| e.to_string())
    } else {
        s.parse().map_err(|e: std::num::ParseIntError| e.to_string())
    }
}

impl Fastq {
    /// Build the optional UMI-in-read-name annotation, validating the tag list
    /// up front. Returns `None` when `--annotate-read-names` is off.
    fn build_umi_header(&self) -> Result<Option<UmiNameAnnotation>> {
        if self.annotate_read_names {
            info!("UMI in read name: from tag(s) {}", self.umi_tag.join(","));
            Ok(Some(UmiNameAnnotation::new(&self.umi_tag, &self.umi_name_delim, &self.umi_sep)?))
        } else {
            Ok(None)
        }
    }

    /// Log the shared conversion configuration once at run start.
    ///
    /// Paired split output always omits the `/1` `/2` suffix (the R1/R2 file
    /// already identifies the mate), so the suffix line reports the *effective*
    /// behavior — not the raw `--no-read-suffix` flag, which paired mode ignores.
    fn log_config(&self) {
        info!("Input: {}", self.input.display());
        info!("Threads: {}", self.threads);
        info!("Exclude flags: 0x{:X}", self.exclude_flags);
        info!("Require flags: 0x{:X}", self.require_flags);
        let suffix = if self.out1.is_some() {
            "omitted (paired split)"
        } else if self.no_suffix {
            "disabled"
        } else {
            "enabled"
        };
        info!("Read name suffix: {suffix}");
    }

    /// Build the per-stage [`FastqOptions`] — the single source of truth for the
    /// flag filters, read-name suffix behavior, and UMI-in-read-name config the
    /// encode step needs, shared by the interleaved and paired-split chain specs.
    fn to_fastq_options(&self) -> Result<FastqOptions> {
        Ok(FastqOptions {
            exclude_flags: self.exclude_flags,
            require_flags: self.require_flags,
            no_suffix: self.no_suffix,
            umi_header: self.build_umi_header()?,
        })
    }

    /// Assemble the typed-step [`ChainSpec`] for a BAM → FASTQ conversion:
    /// `SourceSpec::Bam` → `Stage::Fastq` → the given `sink`.
    ///
    /// [`ChainSpec`]: crate::pipeline::chains::ChainSpec
    fn build_chain_spec(
        &self,
        sink: crate::pipeline::chains::SinkSpec,
        command_line: &str,
    ) -> Result<crate::pipeline::chains::ChainSpec> {
        use crate::pipeline::chains::{ChainSpec, SourceSpec, Stage, StageOptionsBag};

        Ok(ChainSpec {
            stages: vec![Stage::Fastq],
            source: SourceSpec::Bam(self.input.clone()),
            sink,
            stage_opts: StageOptionsBag {
                fastq: Some(self.to_fastq_options()?),
                ..Default::default()
            },
            threading: ThreadingOptions::new(self.threads),
            compression: CompressionOptions { compression_level: FASTQ_GZIP_LEVEL },
            scheduler: self.scheduler.clone(),
            queue_memory: self.resolve_queue_memory(),
            async_reader: false,
            // fastq exposes no read-stream knob; keep the plain sequential reader.
            read_streams: ReadStreams::Fixed(1),
            // Match the default CRC policy every other command uses: verify a file
            // source, skip for stdin (which cannot be re-decoded to re-check).
            verify_crc: !is_stdin_path(&self.input),
            command_line: command_line.to_string(),
        })
    }

    /// Assemble the [`ChainSpec`] for paired-split BAM → FASTQ:
    /// `SourceSpec::Bam` → `Stage::Fastq` → `SinkSpec::FastqPaired`. Only the
    /// sink differs from [`Self::build_chain_spec`]; every other field is shared.
    ///
    /// [`ChainSpec`]: crate::pipeline::chains::ChainSpec
    fn build_paired_chain_spec(
        &self,
        command_line: &str,
    ) -> Result<crate::pipeline::chains::ChainSpec> {
        // clap `requires` guarantees out2 is present whenever out1 is, and
        // execute() only calls this when out1.is_some().
        let out1 = self.out1.clone().expect("out1 present in paired mode");
        let out2 = self.out2.clone().expect("out2 present in paired mode");
        let sink =
            crate::pipeline::chains::SinkSpec::FastqPaired { out1, out2, out0: self.out0.clone() };
        self.build_chain_spec(sink, command_line)
    }

    /// Resolve the queue-memory options for the chain.
    ///
    /// FASTQ blocks are small and the chain is a short linear
    /// source → encode → write pipeline, so the sort/consensus-sized default
    /// (768 MiB *per thread*) is wildly oversized and would inflate RSS as
    /// `--threads` scales. When the flags are left at that default, a lean fixed
    /// total ([`FASTQ_DEFAULT_QUEUE_MEMORY_MB`]) is substituted; any explicit
    /// `--max-memory` / `--memory-per-thread` override is honored untouched. (A
    /// user who explicitly passes the default value is indistinguishable from not
    /// passing it and gets the lean budget — harmless, since 768 MiB/thread is
    /// never what a fastq conversion wants.)
    ///
    /// The "is it the default?" test compares against [`QueueMemoryOptions::default`]
    /// rather than a hardcoded byte count, so it cannot silently stop firing if
    /// that default is ever changed in `common.rs`.
    fn resolve_queue_memory(&self) -> QueueMemoryOptions {
        let mut qm = self.queue_memory.clone();
        let default = QueueMemoryOptions::default();
        let is_default = matches!(
            (&qm.max_memory, &default.max_memory),
            (MemoryLimit::Fixed(a), MemoryLimit::Fixed(b)) if a == b
        ) && qm.memory_per_thread == default.memory_per_thread;
        if is_default {
            qm.max_memory =
                MemoryLimit::Fixed(FASTQ_DEFAULT_QUEUE_MEMORY_MB as usize * 1024 * 1024);
            qm.memory_per_thread = false;
        }
        qm
    }

    /// Reject a write target that resolves to the same file as the `--input`
    /// BAM, which would clobber the file being read (mirrors `retag`/`copy_umi`).
    ///
    /// The input BAM always exists, so it canonicalises; any output that
    /// canonicalises to the same path is rejected before a writer truncates it.
    /// This also catches a symlinked or `./`-spelled output that resolves to the
    /// input. Stdin input (`-`/`/dev/stdin`) has no filesystem entity to
    /// canonicalise and cannot be clobbered, so the guard no-ops there.
    ///
    /// Output-vs-output collisions — two paths naming one destination, and the
    /// stdout-multiplexing case — are handled separately by
    /// [`reject_output_collisions`], which shares the same `(path, flag)` slice.
    fn reject_write_aliasing_input(&self, outputs: &[(&Path, &str)]) -> Result<()> {
        let Ok(input_canon) = std::fs::canonicalize(&self.input) else {
            return Ok(());
        };
        for (path, flag) in outputs {
            if std::fs::canonicalize(path).is_ok_and(|canon| canon == input_canon) {
                anyhow::bail!(
                    "{flag} '{}' is the same file as --input '{}'; choose a different path",
                    path.display(),
                    self.input.display()
                );
            }
        }
        Ok(())
    }
}

impl Command for Fastq {
    fn execute(&self, command_line: &str) -> Result<()> {
        validate_input_exists(&self.input, "Input BAM")?;

        if self.bwa_chunk_size != DEFAULT_BWA_CHUNK_SIZE {
            info!(
                "--bwa-chunk-size ({}) is deprecated and ignored; output batching is now managed \
                 by the pipeline (tune with --max-memory)",
                self.bwa_chunk_size
            );
        }

        // Refuse to clobber the input BAM or route two streams to the same file
        // before opening anything (a sink truncates its path on create). Only
        // user-specified outputs are checked; the default interleaved stdout
        // (`self.output == None`) names no file and is skipped. `reject_output_collisions`
        // handles output-vs-output (including stdout multiplexing and `./`/symlink
        // aliases, with `/dev/null` exempt); `reject_write_aliasing_input` handles
        // output-vs-input.
        let mut outputs: Vec<(&Path, &str)> = Vec::new();
        if let Some(p) = &self.output {
            outputs.push((p.as_path(), "--output"));
        }
        if let Some(p) = &self.out1 {
            outputs.push((p.as_path(), "--out1"));
        }
        if let Some(p) = &self.out2 {
            outputs.push((p.as_path(), "--out2"));
        }
        if let Some(p) = &self.out0 {
            outputs.push((p.as_path(), "--out0"));
        }
        // Paired mode with no `--out0` routes "other" (single-end / ambiguous)
        // reads to stdout (`SinkSpec::FastqPaired`), so register that implicit
        // stdout target; otherwise `--out1 -` (or `--out2 -`) would multiplex
        // onto the same stdout the "other" stream uses and slip past the guard.
        if self.out1.is_some() && self.out0.is_none() {
            outputs.push((Path::new("-"), "--out0 (default: stdout)"));
        }
        reject_output_collisions(&outputs)?;
        self.reject_write_aliasing_input(&outputs)?;

        self.log_config();

        // Paired split output (`-1`/`-2`, optionally `-0`) runs through the
        // typed-step chain: BAM source → Stage::Fastq 3-way encode → three
        // per-file raw writers (BGZF for `.gz`), all sharing one work-stealing
        // pool. clap guarantees `--out1`/`--out2` come together and conflict with
        // `--output`, so `out1.is_some()` is the paired-mode discriminant.
        if let Some(out1) = &self.out1 {
            let out2 = self.out2.as_ref().expect("out2 present in paired mode");
            info!("R1 -> {}", out1.display());
            info!("R2 -> {}", out2.display());
            match &self.out0 {
                Some(out0) => info!("other -> {}", out0.display()),
                None => info!("other -> stdout"),
            }
            return crate::pipeline::chains::build_for(
                self.build_paired_chain_spec(command_line)?,
            )?
            .run();
        }

        // Interleaved output runs through the same chain: BAM source →
        // Stage::Fastq encode (pool-spread) → one raw writer (BGZF for `.gz`).
        // `-` (or an omitted `--output`) writes to stdout.
        let output = self.output.clone().unwrap_or_else(|| PathBuf::from("-"));
        info!("Output -> {}", output.display());
        let sink = crate::pipeline::chains::SinkSpec::Fastq(output);
        crate::pipeline::chains::build_for(self.build_chain_spec(sink, command_line)?)?.run()
    }
}

/// ASCII quality emitted for a base whose quality is entirely absent.
///
/// Per the SAM spec, a record with no quality stores `0xFF` for every base.
/// fgumi previously mapped `0xFF` through [`QUAL_TO_ASCII`] to `~` (Q93), which
/// falsely claims near-perfect quality for a read that has none. Instead emit
/// Q33 (`B`) for an entirely-absent quality string, matching `samtools fastq`
/// (which `fgumi fastq` is modeled on).
const MISSING_QUALITY_ASCII: u8 = b'B';

/// Encode a record's raw quality bytes as Phred+33 ASCII into `out`.
///
/// When the quality is entirely absent (all `0xFF`), fills with
/// [`MISSING_QUALITY_ASCII`] rather than the misleading Q93 that
/// `QUAL_TO_ASCII[0xFF]` would produce. A partially-`0xFF` string is left to the
/// per-byte mapping (a genuinely malformed record, not the "no quality" sentinel).
fn encode_quality_into(quals: &[u8], out: &mut Vec<u8>) {
    out.clear();
    if !quals.is_empty() && quals.iter().all(|&q| q == 0xFF) {
        out.resize(quals.len(), MISSING_QUALITY_ASCII);
    } else {
        out.extend(quals.iter().map(|&s| QUAL_TO_ASCII[s as usize]));
    }
}

/// Default `--bwa-chunk-size`; kept only so a non-default value can be detected
/// and a deprecation notice logged. Output batching is managed by the pipeline
/// now, so the value is otherwise ignored (see [`Fastq::execute`]).
const DEFAULT_BWA_CHUNK_SIZE: u64 = 150_000_000;

/// BGZF compression level for `.gz`/`.bgz`/`.bgzf` FASTQ output.
///
/// 6 is the zlib/bgzip default: a middle trade-off between output size and CPU,
/// appropriate for an intermediate file handed straight to an aligner. Carried
/// on the chain's [`CompressionOptions`] and applied by the sink's `BgzfCompress`
/// step.
const FASTQ_GZIP_LEVEL: u32 = 6;

/// Lean default total in-flight queue-memory budget (MiB) for `fgumi fastq` when
/// the user leaves `--max-memory` at its per-thread default.
///
/// FASTQ blocks are small and the chain is a short linear
/// source → encode → write pipeline, so the sort/consensus-sized default
/// (768 MiB/thread) is wildly oversized; a modest fixed total keeps RSS flat as
/// `--threads` scales while staying well above the pipeline's steady-state
/// working set. See [`Fastq::resolve_queue_memory`].
const FASTQ_DEFAULT_QUEUE_MEMORY_MB: u64 = 256;

/// How to append a record's UMI to its read name.
///
/// Built once per run from the CLI options so the per-record path does no
/// string parsing or allocation beyond the name buffer itself.
#[derive(Debug, Clone)]
pub struct UmiNameAnnotation {
    /// Tags to consult, in priority order; the first present wins.
    tags: Vec<SamTag>,
    /// Delimiter between the read name and the UMI.
    name_delim: Vec<u8>,
    /// Separator written between duplex UMI halves (replaces the stored `-`).
    umi_sep: Vec<u8>,
}

impl UmiNameAnnotation {
    /// Builds the annotation config from raw CLI strings.
    ///
    /// # Errors
    ///
    /// Returns an error if a tag name is not exactly two characters.
    pub fn new(tags: &[String], name_delim: &str, umi_sep: &str) -> Result<Self> {
        let tags = tags
            .iter()
            .map(|t| {
                let bytes = t.as_bytes();
                if bytes.len() != 2 || !SamTag::is_valid_tag_bytes(bytes[0], bytes[1]) {
                    anyhow::bail!(
                        "--umi-tag values must be a valid two-character SAM tag, got '{t}'"
                    );
                }
                Ok(SamTag::new(bytes[0], bytes[1]))
            })
            .collect::<Result<Vec<_>>>()?;
        if tags.is_empty() {
            anyhow::bail!("--umi-tag requires at least one tag");
        }
        Ok(Self {
            tags,
            name_delim: name_delim.as_bytes().to_vec(),
            umi_sep: umi_sep.as_bytes().to_vec(),
        })
    }

    /// Appends `<delim><umi>` to `out` if the record carries one of the configured
    /// tags, rewriting the stored duplex `-` separator to the configured one.
    ///
    /// Records without any of the tags are left unannotated rather than failing —
    /// a BAM can legitimately mix UMI-bearing and UMI-free reads.
    fn append_to(&self, record: &RawRecord, out: &mut Vec<u8>) {
        let aux = aux_data_slice(record);
        let Some(umi) = self.tags.iter().find_map(|&tag| find_string_tag(aux, tag)) else {
            return;
        };
        out.extend_from_slice(&self.name_delim);
        for &byte in umi {
            if byte == b'-' {
                out.extend_from_slice(&self.umi_sep);
            } else {
                out.push(byte);
            }
        }
    }
}

/// Reusable per-record scratch buffers, kept across the conversion loop so the
/// hot path does not allocate.
#[derive(Debug, Default)]
pub(crate) struct FastqRecordBuffers {
    /// Decoded sequence bases.
    seq: Vec<u8>,
    /// Phred+33 encoded qualities.
    qual: Vec<u8>,
    /// Read name, used only when a UMI is appended.
    name: Vec<u8>,
}

impl FastqRecordBuffers {
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            seq: Vec::with_capacity(capacity),
            qual: Vec::with_capacity(capacity),
            name: Vec::with_capacity(256),
        }
    }
}

/// Write a single FASTQ record to the writer.
#[inline]
pub(crate) fn write_fastq_record<W: Write>(
    writer: &mut W,
    record: &RawRecord,
    flags: u16,
    no_suffix: bool,
    buffers: &mut FastqRecordBuffers,
    umi_annotation: Option<&UmiNameAnnotation>,
) -> Result<()> {
    use fgumi_raw_bam::flags as flag_bits;

    // Get read name (without null terminator)
    let name = raw_read_name(record);

    // Determine read suffix (/1 or /2)
    let is_first = (flags & flag_bits::FIRST_SEGMENT) != 0;
    let is_last = (flags & flag_bits::LAST_SEGMENT) != 0;
    let suffix: &[u8] = if no_suffix {
        b""
    } else if is_first && !is_last {
        b"/1"
    } else if is_last && !is_first {
        b"/2"
    } else {
        b"" // Single-end or both flags set
    };

    // Decode sequence from 4-bit BAM encoding to ASCII bases (reuses the scratch buffer).
    extract_sequence_into(record, &mut buffers.seq);

    // Copy quality bytes and transform to Phred+33 ASCII (absent quality → default)
    encode_quality_into(quality_scores_slice(record), &mut buffers.qual);

    if (flags & flag_bits::REVERSE) != 0 {
        // Reverse complement sequence in place using lookup table
        buffers.seq.reverse();
        for base in buffers.seq.iter_mut() {
            *base = COMPLEMENT[*base as usize];
        }
        // Reverse quality in place
        buffers.qual.reverse();
    }

    // Write all parts. The UMI goes between the name and the /1 /2 suffix, matching
    // `samtools fastq -U`.
    writer.write_all(b"@")?;
    if let Some(annotation) = umi_annotation {
        buffers.name.clear();
        buffers.name.extend_from_slice(name);
        annotation.append_to(record, &mut buffers.name);
        writer.write_all(&buffers.name)?;
    } else {
        writer.write_all(name)?;
    }
    writer.write_all(suffix)?;
    writer.write_all(b"\n")?;
    writer.write_all(&buffers.seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&buffers.qual)?;
    writer.write_all(b"\n")?;

    Ok(())
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Which FASTQ stream a record belongs to, based on its segment flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Segment {
    /// First segment of a pair (`FIRST_SEGMENT` set, `LAST_SEGMENT` clear) → R1.
    Read1,
    /// Last segment of a pair (`LAST_SEGMENT` set, `FIRST_SEGMENT` clear) → R2.
    Read2,
    /// Neither or both segment bits set (single-end or ambiguous) → "other".
    Other,
}

/// Classify a record into R1 / R2 / other from its flags, matching how
/// `samtools fastq` routes reads to `-1` / `-2` / `-0`.
pub(crate) fn classify_segment(flags: u16) -> Segment {
    use fgumi_raw_bam::flags as flag_bits;
    let is_first = (flags & flag_bits::FIRST_SEGMENT) != 0;
    let is_last = (flags & flag_bits::LAST_SEGMENT) != 0;
    match (is_first, is_last) {
        (true, false) => Segment::Read1,
        (false, true) => Segment::Read2,
        _ => Segment::Other,
    }
}

/// Returns `true` if `path` names a gzip-family output that must be written as
/// BGZF — extension `gz`/`bgz`/`bgzf`, matched case-insensitively.
///
/// The single gzip-detection function for `fgumi fastq`: it decides both the
/// interleaved `-o` sink and each paired `-1`/`-2`/`-0` sink (through the chain
/// builder's `wire_fastq_output`). Case-insensitive so `OUT.FQ.GZ` is still
/// compressed rather than written as plain text under a `.GZ` name.
pub(crate) fn path_is_gzip(path: &Path) -> bool {
    path.extension().and_then(|e| e.to_str()).is_some_and(|e| {
        e.eq_ignore_ascii_case("gz")
            || e.eq_ignore_ascii_case("bgz")
            || e.eq_ignore_ascii_case("bgzf")
    })
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Per-stage options for [`crate::pipeline::chains::Stage::Fastq`] — the
/// BAM→FASTQ encode. Carries the flag filters, read-name suffix behavior, and
/// the optional UMI-in-read-name config the encode step needs. The output
/// destination(s) and compression are carried by the chain's `SinkSpec`, not
/// here.
#[derive(Clone)]
pub struct FastqOptions {
    /// Exclude reads with any of these flag bits set.
    pub exclude_flags: u16,
    /// Only include reads with all of these flag bits set.
    pub require_flags: u16,
    /// When `true`, omit the `/1` `/2` read-name suffix.
    pub no_suffix: bool,
    /// When `Some`, append the UMI (from the configured tag) to the read name.
    pub umi_header: Option<UmiNameAnnotation>,
}

// ==== impls ported from feat-runall for the chain builder (R2) ====
impl FastqOptions {
    /// Returns `true` if a record with `flags` passes the include/exclude filters.
    #[must_use]
    pub(crate) fn passes_filters(&self, flags: u16) -> bool {
        (flags & self.exclude_flags) == 0 && (flags & self.require_flags) == self.require_flags
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// Build a single-end record carrying an optional UMI tag.
    fn record_with_umi(name: &[u8], umi: Option<(&str, &[u8])>) -> RawRecord {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(name);
        if let Some((tag, value)) = umi {
            let t = tag.as_bytes();
            b.add_string_tag(SamTag::new(t[0], t[1]), value);
        }
        b.build()
    }

    #[rstest]
    // Default delimiters reproduce `samtools fastq -U`: name:UMI, duplex `-` -> `+`.
    #[case::simplex_umi(Some(("RX", &b"ACGT"[..])), ":", "+", "read1:ACGT")]
    #[case::duplex_umi(Some(("RX", &b"ACGT-TTTT"[..])), ":", "+", "read1:ACGT+TTTT")]
    #[case::custom_delims(Some(("RX", &b"ACGT-TTTT"[..])), "_", "|", "read1_ACGT|TTTT")]
    // A record without any of the configured tags is left unannotated rather than failing.
    #[case::no_umi_tag(None, ":", "+", "read1")]
    // OX is consulted when RX is absent (tag priority order).
    #[case::fallback_tag(Some(("OX", &b"GGGG"[..])), ":", "+", "read1:GGGG")]
    fn test_umi_name_annotation(
        #[case] umi: Option<(&str, &[u8])>,
        #[case] name_delim: &str,
        #[case] umi_sep: &str,
        #[case] expected: &str,
    ) {
        let record = record_with_umi(b"read1", umi);
        let annotation =
            UmiNameAnnotation::new(&["RX".to_string(), "OX".to_string()], name_delim, umi_sep)
                .expect("valid annotation config");

        let mut out = b"read1".to_vec();
        annotation.append_to(&record, &mut out);
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[rstest]
    #[case::too_short("R")]
    #[case::too_long("RXX")]
    #[case::invalid_first_char("1X")]
    fn test_umi_name_annotation_rejects_bad_tag(#[case] tag: &str) {
        assert!(UmiNameAnnotation::new(&[tag.to_string()], ":", "+").is_err());
    }

    #[test]
    fn test_umi_name_annotation_rejects_empty_tag_list() {
        assert!(UmiNameAnnotation::new(&[], ":", "+").is_err());
    }

    /// Write reverse complement of sequence bytes to a buffer (test helper).
    fn write_reverse_complement_bytes<W: Write>(writer: &mut W, bases: &[u8]) -> Result<()> {
        for &base in bases.iter().rev() {
            let comp = match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'N' | b'n' => b'N',
                _ => b'N',
            };
            writer.write_all(&[comp])?;
        }
        Ok(())
    }

    /// Write quality scores as Phred+33 ASCII (test helper).
    fn write_quality_bytes<W: Write>(writer: &mut W, quals: &[u8]) -> Result<()> {
        for &score in quals {
            let ascii = score.saturating_add(33).min(126);
            writer.write_all(&[ascii])?;
        }
        Ok(())
    }

    /// Write reversed quality scores as Phred+33 ASCII (test helper).
    fn write_reversed_quality_bytes<W: Write>(writer: &mut W, quals: &[u8]) -> Result<()> {
        for &score in quals.iter().rev() {
            let ascii = score.saturating_add(33).min(126);
            writer.write_all(&[ascii])?;
        }
        Ok(())
    }

    #[test]
    fn test_parse_flags_decimal() {
        assert_eq!(parse_flags("2304").expect("parse decimal '2304' should succeed"), 2304);
        assert_eq!(parse_flags("0").expect("parse decimal '0' should succeed"), 0);
        assert_eq!(parse_flags("65535").expect("parse decimal '65535' should succeed"), 65535);
    }

    #[test]
    fn test_parse_flags_hex() {
        assert_eq!(parse_flags("0x900").expect("parse hex '0x900' should succeed"), 0x900);
        assert_eq!(parse_flags("0X900").expect("parse hex '0X900' should succeed"), 0x900);
        assert_eq!(parse_flags("0xff").expect("parse hex '0xff' should succeed"), 0xff);
        assert_eq!(parse_flags("0xFFFF").expect("parse hex '0xFFFF' should succeed"), 0xFFFF);
    }

    #[test]
    fn test_parse_flags_invalid() {
        assert!(parse_flags("invalid").is_err());
        assert!(parse_flags("0xGGGG").is_err());
        assert!(parse_flags("-1").is_err());
    }

    #[test]
    fn test_qual_to_ascii_lookup_table() {
        // Quality 0 -> ASCII 33 ('!')
        assert_eq!(QUAL_TO_ASCII[0], 33);
        // Quality 30 -> ASCII 63 ('?')
        assert_eq!(QUAL_TO_ASCII[30], 63);
        // Quality 40 -> ASCII 73 ('I')
        assert_eq!(QUAL_TO_ASCII[40], 73);
        // Quality 93 -> ASCII 126 ('~') - max valid
        assert_eq!(QUAL_TO_ASCII[93], 126);
        // Quality > 93 should be clamped to 126
        assert_eq!(QUAL_TO_ASCII[94], 126);
        assert_eq!(QUAL_TO_ASCII[255], 126);
    }

    #[test]
    fn test_complement_lookup_table() {
        // Standard bases
        assert_eq!(COMPLEMENT[b'A' as usize], b'T');
        assert_eq!(COMPLEMENT[b'T' as usize], b'A');
        assert_eq!(COMPLEMENT[b'C' as usize], b'G');
        assert_eq!(COMPLEMENT[b'G' as usize], b'C');
        assert_eq!(COMPLEMENT[b'N' as usize], b'N');

        // Lowercase (case preserved)
        assert_eq!(COMPLEMENT[b'a' as usize], b't');
        assert_eq!(COMPLEMENT[b't' as usize], b'a');
        assert_eq!(COMPLEMENT[b'c' as usize], b'g');
        assert_eq!(COMPLEMENT[b'g' as usize], b'c');
        assert_eq!(COMPLEMENT[b'n' as usize], b'n'); // case preserved (shared table)

        // Unknown bases map to N
        assert_eq!(COMPLEMENT[b'X' as usize], b'N');
        assert_eq!(COMPLEMENT[0], b'N');
    }

    #[test]
    fn test_complement_lookup_table_iupac_and_case() {
        // IUPAC-aware (values shared from fgumi_dna::COMPLEMENT), case-preserving,
        // but keeps the FASTQ-validity policy of folding invalid bytes to N.
        assert_eq!(COMPLEMENT[b'R' as usize], b'Y');
        assert_eq!(COMPLEMENT[b'Y' as usize], b'R');
        assert_eq!(COMPLEMENT[b'K' as usize], b'M');
        assert_eq!(COMPLEMENT[b'S' as usize], b'S');
        assert_eq!(COMPLEMENT[b'a' as usize], b't'); // case preserved
        assert_eq!(COMPLEMENT[b'r' as usize], b'y');
        assert_eq!(COMPLEMENT[b'.' as usize], b'N'); // invalid -> N (FASTQ validity)
    }

    #[test]
    fn test_write_reverse_complement_bytes() {
        let mut output = Vec::new();
        // ACGT reversed = TGCA, complemented = ACGT
        write_reverse_complement_bytes(&mut output, b"ACGT")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"ACGT");

        output.clear();
        write_reverse_complement_bytes(&mut output, b"AAAA")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"TTTT");

        output.clear();
        // ATCG reversed = GCTA, complemented = CGAT
        write_reverse_complement_bytes(&mut output, b"ATCG")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"CGAT");

        output.clear();
        // Test with N
        write_reverse_complement_bytes(&mut output, b"ANCG")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"CGNT");
    }

    /// FASTQ3-03: a record whose quality is entirely absent (all `0xFF`, the SAM
    /// no-quality sentinel) must not be emitted as `~` (Q93, near-perfect); it
    /// gets the `MISSING_QUALITY_ASCII` default (`B`, Q33, matching samtools). A
    /// present quality maps per-byte; a partially-`0xFF` (malformed) string is NOT
    /// the sentinel and is left to the map; empty stays empty.
    #[rstest]
    #[case::absent_all_0xff(&[0xFF, 0xFF, 0xFF, 0xFF], b"BBBB")]
    #[case::present_phred33(&[0, 30, 40, 93], &[33, 63, 73, 126])]
    #[case::partial_0xff_is_not_sentinel(&[30, 0xFF], &[63, 126])]
    #[case::empty(&[], b"")]
    fn encode_quality_maps_present_and_defaults_absent(
        #[case] quals: &[u8],
        #[case] expected: &[u8],
    ) {
        let mut out = Vec::new();
        encode_quality_into(quals, &mut out);
        assert_eq!(out.as_slice(), expected);
    }

    #[test]
    fn test_write_quality_bytes() {
        let mut output = Vec::new();
        // Quality 0 -> ASCII 33 ('!')
        // Quality 30 -> ASCII 63 ('?')
        write_quality_bytes(&mut output, &[0, 30, 40]).expect("write_quality_bytes should succeed");
        assert_eq!(output, vec![33, 63, 73]);
    }

    #[test]
    fn test_write_reversed_quality_bytes() {
        let mut output = Vec::new();
        write_reversed_quality_bytes(&mut output, &[0, 30, 40])
            .expect("write_reversed_quality_bytes should succeed");
        // Reversed: [40, 30, 0] -> [73, 63, 33]
        assert_eq!(output, vec![73, 63, 33]);
    }

    // Output-vs-output collision detection is provided by the reused
    // `reject_output_collisions` helper (tested in `common.rs`), and the
    // output-vs-input clobber guard (`reject_write_aliasing_input`) is exercised
    // end-to-end by the integration tests (`test_fastq_output_same_as_input_rejected`,
    // `_symlink_to_input_rejected`, `_paired_duplicate_output_rejected`), so no
    // unit test duplicates them here.

    #[test]
    fn test_quality_encoding_edge_cases() {
        let mut output = Vec::new();
        // Test max valid quality (93)
        write_quality_bytes(&mut output, &[93]).expect("write_quality_bytes should succeed");
        assert_eq!(output, vec![126]); // '~'

        output.clear();
        // Test overflow clamping (94+ should clamp to 126)
        write_quality_bytes(&mut output, &[94, 100, 255])
            .expect("write_quality_bytes should succeed");
        assert_eq!(output, vec![126, 126, 126]);
    }

    /// `classify_segment` routes by the FIRST/LAST segment bits exactly as
    /// `samtools fastq` does: first-only → R1, last-only → R2, neither/both → other.
    #[rstest]
    #[case::read1(fgumi_raw_bam::flags::FIRST_SEGMENT, Segment::Read1)]
    #[case::read2(fgumi_raw_bam::flags::LAST_SEGMENT, Segment::Read2)]
    #[case::both_bits(
        fgumi_raw_bam::flags::FIRST_SEGMENT | fgumi_raw_bam::flags::LAST_SEGMENT,
        Segment::Other
    )]
    #[case::neither(0, Segment::Other)]
    fn classify_segment_routes_by_first_last_bits(#[case] flags: u16, #[case] expected: Segment) {
        assert_eq!(classify_segment(flags), expected);
    }

    /// `passes_filters` requires ALL `require_flags` set and NONE of
    /// `exclude_flags` set (an inverted bit-mask here would pass silently).
    #[rstest]
    // exclude FIRST_SEGMENT: an R1 read is excluded, others pass.
    #[case::excluded(
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        0,
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        false
    )]
    #[case::not_excluded(fgumi_raw_bam::flags::FIRST_SEGMENT, 0, 0, true)]
    // require PAIRED: only paired reads pass.
    #[case::required_present(0, fgumi_raw_bam::flags::PAIRED, fgumi_raw_bam::flags::PAIRED, true)]
    #[case::required_absent(0, fgumi_raw_bam::flags::PAIRED, 0, false)]
    // both: must satisfy require AND avoid exclude.
    #[case::require_met_exclude_hit(
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        fgumi_raw_bam::flags::PAIRED,
        fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
        false
    )]
    #[case::require_met_exclude_clear(
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        fgumi_raw_bam::flags::PAIRED,
        fgumi_raw_bam::flags::PAIRED,
        true
    )]
    fn passes_filters_applies_exclude_and_require(
        #[case] exclude_flags: u16,
        #[case] require_flags: u16,
        #[case] flags: u16,
        #[case] expected: bool,
    ) {
        let opts =
            FastqOptions { exclude_flags, require_flags, no_suffix: false, umi_header: None };
        assert_eq!(opts.passes_filters(flags), expected);
    }

    /// `path_is_gzip` matches `gz`/`bgz`/`bgzf` as the final extension,
    /// case-insensitively (so `OUT.FQ.GZ` is still compressed), and does not
    /// match `.gz` mid-stem.
    #[rstest]
    #[case::gz("out.fq.gz", true)]
    #[case::gz_upper("OUT.FQ.GZ", true)]
    #[case::bgz("out.fq.bgz", true)]
    #[case::bgzf("out.fq.bgzf", true)]
    #[case::bgzf_mixed_case("out.fq.BgzF", true)]
    #[case::plain_fq("out.fq", false)]
    #[case::no_extension("out", false)]
    #[case::gz_in_stem("out.gz.fq", false)]
    fn path_is_gzip_matches_gz_bgz_bgzf(#[case] path: &str, #[case] expected: bool) {
        assert_eq!(path_is_gzip(std::path::Path::new(path)), expected);
    }
}
