//! # CODEC Consensus Calling Command
//!
//! Implementation of the `codec` command for calling consensus from CODEC sequencing data.
//!
//! CODEC (Bae et al 2023) is a sequencing protocol where each read-pair sequences both strands
//! of the original duplex molecule. R1 comes from one strand, R2 from the opposite strand,
//! allowing even a single read-pair to generate duplex consensus.

use crate::commands::command::Command;
use anyhow::{Result, bail};
use clap::Parser;

use super::common::{
    BamIoOptions, CompressionOptions, QueueMemoryOptions, ReadGroupOptions, RejectsOptions,
    SchedulerOptions, StatsOptions, ThreadingOptions,
};

/// Codec-specific tuning, flattened into both the standalone `Codec`
/// command (as bare `--min-duplex-length`, `--outer-bases-qual`, …) and
/// the `RunAll` command (as `--codec::min-duplex-length`, etc. via the
/// `MultiCodecOptions` companion struct).
///
/// `Default` matches each field's clap `default_value` so the macro's
/// generated `default_value_t = CodecOptions::default().field` reads
/// the same value the standalone command's CLI default would yield.
#[fgumi_cli_macros::multi_options("codec", "Codec Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct CodecOptions {
    // ── Consensus-calling options (inlined from ConsensusCallingOptions). ──
    /// Phred-scaled error rate prior to UMI integration
    #[arg(short = '1', long = "error-rate-pre-umi", default_value = "45")]
    pub error_rate_pre_umi: u8,

    /// Phred-scaled error rate post UMI integration
    #[arg(short = '2', long = "error-rate-post-umi", default_value = "40")]
    pub error_rate_post_umi: u8,

    /// Minimum base quality in raw reads to use for consensus
    #[arg(short = 'm', long = "min-input-base-quality", default_value = "10")]
    pub min_input_base_quality: u8,

    /// Produce per-base tags (ad/bd, ae/be, ac/bc, aq/bq) in addition to per-read tags
    #[arg(short = 'B', long = "output-per-base-tags", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = super::common::parse_bool)]
    pub output_per_base_tags: bool,

    /// Quality-trim reads before consensus calling (removes low-quality bases from ends)
    #[arg(long = "trim", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = super::common::parse_bool)]
    pub trim: bool,

    /// Minimum consensus base quality (output consensus bases below this are masked to N)
    #[arg(long = "min-consensus-base-quality", default_value = "2")]
    pub min_consensus_base_quality: u8,

    // ── Codec-specific tuning. ──
    /// Minimum read pairs per strand to form consensus (same as --min-reads)
    #[arg(short = 'M', long = "min-reads", default_value = "1")]
    pub min_reads: usize,

    /// Maximum read pairs per strand (downsample if exceeded)
    #[arg(long = "max-reads")]
    pub max_reads: Option<usize>,

    /// Minimum duplex overlap length in bases.
    #[arg(short = 'd', long = "min-duplex-length", default_value = "1")]
    pub min_duplex_length: usize,

    /// Reduce single-strand region quality to this value (0-93).
    #[arg(long = "single-strand-qual")]
    pub single_strand_qual: Option<u8>,

    /// Reduce outer bases quality to this value (0-93).
    #[arg(short = 'Q', long = "outer-bases-qual")]
    pub outer_bases_qual: Option<u8>,

    /// Number of outer bases to reduce quality for.
    #[arg(short = 'O', long = "outer-bases-length", default_value = "5")]
    pub outer_bases_length: usize,

    /// Maximum duplex disagreement rate (0.0-1.0).
    #[arg(short = 'x', long = "max-duplex-disagreement-rate", default_value = "1.0")]
    pub max_duplex_disagreement_rate: f64,

    /// Maximum number of duplex disagreements.
    #[arg(short = 'X', long = "max-duplex-disagreements")]
    pub max_duplex_disagreements: Option<usize>,

    // ── Chain-builder slots (populated by Codec::execute, not by clap).
    // Invisible to the CLI; carried so build_codec_chain can access the
    // full command configuration without holding a reference back to the
    // Codec struct. ────────────────────────────────────────────────────
    /// Input / output BAM paths. Populated by `Codec::execute`.
    #[arg(skip)]
    pub io: super::common::BamIoOptions,

    /// Rejects output options. Populated by `Codec::execute`.
    #[arg(skip)]
    pub rejects_opts: super::common::RejectsOptions,

    /// Statistics output options. Populated by `Codec::execute`.
    #[arg(skip)]
    pub stats_opts: super::common::StatsOptions,

    /// Read group / read name prefix options. Populated by `Codec::execute`.
    #[arg(skip)]
    pub read_group: super::common::ReadGroupOptions,
}

impl Default for CodecOptions {
    fn default() -> Self {
        let consensus = super::common::ConsensusCallingOptions::default();
        Self {
            error_rate_pre_umi: consensus.error_rate_pre_umi,
            error_rate_post_umi: consensus.error_rate_post_umi,
            min_input_base_quality: consensus.min_input_base_quality,
            output_per_base_tags: consensus.output_per_base_tags,
            trim: consensus.trim,
            min_consensus_base_quality: consensus.min_consensus_base_quality,
            min_reads: 1,
            max_reads: None,
            min_duplex_length: 1,
            single_strand_qual: None,
            outer_bases_qual: None,
            outer_bases_length: 5,
            max_duplex_disagreement_rate: 1.0,
            max_duplex_disagreements: None,
            io: super::common::BamIoOptions::default(),
            rejects_opts: super::common::RejectsOptions::default(),
            stats_opts: super::common::StatsOptions::default(),
            read_group: super::common::ReadGroupOptions::default(),
        }
    }
}

impl CodecOptions {
    /// Reconstruct the shared [`ConsensusCallingOptions`] from the inlined
    /// flat fields.
    #[must_use]
    pub fn consensus(&self) -> super::common::ConsensusCallingOptions {
        super::common::ConsensusCallingOptions {
            error_rate_pre_umi: self.error_rate_pre_umi,
            error_rate_post_umi: self.error_rate_post_umi,
            min_input_base_quality: self.min_input_base_quality,
            output_per_base_tags: self.output_per_base_tags,
            trim: self.trim,
            min_consensus_base_quality: self.min_consensus_base_quality,
        }
    }

    /// Validate the numeric/semantic codec parameters (everything except I/O
    /// path existence). Shared by the standalone `Codec::validate` and the
    /// `runall` `Stage::Codec` arm so both reject degenerate configs
    /// identically (S5c2-001).
    pub(crate) fn validate(&self) -> Result<()> {
        // Validate error rates.
        if self.error_rate_pre_umi == 0 {
            bail!("error-rate-pre-umi must be > 0");
        }
        if self.error_rate_post_umi == 0 {
            bail!("error-rate-post-umi must be > 0");
        }

        // Validate optional quality ceilings (0-93 Phred range).
        const MAX_PHRED: u8 = 93;
        if let Some(qual) = self.single_strand_qual {
            if qual > MAX_PHRED {
                bail!("single-strand-qual ({qual}) exceeds maximum Phred score ({MAX_PHRED})");
            }
        }
        if let Some(qual) = self.outer_bases_qual {
            if qual > MAX_PHRED {
                bail!("outer-bases-qual ({qual}) exceeds maximum Phred score ({MAX_PHRED})");
            }
        }

        // Validate min/max reads.
        if self.min_reads == 0 {
            bail!("min-reads must be >= 1");
        }
        if let Some(max) = self.max_reads {
            if max < self.min_reads {
                bail!("max-reads ({}) must be >= min-reads ({})", max, self.min_reads);
            }
        }

        // Validate duplex length.
        if self.min_duplex_length == 0 {
            bail!("min-duplex-length must be >= 1");
        }

        // Validate disagreement rate. Reject non-finite values (NaN/inf) first,
        // since NaN comparisons are always false and would otherwise slip past
        // the `0.0..=1.0` range test below.
        if !self.max_duplex_disagreement_rate.is_finite()
            || self.max_duplex_disagreement_rate < 0.0
            || self.max_duplex_disagreement_rate > 1.0
        {
            bail!("max-duplex-disagreement-rate must be between 0.0 and 1.0");
        }

        Ok(())
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "codec",
    about = "\x1b[38;5;180m[CONSENSUS]\x1b[0m      \x1b[36mCall CODEC consensus reads from grouped BAM\x1b[0m",
    long_about = r#"
Calls duplex consensus reads from CODEC sequencing data.

CODEC (Bae et al. 2023) is a protocol in which a single read pair sequences both strands of the
original duplex molecule: R1 comes from one strand and R2 from the opposite strand, so even a
single read pair can yield a duplex consensus. Prior to running this tool, group reads with
`group` using the `adjacency` (or `identity`/`edit`) strategy — NOT `paired`; MI tags must not
carry `/A` or `/B` suffixes (CODEC pairs the two strands by position within each MI group).

For each molecule the tool clips read pairs that extend past their mate, filters R1s and R2s for
compatible alignments, requires sufficient duplex overlap, calls a single-strand consensus from
the R1s and from the R2s separately, then combines the two into a duplex consensus and applies
quality masking. The consensus reads produced are unmapped single fragments (not paired-end) and
should be aligned afterwards.

Input must be template-coordinate sorted (or queryname sorted) and carry `RX` (raw UMI) and `MI`
(molecule ID, from `group`) tags; `CB` (cell barcode) is used when present.

Methylation-aware consensus calling (`--methylation-mode`) is not supported for CODEC; use
`simplex` or `duplex` for EM-seq / TAPs workflows.

Quality-reduction and duplex-disagreement knobs:

* `--single-strand-qual` caps the reported quality of bases supported by only one strand
  (single-strand regions), since those bases lack duplex confirmation.
* `--outer-bases-qual` / `--outer-bases-length` cap the quality of the outermost N bases of each
  read, which are more error-prone (adapter/ligation artifacts).
* `--min-duplex-length` sets the minimum overlap (in bases) required for a region to be called as
  duplex.
* `--max-duplex-disagreements` / `--max-duplex-disagreement-rate` bound how much the two strands
  may disagree before the molecule is rejected.

Reads that fail these thresholds can be captured with `--rejects`.

Consensus reads carry optional per-read and per-base tags. The first letter (a, b, c) marks the
first single-strand consensus (a), the second single-strand consensus (b), or the final duplex
consensus (c); the second letter is the quantity, upper case for per-read values and lower case
for per-base values:

  consensus depth      [aD,bD,cD] (int)  : maximum depth of raw reads
  consensus min depth  [aM,bM,cM] (int)  : minimum depth of raw reads
  consensus error rate [aE,bE,cE] (float): fraction of raw bases disagreeing with the consensus
  consensus depth      [ad,bd] (short[]) : per-base count of contributing reads
  consensus errors     [ae,be] (short[]) : per-base count of disagreeing reads
  consensus bases      [ac,bc] (string)  : single-strand consensus bases
  consensus quals      [aq,bq] (string)  : single-strand consensus qualities

Example:

  fgumi codec -i grouped.bam -o codec_consensus.bam --min-reads 1
"#
)]
pub struct Codec {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Options for writing rejected reads
    #[command(flatten)]
    pub rejects_opts: RejectsOptions,

    /// Options for writing statistics
    #[command(flatten)]
    pub stats_opts: StatsOptions,

    /// Read group and name prefix options
    #[command(flatten)]
    pub read_group: ReadGroupOptions,

    /// Per-stage codec tuning. Flattened so the standalone command
    /// exposes `--error-rate-pre-umi`, `--min-reads`,
    /// `--min-duplex-length`, `--outer-bases-qual`, etc. directly and
    /// `runall` exposes them as `--codec::error-rate-pre-umi`,
    /// `--codec::min-reads`, `--codec::min-duplex-length`,
    /// `--codec::outer-bases-qual`, etc.
    #[command(flatten)]
    pub options: CodecOptions,

    /// Threading options for parallel processing
    #[command(flatten)]
    pub threading: ThreadingOptions,

    /// Compression options for output
    #[command(flatten)]
    pub compression: CompressionOptions,

    /// Scheduler and pipeline statistics options.
    #[command(flatten)]
    pub scheduler_opts: SchedulerOptions,

    /// Queue memory options.
    #[command(flatten)]
    pub queue_memory: QueueMemoryOptions,
}

impl Command for Codec {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs
        self.validate()?;

        // Route through the typed-step chain for every thread count. `--threads`
        // unset resolves to a one-worker chain (`num_threads() == 1`) — output
        // identical to `--threads 1`, with no separate single-threaded path. The
        // chain opens the source once (BAM or SAM) and runs CODEC consensus in
        // `add_codec`. Populate the #[arg(skip)] chain-builder slots on a clone of
        // CodecOptions before handing it to the ChainSpec.
        use crate::pipeline::chains::{ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag};
        let mut options = self.options.clone();
        options.io.clone_from(&self.io);
        options.rejects_opts.clone_from(&self.rejects_opts);
        options.stats_opts.clone_from(&self.stats_opts);
        options.read_group.clone_from(&self.read_group);

        let spec = ChainSpec {
            async_reader: self.io.async_reader,
            stages: vec![Stage::Codec],
            source: SourceSpec::Bam(self.io.input.clone()),
            sink: SinkSpec::Bam(self.io.output.clone()),
            stage_opts: StageOptionsBag { codec: Some(options), ..Default::default() },
            threading: self.threading.clone(),
            compression: self.compression.clone(),
            scheduler: self.scheduler_opts.clone(),
            queue_memory: self.queue_memory.clone(),
            command_line: command_line.to_string(),
        };
        crate::pipeline::chains::build_for(spec)?.run()
    }
}
impl Codec {
    /// Validates command-line arguments
    fn validate(&self) -> Result<()> {
        // Shared validator: input path existence (skipped for stdin),
        // matching the sibling consensus commands (e.g. `Simplex`).
        self.io.validate()?;

        // Numeric/semantic checks shared with the runall Stage::Codec arm.
        self.options.validate()?;

        Ok(())
    }
}

// Helper methods moved to consensus_runner module

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::consensus_runner::ConsensusStatsOps;
    use crate::consensus::codec_caller::CodecConsensusStats;
    use crate::sam::SamTag;
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use rstest::rstest;
    use std::path::PathBuf;

    /// Helper to create a Codec with specified input/output paths
    fn create_codec_with_paths(input: PathBuf, output: PathBuf) -> Codec {
        Codec {
            io: BamIoOptions { input, output, ..Default::default() },
            rejects_opts: RejectsOptions::default(),
            stats_opts: StatsOptions::default(),
            read_group: ReadGroupOptions::default(),
            options: CodecOptions {
                output_per_base_tags: false,
                min_consensus_base_quality: 0,
                min_reads: 1,
                max_reads: None,
                ..CodecOptions::default()
            },
            threading: ThreadingOptions::none(),
            compression: CompressionOptions { compression_level: 1 },
            scheduler_opts: SchedulerOptions::default(),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    fn create_test_codec() -> Codec {
        create_codec_with_paths(PathBuf::from("input.bam"), PathBuf::from("output.bam"))
    }

    #[test]
    fn test_validation() {
        // `validate()` now checks input-path existence via `io.validate()`,
        // so point the command at a real (empty) temp file.
        let tmp = TempDir::new().expect("create temp dir");
        let input_path = tmp.path().join("input.bam");
        std::fs::write(&input_path, b"").expect("write temp input");
        let mut cmd = create_codec_with_paths(input_path, tmp.path().join("output.bam"));

        assert!(cmd.validate().is_ok());

        // Test invalid error rate
        cmd.options.error_rate_pre_umi = 0;
        assert!(cmd.validate().is_err());
        cmd.options.error_rate_pre_umi = 45;

        // Test invalid min reads
        cmd.options.min_reads = 0;
        assert!(cmd.validate().is_err());
        cmd.options.min_reads = 1;

        // Test invalid max < min
        cmd.options.max_reads = Some(0);
        assert!(cmd.validate().is_err());
        cmd.options.max_reads = None;

        // Test invalid duplex length
        cmd.options.min_duplex_length = 0;
        assert!(cmd.validate().is_err());
        cmd.options.min_duplex_length = 1;

        // Test invalid disagreement rate
        cmd.options.max_duplex_disagreement_rate = 1.5;
        assert!(cmd.validate().is_err());
        cmd.options.max_duplex_disagreement_rate = 1.0;

        assert!(cmd.validate().is_ok());
    }

    #[test]
    fn test_stats_to_metrics() {
        let stats = CodecConsensusStats {
            total_input_reads: 1000,
            consensus_reads_generated: 500,
            reads_filtered: 500,
            ..Default::default()
        };

        let metrics = stats.to_metrics();
        assert_eq!(metrics.total_input_reads, 1000);
        assert_eq!(metrics.consensus_reads, 500);
        assert_eq!(metrics.filtered_reads, 500);
    }

    // Integration tests
    use fgumi_raw_bam::{
        SamBuilder as RawSamBuilder, flags, raw_record_to_record_buf, testutil::encode_op,
    };
    use noodles::sam::Header;
    use noodles::sam::alignment::record_buf::RecordBuf;
    use tempfile::TempDir;

    fn to_record_buf(raw: fgumi_raw_bam::RawRecord) -> RecordBuf {
        raw_record_to_record_buf(&raw, &Header::default())
            .expect("raw_record_to_record_buf failed in test")
    }

    /// Helper to create a CODEC-style FR read pair with proper overlap
    /// R1 forward at start1, R2 reverse at start2, with `read_len` bases each
    #[allow(clippy::cast_possible_truncation)]
    fn create_codec_fr_pair_overlapping(
        mi_value: &str,
        start1: usize,
        start2: usize,
        read_len: usize,
        quals: &[u8],
    ) -> (RecordBuf, RecordBuf) {
        // Use a simple reference-matching sequence
        let seq_forward = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seq = &seq_forward[..read_len];

        // For R2 reverse, we need the reverse complement
        let seq_rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect();

        let insert_size: i32 = (start2 + read_len) as i32 - start1 as i32;
        let cigar = encode_op(0, read_len);
        // R1: PAIRED | PROPERLY_SEGMENTED | FIRST_SEGMENT | MATE_REVERSE = 0x1|0x2|0x40|0x20 = 99
        const PROPERLY_PAIRED: u16 = 0x2;
        let r1_flags = flags::PAIRED | PROPERLY_PAIRED | flags::FIRST_SEGMENT | flags::MATE_REVERSE;
        // R2: PAIRED | PROPERLY_SEGMENTED | LAST_SEGMENT | REVERSE = 0x1|0x2|0x80|0x10 = 147
        let r2_flags = flags::PAIRED | PROPERLY_PAIRED | flags::LAST_SEGMENT | flags::REVERSE;

        let mut b1 = RawSamBuilder::new();
        b1.read_name(format!("read_{mi_value}").as_bytes())
            .flags(r1_flags)
            .ref_id(0)
            .pos(start1 as i32 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(seq)
            .qualities(quals)
            .mate_ref_id(0)
            .mate_pos(start2 as i32 - 1)
            .template_length(insert_size);
        b1.add_string_tag(SamTag::MI, mi_value.as_bytes());
        b1.add_string_tag(SamTag::RG, b"A");
        let r1 = to_record_buf(b1.build());

        let mut b2 = RawSamBuilder::new();
        b2.read_name(format!("read_{mi_value}").as_bytes())
            .flags(r2_flags)
            .ref_id(0)
            .pos(start2 as i32 - 1)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq_rc)
            .qualities(quals)
            .mate_ref_id(0)
            .mate_pos(start1 as i32 - 1)
            .template_length(-insert_size);
        b2.add_string_tag(SamTag::MI, mi_value.as_bytes());
        b2.add_string_tag(SamTag::RG, b"A");
        let r2 = to_record_buf(b2.build());

        (r1, r2)
    }

    fn write_codec_bam(path: &std::path::Path, records: Vec<RecordBuf>) -> Result<()> {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::{ReferenceSequence, header::Version};
        use std::num::NonZeroUsize;

        let mut header = Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(Version::new(
                1, 6,
            )))
            .build();

        // Add reference sequence
        let rs = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(1000).expect("non-zero chromosome length"),
        );
        header.reference_sequences_mut().insert(bstr::BString::from("chr1"), rs);

        // Add read group
        let rg = Map::<noodles::sam::header::record::value::map::ReadGroup>::default();
        header.read_groups_mut().insert(bstr::BString::from("A"), rg);

        let mut writer = noodles::bam::io::writer::Builder.build_from_path(path)?;
        writer.write_header(&header)?;

        for record in &records {
            writer.write_alignment_record(&header, record)?;
        }

        Ok(())
    }

    fn read_bam_records(path: &std::path::Path) -> Result<Vec<RecordBuf>> {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(path)?;
        let header = reader.read_header()?;
        let records: Vec<_> = reader.record_bufs(&header).collect::<std::io::Result<Vec<_>>>()?;
        Ok(records)
    }

    #[test]
    fn test_codec_execute_basic() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create FR pairs with same MI (molecule ID)
        // For proper overlap: R1 at pos 100 forward, R2 at pos 105 reverse, both 20bp
        // This gives overlap from pos 105-119 (15 bases overlap)
        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;

        // Execute should complete without errors (coverage is the goal)
        cmd.execute("test")?;

        // Output file should be created
        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_rejects() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let rejects_path = dir.path().join("rejects.bam");

        let mut records = Vec::new();
        // Create multiple FR pairs
        for i in 0..3 {
            let (r1, r2) = create_codec_fr_pair_overlapping(
                &format!("UMI{i:03}"),
                100 + i * 50,
                105 + i * 50,
                20,
                &[30; 20],
            );
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.rejects_opts.rejects = Some(rejects_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;

        cmd.execute("test")?;

        // Rejects file should exist
        assert!(rejects_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_stats() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");
        let stats_path = dir.path().join("stats.txt");

        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.stats_opts.stats = Some(stats_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;

        cmd.execute("test")?;

        // Stats file should exist and contain data
        assert!(stats_path.exists());
        let stats_content = std::fs::read_to_string(&stats_path)?;
        assert!(!stats_content.is_empty());

        Ok(())
    }

    #[test]
    fn test_codec_execute_multithreaded() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_single = dir.path().join("output_single.bam");
        let output_multi = dir.path().join("output_multi.bam");
        let rejects_path = dir.path().join("rejects.bam");

        // Create many UMI groups to test parallel processing (25+ to trigger multiple batches with BATCH_SIZE=10)
        let mut records = Vec::new();
        for i in 0..25 {
            let (r1, r2) = create_codec_fr_pair_overlapping(
                &format!("UMI{i:03}"),
                100 + i * 50,
                105 + i * 50,
                20,
                &[30; 20],
            );
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        // Single-threaded
        let mut cmd_single = create_codec_with_paths(input_path.clone(), output_single.clone());
        cmd_single.read_group.read_name_prefix = Some("codec".to_string());
        cmd_single.options.outer_bases_length = 0;

        cmd_single.execute("test")?;

        // Multi-threaded with rejects
        let mut cmd_multi = create_codec_with_paths(input_path, output_multi.clone());
        cmd_multi.rejects_opts.rejects = Some(rejects_path.clone());
        cmd_multi.read_group.read_name_prefix = Some("codec".to_string());
        cmd_multi.options.outer_bases_length = 0;
        cmd_multi.threading = ThreadingOptions::new(4);

        cmd_multi.execute("test")?;

        let records_single = read_bam_records(&output_single)?;
        let records_multi = read_bam_records(&output_multi)?;

        // Both should produce the same number of consensus reads
        assert_eq!(
            records_single.len(),
            records_multi.len(),
            "Single and multi-threaded should produce same number of reads"
        );

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_per_base_tags() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;
        cmd.options.output_per_base_tags = true; // Enable per-base tags

        cmd.execute("test")?;

        let output_records = read_bam_records(&output_path)?;
        assert!(!output_records.is_empty());

        Ok(())
    }

    #[test]
    fn test_codec_execute_with_trim() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        let mut records = Vec::new();
        // Create reads with low quality at ends
        let quals = [5, 5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5, 5, 5, 5]; // Low quality at ends
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &quals);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;
        cmd.options.trim = true; // Enable trimming

        cmd.execute("test")?;

        // Should complete successfully (may or may not produce output depending on overlap after trim)
        assert!(output_path.exists());

        Ok(())
    }

    // Note: Overlapping consensus tests removed - CODEC does not support overlapping consensus
    // (matching fgbio's CallCodecConsensusReads which has no such option).

    #[test]
    fn test_codec_execute_with_max_reads() -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create multiple read pairs for same UMI (to test max_reads downsampling)
        let mut records = Vec::new();
        for i in 0..5 {
            let (r1, r2) =
                create_codec_fr_pair_overlapping("UMI001", 100 + i, 105 + i, 20, &[30; 20]);
            records.push(r1);
            records.push(r2);
        }

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;
        cmd.options.max_reads = Some(2); // Limit to 2 reads per strand

        cmd.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    #[test]
    fn test_codec_execute_validation_errors() {
        // Test various validation errors
        let mut cmd = create_test_codec();

        // Test post-UMI error rate = 0
        cmd.options.error_rate_post_umi = 0;
        assert!(cmd.validate().is_err());
        cmd.options.error_rate_post_umi = 40;

        // Test max < min reads
        cmd.options.min_reads = 5;
        cmd.options.max_reads = Some(2);
        assert!(cmd.validate().is_err());
        cmd.options.min_reads = 1;
        cmd.options.max_reads = None;

        // Test invalid disagreement rate
        cmd.options.max_duplex_disagreement_rate = -0.1;
        assert!(cmd.validate().is_err());
    }

    /// Parameterized test for all threading modes.
    ///
    /// Tests:
    /// - `None`: Single-threaded fast path, no pipeline
    /// - `Some(1)`: Pipeline with 1 thread
    /// - `Some(2)`: Pipeline with 2 threads
    #[rstest]
    #[case::fast_path(ThreadingOptions::none())]
    #[case::pipeline_1(ThreadingOptions::new(1))]
    #[case::pipeline_2(ThreadingOptions::new(2))]
    fn test_threading_modes(#[case] threading: ThreadingOptions) -> Result<()> {
        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");
        let output_path = dir.path().join("output.bam");

        // Create FR pairs with same MI (molecule ID)
        let mut records = Vec::new();
        let (r1, r2) = create_codec_fr_pair_overlapping("UMI001", 100, 105, 20, &[30; 20]);
        records.push(r1);
        records.push(r2);

        write_codec_bam(&input_path, records)?;

        let mut cmd = create_codec_with_paths(input_path, output_path.clone());
        cmd.read_group.read_name_prefix = Some("codec".to_string());
        cmd.options.outer_bases_length = 0;
        cmd.threading = threading;
        cmd.execute("test")?;

        assert!(output_path.exists());

        Ok(())
    }

    /// Asserts that single-threaded and multi-threaded codec produce the same number of
    /// consensus records and identical CB tag presence when some groups have CB and some do not.
    #[test]
    fn test_threading_parity_mixed_cb() -> Result<()> {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::value::Value;

        let dir = TempDir::new()?;
        let input_path = dir.path().join("input.bam");

        // Group 1: MI001 with CB=CELL1
        let (mut r1, mut r2) = create_codec_fr_pair_overlapping("MI001", 100, 105, 20, &[30; 20]);
        let cb_tag = Tag::new(b'C', b'B');
        r1.data_mut().insert(cb_tag, Value::String("CELL1".into()));
        r2.data_mut().insert(cb_tag, Value::String("CELL1".into()));

        // Group 2: MI002 without CB
        let (r3, r4) = create_codec_fr_pair_overlapping("MI002", 200, 205, 20, &[30; 20]);

        write_codec_bam(&input_path, vec![r1, r2, r3, r4])?;

        // Run single-threaded
        let out_st = dir.path().join("out_st.bam");
        let mut cmd_st = create_codec_with_paths(input_path.clone(), out_st.clone());
        cmd_st.options.outer_bases_length = 0;
        cmd_st.threading = ThreadingOptions::none();
        cmd_st.execute("test")?;
        let records_st = read_bam_records(&out_st)?;

        // Run multi-threaded
        let out_mt = dir.path().join("out_mt.bam");
        let mut cmd_mt = create_codec_with_paths(input_path, out_mt.clone());
        cmd_mt.options.outer_bases_length = 0;
        cmd_mt.threading = ThreadingOptions::new(2);
        cmd_mt.execute("test")?;
        let records_mt = read_bam_records(&out_mt)?;

        assert_eq!(
            records_st.len(),
            records_mt.len(),
            "single-threaded and multi-threaded should produce the same number of consensus records"
        );

        // Both runs should have the same CB tag presence on each record
        let cb_presence_st: Vec<bool> =
            records_st.iter().map(|r| r.data().get(&cb_tag).is_some()).collect();
        let cb_presence_mt: Vec<bool> =
            records_mt.iter().map(|r| r.data().get(&cb_tag).is_some()).collect();
        assert_eq!(
            cb_presence_st, cb_presence_mt,
            "CB tag presence should match between single-threaded and multi-threaded modes"
        );

        Ok(())
    }

    // ── S5c2-001: numeric validation shared by standalone codec and runall ──

    #[test]
    fn codec_options_validate_accepts_defaults() {
        // The Default config (min_reads = 1, valid error rates, etc.) is valid.
        CodecOptions::default().validate().expect("default CodecOptions should validate");
    }

    #[rstest]
    #[case::min_reads_zero(CodecOptions { min_reads: 0, ..Default::default() })]
    #[case::max_lt_min(CodecOptions { min_reads: 3, max_reads: Some(2), ..Default::default() })]
    #[case::error_pre_zero(CodecOptions { error_rate_pre_umi: 0, ..Default::default() })]
    #[case::error_post_zero(CodecOptions { error_rate_post_umi: 0, ..Default::default() })]
    #[case::ss_qual_too_high(CodecOptions { single_strand_qual: Some(94), ..Default::default() })]
    #[case::outer_qual_too_high(CodecOptions { outer_bases_qual: Some(94), ..Default::default() })]
    #[case::min_duplex_zero(CodecOptions { min_duplex_length: 0, ..Default::default() })]
    #[case::disagreement_rate_high(
        CodecOptions { max_duplex_disagreement_rate: 1.5, ..Default::default() }
    )]
    #[case::disagreement_rate_nan(
        CodecOptions { max_duplex_disagreement_rate: f64::NAN, ..Default::default() }
    )]
    #[case::disagreement_rate_negative(
        CodecOptions { max_duplex_disagreement_rate: -0.1, ..Default::default() }
    )]
    #[case::disagreement_rate_inf(
        CodecOptions { max_duplex_disagreement_rate: f64::INFINITY, ..Default::default() }
    )]
    fn codec_options_validate_rejects_degenerate(#[case] options: CodecOptions) {
        assert!(
            options.validate().is_err(),
            "expected degenerate CodecOptions to be rejected: {options:?}"
        );
    }
}
