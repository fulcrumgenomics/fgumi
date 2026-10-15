//! A library entry point for sorting a BAM: [`sort_bam`] builds the same
//! chain-native sort `fgumi sort` runs, from a plain [`SortConfig`].

use std::path::PathBuf;

use anyhow::{Result, bail};

use crate::commands::common::{CompressionOptions, MemoryLimit, ThreadingOptions};
use crate::commands::sort::{SortOptions, SortOrderArg};
use crate::pipeline::chains::staged_output::StagedSortOutput;
use crate::pipeline::chains::{ChainSpec, SourceSpec, Stage, StageOptionsBag, build_for};

pub use fgumi_sort::{KeyTypesSpec, SortStats};

/// What [`sort_bam`] sorts, and how.
///
/// Template-coordinate keys the `CB` cell tag, as `fgumi sort` does; whether
/// its lane is kept is `key_types`'.
#[derive(Debug, Clone)]
pub struct SortConfig {
    /// The input BAM (or SAM).
    pub input: PathBuf,
    /// The output BAM.
    pub output: PathBuf,
    /// The sort order.
    pub order: SortOrderArg,
    /// Worker threads (`fgumi sort --threads`).
    pub threads: usize,
    /// The in-memory budget in bytes, for the whole run (not per thread);
    /// `None` is `fgumi sort`'s default (768 MiB per thread).
    pub max_memory: Option<u64>,
    /// Directories for spill files; empty is the system temp directory.
    pub tmp_dirs: Vec<PathBuf>,
    /// The output's BGZF compression level.
    pub compression: u32,
    /// Write a BAI beside the output (coordinate order only).
    pub write_index: bool,
    /// The template-coordinate key lanes to keep (`--key-types`).
    pub key_types: KeyTypesSpec,
}

/// Sort a BAM with the same chain `fgumi sort` builds.
///
/// The output (and its BAI) is staged beside the destination and published
/// only after a successful run that wrote every record it read, so a failed
/// sort leaves the destination as it was.
///
/// # Errors
///
/// An invalid configuration (`write_index` with a non-coordinate order), any
/// failure of the sort, or a record-count mismatch (nothing is published).
pub fn sort_bam(config: &SortConfig) -> Result<SortStats> {
    if config.write_index && !matches!(config.order, SortOrderArg::Coordinate) {
        bail!("write_index is only valid for coordinate sort");
    }
    let mut sort = SortOptions {
        order: config.order,
        key_types: Some(config.key_types),
        tmp_dirs: config.tmp_dirs.clone(),
        ..SortOptions::default()
    };
    if let Some(bytes) = config.max_memory {
        sort.max_memory = MemoryLimit::Fixed(usize::try_from(bytes).unwrap_or(usize::MAX));
        sort.memory_per_thread = false;
    }
    let (scheduler, queue_memory) = crate::commands::common::sort_chain_runtime(
        sort.max_memory,
        sort.memory_reserve,
        sort.memory_per_thread,
    );
    let staged = StagedSortOutput::create(&config.output, config.write_index, "sort")?;
    let spec = ChainSpec {
        stages: vec![Stage::Sort],
        source: SourceSpec::Bam(config.input.clone()),
        sink: staged.sink(),
        stage_opts: StageOptionsBag { sort: Some(sort), ..Default::default() },
        threading: ThreadingOptions { threads: Some(config.threads) },
        compression: CompressionOptions { compression_level: config.compression },
        scheduler,
        queue_memory,
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Auto,
        verify_crc: true,
        command_line: "fgumi_lib::sort::sort_bam".to_string(),
    };
    let built = build_for(spec)?;
    let slot = built.sort_stats_slot.clone();
    built.run()?;
    let stats = slot.and_then(|s| *s.lock()).unwrap_or_default();
    staged.finish(&stats)?;
    Ok(stats)
}
