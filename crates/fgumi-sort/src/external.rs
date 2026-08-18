//! Raw-bytes sorting implementation for BAM files.
//!
//! Items in this module are reachable only via the crate-root re-exports
//! ([`crate::RawExternalSorter`] etc.) and via siblings inside `fgumi-sort`.
//! Some helpers (legacy spill-chunk writer, alternative readers) are retained
//! for benchmarking and as fall-back paths.

#![allow(dead_code)]
// Pre-existing module body follows.
//! Instead of fully decoding each BAM record into `RecordBuf`, it uses noodles'
//! lazy `Record` type that stores raw bytes and only parses fields on-demand.
//!
//! # Performance Benefits
//!
//! - **3-4x lower memory usage**: Raw bytes are ~200-400 bytes vs ~800-1200 bytes decoded
//! - **No re-encoding overhead**: Records are written back as raw bytes
//! - **Lazy field access**: Only sort-key fields are parsed
//!
//! # Algorithm
//!
//! 1. Read BAM records as lazy `Record` objects (raw bytes)
//! 2. Extract only the fields needed for sort keys (tid, pos, flags, name, MI)
//! 3. Sort by keys while keeping raw records
//! 4. Write raw record bytes to output

use crate::inline::{
    CbKey32, InMemoryChunk, ProbeableBuffer, RecordBuffer, TEMPLATE_HEADER_SIZE, TemplateKey,
    TemplateKey24, TemplateKey40, TemplateLaneKey, TemplateRecordBuffer, TemplateRecordRef,
    TertKey32,
};
use crate::keys::{QuerynameComparator, RawSortKey, SortOrder};
use crate::memory_probe::{
    BufferProbeStats, ConsumerProbeStats, MergeProbe, SpillProbe, force_mi_collect, log_snapshot,
};
use crate::pooled_chunk_writer::PooledChunkWriter;
use crate::read_ahead::{PooledInputStream, RawReadAheadReader, RecordSource};
use crate::tmp_dir_alloc::TmpDirAllocator;
use crate::worker_pool::SortWorkerPool;
use anyhow::{Context as _, Result};
use crossbeam_channel::{Receiver, Sender, bounded};
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_bam_io::{is_stdin_path, is_stdout_path};
use fgumi_raw_bam::SamTag;
use log::{debug, info, warn};
use noodles::sam::Header;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles_bgzf::io::{
    MultithreadedWriter, Reader as BgzfReader, Writer as BgzfWriter, multithreaded_writer,
    writer::CompressionLevel,
};
use std::collections::{HashMap, HashSet};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};
use tempfile::TempDir;

/// Maximum number of records held in one in-memory chunk.
///
/// Each key carries its ingest position within the chunk (see
/// [`RawSortKey::set_position`](crate::keys::RawSortKey::set_position)) so that
/// the key is a *total* order and the chunk sort can be unstable. That position
/// is a `u32`, so a chunk that grew past `u32::MAX` records would have to stamp
/// two records with the same position — exact name/flag ties would stop being
/// totally ordered and the unstable sort could reorder them, breaking output
/// identity with `samtools sort -n`.
///
/// A chunk spills on `memory_limit` long before this in any real run (holding
/// this many records needs hundreds of GB). This is the backstop that makes the
/// `usize` -> `u32` narrowing at the stamp site provably lossless rather than
/// merely unlikely to overflow.
const MAX_CHUNK_RECORDS: usize = u32::MAX as usize;

/// Buffer wrapped around sort's input when it is read synchronously.
///
/// The worker pool's block reader pulls in large gulps, so the input is
/// buffered well above the 8 KiB `io::Stdin` and `File` default. stdin always
/// gets this buffer; a regular file gets it only on the synchronous path,
/// because `--async-reader` hands the file to `PrefetchReader`, which already
/// reads ahead into its own bounded queue.
const SORT_INPUT_BUFFER_SIZE: usize = 2 * 1024 * 1024;

/// Whether the current chunk must be spilled before accepting another record.
///
/// Spilling is normally driven by `memory_limit`; `MAX_CHUNK_RECORDS` is the
/// second, independent trigger that keeps ingest positions unique within a
/// chunk (see [`MAX_CHUNK_RECORDS`]).
fn should_spill_chunk(memory_used: usize, memory_limit: usize, chunk_records: usize) -> bool {
    memory_used >= memory_limit || chunk_records >= MAX_CHUNK_RECORDS
}

/// Render the effective per-phase worker counts for a log line.
///
/// Collapses to a single number when the phases agree (the common case, where
/// neither `--sort-threads` nor `--merge-threads` was overridden) and spells
/// both out otherwise. Shared by the engine's end-of-run summary and the `sort`
/// command's config dump so the two emitters cannot drift.
#[must_use]
pub fn format_thread_counts(sort_threads: usize, merge_threads: usize) -> String {
    if sort_threads == merge_threads {
        format!("{sort_threads}")
    } else {
        format!("sort {sort_threads}, merge {merge_threads}")
    }
}

// ============================================================================
// Per-Phase Timing for Sort Pipeline
// ============================================================================

/// Tracks wall-clock time spent in each phase of the sort pipeline.
///
/// Used to identify bottlenecks and validate thread architecture changes.
/// All times are cumulative (multiple spill cycles accumulate).
#[derive(Debug, Default)]
struct SortPhaseTimer {
    /// Time reading records from input BAM (includes BGZF decompression).
    read_secs: f64,
    /// Time sorting in-memory buffers (rayon parallel sort or single-threaded).
    sort_secs: f64,
    /// Time writing sorted chunks to temp files (BGZF compression).
    spill_write_secs: f64,
    /// Time consolidating temp files when limit exceeded.
    consolidate_secs: f64,
    /// Time in the final k-way merge phase (includes reader decompression + writer compression).
    merge_secs: f64,
    /// Time writing in-memory-only output (no merge needed).
    write_output_secs: f64,
    /// Number of spill cycles (sort + write).
    spill_count: usize,
    /// Number of consolidation merges.
    consolidate_count: usize,
    /// Total bytes written to spill files.
    total_spill_bytes: u64,
    /// Wall-clock start of the entire sort operation.
    overall_start: Option<Instant>,
    /// Tracks the start of the current read span (between spills).
    read_span_start: Option<Instant>,
}

impl SortPhaseTimer {
    fn new() -> Self {
        Self {
            overall_start: Some(Instant::now()),
            read_span_start: Some(Instant::now()),
            ..Default::default()
        }
    }

    /// End a read span (call before sort/spill). Returns elapsed read time.
    fn end_read_span(&mut self) -> Duration {
        if let Some(start) = self.read_span_start.take() {
            let elapsed = start.elapsed();
            self.read_secs += elapsed.as_secs_f64();
            elapsed
        } else {
            Duration::ZERO
        }
    }

    /// Start a new read span (call after spill write completes).
    fn begin_read_span(&mut self) {
        self.read_span_start = Some(Instant::now());
    }

    /// Time a closure and accumulate elapsed seconds into `field`.
    fn time<T>(field: &mut f64, f: impl FnOnce() -> T) -> T {
        let start = Instant::now();
        let result = f();
        *field += start.elapsed().as_secs_f64();
        result
    }

    /// Time a sort operation.
    fn time_sort<T>(&mut self, f: impl FnOnce() -> T) -> T {
        Self::time(&mut self.sort_secs, f)
    }

    /// Time a spill write operation.
    fn time_spill_write<T>(&mut self, f: impl FnOnce() -> Result<T>) -> Result<T> {
        let result = Self::time(&mut self.spill_write_secs, f);
        self.spill_count += 1;
        result
    }

    /// How many chunks have been spilled, counted where the spill is timed.
    fn spill_count(&self) -> usize {
        self.spill_count
    }

    /// Record the size of a spill file.
    fn record_spill_size(&mut self, path: &Path) {
        if let Ok(meta) = std::fs::metadata(path) {
            self.total_spill_bytes += meta.len();
        }
    }

    /// Add the bytes a spill run grew by, given its size before the chunk was
    /// written.
    ///
    /// Run formation appends several chunks to one file, so re-adding the file's
    /// whole length after each chunk would sum `C + 2C + ... + NC` and report a
    /// spill volume quadratic in the chunk count. Spill volume is the figure used
    /// to judge whether run formation reduced spill I/O at all, so it has to be
    /// the delta.
    fn record_spill_growth(&mut self, path: &Path, size_before: u64) {
        if let Ok(meta) = std::fs::metadata(path) {
            self.total_spill_bytes += meta.len().saturating_sub(size_before);
        }
    }

    /// Time a consolidation operation.
    fn time_consolidate(&mut self, f: impl FnOnce() -> Result<()>) -> Result<()> {
        let start = Instant::now();
        let result = f();
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed > 0.001 {
            // Only count if consolidation actually happened
            self.consolidate_secs += elapsed;
            self.consolidate_count += 1;
        }
        result
    }

    /// Time the merge phase.
    fn time_merge<T>(&mut self, f: impl FnOnce() -> Result<T>) -> Result<T> {
        Self::time(&mut self.merge_secs, f)
    }

    /// Time writing in-memory-only output (no merge needed).
    fn time_write_output<T>(&mut self, f: impl FnOnce() -> Result<T>) -> Result<T> {
        Self::time(&mut self.write_output_secs, f)
    }

    /// Log the final summary.
    ///
    /// `sort_threads` and `merge_threads` are the *effective* per-phase worker
    /// counts ([`RawExternalSorter::phase1_threads`] /
    /// [`RawExternalSorter::phase2_threads`]), not the configured `threads`
    /// they default from — the summary reports what the phases ran with.
    ///
    /// `max_temp_files` is reported so a run that consolidated says which limit
    /// it consolidated against.
    #[allow(clippy::cast_precision_loss)]
    fn log_summary(&self, sort_threads: usize, merge_threads: usize, max_temp_files: usize) {
        let overall = self.overall_start.map_or(0.0, |s| s.elapsed().as_secs_f64());
        // Guard against division by zero when sort completes in negligible time.
        let overall_nonzero = if overall > 0.0 { overall } else { f64::EPSILON };
        let read_pct = 100.0 * self.read_secs / overall_nonzero;
        let sort_pct = 100.0 * self.sort_secs / overall_nonzero;
        let spill_pct = 100.0 * self.spill_write_secs / overall_nonzero;
        let spill_count = self.spill_count;
        let read_secs = self.read_secs;
        let sort_secs = self.sort_secs;
        let spill_secs = self.spill_write_secs;

        info!("=== Sort Phase Timing ===");
        info!("  Read + decompress: {read_secs:.1}s ({read_pct:.0}%)");
        info!("  In-memory sort:    {sort_secs:.1}s ({sort_pct:.0}%) [{spill_count} spills]");
        let spill_mb = self.total_spill_bytes as f64 / (1024.0 * 1024.0);
        info!(
            "  Spill write:       {spill_secs:.1}s ({spill_pct:.0}%) [{spill_count} writes, {spill_mb:.1} MB total]"
        );
        if self.consolidate_count > 0 {
            let cons_secs = self.consolidate_secs;
            let cons_pct = 100.0 * cons_secs / overall_nonzero;
            let cons_count = self.consolidate_count;
            // Report the limit that triggered it: consolidation rewrites data
            // that is already sorted, so a run that consolidated wants to know
            // what it was up against. Reporting the number and not what to do
            // about it keeps the recommendation with the caller that chose the
            // limit -- the engine cannot know whether it was requested.
            info!(
                "  Consolidation:     {cons_secs:.1}s ({cons_pct:.0}%) [{cons_count} merges, limit {max_temp_files}]"
            );
        }
        if self.merge_secs > 0.0 {
            let merge_secs = self.merge_secs;
            let merge_pct = 100.0 * merge_secs / overall_nonzero;
            info!("  K-way merge:       {merge_secs:.1}s ({merge_pct:.0}%)");
        }
        if self.write_output_secs > 0.0 {
            let write_secs = self.write_output_secs;
            let write_pct = 100.0 * write_secs / overall_nonzero;
            info!("  Write output:      {write_secs:.1}s ({write_pct:.0}%)");
        }
        info!("  Total wall clock:  {overall:.1}s");
        info!("  Threads: {}", format_thread_counts(sort_threads, merge_threads));
        info!("=========================");
    }
}

// ============================================================================
// Library Lookup for Template-Coordinate Sort
// ============================================================================

/// Deterministic hasher for cell barcode hashing in template-coordinate sort.
///
/// Uses arbitrary fixed seeds so that hash values are reproducible across runs.
#[must_use]
pub fn cb_hasher() -> ahash::RandomState {
    // Arbitrary fixed seeds — chosen for uniqueness, not cryptographic strength.
    ahash::RandomState::with_seeds(
        0xa1b2_c3d4_e5f6_0718,
        0x9182_7364_5546_3728,
        0xfede_dcba_0987_6543,
        0x0011_2233_4455_6677,
    )
}

/// Where a merge writes its output.
///
/// For a regular file, the merge writes to an exclusive, uniquely-named
/// [`tempfile::NamedTempFile`] in the output's directory and atomically
/// [`persist`](MergeOutputTarget::persist)s it into place on success. Because the
/// temp is RAII-owned, *any* error path (a sort-order violation, a write failure,
/// `finish` failing) drops it and removes the file — a rejected or failed merge
/// never leaves a partial/corrupt output or a stray temp behind. The
/// same-directory temp keeps the rename atomic (same filesystem), and exclusive
/// creation avoids clobbering a concurrent run or a stale leftover.
///
/// Stdout can't be renamed, so it is written directly (a mid-merge failure there
/// leaves whatever was already streamed, which is unavoidable for a pipe).
///
/// A `NamedTempFile` is created `0600`, so persisting it verbatim would silently
/// make merged output owner-only. To preserve the semantics of the pre-atomic-temp
/// direct write (`File::create`), the `Temp` variant carries the mode the final
/// file must end up with (see [`target_file_mode`]) and applies it before the
/// rename.
enum MergeOutputTarget {
    /// A regular-file output, staged in a same-directory temp. `dest` is the
    /// resolved final path the temp is atomically renamed onto (a symlinked
    /// output is followed to its real target, so the temp is staged next to —
    /// and renamed onto — the linked file rather than the link itself). `mode`
    /// is the Unix mode to stamp onto the temp before persisting (`None` on
    /// non-Unix, where file modes are managed by the platform).
    Temp {
        temp: tempfile::NamedTempFile,
        dest: PathBuf,
        mode: Option<u32>,
    },
    Stdout(PathBuf),
}

impl MergeOutputTarget {
    fn create(output: &Path) -> Result<Self> {
        if is_stdout_path(output) {
            return Ok(Self::Stdout(output.to_path_buf()));
        }
        // Follow a symlinked destination to the real file it points at, so the
        // atomic rename updates the linked file in place (matching `File::create`,
        // which follows symlinks) instead of replacing the symlink with a regular
        // file. The temp is then staged next to — and renamed onto — that real
        // target, keeping the rename same-directory (hence same-filesystem/atomic).
        let dest = resolve_symlink_output(output)?;
        // The temp+rename `persist` replaces `dest` with a regular file, so an
        // existing FIFO, device, socket, or directory would be silently clobbered
        // (a plain `File::create` would instead write *through* a FIFO/device).
        // Only a missing path or an existing regular file can be staged this way;
        // reject anything else with a clear error rather than corrupt it.
        if let Ok(metadata) = std::fs::metadata(&dest) {
            anyhow::ensure!(
                metadata.file_type().is_file(),
                "merge output must be a regular file or stdout: {}",
                dest.display()
            );
        }
        let dir = dest
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
        let temp = tempfile::Builder::new()
            .prefix(".fgumi-merge-")
            .suffix(".tmp")
            .tempfile_in(&dir)
            .with_context(|| format!("failed to create a merge temp file in {}", dir.display()))?;
        // Resolve the final mode now (before any BGZF writer threads spawn) so the
        // umask read below is single-threaded and cannot race a concurrent create.
        let mode = target_file_mode(&dest);
        Ok(Self::Temp { temp, dest, mode })
    }

    /// The path the merged BAM is written to.
    fn path(&self) -> &Path {
        match self {
            Self::Temp { temp, .. } => temp.path(),
            Self::Stdout(path) => path,
        }
    }

    /// Atomically moves the finished temp onto its resolved destination (a no-op
    /// for stdout), first stamping the resolved mode so the output is not left
    /// temp-private (`0600`). Consumes `self`, disarming the RAII auto-remove.
    fn persist(self) -> Result<()> {
        if let Self::Temp { temp, dest, mode } = self {
            if let Some(mode) = mode {
                #[cfg(unix)]
                {
                    use std::os::unix::fs::PermissionsExt;
                    temp.as_file()
                        .set_permissions(std::fs::Permissions::from_mode(mode))
                        .with_context(|| {
                            format!("failed to set mode on merge temp for {}", dest.display())
                        })?;
                }
                #[cfg(not(unix))]
                let _ = mode;
            }
            temp.persist(&dest).map_err(|e| e.error).with_context(|| {
                format!("failed to finalize merged output at {}", dest.display())
            })?;
        }
        Ok(())
    }
}

/// Resolves a symlinked output path to the real file the atomic rename must
/// target, so a `latest.bam -> run.bam` symlink is followed and updated in
/// place (matching `File::create`, which follows symlinks) instead of being
/// replaced by a regular file. Non-symlink paths — including ones that do not
/// exist yet — are returned unchanged. Bails on a symlink cycle (or an
/// unreasonably long chain) rather than looping forever.
fn resolve_symlink_output(output: &Path) -> Result<PathBuf> {
    // POSIX ELOOP triggers around 40 levels; cap here at the same bound so a
    // cyclic/pathological chain fails fast instead of spinning.
    const MAX_LINKS: usize = 40;
    let mut path = output.to_path_buf();
    for _ in 0..MAX_LINKS {
        // `is_symlink` returns false for a nonexistent path, so a brand-new
        // output (or its final component) short-circuits here unchanged.
        if !path.is_symlink() {
            return Ok(path);
        }
        let target = std::fs::read_link(&path)
            .with_context(|| format!("failed to read symlink output {}", path.display()))?;
        // Relative link targets resolve against the link's own directory.
        path = if target.is_absolute() {
            target
        } else {
            path.parent().map_or_else(|| target.clone(), |parent| parent.join(&target))
        };
    }
    anyhow::bail!("too many levels of symbolic links while resolving output {}", output.display())
}

/// The mode the merged output file must carry, matching the pre-atomic-temp write
/// path (`File::create`): an existing destination keeps its current mode; a new
/// file gets `0o666 & !umask`. Returns `None` on non-Unix, where the temp's
/// platform-managed permissions are left as-is.
#[cfg(unix)]
// The `not(unix)` sibling returns `None`, so the `Option` is load-bearing there.
#[allow(clippy::unnecessary_wraps, reason = "non-unix sibling returns None")]
fn target_file_mode(output: &Path) -> Option<u32> {
    use std::os::unix::fs::PermissionsExt;
    let mode = match std::fs::metadata(output) {
        // Overwriting: `File::create` never changes an existing file's mode, so keep it.
        Ok(meta) => meta.permissions().mode() & 0o777,
        // New file: `File::create` opens `0o666`, then the kernel masks it with umask.
        Err(_) => 0o666 & !process_umask(),
    };
    Some(mode)
}

#[cfg(not(unix))]
fn target_file_mode(_output: &Path) -> Option<u32> {
    None
}

/// Reads the process file-creation mask (`umask`) without leaving it changed.
///
/// `umask(2)` can only *set* the mask (returning the previous value), so reading it
/// means setting it to `0` and immediately restoring it. That read-modify-restore is
/// not atomic: two concurrent probes can interleave so the second reads the first's
/// transient `0` and restores `0`, permanently clearing the process mask (and mis-
/// computing every subsequent new-file mode). A process-global lock serializes the
/// probes so each is atomic with respect to every other — needed because
/// `fgumi_lib` is a library and callers may run merges concurrently in one process.
#[cfg(unix)]
fn process_umask() -> u32 {
    use std::sync::Mutex;

    // Serializes the non-atomic umask read-restore against other concurrent probes.
    static UMASK_LOCK: Mutex<()> = Mutex::new(());
    let _guard = UMASK_LOCK.lock().unwrap_or_else(std::sync::PoisonError::into_inner);

    // SAFETY: `umask` has no preconditions and cannot fail; it is `unsafe` only
    // because it is a raw libc binding. We restore the original mask on the very
    // next call under `UMASK_LOCK`, so the process-wide umask is unchanged on return.
    #[allow(unsafe_code)]
    let previous = unsafe {
        let previous = libc::umask(0);
        libc::umask(previous);
        previous
    };
    // `libc::mode_t` is `u16` on macOS (the conversion widens) and `u32` on Linux
    // (the conversion is a no-op), so `useless_conversion` fires only on Linux.
    #[allow(
        clippy::useless_conversion,
        reason = "libc::mode_t is u16 on macOS and u32 on Linux; the From keeps this portable"
    )]
    u32::from(previous)
}

/// Maps read group ID -> library ordinal for O(1) comparison.
///
/// Pre-computes ordinals by sorting library names alphabetically.
/// Empty/unknown library sorts first (ordinal 0).
pub struct LibraryLookup {
    /// RG ID -> library ordinal
    rg_to_ordinal: HashMap<Vec<u8>, u32>,
    /// Deterministic hasher for read name hashing, constructed once for reuse.
    hasher: ahash::RandomState,
}

impl LibraryLookup {
    /// Build lookup from BAM header.
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn from_header(header: &Header) -> Self {
        // Collect all unique library names from read groups
        let mut libraries: Vec<String> = header
            .read_groups()
            .iter()
            .filter_map(|(_, rg)| {
                rg.other_fields().get(&rg_tag::LIBRARY).map(std::string::ToString::to_string)
            })
            .collect();

        // Sort alphabetically and deduplicate
        libraries.sort();
        libraries.dedup();

        // Build library name -> ordinal mapping
        // Empty string gets ordinal 0, then libraries in sorted order
        let mut lib_to_ordinal: HashMap<String, u32> = HashMap::new();
        lib_to_ordinal.insert(String::new(), 0);
        for (i, lib) in libraries.iter().enumerate() {
            lib_to_ordinal.insert(lib.clone(), (i + 1) as u32);
        }

        // Build RG ID -> ordinal mapping
        let rg_to_ordinal: HashMap<Vec<u8>, u32> = header
            .read_groups()
            .iter()
            .map(|(id, rg)| {
                let lib = rg
                    .other_fields()
                    .get(&rg_tag::LIBRARY)
                    .map(std::string::ToString::to_string)
                    .unwrap_or_default();
                let ordinal = *lib_to_ordinal.get(&lib).unwrap_or(&0);
                (id.to_vec(), ordinal)
            })
            .collect();

        // Arbitrary fixed seeds — chosen for uniqueness, not cryptographic strength.
        let hasher = ahash::RandomState::with_seeds(
            0x517c_c1b7_2722_0a95,
            0x1234_5678_90ab_cdef,
            0xfedc_ba98_7654_3210,
            0x0123_4567_89ab_cdef,
        );

        Self { rg_to_ordinal, hasher }
    }

    /// Hash a read name deterministically.
    #[inline]
    #[must_use]
    pub fn hash_name(&self, name: &[u8]) -> u64 {
        self.hasher.hash_one(name)
    }

    /// Get library ordinal for a record (from RG tag in aux data).
    ///
    /// Only used in tests; production code uses [`ordinal_from_rg`](Self::ordinal_from_rg)
    /// with pre-extracted RG bytes from the single-pass aux scan.
    #[cfg(test)]
    #[must_use]
    pub fn get_ordinal(&self, bam: &[u8]) -> u32 {
        fgumi_raw_bam::RawRecordView::new(bam)
            .tags()
            .find_string(SamTag::RG)
            .and_then(|rg| self.rg_to_ordinal.get(rg))
            .copied()
            .unwrap_or(0)
    }

    /// Get library ordinal from pre-extracted RG tag bytes.
    #[inline]
    #[must_use]
    pub fn ordinal_from_rg(&self, rg: Option<&[u8]>) -> u32 {
        rg.and_then(|rg| self.rg_to_ordinal.get(rg)).copied().unwrap_or(0)
    }

    /// Number of distinct library ordinals realized by the header's read groups.
    ///
    /// Used to decide whether the template-key `tertiary` library lane can vary:
    /// `<= 1` means every header read group maps to the same library ordinal, so
    /// the lane is provisionally constant (decode-time verify still backstops the
    /// case where a record carries an RG absent from the header → ordinal 0).
    /// Returns 0 when the header declares no read groups.
    #[must_use]
    pub fn distinct_header_ordinals(&self) -> usize {
        self.rg_to_ordinal.values().copied().collect::<HashSet<u32>>().len()
    }
}

/// Number of records to prefetch per chunk during merge.
/// Larger buffer reduces I/O latency impact during merge.
const MERGE_PREFETCH_SIZE: usize = 1024;

// A bare `RawExternalSorter` keeps the fixed, portable `FALLBACK_MAX_TEMP_FILES`
// as its consolidation limit. Sizing that limit to the process's soft
// `RLIMIT_NOFILE` is a policy decision left to the caller, which opts in through
// [`crate::fd_limit`] — `fgumi sort` does, `fgumi merge` and `fgumi simulate` do not.

/// Working estimate of the raw BAM bytes per template-coordinate record (excluding the
/// inline header and the `TemplateRecordRef<K>` index entry).  Used by the capacity
/// estimator in `sort_template_coordinate_impl` to split `effective_initial_capacity`
/// between the data arena and the ref index.
const EST_BAM_BYTES_PER_TEMPLATE_RECORD: usize = 250;

/// Counting semaphore for limiting concurrent chunk reader I/O.
/// Pre-filled with N tokens; readers acquire before decompressing, release after.
pub(crate) type ChunkReaderSemaphore = (Sender<()>, Receiver<()>);

/// Create a counting semaphore that allows `threads` concurrent readers.
pub(crate) fn make_reader_semaphore(threads: usize) -> Arc<ChunkReaderSemaphore> {
    let limit = threads.max(1);
    let (tx, rx) = bounded(limit);
    for _ in 0..limit {
        tx.send(()).expect("semaphore channel must not be disconnected during initialization");
    }
    Arc::new((tx, rx))
}

// ============================================================================
// Generic Keyed Temp File I/O (works with any RawSortKey)
// ============================================================================
//
// Stores pre-computed sort keys alongside each record for O(1) merge comparisons.
// Format: [key: serialized][len: 4 bytes][record: len bytes] per record

use std::marker::PhantomData;

/// Wrapper for temp chunk writers supporting both raw and compressed output.
enum ChunkWriterInner {
    /// Uncompressed raw output (fastest).
    Raw(BufWriter<std::fs::File>),
    /// Single-threaded BGZF-compressed output.
    SingleThreaded(BgzfWriter<BufWriter<std::fs::File>>),
    /// Multi-threaded BGZF-compressed output (faster for large chunks).
    MultiThreaded(MultithreadedWriter<BufWriter<std::fs::File>>),
}

impl Write for ChunkWriterInner {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            ChunkWriterInner::Raw(w) => w.write(buf),
            ChunkWriterInner::SingleThreaded(w) => w.write(buf),
            ChunkWriterInner::MultiThreaded(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            ChunkWriterInner::Raw(w) => w.flush(),
            ChunkWriterInner::SingleThreaded(w) => w.flush(),
            ChunkWriterInner::MultiThreaded(w) => w.flush(),
        }
    }
}

impl ChunkWriterInner {
    fn finish(self) -> Result<()> {
        match self {
            ChunkWriterInner::Raw(mut w) => {
                w.flush()?;
                Ok(())
            }
            ChunkWriterInner::SingleThreaded(w) => {
                w.finish()?;
                Ok(())
            }
            ChunkWriterInner::MultiThreaded(mut w) => {
                w.finish()?;
                Ok(())
            }
        }
    }
}

/// Generic writer for keyed temp chunks with pre-computed sort keys.
///
/// Works with any type implementing `RawSortKey`.
/// Supports optional BGZF compression for reduced disk usage.
pub(crate) struct GenericKeyedChunkWriter<K: RawSortKey> {
    writer: ChunkWriterInner,
    _marker: PhantomData<K>,
}

impl<K: RawSortKey> GenericKeyedChunkWriter<K> {
    /// Create a new keyed chunk writer with optional compression.
    ///
    /// - `compression_level` 0 = uncompressed (fastest, uses most disk).
    /// - `compression_level` > 0 = BGZF compression at specified level.
    /// - `threads` > 1 enables multi-threaded compression.
    ///
    /// # Errors
    ///
    /// Returns an error if the output file cannot be created.
    ///
    /// # Panics
    ///
    /// Panics if `threads` is greater than 1 but `NonZero::new` receives zero.
    pub fn create(path: &Path, compression_level: u32, threads: usize) -> Result<Self> {
        let file = std::fs::File::create(path)?;
        let buf = BufWriter::with_capacity(256 * 1024, file);

        let writer = if compression_level == 0 {
            ChunkWriterInner::Raw(buf)
        } else if threads > 1 {
            // Use multi-threaded BGZF for faster compression
            let worker_count = NonZero::new(threads).expect("threads > 1");
            let mut builder =
                multithreaded_writer::Builder::default().set_worker_count(worker_count);
            #[allow(clippy::cast_possible_truncation)]
            if let Some(level) = CompressionLevel::new(compression_level as u8) {
                builder = builder.set_compression_level(level);
            }
            ChunkWriterInner::MultiThreaded(builder.build_from_writer(buf))
        } else {
            // Single-threaded BGZF with specified compression level
            #[allow(clippy::cast_possible_truncation)]
            let level = CompressionLevel::new(compression_level as u8).unwrap_or_else(|| {
                CompressionLevel::new(6).expect("compression level 6 is always valid")
            });
            let writer = noodles_bgzf::io::writer::Builder::default()
                .set_compression_level(level)
                .build_from_writer(buf);
            ChunkWriterInner::SingleThreaded(writer)
        };

        Ok(Self { writer, _marker: PhantomData })
    }

    /// Write a keyed record.
    ///
    /// When `K::EMBEDDED_IN_RECORD` is true, the key is embedded in the record
    /// bytes so only the record is written. Otherwise, the key prefix is written
    /// followed by the record.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the underlying writer fails.
    #[inline]
    #[allow(clippy::cast_possible_truncation)]
    pub fn write_record(&mut self, key: &K, record: &[u8]) -> Result<()> {
        if !K::EMBEDDED_IN_RECORD {
            key.write_to(&mut self.writer)?;
        }
        self.writer.write_all(&(record.len() as u32).to_le_bytes())?;
        self.writer.write_all(record)?;
        Ok(())
    }

    /// Finish writing and flush.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing the writer fails.
    pub fn finish(self) -> Result<()> {
        self.writer.finish()
    }
}

/// Result of reading a keyed record from a chunk: `Ok(Some(...))` for a record,
/// `Ok(None)` for EOF, or `Err(...)` for an I/O error.
type ChunkReadResult<K> = Result<Option<(K, Vec<u8>)>>;

/// Read exactly `buf.len()` bytes, distinguishing clean EOF from truncation.
///
/// Returns `Ok(true)` when `buf` is fully filled, `Ok(false)` when zero bytes are
/// available (clean EOF), and `Err` for any partial read or I/O error.
pub(crate) fn read_exact_or_eof<R: Read>(reader: &mut R, buf: &mut [u8]) -> std::io::Result<bool> {
    let mut offset = 0;
    while offset < buf.len() {
        match reader.read(&mut buf[offset..]) {
            Ok(0) => {
                return if offset == 0 {
                    Ok(false) // clean EOF
                } else {
                    Err(std::io::Error::new(
                        std::io::ErrorKind::UnexpectedEof,
                        format!("truncated chunk: read {} of {} bytes", offset, buf.len()),
                    ))
                };
            }
            Ok(n) => offset += n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(true)
}

/// Generic reader for keyed temp chunks with background prefetching.
///
/// Works with any type implementing `RawSortKey`.
/// Auto-detects BGZF compression via magic bytes.
pub(crate) struct GenericKeyedChunkReader<K: RawSortKey + 'static> {
    receiver: Receiver<ChunkReadResult<K>>,
    /// Return channel for empty buffers — the consumer sends its old buffer
    /// back so the producer can reuse the allocation instead of allocating.
    buf_return: Sender<Vec<u8>>,
    _handle: JoinHandle<()>,
}

impl<K: RawSortKey + 'static> GenericKeyedChunkReader<K> {
    /// Open a keyed chunk file for reading with background prefetching.
    /// Auto-detects BGZF/gzip compression via magic bytes (0x1f 0x8b).
    ///
    /// An optional `concurrency_limit` semaphore can be provided to cap the number
    /// of reader threads actively performing I/O. Readers acquire a token before
    /// reading a batch of records and release it after sending.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened.
    #[allow(clippy::unnecessary_wraps)]
    pub fn open(path: &Path, concurrency_limit: Option<Arc<ChunkReaderSemaphore>>) -> Result<Self> {
        let (tx, rx) = bounded(MERGE_PREFETCH_SIZE);
        let (buf_tx, buf_rx) = bounded::<Vec<u8>>(MERGE_PREFETCH_SIZE);
        let path = path.to_path_buf();

        let handle = thread::spawn(move || {
            let file = match std::fs::File::open(&path) {
                Ok(f) => f,
                Err(e) => {
                    let _ = tx.send(Err(anyhow::anyhow!(
                        "Failed to open keyed chunk {}: {e}",
                        path.display()
                    )));
                    return;
                }
            };
            let mut buf_reader = BufReader::with_capacity(2 * 1024 * 1024, file);

            // Detect magic bytes at file start:
            //   BGZF/gzip  : 0x1f 0x8b
            //   ZSP1 spill : ASCII "ZSP1" followed by [u32 LE len][zstd frame]+
            //
            // Use `read_exact_or_eof` so a short `read()` (legal for `Read`)
            // can't truncate `ZSPILL_MAGIC` into a `None` and silently route a
            // zstd spill through the uncompressed path. Clean EOF (empty file)
            // is preserved as `read_n == 0` — same behavior as before.
            let mut magic = [0u8; 4];
            let read_n = match read_exact_or_eof(&mut buf_reader, &mut magic) {
                Ok(true) => 4,
                Ok(false) => 0,
                Err(e) => {
                    let _ = tx.send(Err(anyhow::anyhow!(
                        "Failed to read keyed chunk magic {}: {e}",
                        path.display()
                    )));
                    return;
                }
            };
            let codec = crate::codec::SpillCodec::from_magic(&magic[..read_n]);

            // Seek back to start so each decoder consumes the magic itself.
            if buf_reader.seek(SeekFrom::Start(0)).is_err() {
                let _ = tx
                    .send(Err(anyhow::anyhow!("Failed to seek in keyed chunk {}", path.display())));
                return;
            }

            match codec {
                Some(crate::codec::SpillCodec::Zstd) => {
                    match crate::zspill_stream::ZspillStreamReader::new(buf_reader) {
                        Ok(rdr) => Self::read_records(rdr, tx, buf_rx, concurrency_limit),
                        Err(e) => {
                            let _ =
                                tx.send(Err(anyhow::anyhow!("zstd spill reader open failed: {e}")));
                        }
                    }
                }
                Some(crate::codec::SpillCodec::Bgzf) => {
                    let bgzf_reader = BgzfReader::new(buf_reader);
                    Self::read_records(bgzf_reader, tx, buf_rx, concurrency_limit);
                }
                None => {
                    // Treat as uncompressed (existing test behavior).
                    Self::read_records(buf_reader, tx, buf_rx, concurrency_limit);
                }
            }
        });

        Ok(Self { receiver: rx, buf_return: buf_tx, _handle: handle })
    }

    /// Open a keyed chunk file with pool-based BGZF decompression.
    ///
    /// Instead of decompressing BGZF blocks on the reader thread, raw blocks are
    /// submitted to the shared `SortWorkerPool` for decompression. This matches
    /// Read records from a reader and send them through the channel.
    ///
    /// When a semaphore is provided, reads records in batches of 64: acquires
    /// a token, reads the batch (I/O + decompression), releases the token,
    /// then sends the batch through the channel. This prevents deadlock — the
    /// token is never held during a blocking `tx.send()`.
    #[allow(clippy::needless_pass_by_value)]
    fn read_records<R: Read>(
        mut reader: R,
        tx: crossbeam_channel::Sender<ChunkReadResult<K>>,
        buf_pool: crossbeam_channel::Receiver<Vec<u8>>,
        semaphore: Option<Arc<ChunkReaderSemaphore>>,
    ) {
        const BATCH_SIZE: usize = 64;

        loop {
            // Phase 1: Acquire token and read a batch of records from disk.
            if let Some(ref sem) = semaphore {
                let _ = sem.1.recv();
            }

            let mut batch: Vec<(K, Vec<u8>)> = Vec::with_capacity(BATCH_SIZE);
            let mut eof = false;
            let mut read_error: Option<String> = None;

            for _ in 0..BATCH_SIZE {
                if K::EMBEDDED_IN_RECORD {
                    // Keyless format: read record, extract key from BAM bytes.
                    let mut len_buf = [0u8; 4];
                    match read_exact_or_eof(&mut reader, &mut len_buf) {
                        Ok(true) => {}
                        Ok(false) => {
                            eof = true;
                            break;
                        }
                        Err(e) => {
                            read_error = Some(format!("Error reading chunk record length: {e}"));
                            break;
                        }
                    }
                    let len = u32::from_le_bytes(len_buf) as usize;
                    let mut record = buf_pool.try_recv().unwrap_or_default();
                    record.clear();
                    record.resize(len, 0);
                    if let Err(e) = reader.read_exact(&mut record) {
                        read_error = Some(format!("Error reading chunk record: {e}"));
                        break;
                    }
                    let key = K::extract_from_record(&record);
                    batch.push((key, record));
                } else {
                    // Keyed format: read key prefix, then record.
                    let key = match K::read_from(&mut reader) {
                        Ok(k) => k,
                        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                            eof = true;
                            break;
                        }
                        Err(e) => {
                            read_error = Some(format!("Error reading keyed chunk key: {e}"));
                            break;
                        }
                    };

                    let mut len_buf = [0u8; 4];
                    match reader.read_exact(&mut len_buf) {
                        Ok(()) => {}
                        Err(e) => {
                            read_error = Some(format!("Error reading keyed chunk length: {e}"));
                            break;
                        }
                    }
                    let len = u32::from_le_bytes(len_buf) as usize;

                    let mut record = buf_pool.try_recv().unwrap_or_default();
                    record.clear();
                    record.resize(len, 0);
                    if let Err(e) = reader.read_exact(&mut record) {
                        read_error = Some(format!("Error reading keyed chunk record: {e}"));
                        break;
                    }

                    batch.push((key, record));
                }
            }

            // Phase 2: Release token before any blocking channel sends.
            if let Some(ref sem) = semaphore {
                let _ = sem.0.send(());
            }

            // Phase 3: Send the batch through the channel (may block).
            for record in batch {
                if tx.send(Ok(Some(record))).is_err() {
                    return; // Receiver dropped
                }
            }

            if let Some(msg) = read_error {
                let _ = tx.send(Err(anyhow::anyhow!("{msg}")));
                break;
            }

            if eof {
                let _ = tx.send(Ok(None));
                break;
            }
        }
    }

    /// Read the next keyed record from the prefetch buffer into `buf`.
    ///
    /// On success the record bytes are swapped into `buf` and the sort key is
    /// returned. The old contents of `buf` are returned to the producer thread
    /// for reuse, avoiding per-record allocation on the disk path.
    ///
    /// # Errors
    ///
    /// Returns an error if the background reader encountered an I/O error.
    pub fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<K>> {
        match self.receiver.recv() {
            Ok(Ok(Some((key, mut data)))) => {
                std::mem::swap(buf, &mut data);
                // Return the old buffer to the producer for reuse.
                let _ = self.buf_return.try_send(data);
                Ok(Some(key))
            }
            Ok(Ok(None)) => Ok(None),
            Ok(Err(e)) => Err(e),
            // Channel disconnected — the producer thread panicked or was dropped
            // without sending an EOF sentinel. Treat as an error, not clean EOF.
            Err(_) => Err(anyhow::anyhow!("chunk reader thread terminated unexpectedly")),
        }
    }

    /// Try to read the next record without blocking.
    ///
    /// Returns `Some(Ok(Some(...)))` if a record is available, `Some(Ok(None))` if the
    /// stream is exhausted, `Some(Err(...))` on read error, or `None` if no record is
    /// currently available (channel empty).
    pub fn try_next_record(&mut self) -> Option<ChunkReadResult<K>> {
        match self.receiver.try_recv() {
            Ok(result) => Some(result),
            Err(crossbeam_channel::TryRecvError::Disconnected) => {
                Some(Err(anyhow::anyhow!("chunk reader thread terminated unexpectedly")))
            }
            Err(crossbeam_channel::TryRecvError::Empty) => None,
        }
    }
}

/// Container for the in-memory chunks passed into the merge.
///
/// Inline-buffer sorts (coordinate, template) produce
/// `Shared` chunks that share an `Arc<SegmentedBuf>` backing store —
/// zero per-record allocation, one memcpy per record at merge time.
///
/// The queryname streaming path accumulates records upstream as
/// individual `RawRecord`s and produces `Owned` chunks — the
/// per-record allocation is sunk cost, so we preserve the original
/// zero-copy `mem::swap` merge bridge.
pub(crate) enum MemorySources<K: RawSortKey + Default + 'static> {
    Shared(Vec<InMemoryChunk<K>>),
    Owned(Vec<Vec<(K, fgumi_raw_bam::RawRecord)>>),
}

impl<K: RawSortKey + Default + 'static> MemorySources<K> {
    fn num_non_empty(&self) -> usize {
        match self {
            Self::Shared(chunks) => chunks.iter().filter(|c| !c.is_empty()).count(),
            Self::Owned(chunks) => chunks.iter().filter(|c| !c.is_empty()).count(),
        }
    }

    fn push_into(self, sources: &mut Vec<ChunkSource<K>>) {
        match self {
            Self::Shared(chunks) => {
                for chunk in chunks {
                    if !chunk.is_empty() {
                        sources.push(ChunkSource::Memory { chunk, idx: usize::MAX });
                    }
                }
            }
            Self::Owned(chunks) => {
                for records in chunks {
                    if !records.is_empty() {
                        sources.push(ChunkSource::MemoryOwned { records, idx: usize::MAX });
                    }
                }
            }
        }
    }
}

/// Source for keyed chunks during merge (disk or in-memory).
///
/// Each variant owns its "current record" state internally so the merge loop
/// can borrow the current record's bytes via [`Self::current_bytes`] without
/// copying. The `PoolDisk` variant owns a `scratch: Vec<u8>` holding the most
/// recently read record's bytes; `Memory`/`MemoryOwned` borrow directly from
/// their backing store, eliminating the per-record memcpy on the in-memory path.
enum ChunkSource<K: RawSortKey + Default + 'static> {
    /// In-memory chunk from an inline-buffer sort (coordinate / template).
    /// All records borrow their bytes from a shared `Arc<SegmentedBuf>`;
    /// `current_bytes` borrows them directly (zero-copy).
    ///
    /// `idx` is the index of the CURRENT record. Sentinel `usize::MAX` means
    /// "no record loaded yet"; the first `advance` lands on 0.
    Memory { chunk: InMemoryChunk<K>, idx: usize },
    /// In-memory chunk from a queryname-style sort (each record an owned
    /// `RawRecord`). `current_bytes` borrows `records[idx].1` directly.
    ///
    /// `idx` uses the same `usize::MAX` sentinel as `Memory`.
    MemoryOwned { records: Vec<(K, fgumi_raw_bam::RawRecord)>, idx: usize },
    /// Pool-integrated disk source — workers read and decompress, main thread parses.
    PoolDisk {
        source_id: usize,
        /// Bytes for the current record (filled by the most recent `advance`).
        scratch: Vec<u8>,
    },
}

impl<K: RawSortKey + Default + 'static> ChunkSource<K> {
    /// Advance to the next record. Returns its key, or `None` at EOF.
    ///
    /// For `PoolDisk` sources, `consumer` must be `Some`. After this returns
    /// `Some(_)`, callers may invoke [`Self::current_bytes`] to borrow the
    /// record's bytes without copying.
    fn advance(&mut self, consumer: Option<&mut MainThreadChunkConsumer<K>>) -> Result<Option<K>> {
        match self {
            ChunkSource::PoolDisk { source_id, scratch } => {
                let c = consumer.ok_or_else(|| {
                    anyhow::anyhow!(
                        "PoolDisk source (id {source_id}) requires a MainThreadChunkConsumer \
                         but none was provided — this is a bug in the sort pipeline"
                    )
                })?;
                c.next_record(*source_id, scratch)
            }
            ChunkSource::Memory { chunk, idx } => {
                // Sentinel `usize::MAX` → first advance lands on 0.
                let next = if *idx == usize::MAX { 0 } else { *idx + 1 };
                if next < chunk.len() {
                    let key = chunk.take_key(next);
                    *idx = next;
                    Ok(Some(key))
                } else {
                    Ok(None)
                }
            }
            ChunkSource::MemoryOwned { records, idx } => {
                let next = if *idx == usize::MAX { 0 } else { *idx + 1 };
                if next < records.len() {
                    let key = std::mem::take(&mut records[next].0);
                    *idx = next;
                    Ok(Some(key))
                } else {
                    Ok(None)
                }
            }
        }
    }

    /// Bytes for the current record (the one whose key was returned by the
    /// most recent successful [`Self::advance`]).
    ///
    /// Caller MUST have received `Some(_)` from a prior `advance`. The merge
    /// loops uphold this: the init loop only retains a source whose first
    /// `advance` returned `Some`, and the main loop only calls this on the
    /// active tree winner — so the `PoolDisk` empty-scratch case is
    /// unreachable.
    fn current_bytes(&self) -> &[u8] {
        match self {
            ChunkSource::PoolDisk { scratch, .. } => scratch,
            ChunkSource::Memory { chunk, idx } => {
                debug_assert!(*idx != usize::MAX, "current_bytes called before advance");
                chunk.record_bytes(*idx)
            }
            ChunkSource::MemoryOwned { records, idx } => {
                debug_assert!(*idx != usize::MAX, "current_bytes called before advance");
                &records[*idx].1[..]
            }
        }
    }
}

/// Bytes of `source`'s current record to hand to the merge writer.
///
/// For `PoolDisk` sources this borrows directly from the pooled decompressed
/// block when the record fit in one block (the zero-copy fast path), falling
/// back to the source's scratch buffer for a record that straddled a block
/// boundary. All other source variants borrow from their own backing store via
/// [`ChunkSource::current_bytes`] and ignore `consumer`.
///
/// `consumer` must be `Some` whenever `source` is `PoolDisk` (a `PoolDisk`
/// source only exists when Phase 2 is active, which is exactly when the merge
/// loops set the consumer); a missing consumer is a pipeline bug and returns an
/// error here — matching [`ChunkSource::advance`]'s handling of the same
/// invariant — rather than silently writing empty/stale scratch bytes.
/// How the merge presented each record to the writer: borrowed straight from
/// the decompressed block, or reassembled into scratch because it straddled a
/// block boundary.
///
/// Process-wide rather than per-merge because the merge loop hands out `&[u8]`
/// and threading a counter through it would change the signature of the hot
/// path. These are never reset — concurrent sorts sharing the process would
/// race a reset — so the merge loop snapshots them before and after its run and
/// reports only the delta (see `log_merge_sub_phases`) rather than the running
/// total. The delta is exact for sequential sorts; if another sort runs
/// concurrently in the same process its records fall inside the window and
/// inflate the delta, but that is a diagnostic-only skew and is preferable to a
/// race-prone reset.
static RECORD_BORROWED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
static RECORD_REASSEMBLED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Consumer-side per-merge diagnostics for the sub-phase log.
///
/// Grouped into one value so they travel together and keep
/// [`RawExternalSorter::log_merge_sub_phases`] within its argument budget. All
/// four are scoped to a single merge: the backpressure fields are measured over
/// this merge's loop, and the presentation counts are a delta of the
/// process-wide [`RECORD_BORROWED`]/[`RECORD_REASSEMBLED`] atomics taken across
/// the loop.
#[derive(Clone, Copy)]
struct MergeConsumerDiag {
    /// Seconds the consumer blocked waiting for an output permit (exact).
    backpressure_secs: f64,
    /// Number of output-permit waits.
    backpressure_waits: u64,
    /// Records handed to the writer borrowed zero-copy from a decompressed block.
    borrowed: u64,
    /// Records reassembled into scratch because they straddled a block boundary.
    reassembled: u64,
}

#[inline]
fn winner_record_bytes<'a, K: RawSortKey + Default + 'static>(
    source: &'a ChunkSource<K>,
    consumer: Option<&'a MainThreadChunkConsumer<K>>,
) -> Result<&'a [u8]> {
    match source {
        ChunkSource::PoolDisk { source_id, scratch } => {
            let consumer = consumer.ok_or_else(|| {
                anyhow::anyhow!(
                    "PoolDisk source (id {source_id}) requires a MainThreadChunkConsumer \
                     during merge but none was provided — this is a bug in the sort pipeline"
                )
            })?;
            // `Some` on the zero-copy fast path; `None` means the record
            // straddled a block boundary and was reassembled into `scratch`.
            //
            // Counted, not timed: this runs once per record on the one thread
            // that touches every record, so a clock read here would cost more
            // than the step. What was unknown is how often the copy path fires
            // at all -- a frequency answers that, and the per-record cost is a
            // memcpy of a ~100-byte record either way.
            if let Some(borrowed) = consumer.current_record_bytes(*source_id) {
                RECORD_BORROWED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                Ok(borrowed)
            } else {
                RECORD_REASSEMBLED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                Ok(scratch.as_slice())
            }
        }
        other => Ok(other.current_bytes()),
    }
}

// ============================================================================
// MainThreadChunkConsumer — pool-integrated merge reader
// ============================================================================
//
// The pool's worker threads read raw BGZF blocks from disk, decompress them,
// and insert the decompressed blocks into per-file `ReorderBuffer`s held in
// each `Phase2FileState`. The main thread (the merge loop) holds a snapshot of
// those file states and pulls one decompressed block at a time per source as
// the loser tree advances. There is no global queue and no per-source
// reorder/buffering layer in the consumer itself — the per-file
// `Phase2FileState.decompressed` reorder buffer IS the per-source buffer.
//
// Backpressure is enforced inside the worker pool: when a per-file reorder
// buffer reaches `PHASE2_DECOMP_CAP`, workers stop pulling new raw blocks for
// that file (with a deadlock-free escape hatch for the gap-filling serial).

/// Where the current record's bytes live for zero-copy presentation to the
/// merge loop.
///
/// The common case (`Borrowed`) is a record that lies contiguously inside the
/// source's current decompressed block, so it can be handed to the writer as a
/// slice with no intermediate copy. A record that straddles a block boundary is
/// reassembled into the source's scratch buffer (`Scratch`) as before.
#[derive(Clone, Copy)]
enum CurRecord {
    /// Record occupies `current_buf[start..end]` — borrow it directly.
    Borrowed { start: usize, end: usize },
    /// Record was reassembled into the `ChunkSource`'s scratch buffer.
    Scratch,
}

/// Per-source byte-stream parser state owned by the main thread.
///
/// As blocks are pulled from `Phase2FileState.decompressed`, the bytes are
/// stashed in `current_buf` and consumed left-to-right by the record parser.
struct SourceParserState {
    /// Current decompressed block being consumed.
    current_buf: Vec<u8>,
    /// Read position within `current_buf`.
    current_pos: usize,
    /// Location of the most recently parsed record's bytes.
    cur_record: CurRecord,
}

impl SourceParserState {
    fn new() -> Self {
        Self { current_buf: Vec::new(), current_pos: 0, cur_record: CurRecord::Scratch }
    }

    fn remaining(&self) -> usize {
        self.current_buf.len() - self.current_pos
    }
}

/// Reads from per-file decompressed-block reorder buffers and presents records
/// to the merge loop.
///
/// The main thread drives all record consumption; no threads are spawned here.
/// Sort-pool workers do the disk reads and BGZF decompression in parallel and
/// publish results into per-file `Phase2FileState.decompressed` reorder
/// buffers, which this consumer drains in serial order.
///
/// # Type Parameter
///
/// `K` is the sort key type (`RawCoordinateKey`, `TemplateKey`, etc.).
pub(crate) struct MainThreadChunkConsumer<K: RawSortKey + 'static> {
    /// Snapshot of the pool's Phase 2 file vector. Indexed by `source_id`.
    files: Arc<Vec<crate::worker_pool::Phase2FileState>>,
    /// Per-source parser state.
    parser_state: Vec<SourceParserState>,
    /// Set by a pool worker when BGZF decompression of a chunk block fails.
    decompression_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Set by a pool worker when a chunk file I/O read fails.
    chunk_read_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Set by `do_shutdown` when a worker thread panicked unexpectedly.
    worker_panicked: std::sync::Arc<std::sync::atomic::AtomicBool>,
    /// Where and for how long the merge loop blocked. Only this thread touches
    /// it -- see [`crate::merge_stalls::ConsumerStallTracker`].
    stalls: crate::merge_stalls::ConsumerStallTracker,
    /// Shared pool state, for the epoch clock and the trace counters that
    /// record both halves of a producer/consumer handoff.
    shared: Arc<crate::worker_pool::SharedPipelineState>,
    /// Source the previous block came from, and how many consecutive blocks
    /// have now come from it. Says whether the merge dwells on one run at a
    /// time -- in which case lookahead on that run would pay -- or hops.
    current_run: Option<(usize, u64)>,
    _phantom: std::marker::PhantomData<K>,
}

/// Enforce the permutation invariant, removing the output when it fails.
///
/// A sort is a permutation: every record read must be written. By the time this
/// runs the writer has finalized `output`, so a mismatch leaves a
/// truncated-but-*valid* BAM (and, under `write_index`, a `.bai` describing it)
/// on disk. Saying the output "must not be used" does not stop a wrapper that
/// checks only an exit code, or a user who re-runs and reads the stale file, so
/// a file output is removed rather than merely described.
///
/// That also matches the rest of the crate rather than making this the one
/// exception: the merge writes through [`MergeOutputTarget`] and drops its temp
/// on error, and `test_sort_records_producer_error_aborts` pins that an aborted
/// sort leaves no output behind.
///
/// What is *not* removed is as deliberate as what is, and mirrors the guards
/// [`MergeOutputTarget::create`] already applies: stdout cannot be un-written,
/// and a FIFO, device or directory was written *through* rather than created, so
/// unlinking it would destroy something this sort does not own. Those cases keep
/// the "must not be used" wording, which is all that is true of them. A removal
/// that fails is reported rather than swallowed -- the caller is being told the
/// output is dangerous, so whether it is still there is exactly what they need.
fn enforce_record_count(stats: &RawSortStats, output: &Path, write_index: bool) -> Result<()> {
    if stats.total_records == stats.output_records {
        return Ok(());
    }

    let disposition = if is_stdout_path(output) {
        // The records are already past this process; there is nothing to unlink,
        // and `/dev/stdout` itself must certainly not be.
        "The records written to stdout are incomplete and must not be used.".to_string()
    } else {
        let mut targets = vec![output.to_path_buf()];
        if write_index {
            targets.push(fgumi_bam_io::bai_sidecar_path(output));
        }
        let left_behind: Vec<PathBuf> =
            targets.into_iter().filter(|path| !remove_incomplete_output(path)).collect();
        if left_behind.is_empty() {
            format!("The incomplete output at {} has been removed.", output.display())
        } else {
            let paths =
                left_behind.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ");
            format!("The incomplete output at {paths} could not be removed and must not be used.")
        }
    };

    anyhow::bail!(
        "sort lost records: read {} but wrote {} (differ by {}). {}",
        stats.total_records,
        stats.output_records,
        stats.total_records.abs_diff(stats.output_records),
        disposition
    )
}

/// Remove an incomplete sort output, reporting whether it is now gone.
///
/// Follows a symlink to the file the sort actually wrote through, so unlinking
/// the link would not leave the real short BAM in place -- the same resolution
/// [`MergeOutputTarget::create`] does for the opposite reason. Only regular
/// files are removed; anything else was written through, not created, and a
/// path that cannot be resolved or stat-ed is reported as still present rather
/// than assumed gone.
fn remove_incomplete_output(path: &Path) -> bool {
    let Ok(resolved) = resolve_symlink_output(path) else {
        return false;
    };
    match std::fs::metadata(&resolved) {
        // Already absent is the outcome we wanted, not a failure.
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => true,
        Err(_) => false,
        Ok(metadata) if !metadata.is_file() => false,
        Ok(_) => std::fs::remove_file(&resolved).is_ok(),
    }
}

/// Whether the merge block-lifecycle block has nothing to print.
///
/// Every report that block prints must appear here. A report missing from the
/// gate is dropped whenever it is the only one populated -- which is exactly
/// when it is the whole story. `scans` is the case that bit: the fruitless-scan
/// cost is the figure the `WorkUnclaimed` verdict tells the reader to compare
/// against, so a scans-only run lost the number its own verdict pointed at.
///
/// Pure so the coverage is testable without capturing log output, following
/// [`classify_scan`](crate::merge_stalls::classify_scan).
fn block_lifecycle_is_silent(
    life: &crate::merge_trace::BlockLifecycleReport,
    refill: &crate::merge_trace::RefillReport,
    consumer: &crate::merge_trace::ConsumerTraceReport,
    scans: &crate::merge_trace::HistogramReport,
) -> bool {
    life.is_empty() && refill.is_empty() && consumer.is_empty() && scans.is_empty()
}

/// Whether the merge stall block has nothing to print.
///
/// Deliberately separate from [`block_lifecycle_is_silent`], and never a gate on
/// it: these three reports are populated by the pipeline *stalling*, the
/// lifecycle reports by it *working*. A merge that flowed cleanly leaves this
/// gate silent and the lifecycle gate loud, so letting this one decide whether
/// the lifecycle block prints would drop the whole trace on exactly the runs
/// that produced a complete one.
///
/// Pure so the coverage is testable without capturing log output, following
/// [`classify_scan`](crate::merge_stalls::classify_scan). The emptiness check on
/// `stalls` is deliberately kept here even though the caller already drops an
/// empty report: the caller drops it for its own reason -- an empty report must
/// not reach the block that prints it -- and a gate that answers "is there
/// anything to print" only when someone else has pre-filtered its input is a
/// gate whose test proves nothing.
fn merge_stalls_are_silent(
    stalls: Option<&crate::merge_stalls::ConsumerStallReport>,
    scans: &crate::merge_stalls::Phase2ScanReport,
    wake: &crate::merge_stalls::WakeLatencyReport,
) -> bool {
    stalls.is_none_or(|s| s.is_empty()) && scans.is_empty() && wake.is_empty()
}

impl<K: RawSortKey + 'static> MainThreadChunkConsumer<K> {
    /// Create a new consumer for the given pool file snapshot.
    #[must_use]
    pub(crate) fn new(
        files: Arc<Vec<crate::worker_pool::Phase2FileState>>,
        decompression_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
        chunk_read_error: std::sync::Arc<std::sync::atomic::AtomicBool>,
        worker_panicked: std::sync::Arc<std::sync::atomic::AtomicBool>,
        shared: Arc<crate::worker_pool::SharedPipelineState>,
    ) -> Self {
        let parser_state = (0..files.len()).map(|_| SourceParserState::new()).collect();
        let stalls = crate::merge_stalls::ConsumerStallTracker::new(files.len());
        Self {
            files,
            parser_state,
            decompression_error,
            chunk_read_error,
            worker_panicked,
            stalls,
            shared,
            current_run: None,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Note that the merge is taking a block from `source_id`, closing the
    /// previous source's run if it switched.
    fn note_source_run(&mut self, source_id: usize) {
        match self.current_run {
            Some((prev, blocks)) if prev == source_id => {
                self.current_run = Some((prev, blocks + 1));
            }
            Some((_, blocks)) => {
                self.shared.consumer_trace.record_source_run(blocks);
                self.current_run = Some((source_id, 1));
            }
            None => self.current_run = Some((source_id, 1)),
        }
    }

    /// Flush the run in progress, so the last one is not lost at end of merge.
    pub(crate) fn finish_source_run(&mut self) {
        if let Some((_, blocks)) = self.current_run.take() {
            self.shared.consumer_trace.record_source_run(blocks);
        }
    }

    /// Start the stall accounting from now, discarding what came before.
    ///
    /// Called once the loser tree is seeded, so the stall figures cover the same
    /// interval as `loop_total` -- see
    /// [`ConsumerStallTracker::restart`](crate::merge_stalls::ConsumerStallTracker::restart).
    pub(crate) fn restart_stalls(&mut self) {
        self.stalls.restart();
        // The seeding pulls opened a run that belongs to the same discarded
        // interval; drop it rather than let it extend the merge's first run.
        self.current_run = None;
    }

    /// Where and for how long the merge loop blocked waiting for blocks.
    pub(crate) fn stall_report(&self) -> crate::merge_stalls::ConsumerStallReport {
        self.stalls.snapshot()
    }

    /// Observe every file's buffering state, for attributing one park.
    ///
    /// `try_lock` throughout: this runs while the merge is stalled and must not
    /// add to the stall it is measuring, nor perturb the very contention it is
    /// trying to detect. A file whose locks are both held is reported as
    /// `contended` rather than guessed at. The two locks are taken and released
    /// one at a time, so this introduces no lock ordering of its own.
    fn census(&self, waited_on: usize) -> crate::merge_stalls::ParkCensus {
        use crate::merge_stalls::AwaitedState;
        use crate::worker_pool::PHASE2_DECOMP_CAP;

        let mut census = crate::merge_stalls::ParkCensus::default();
        for (idx, file) in self.files.iter().enumerate() {
            let Ok(raw_len) = file.raw_blocks.try_lock().map(|g| g.len()) else {
                census.contended += 1;
                continue;
            };
            let Ok(decomp_len) = file.decompressed.try_lock().map(|g| g.len()) else {
                census.contended += 1;
                continue;
            };
            let in_flight = file.decomp_in_flight.load(std::sync::atomic::Ordering::Relaxed);
            let at_eof = file.reader_eof.load(std::sync::atomic::Ordering::Relaxed);
            let empty = raw_len == 0 && decomp_len == 0 && in_flight == 0;

            if idx == waited_on {
                // The caller has just failed `try_pop_next` on this file, so a
                // non-empty reorder buffer is proof it holds serials other than
                // the one wanted -- a gap, not merely "some data". That is what
                // makes these five states a decision rather than an inference.
                census.awaited = Some(AwaitedState::classify(decomp_len, in_flight, raw_len));
            }

            if decomp_len >= PHASE2_DECOMP_CAP {
                census.capped += 1;
            } else if empty && at_eof {
                census.drained += 1;
            } else if empty {
                census.starved += 1;
            } else {
                census.working += 1;
            }
        }
        census
    }

    /// Get the next record from a specific source.
    ///
    /// Pulls decompressed blocks from the source's per-file reorder buffer as
    /// needed. Parses the next record from the source's byte stream.
    ///
    /// Returns `Ok(Some(key))` with record bytes swapped into `buf`, `Ok(None)` at EOF.
    ///
    /// # Errors
    ///
    /// Returns an error if a record is truncated or if a worker reported an error.
    pub fn next_record(&mut self, source_id: usize, buf: &mut Vec<u8>) -> Result<Option<K>> {
        // Make sure we have data (or detect EOF) before parsing.
        if self.parser_state[source_id].remaining() == 0
            && !self.advance_to_next_block(source_id)?
        {
            return Ok(None);
        }
        self.parse_next_record(source_id, buf)
    }

    /// Bytes of `source_id`'s current record when they were borrowed directly
    /// from the decompressed block (the zero-copy fast path). Returns `None`
    /// when the record straddled a block boundary and was reassembled into the
    /// source's scratch buffer, in which case the caller reads it from there.
    ///
    /// The returned slice is valid until the next `advance` of this source
    /// (which is what may replace `current_buf`); the merge loop always writes
    /// the current record before advancing its source, so this holds.
    #[inline]
    fn current_record_bytes(&self, source_id: usize) -> Option<&[u8]> {
        let st = &self.parser_state[source_id];
        match st.cur_record {
            CurRecord::Borrowed { start, end } => Some(&st.current_buf[start..end]),
            CurRecord::Scratch => None,
        }
    }

    /// Pull the next decompressed block for `source_id` from the per-file
    /// reorder buffer, blocking the main thread (`std::thread::park`) until
    /// either a block becomes available, the source drains, or a worker error
    /// is reported.
    ///
    /// Returns `Ok(true)` if a new block was loaded into `current_buf`,
    /// `Ok(false)` if the source has produced all its data, or an error if a
    /// worker reported a fatal failure.
    fn advance_to_next_block(&mut self, source_id: usize) -> Result<bool> {
        self.stalls.note_block_pull();
        let mut parked_yet = false;
        // The file is re-borrowed per iteration rather than held: the stall
        // tracker below needs `&mut self`, and a long-lived `&self.files[..]`
        // would keep `self` immutably borrowed across the whole loop.
        loop {
            // Try to pop the next-in-order decompressed block.
            //
            // Whether a block was actually taken is carried out of the block
            // scope rather than returned from inside it, because crediting the
            // source's run needs `&mut self` and the borrow of `self.files`
            // taken below is still live in there.
            let mut popped = false;
            {
                let file = &self.files[source_id];
                let mut guard =
                    file.decompressed.lock().expect("phase2 decompressed mutex poisoned");
                if let Some(block) = guard.try_pop_next() {
                    let remaining = guard.len();
                    file.decomp_len.store(remaining, std::sync::atomic::Ordering::Relaxed);
                    // This pop may have drained the file, which opens a refill
                    // cycle. Opening it here, at the instant the buffer hits
                    // zero rather than when the consumer next comes back for a
                    // block, measures the pipeline's response time rather than
                    // the consumer's round trip.
                    //
                    // Opened while this mutex is still held, because a worker
                    // closes the cycle under the same mutex and the two must not
                    // interleave: otherwise an insert could observe a cycle,
                    // have this pop close it and open the next, and then record
                    // its latency against the new cycle and clear it -- losing
                    // that cycle's real measurement and booking a zero for it.
                    //
                    // Only a drain pays for the extra work in the critical
                    // section, so the cost falls once per refill cycle rather
                    // than once per block, and the clock read and the histogram
                    // record both stay outside the lock on the common path.
                    let emptied_cause = if remaining == 0 {
                        let cause = crate::merge_trace::EmptyCause::classify(
                            file.raw_len.load(std::sync::atomic::Ordering::Relaxed),
                            file.decomp_in_flight.load(std::sync::atomic::Ordering::Relaxed),
                        );
                        file.mark_emptied(self.shared.now_nanos(), cause);
                        Some(cause)
                    } else {
                        None
                    };
                    drop(guard);
                    if emptied_cause.is_some() {
                        // This file needs a worker now. Waking one here is the
                        // difference between the pool learning that within a
                        // microsecond and learning it whenever some worker's
                        // backoff happens to expire. After the drop, so the
                        // woken worker does not immediately block on the mutex
                        // this thread would still be holding.
                        self.shared.wake_one_worker();
                    }
                    if let Some(cause) = emptied_cause {
                        self.shared.refill.record_empty(cause);
                    }
                    let now = self.shared.now_nanos();
                    self.shared
                        .block_lifecycle
                        .reorder_dwell
                        .record(now.saturating_sub(block.inserted_nanos));
                    let st = &mut self.parser_state[source_id];
                    st.current_buf = block.data;
                    st.current_pos = 0;
                    popped = true;
                }
            }
            if popped {
                // Only a block that was actually consumed extends the source's
                // run. Crediting the call instead adds one phantom block to
                // every completed run, because a source at EOF is pulled from
                // once more and comes back empty.
                self.note_source_run(source_id);
                return Ok(true);
            }

            // No block ready. Check error flags first — they take precedence
            // over EOF detection so a single-source sort that fails on its
            // last block surfaces the error rather than silently truncating.
            if self.decompression_error.load(std::sync::atomic::Ordering::Acquire) {
                return Err(anyhow::anyhow!(
                    "BGZF decompression error on chunk blocks (see log for details)"
                ));
            }
            if self.chunk_read_error.load(std::sync::atomic::Ordering::Acquire) {
                return Err(anyhow::anyhow!("I/O error reading chunk file (see log for details)"));
            }
            if self.worker_panicked.load(std::sync::atomic::Ordering::Acquire) {
                return Err(anyhow::anyhow!(
                    "a sort worker thread panicked unexpectedly (see log for details)"
                ));
            }

            // Source produced everything it ever will?
            if self.files[source_id].is_drained() {
                // Move the drain frontier past it, so workers scanning for the
                // next hungry file start at one that still has records.
                self.shared.advance_phase2_frontier(&self.files);
                return Ok(false);
            }

            // Park until a worker unparks us. Workers unpark after pushing a
            // decompressed block, after setting reader.eof, and after error
            // flags. The loop re-checks all conditions on wake-up so spurious
            // wake-ups are harmless.
            //
            // This is where the merge loop actually blocks, and it is timed
            // exactly rather than sampled: parks happen once per ~64 KB block
            // instead of once per ~100-byte record, so an `Instant` pair here is
            // three orders of magnitude cheaper than one per record -- and it
            // separates *waiting* from the record copy that the sampled "fetch
            // next record" figure folds in with it.
            if self.stalls.should_census() {
                let census = self.census(source_id);
                self.stalls.record_census(census);
            }
            // Classify the awaited file on *every* park, not just the censused
            // ones: the depths are plain atomics, so unlike the pool-wide
            // census this costs three relaxed loads rather than 86 `try_lock`s.
            // Park durations are heavy-tailed, so the per-state distribution is
            // exactly the thing a sampled version would get wrong.
            let (raw_len, decomp_len, in_flight) = self.files[source_id].depths();
            let state = crate::merge_stalls::AwaitedState::classify(decomp_len, in_flight, raw_len);
            // Tell the pool which file the merge is blocked on, so it reads
            // ahead deeply there. The drain frontier is the same file only when
            // sources drain in order; on a partially-correlated input the merge
            // parks on a source that is not the lowest undrained one, and that
            // source was reading at the shallow allowance.
            self.shared
                .phase2_awaited_source
                .store(source_id, std::sync::atomic::Ordering::Relaxed);
            // Arm the park handshake. Workers stamp when the awaited block is
            // first claimed and when the consumer is woken; the difference from
            // `park_start` splits this park into additive stages. Cleared here
            // rather than after so a stamp can only ever be attributed to a park
            // that had already begun.
            let ord = std::sync::atomic::Ordering::Relaxed;
            self.shared.awaited_claim_nanos.store(0, ord);
            self.shared.awaited_publish_nanos.store(0, ord);
            let park_started_at = self.shared.now_nanos();
            // Sampled BEFORE the park, not after: the question is why nobody had
            // already started on this block, and the pool's state on *resume* is
            // the answer to a different question -- by then somebody has.
            let supply = self.shared.park_supply_now();

            let park_start = Instant::now();
            std::thread::park();
            let parked_ns = crate::merge_trace::elapsed_nanos(park_start);
            self.shared.park_supply.record(supply, parked_ns);
            self.stalls.record_park(source_id, parked_ns, !parked_yet);
            self.shared.consumer_trace.record_park(state, parked_ns, in_flight);

            // Split the park. `now_nanos` is the same clock the workers stamp
            // with, so the segments are comparable; `parked_ns` comes from
            // `Instant` and is only used for the existing counters.
            let resumed_at = self.shared.now_nanos();
            let claim = self.shared.awaited_claim_nanos.load(ord);
            let publish = self.shared.awaited_publish_nanos.load(ord);
            let segments = crate::merge_stalls::split_park(
                park_started_at,
                (claim != 0).then_some(claim),
                (publish != 0).then_some(publish),
                resumed_at,
            );
            // Depth on the critical path: how many blocks the awaited file had
            // ready the moment the consumer could run again. A mean near 1 means
            // every block is fetched on demand, which is a different problem
            // from a slow fetch.
            let (_, ready, _) = self.files[source_id].depths();
            self.shared.park_attribution.record(segments, claim != 0, ready as u64);
            parked_yet = true;
        }
    }

    /// Parse the next record from a source's byte stream.
    ///
    /// Handles the format: for `EMBEDDED_IN_RECORD` keys, reads [len(4)][record(len)].
    /// For keyed format, reads [key][len(4)][record(len)].
    fn parse_next_record(&mut self, source_id: usize, buf: &mut Vec<u8>) -> Result<Option<K>> {
        let mut len_buf = [0u8; 4];

        if K::EMBEDDED_IN_RECORD {
            if !self.read_exact_from_source(source_id, &mut len_buf)? {
                return Ok(None);
            }
            let len = u32::from_le_bytes(len_buf) as usize;

            // Fast path: the whole record is contiguous in the current block, so
            // borrow it in place rather than copying it into `buf` (the source's
            // scratch). A record only straddles a block boundary right at the
            // ~64KB block edges, so this hits for the vast majority of records
            // and removes one full-record memcpy per record from the merge.
            if self.parser_state[source_id].remaining() >= len {
                let start = self.parser_state[source_id].current_pos;
                let end = start + len;
                {
                    let st = &mut self.parser_state[source_id];
                    st.current_pos = end;
                    st.cur_record = CurRecord::Borrowed { start, end };
                }
                let key =
                    K::extract_from_record(&self.parser_state[source_id].current_buf[start..end]);
                return Ok(Some(key));
            }

            // Slow path: record spans a block boundary — reassemble into scratch.
            buf.clear();
            buf.resize(len, 0);
            if !self.read_exact_from_source(source_id, buf)? {
                return Err(anyhow::anyhow!("truncated record in chunk source {source_id}"));
            }
            self.parser_state[source_id].cur_record = CurRecord::Scratch;
            let key = K::extract_from_record(buf);
            Ok(Some(key))
        } else {
            let Some(key) = self.read_key_from_source::<K>(source_id)? else {
                return Ok(None);
            };

            if !self.read_exact_from_source(source_id, &mut len_buf)? {
                return Err(anyhow::anyhow!("truncated record length in chunk source {source_id}"));
            }
            let len = u32::from_le_bytes(len_buf) as usize;

            // Fast path, same as the embedded case: the key and length have been
            // consumed, so if the record itself is contiguous in the current
            // block it can be borrowed in place. Previously every keyed record
            // was copied into scratch, which on a template-coordinate merge is a
            // full-record memcpy per record on the one thread that touches every
            // record -- measured at 776,167,078 of 776,167,078 records (100%),
            // against a comment claiming the borrow path took "the vast
            // majority". Nothing about a keyed record makes it uncopyable; the
            // borrow simply had not been extended past the embedded path.
            if self.parser_state[source_id].remaining() >= len {
                let st = &mut self.parser_state[source_id];
                let start = st.current_pos;
                let end = start + len;
                st.current_pos = end;
                st.cur_record = CurRecord::Borrowed { start, end };
                return Ok(Some(key));
            }

            // Slow path: the record spans a block boundary, so it has to be
            // reassembled into the source's scratch buffer.
            buf.clear();
            buf.resize(len, 0);
            if !self.read_exact_from_source(source_id, buf)? {
                return Err(anyhow::anyhow!("truncated record in chunk source {source_id}"));
            }
            self.parser_state[source_id].cur_record = CurRecord::Scratch;
            Ok(Some(key))
        }
    }

    /// Read exactly `out.len()` bytes from a source into `out`, pulling more
    /// blocks from the per-file reorder buffer as needed.
    ///
    /// Returns `Ok(false)` at clean EOF (zero bytes available), `Ok(true)` on success.
    fn read_exact_from_source(&mut self, source_id: usize, out: &mut [u8]) -> Result<bool> {
        let n = out.len();
        let mut filled = 0;

        while filled < n {
            if self.parser_state[source_id].remaining() == 0
                && !self.advance_to_next_block(source_id)?
            {
                if filled == 0 {
                    return Ok(false);
                }
                return Err(anyhow::anyhow!(
                    "truncated data in chunk source {source_id}: got {filled} of {n} bytes",
                ));
            }

            let st = &mut self.parser_state[source_id];
            let take = (n - filled).min(st.remaining());
            out[filled..filled + take]
                .copy_from_slice(&st.current_buf[st.current_pos..st.current_pos + take]);
            st.current_pos += take;
            filled += take;
        }

        Ok(true)
    }

    /// Read a sort key from a source's byte stream.
    ///
    /// Returns `Ok(None)` at clean EOF.
    fn read_key_from_source<KK: RawSortKey>(&mut self, source_id: usize) -> Result<Option<KK>> {
        let mut adapter = SourceReadAdapter { consumer: self, source_id, bytes_read: 0 };
        match KK::read_from(&mut adapter) {
            Ok(key) => Ok(Some(key)),
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
                if adapter.bytes_read == 0 {
                    Ok(None)
                } else {
                    Err(anyhow::anyhow!(
                        "truncated key in chunk source {source_id}: \
                         got {n} bytes then EOF",
                        n = adapter.bytes_read
                    ))
                }
            }
            Err(e) => Err(anyhow::anyhow!("error reading key from source {source_id}: {e}")),
        }
    }

    /// Gather probe statistics from the per-file Phase 2 state.
    fn probe_consumer_stats(&self) -> ConsumerProbeStats {
        let mut pending_blocks: u64 = 0;
        let mut pending_bytes: u64 = 0;
        let mut active_sources: u64 = 0;

        for file in self.files.iter() {
            let (blocks, bytes, active) = file.probe_stats();
            pending_blocks += blocks;
            pending_bytes += bytes;
            if active {
                active_sources += 1;
            }
        }

        ConsumerProbeStats {
            current_bytes: 0,
            current_capacity: 0,
            pending_blocks,
            pending_bytes,
            active_sources,
        }
    }
}

/// Adapter that implements `std::io::Read` over a `MainThreadChunkConsumer` source.
///
/// This allows `K::read_from(&mut reader)` to read from the pool-based byte stream.
/// `bytes_read` tracks how many bytes have been consumed so `read_key_from_source` can
/// distinguish a clean EOF (zero bytes seen) from a truncated read (some bytes then EOF).
struct SourceReadAdapter<'a, K: RawSortKey + 'static> {
    consumer: &'a mut MainThreadChunkConsumer<K>,
    source_id: usize,
    bytes_read: usize,
}

impl<K: RawSortKey + 'static> Read for SourceReadAdapter<'_, K> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        // Pull the next block if the current one is exhausted.
        if self.consumer.parser_state[self.source_id].remaining() == 0 {
            match self.consumer.advance_to_next_block(self.source_id) {
                Ok(true) => {}
                Ok(false) => return Ok(0), // clean EOF
                Err(e) => return Err(std::io::Error::other(e.to_string())),
            }
        }

        let st = &mut self.consumer.parser_state[self.source_id];
        let take = buf.len().min(st.remaining());
        buf[..take].copy_from_slice(&st.current_buf[st.current_pos..st.current_pos + take]);
        st.current_pos += take;
        self.bytes_read += take;
        Ok(take)
    }
}

/// Generates unique file paths for chunk and merged temp files.
///
/// Maintains monotonic counters for both chunk files (`chunk_0000.keyed`, ...)
/// and merged files (`merged_0000.keyed`, ...) to prevent naming collisions
/// after consolidation drains entries from the chunk file list.
///
/// When multiple temp directories are supplied via [`TmpDirAllocator`], chunk
/// and merged files are distributed across them in round-robin order.
struct ChunkNamer<'a> {
    alloc: &'a mut TmpDirAllocator,
    chunk_count: usize,
    merge_count: usize,
}

impl<'a> ChunkNamer<'a> {
    fn new(alloc: &'a mut TmpDirAllocator) -> Self {
        Self { alloc, chunk_count: 0, merge_count: 0 }
    }

    /// Returns the next unique chunk file path, drawing from the allocator.
    fn next_chunk_path(&mut self) -> Result<PathBuf> {
        let base = self.alloc.next()?;
        let path = base.join(format!("chunk_{:04}.keyed", self.chunk_count));
        self.chunk_count += 1;
        Ok(path)
    }

    /// Returns the next unique merged file path, drawing from the allocator.
    fn next_merged_path(&mut self) -> Result<PathBuf> {
        let base = self.alloc.next()?;
        let path = base.join(format!("merged_{:04}.keyed", self.merge_count));
        self.merge_count += 1;
        Ok(path)
    }
}

/// A spill write that is finishing in the background.
///
/// Used for pipelining: the I/O thread continues writing while the main thread
/// reads the next batch. The chunk path is stored alongside the handle so it
/// can be pushed to `chunk_files` only after the write completes.
struct PendingSpill {
    handle: crate::pooled_chunk_writer::SpillWriteHandle,
    chunk_path: PathBuf,
    /// Whether this chunk extended an existing run rather than starting one.
    /// An extended run is already in `chunk_files`, so it must not be added twice.
    appended: bool,
    /// The run file's length before this chunk was written, so the drain can
    /// account for the bytes this chunk added rather than the file's whole size.
    size_before: u64,
}

/// Decides whether the next sorted chunk extends the open spill run or starts a
/// new one — "natural run formation".
///
/// Phase 1 chunks its input sequentially, so an input already in the requested
/// order produces chunks that are not merely disjoint but *totally ordered*:
/// chunk 0 holds the earliest records, chunk 1 the next. Writing each to its own
/// file manufactures a k-way merge for a sequence that was never interleaved.
/// Appending instead collapses sorted input to a single run, and almost-sorted
/// input to its natural runs.
///
/// # Why only the most recent run
///
/// Considering older runs too (best fit) yields fewer runs on adversarial input
/// but breaks stability. Suppose run 0 ended at key 100 and stopped receiving,
/// run 1 is at 200, and a chunk arrives with minimum 150: best fit appends it to
/// run 0, so run 0 holds chunks {0,1,2,9} while run 1 holds {3,4}. Source index
/// order no longer matches ingest order, and for a key present in both chunk 3
/// and chunk 9 the merge would emit chunk 9's copy first.
///
/// Restricting appends to the most recent run keeps every run a *consecutive*
/// group of chunks, so source order equals ingest order. Combined with the loser
/// tree's tie-break on source index (`LoserTree::is_greater`), output is
/// byte-identical to writing one file per chunk.
#[derive(Default)]
struct RunFormer<K> {
    /// The open run: its path and the key of its last record. `None` before the
    /// first spill, and after a consolidation that swallowed the run.
    open: Option<(PathBuf, K)>,
}

/// Report how many spilled chunks became how many runs.
///
/// Sorted input collapses to one run, almost-sorted input to its natural runs,
/// shuffled input to one run per chunk. Logging the ratio means every sort
/// reports how ordered its own input was, which is the cheapest way to learn how
/// common near-sorted inputs are in practice.
///
/// `runs` must be the number of runs *formed* (`RawSortStats::runs_written`),
/// not the number of spill files left at merge time. Consolidation merges the
/// oldest half of the files whenever their count reaches `--max-temp-files`, so
/// the surviving file count can be far lower than the run count for reasons that
/// have nothing to do with how ordered the input was -- which is the one thing
/// this line exists to report.
fn log_run_formation(chunks_spilled: usize, runs: usize) {
    info!(
        "Spill runs: {runs} from {chunks_spilled} chunks ({} extended an existing run)",
        chunks_spilled.saturating_sub(runs)
    );
}

/// The smallest and largest keys in an already-sorted slice.
///
/// Read off the ends rather than scanned, which is what makes run formation free:
/// the caller has just sorted the buffer, so its extremes are its first and last
/// elements. `None` only for an empty slice, which a spill never produces.
fn chunk_bounds<T, K>(items: &[T], key: impl Fn(&T) -> K) -> Option<(K, K)> {
    Some((key(items.first()?), key(items.last()?)))
}

impl<K: Clone + Ord> RunFormer<K> {
    /// Choose where a sorted chunk should be written.
    ///
    /// Returns the path and whether it extends an existing run. The chunk extends
    /// the open run when every key in it sorts at or after that run's last key;
    /// otherwise a fresh path is allocated and a new run begins.
    ///
    /// Both pre-write steps in one method rather than two, because they have a
    /// mandatory order — forget a consolidated run, *then* test appendability —
    /// and nothing in the types enforced it when callers drove it themselves.
    /// Recording the new boundary necessarily follows the write and stays in
    /// [`Self::extended`], which [`RawExternalSorter::spill_chunk`] — the only
    /// caller of either — invokes once the path is known.
    fn place(
        &mut self,
        chunk_files: &[PathBuf],
        bounds: Option<&(K, K)>,
        namer: &mut ChunkNamer<'_>,
    ) -> Result<(PathBuf, bool)> {
        // Consolidation may have merged the open run away while the previous
        // spill was draining, leaving a path that no longer holds what we think.
        if let Some((path, _)) = &self.open
            && !chunk_files.contains(path)
        {
            self.open = None;
        }

        let appendable = match (&self.open, bounds) {
            (Some((path, last_key)), Some((chunk_min, _))) if last_key <= chunk_min => {
                Some(path.clone())
            }
            _ => None,
        };

        match appendable {
            Some(path) => Ok((path, true)),
            // Always allocate through the namer on this branch: it is what
            // advances the temp-directory round-robin and periodically re-checks
            // free space (`TmpDirAllocator::next`).
            None => Ok((namer.next_chunk_path()?, false)),
        }
    }

    /// Record that a chunk ending at `chunk_max` was written to `path`.
    fn extended(&mut self, path: PathBuf, chunk_max: K) {
        self.open = Some((path, chunk_max));
    }
}

/// Build `BufferProbeStats` from any buffer implementing `ProbeableBuffer`.
#[allow(clippy::cast_possible_truncation)]
fn probe_stats(buf: &impl ProbeableBuffer) -> BufferProbeStats {
    BufferProbeStats {
        usage: buf.memory_usage() as u64,
        capacity: buf.allocated_capacity() as u64,
        records: buf.len() as u64,
        segments: buf.num_segments() as u64,
    }
}

/// Raw-bytes external sorter for BAM files.
///
/// This sorter uses lazy record parsing to minimize memory usage and avoid
/// re-encoding overhead. It's significantly faster than the RecordBuf-based
/// sorter for large files.
pub struct RawExternalSorter {
    /// Sort order to use.
    sort_order: SortOrder,
    /// Maximum memory to use for in-memory sorting.
    memory_limit: usize,
    /// Temporary directories for spill files.
    ///
    /// When empty, a single directory is created under the system default
    /// temp location. When one or more paths are given, spill files are
    /// distributed across them in free-space-aware round-robin order via
    /// [`TmpDirAllocator`].
    temp_dirs: Vec<PathBuf>,
    /// Number of threads for parallel operations. Used by both phases unless
    /// overridden by `sort_threads` / `merge_threads`.
    threads: usize,
    /// Phase-1 (accumulate/sort/spill) worker count. `None` means "use
    /// `threads`". Lets a caller cede cores to an upstream producer during
    /// ingest while keeping the merge wide.
    sort_threads: Option<usize>,
    /// Phase-2 (merge/write) worker count. `None` means "use `threads`".
    merge_threads: Option<usize>,
    /// Compression level for output.
    output_compression: u32,
    /// Compression level for temporary chunk files (0 = uncompressed).
    temp_compression: u32,
    /// Codec used for temporary chunk files (Bgzf or Zstd).
    spill_codec: crate::codec::SpillCodec,
    /// Whether to write BAM index alongside output (coordinate sort only).
    write_index: bool,
    /// Program record info (version, `command_line`) for @PG header.
    pg_info: Option<(String, String)>,
    /// Maximum temp files before consolidation (0 = unlimited).
    max_temp_files: usize,
    /// Cell barcode tag for template-coordinate sort (e.g., `SamTag::CB`).
    /// When `Some`, CB hash is included in sort key for single-cell data.
    cell_tag: Option<SamTag>,
    /// Initial buffer capacity hint (bytes) for pre-allocation.
    ///
    /// Decoupled from `memory_limit` so that auto-detected limits can start with
    /// a modest allocation and let `Vec` grow on demand, while explicit limits
    /// pre-allocate the full budget upfront (preserving prior behavior).
    initial_capacity: Option<usize>,
    /// When true, wrap input in a `PrefetchReader` for async I/O.
    async_reader: bool,
    /// Which optional template-key lanes to retain (template-coordinate only).
    ///
    /// Defaults to [`KeyTypesSpec::Auto`], which provisions the narrowest key
    /// that fits the first record's optional lanes.
    key_types: KeyTypesSpec,
}

/// RAII guard that ensures Phase 2 teardown runs on every exit path between
/// `pool.set_phase(PHASE2)` and the end of the merge. Without this, any `?`
/// early-return would leave the pool stuck in PHASE2 with `phase2_files` still
/// published.
///
/// Merges finish by handing the guard to
/// [`finish_output`](Self::finish_output), which owns the whole teardown order —
/// the merge is done with its *sources* well before it is done with the *pool*,
/// and getting that order wrong is silent. `Drop` still deactivates, so a merge
/// that fails before it gets that far tears down completely.
struct Phase2Guard<'a, K: RawSortKey + 'static> {
    pool: &'a Arc<SortWorkerPool>,
    consumer: Option<MainThreadChunkConsumer<K>>,
    active: bool,
    sources_released: bool,
}

impl<K: RawSortKey + 'static> Phase2Guard<'_, K> {
    fn consumer_mut(&mut self) -> Option<&mut MainThreadChunkConsumer<K>> {
        self.consumer.as_mut()
    }

    fn consumer_ref(&self) -> Option<&MainThreadChunkConsumer<K>> {
        self.consumer.as_ref()
    }

    /// Release the drained merge sources, finalize the output, then leave Phase 2.
    ///
    /// Releasing the sources first is the point: it keeps a slow `finish` from
    /// holding the merge's file descriptors and per-file reorder buffers (a 2 MiB
    /// `BufReader` per spill file) alive while it drains. Leaving Phase 2 last
    /// then costs nothing, and keeps the pool's phase an accurate description of
    /// what the pool is doing for as long as it is doing it.
    ///
    /// Note what this ordering is *not* responsible for. The output's compression
    /// level does not depend on it: each block carries its own
    /// [`CompressTarget`](crate::worker_pool::CompressTarget), stamped by the
    /// writer that produced it, so trailing output blocks are compressed at
    /// `output_compression` whatever phase the pool is in when a worker pops
    /// them. That was not always true — the compressor used to be selected from
    /// the pool's phase at pop time, which is why finalizing after the return to
    /// `LEGACY` silently wrote the tail of the BAM at `temp_compression`.
    ///
    /// Takes `self` by value: finalizing the output is the guard's terminal
    /// operation, so a second call is a compile error rather than a runtime
    /// check. The `Drop` at the end of this method is a no-op, since
    /// `deactivate` has already cleared `active`.
    fn finish_output<T>(mut self, finish: impl FnOnce() -> Result<T>) -> Result<T> {
        self.release_sources();
        let finished = finish();
        self.deactivate();
        finished
    }

    /// Drop the consumer (releasing the Arc snapshot of the per-file vector) and
    /// unpublish the spill files, while leaving the pool in Phase 2.
    fn release_sources(&mut self) {
        if self.active && !self.sources_released {
            drop(self.consumer.take());
            self.pool.clear_phase2_files();
            self.sources_released = true;
        }
    }

    /// Leave Phase 2 entirely: release the sources if that has not happened yet,
    /// then return the pool to `LEGACY`.
    fn deactivate(&mut self) {
        if self.active {
            self.release_sources();
            self.pool.set_phase(crate::worker_pool::phase::LEGACY);
            self.active = false;
        }
    }
}

impl<K: RawSortKey + 'static> Drop for Phase2Guard<'_, K> {
    fn drop(&mut self) {
        self.deactivate();
    }
}

impl RawExternalSorter {
    /// Create a new raw external sorter with the given sort order.
    #[must_use]
    pub fn new(sort_order: SortOrder) -> Self {
        Self {
            sort_order,
            memory_limit: 512 * 1024 * 1024, // 512 MB default
            temp_dirs: Vec::new(),
            threads: 1,
            sort_threads: None,
            merge_threads: None,
            output_compression: 6,
            temp_compression: 1, // Default: fast compression
            spill_codec: crate::codec::SpillCodec::default(),
            write_index: false,
            pg_info: None,
            max_temp_files: crate::fd_limit::FALLBACK_MAX_TEMP_FILES,
            cell_tag: None,
            initial_capacity: None,
            async_reader: false,
            key_types: KeyTypesSpec::default(),
        }
    }

    /// Set the memory limit for in-memory sorting.
    #[must_use]
    pub fn memory_limit(mut self, limit: usize) -> Self {
        self.memory_limit = limit;
        self
    }

    /// Set a single temporary directory for spill files.
    ///
    /// Equivalent to calling [`Self::temp_dirs`] with a single-element vector.
    #[must_use]
    pub fn temp_dir(mut self, path: PathBuf) -> Self {
        self.temp_dirs = vec![path];
        self
    }

    /// Set multiple temporary directories for spill files.
    ///
    /// Spill files are distributed across the supplied directories in
    /// free-space-aware round-robin order. Passing an empty vector falls
    /// back to a single directory under the system temp location.
    #[must_use]
    pub fn temp_dirs(mut self, paths: Vec<PathBuf>) -> Self {
        self.temp_dirs = paths;
        self
    }

    /// Set the number of threads.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Set the Phase-1 (accumulate/sort/spill) worker count, overriding
    /// [`threads`](Self::threads) for that phase only.
    ///
    /// Phase 1 competes with whatever is feeding the sort, so a pipeline like
    /// `bwa mem -t 32 | fgumi sort` can lower this to cede cores to the producer
    /// during ingest while leaving the merge at full width.
    ///
    /// This is purely a scheduling knob: output is byte-identical regardless.
    #[must_use]
    pub fn sort_threads(mut self, n: usize) -> Self {
        self.sort_threads = Some(n);
        self
    }

    /// Set the Phase-2 (merge/write) worker count, overriding
    /// [`threads`](Self::threads) for that phase only.
    ///
    /// This is purely a scheduling knob: output is byte-identical regardless.
    #[must_use]
    pub fn merge_threads(mut self, n: usize) -> Self {
        self.merge_threads = Some(n);
        self
    }

    /// Effective Phase-1 worker count: `sort_threads` if set, else `threads`.
    #[must_use]
    pub fn phase1_threads(&self) -> usize {
        self.sort_threads.unwrap_or(self.threads).max(1)
    }

    /// Effective Phase-2 worker count: `merge_threads` if set, else `threads`.
    #[must_use]
    pub fn phase2_threads(&self) -> usize {
        self.merge_threads.unwrap_or(self.threads).max(1)
    }

    /// Transition the worker pool from ingest (Phase 1) to output (Phase 2).
    ///
    /// Every sort path calls this exactly once, immediately after ingest
    /// completes, to hand the pool over to Phase 2 and lift the active-worker cap
    /// to `merge_threads`. Centralized here so the transition cannot drift
    /// between sort orders.
    ///
    /// Workers idled by the Phase-1 cap re-check within `MAX_BACKOFF_US`, so no
    /// explicit wake is needed; capped workers always drained their held items, so
    /// raising the cap can never be required for correctness, only for width.
    ///
    /// Omitting this call once left in-memory sorts writing their output at
    /// `temp_compression`, because the phase used to decide which compressor a
    /// block went through. It no longer does — blocks carry their own
    /// [`CompressTarget`](crate::worker_pool::CompressTarget) — so what the
    /// phase buys now is Phase 2 scheduling and the wider worker cap. See
    /// [`SortWorkerPool::begin_phase2`].
    fn enter_output_phase(&self, pool: &SortWorkerPool) {
        pool.begin_phase2(self.phase2_threads());
    }

    /// Configured maximum number of temporary spill files kept before the
    /// oldest runs are consolidated into one (`0` means unlimited). Exposed so
    /// callers can assert how their `max_temp_files` setting resolved.
    #[must_use]
    pub fn temp_file_limit(&self) -> usize {
        self.max_temp_files
    }

    /// Set the output compression level.
    #[must_use]
    pub fn output_compression(mut self, level: u32) -> Self {
        self.output_compression = level;
        self
    }

    /// Set compression level for temporary chunk files.
    ///
    /// For [`SpillCodec::Bgzf`], level 0 disables compression (fastest, uses
    /// most disk space). For [`SpillCodec::Zstd`] there is no level-0 "stored"
    /// mode; pass a level >= 1 or [`Self::sort`] will reject the combination.
    /// Level 1 (default) provides fast compression with reasonable space savings.
    /// Higher levels provide better compression but are slower.
    ///
    /// [`SpillCodec::Bgzf`]: crate::codec::SpillCodec::Bgzf
    /// [`SpillCodec::Zstd`]: crate::codec::SpillCodec::Zstd
    #[must_use]
    pub fn temp_compression(mut self, level: u32) -> Self {
        self.temp_compression = level;
        self
    }

    /// Set the codec used for temporary chunk files.
    ///
    /// Defaults to [`SpillCodec::Zstd`](crate::codec::SpillCodec::Zstd) which is significantly faster than
    /// BGZF at comparable ratios for BAM-record data.
    #[must_use]
    pub fn spill_codec(mut self, codec: crate::codec::SpillCodec) -> Self {
        self.spill_codec = codec;
        self
    }

    /// Enable writing BAM index alongside output.
    ///
    /// Only valid for coordinate sort. When enabled, writes `<output>.bai`
    /// alongside the output BAM file. Output BGZF compression stays
    /// multi-threaded (scales with the configured thread count); BAI virtual
    /// offsets are recovered from each BGZF block as it finalizes, so indexing
    /// does not serialize compression.
    #[must_use]
    pub fn write_index(mut self, enabled: bool) -> Self {
        self.write_index = enabled;
        self
    }

    /// Set program record info for @PG header entry.
    #[must_use]
    pub fn pg_info(mut self, version: String, command_line: String) -> Self {
        self.pg_info = Some((version, command_line));
        self
    }

    /// Set maximum temp files before consolidation.
    ///
    /// When the number of temp files exceeds this limit, the oldest files
    /// are merged together to reduce the count. Set to 0 for unlimited.
    ///
    /// The default is a fixed, portable value. To size it to the host's
    /// descriptor budget instead — which is what `fgumi sort` does — pass
    /// [`resolve_temp_file_limit`](crate::resolve_temp_file_limit).
    #[must_use]
    pub fn max_temp_files(mut self, max: usize) -> Self {
        self.max_temp_files = max;
        self
    }

    /// Set the cell barcode tag for template-coordinate sort.
    ///
    /// When set, the CB hash is included in the sort key so that templates
    /// from different cells at the same locus are not interleaved.
    #[must_use]
    pub fn cell_tag(mut self, tag: SamTag) -> Self {
        self.cell_tag = Some(tag);
        self
    }

    /// Set the initial buffer capacity hint (bytes).
    ///
    /// When set, buffer pre-allocation uses this value instead of `memory_limit`.
    /// This avoids huge upfront allocations when auto-detecting memory, while
    /// still allowing the buffer to grow up to `memory_limit` before spilling.
    #[must_use]
    pub fn initial_capacity(mut self, bytes: usize) -> Self {
        self.initial_capacity = Some(bytes);
        self
    }

    /// Set which optional template-key lanes to retain.
    ///
    /// Only affects the template-coordinate sort order. The default
    /// ([`KeyTypesSpec::Auto`]) provisions the narrowest key consistent with the
    /// first record's optional lanes; the decode-time verify guarantees the
    /// dropped lanes are constant across all records.
    #[must_use]
    pub fn key_types(mut self, spec: KeyTypesSpec) -> Self {
        self.key_types = spec;
        self
    }

    /// Enable/disable the async prefetch reader on input.
    ///
    /// When enabled, the input BAM is wrapped in a `PrefetchReader` before the
    /// BGZF layer, which overlaps block I/O with decompression.
    #[must_use]
    pub fn async_reader(mut self, enabled: bool) -> Self {
        self.async_reader = enabled;
        self
    }

    /// Returns the effective initial capacity for buffer pre-allocation.
    ///
    /// Uses `initial_capacity` if set, otherwise falls back to `memory_limit`.
    fn effective_initial_capacity(&self) -> usize {
        self.initial_capacity.unwrap_or(self.memory_limit).min(self.memory_limit)
    }

    /// Build a rayon thread pool sized to `self.threads`.
    ///
    /// The sort path uses `par_sort` and friends at several points. Rayon's
    /// global pool defaults to `num_cpus::get()`, which silently violates the
    /// user's `--threads` contract on machines where more physical cores are
    /// available. Every rayon call site is wrapped with `pool.install(...)`
    /// so that `rayon::current_num_threads()` returns `self.threads` and
    /// fan-out is bounded to the requested thread count.
    ///
    /// Oversubscription with the `SortWorkerPool` is not a concern because
    /// every call site is preceded by `drain_pending_spill`, which joins the
    /// prior chunk's I/O thread and therefore guarantees all sort workers are
    /// idle at the moment rayon fans out.
    fn build_sort_rayon_pool(&self) -> Result<rayon::ThreadPool> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.phase1_threads())
            .thread_name(|i| format!("fgumi-sort-rayon-{i}"))
            .build()
            .map_err(|e| anyhow::anyhow!("failed to build rayon sort pool: {e}"))
    }

    /// Write one sorted chunk, extending the open spill run when it can.
    ///
    /// Shared by all four sort functions: identical aside from the key type and
    /// the per-record write, which the caller supplies as a closure. Keeping the
    /// run-formation bookkeeping here means the append/create choice, the
    /// pre-write size needed for spill accounting, and the boundary key recorded
    /// afterwards cannot drift between the four copies.
    fn spill_chunk<K: RawSortKey + Clone + Ord + Default + 'static>(
        run_former: &mut RunFormer<K>,
        chunk_files: &[PathBuf],
        bounds: Option<(K, K)>,
        namer: &mut ChunkNamer<'_>,
        pool: &Arc<SortWorkerPool>,
        timer: &mut SortPhaseTimer,
        write: impl FnOnce(&mut PooledChunkWriter<K>) -> Result<()>,
    ) -> Result<PendingSpill> {
        let (chunk_path, appended) = run_former.place(chunk_files, bounds.as_ref(), namer)?;

        // Captured before the write so the drain can charge only the bytes this
        // chunk adds; an appended run's file already holds its predecessors.
        let size_before =
            if appended { std::fs::metadata(&chunk_path).map_or(0, |m| m.len()) } else { 0 };

        let handle = timer.time_spill_write(|| {
            let mut writer = PooledChunkWriter::<K>::open(
                Arc::clone(pool),
                &chunk_path,
                pool.spill_codec(),
                appended,
            )?;
            write(&mut writer)?;
            writer.start_finish()
        })?;

        if let Some((_, chunk_max)) = bounds {
            run_former.extended(chunk_path.clone(), chunk_max);
        }
        Ok(PendingSpill { handle, chunk_path, appended, size_before })
    }

    /// Consolidate temp files if we've exceeded the limit.
    /// Wait for a pending spill to complete and, if one was present, run consolidation.
    ///
    /// Shared by all four sort functions: the in-loop "wait before next spill" drain and
    /// the post-loop "drain before merge" drain are identical aside from the key type.
    fn drain_pending_spill<K: RawSortKey + Default + 'static>(
        &self,
        pending: &mut Option<PendingSpill>,
        chunk_files: &mut Vec<PathBuf>,
        stats: &mut RawSortStats,
        timer: &mut SortPhaseTimer,
        namer: &mut ChunkNamer<'_>,
        pool: &std::sync::Arc<crate::worker_pool::SortWorkerPool>,
    ) -> Result<()> {
        if let Some(prev) = pending.take() {
            prev.handle.wait()?;
            timer.record_spill_growth(&prev.chunk_path, prev.size_before);
            // A chunk that extended an existing run is already represented in
            // `chunk_files`; only a chunk that started a run adds a merge source.
            if !prev.appended {
                chunk_files.push(prev.chunk_path);
                stats.runs_written += 1;
            }

            timer.time_consolidate(|| {
                self.maybe_consolidate_temp_files::<K>(chunk_files, namer, pool)
            })?;
        }
        Ok(())
    }

    ///
    /// Merges the oldest half of temp files into a single new file to reduce
    /// the total count while maintaining sort order.
    fn maybe_consolidate_temp_files<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &mut Vec<PathBuf>,
        namer: &mut ChunkNamer<'_>,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<()> {
        use crate::loser_tree::LoserTree;

        if self.max_temp_files == 0 || chunk_files.len() < self.max_temp_files {
            return Ok(());
        }

        // Need at least 2 files to consolidate meaningfully
        if self.max_temp_files < 2 {
            return Ok(());
        }

        // Merge oldest half of files into one (at least 2)
        let n_to_merge = (self.max_temp_files / 2).max(2).min(chunk_files.len());
        let files_to_merge: Vec<PathBuf> = chunk_files.drain(..n_to_merge).collect();

        debug!(
            "Consolidating {} temp files into 1 (total was {})...",
            n_to_merge,
            n_to_merge + chunk_files.len()
        );

        // Create merged output file
        let merged_path = namer.next_merged_path()?;

        // Open readers with semaphore to cap concurrent I/O.
        let sem = make_reader_semaphore(self.phase1_threads());
        let mut readers: Vec<GenericKeyedChunkReader<K>> = files_to_merge
            .iter()
            .map(|p| GenericKeyedChunkReader::<K>::open(p, Some(Arc::clone(&sem))))
            .collect::<Result<Vec<_>>>()?;

        // Use pooled writer for parallel compression during consolidation.
        // Match the pool's spill codec so the merged file is the same format
        // as the per-chunk spill files it consolidates.
        let mut writer =
            PooledChunkWriter::<K>::new(Arc::clone(pool), &merged_path, pool.spill_codec())?;

        // Initialize loser tree with first record from each reader
        let mut initial_keys: Vec<K> = Vec::with_capacity(readers.len());
        let mut records: Vec<Vec<u8>> = Vec::with_capacity(readers.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(readers.len());

        for (reader_idx, reader) in readers.iter_mut().enumerate() {
            let mut record = Vec::new();
            if let Some(key) = reader.next_record(&mut record)? {
                initial_keys.push(key);
                records.push(record);
                source_map.push(reader_idx);
            }
        }

        if initial_keys.is_empty() {
            writer.finish()?;
            // Insert at beginning to preserve stable order
            chunk_files.insert(0, merged_path);
            // Clean up old files
            for path in &files_to_merge {
                let _ = std::fs::remove_file(path);
            }
            return Ok(());
        }

        let mut tree = LoserTree::new(initial_keys);

        while tree.winner_is_active() {
            let winner = tree.winner();
            let reader_idx = source_map[winner];
            writer.write_record(tree.winner_key(), &records[winner])?;

            if let Some(next_key) = readers[reader_idx].next_record(&mut records[winner])? {
                tree.replace_winner(next_key);
            } else {
                tree.remove_winner();
            }
        }

        writer.finish()?;

        // Insert merged file at the beginning to preserve stable order for equal keys.
        // The merged file contains the oldest records, so it should be processed first.
        chunk_files.insert(0, merged_path);

        // Clean up old files
        for path in &files_to_merge {
            let _ = std::fs::remove_file(path);
        }

        debug!("Consolidation complete, {} temp files remain", chunk_files.len());

        Ok(())
    }

    /// Sort a BAM file using raw-bytes approach.
    ///
    /// # Errors
    ///
    /// Returns an error if reading, sorting, or writing the BAM file fails.
    pub fn sort(&self, input: &Path, output: &Path) -> Result<RawSortStats> {
        let pool = self.create_worker_pool()?;

        // Open input BAM and create record source
        // N+2 model: workers do ReadInputBlocks + DecompressInput,
        // main thread reads records directly from PooledInputStream.
        debug!(
            "Phase 1: Pool-integrated input reading ({} workers, N+2 model)",
            pool.num_workers()
        );
        let (record_source, header) = {
            let (reader, header) =
                create_raw_bam_reader_pool_integrated(input, &pool, self.async_reader)?;
            (RecordSource::direct(reader), header)
        };

        self.sort_from_source(record_source, &header, output, pool)
    }

    /// Sort records produced in-process, without a BAM input file.
    ///
    /// This is the streaming counterpart to [`Self::sort`]: instead of reading
    /// an input BAM, records are pulled from `records` as the sort's accumulate
    /// phase consumes them. Producers that would otherwise write an unsorted BAM
    /// only to hand it straight back for sorting should use this — it removes a
    /// full compress/write/read/decompress round-trip, and removes the
    /// intermediate file's disk footprint entirely (which for inputs that exceed
    /// `memory_limit` would otherwise sit alongside the spill chunks).
    ///
    /// Spill and merge behave exactly as in [`Self::sort`]: when everything fits
    /// within `memory_limit` nothing touches disk except the output; past that
    /// limit, sorted chunks spill to the configured temp directories and are
    /// k-way merged.
    ///
    /// `header` is the input header — the output header's sort-order tags are
    /// derived from it the same way [`Self::sort`] derives them from the input
    /// BAM's header, and `pg_info` (if set) is applied to it.
    ///
    /// An `Err` item from `records` aborts the sort and is propagated; the
    /// output file is not written, so a producer failure can never leave a
    /// silently truncated result behind.
    ///
    /// # Errors
    ///
    /// Returns an error if the producer yields one, or if sorting, spilling, or
    /// writing the output BAM fails.
    pub fn sort_records<I>(
        &self,
        records: I,
        header: &Header,
        output: &Path,
    ) -> Result<RawSortStats>
    where
        I: IntoIterator<Item = Result<fgumi_raw_bam::RawRecord>>,
        I::IntoIter: Send + 'static,
    {
        let pool = self.create_worker_pool()?;
        debug!("Phase 1: In-process record stream (no input BAM)");
        self.sort_from_source(RecordSource::stream(records), header, output, pool)
    }

    /// Build the shared worker pool used by every phase of a sort run.
    ///
    /// Also enforces the spill-codec/compression-level invariant, so no entry
    /// point can reach a misconfigured pool.
    fn create_worker_pool(&self) -> Result<Arc<SortWorkerPool>> {
        // zstd has no level-0 "stored" mode; silently remapping to 1 would
        // surprise API callers who set temp_compression(0) expecting an
        // uncompressed spill (which works for BGZF). Mirror the CLI guard so
        // direct RawExternalSorter callers cannot reach a misconfigured pool.
        anyhow::ensure!(
            !(self.temp_compression == 0
                && matches!(self.spill_codec, crate::codec::SpillCodec::Zstd)),
            "temp_compression=0 is only supported with SpillCodec::Bgzf; \
             zstd does not have an uncompressed mode. Pass \
             spill_codec(SpillCodec::Bgzf) to keep level-0 spill, or pick a \
             zstd level >= 1."
        );

        debug!("Starting raw-bytes sort with order: {:?}", self.sort_order);
        debug!("Memory limit: {} MB", self.memory_limit / (1024 * 1024));
        debug!(
            "Threads: {} (default for both phases; phase 1: {}, phase 2: {})",
            self.threads,
            self.phase1_threads(),
            self.phase2_threads()
        );

        // One worker pool spans both phases, so size it to the wider of the two
        // and cap the active count per phase via `set_active_workers`.
        let pool = SortWorkerPool::new(
            self.phase1_threads().max(self.phase2_threads()),
            self.temp_compression,
            self.output_compression,
            self.spill_codec,
        );
        // Phase 1 runs first; the merge raises this to `phase2_threads()`.
        pool.set_active_workers(self.phase1_threads());
        Ok(Arc::new(pool))
    }

    /// Run the sort over an already-constructed record source.
    ///
    /// Shared tail of [`Self::sort`] and [`Self::sort_records`]: applies
    /// `pg_info` to the header, provisions temp directories, and dispatches to
    /// the per-order implementation.
    fn sort_from_source(
        &self,
        record_source: RecordSource,
        header: &Header,
        output: &Path,
        pool: Arc<SortWorkerPool>,
    ) -> Result<RawSortStats> {
        // Say so when the input claims to be in the order being requested. This
        // is checked before @PG is added so it reflects the input as given.
        if crate::header_declares_order(header, self.sort_order) {
            log::warn!(
                "Input header already declares this sort order ({}); sorting it again. \
                 Headers are not verified, so the sort still runs -- use `--verify` to \
                 check an input's order without rewriting it.",
                self.sort_order_flag_value()
            );
        }

        // Add @PG record if pg_info was provided
        let header = if let Some((ref version, ref command_line)) = self.pg_info {
            fgumi_bam_io::header::add_pg_record(header.clone(), version, command_line)?
        } else {
            header.clone()
        };

        // _temp_dirs: RAII handles; kept alive until sort returns.
        let (_temp_dirs, mut alloc) = self.create_temp_dirs()?;

        // Sort based on order
        let stats = match self.sort_order {
            SortOrder::Coordinate => {
                self.sort_coordinate(record_source, pool, &header, output, &mut alloc)
            }
            SortOrder::Queryname(comparator) => {
                self.sort_queryname(record_source, pool, &header, output, &mut alloc, comparator)
            }
            SortOrder::TemplateCoordinate => {
                self.sort_template_coordinate(record_source, pool, &header, output, &mut alloc)
            }
        }?;

        // A sort is a permutation: every record read must be written. Checked here
        // rather than in the CLI so it covers every entry point -- `sort_records`
        // has callers (`fgumi simulate`) that discard the stats entirely, and a
        // library consumer should not have to opt in to the engine's own
        // guarantee. Both counts are measured independently: `total_records` by
        // the ingest loop, `output_records` by whichever writer produced the
        // output.
        enforce_record_count(&stats, output, self.write_index)?;
        Ok(stats)
    }

    /// Warn when a merge is about to ask for more descriptors than the process
    /// has.
    ///
    /// A merge opens every input at once, so its descriptor cost is the input
    /// count — user-controlled, and not bounded by `max_temp_files`, which the
    /// merge path never consults. Without this the failure surfaces as a bare
    /// "Too many open files" naming one path, which reads like a problem with
    /// that file rather than with the size of the merge.
    ///
    /// Warns rather than refuses: the reserve is deliberately conservative, so a
    /// merge slightly over the estimate may still succeed, and refusing it would
    /// turn a working command into a broken one. The subsequent open fails
    /// immediately anyway — every reader is opened before any merging starts —
    /// so nothing long-running is wasted either way.
    ///
    /// Returns the message rather than logging it so the threshold and the
    /// wording are assertable without a log capture; `merge_bams` emits it.
    fn fd_budget_warning(num_inputs: usize) -> Option<String> {
        // One reading of the budget answers both questions in the helper.
        // Reading it per question would let the warning name a budget other
        // than the one that tripped it.
        Self::fd_budget_warning_from_nofile(crate::fd_limit::soft_nofile(), num_inputs)
    }

    /// [`Self::fd_budget_warning`] against a supplied budget rather than the
    /// process's own.
    ///
    /// Split out for the reason `fd_limit::temp_file_limit_from_soft_nofile`
    /// is: the arithmetic is then testable without the host's real
    /// `RLIMIT_NOFILE` deciding the answer. Reading it inside made both tests
    /// host-dependent — the silent case needs a soft limit of at least
    /// `8 + FD_RESERVE`, so it failed under a `ulimit -n` below 40, and the
    /// warning case only reached the branch because `usize::MAX` saturates the
    /// comparison on a 64-bit target. The CLI's own `fd_budget_warning` already
    /// takes the budget as a parameter; this makes the engine's match.
    fn fd_budget_warning_from_nofile(soft: Option<u64>, num_inputs: usize) -> Option<String> {
        if crate::fd_limit::fits_nofile_budget(soft, num_inputs) {
            return None;
        }
        let budget = soft.map_or_else(
            || "unknown".to_string(),
            |soft| soft.saturating_sub(crate::fd_limit::FD_RESERVE).to_string(),
        );
        Some(format!(
            "Merging {num_inputs} inputs opens {num_inputs} files at once, which exceeds this \
             process's open file budget of about {budget}. Raise it with `ulimit -n`, or merge in \
             stages, if this fails with \"Too many open files\"."
        ))
    }

    /// Merge multiple pre-sorted BAM files into a single sorted BAM.
    ///
    /// Each input BAM must already be sorted in the order specified by
    /// `self.sort_order`. The output preserves the sort order.
    ///
    /// # Errors
    ///
    /// Returns an error if any input cannot be opened, or writing fails.
    pub fn merge_bams(&self, inputs: &[PathBuf], header: &Header, output: &Path) -> Result<u64> {
        use crate::inline::extract_coordinate_key_inline;
        use crate::keys::{
            QuerynameComparator, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
            SortContext,
        };

        debug!("Starting k-way merge of {} BAM files", inputs.len());
        if let Some(warning) = Self::fd_budget_warning(inputs.len()) {
            warn!("{warning}");
        }

        let mut readers = Self::open_bam_prefetch_readers(inputs)?;
        let output_header = self.create_output_header(header);

        match self.sort_order {
            SortOrder::TemplateCoordinate => {
                let lib_lookup = LibraryLookup::from_header(header);
                let cell_tag = self.cell_tag;
                let hasher = cb_hasher();
                self.run_merge_loop(&mut readers, inputs, &output_header, output, |bam| {
                    extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher)
                })
            }
            SortOrder::Coordinate => {
                #[allow(clippy::cast_possible_truncation)]
                let nref = header.reference_sequences().len() as u32;
                self.run_merge_loop(&mut readers, inputs, &output_header, output, |bam| {
                    RawCoordinateKey { sort_key: extract_coordinate_key_inline(bam, nref) }
                })
            }
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                let ctx = SortContext::from_header(header);
                self.run_merge_loop(&mut readers, inputs, &output_header, output, |bam| {
                    RawQuerynameLexKey::extract(bam, &ctx)
                })
            }
            SortOrder::Queryname(QuerynameComparator::Natural) => {
                let ctx = SortContext::from_header(header);
                self.run_merge_loop(&mut readers, inputs, &output_header, output, |bam| {
                    RawQuerynameKey::extract(bam, &ctx)
                })
            }
        }
    }

    /// The `fgumi sort --order` value for this sorter's order, for error hints.
    fn sort_order_flag_value(&self) -> &'static str {
        match self.sort_order {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => "queryname",
            SortOrder::Queryname(QuerynameComparator::Natural) => "queryname::natural",
            SortOrder::TemplateCoordinate => "template-coordinate",
        }
    }

    /// Open background prefetch readers for multiple BAM files.
    fn open_bam_prefetch_readers(inputs: &[PathBuf]) -> Result<Vec<RawReadAheadReader>> {
        inputs
            .iter()
            .map(|path| {
                let (reader, _header) = create_raw_bam_reader(path, 1)?;
                Ok(RawReadAheadReader::new(reader))
            })
            .collect()
    }

    /// K-way merge loop: extract keys on the merge thread, write to output.
    ///
    /// Uses reusable per-source record buffers to avoid per-record heap
    /// allocations during the merge.
    fn run_merge_loop<K: Ord>(
        &self,
        readers: &mut [RawReadAheadReader],
        inputs: &[PathBuf],
        output_header: &Header,
        output: &Path,
        extract_key: impl Fn(&[u8]) -> K,
    ) -> Result<u64> {
        use crate::loser_tree::LoserTree;

        // Initialize: collect first record + key from each reader using
        // reusable per-source buffers (one allocation per source, not per record).
        let mut initial_keys: Vec<K> = Vec::with_capacity(readers.len());
        let mut records: Vec<Vec<u8>> = Vec::with_capacity(readers.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(readers.len());

        for (idx, reader) in readers.iter_mut().enumerate() {
            if let Some(raw_record) = reader.next() {
                let mut buf = Vec::with_capacity(raw_record.as_ref().len());
                buf.extend_from_slice(raw_record.as_ref());
                initial_keys.push(extract_key(&buf));
                records.push(buf);
                source_map.push(idx);
            }
        }

        // Write to an exclusive temp and atomically persist on success, so a
        // mid-merge sort-order violation (below) — or any other error — never
        // leaves a partial/corrupt output or a stray temp behind. Chosen before
        // the empty-input branch so header-only outputs are finalized through the
        // same atomic path as non-empty merges.
        let out_target = MergeOutputTarget::create(output)?;

        if initial_keys.is_empty() {
            debug!("Merge complete: 0 records merged");
            let writer = fgumi_bam_io::create_raw_bam_writer(
                out_target.path(),
                output_header,
                self.threads,
                self.output_compression,
            )?;
            writer.finish()?;
            out_target.persist()?;
            return Ok(0);
        }

        let mut tree = LoserTree::new(initial_keys);

        let mut writer = fgumi_bam_io::create_raw_bam_writer(
            out_target.path(),
            output_header,
            self.threads,
            self.output_compression,
        )?;

        let mut records_merged = 0u64;
        let merge_progress = ProgressTracker::new("Merged records").with_interval(1_000_000);
        // Set when an input is found not to be monotonic in the merge order:
        // (source index into `inputs`, 1-based record number within that input).
        let mut violation: Option<(usize, u64)> = None;
        // Per-source count of records pulled so far (each source starts with its
        // first record already loaded), used to report the offending record's
        // 1-based position within its input.
        let mut pulled_per_source = vec![1u64; records.len()];

        while tree.winner_is_active() {
            let winner = tree.winner();

            writer.write_raw_record(&records[winner])?;
            records_merged += 1;
            merge_progress.log_if_needed(1);

            let reader_idx = source_map[winner];
            if let Some(raw_record) = readers[reader_idx].next() {
                let buf = &mut records[winner];
                buf.clear();
                buf.extend_from_slice(raw_record.as_ref());
                let new_key = extract_key(buf);
                // MERGE3-01 streaming verify: the merge assumes each input is
                // already sorted in the merge order. `tree.winner_key()` still
                // holds the key of the record we just emitted from this source;
                // if the next record sorts before it, the input isn't monotonic
                // and the k-way merge would silently corrupt the output.
                if new_key < *tree.winner_key() {
                    violation = Some((reader_idx, pulled_per_source[winner] + 1));
                    break;
                }
                pulled_per_source[winner] += 1;
                tree.replace_winner(new_key);
            } else {
                tree.remove_winner();
            }
        }

        if let Some((reader_idx, record_no)) = violation {
            // Close the writer's handle to the partial output. `out_target` is
            // still owned here, so returning the error below drops it and removes
            // the temp — nothing partial is left behind. Not finalizing a doomed
            // output also avoids `finish` masking the real (sort-order) error.
            drop(writer);
            let source = inputs[reader_idx].display();
            let order = self.sort_order_flag_value();
            // The path is named for identification only; the remediation commands
            // use a `<input>` placeholder so a path with shell metacharacters
            // can't turn the copy-pasteable hint into something unexpected.
            anyhow::bail!(
                "Input '{source}' is not sorted in {order} order: record {record_no} sorts \
                 before a preceding record from the same input, so the k-way merge would \
                 corrupt the output. Sort that input in {order} order first \
                 (`fgumi sort -i <input> -o sorted.bam --order {order}`), or locate the \
                 violation with `fgumi sort -i <input> --verify --order {order}`."
            );
        }

        writer.finish()?;
        out_target.persist()?;
        merge_progress.log_final();

        Ok(records_merged)
    }

    /// Sort by coordinate order using optimized radix sort for large arrays.
    fn sort_coordinate(
        &self,
        record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        if self.write_index {
            self.sort_coordinate_with_index(record_source, pool, header, output, alloc)
        } else {
            self.sort_coordinate_optimized(record_source, pool, header, output, alloc)
        }
    }

    /// Optimized coordinate sort using inline buffer for reduced memory overhead.
    ///
    /// Uses `RecordBuffer` which stores records in a single contiguous allocation
    /// with pre-computed sort keys, eliminating per-record heap allocations.
    #[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
    fn sort_coordinate_optimized(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        use crate::keys::RawCoordinateKey;

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        // Get number of references (unmapped reads map to nref)
        let nref = header.reference_sequences().len() as u32;

        // Estimate capacity from initial_capacity (not memory_limit) to avoid
        // huge upfront allocations when auto-detecting memory.
        let init_cap = self.effective_initial_capacity();
        // Per-record footprint: ~200 bytes BAM + 8 header + 24 ref ≈ 232 bytes (rounded to 240 for headroom)
        let estimated_records = init_cap / 240;
        // Data bytes = init_cap minus ref overhead (24 bytes/record)
        let estimated_data_bytes = init_cap.saturating_sub(estimated_records * 24);

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, estimated_data_bytes, nref);
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        // Natural run formation: consecutive chunks that are already in order
        // extend one run instead of each becoming its own merge source.
        let mut run_former: RunFormer<crate::keys::RawCoordinateKey> = RunFormer::default();
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        debug!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        // Borrow each record's bytes straight out of the decompressed block and
        // push them into the arena, skipping the intermediate `RawRecord` copy
        // (the borrowed slice is invalidated by the next `next_record_borrowed`
        // call, which is fine — `push_coordinate` copies the bytes into the buffer).
        while let Some(record) = record_source.next_record_borrowed()? {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Push directly to buffer - key extracted inline from raw bytes
            buffer.push_coordinate(record)?;

            if probe.should_sample_read(stats.total_records) {
                probe.log_mid_read(probe_stats(&buffer), Some(pool.phase1_queue_depths()));
            }

            // Check memory usage
            if buffer.memory_usage() >= self.memory_limit {
                timer.end_read_span();
                let bstats = probe_stats(&buffer);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete before starting a new one
                self.drain_pending_spill::<RawCoordinateKey>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                // The buffer is sorted, so the chunk's extreme keys are its ends.
                let bounds =
                    chunk_bounds(buffer.refs(), |r| RawCoordinateKey { sort_key: r.sort_key });

                // Pipelined: `spill_chunk` returns once the write is submitted, so
                // I/O continues in the background while we read the next batch.
                pending_spill = Some(Self::spill_chunk::<RawCoordinateKey>(
                    &mut run_former,
                    &chunk_files,
                    bounds,
                    &mut namer,
                    &pool,
                    &mut timer,
                    |writer| {
                        for r in buffer.refs() {
                            let key = RawCoordinateKey { sort_key: r.sort_key };
                            writer.write_record(&key, buffer.get_record(r))?;
                        }
                        Ok(())
                    },
                )?);

                buffer.clear();
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        progress.log_final();
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<RawCoordinateKey>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(buffer.memory_usage() as u64);

        // Ingest is done: everything from here is Phase 2 (merge/write).
        self.enter_output_phase(&pool);

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            stats.output_records = timer.time_write_output(|| {
                use crate::pooled_bam_writer::PooledBamWriter;
                let output_header = self.create_output_header(header);
                let mut writer = PooledBamWriter::new(Arc::clone(&pool), output, &output_header)?;

                let mut written: u64 = 0;
                for record_bytes in buffer.iter_sorted() {
                    writer.write_raw_record(record_bytes)?;
                    written += 1;
                }
                writer.finish()?;
                Ok(written)
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks (avoids
            // intermediate merge back into a single sorted buffer); each
            // chunk becomes its own in-memory merge source.
            let memory_chunks: Vec<InMemoryChunk<RawCoordinateKey>> = if buffer.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort_into_chunks(self.phase1_threads()))
                })
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                vec![buffer.drain_into_single_chunk()]
            };

            let memory_chunks = MemorySources::Shared(memory_chunks);
            let n_memory = memory_chunks.num_non_empty();
            log_run_formation(timer.spill_count(), stats.runs_written);
            debug!(
                "Phase 2: Merging {} chunks (keyed O(1) comparisons)...",
                chunk_files.len() + n_memory
            );

            // Merge disk chunks + in-memory chunks using O(1) key comparisons
            stats.output_records = timer.time_merge(|| {
                self.merge_chunks_generic::<RawCoordinateKey>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )
            })?;
        }

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.phase1_threads(), self.phase2_threads(), self.max_temp_files);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Coordinate sort with BAM index generation.
    ///
    /// Similar to `sort_coordinate_optimized` but uses `IndexingBamWriter` to
    /// build the BAI index incrementally during write. Output BGZF compression
    /// stays multi-threaded (scales with the configured thread count); BAI
    /// virtual offsets are recovered from each BGZF block as it finalizes, so
    /// indexing does not serialize compression.
    #[allow(clippy::cast_possible_truncation, clippy::too_many_lines)]
    fn sort_coordinate_with_index(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        use crate::keys::RawCoordinateKey;
        use fgumi_bam_io::{bai_sidecar_path, create_indexing_bam_writer, write_bai_index};

        debug!("Indexing enabled: will write BAM index alongside output");

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        let nref = header.reference_sequences().len() as u32;
        let init_cap = self.effective_initial_capacity();
        // Per-record footprint: ~200 bytes BAM + 8 header + 24 ref ≈ 232 bytes (rounded to 240 for headroom)
        let estimated_records = init_cap / 240;
        let estimated_data_bytes = init_cap.saturating_sub(estimated_records * 24);

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer = RecordBuffer::with_capacity(estimated_records, estimated_data_bytes, nref);
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        // Natural run formation: consecutive chunks that are already in order
        // extend one run instead of each becoming its own merge source.
        let mut run_former: RunFormer<crate::keys::RawCoordinateKey> = RunFormer::default();
        let rayon_pool = self.build_sort_rayon_pool()?;

        debug!("Phase 1: Reading and sorting chunks (inline buffer, keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            buffer.push_coordinate(record.as_ref())?;

            if probe.should_sample_read(stats.total_records) {
                probe.log_mid_read(probe_stats(&buffer), Some(pool.phase1_queue_depths()));
            }

            if buffer.memory_usage() >= self.memory_limit {
                timer.end_read_span();
                let bstats = probe_stats(&buffer);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete
                self.drain_pending_spill::<RawCoordinateKey>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                let bounds =
                    chunk_bounds(buffer.refs(), |r| RawCoordinateKey { sort_key: r.sort_key });

                // Pipelined: `spill_chunk` returns once the write is submitted, so
                // I/O continues in the background while we read the next batch.
                pending_spill = Some(Self::spill_chunk::<RawCoordinateKey>(
                    &mut run_former,
                    &chunk_files,
                    bounds,
                    &mut namer,
                    &pool,
                    &mut timer,
                    |writer| {
                        for r in buffer.refs() {
                            let key = RawCoordinateKey { sort_key: r.sort_key };
                            writer.write_record(&key, buffer.get_record(r))?;
                        }
                        Ok(())
                    },
                )?);

                buffer.clear();
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        debug!("Read {} records total", stats.total_records);
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<RawCoordinateKey>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(buffer.memory_usage() as u64);

        let output_header = self.create_output_header(header);

        // Ingest is done: everything from here is Phase 2 (merge/write).
        self.enter_output_phase(&pool);

        if chunk_files.is_empty() {
            // All records fit in memory - no merge needed
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            stats.output_records = timer.time_write_output(|| {
                let mut writer = create_indexing_bam_writer(
                    output,
                    &output_header,
                    self.output_compression,
                    self.phase2_threads(),
                )?;

                let mut written: u64 = 0;
                for record_bytes in buffer.iter_sorted() {
                    writer.write_raw_record(record_bytes)?;
                    written += 1;
                }

                let index = writer.finish()?;

                let index_path = bai_sidecar_path(output);
                write_bai_index(&index_path, &index)?;
                info!("Wrote BAM index: {}", index_path.display());
                Ok(written)
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks
            let memory_chunks: Vec<InMemoryChunk<RawCoordinateKey>> = if buffer.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort_into_chunks(self.phase1_threads()))
                })
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                vec![buffer.drain_into_single_chunk()]
            };

            let memory_chunks = MemorySources::Shared(memory_chunks);
            let n_memory = memory_chunks.num_non_empty();
            log_run_formation(timer.spill_count(), stats.runs_written);
            debug!(
                "Phase 2: Merging {} chunks with index generation...",
                chunk_files.len() + n_memory
            );

            stats.output_records = timer.time_merge(|| {
                let (index, records_merged) = self.merge_chunks_with_index::<RawCoordinateKey>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )?;

                let index_path = bai_sidecar_path(output);
                write_bai_index(&index_path, &index)?;
                info!("Wrote BAM index: {}", index_path.display());
                Ok(records_merged)
            })?;
        }

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.phase1_threads(), self.phase2_threads(), self.max_temp_files);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort by queryname order using keyed temp files for O(1) merge comparisons.
    ///
    /// Dispatches to the appropriate key type based on the comparator:
    /// - `Lexicographic`: uses `RawQuerynameLexKey` (byte comparison)
    /// - `Natural`: uses `RawQuerynameKey` (natural numeric comparison)
    fn sort_queryname(
        &self,
        record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
        comparator: QuerynameComparator,
    ) -> Result<RawSortStats> {
        use crate::keys::{RawQuerynameKey, RawQuerynameLexKey};
        debug!("Using queryname sort with {comparator} comparator");
        match comparator {
            QuerynameComparator::Lexicographic => self.sort_queryname_keyed::<RawQuerynameLexKey>(
                record_source,
                pool,
                header,
                output,
                alloc,
            ),
            QuerynameComparator::Natural => self.sort_queryname_keyed::<RawQuerynameKey>(
                record_source,
                pool,
                header,
                output,
                alloc,
            ),
        }
    }

    /// Generic queryname sort using a specific key type.
    #[allow(clippy::too_many_lines)]
    fn sort_queryname_keyed<K: RawSortKey + Default + 'static>(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        use crate::keys::SortContext;

        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        let ctx = SortContext::from_header(header);

        // Estimate capacity from initial_capacity to avoid huge upfront allocations.
        let init_cap = self.effective_initial_capacity();
        let estimated_records = init_cap / 300;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut entries: Vec<(K, fgumi_raw_bam::RawRecord)> = Vec::with_capacity(estimated_records);
        let mut memory_used = 0usize;
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        // Natural run formation: consecutive chunks already in queryname order
        // extend one run instead of each becoming its own merge source.
        let mut run_former: RunFormer<K> = RunFormer::default();
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        debug!("Phase 1: Reading and sorting chunks (keyed output)...");
        let mut probe = SpillProbe::new("phase1");

        for record in record_source.by_ref() {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Extract key from raw bytes. Stamp the ingest position within this
            // chunk so the key is totally ordered: read name + flags alone is not
            // (secondary alignments of one read collide), and without a tiebreak
            // the chunk sort would have to be stable to preserve ingest order.
            let mut key = K::extract(record.as_ref(), &ctx);
            key.set_position(
                u32::try_from(entries.len()).expect("chunk length is capped at MAX_CHUNK_RECORDS"),
            );

            // Estimate memory: record bytes + key overhead
            let record_size = record.as_ref().len() + 50; // approximate key size
            memory_used += record_size;

            entries.push((key, record));

            if probe.should_sample_read(stats.total_records) {
                let bstats = BufferProbeStats::simple(memory_used as u64, entries.len() as u64);
                probe.log_mid_read(bstats, Some(pool.phase1_queue_depths()));
            }

            // Check if we need to spill to disk
            if should_spill_chunk(memory_used, self.memory_limit, entries.len()) {
                timer.end_read_span();
                let bstats = BufferProbeStats::simple(memory_used as u64, entries.len() as u64);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete
                self.drain_pending_spill::<K>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(bstats, Some(pool.phase1_queue_depths()));

                timer.time_sort(|| {
                    use rayon::prelude::*;
                    // `entries` is in ingest order, and exact-queryname-key ties must keep
                    // that order — matching `samtools sort -n`, fgbio, and the arena/runall
                    // path. The key carries its ingest position (see
                    // `RawSortKey::set_position`), so it is a total order and an unstable
                    // sort reproduces ingest order for ties without paying for a stable
                    // sort. The loser-tree merge breaks cross-chunk ties by chunk (ingest)
                    // order, so the whole queryname sort stays deterministic.
                    rayon_pool.install(|| entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0)));
                });

                let bounds = chunk_bounds(&entries, |(k, _)| k.clone());

                pending_spill = Some(Self::spill_chunk::<K>(
                    &mut run_former,
                    &chunk_files,
                    bounds,
                    &mut namer,
                    &pool,
                    &mut timer,
                    |writer| {
                        for (key, record) in entries.drain(..) {
                            writer.write_record(&key, record.as_ref())?;
                        }
                        Ok(())
                    },
                )?);

                memory_used = 0;
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        progress.log_final();
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<K>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(memory_used as u64);

        // Ingest is done: everything from here is Phase 2 (merge/write).
        self.enter_output_phase(&pool);

        if chunk_files.is_empty() {
            // All records fit in memory
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                use rayon::prelude::*;
                // Unstable sort over a totally ordered key: position is part of the key,
                // so ingest order is preserved for exact-key ties (see the per-chunk sort
                // above for the full rationale).
                rayon_pool.install(|| entries.par_sort_unstable_by(|a, b| a.0.cmp(&b.0)));
            });

            stats.output_records = timer.time_write_output(|| {
                use crate::pooled_bam_writer::PooledBamWriter;
                let output_header = self.create_output_header(header);
                let mut writer = PooledBamWriter::new(Arc::clone(&pool), output, &output_header)?;

                let mut written: u64 = 0;
                for (_key, record) in entries {
                    writer.write_raw_record(&record)?;
                    written += 1;
                }
                writer.finish()?;
                Ok(written)
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks (avoids
            // intermediate merge back into a single sorted buffer). The
            // queryname path accumulates records as individual `RawRecord`
            // allocations upstream (per-record alloc already paid), so we
            // sort the keyed `Vec`s in place and pack the owned `RawRecord`s
            // into `MemorySources::Owned` for the merge interface — unlike
            // the inline-buffer paths (coordinate, template), which produce
            // `InMemoryChunk`s sharing an `Arc<SegmentedBuf>`.
            let keyed_chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>> = if entries.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    use rayon::prelude::*;
                    let chunk_size = entries.len().div_ceil(self.phase1_threads());
                    rayon_pool.install(|| {
                        entries.par_chunks_mut(chunk_size).for_each(|chunk| {
                            // Unstable sort per chunk: the key's ingest position totally
                            // orders records within a chunk, so ties keep ingest order
                            // without stability; the loser-tree merge preserves it across
                            // chunks.
                            chunk.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                        });
                    });
                    // Carve sub-chunks aligned with par_chunks_mut boundaries.
                    // par_chunks_mut produces [0..cs), [cs..2cs), ..., [n-tail..n).
                    // We peel the short tail first, then full chunks from the end,
                    // and reverse. Each split_off is O(chunk_size) → O(n) total.
                    let mut remaining = std::mem::take(&mut entries);
                    let num_chunks = remaining.len().div_ceil(chunk_size);
                    let mut chunks: Vec<Vec<(K, fgumi_raw_bam::RawRecord)>> =
                        Vec::with_capacity(num_chunks);
                    let tail_len = remaining.len() % chunk_size;
                    if tail_len != 0 {
                        let split_at = remaining.len() - tail_len;
                        chunks.push(remaining.split_off(split_at));
                    }
                    while !remaining.is_empty() {
                        let split_at = remaining.len().saturating_sub(chunk_size);
                        chunks.push(remaining.split_off(split_at));
                    }
                    chunks.reverse();
                    chunks
                })
            } else {
                timer.time_sort(|| {
                    // Unstable sort: the key's ingest position preserves ingest order
                    // for exact-key ties.
                    entries.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                });
                vec![entries]
            };
            let memory_chunks = MemorySources::Owned(keyed_chunks);

            let n_memory = memory_chunks.num_non_empty();
            log_run_formation(timer.spill_count(), stats.runs_written);
            debug!(
                "Phase 2: Merging {} chunks (keyed comparisons)...",
                chunk_files.len() + n_memory
            );

            // Merge disk chunks + in-memory records using keyed comparisons
            stats.output_records = timer.time_merge(|| {
                self.merge_chunks_generic::<K>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )
            })?;
        }

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.phase1_threads(), self.phase2_threads(), self.max_temp_files);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Sort by template-coordinate order using inline buffer for reduced memory.
    ///
    /// Uses `TemplateRecordBuffer` which stores records in a single contiguous allocation
    /// with packed sort keys, eliminating per-record heap allocations for names.
    ///
    /// Writes keyed temp chunks that preserve pre-computed sort keys, enabling O(1)
    /// comparisons during merge (instead of expensive CIGAR/aux parsing).
    fn sort_template_coordinate(
        &self,
        mut record_source: RecordSource,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
    ) -> Result<RawSortStats> {
        let lib_lookup = LibraryLookup::from_header(header);
        let cb_hasher = cb_hasher();

        // Peek the first record to provision the variant (extraction stays
        // full-width). The first record is handed (by reference) to the impl,
        // which processes it before draining the rest of `record_source`; the
        // wrapper only extracts `first_key` to choose `K` — it does not re-chain
        // the record into the iterator.
        let first_record = record_source.by_ref().next();
        let first_key = first_record.as_ref().map(|r| {
            extract_template_key_inline(r.as_ref(), &lib_lookup, self.cell_tag, &cb_hasher)
        });
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        let variant =
            select_template_variant(first_key.as_ref(), self.key_types, header_library_varies);

        debug!(
            "template-coordinate variant: {} lanes (cb={}, tertiary={})",
            variant.lanes(),
            variant.cb,
            variant.tertiary
        );

        match (variant.cb, variant.tertiary) {
            (false, false) => self.sort_template_coordinate_impl::<TemplateKey24>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
            (true, false) => self.sort_template_coordinate_impl::<CbKey32>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
            (false, true) => self.sort_template_coordinate_impl::<TertKey32>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
            (true, true) => self.sort_template_coordinate_impl::<TemplateKey40>(
                record_source,
                first_record.as_ref(),
                first_key,
                variant,
                pool,
                header,
                output,
                alloc,
                &lib_lookup,
                &cb_hasher,
            ),
        }
    }

    /// Generic template-coordinate sort body, monomorphized over the chosen
    /// lane-key type `K`.
    ///
    /// The dispatch wrapper [`Self::sort_template_coordinate`] peeks the first
    /// record to choose `K`, then hands the captured `first_record` (already
    /// pulled out of `record_source`, borrowed here) plus the rest of the
    /// iterator. This impl processes `first_record` in a pre-loop block and
    /// drains the remainder in the read loop. `first_key` is the full key of the
    /// first record; the decode-time verify compares every subsequent record's
    /// dropped lanes against it.
    #[allow(clippy::too_many_lines)]
    // The generic body keeps the same parameter shape as the legacy function
    // plus the peeked-first-record handoff; the existing code already allows
    // this pattern for the merge helpers.
    #[allow(clippy::too_many_arguments)]
    fn sort_template_coordinate_impl<K: TemplateLaneKey + RawSortKey + Default + 'static>(
        &self,
        mut record_source: RecordSource,
        first_record: Option<&fgumi_raw_bam::RawRecord>,
        first_key: Option<TemplateKey>,
        variant: TemplateKeyVariant,
        pool: Arc<SortWorkerPool>,
        header: &Header,
        output: &Path,
        alloc: &mut TmpDirAllocator,
        lib_lookup: &LibraryLookup,
        cb_hasher: &ahash::RandomState,
    ) -> Result<RawSortStats> {
        let mut stats = RawSortStats::default();
        let mut timer = SortPhaseTimer::new();

        // The full key the decode-time verify compares dropped lanes against.
        // On empty input, `first_record` is `None` so the pre-loop push is
        // skipped and the read loop runs zero times — `first` is never read in
        // that case, so a default (all-zero) key is a harmless placeholder that
        // lets the function fall through to the post-loop empty-output path
        // (header-only BAM), matching the coordinate and queryname orders.
        let first = first_key.unwrap_or_default();

        // Per-record arena footprint (data side: inline header + BAM bytes; ref side:
        // the cached key + offset/len/pad). Ref size tracks the chosen key width K:
        //   K=TemplateKey24 → ref=40 B, bpr=298 B, data%≈86.6%
        //   K=TemplateKey32 → ref=48 B, bpr=306 B, data%≈84.3%
        //   K=TemplateKey40 → ref=56 B, bpr=314 B, data%≈82.2%
        let ref_bytes = std::mem::size_of::<TemplateRecordRef<K>>();
        let data_bytes_per_record = TEMPLATE_HEADER_SIZE + EST_BAM_BYTES_PER_TEMPLATE_RECORD;
        let bytes_per_record = data_bytes_per_record + ref_bytes;

        let init_cap = self.effective_initial_capacity();
        let estimated_records = (init_cap / bytes_per_record).max(1);
        let estimated_data_bytes = init_cap * data_bytes_per_record / bytes_per_record;

        let mut chunk_files: Vec<PathBuf> = Vec::new();
        let mut buffer =
            TemplateRecordBuffer::<K>::with_capacity(estimated_records, estimated_data_bytes);
        let mut namer = ChunkNamer::new(alloc);
        let mut pending_spill: Option<PendingSpill> = None;
        // Natural run formation: consecutive chunks already in template order
        // extend one run instead of each becoming its own merge source.
        let mut run_former: RunFormer<K> = RunFormer::default();
        let rayon_pool = self.build_sort_rayon_pool()?;

        let progress = ProgressTracker::new("Read records").with_interval(1_000_000);
        debug!("Phase 1: Reading and sorting chunks (inline buffer)...");
        let mut probe = SpillProbe::new("phase1");

        // Process the captured first record before draining the rest: extract its
        // full key, verify the dropped lanes match `first` (trivially true for the
        // first record itself), and push the narrowed key. A single record cannot
        // exceed the memory limit, so no spill check is needed here.
        if let Some(record) = first_record {
            stats.total_records += 1;
            progress.log_if_needed(1);

            let bam_bytes = record.as_ref();
            let full = extract_template_key_inline(bam_bytes, lib_lookup, self.cell_tag, cb_hasher);
            if let Some(violation) = verify_dropped_lanes(&first, &full, variant) {
                let name = String::from_utf8_lossy(
                    fgumi_raw_bam::RawRecordView::new(bam_bytes).read_name(),
                )
                .into_owned();
                return Err(dropped_lane_error(&name, violation));
            }
            buffer.push(bam_bytes, K::from_full(&full))?;
        }

        // Borrow each record's bytes in place (see the coordinate ingest loop);
        // the key is extracted and the bytes copied into the buffer before the
        // borrow ends, so no owned `RawRecord` is needed here.
        while let Some(bam_bytes) = record_source.next_record_borrowed()? {
            stats.total_records += 1;
            progress.log_if_needed(1);

            // Extract the full template key, verify the lanes the chosen variant
            // dropped are constant relative to the first record, then push the
            // narrowed key.
            let full = extract_template_key_inline(bam_bytes, lib_lookup, self.cell_tag, cb_hasher);
            if let Some(violation) = verify_dropped_lanes(&first, &full, variant) {
                let name = String::from_utf8_lossy(
                    fgumi_raw_bam::RawRecordView::new(bam_bytes).read_name(),
                )
                .into_owned();
                return Err(dropped_lane_error(&name, violation));
            }
            buffer.push(bam_bytes, K::from_full(&full))?;

            if probe.should_sample_read(stats.total_records) {
                probe.log_mid_read(probe_stats(&buffer), Some(pool.phase1_queue_depths()));
            }

            // Check memory usage
            if buffer.memory_usage() >= self.memory_limit {
                timer.end_read_span();
                let bstats = probe_stats(&buffer);
                let depths = Some(pool.phase1_queue_depths());
                probe.pre_spill(bstats, depths);

                // Wait for any previous spill to complete
                self.drain_pending_spill::<K>(
                    &mut pending_spill,
                    &mut chunk_files,
                    &mut stats,
                    &mut timer,
                    &mut namer,
                    &pool,
                )?;
                probe.post_drain(probe_stats(&buffer), Some(pool.phase1_queue_depths()));

                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });

                let bounds = chunk_bounds(buffer.refs(), |r| r.key);

                pending_spill = Some(Self::spill_chunk::<K>(
                    &mut run_former,
                    &chunk_files,
                    bounds,
                    &mut namer,
                    &pool,
                    &mut timer,
                    |writer| {
                        for (key, record) in buffer.iter_sorted_keyed() {
                            writer.write_record(&key, record)?;
                        }
                        Ok(())
                    },
                )?);

                buffer.clear();
                force_mi_collect();
                probe.post_spill(Some(pool.phase1_queue_depths()));
                timer.begin_read_span();
            }
        }

        timer.end_read_span();
        progress.log_final();
        if let Some(err) = record_source.take_error() {
            return Err(anyhow::Error::from(err));
        }

        // Drain any pending spill before merge
        self.drain_pending_spill::<K>(
            &mut pending_spill,
            &mut chunk_files,
            &mut stats,
            &mut timer,
            &mut namer,
            &pool,
        )?;
        probe.phase1_end(buffer.memory_usage() as u64);

        // Ingest is done: everything from here is Phase 2 (merge/write).
        self.enter_output_phase(&pool);

        if chunk_files.is_empty() {
            // All records fit in memory
            debug!("All records fit in memory, performing in-memory sort");

            timer.time_sort(|| {
                rayon_pool.install(|| buffer.par_sort());
            });

            stats.output_records = timer.time_write_output(|| {
                use crate::pooled_bam_writer::PooledBamWriter;
                let output_header = self.create_output_header(header);
                let mut writer = PooledBamWriter::new(Arc::clone(&pool), output, &output_header)?;

                let mut written: u64 = 0;
                for record_bytes in buffer.iter_sorted() {
                    writer.write_raw_record(record_bytes)?;
                    written += 1;
                }
                writer.finish()?;
                Ok(written)
            })?;
        } else {
            // Sort remaining records into separate sub-array chunks (avoids
            // intermediate merge back into a single sorted buffer)
            let memory_chunks: Vec<InMemoryChunk<K>> = if buffer.is_empty() {
                Vec::new()
            } else if self.phase1_threads() > 1 {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort_into_chunks(self.phase1_threads()))
                })
            } else {
                timer.time_sort(|| {
                    rayon_pool.install(|| buffer.par_sort());
                });
                vec![buffer.drain_into_single_chunk()]
            };

            let memory_chunks = MemorySources::Shared(memory_chunks);
            let n_memory = memory_chunks.num_non_empty();
            log_run_formation(timer.spill_count(), stats.runs_written);
            debug!("Phase 2: Merging {} chunks...", chunk_files.len() + n_memory);

            // Merge using O(1) key comparisons
            stats.output_records = timer.time_merge(|| {
                self.merge_chunks_generic::<K>(
                    &chunk_files,
                    memory_chunks,
                    header,
                    output,
                    stats.total_records,
                    &pool,
                )
            })?;
        }

        if let Ok(pool) = Arc::try_unwrap(pool) {
            pool.shutdown();
        }
        timer.log_summary(self.phase1_threads(), self.phase2_threads(), self.max_temp_files);
        debug!("Sort complete: {} records processed", stats.total_records);

        Ok(stats)
    }

    /// Build chunk sources from disk files and in-memory chunks.
    ///
    /// Disk chunks become `PoolDisk` sources: the shared worker pool reads and
    /// decompresses them in the background while the main thread parses records,
    /// so no per-source threads are spawned. Both the plain
    /// ([`merge_chunks_generic`]) and indexed ([`merge_chunks_with_index`])
    /// merges go through this pool-integrated path.
    fn build_chunk_sources<K: RawSortKey + Default + 'static>(
        chunk_files: &[PathBuf],
        memory_chunks: MemorySources<K>,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<Vec<ChunkSource<K>>> {
        let num_disk = chunk_files.len();
        let num_memory = memory_chunks.num_non_empty();
        let mut sources: Vec<ChunkSource<K>> = Vec::with_capacity(num_disk + num_memory);

        if !chunk_files.is_empty() {
            // Pool-integrated path: install per-file Phase 2 state on the
            // pool, then create one PoolDisk source per file. Workers
            // cooperatively read+decompress all files via work-stealing.
            pool.set_phase2_files(chunk_files)?;

            for source_id in 0..num_disk {
                sources.push(ChunkSource::PoolDisk { source_id, scratch: Vec::new() });
            }
        }

        memory_chunks.push_into(&mut sources);

        Ok(sources)
    }

    /// Build the pooled merge sources and activate Phase 2.
    ///
    /// Shared by both the plain ([`merge_chunks_generic`]) and indexed
    /// ([`merge_chunks_with_index`]) merges so the Phase-2 lifecycle is
    /// single-sourced and the two paths cannot drift. Returns the sources plus
    /// an RAII [`Phase2Guard`] (borrowing `pool`); callers finish through
    /// [`Phase2Guard::finish_output`], and `Drop` resets Phase 2 on any path
    /// that does not get that far.
    ///
    /// Phase 2 is armed whether or not there are disk chunks: it is what lets
    /// workers schedule `SortStep::Phase2FileWork` and reach the merge's
    /// per-file reorder buffers. With no spill files `try_phase2_file_work` sees
    /// an empty file vector and returns immediately, so the only work left is
    /// compression — which is correct in any phase, since each block carries its
    /// own [`CompressTarget`](crate::worker_pool::CompressTarget). The consumer
    /// is created only for disk sources — it drains the pool's per-file reorder
    /// buffers, and there is nothing to drain when everything is already in
    /// memory.
    fn setup_phase2_merge<'a, K: RawSortKey + Default + 'static>(
        chunk_files: &[PathBuf],
        memory_chunks: MemorySources<K>,
        pool: &'a Arc<SortWorkerPool>,
    ) -> Result<(Vec<ChunkSource<K>>, Phase2Guard<'a, K>)> {
        let num_disk = chunk_files.len();
        if num_disk > 0 {
            debug!(
                "Pool-integrated merge: {} disk sources, {} pool workers (N+2 model)",
                num_disk,
                pool.num_workers()
            );
        }

        let sources = Self::build_chunk_sources::<K>(chunk_files, memory_chunks, pool)?;
        debug!("Merging from {} sources...", sources.len());

        // Activate Phase 2 for the whole merge, spill files or not, so workers can
        // schedule `SortStep::Phase2FileWork` and reach the per-file reorder
        // buffers. Compression level is no longer a reason to arm it — each block
        // carries its own `CompressTarget` and is compressed correctly in any phase.
        // The consumer is created only for disk sources: it holds an Arc snapshot of
        // the pool's per-file Phase 2 state — workers populate per-file reorder
        // buffers, the consumer pops from them. There is nothing to consume when
        // everything is already in memory.
        let consumer = (num_disk > 0).then(|| {
            MainThreadChunkConsumer::new(
                pool.phase2_files(),
                pool.decompress_error_flag(),
                pool.chunk_read_error_flag(),
                pool.worker_panicked_flag(),
                pool.shared_state(),
            )
        });
        pool.set_phase(crate::worker_pool::phase::PHASE2);
        let guard = Phase2Guard { pool, consumer, active: true, sources_released: false };

        Ok((sources, guard))
    }

    /// Say, in one line, what limited the merge.
    ///
    /// Two numbers decide it. `utilization` is the worker pool's share of its
    /// own capacity; `fetch_fraction` is how much of the consumer's loop went
    /// into fetching the next record, which is where it blocks waiting for a
    /// decompressed block. The interesting case is *both* being unfavourable:
    /// if the pool is idle **and** the consumer is waiting, then neither side is
    /// the bottleneck and the merge is waiting on storage.
    ///
    /// Deliberately hedged ("suggests", "likely"). These thresholds come from a
    /// handful of measured runs, not a model, and the cost of a confident wrong
    /// diagnosis -- someone raising `--threads` on an I/O-bound sort -- is worse
    /// than the cost of a vague right one. The numbers above the verdict are the
    /// evidence; this line only points at them.
    fn log_merge_verdict(utilization: Option<f64>, fetch_fraction: f64) {
        use crate::merge_phases::{MergeVerdict, classify_merge};

        let Some(utilization) = utilization else { return };
        let (util_pct, fetch_pct) = (100.0 * utilization, 100.0 * fetch_fraction);

        match classify_merge(utilization, fetch_fraction) {
            MergeVerdict::CpuBound => info!(
                "  Verdict: worker pool saturated ({util_pct:.0}%); the merge is CPU-bound, so \
                 more threads may help"
            ),
            MergeVerdict::IoBound => info!(
                "  Verdict: workers idle ({util_pct:.0}% of capacity) while the consumer spent \
                 {fetch_pct:.0}% of its loop waiting for data -- neither side is the constraint. \
                 More threads are unlikely to help. Storage is one candidate; the pipeline's own \
                 coordination -- per-file read-ahead depth, lock contention, how late idle \
                 workers notice new work -- is another, and the Merge Stalls block below is what \
                 separates them."
            ),
            MergeVerdict::Mixed => info!(
                "  Verdict: worker pool {util_pct:.0}% utilized, consumer {fetch_pct:.0}% \
                 waiting on data -- neither clearly saturated"
            ),
        }
    }

    /// Log the merge's component breakdown: main-thread work and worker work.
    ///
    /// The two halves are measured differently and must be read differently.
    /// Consumer figures are sampled 1-in-N and scaled, so they estimate that one
    /// thread's wall time and do partition `loop wall`. Worker figures are exact
    /// busy-time sums across every pool thread, so they overlap each other AND
    /// the consumer, and routinely exceed `loop wall` -- see
    /// [`crate::merge_phases`]. Comparing a worker figure to `loop wall` is
    /// meaningless; comparing worker figures to each other is the point.
    ///
    /// The two wall clocks are therefore separate arguments, and each denominator
    /// takes the one it belongs to. `loop_total` covers the consumer loop alone,
    /// which is what the sampled consumer rows partition. `merge_total` runs
    /// through output finalization, which is where the queued tail of
    /// `output_compress` is actually done, so it is what worker utilization
    /// divides by; using `loop_total` there would charge worker busy time to a
    /// window that ended before some of it happened.
    #[allow(clippy::cast_precision_loss)]
    /// Log the consumer's CPU cost, partitioned from exact totals.
    ///
    /// The sampled sub-phase rows are magnitude-biased -- they routinely sum past
    /// loop wall, and `log_merge_sub_phases` says so -- but their *ratios* hold up
    /// far better than their absolute values, because the bias comes from scaling
    /// whichever samples happened to land on a stall and inflates the buckets
    /// roughly together. So this takes the split from sampling and the total from
    /// quantities that are exact: loop wall minus the two stalls that are timed
    /// exactly. The breakdown then partitions the loop by construction rather
    /// than overshooting it by 30-50%.
    fn log_consumer_cpu(
        sampled: (f64, f64, f64, f64),
        records: (u64, f64),
        stalls: Option<&crate::merge_stalls::ConsumerStallReport>,
    ) {
        let (loop_total, est_read, est_tree, est_write) = sampled;
        let (records_merged, blocked_secs) = records;
        if records_merged == 0 {
            return;
        }
        let park_secs = stalls.map_or(0.0, |s| s.park_secs);
        let consumer_cpu = (loop_total - park_secs - blocked_secs).max(0.0);
        #[allow(clippy::cast_precision_loss, reason = "record counts stay far below 2^52")]
        let records = records_merged as f64;
        // Park is inside the sampled "fetch" bucket, and output backpressure is
        // inside the sampled "write" bucket. `consumer_cpu` already excludes
        // both, so remove them from their buckets before taking ratios or the
        // waits would be counted on both sides and inflate their shares.
        let fetch_cpu = (est_read - park_secs).max(0.0);
        let write_cpu = (est_write - blocked_secs).max(0.0);
        let ratio_total = fetch_cpu + est_tree + write_cpu;
        info!(
            "  Consumer CPU: {consumer_cpu:.1}s exact (loop {loop_total:.1}s - park \
             {park_secs:.1}s - backpressure {blocked_secs:.1}s) = {:.0} ns/record",
            consumer_cpu * 1e9 / records
        );
        if ratio_total <= 0.0 {
            return;
        }
        let share = |part: f64| consumer_cpu * part / ratio_total;
        let ns = |part: f64| share(part) * 1e9 / records;
        info!("    parse + key extract: {:.1}s ({:.0} ns/rec)", share(fetch_cpu), ns(fetch_cpu));
        info!("    loser tree:          {:.1}s ({:.0} ns/rec)", share(est_tree), ns(est_tree));
        info!("    enqueue write:       {:.1}s ({:.0} ns/rec)", share(write_cpu), ns(write_cpu));
        info!(
            "    (totals exact; the split between them is sampled, so read the shares as \
             proportions rather than to the tenth of a second)"
        );
    }

    fn log_merge_sub_phases(
        walls: (f64, f64),
        consumer: (f64, f64, f64),
        sampling: (u64, u64),
        active_workers: usize,
        pool: &Arc<SortWorkerPool>,
        stalls: Option<crate::merge_stalls::ConsumerStallReport>,
        consumer_diag: MergeConsumerDiag,
    ) {
        let (loop_total, merge_total) = walls;
        let (write_secs, read_secs, tree_secs) = consumer;
        let (samples_taken, records_merged) = sampling;
        if records_merged == 0 {
            return;
        }
        #[allow(clippy::cast_precision_loss, reason = "sample counts stay far below 2^52")]
        let scale =
            if samples_taken > 0 { records_merged as f64 / samples_taken as f64 } else { 1.0 };
        let (est_write, est_read, est_tree) =
            (write_secs * scale, read_secs * scale, tree_secs * scale);

        info!("=== Merge Sub-Phase Timing ===");
        info!(
            "  Consumer (main thread; {samples_taken} samples of {records_merged} records, scaled {scale:.0}x)"
        );
        info!("    Fetch next record: {est_read:.1}s  (includes waiting on decompressed blocks)");
        let MergeConsumerDiag {
            backpressure_secs: blocked_secs,
            backpressure_waits: blocked_waits,
            borrowed,
            reassembled,
        } = consumer_diag;
        if blocked_waits > 0 {
            #[allow(clippy::cast_precision_loss, reason = "wait counts stay far below 2^52")]
            let mean_us = blocked_secs * 1e6 / blocked_waits as f64;
            info!(
                "    Output backpressure: {blocked_secs:.1}s over {blocked_waits} waits \
                 (mean {mean_us:.0} us, exact) -- the consumer blocked for an output permit"
            );
        } else {
            info!("    Output backpressure: none -- the compressors always had a permit ready");
        }
        if borrowed + reassembled > 0 {
            #[allow(clippy::cast_precision_loss, reason = "record counts stay far below 2^52")]
            let pct = 100.0 * reassembled as f64 / (borrowed + reassembled) as f64;
            info!(
                "    Record presentation: {borrowed} borrowed zero-copy, {reassembled} \
                 reassembled across a block boundary ({pct:.2}%, exact)"
            );
        }
        info!("    Loser tree:        {est_tree:.1}s");

        Self::log_consumer_cpu(
            (loop_total, est_read, est_tree, est_write),
            (records_merged, blocked_secs),
            stalls.as_ref(),
        );

        // Parking is a subset of "fetch next record", so exact park time cannot
        // legitimately exceed the sampled estimate of it. When it does, the
        // sample missed stalls: parks are rare per record and expensive when
        // they happen, and a 1-in-1021 record sample of a heavy-tailed event has
        // enough variance to land far off in either direction. The exact figure
        // in the stall block below supersedes this row whenever they disagree.
        if let Some(s) = stalls
            && s.park_secs > est_read * 1.05
        {
            info!(
                "    NOTE: exact park time is {:.1}s, above the sampled fetch estimate of \
                 {est_read:.1}s -- the sample under-caught stalls; prefer the Merge Stalls block",
                s.park_secs
            );
        }

        let workers = pool.merge_phase_breakdown();
        if !workers.is_empty() {
            let busy = workers.total_busy_secs();
            let pct = |secs: f64| if busy > 0.0 { 100.0 * secs / busy } else { 0.0 };
            let (read_s, read_n) = workers.read;
            let (dec_s, dec_n) = workers.decompress;
            let (comp_s, comp_n) = workers.output_compress;
            let (spill_s, spill_n) = workers.spill_compress;
            let workers_n = active_workers;
            info!(
                "  Workers ({workers_n} active threads; busy time, overlaps the above and itself)"
            );
            info!("    Spill disk read:   {read_s:.1}s ({:.0}%) [{read_n} batches]", pct(read_s));
            // Blocks-per-batch, split by allowance. A deep share near zero means
            // the deep-read gate (drain frontier or awaited source) is not
            // firing, which reads identically to "the deeper read-ahead did not
            // help" unless it is reported separately.
            let (deep_b, deep_blk, shal_b, shal_blk) = pool.read_batch_split();
            #[allow(clippy::cast_precision_loss, reason = "block counts stay far below 2^52")]
            let per = |blk: u64, b: u64| if b == 0 { 0.0 } else { blk as f64 / b as f64 };
            info!(
                "      deep {deep_b} batches / {deep_blk} blocks ({:.1} per batch), \
                 other {shal_b} batches / {shal_blk} blocks ({:.1} per batch)",
                per(deep_blk, deep_b),
                per(shal_blk, shal_b)
            );
            info!("    Spill decompress:  {dec_s:.1}s ({:.0}%) [{dec_n} blocks]", pct(dec_s));
            info!("    Output compress:   {comp_s:.1}s ({:.0}%) [{comp_n} blocks]", pct(comp_s));
            info!("    Total worker busy: {busy:.1}s  (NOT comparable to loop wall clock)");
            // Utilization is the thread-efficiency question: well below 100%
            // means the pool idled, the merge was bound by something other than
            // worker CPU, and adding compression threads cannot help.
            if let Some(util) = workers.worker_utilization(merge_total, workers_n) {
                info!(
                    "    Worker utilization: {:.0}% of {workers_n} active threads x \
                     {merge_total:.1}s",
                    100.0 * util
                );
            }
            // Phase 1 spill compression rides the same worker step, so it is
            // reported for context but excluded from the merge totals above.
            info!("  Phase 1 (not part of the merge)");
            info!("    Spill compress:    {spill_s:.1}s [{spill_n} blocks]");

            // Utilization over the full merge window; the fetch fraction stays
            // relative to the consumer loop, which is the only thing it is a
            // fraction of.
            Self::log_merge_verdict(
                workers.worker_utilization(merge_total, workers_n),
                if loop_total > 0.0 { est_read / loop_total } else { 0.0 },
            );
        }
        Self::log_merge_stalls(loop_total, merge_total, active_workers, stalls, pool);
        info!("==============================");
    }

    /// Log why the merge stalled, as opposed to where its time went.
    ///
    /// The three blocks answer one question each, and only together: the
    /// consumer's exact park time and the state of the other files when it
    /// parked; why worker file-scans came back empty; and how much of the
    /// workers' idle time was spent asleep *through* work becoming available.
    /// See [`crate::merge_stalls`].
    ///
    /// `merge_total` and `active_workers` are the utilization denominators and
    /// must be the same pair [`Self::log_merge_sub_phases`] divides by, because
    /// `SATURATED_UTILIZATION` is shared between `classify_merge` and
    /// `classify_stall` so the two verdicts agree about whether the pool was
    /// busy. `loop_total` stays the denominator for the consumer's park
    /// fraction, which is a fraction of the merge loop alone.
    #[allow(clippy::cast_precision_loss)]
    fn log_merge_stalls(
        loop_total: f64,
        merge_total: f64,
        active_workers: usize,
        stalls: Option<crate::merge_stalls::ConsumerStallReport>,
        pool: &Arc<SortWorkerPool>,
    ) {
        use crate::merge_stalls::{Phase2Skip, ScanVerdict};

        // Utilization gates the stall shape: "no worker on it" means a
        // scheduling defect on an idle pool and plain saturation on a busy one.
        // The active cap, not the pool width -- a Phase 2 narrower than Phase 1
        // otherwise charges the merge for workers that were never allowed to
        // help, and reports a saturated pool as idle.
        let utilization = pool
            .merge_phase_breakdown()
            .worker_utilization(merge_total, active_workers)
            .unwrap_or(0.0);

        let scans = pool.phase2_scan_report();
        let wake = pool.wake_latency_report();
        let stalls = stalls.filter(|s| !s.is_empty());
        if !merge_stalls_are_silent(stalls.as_ref(), &scans, &wake) {
            info!("=== Merge Stalls ===");

            if let Some(s) = stalls {
                Self::log_consumer_stalls(loop_total, utilization, s);
            }

            if !scans.is_empty() {
                let reasons = Phase2Skip::ALL
                    .iter()
                    .map(|&r| format!("{}={}", r.label(), scans.skips[r as usize]))
                    .collect::<Vec<_>>()
                    .join(" ");
                info!("  Worker scans finding no work: {} ({reasons})", scans.scans);
                info!(
                    "    Of those: {:.0}% backpressured, {:.0}% waiting on a peer's read, {:.0}% \
                     contended",
                    100.0 * scans.verdict_share(ScanVerdict::Backpressured),
                    100.0 * scans.verdict_share(ScanVerdict::IoWait),
                    100.0 * scans.verdict_share(ScanVerdict::Contended)
                );
            }

            if !wake.is_empty() {
                info!(
                    "  Worker discovery lag: ~{:.1}s of {:.1}s deep-sleep worker-seconds; {} \
                     sleeps ended in a find, {:.0}% of them after >=320us",
                    wake.estimated_discovery_lag_secs(),
                    wake.deep_sleep_idle_secs(),
                    wake.productive_sleep_count(),
                    100.0 * wake.deep_sleep_wake_share()
                );
                info!(
                    "    (the consumer unparks one worker when a reorder buffer drains, so this \
                     bounds how late work arriving any other way is noticed; it delays the merge \
                     only when every worker is asleep at once)"
                );
                // Why the discovery lag above is as large as it is. A wake aimed at a
                // worker that is already running does nothing, and the consumer then
                // waits for some other worker's backoff to expire -- bounded by
                // MAX_BACKOFF_US, not by unpark cost. `recoverable` is the subset where a
                // parked worker existed and the rotating target passed it over; the rest
                // is a genuinely saturated pool and is nobody's fault.
                let wakes = pool.wake_accounting();
                if wakes.issued > 0 {
                    info!(
                        "  Wake targeting: {} of {} wakes hit an already-running worker ({:.0}%), {} of \
                 those had a parked worker available ({:.0}% recoverable)",
                        wakes.on_running,
                        wakes.issued,
                        Self::percent(wakes.on_running, wakes.issued),
                        wakes.recoverable,
                        Self::percent(wakes.recoverable, wakes.issued)
                    );
                }
                info!(
                    "    Wakes issued: {} (backoff 10us doubling to {}us ceiling; one worker woken \
                     per wake)",
                    pool.wakes_issued(),
                    crate::worker_pool::MAX_BACKOFF_US
                );
            }

            // What the pool looked like when the consumer parked, which is the
            // question "why did nobody start on this block yet" reduced to three
            // mutually exclusive answers with three different fixes. Splitting park
            // *time* as well as park counts matters: a rare-but-long class is
            // invisible in counts alone. Outside the wake gate above so it stands
            // on its own park count: a heavy-parking consumer whose wake report is
            // empty would otherwise collect this census and never print it.
            let supply = pool.park_supply_report();
            if supply.total_parks() > 0 {
                const LABELS: [&str; crate::merge_stalls::ParkSupply::COUNT] =
                    ["a worker was asleep", "all busy, compress queued", "all busy merging"];
                info!("  Why nobody had started on the awaited block, at the moment of the park:");
                for (i, label) in LABELS.iter().enumerate() {
                    info!(
                        "    {label:<26} {:>10} parks ({:>4.1}%)  {:>7.1}s ({:>4.1}% of park time)",
                        supply.counts[i],
                        Self::percent(supply.counts[i], supply.total_parks()),
                        Self::secs(supply.nanos[i]),
                        Self::percent(supply.nanos[i], supply.total_nanos())
                    );
                }
            }

            Self::log_park_attribution(pool, loop_total);

            // Close this block before delegating: `log_block_lifecycle` opens and
            // closes its own, so without a terminator here the lifecycle block
            // reads as nested inside the stall block rather than following it.
            info!("====================");
        }

        // Outside the gate above, not inside it. The lifecycle reports record the
        // pipeline working and these three record it stalling, so a merge that
        // never stalled is precisely the run with a complete lifecycle trace and
        // nothing to say above. `log_block_lifecycle` carries its own gate.
        Self::log_block_lifecycle(pool);
    }

    /// Nanoseconds as seconds.
    #[expect(clippy::cast_precision_loss, reason = "nanosecond totals stay below 2^52")]
    fn secs(nanos: u64) -> f64 {
        nanos as f64 / 1e9
    }

    /// `part` as a percentage of `whole`, or 0.0 when `whole` is zero.
    #[expect(clippy::cast_precision_loss, reason = "counts stay far below 2^52")]
    fn percent(part: u64, whole: u64) -> f64 {
        if whole == 0 {
            return 0.0;
        }
        100.0 * part as f64 / whole as f64
    }

    /// Log a latency distribution for every stage a block passes through.
    ///
    /// [`crate::merge_phases::MergePhaseBreakdown`] already gives each stage a
    /// busy total and a block count, which yields a mean. A mean cannot settle
    /// an argument: output compression averages ~187 us over five million blocks
    /// and a uniform 187 us behaves nothing like a bimodal mix with a long tail,
    /// yet only the tail explains a worker being unavailable when the consumer
    /// needs one. So report count, total, mean and three percentiles for each,
    /// including the writer -- previously visible only as the consumer's
    /// backpressure wait, which jumped 0.0s to 29.9s across one sweep with no
    /// way to attribute it.
    ///
    /// `wasted visits` is the other half: files a worker passed over before one
    /// gave it work. The scan tally is published only when a scan finds nothing,
    /// so the walk to a *successful* claim was never counted -- and on an 89-way
    /// merge that walk is most of what an "idle" worker is doing.
    ///
    /// The distributions log at **debug** and the wasted-visit line at **info**,
    /// the same split [`Self::log_block_lifecycle`] uses: a percentile table is
    /// for an investigation, not for someone who just sorted a BAM. Collection
    /// is unconditional either way.
    fn log_stage_latency(
        pool: &SortWorkerPool,
        writer: &(
            crate::merge_trace::HistogramReport,
            crate::merge_trace::HistogramReport,
            crate::merge_trace::HistogramReport,
        ),
    ) {
        let stage = pool.stage_latency();
        let &(write_dur, reorder_wait, reorder_depth) = writer;
        let rows: [(&str, crate::merge_trace::HistogramReport); 6] = [
            ("read (batched)", stage.read.snapshot()),
            ("decompress spill", stage.decompress.snapshot()),
            ("compress output", stage.output_compress.snapshot()),
            ("compress spill (ph1)", stage.spill_compress.snapshot()),
            ("write block", write_dur),
            ("write reorder wait", reorder_wait),
        ];
        if rows.iter().all(|(_, r)| r.is_empty()) {
            return;
        }
        // Distributions go to debug and decision-changing numbers stay at info,
        // matching `log_block_lifecycle` above: a six-row percentile table is
        // what an investigation needs and far more than someone who just sorted
        // a BAM should have to scroll past. Collection stays unconditional --
        // gating collection is what made these questions unanswerable from logs
        // already in hand.
        debug!("=== Stage Latency ===");
        debug!(
            "  {:<22} {:>10} {:>9} {:>9} {:>9} {:>9} {:>9}",
            "stage", "count", "total", "mean", "p50", "p90", "p99"
        );
        for (label, r) in rows {
            if r.is_empty() {
                continue;
            }
            debug!(
                "  {label:<22} {:>10} {:>8.1}s {:>8.0}us {:>8}us {:>8}us {:>8}us",
                r.count,
                r.total_secs(),
                r.mean_micros(),
                r.percentile_micros(0.50),
                r.percentile_micros(0.90),
                r.percentile_micros(0.99)
            );
        }
        if !reorder_depth.is_empty() {
            debug!(
                "  Writer reorder depth: mean {:.1} blocks, p90 {}, p99 {} (blocks held waiting \
                 for an earlier serial)",
                reorder_depth.mean_micros(),
                reorder_depth.percentile_micros(0.90),
                reorder_depth.percentile_micros(0.99)
            );
        }
        // Closes the debug block above, so the separator does not outlive its
        // header when only info is enabled.
        debug!("=====================");
        // Stays at info: this one is a decision-changer, not a distribution. It
        // is how much of an "idle" worker is scanning rather than waiting, and
        // reading it wrong sends an investigation at the scan loop -- filling
        // the walk with more work per claim was measured at +36% (see
        // `PHASE2_DECOMP_CAP`).
        let claims = stage.useful_claims.load(std::sync::atomic::Ordering::Relaxed);
        if claims > 0 {
            info!(
                "  Wasted file visits: {:.1} per claim over {} claims ({} visits produced \
                 nothing) -- what an idle worker is actually doing",
                stage.wasted_visits_per_claim(),
                claims,
                stage.wasted_visits.load(std::sync::atomic::Ordering::Relaxed)
            );
        }
    }

    /// Log where consumer park time went, and what each worker did.
    ///
    /// The park table is the one figure in this report that partitions an exact
    /// total by construction, which is why it is the one to read first. Earlier
    /// attempts to explain park were shares of park *events* (blind to a
    /// rare-but-long cause) or sums over workers (which overcount, because the
    /// consumer waits for whichever worker arrives first). `unattributed` is the
    /// honesty term: a large value means the model is incomplete, not that the
    /// merge is idle for no reason.
    #[expect(
        clippy::cast_precision_loss,
        reason = "nanosecond and count totals are within f64's exact-integer range here"
    )]
    fn log_park_attribution(pool: &SortWorkerPool, loop_total: f64) {
        let park = pool.park_attribution_report();
        if park.is_empty() {
            return;
        }
        let total = park.total_nanos();
        if total == 0 {
            return;
        }
        let secs = |ns: u64| ns as f64 / 1e9;
        let share = |ns: u64| 100.0 * ns as f64 / total as f64;
        let per_park = |ns: u64| ns as f64 / park.parks as f64 / 1e3;

        info!("  Consumer park, by stage (exact, partitions the park):");
        info!("    {:<28} {:>8} {:>7} {:>10}", "stage", "time", "share", "per park");
        for (label, ns) in [
            ("waiting for a worker", park.to_claim_nanos),
            ("read + decompress work", park.work_nanos),
            ("waiting for its own wake", park.to_resume_nanos),
            ("unattributed", park.unattributed_nanos),
        ] {
            info!("    {label:<28} {:>7.1}s {:>6.0}% {:>9.0}us", secs(ns), share(ns), per_park(ns));
        }
        info!(
            "    {:<28} {:>7.1}s {:>6.0}% {:>9.0}us  over {} parks ({:.0}% of loop wall)",
            "TOTAL",
            secs(total),
            100.0,
            per_park(total),
            park.parks,
            if loop_total > 0.0 { 100.0 * secs(total) / loop_total } else { 0.0 }
        );
        info!(
            "    Blocks ready on the awaited file at resume: mean {:.2} (1.0 = every block \
             fetched on demand, one round trip per block)",
            park.mean_ready_on_resume()
        );
        info!(
            "    Parks with no claim during them: {} of {} ({:.0}%) -- the block was already \
             in flight or already done",
            park.unclaimed_parks,
            park.parks,
            100.0 * park.unclaimed_parks as f64 / park.parks as f64
        );

        let threads = pool.per_thread_report();
        let claims_total: u64 = threads.iter().map(|&(_, _, c)| c).sum();
        if claims_total == 0 {
            return;
        }
        info!("  Per worker (merge + phase 1 combined):");
        info!(
            "    {:>3} {:>9} {:>9} {:>7} {:>10} {:>5}",
            "wid", "busy", "idle", "busy%", "claims", "share"
        );
        for (w, &(busy, idle, claims)) in threads.iter().enumerate() {
            let denom = busy + idle;
            info!(
                "    {w:>3} {:>8.1}s {:>8.1}s {:>6.0}% {:>10} {:>4.0}%",
                secs(busy),
                secs(idle),
                if denom > 0 { 100.0 * busy as f64 / denom as f64 } else { 0.0 },
                claims,
                100.0 * claims as f64 / claims_total as f64
            );
        }
    }

    /// Log every stage of a spill block's journey, and the refill cycle.
    ///
    /// When the stall block above printed, it said the consumer waits for a
    /// block that is being produced; this says how long each step of producing
    /// it takes, and -- through the refill numbers -- how much of the wait is the
    /// pipeline working versus the pipeline not having started. When it stayed
    /// silent, these figures are the whole report, which is why they are not
    /// gated on it. See [`crate::merge_trace`].
    #[allow(clippy::cast_precision_loss)]
    /// Log how deep the pool was on the file the consumer was blocked on.
    ///
    /// Separates the two explanations for a starved consumer, which the park
    /// totals alone cannot: a pool that is barely on the awaited file is
    /// scheduling badly and can be steered, while one running near its tracked
    /// depth is paying the head block's decompress latency and cannot be helped
    /// by putting more workers on that same file.
    fn log_awaited_file_depth(
        consumer: &crate::merge_trace::ConsumerTraceReport,
        pool: &Arc<SortWorkerPool>,
    ) {
        use crate::merge_stalls::AwaitedState;
        use crate::merge_trace::MAX_TRACKED_IN_FLIGHT;

        // A merge that never stalled still has source runs to report, so this
        // block stands on its own park count rather than on the report being
        // non-empty -- otherwise it prints a header with nothing under it and
        // divides by zero to fill the line below.
        let parks = consumer.parks();
        if parks == 0 {
            return;
        }
        debug!("  Park duration by what the awaited file was doing");
        for state in AwaitedState::ALL {
            let hist = consumer.park_by_state[state as usize];
            if !hist.is_empty() {
                debug!("    {:<14} {}", state.label(), hist.summary());
            }
        }
        // The same five states weighted by TIME, not by park count, and at info
        // because the two disagree and only one of them is a cost.
        //
        // The park-supply census made this concrete: by count its largest class
        // was "all busy, compress queued" at 41%, by time that class was 13%
        // while "a worker was asleep" was 81%. Counts named the wrong fix,
        // confidently. These histograms have always carried the nanoseconds --
        // they were simply reported one level down, so every arm collected the
        // answer and none printed it.
        let state_secs: [f64; AwaitedState::COUNT] =
            std::array::from_fn(|i| consumer.park_by_state[i].total_secs());
        let state_total: f64 = state_secs.iter().sum();
        if state_total > 0.0 {
            let share = |i: usize| 100.0 * state_secs[i] / state_total;
            info!(
                "    The awaited file by park TIME ({state_total:.1}s): {:.0}% gap-filling, \
                 {:.0}% gap-stalled, {:.0}% decompressing, {:.0}% raw-queued, {:.0}% starved",
                share(AwaitedState::ReorderGapFilling as usize),
                share(AwaitedState::ReorderGapStalled as usize),
                share(AwaitedState::Decompressing as usize),
                share(AwaitedState::RawQueued as usize),
                share(AwaitedState::Starved as usize)
            );
        }

        let parks_pct = |count: u64| {
            #[allow(clippy::cast_precision_loss, reason = "park counts stay far below 2^52")]
            let pct = 100.0 * count as f64 / parks as f64;
            pct
        };
        let mean_depth = consumer.mean_in_flight();
        info!(
            "    Workers on the awaited file at a park: none {:.0}%, exactly one {:.0}%, \
             two or more {:.0}% (mean {mean_depth:.1}, tracked to {MAX_TRACKED_IN_FLIGHT})",
            parks_pct(consumer.idle_file_parks()),
            parks_pct(consumer.single_worker_parks()),
            parks_pct(consumer.multi_worker_parks()),
        );

        // Keyed on depth against the cap rather than on which bucket is
        // largest. "Two or more" becomes the majority long before the pool is
        // actually deep, so a bucket comparison would call a pool running two
        // deep saturated and retire a question the data has not answered.
        #[allow(clippy::cast_precision_loss, reason = "the cap is a small constant")]
        let cap = MAX_TRACKED_IN_FLIGHT as f64;
        // Which admission gate holds the depth where it is. Rendered next to
        // the depth itself: "only 1.9 deep of 8" is the symptom, and without
        // the gate that produced it the next step is a guess.
        // What the byte target actually derived. Without this a sizing
        // regression is invisible: the wall clock moves and nothing says which
        // batch produced it.
        let (mean_block_bytes, derived_cap, derived_batch) = pool.awaited_sizing();
        if mean_block_bytes > 0 {
            info!(
                "    Hot-file refill sized from measured blocks: {mean_block_bytes} B/block \
                 -> batch {derived_batch}, cap {derived_cap}"
            );
        }
        let skips = pool.awaited_skip_counts();
        let skip_total: u64 = skips.iter().sum();
        if skip_total > 0 {
            #[allow(clippy::cast_precision_loss, reason = "skip counts stay far below 2^52")]
            let pct = |n: u64| 100.0 * n as f64 / skip_total as f64;
            info!(
                "    Why the pool passed over the awaited file: raw-lock {:.0}%, raw-empty \
                 {:.0}%, decomp-lock {:.0}%, decomp-capped {:.0}% (of {skip_total})",
                pct(skips[0]),
                pct(skips[1]),
                pct(skips[2]),
                pct(skips[3]),
            );
        }
        if mean_depth >= cap / 2.0 {
            info!(
                "      -> the pool is running near its tracked depth on the file the merge is \
                 blocked on, so these parks are the head block's decompress latency; more \
                 concurrency on that file cannot shorten them"
            );
        } else if consumer.multi_worker_parks() > 0 {
            info!(
                "      -> the pool is on the awaited file but only {mean_depth:.1} deep of \
                 {MAX_TRACKED_IN_FLIGHT} tracked, so capacity is going unused on the file the \
                 merge is blocked on -- supply to that file, not decompress latency, is the \
                 first thing to check"
            );
        }
    }

    fn log_block_lifecycle(pool: &Arc<SortWorkerPool>) {
        use crate::merge_trace::EmptyCause;

        let life = pool.block_lifecycle_report();
        let refill = pool.refill_report();
        let consumer = pool.consumer_trace_report();
        let scans = pool.fruitless_scan_report();
        if block_lifecycle_is_silent(&life, &refill, &consumer, &scans) {
            return;
        }

        // Two audiences, one measurement. The distributions are what an
        // investigation needs and are far too much for someone who just sorted a
        // BAM, so the per-stage histograms go to debug while the numbers that
        // change a decision stay at info. Collection is unconditional either
        // way -- gating collection is what made the original question
        // unanswerable from logs we had already collected.
        info!("=== Merge Block Lifecycle ===");
        debug!("  disk read   -> {}", life.read_batch.summary());
        debug!("  raw dwell   -> {}   (queued, waiting for a worker)", life.raw_dwell.summary());
        debug!("  decompress  -> {}", life.decompress.summary());
        debug!(
            "  reorder     -> {}   (decompressed, waiting for the consumer)",
            life.reorder_dwell.summary()
        );
        info!(
            "  Per block: {:.0}us in the raw FIFO unclaimed, {:.0}us decompressing, {:.0}us \
             buffered before use (p50)",
            life.raw_dwell.percentile_micros(0.50),
            life.decompress.percentile_micros(0.50),
            life.reorder_dwell.percentile_micros(0.50)
        );
        if life.reorder_is_pass_through() {
            info!(
                "    NOTE: blocks are consumed almost as fast as they are inserted, so the \
                 reorder buffer is a pass-through and PHASE2_DECOMP_CAP is not the binding \
                 constraint -- however full the other files look"
            );
        }

        if !refill.is_empty() {
            info!("  Refill cycle ({} times a file's buffer ran dry)", refill.empties());
            info!(
                "    At the moment it emptied: {:.0}% had raw blocks unclaimed, {:.0}% already \
                 decompressing, {:.0}% nothing at all",
                100.0 * refill.cause_share(EmptyCause::RawReady),
                100.0 * refill.cause_share(EmptyCause::Decompressing),
                100.0 * refill.cause_share(EmptyCause::Dry)
            );
            debug!("    empty -> claimed  {}", refill.claim_lag.summary());
            debug!("    empty -> inserted {}", refill.insert_lag.summary());
            if !refill.read_lag.is_empty() {
                debug!("    empty -> read     {}", refill.read_lag.summary());
            }
            info!(
                "    -> {:.0}% of refill latency is spent waiting for a worker to START, {:.0}% \
                 doing the work",
                100.0 * refill.claim_share(),
                100.0 * (1.0 - refill.claim_share())
            );
        }

        if !consumer.is_empty() {
            // A merge that never stalled still has source runs to report, so the
            // park block has to stand on its own count rather than on the report
            // being non-empty -- otherwise it prints a header with nothing under
            // it and divides by zero to fill the line below.
            Self::log_awaited_file_depth(&consumer, pool);
            if !consumer.source_run_length.is_empty() {
                info!(
                    "  Consecutive blocks per source: {}",
                    consumer.source_run_length.summary_blocks()
                );
                if consumer.source_run_length.percentile_micros(0.90) <= 1 {
                    info!(
                        "    -> the merge switches source almost every block, so there is no hot \
                         file to prioritise; demand is spread across all runs at once"
                    );
                }
            }
        }

        if !scans.is_empty() {
            debug!("  Fruitless worker scan cost: {}", scans.summary());
        }
        info!("=============================");
    }

    /// Log where the merge loop blocked and what the other files were doing.
    ///
    /// Exact, not sampled -- see [`crate::merge_stalls`]. When this disagrees
    /// with the sampled "fetch next record" row above, this is the one to
    /// believe.
    #[allow(clippy::cast_precision_loss)]
    fn log_consumer_stalls(
        loop_total: f64,
        utilization: f64,
        s: crate::merge_stalls::ConsumerStallReport,
    ) {
        use crate::merge_stalls::{StallShape, classify_stall};

        let park_fraction = if loop_total > 0.0 { s.park_secs / loop_total } else { 0.0 };
        info!(
            "  Consumer parked: {:.1}s ({:.0}% of loop wall, exact) over {} parks",
            s.park_secs,
            100.0 * park_fraction,
            s.parks
        );
        info!(
            "    Block pulls that had to wait: {}/{} ({:.0}%), {:.1} parks each (1.0 = no wasted \
             wake-ups)",
            s.stalled_pulls,
            s.block_pulls,
            100.0 * s.stall_rate(),
            s.parks_per_stall()
        );
        info!(
            "    Worst source: #{} at {:.1}s ({:.0}% of park time; {} sources parked on)",
            s.top_source,
            s.top_source_park_secs,
            100.0 * s.top_source_share(),
            s.sources_parked_on
        );
        if s.censuses > 0 {
            info!(
                "    Other files at a park ({} parks sampled): {:.0}% at cap, {:.0}% starved, \
                 {:.0}% unreadable",
                s.censuses,
                100.0 * s.capped_share,
                100.0 * s.starved_share,
                100.0 * s.contended_share
            );
            info!(
                "    The awaited file: {:.0}% gap-filling, {:.0}% gap-stalled, {:.0}% \
                 decompressing, {:.0}% raw-queued, {:.0}% starved",
                100.0 * s.awaited.reorder_gap_filling,
                100.0 * s.awaited.reorder_gap_stalled,
                100.0 * s.awaited.decompressing,
                100.0 * s.awaited.raw_queued,
                100.0 * s.awaited.starved
            );
            info!(
                "      -> block not read yet {:.0}%, exists but unclaimed {:.0}%, being produced \
                 {:.0}%",
                100.0 * s.awaited.starved,
                100.0 * s.awaited.unclaimed(),
                100.0 * s.awaited.in_progress()
            );
        }
        match classify_stall(park_fraction, utilization, s.contended_share, s.awaited) {
            StallShape::NotStalled => info!("    Shape: the consumer is not waiting on blocks"),
            StallShape::PoolSaturated => info!(
                "    Shape: pool saturated -- the consumer waits because every worker is busy, \
                 which is what a healthy CPU-bound merge looks like. Fewer bytes to compress or \
                 more threads would help; nothing here is misscheduled"
            ),
            StallShape::HeadOfLine => info!(
                "    Shape: head-of-line -- the awaited file has nothing anywhere in its \
                 pipeline, so the block has not been read from disk yet. The constraint is \
                 upstream of the pool: storage, or read concurrency"
            ),
            StallShape::WorkUnclaimed => info!(
                "    Shape: work unclaimed -- the block the consumer needs already exists and no \
                 worker is on it. Capacity is not the problem; scheduling and wake latency are. \
                 Compare the discovery-lag line below"
            ),
            StallShape::DecompressLatency => info!(
                "    Shape: decompression latency -- a worker is already producing the needed \
                 block, so the consumer is paying the per-block cost serially. Check the reorder \
                 dwell below before reaching for a deeper cap: if blocks are consumed as fast as \
                 they are inserted, the buffer is not what the pipeline is running into"
            ),
            StallShape::Contended => info!(
                "    Shape: lock contention -- a large share of file state could not be read \
                 without blocking, so the shares above understate what was available"
            ),
            StallShape::Mixed => info!("    Shape: no single candidate dominates"),
        }
    }

    /// Generic merge for keyed chunks using `O(1)` key comparisons.
    ///
    /// This is the unified merge function that works with any `RawSortKey` type.
    /// It provides `O(1)` comparisons during merge for fixed-size keys (coordinate, template)
    /// and `O(name_len)` for variable-size keys (queryname).
    #[allow(clippy::too_many_lines)]
    fn merge_chunks_generic<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_chunks: MemorySources<K>,
        header: &Header,
        output: &Path,
        total_records: u64,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<u64> {
        use crate::loser_tree::LoserTree;
        use crate::pooled_bam_writer::PooledBamWriter;

        let (mut sources, mut guard) =
            Self::setup_phase2_merge::<K>(chunk_files, memory_chunks, pool)?;

        let output_header = self.create_output_header(header);

        // Initialize loser tree with first record from each source.
        // Each source owns its current-record state internally; no
        // intermediate `records: Vec<Vec<u8>>` is needed.
        let mut initial_keys: Vec<K> = Vec::with_capacity(sources.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(sources.len());

        for (idx, source) in sources.iter_mut().enumerate() {
            if let Some(key) = source.advance(guard.consumer_mut())? {
                initial_keys.push(key);
                source_map.push(idx);
            }
        }

        if initial_keys.is_empty() {
            debug!("Merge complete: 0 records merged");
            guard.finish_output(|| {
                PooledBamWriter::new(Arc::clone(pool), output, &output_header)?.finish()
            })?;
            return Ok(0);
        }

        let mut tree = LoserTree::new(initial_keys);

        debug!("Merge thread budget: {} pool workers + 1 I/O + 1 main (N+2)", pool.num_workers());
        let mut writer = PooledBamWriter::new(Arc::clone(pool), output, &output_header)?;

        let mut records_merged = 0u64;
        let merge_progress = ProgressTracker::new("Merged records")
            .with_interval(1_000_000)
            .with_total(total_records);

        let mut merge_probe = MergeProbe::new();

        // Sub-phase timing. Sampled 1-in-`merge_sample_interval`, so on a
        // billion-record merge this is a few million `Instant::now()` pairs
        // rather than two billion -- cheap enough to run unconditionally. It
        // used to be gated on debug logging, which meant the one breakdown that
        // explains where a merge-dominated sort spends its time was absent from
        // every benchmark log we actually collect.
        // Prime, and deliberately NOT 1024. Records are ~100 bytes and spill
        // blocks ~64 KB, so a power-of-two interval aliases against block
        // boundaries: the sampled record keeps landing on the refill that blocks
        // waiting for the next block, and scaling those up overestimated the
        // consumer's cost by ~2x -- enough that the sub-phases summed to more
        // than the loop wall clock they are supposed to partition. A prime
        // interval decorrelates the sample from any block-size-derived period.
        let merge_sample_interval: u64 = 1021;
        let mut merge_write_secs = 0.0f64;
        let mut merge_read_secs = 0.0f64;
        let mut merge_tree_secs = 0.0f64;
        let mut samples_taken: u64 = 0;
        // Countdown rather than `records_merged % interval`. A non-power-of-two
        // modulo compiles to a multiply-shift, which is a few cycles per record
        // -- ~0.1% of a billion-record merge. A decrement and a
        // perfectly-predicted branch is not measurable at all, which is what
        // lets this run unconditionally instead of behind a flag nobody
        // remembers to set.
        let mut sample_countdown: u64 = 0;
        // Seeding the loser tree pulled one block per source, and those pulls
        // are the likeliest of the whole merge to park -- every file is cold.
        // They precede this clock, so they must precede the stall counters too,
        // or `park_fraction` covers a longer interval than the wall time it is
        // divided by and `classify_stall` reads a merge that never stalled as
        // one that did.
        if let Some(consumer) = guard.consumer_mut() {
            consumer.restart_stalls();
        }
        // Snapshot the process-wide presentation counters so the sub-phase log
        // reports only this merge's delta, not totals accumulated by any prior
        // (sequential or concurrent) sort sharing the process.
        let borrowed_before = RECORD_BORROWED.load(std::sync::atomic::Ordering::Relaxed);
        let reassembled_before = RECORD_REASSEMBLED.load(std::sync::atomic::Ordering::Relaxed);
        let loop_start = Instant::now();

        while tree.winner_is_active() {
            let winner = tree.winner();
            let src_idx = source_map[winner];
            let sample_this = sample_countdown == 0;
            if sample_this {
                sample_countdown = merge_sample_interval - 1;
            } else {
                sample_countdown -= 1;
            }

            let record_bytes = winner_record_bytes(&sources[src_idx], guard.consumer_ref())?;
            if sample_this {
                let t0 = Instant::now();
                writer.write_raw_record(record_bytes)?;
                merge_write_secs += t0.elapsed().as_secs_f64();
            } else {
                writer.write_raw_record(record_bytes)?;
            }

            records_merged += 1;
            merge_progress.log_if_needed(1);

            if merge_probe.should_sample(records_merged) {
                let depths = pool.phase1_queue_depths();
                let consumer_stats = guard.consumer_mut().map(|c| c.probe_consumer_stats());
                merge_probe.log_mid_with_depths(depths, consumer_stats);
            }

            if sample_this {
                let t0 = Instant::now();
                let next = sources[src_idx].advance(guard.consumer_mut())?;
                merge_read_secs += t0.elapsed().as_secs_f64();

                let t0 = Instant::now();
                if let Some(key) = next {
                    tree.replace_winner(key);
                } else {
                    tree.remove_winner();
                }
                merge_tree_secs += t0.elapsed().as_secs_f64();
                samples_taken += 1;
            } else {
                let next = sources[src_idx].advance(guard.consumer_mut())?;
                if let Some(key) = next {
                    tree.replace_winner(key);
                } else {
                    tree.remove_winner();
                }
            }
        }

        let loop_total = loop_start.elapsed().as_secs_f64();
        let borrowed_this_merge =
            RECORD_BORROWED.load(std::sync::atomic::Ordering::Relaxed) - borrowed_before;
        let reassembled_this_merge =
            RECORD_REASSEMBLED.load(std::sync::atomic::Ordering::Relaxed) - reassembled_before;
        // The active limit, not the pool width: Phase 2 caps the pool to
        // `phase2_threads`, so on a run with a wider Phase 1 the extra threads
        // cannot take merge work and must not sit in the utilization
        // denominator -- a saturated pool would read as idle, and the verdict
        // would call a CPU-bound merge I/O-bound. Read while Phase 2 is still
        // active, so the number cannot depend on what teardown does to the cap.
        let active_workers = pool.active_workers();

        // Output backpressure, harvested before finalize for the same reason as
        // the stall report below: `finish` drains the output queue, and permits
        // taken during that drain belong to the drain, not to the merge loop.
        //
        // This is the merge's *second* consumer stall, and until now the only
        // unmeasured one. The consumer blocks here once per output block when the
        // compressors are behind, which the sampled breakdown charged to "enqueue
        // write" -- a bucket documented as a handoff that excludes compression.
        let (write_blocked_secs, write_blocked_waits) = writer.write_backpressure();
        // Retain the permit pool so the writer histograms can be harvested *after*
        // `finish` drains the output queue: the pool outlives the writer, and the
        // block writes and reorder waits incurred during that drain belong in the
        // "write block" and "write reorder wait" rows. A snapshot taken here,
        // before the drain, would omit the tail -- the same reason the stall
        // report below is harvested pre-finalize but the histograms are not.
        let writer_permit_pool = writer.permit_pool();

        // Harvest the consumer's stall report before finalizing: `finish_output`
        // releases the merge sources and with them the consumer, and the report
        // describes the loop that has just ended rather than the output drain
        // that follows.
        let stalls = {
            // Close the run in progress first, or the last (and often longest)
            // stretch on one source is dropped from the histogram.
            if let Some(consumer) = guard.consumer_mut() {
                consumer.finish_source_run();
            }
            guard.consumer_ref().map(MainThreadChunkConsumer::stall_report)
        };

        // Finalize before logging. `finish` drains the output queue, and every
        // block still in it is compressed by the same workers this breakdown
        // reports -- so a snapshot taken first omits the tail of
        // `output_compress` and divides the rest by a wall clock that stops
        // before the drain.
        guard.finish_output(|| writer.finish())?;

        // Harvest now, after the drain, from the retained pool. Empty only if the
        // writer was already finalized when the pool was retained, which cannot
        // happen on this path.
        let writer_stats = writer_permit_pool
            .as_deref()
            .map_or_else(Default::default, crate::worker_pool::PermitPool::writer_stats);

        Self::log_merge_sub_phases(
            (loop_total, loop_start.elapsed().as_secs_f64()),
            (merge_write_secs, merge_read_secs, merge_tree_secs),
            (samples_taken, records_merged),
            active_workers,
            pool,
            stalls,
            MergeConsumerDiag {
                backpressure_secs: write_blocked_secs,
                backpressure_waits: write_blocked_waits,
                borrowed: borrowed_this_merge,
                reassembled: reassembled_this_merge,
            },
        );
        Self::log_stage_latency(pool, &writer_stats);

        merge_progress.log_final();
        log_snapshot("phase2.end", 0);

        Ok(records_merged)
    }

    /// Merge keyed chunks with BAM index generation.
    ///
    /// Identical input/output pipeline to `merge_chunks_generic` — the shared
    /// worker pool decompresses the input runs (`PoolDisk` sources) and
    /// compresses the output blocks — but the output goes through
    /// [`PooledBamWriter::new_indexing`], which tracks each record's virtual
    /// file offset and returns the generated BAI index. A single pool of
    /// workers therefore serves both decompression and compression; there is no
    /// separate writer thread pool and no serialized single-reader input.
    fn merge_chunks_with_index<K: RawSortKey + Default + 'static>(
        &self,
        chunk_files: &[PathBuf],
        memory_chunks: MemorySources<K>,
        header: &Header,
        output: &Path,
        total_records: u64,
        pool: &Arc<SortWorkerPool>,
    ) -> Result<(noodles::bam::bai::Index, u64)> {
        use crate::loser_tree::LoserTree;
        use crate::pooled_bam_writer::PooledBamWriter;

        let (mut sources, mut guard) =
            Self::setup_phase2_merge::<K>(chunk_files, memory_chunks, pool)?;

        let output_header = self.create_output_header(header);

        // Initialize loser tree with first record from each source.
        let mut initial_keys: Vec<K> = Vec::with_capacity(sources.len());
        let mut source_map: Vec<usize> = Vec::with_capacity(sources.len());

        for (idx, source) in sources.iter_mut().enumerate() {
            if let Some(key) = source.advance(guard.consumer_mut())? {
                initial_keys.push(key);
                source_map.push(idx);
            }
        }

        if initial_keys.is_empty() {
            let index = guard.finish_output(|| {
                PooledBamWriter::new_indexing(Arc::clone(pool), output, &output_header)?
                    .finish_index()
            })?;
            debug!("Merge complete: 0 records merged");
            return Ok((index, 0));
        }

        let mut tree = LoserTree::new(initial_keys);

        let mut writer = PooledBamWriter::new_indexing(Arc::clone(pool), output, &output_header)?;
        let merge_progress = ProgressTracker::new("Merged records")
            .with_interval(1_000_000)
            .with_total(total_records);

        // Seeding the loser tree pulled one block per source, and those pulls
        // are the likeliest of the whole merge to park -- every file is cold.
        // They precede this clock, so they must precede the stall counters too,
        // or `park_fraction` covers a longer interval than the wall time it is
        // divided by and `classify_stall` reads a merge that never stalled as
        // one that did.
        if let Some(consumer) = guard.consumer_mut() {
            consumer.restart_stalls();
        }
        let loop_start = Instant::now();
        let mut records_merged: u64 = 0;
        while tree.winner_is_active() {
            let winner = tree.winner();
            let src_idx = source_map[winner];
            let record_bytes = winner_record_bytes(&sources[src_idx], guard.consumer_ref())?;
            writer.write_raw_record(record_bytes)?;
            records_merged += 1;
            merge_progress.log_if_needed(1);

            if let Some(key) = sources[src_idx].advance(guard.consumer_mut())? {
                tree.replace_winner(key);
            } else {
                tree.remove_winner();
            }
        }

        // This path has no sampled sub-phase breakdown -- adding the per-record
        // sampling here would duplicate the merge loop rather than share it --
        // but the stall counters need no sampling at all, so an indexed sort is
        // not left with nothing. `--write-index` is the default for a
        // coordinate sort, so that gap covers most production merges.
        //
        // The active cap is read while Phase 2 is still active, so teardown
        // cannot change it, and the consumer's report is harvested before
        // `finish_output` releases the merge sources and with them the
        // consumer. Both describe the loop that has just ended.
        let loop_total = loop_start.elapsed().as_secs_f64();
        let active_workers = pool.active_workers();
        let stalls = {
            // Close the run in progress first, or the last (and often longest)
            // stretch on one source is dropped from the histogram.
            if let Some(consumer) = guard.consumer_mut() {
                consumer.finish_source_run();
            }
            guard.consumer_ref().map(MainThreadChunkConsumer::stall_report)
        };

        // Retain the permit pool so the writer histograms can be harvested *after*
        // `finish_index` drains the output queue, exactly as the generic path
        // does: the pool outlives the writer, and the block writes and reorder
        // waits incurred during that drain belong in the "write block" and
        // "write reorder wait" rows. A snapshot taken here, before the drain,
        // would omit the tail.
        let writer_permit_pool = writer.permit_pool();

        // Finalize before logging, for the reason the generic path does: `finish`
        // drains the output queue, and every block still in it is compressed by
        // the same workers `log_merge_stalls` divides into `merge_total`. Logging
        // first omits that tail of `output_compress` while still charging the
        // window it was queued in, understating utilization -- and utilization is
        // the input `classify_stall` uses to tell a saturated pool from a
        // scheduling defect, so the bias flips a verdict rather than shading a
        // number. `loop_total` stays the consumer's park-fraction denominator,
        // which is a fraction of the merge loop alone.
        let index = guard.finish_output(|| writer.finish_index())?;

        // Harvest now, after the drain, from the retained pool. Empty only if the
        // writer was already finalized when the pool was retained, which cannot
        // happen on this path.
        let writer_stats = writer_permit_pool
            .as_deref()
            .map_or_else(Default::default, crate::worker_pool::PermitPool::writer_stats);

        Self::log_merge_stalls(
            loop_total,
            loop_start.elapsed().as_secs_f64(),
            active_workers,
            stalls,
            pool,
        );
        // `--write-index` is the default for a coordinate sort, so wiring the
        // stage-latency table here -- matching `merge_chunks_generic` -- is what
        // gives most production merges any stage-latency and writer-histogram
        // rows at all.
        Self::log_stage_latency(pool, &writer_stats);

        merge_progress.log_final();
        Ok((index, records_merged))
    }

    /// Create output header with appropriate sort order tags.
    fn create_output_header(&self, header: &Header) -> Header {
        super::create_output_header(self.sort_order, header)
    }

    /// Create per-base temp directories and an allocator over their subdirs.
    ///
    /// For each user-supplied base directory, a fresh sort-run subdirectory is
    /// created (via `tempfile::TempDir`). The returned [`Vec<TempDir>`] owns
    /// those handles so subdirs are removed on drop; the allocator hands out
    /// the corresponding subdir paths for chunk/merged file placement.
    ///
    /// When `temp_dirs` is empty, a single subdirectory is created under the
    /// system default temp location.
    fn create_temp_dirs(&self) -> Result<(Vec<TempDir>, TmpDirAllocator)> {
        use super::create_temp_dir;

        if self.temp_dirs.is_empty() {
            let td = create_temp_dir(None)?;
            let base = td.path().to_path_buf();
            let alloc = TmpDirAllocator::new(vec![base])?;
            return Ok((vec![td], alloc));
        }

        let mut handles = Vec::with_capacity(self.temp_dirs.len());
        let mut subdirs = Vec::with_capacity(self.temp_dirs.len());
        for base in &self.temp_dirs {
            let td = create_temp_dir(Some(base))?;
            subdirs.push(td.path().to_path_buf());
            handles.push(td);
        }
        let alloc = TmpDirAllocator::new(subdirs)?;
        Ok((handles, alloc))
    }
}

/// Extract a packed `TemplateKey` directly from BAM record bytes.
///
/// This function computes the template-coordinate sort key inline, avoiding
/// heap allocations for the read name by using a hash instead.
///
/// When `cell_tag` is `Some`, the CB (cellular barcode) tag value is hashed
/// and included in the sort key between neg2 and MI, matching fgbio's order.
#[must_use]
pub fn extract_template_key_inline(
    bam_bytes: &[u8],
    lib_lookup: &LibraryLookup,
    cell_tag: Option<SamTag>,
    cb_hasher: &ahash::RandomState,
) -> TemplateKey {
    use fgumi_raw_bam;
    use fgumi_raw_bam::{flags, mate_unclipped_5prime, unclipped_5prime_raw};

    // Single-pass extraction of all aux tags (MI, RG, cell barcode, MC)
    let aux = fgumi_raw_bam::extract_template_aux_tags(bam_bytes, cell_tag);
    let mi = aux.mi;
    let library = lib_lookup.ordinal_from_rg(aux.rg);
    let cb_hash = aux.cell.map_or(0u64, |cb_bytes| cb_hasher.hash_one(cb_bytes));

    // Extract fields from raw bytes
    let v = fgumi_raw_bam::RawRecordView::new(bam_bytes);
    let tid = v.ref_id();
    let pos = v.pos();
    let l_read_name = v.l_read_name() as usize;
    let flag = v.flags();
    let mate_tid = v.mate_ref_id();
    let mate_pos = v.mate_pos();

    // Extract flags
    let is_unmapped = (flag & flags::UNMAPPED) != 0;
    let mate_unmapped = (flag & flags::MATE_UNMAPPED) != 0;
    let is_reverse = (flag & flags::REVERSE) != 0;
    let mate_reverse = (flag & flags::MATE_REVERSE) != 0;
    let is_paired = (flag & flags::PAIRED) != 0;

    // Hash read name (exclude null terminator)
    let name_len = l_read_name.saturating_sub(1);
    let name = if name_len > 0 && 32 + name_len <= bam_bytes.len() {
        &bam_bytes[32..32 + name_len]
    } else {
        &[]
    };
    let name_hash = lib_lookup.hash_name(name);

    // Secondary/supplementary reads cannot reconstruct their template coordinate
    // from their own record (they carry their own and their mate's position, but
    // not their own primary's). When `fgumi zipper` has stamped the exact
    // template coordinate into the `tc` tag, key on it so the read sorts at its
    // primary pair's coordinate — a placement per-record keying (samtools/fgbio)
    // cannot achieve. Falls through to the position-based computation when absent.
    let is_secondary = (flag & flags::SECONDARY) != 0;
    let is_supplementary = (flag & flags::SUPPLEMENTARY) != 0;
    if is_secondary || is_supplementary {
        let aux_slice = fgumi_raw_bam::aux_data_slice(bam_bytes);
        if let Some([tid1, pos1, neg1, tid2, pos2, neg2]) =
            fgumi_raw_bam::read_tc_template_coordinate(aux_slice)
        {
            // `tc` is canonicalized by (tid, pos) only (see zipper's
            // `add_template_coordinate_tags_raw`), whereas the per-record path
            // additionally tiebreaks on strand for equal positions. In the
            // degenerate case where both primaries share an identical tid and
            // unclipped-5' position but differ in strand, the supplementary's
            // secondary-lane strand order can differ from its primaries'. This is
            // harmless: the primary lane (tid1, tid2, pos1) still matches, so the
            // read co-locates with its template; only a same-coordinate tiebreak
            // differs.
            //
            // `TemplateKey::new`'s final argument is the `is_upper` bit, packed
            // into `name_hash_upper` as the *last* tiebreak (after position,
            // strand, CB, library, MI, and name). For a mapped read it answers
            // "is this read's own primary the upper (lane-2) end?", derived from
            // this read's position vs. its mate. A secondary/supplementary read
            // cannot recover that: `tc` is canonicalized by (tid, pos) and does
            // not record which segment (R1/R2) landed on which lane, and the
            // read's own alignment position is exactly the misleading coordinate
            // this branch exists to avoid. We use `is_read2` because it mirrors
            // the primary's `is_upper` for the standard FR pair (R1 is the lower
            // lane -> is_upper=false; R2 is the upper lane -> is_upper=true) and
            // is fully deterministic; see
            // `test_extract_template_key_secondary_supplementary_matches_primary_via_tc`,
            // which asserts the tc-keyed key is byte-identical to the primary's.
            // Any residual disagreement (an exotic pair whose R2 is the lower
            // lane) only reorders same-name records at an identical coordinate,
            // so it cannot affect co-location or dedup grouping.
            let is_read2 = (flag & 0x80) != 0;
            return TemplateKey::new(
                tid1,
                pos1,
                neg1 != 0,
                tid2,
                pos2,
                neg2 != 0,
                cb_hash,
                library,
                mi,
                name_hash,
                is_read2,
            );
        }
    }

    // Handle unmapped reads
    if is_unmapped {
        if is_paired && !mate_unmapped {
            // Unmapped read with mapped mate - use mate's position as primary key
            let mate_unclipped =
                aux.mc.map_or(mate_pos, |mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc));

            return TemplateKey::new(
                mate_tid,
                mate_unclipped,
                mate_reverse,
                i32::MAX,
                i32::MAX,
                false,
                cb_hash,
                library,
                mi,
                name_hash,
                true, // Unmapped read is always "upper" relative to mapped mate
            );
        }

        // Completely unmapped - sort to end. Still carry the read's library and
        // MI lanes so a fully-unmapped read realizes the same tertiary as its
        // mapped, same-library peers (otherwise the dropped-lane verify treats a
        // single-library file as if the library varied; see #375).
        let is_read2 = (flag & 0x80) != 0; // is_last_segment flag
        return TemplateKey::unmapped(name_hash, cb_hash, library, mi, is_read2);
    }

    // Calculate unclipped 5' position for this read (zero-alloc: reads cigar directly)
    let this_pos = unclipped_5prime_raw(bam_bytes, pos, is_reverse);

    // Calculate mate's unclipped 5' position
    let mate_unclipped = if is_paired && !mate_unmapped {
        aux.mc.map_or(mate_pos, |mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc))
    } else {
        mate_pos
    };

    // Determine canonical ordering
    let (tid1, tid2, pos1, pos2, neg1, neg2, is_upper) = if is_paired && !mate_unmapped {
        // Samtools logic: is_upper if pos > mate_pos, or (pos == mate_pos && this read is reverse)
        let is_upper = (tid, this_pos) > (mate_tid, mate_unclipped)
            || ((tid, this_pos) == (mate_tid, mate_unclipped) && is_reverse);

        if is_upper {
            // Swap: mate's position comes first
            (mate_tid, tid, mate_unclipped, this_pos, mate_reverse, is_reverse, true)
        } else {
            // No swap: this read's position comes first
            (tid, mate_tid, this_pos, mate_unclipped, is_reverse, mate_reverse, false)
        }
    } else {
        // Unpaired or mate unmapped - use MAX for tid2/pos2
        (tid, i32::MAX, this_pos, i32::MAX, is_reverse, false, false)
    };

    TemplateKey::new(tid1, pos1, neg1, tid2, pos2, neg2, cb_hash, library, mi, name_hash, is_upper)
}

// SortStats is defined in crate root (lib.rs); alias it here for internal use.
pub(crate) use crate::SortStats as RawSortStats;

// ============================================================================
// Pool-integrated BAM reader construction
//
// This function depends on sort-internal types (`SortWorkerPool`,
// `PooledInputStream`, `phase`), which is why it lives here in `fgumi-sort`
// rather than in `fgumi-bam-io`.
// ============================================================================

/// Create a raw BAM reader using the pool's Phase 1 integrated reading.
///
/// Workers in the pool do `ReadInputBlocks` + `DecompressInput`. The main
/// thread consumes decompressed bytes via `PooledInputStream`.
///
/// When `async_reader` is false, no extra threads are spawned: the pool's
/// block reader reads directly from the input file. When `async_reader` is
/// true, the input file is wrapped in a `PrefetchReader`, which spawns one
/// dedicated OS thread (`fgumi-prefetch`) that reads raw bytes ahead into a
/// bounded queue so the pool's block reader never blocks on disk I/O.
///
/// # Flow
///
/// 1. Open the input once, transcoding it if it is uncompressed SAM
/// 2. Parse the header, capturing the bytes it consumed so they can be replayed
/// 3. Set the replaying stream as the pool's input file
/// 4. Set pool to PHASE1 — workers start reading/decompressing
/// 5. Main thread skips header bytes from decompressed stream
/// 6. Returns `RawBamReader<PooledInputStream>` for direct record iteration
///
/// # Errors
///
/// Returns an error if the input cannot be opened, the header cannot be
/// parsed, or header bytes cannot be skipped from the decompressed stream.
fn create_raw_bam_reader_pool_integrated<P: AsRef<Path>>(
    path: P,
    pool: &Arc<SortWorkerPool>,
    async_reader: bool,
) -> Result<(fgumi_raw_bam::RawBamReader<PooledInputStream>, Header)> {
    use crate::worker_pool::phase;
    use std::io;

    let path_ref = path.as_ref();

    let opened: Box<dyn io::Read + Send> = if is_stdin_path(path_ref) {
        if async_reader {
            // `--async-reader` is about decoupling the read from the block
            // reader, which stdin needs at least as much as a file does: the
            // prefetch thread also subsumes the buffering below, reading ahead
            // in chunks into a bounded queue.
            log::debug!("async sort reader enabled: spawning fgumi-prefetch thread for stdin");
            Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::new(io::stdin()))
        } else {
            // `io::Stdin` re-acquires a mutex and reads through an 8 KiB buffer
            // on every call; the pool's block reader wants far bigger gulps than
            // that, so give the stdin path the same 2 MiB buffer the file path
            // gets.
            Box::new(io::BufReader::with_capacity(SORT_INPUT_BUFFER_SIZE, io::stdin()))
        }
    } else {
        let file = std::fs::File::open(path_ref)
            .with_context(|| format!("Failed to open input BAM: {}", path_ref.display()))?;

        // Grow the per-fd readahead window. This is the plain sequential hint
        // and applies however the bytes are subsequently read; the WILLNEED
        // hints that `PrefetchReader` issues are a separate, async-only extra.
        fgumi_bam_io::os_hints::advise_sequential(&file);

        if async_reader {
            log::debug!(
                "async sort reader enabled: spawning fgumi-prefetch thread for {}",
                path_ref.display()
            );
            Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::from_file(file))
        } else {
            Box::new(io::BufReader::with_capacity(SORT_INPUT_BUFFER_SIZE, file))
        }
    };

    // Uncompressed SAM becomes a BGZF stream here, so the worker pool below
    // only ever decompresses BGZF blocks.
    let opened = fgumi_bam_io::sam_input::normalize_to_bgzf(opened, path_ref)?;

    // Parse the header through a tee and replay the bytes it consumed, rather
    // than rewinding: the input is opened exactly once, which keeps stdin,
    // FIFOs and transcoded SAM on the same path as a plain file.
    let (header, reader) = fgumi_bam_io::read_header_and_replay(opened, path_ref)?;

    pool.set_input_file(reader);
    pool.set_phase(phase::PHASE1);

    let mut pooled_input = PooledInputStream::new(
        pool.decompressed_input_queue(),
        pool.decompressed_input_done_flag(),
        pool.input_read_error_flag(),
        pool.decompress_error_flag(),
    );

    // Deliberately not phrased as a header failure. The header was parsed
    // successfully above; by the time the main thread skips its bytes here the
    // pool's workers have already read and decompressed well past it, so a
    // fault surfacing at this point is usually a *record* fault (a malformed
    // SAM line, a truncated block) rather than anything to do with the header.
    // The worker logs the specific cause before setting the flag this reads.
    skip_bam_header(&mut pooled_input)
        .with_context(|| format!("Failed to read input: {}", path_ref.display()))?;

    let raw_reader = fgumi_raw_bam::RawBamReader::new(pooled_input);
    Ok((raw_reader, header))
}

/// Skip the BAM header from a reader positioned at the start of a BAM stream.
///
/// Reads and discards: magic (4 bytes), header text length + text,
/// `n_ref` + reference entries.
fn skip_bam_header<R: Read>(reader: &mut R) -> Result<()> {
    use std::io;
    let mut buf4 = [0u8; 4];

    reader.read_exact(&mut buf4)?;
    anyhow::ensure!(&buf4 == b"BAM\x01", "Not a BAM file (bad magic)");

    reader.read_exact(&mut buf4)?;
    let l_text = u32::from_le_bytes(buf4) as usize;
    let copied = io::copy(&mut reader.take(l_text as u64), &mut io::sink())?;
    anyhow::ensure!(
        copied == l_text as u64,
        "BAM header text truncated: expected {l_text} bytes, got {copied}"
    );

    reader.read_exact(&mut buf4)?;
    let n_ref = u32::from_le_bytes(buf4) as usize;

    for _ in 0..n_ref {
        reader.read_exact(&mut buf4)?;
        let l_name = u32::from_le_bytes(buf4) as usize;
        let copied = io::copy(&mut reader.take(l_name as u64), &mut io::sink())?;
        anyhow::ensure!(
            copied == l_name as u64,
            "BAM reference name truncated: expected {l_name} bytes, got {copied}"
        );
        reader.read_exact(&mut buf4)?;
    }

    Ok(())
}

// ============================================================================
// Template-key variant selection
// ============================================================================

/// Which optional lanes a chosen template key retains.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct TemplateKeyVariant {
    /// Retain the `cb_hash` lane.
    pub cb: bool,
    /// Retain the tertiary (library|mi) lane.
    pub tertiary: bool,
}

impl TemplateKeyVariant {
    /// Number of u64 lanes: 3 + cb + tertiary.
    #[must_use]
    pub fn lanes(self) -> usize {
        3 + usize::from(self.cb) + usize::from(self.tertiary)
    }
}

/// Parsed `--key-types` spec controlling which optional sort-key lanes are kept.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum KeyTypesSpec {
    /// Auto-detect from the first record (default).
    #[default]
    Auto,
    /// Force all lanes (N=5).
    Full,
    /// Force no optional lanes (N=3).
    None,
    /// Force a specific optional set.
    Explicit {
        /// Retain the `cb_hash` lane.
        cb: bool,
        /// Retain the tertiary (library|mi) lane.
        tertiary: bool,
    },
}

/// Mask for the MI component of the `tertiary` key lane (low 48 bits).
/// `tertiary` packs `library_ordinal << 48 | (mi_value << 1 | !mi_suffix)`, so the
/// high 16 bits are the library ordinal and the low 48 bits encode MI. A nonzero
/// low-48 region means the first record carries an MI tag.
const TERTIARY_MI_MASK: u64 = (1u64 << 48) - 1;

/// Choose the template-key variant for this sort.
///
/// `first_key` is the full key extracted from the first record (`None` for an
/// empty input — then the variant is irrelevant; default to lite). `spec` is the
/// parsed `--key-types` override. `header_library_varies` is `true` when the
/// header realizes more than one distinct library ordinal (see
/// [`LibraryLookup::distinct_header_ordinals`]).
///
/// Selection only *provisions* the variant; the decode-time verify (fill loop)
/// still *guarantees* the dropped lanes are constant. Under `Auto`, the tertiary
/// library lane is kept only when it can actually vary: the first record carries
/// an MI tag, or the header realizes more than one library ordinal. A single,
/// constant library ordinal in the high bits is droppable — the header informs
/// provisioning, while verify backstops the rare case of a record whose RG is
/// absent from the header (which realizes ordinal 0).
#[must_use]
pub fn select_template_variant(
    first_key: Option<&TemplateKey>,
    spec: KeyTypesSpec,
    header_library_varies: bool,
) -> TemplateKeyVariant {
    match spec {
        KeyTypesSpec::Full => TemplateKeyVariant { cb: true, tertiary: true },
        KeyTypesSpec::None => TemplateKeyVariant { cb: false, tertiary: false },
        KeyTypesSpec::Explicit { cb, tertiary } => TemplateKeyVariant { cb, tertiary },
        KeyTypesSpec::Auto => match first_key {
            None => TemplateKeyVariant { cb: false, tertiary: false },
            Some(k) => {
                // Keep the tertiary lane only if it can actually vary: MI present on
                // the first record (low-48 bits), or the header realizes more than
                // one library ordinal. A single (constant) library ordinal in the
                // high bits is droppable — decode-time verify backstops the rare
                // case of a record whose RG is absent from the header (ordinal 0).
                let mi_present = (k.tertiary & TERTIARY_MI_MASK) != 0;
                TemplateKeyVariant {
                    cb: k.cb_hash != 0,
                    tertiary: mi_present || header_library_varies,
                }
            }
        },
    }
}

// ============================================================================
// Dropped-lane violation detection
// ============================================================================

/// A dropped-lane violation discovered during decode-time verify.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DroppedLaneViolation {
    /// A record carried a CB tag absent from the first record.
    Cb,
    /// A record carried an MI value differing from the first record.
    Mi,
    /// A record's realized library ordinal differed from the first record.
    Library,
}

impl DroppedLaneViolation {
    /// The `--key-types` token that re-includes this lane.
    #[must_use]
    pub fn key_types_token(self) -> &'static str {
        match self {
            DroppedLaneViolation::Cb => "cb",
            DroppedLaneViolation::Mi => "mi",
            DroppedLaneViolation::Library => "library",
        }
    }
}

/// Verify that the lanes the chosen `variant` dropped are constant for `cur`
/// relative to the sort's `first` full key. Returns the first violation found,
/// or `None` if all dropped lanes match.
///
/// Tertiary sub-fields are decoded only to attribute the violation: the high 16
/// bits are the library ordinal; the low 48 bits are `(mi_value << 1) | !mi_suffix`
/// (`mi_value` occupies bits 47..1; `!mi_suffix` is bit 0).
#[must_use]
pub fn verify_dropped_lanes(
    first: &TemplateKey,
    cur: &TemplateKey,
    variant: TemplateKeyVariant,
) -> Option<DroppedLaneViolation> {
    if !variant.cb && cur.cb_hash != first.cb_hash {
        return Some(DroppedLaneViolation::Cb);
    }
    if !variant.tertiary && cur.tertiary != first.tertiary {
        // Attribute: library is the high 16 bits, mi the low 48.
        let lib_changed = (cur.tertiary >> 48) != (first.tertiary >> 48);
        return Some(if lib_changed {
            DroppedLaneViolation::Library
        } else {
            DroppedLaneViolation::Mi
        });
    }
    None
}

/// Build the actionable error message for a dropped-lane violation.
fn dropped_lane_error(name: &str, v: DroppedLaneViolation) -> anyhow::Error {
    let field = match v {
        DroppedLaneViolation::Cb => "CB",
        DroppedLaneViolation::Mi => "MI",
        DroppedLaneViolation::Library => "library",
    };
    anyhow::anyhow!(
        "record {name} carries a {field} value absent from the input's first record; \
         re-run with --key-types {}",
        v.key_types_token(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use fgumi_raw_bam::flags;
    use fgumi_sam::{PairBuilder, SamBuilder};
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;
    use rstest::rstest;

    // ========================================================================
    // Record-count invariant
    // ========================================================================

    /// A failed permutation check must not leave a readable BAM behind.
    ///
    /// The output is finalized by the time the check runs, so it is a valid,
    /// short BAM -- the most dangerous kind of failure output, because every
    /// existence check and every reader accepts it. The rest of the crate
    /// deletes on abort (`test_sort_records_producer_error_aborts`), and this
    /// path has to as well or the error's own claim is false.
    #[rstest]
    #[case::index_written(true)]
    #[case::no_index(false)]
    fn record_count_mismatch_removes_the_output(#[case] write_index: bool) {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let output = dir.path().join("out.bam");
        let bai = fgumi_bam_io::bai_sidecar_path(&output);
        std::fs::write(&output, b"a short but perfectly valid BAM").expect("write output");
        std::fs::write(&bai, b"an index describing it").expect("write index");

        let stats = RawSortStats { total_records: 100, output_records: 99, runs_written: 1 };
        let err = enforce_record_count(&stats, &output, write_index)
            .expect_err("a record-count mismatch must fail the sort");
        let msg = format!("{err:#}");

        assert!(!output.exists(), "the incomplete output must not survive the error");
        assert!(msg.contains("read 100 but wrote 99"), "got: {msg}");
        assert!(msg.contains("differ by 1"), "got: {msg}");
        assert!(
            msg.contains("has been removed"),
            "the message must state what was done, not just what not to do: {msg}"
        );
        // The index describes a file that no longer exists, so it goes too --
        // but only when this sort is the one that wrote it.
        assert_eq!(
            bai.exists(),
            !write_index,
            "index removal must follow write_index (write_index={write_index})"
        );
    }

    /// stdout cannot be un-written, and `/dev/stdout` must certainly not be
    /// unlinked. The error still fires; it just cannot claim a removal.
    #[rstest]
    #[case::dash("-")]
    #[case::dev_stdout("/dev/stdout")]
    fn record_count_mismatch_never_unlinks_stdout(#[case] output: &str) {
        let stats = RawSortStats { total_records: 10, output_records: 4, runs_written: 1 };
        let err = enforce_record_count(&stats, Path::new(output), false)
            .expect_err("a record-count mismatch must fail the sort");
        let msg = format!("{err:#}");

        assert!(Path::new("/dev/stdout").exists(), "/dev/stdout must survive");
        assert!(msg.contains("written to stdout"), "got: {msg}");
        assert!(
            !msg.contains("has been removed"),
            "nothing was removed, so the message must not say so: {msg}"
        );
    }

    /// A symlinked output is written *through* to the real file, so unlinking
    /// the link would leave the short BAM exactly where a reader would find it.
    #[test]
    fn record_count_mismatch_removes_the_symlink_target() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let real = dir.path().join("real.bam");
        let link = dir.path().join("link.bam");
        std::fs::write(&real, b"a short but perfectly valid BAM").expect("write output");
        std::os::unix::fs::symlink(&real, &link).expect("symlink");

        let stats = RawSortStats { total_records: 10, output_records: 9, runs_written: 1 };
        enforce_record_count(&stats, &link, false).expect_err("mismatch must fail the sort");

        assert!(!real.exists(), "the file actually written must be the one removed");
    }

    /// A FIFO was written through, not created, so it is not this sort's to
    /// unlink -- and the message must not claim otherwise.
    #[test]
    fn record_count_mismatch_leaves_a_non_regular_output_alone() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let fifo = dir.path().join("out.fifo");
        let mkfifo =
            std::process::Command::new("mkfifo").arg(&fifo).status().expect("failed to run mkfifo");
        assert!(mkfifo.success(), "mkfifo failed");

        let stats = RawSortStats { total_records: 10, output_records: 9, runs_written: 1 };
        let err = enforce_record_count(&stats, &fifo, false).expect_err("mismatch must fail");
        let msg = format!("{err:#}");

        assert!(fifo.exists(), "a FIFO must not be unlinked by the sort that wrote through it");
        assert!(msg.contains("could not be removed"), "got: {msg}");
    }

    /// The check must be silent and side-effect-free when the counts agree --
    /// otherwise every successful sort would delete its own output.
    #[test]
    fn matching_record_counts_leave_the_output_alone() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let output = dir.path().join("out.bam");
        std::fs::write(&output, b"the real output").expect("write output");

        let stats = RawSortStats { total_records: 100, output_records: 100, runs_written: 1 };
        enforce_record_count(&stats, &output, true).expect("matching counts must succeed");
        assert!(output.exists(), "a successful sort must keep its output");
    }

    /// An output that is already gone is the state the removal wanted, so it
    /// must not be reported as a removal failure on top of the count error.
    #[test]
    fn record_count_mismatch_tolerates_a_missing_output() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let output = dir.path().join("never-written.bam");

        let stats = RawSortStats { total_records: 5, output_records: 0, runs_written: 0 };
        let err = enforce_record_count(&stats, &output, true)
            .expect_err("a record-count mismatch must fail the sort");
        let msg = format!("{err:#}");
        assert!(msg.contains("has been removed"), "got: {msg}");
        assert!(!msg.contains("could not be removed"), "absent is not a failure: {msg}");
    }

    // ========================================================================
    // Merge diagnostics reporting
    // ========================================================================

    /// The block-lifecycle gate must name every report the block prints.
    ///
    /// Each case populates exactly one report and asserts the block still
    /// speaks. A report left out of the gate passes every case but its own,
    /// which is what makes one case per report the point rather than a
    /// formality -- `scans` was the one that had been omitted.
    #[rstest]
    #[case::nothing_recorded_stays_silent("none", true)]
    #[case::block_lifecycle_alone_speaks("life", false)]
    #[case::refill_cycle_alone_speaks("refill", false)]
    #[case::consumer_trace_alone_speaks("consumer", false)]
    #[case::fruitless_scans_alone_speak("scans", false)]
    fn block_lifecycle_gate_covers_every_report(#[case] populate: &str, #[case] silent: bool) {
        use crate::merge_trace::{
            BlockLifecycleStats, ConsumerTraceStats, DurationHistogram, EmptyCause, RefillStats,
        };

        let life = {
            let stats = BlockLifecycleStats::default();
            if populate == "life" {
                stats.raw_dwell.record(1_000);
            }
            let mut report = stats.snapshot();
            // `is_empty` keys off the working stages, which the pool fills in.
            if populate == "life" {
                report.decompress = {
                    let h = DurationHistogram::default();
                    h.record(1_000);
                    h.snapshot()
                };
            }
            report
        };
        let refill = {
            let stats = RefillStats::default();
            if populate == "refill" {
                stats.record_empty(EmptyCause::Dry);
            }
            stats.snapshot()
        };
        let consumer = {
            let stats = ConsumerTraceStats::default();
            if populate == "consumer" {
                stats.record_source_run(4);
            }
            stats.snapshot()
        };
        let scans = {
            let hist = DurationHistogram::default();
            if populate == "scans" {
                hist.record(1_000);
            }
            hist.snapshot()
        };

        assert_eq!(
            block_lifecycle_is_silent(&life, &refill, &consumer, &scans),
            silent,
            "gate disagreed with the only populated report ({populate})"
        );
    }

    /// The stall gate must name every report the stall block prints.
    ///
    /// Same shape, and same reason, as
    /// [`block_lifecycle_gate_covers_every_report`]: one case per report, each
    /// populating only that report, so a report dropped from the gate fails its
    /// own case rather than passing on a neighbour's.
    #[rstest]
    #[case::nothing_recorded_stays_silent("none", true)]
    #[case::consumer_stalls_alone_speak("stalls", false)]
    #[case::fruitless_scans_alone_speak("scans", false)]
    #[case::wake_latency_alone_speaks("wake", false)]
    fn merge_stall_gate_covers_every_report(#[case] populate: &str, #[case] silent: bool) {
        use crate::merge_stalls::{
            ConsumerStallTracker, Phase2ScanStats, Phase2ScanTally, Phase2Skip, WakeLatencyStats,
            WakePhase,
        };

        let stalls = {
            let mut tracker = ConsumerStallTracker::new(1);
            if populate == "stalls" {
                tracker.note_block_pull();
            }
            tracker.snapshot()
        };
        let scans = {
            let stats = Phase2ScanStats::default();
            if populate == "scans" {
                let mut tally = Phase2ScanTally::default();
                tally.note(Phase2Skip::RawFull);
                stats.record_fruitless_scan(tally);
            }
            stats.snapshot()
        };
        let wake = {
            let stats = WakeLatencyStats::default();
            if populate == "wake" {
                stats.record_sleep(WakePhase::Merge, 320, 1_000);
            }
            stats.snapshot(WakePhase::Merge)
        };

        // Mirrors the caller, which drops an all-zero stall report before the
        // gate sees it -- the tracker always exists, so `Some(empty)` is the
        // shape a merge with no stalls actually produces.
        let stalls = Some(stalls).filter(|s| !s.is_empty());
        assert_eq!(
            merge_stalls_are_silent(stalls.as_ref(), &scans, &wake),
            silent,
            "gate disagreed with the only populated report ({populate})"
        );
    }

    /// The stall gate must not decide whether the lifecycle block prints.
    ///
    /// The two are populated by opposite conditions: a merge that flowed without
    /// ever stalling records a full lifecycle trace and nothing at all in the
    /// three stall reports. Gating the lifecycle block on the stall gate --
    /// which `log_merge_stalls` did by early-returning -- drops the entire trace
    /// on exactly those runs.
    #[test]
    fn silent_stall_gate_does_not_silence_the_lifecycle_gate() {
        use crate::merge_stalls::{
            ConsumerStallTracker, Phase2ScanStats, WakeLatencyStats, WakePhase,
        };
        use crate::merge_trace::{
            BlockLifecycleStats, ConsumerTraceStats, DurationHistogram, RefillStats,
        };

        // A merge that never stalled: no block pull waited, no scan came back
        // empty, no worker slept.
        let stalls = ConsumerStallTracker::new(1).snapshot();
        let scans = Phase2ScanStats::default().snapshot();
        let wake = WakeLatencyStats::default().snapshot(WakePhase::Merge);
        assert!(
            merge_stalls_are_silent(Some(&stalls), &scans, &wake),
            "a merge with no stalls should leave the stall block silent"
        );

        // The same merge still moved every block through the pipeline.
        let life = {
            let stats = BlockLifecycleStats::default();
            stats.raw_dwell.record(1_000);
            let mut report = stats.snapshot();
            report.decompress = {
                let h = DurationHistogram::default();
                h.record(1_000);
                h.snapshot()
            };
            report
        };
        assert!(
            !block_lifecycle_is_silent(
                &life,
                &RefillStats::default().snapshot(),
                &ConsumerTraceStats::default().snapshot(),
                &DurationHistogram::default().snapshot(),
            ),
            "the lifecycle block has a full trace to print and must not be gated on the stalls"
        );
    }

    // ========================================================================
    // Chunk spill triggers
    // ========================================================================

    /// The unstable chunk sort is only correct because each key carries a unique
    /// ingest position, and that position is a `u32`. `MAX_CHUNK_RECORDS` is
    /// what keeps the stamp site's `usize` -> `u32` narrowing lossless, so it
    /// must never exceed what a `u32` can hold.
    ///
    /// The stamp happens *before* the push, so the largest position a chunk can
    /// produce is `MAX_CHUNK_RECORDS - 1`.
    #[test]
    fn test_max_chunk_records_keeps_ingest_positions_in_u32() {
        let largest_position = MAX_CHUNK_RECORDS - 1;

        assert!(
            u32::try_from(largest_position).is_ok(),
            "a full chunk would stamp position {largest_position}, which must fit in u32",
        );
    }

    /// Memory pressure and the record cap are independent spill triggers: either
    /// one alone must force the spill.
    #[rstest]
    #[case::empty(0, 1_000, 0, false)]
    #[case::under_both(999, 1_000, 10, false)]
    #[case::at_memory_limit(1_000, 1_000, 10, true)]
    #[case::over_memory_limit(1_001, 1_000, 10, true)]
    #[case::one_under_record_cap(0, usize::MAX, MAX_CHUNK_RECORDS - 1, false)]
    #[case::at_record_cap(0, usize::MAX, MAX_CHUNK_RECORDS, true)]
    // `MAX_CHUNK_RECORDS` is `u32::MAX`, which on a 32-bit target is `usize::MAX`
    // — `+ 1` would not be representable there. Saturating leaves the case at the
    // cap on such a target, which is still a spill.
    #[case::over_record_cap(0, usize::MAX, MAX_CHUNK_RECORDS.saturating_add(1), true)]
    fn test_should_spill_chunk(
        #[case] memory_used: usize,
        #[case] memory_limit: usize,
        #[case] chunk_records: usize,
        #[case] expected: bool,
    ) {
        assert_eq!(should_spill_chunk(memory_used, memory_limit, chunk_records), expected);
    }

    // ========================================================================
    // process_umask concurrency
    // ========================================================================

    /// `process_umask` reads the mask via a non-atomic set-0-then-restore. Two
    /// interleaved probes can leave the process mask permanently `0` (and every
    /// probe observe the wrong value); `UMASK_LOCK` must serialize them. Read the
    /// mask once, hammer it from many threads, and assert every probe — and the
    /// final mask — matches the initial value. Without the lock this corrupts to
    /// `0` reliably under contention. (Uses only `process_umask` so it introduces
    /// no new `unsafe` site.)
    #[cfg(unix)]
    #[test]
    fn test_process_umask_is_concurrency_safe() {
        use std::thread;

        let expected = process_umask();
        let threads: Vec<_> = (0..8)
            .map(|_| {
                thread::spawn(move || {
                    for _ in 0..2000 {
                        assert_eq!(
                            process_umask(),
                            expected,
                            "a concurrent probe observed a corrupted umask"
                        );
                    }
                })
            })
            .collect();
        for t in threads {
            t.join().expect("umask probe thread panicked");
        }
        assert_eq!(
            process_umask(),
            expected,
            "concurrent probes must leave the process umask unchanged"
        );
    }

    // ========================================================================
    // resolve_symlink_output / target_file_mode
    // ========================================================================

    /// A non-symlink path — even one that does not exist yet — is returned
    /// unchanged, so a brand-new merge output is not perturbed by the resolver.
    #[cfg(unix)]
    #[test]
    fn test_resolve_symlink_output_passthrough_for_nonexistent() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("merged.bam");
        assert_eq!(resolve_symlink_output(&out).unwrap(), out);
    }

    /// A relative symlink target resolves against the link's own directory, so a
    /// `latest.bam -> run.bam` link is followed to the sibling it points at.
    #[cfg(unix)]
    #[test]
    fn test_resolve_symlink_output_follows_relative_target() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("run.bam");
        std::fs::write(&target, b"x").unwrap();
        let link = dir.path().join("latest.bam");
        // Relative target ("run.bam") must resolve against the link's directory.
        std::os::unix::fs::symlink("run.bam", &link).unwrap();
        assert_eq!(resolve_symlink_output(&link).unwrap(), target);
    }

    /// A cyclic symlink chain bails at `MAX_LINKS` with an error instead of
    /// looping forever.
    #[cfg(unix)]
    #[test]
    fn test_resolve_symlink_output_detects_cycle() {
        let dir = tempfile::tempdir().unwrap();
        let a = dir.path().join("a");
        let b = dir.path().join("b");
        std::os::unix::fs::symlink(&b, &a).unwrap();
        std::os::unix::fs::symlink(&a, &b).unwrap();
        assert!(resolve_symlink_output(&a).is_err());
    }

    /// Overwriting an existing destination keeps its current mode (matching
    /// `File::create`, which never re-chmods an existing file).
    #[cfg(unix)]
    #[test]
    fn test_target_file_mode_keeps_existing_file_mode() {
        use std::os::unix::fs::PermissionsExt;
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bam");
        std::fs::write(&path, b"x").unwrap();
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o640)).unwrap();
        assert_eq!(target_file_mode(&path), Some(0o640));
    }

    /// A new destination gets `0o666 & !umask` — the mode `File::create` would
    /// have produced — not the temp file's private `0600`.
    #[cfg(unix)]
    #[test]
    fn test_target_file_mode_new_file_uses_umask() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("does-not-exist.bam");
        assert_eq!(target_file_mode(&path), Some(0o666 & !process_umask()));
    }

    // ========================================================================
    // LibraryLookup tests
    // ========================================================================

    #[test]
    fn test_library_lookup_empty_header() {
        let header = Header::builder().build();
        let lookup = LibraryLookup::from_header(&header);
        assert!(lookup.rg_to_ordinal.is_empty());
    }

    #[test]
    fn test_library_lookup_single_rg() {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let header = Header::builder().add_read_group(BString::from("rg1"), rg).build();

        let lookup = LibraryLookup::from_header(&header);
        assert_eq!(lookup.rg_to_ordinal.len(), 1);
        // LibA is the only library, so it gets ordinal 1 (0 is reserved for empty/unknown)
        assert_eq!(
            *lookup.rg_to_ordinal.get(b"rg1".as_slice()).expect("rg1 should be in ordinal map"),
            1
        );
    }

    #[test]
    fn test_library_lookup_multiple_libraries() {
        let rg_a = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibC"))
            .build()
            .expect("valid");
        let rg_b = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_c = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .expect("valid");

        let header = Header::builder()
            .add_read_group(BString::from("rg1"), rg_a)
            .add_read_group(BString::from("rg2"), rg_b)
            .add_read_group(BString::from("rg3"), rg_c)
            .build();

        let lookup = LibraryLookup::from_header(&header);
        assert_eq!(lookup.rg_to_ordinal.len(), 3);

        // Libraries sorted alphabetically: LibA=1, LibB=2, LibC=3
        let rg2 = *lookup.rg_to_ordinal.get(b"rg2".as_slice()).expect("rg2");
        let rg3 = *lookup.rg_to_ordinal.get(b"rg3".as_slice()).expect("rg3");
        let rg1 = *lookup.rg_to_ordinal.get(b"rg1".as_slice()).expect("rg1");
        assert_eq!(rg2, 1); // LibA
        assert_eq!(rg3, 2); // LibB
        assert_eq!(rg1, 3); // LibC
    }

    #[test]
    fn test_distinct_header_ordinals() {
        // (i) no read groups -> 0.
        let empty = LibraryLookup::from_header(&Header::builder().build());
        assert_eq!(empty.distinct_header_ordinals(), 0);

        // (ii) one RG with an LB -> 1.
        let rg_a = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let one = LibraryLookup::from_header(
            &Header::builder().add_read_group(BString::from("rg1"), rg_a).build(),
        );
        assert_eq!(one.distinct_header_ordinals(), 1);

        // (iii) two RGs with two different LBs -> 2.
        let rg_b = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_c = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .expect("valid");
        let two = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_b)
                .add_read_group(BString::from("rg2"), rg_c)
                .build(),
        );
        assert_eq!(two.distinct_header_ordinals(), 2);

        // (iv) two RGs with the SAME LB -> 1.
        let rg_d = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_e = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let same = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_d)
                .add_read_group(BString::from("rg2"), rg_e)
                .build(),
        );
        assert_eq!(same.distinct_header_ordinals(), 1);

        // (v) one RG with an LB + one RG WITHOUT an LB -> 2 (ordinals {1, 0}).
        let rg_f = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_no_lb = Map::<ReadGroup>::builder().build().expect("valid");
        let mixed = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_f)
                .add_read_group(BString::from("rg2"), rg_no_lb)
                .build(),
        );
        assert_eq!(mixed.distinct_header_ordinals(), 2);
    }

    #[test]
    fn test_library_lookup_unknown_rg_returns_zero() {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let header = Header::builder().add_read_group(BString::from("rg1"), rg).build();

        let lookup = LibraryLookup::from_header(&header);
        // A BAM record with no RG tag or unknown RG should return ordinal 0
        // We can test get_ordinal with a minimal BAM record that has no aux data
        let mut bam = vec![0u8; 36];
        bam[8] = 4; // l_read_name = 4
        bam[32..36].copy_from_slice(b"rea\0");
        assert_eq!(lookup.get_ordinal(&bam), 0);
    }

    // ========================================================================
    // RawExternalSorter builder tests
    // ========================================================================

    #[test]
    fn test_raw_sorter_defaults() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate);
        assert_eq!(sorter.memory_limit, 512 * 1024 * 1024);
        assert!(sorter.temp_dirs.is_empty());
        assert_eq!(sorter.threads, 1);
        assert_eq!(sorter.output_compression, 6);
        assert_eq!(sorter.temp_compression, 1);
        assert!(!sorter.write_index);
        assert!(sorter.pg_info.is_none());
        // Fixed, not derived: a bare sorter must behave the same on every host.
        // Deriving from `ulimit -n` is a decision the command line makes, not
        // one the engine imposes on every embedder.
        assert_eq!(sorter.max_temp_files, crate::fd_limit::FALLBACK_MAX_TEMP_FILES);
    }

    #[test]
    fn test_raw_sorter_builder_chain() {
        let sorter = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .memory_limit(1024)
            .temp_dir(PathBuf::from("/tmp/test"))
            .threads(8)
            .output_compression(9)
            .temp_compression(3)
            .write_index(true)
            .pg_info("1.0".to_string(), "fgumi sort".to_string())
            .max_temp_files(128);

        assert_eq!(sorter.memory_limit, 1024);
        assert_eq!(sorter.temp_dirs, vec![PathBuf::from("/tmp/test")]);
        assert_eq!(sorter.threads, 8);
        assert_eq!(sorter.output_compression, 9);
        assert_eq!(sorter.temp_compression, 3);
        assert!(sorter.write_index);
        assert_eq!(sorter.pg_info, Some(("1.0".to_string(), "fgumi sort".to_string())));
        assert_eq!(sorter.max_temp_files, 128);
    }

    #[test]
    fn test_raw_sorter_memory_limit() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate).memory_limit(256 * 1024 * 1024);
        assert_eq!(sorter.memory_limit, 256 * 1024 * 1024);
    }

    #[test]
    fn test_raw_sorter_temp_compression() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate).temp_compression(0);
        assert_eq!(sorter.temp_compression, 0);
    }

    /// `RawExternalSorter::sort` rejects `temp_compression=0` + `SpillCodec::Zstd`
    /// before doing any work, mirroring the CLI guard in `commands::sort`. zstd
    /// has no level-0 "stored" mode, and silently remapping to 1 would surprise
    /// API callers who pass 0 expecting uncompressed spills.
    #[test]
    fn test_raw_sorter_sort_rejects_temp_compression_zero_with_zstd() {
        use fgumi_sam::SamBuilder;

        let mut builder = SamBuilder::new();
        let _ = builder.add_pair().name("read0").start1(1).start2(101).build();

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        let err = RawExternalSorter::new(SortOrder::Coordinate)
            .spill_codec(crate::codec::SpillCodec::Zstd)
            .temp_compression(0)
            .sort(&input, &output)
            .expect_err("sort should reject temp_compression=0 + zstd");
        let msg = err.to_string();
        assert!(
            msg.contains("temp_compression=0 is only supported with SpillCodec::Bgzf"),
            "unexpected error: {msg}"
        );
        assert!(!output.exists(), "no output should be produced on validation failure");
    }

    #[test]
    fn test_raw_sorter_max_temp_files() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate).max_temp_files(0);
        assert_eq!(sorter.max_temp_files, 0);
    }

    /// An ordinary merge width must not warn — every `fgumi merge` opens a
    /// handful of inputs, so a warning there would be noise on every run —
    /// while a merge past the budget must carry the two things that make the
    /// warning actionable: how many inputs are being opened, and the
    /// `ulimit -n` lever.
    ///
    /// The budget is supplied rather than read from the host, so every branch
    /// is exercised on every target. Reading it made both assertions depend on
    /// the machine: the silent case needs a soft limit of at least
    /// `8 + FD_RESERVE` and failed under a `ulimit -n` below 40, and the
    /// warning case reached its branch only because `usize::MAX` saturates the
    /// comparison on a 64-bit target.
    #[rstest]
    #[case::ordinary_merge_fits(Some(1024), 8, None)]
    // Exactly at the budget: 64 - FD_RESERVE is 32, and 32 inputs still fit.
    #[case::exactly_at_the_budget(Some(64), 32, None)]
    #[case::one_past_the_budget(Some(64), 33, Some(32))]
    // A budget below the reserve saturates to zero, so any merge overruns it.
    #[case::budget_below_the_reserve(Some(8), 1, Some(0))]
    // Nothing to exceed where the limit cannot be read; silent by design.
    #[case::unreadable_budget_never_warns(None, usize::MAX, None)]
    fn test_fd_budget_warning_names_the_width_and_the_lever(
        #[case] soft: Option<u64>,
        #[case] num_inputs: usize,
        #[case] expected_budget: Option<u64>,
    ) {
        let warning = RawExternalSorter::fd_budget_warning_from_nofile(soft, num_inputs);
        match expected_budget {
            None => assert_eq!(warning, None, "a merge within the budget must not warn"),
            Some(budget) => {
                let warning = warning.expect("a merge past the budget must warn");
                assert!(
                    warning.contains(&num_inputs.to_string()),
                    "input count missing: {warning}"
                );
                assert!(
                    warning.contains(&format!("budget of about {budget}")),
                    "warning must name the budget it tripped: {warning}"
                );
                assert!(warning.contains("ulimit -n"), "remedy missing: {warning}");
            }
        }
    }

    /// The production wrapper must be exactly the helper applied to the host's
    /// own budget — the seam exists to make the arithmetic testable, not to let
    /// the two paths answer differently.
    ///
    /// Host-independent despite reading the real limit, because both sides read
    /// the same one: whatever this machine reports, the answers must agree.
    #[test]
    fn test_fd_budget_warning_reads_the_process_budget() {
        let soft = crate::fd_limit::soft_nofile();
        for num_inputs in [1_usize, 8, 1024, usize::MAX] {
            assert_eq!(
                RawExternalSorter::fd_budget_warning(num_inputs),
                RawExternalSorter::fd_budget_warning_from_nofile(soft, num_inputs),
                "wrapper and helper disagree at {num_inputs} inputs"
            );
        }
    }

    // ========================================================================
    // create_output_header tests
    // ========================================================================

    #[test]
    fn test_create_output_header_coordinate() {
        let sorter = RawExternalSorter::new(SortOrder::Coordinate);
        let header = Header::builder().build();
        let output_header = sorter.create_output_header(&header);

        let hd = output_header.header().expect("header should have HD record");
        let so = hd.other_fields().get(b"SO").expect("should have SO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"coordinate");
    }

    #[test]
    fn test_create_output_header_queryname() {
        let sorter = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()));
        let header = Header::builder().build();
        let output_header = sorter.create_output_header(&header);

        let hd = output_header.header().expect("header should have HD record");
        let so = hd.other_fields().get(b"SO").expect("should have SO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"queryname");
    }

    #[test]
    fn test_create_output_header_template_coordinate() {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate);
        let header = Header::builder().build();
        let output_header = sorter.create_output_header(&header);

        let hd = output_header.header().expect("header should have HD record");
        let fields = hd.other_fields();

        let so = fields.get(b"SO").expect("should have SO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(so), b"unsorted");

        let go = fields.get(b"GO").expect("should have GO tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(go), b"query");

        let ss = fields.get(b"SS").expect("should have SS tag");
        assert_eq!(<_ as AsRef<[u8]>>::as_ref(ss), b"unsorted:template-coordinate");
    }

    // ========================================================================
    // find_string_tag_in_record tests (via fgumi_raw_bam)
    // ========================================================================

    /// Helper to build minimal BAM bytes with aux data appended.
    /// Fixed 32-byte header + read name "rea\0" (`l_read_name=4`) + no cigar + no seq + aux.
    fn build_bam_with_aux(aux_data: &[u8]) -> Vec<u8> {
        let l_read_name: u8 = 4; // "rea" + null
        let mut bam = vec![0u8; 32];
        bam[8] = l_read_name;
        // n_cigar_op = 0 (bytes 12-13 already zero)
        // l_seq = 0 (bytes 16-19 already zero)
        // read name
        bam.extend_from_slice(b"rea\0");
        // no cigar, no seq, no qual
        // append aux data
        bam.extend_from_slice(aux_data);
        bam
    }

    #[rstest::rstest]
    #[case::present(b"RGZgroup1\0".as_slice(), Some(b"group1".as_slice()))]
    #[case::absent(b"".as_slice(), None)]
    fn test_find_rg_tag(#[case] aux_data: &[u8], #[case] expected: Option<&[u8]>) {
        let bam = build_bam_with_aux(aux_data);
        assert_eq!(
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(SamTag::RG),
            expected
        );
    }

    #[test]
    fn test_find_rg_tag_after_other_tags() {
        // XY:i:42 followed by RG:Z:mygroup — aux built dynamically for the integer tag
        let mut aux = Vec::new();
        aux.extend_from_slice(b"XYi");
        aux.extend_from_slice(&42i32.to_le_bytes());
        aux.extend_from_slice(b"RGZmygroup\0");
        let bam = build_bam_with_aux(&aux);
        assert_eq!(
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(SamTag::RG),
            Some(b"mygroup".as_slice())
        );
    }

    // ========================================================================
    // RawExternalSorter cell_tag builder tests
    // ========================================================================

    #[test]
    fn test_raw_sorter_cell_tag_default_is_none() {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate);
        assert!(sorter.cell_tag.is_none());
    }

    #[test]
    fn test_raw_sorter_cell_tag_builder() {
        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate).cell_tag(SamTag::CB);
        assert_eq!(sorter.cell_tag, Some(SamTag::CB));
    }

    /// Smoke test: an auto-detected lite (no-CB, no-tertiary) template-coordinate
    /// sort runs end-to-end and preserves the record count. With the default
    /// `KeyTypesSpec::Auto`, records lacking CB/MI/multi-library auxiliary data
    /// select the narrow `TemplateKey24` lane key; this exercises that path
    /// through Phase 1 and Phase 2 and confirms output completeness.
    #[test]
    fn template_sort_auto_lite_roundtrips_record_count() {
        let dir = tempfile::tempdir().unwrap();
        let (sorted, names) =
            create_sorted_bam(dir.path(), "lite", 50, 0, SortOrder::TemplateCoordinate);
        assert_eq!(collect_read_names(&sorted).len(), names.len() * 2);
    }

    /// Regression: a 0-record template-coordinate sort must still write a valid
    /// header-only output BAM (matching the coordinate and queryname orders),
    /// not skip the output entirely. The empty-input path skips the pre-loop
    /// first-record push, runs the read loop zero times, and falls through to
    /// the `chunk_files.is_empty()` fast path that writes the header-only BAM.
    #[test]
    fn template_sort_empty_input_writes_header_only_bam() {
        let dir = tempfile::tempdir().unwrap();
        let (sorted, names) =
            create_sorted_bam(dir.path(), "empty", 0, 0, SortOrder::TemplateCoordinate);
        assert!(names.is_empty());
        assert!(sorted.exists(), "empty-input sort must still produce an output file");
        assert_eq!(count_bam_records(&sorted), 0, "output must contain only the header");
    }

    // ========================================================================
    // Phase 2 teardown: Phase2Guard::finish_output
    // ========================================================================

    /// Everything an empty-source merge needs: a sorter, its pool (already
    /// handed over to Phase 2), the output header, and a temp dir holding one
    /// spill chunk that decodes to zero records plus the output path.
    ///
    /// `_dir` is returned only to keep the temp dir alive for the caller's
    /// lifetime; dropping it deletes both paths.
    struct EmptyMergeFixture {
        sorter: RawExternalSorter,
        pool: Arc<SortWorkerPool>,
        header: Header,
        chunk: PathBuf,
        output: PathBuf,
        _dir: tempfile::TempDir,
    }

    /// Builds an [`EmptyMergeFixture`].
    ///
    /// A real Phase 1 checks its spill trigger only after pushing a record, so
    /// it never writes an empty chunk and the "every source is empty" merge
    /// branch is unreachable from `sort()`. Handing the merge a chunk that holds
    /// nothing but a BGZF EOF block is the way to drive that branch directly,
    /// and it is the shape `consolidate_chunks` emits for a zero-record run.
    ///
    /// No compression setting is configured, because none affects these tests:
    /// the read side derives the chunk's codec from its magic bytes
    /// (`SortWorkerPool::set_phase2_files`) rather than from `spill_codec`, and
    /// the assertions are on the pipeline phase and the record count, not on
    /// output size. `threads` is set because it sizes the pool.
    fn empty_merge_fixture() -> EmptyMergeFixture {
        use crate::keys::RawCoordinateKey;

        let dir = tempfile::tempdir().expect("tempdir");
        let chunk = dir.path().join("empty.chunk");
        GenericKeyedChunkWriter::<RawCoordinateKey>::create(
            &chunk, /* compression_level */ 1, /* num_threads */ 1,
        )
        .expect("create empty chunk")
        .finish()
        .expect("finish empty chunk");

        let sorter = RawExternalSorter::new(SortOrder::Coordinate).threads(2);
        let pool = sorter.create_worker_pool().expect("worker pool");
        sorter.enter_output_phase(&pool);

        EmptyMergeFixture {
            sorter,
            pool,
            header: default_merge_header(),
            output: dir.path().join("output.bam"),
            chunk,
            _dir: dir,
        }
    }

    /// The merge sources are released *before* the output is finalized, and the
    /// pool leaves Phase 2 only *after*.
    ///
    /// [`Phase2Guard::finish_output`] exists for the first half: a `finish` that
    /// blocks draining the writer must not still be pinning the merge's file
    /// descriptors and 2 MiB-per-chunk reorder buffers. So the observation is
    /// made from *inside* the finalize closure, where the real writer's
    /// `finish()` runs — not merely before and after it.
    ///
    /// The phase is checked in the same place, but note it is no longer
    /// load-bearing for compression: blocks carry their own `CompressTarget`, so
    /// the tail is written at `output_compression` in any phase. What it pins is
    /// that the pool's phase still describes what the pool is doing while it is
    /// doing it.
    #[test]
    fn test_phase2_guard_releases_sources_before_finalizing_output() {
        use crate::keys::RawCoordinateKey;
        use crate::worker_pool::phase;

        let fixture = empty_merge_fixture();
        let pool = &fixture.pool;
        let (_sources, guard) = RawExternalSorter::setup_phase2_merge::<RawCoordinateKey>(
            std::slice::from_ref(&fixture.chunk),
            MemorySources::Shared(Vec::new()),
            pool,
        )
        .expect("setup phase 2");

        assert_eq!(pool.current_phase(), phase::PHASE2, "setup must arm Phase 2");
        assert_eq!(pool.phase2_files().len(), 1, "the disk source must be published");

        let observed = guard
            .finish_output(|| Ok((pool.current_phase(), pool.phase2_files().len())))
            .expect("finish_output must propagate the closure's Ok");

        assert_eq!(
            observed,
            (phase::PHASE2, 0),
            "the output must be finalized with the drained spill files already unpublished, and \
             before the pool leaves Phase 2"
        );
        // `finish_output` consumed the guard, so its `Drop` has already run — a
        // second deactivate must not have disturbed the phase it just set.
        assert_eq!(
            pool.current_phase(),
            phase::LEGACY,
            "finish_output must return the pool to LEGACY once the output is done"
        );
    }

    /// A failing finalizer still hands the pool back in `LEGACY` with the spill
    /// files unpublished, and the error still reaches the caller unwrapped.
    ///
    /// A writer that fails mid-`finish` is the one exit the merge cannot retry,
    /// so the pool it leaves behind is what the rest of the sort — and the
    /// pool's own shutdown — has to cope with. Because
    /// [`Phase2Guard::finish_output`] consumes the guard, the reset here is
    /// guaranteed twice over (explicitly, then by `Drop` as the method returns);
    /// this pins the resulting contract rather than either mechanism.
    #[test]
    fn test_phase2_guard_finish_output_resets_phase_on_error() {
        use crate::keys::RawCoordinateKey;
        use crate::worker_pool::phase;

        let fixture = empty_merge_fixture();
        let pool = &fixture.pool;
        let (_sources, guard) = RawExternalSorter::setup_phase2_merge::<RawCoordinateKey>(
            std::slice::from_ref(&fixture.chunk),
            MemorySources::Shared(Vec::new()),
            pool,
        )
        .expect("setup phase 2");

        let failed =
            guard.finish_output(|| -> Result<()> { Err(anyhow::anyhow!("writer blew up")) });

        assert!(failed.is_err(), "finish_output must propagate the finalizer's error");
        assert_eq!(
            failed.unwrap_err().to_string(),
            "writer blew up",
            "the finalizer's own error must reach the caller, not a wrapped one"
        );
        assert_eq!(
            pool.current_phase(),
            phase::LEGACY,
            "a failed finalize must still return the pool to LEGACY"
        );
        assert!(
            pool.phase2_files().is_empty(),
            "a failed finalize must still unpublish the spill files"
        );
    }

    /// A merge whose every source is empty writes a header-only BAM and still
    /// hands the pool back in `LEGACY`.
    ///
    /// This early-return branch finalizes its writer through the same
    /// [`Phase2Guard::finish_output`] as the main loop, so it has the same way
    /// to go wrong.
    #[test]
    fn test_merge_chunks_generic_with_only_empty_sources() {
        use crate::keys::RawCoordinateKey;
        use crate::worker_pool::phase;

        let fixture = empty_merge_fixture();
        let merged = fixture
            .sorter
            .merge_chunks_generic::<RawCoordinateKey>(
                std::slice::from_ref(&fixture.chunk),
                MemorySources::Shared(Vec::new()),
                &fixture.header,
                &fixture.output,
                0,
                &fixture.pool,
            )
            .expect("empty merge should succeed");

        assert_eq!(merged, 0, "an all-empty merge must report zero records merged");
        assert_eq!(
            count_bam_records(&fixture.output),
            0,
            "output must hold the header and nothing else"
        );
        assert_eq!(
            fixture.pool.current_phase(),
            phase::LEGACY,
            "the merge must hand the pool back in LEGACY"
        );
        assert!(fixture.pool.phase2_files().is_empty(), "the merge must unpublish the spill files");
    }

    /// The indexed counterpart of
    /// [`test_merge_chunks_generic_with_only_empty_sources`]: the BAI still has
    /// to be built and cover every reference, even with no records to index.
    #[test]
    fn test_merge_chunks_with_index_with_only_empty_sources() {
        use crate::keys::RawCoordinateKey;
        use crate::worker_pool::phase;

        let fixture = empty_merge_fixture();
        let (index, records_merged) = fixture
            .sorter
            .merge_chunks_with_index::<RawCoordinateKey>(
                std::slice::from_ref(&fixture.chunk),
                MemorySources::Shared(Vec::new()),
                &fixture.header,
                &fixture.output,
                0,
                &fixture.pool,
            )
            .expect("empty indexed merge should succeed");

        assert_eq!(records_merged, 0, "an empty merge must report zero records written");
        assert_eq!(
            count_bam_records(&fixture.output),
            0,
            "output must hold the header and nothing else"
        );
        assert_eq!(
            index.reference_sequences().len(),
            fixture.header.reference_sequences().len(),
            "the index must cover every reference sequence even with no records"
        );
        assert_eq!(
            fixture.pool.current_phase(),
            phase::LEGACY,
            "the indexed merge must hand the pool back in LEGACY"
        );
        assert!(fixture.pool.phase2_files().is_empty(), "the merge must unpublish the spill files");
    }

    // ========================================================================
    // extract_template_key_inline cell_tag tests
    // ========================================================================

    fn test_cb_hasher() -> ahash::RandomState {
        cb_hasher()
    }

    /// Build minimal BAM bytes with mapped read at (tid, pos) on the forward strand,
    /// with optional aux data appended.
    #[allow(clippy::cast_possible_truncation)]
    fn build_mapped_bam(tid: i32, pos: i32, name: &[u8], aux: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let mut bam = vec![0u8; 32];
        // ref_id at offset 0..4
        bam[0..4].copy_from_slice(&tid.to_le_bytes());
        // pos at offset 4..8
        bam[4..8].copy_from_slice(&pos.to_le_bytes());
        // l_read_name at offset 8
        bam[8] = l_read_name;
        // flags at offset 14..16: paired + proper pair = 0x03
        bam[14..16].copy_from_slice(&3u16.to_le_bytes());
        // mate_ref_id at offset 20..24: same tid
        bam[20..24].copy_from_slice(&tid.to_le_bytes());
        // mate_pos at offset 24..28: same pos
        bam[24..28].copy_from_slice(&pos.to_le_bytes());
        // read name
        bam.extend_from_slice(name);
        bam.push(0); // null terminator
        // no cigar, no seq, no qual
        // aux data
        bam.extend_from_slice(aux);
        bam
    }

    /// Build CB:Z: aux tag bytes.
    fn cb_aux(value: &[u8]) -> Vec<u8> {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"CBZ");
        aux.extend_from_slice(value);
        aux.push(0); // null terminator
        aux
    }

    #[test]
    fn test_extract_template_key_cb_present_has_nonzero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "CB present should produce non-zero cb_hash");
    }

    #[test]
    fn test_extract_template_key_cb_absent_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        // No CB tag in aux data
        let bam = build_mapped_bam(0, 100, b"read1", &[]);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_eq!(key.cb_hash, 0, "missing CB tag should produce cb_hash=0");
    }

    #[test]
    fn test_extract_template_key_cell_tag_none_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key = extract_template_key_inline(&bam, &lib_lookup, None, &test_cb_hasher());
        assert_eq!(key.cb_hash, 0, "cell_tag=None should produce cb_hash=0");
    }

    #[test]
    fn test_extract_template_key_different_cb_values_differ() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        let aux1 = cb_aux(b"ACGTACGT");
        let bam1 = build_mapped_bam(0, 100, b"read1", &aux1);
        let key1 =
            extract_template_key_inline(&bam1, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());

        let aux2 = cb_aux(b"TGCATGCA");
        let bam2 = build_mapped_bam(0, 100, b"read1", &aux2);
        let key2 =
            extract_template_key_inline(&bam2, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());

        assert_ne!(
            key1.cb_hash, key2.cb_hash,
            "different CB values should produce different hashes"
        );
    }

    #[test]
    fn test_extract_template_key_cb_hash_is_deterministic() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key1 =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        let key2 =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_eq!(key1.cb_hash, key2.cb_hash, "same input should produce same cb_hash");
    }

    /// Encode a `tc:B:i` aux tag carrying the 6-element template-coordinate
    /// array `[tid1, pos1, neg1, tid2, pos2, neg2]` that `fgumi zipper` stamps
    /// onto secondary/supplementary reads.
    fn tc_aux(vals: &[i32; 6]) -> Vec<u8> {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"tcB"); // tag "tc", array type 'B'
        aux.push(b'i'); // subtype: int32
        #[allow(clippy::cast_possible_truncation)]
        aux.extend_from_slice(&(vals.len() as u32).to_le_bytes());
        for v in vals {
            aux.extend_from_slice(&v.to_le_bytes());
        }
        aux
    }

    /// Build minimal BAM bytes for a mapped secondary/supplementary alignment:
    /// paired, mate mapped, at `own_pos`, with its mate at `mate_pos`, carrying the
    /// caller-supplied `sec_supp_flag` (`flags::SECONDARY` or `flags::SUPPLEMENTARY`),
    /// plus optional aux. Both flags route through the same `tc` keying branch.
    #[allow(clippy::cast_possible_truncation)]
    fn build_secondary_or_supplementary_bam(
        own_pos: i32,
        mate_pos: i32,
        sec_supp_flag: u16,
        name: &[u8],
        aux: &[u8],
    ) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8;
        let mut bam = vec![0u8; 32];
        bam[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id = chr0
        bam[4..8].copy_from_slice(&own_pos.to_le_bytes());
        bam[8] = l_read_name;
        // flags: PAIRED | (SECONDARY or SUPPLEMENTARY) -- mate is mapped.
        bam[14..16].copy_from_slice(&(flags::PAIRED | sec_supp_flag).to_le_bytes());
        bam[20..24].copy_from_slice(&0i32.to_le_bytes()); // mate ref_id = chr0
        bam[24..28].copy_from_slice(&mate_pos.to_le_bytes());
        bam.extend_from_slice(name);
        bam.push(0);
        bam.extend_from_slice(aux);
        bam
    }

    /// Build minimal BAM bytes for a mapped primary read on `chr0` with the
    /// caller-supplied `flag`, aligned at `own_pos` with its mate at `mate_pos`.
    /// No CIGAR, so the unclipped 5' position equals `own_pos`/`mate_pos`. Used to
    /// key the two primaries of a pair through the per-record (mapped) path so
    /// their `TemplateKey`s can be compared against the tc-keyed sec/supp keys.
    #[allow(clippy::cast_possible_truncation)]
    fn build_primary_bam(own_pos: i32, mate_pos: i32, flag: u16, name: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8;
        let mut bam = vec![0u8; 32];
        bam[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id = chr0
        bam[4..8].copy_from_slice(&own_pos.to_le_bytes());
        bam[8] = l_read_name;
        bam[14..16].copy_from_slice(&flag.to_le_bytes());
        bam[20..24].copy_from_slice(&0i32.to_le_bytes()); // mate ref_id = chr0
        bam[24..28].copy_from_slice(&mate_pos.to_le_bytes());
        bam.extend_from_slice(name);
        bam.push(0);
        bam
    }

    /// A secondary/supplementary alignment carries its own position and its mate's
    /// position, but **not** its own primary's position. Per-record keying (samtools
    /// / fgbio) can therefore only place it by the mate coordinate, which can drift a
    /// fragment-length from the template's true coordinate. Because `fgumi zipper`
    /// pre-computes the exact template coordinate into the `tc` tag, fgumi's
    /// template-coordinate sort can place the read **exactly** at its primary pair's
    /// coordinate — a placement samtools/fgbio cannot achieve. This is the novel
    /// fgumi capability; both SECONDARY and SUPPLEMENTARY take the tc branch.
    #[rstest]
    #[case::secondary(flags::SECONDARY)]
    #[case::supplementary(flags::SUPPLEMENTARY)]
    fn test_extract_template_key_secondary_supplementary_uses_tc_tag(#[case] sec_supp_flag: u16) {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        // Sec/supp of a left-anchored R1: the primary is at pos 1000, its mate (R2
        // primary) at 1400, but this alignment is at 9000. The exact template
        // coordinate (pos1) is 1000, carried in `tc`.
        let tc = tc_aux(&[0, 1000, 0, 0, 1400, 1]);
        let supp = build_secondary_or_supplementary_bam(9000, 1400, sec_supp_flag, b"tmpl", &tc);
        let supp_key = extract_template_key_inline(&supp, &lib_lookup, None, &test_cb_hasher());

        // A different template whose true coordinate (pos1) is 1200 — between the
        // read's exact coordinate (1000, from tc) and the mate-only approximation
        // (1400).
        let filler = build_mapped_bam(0, 1200, b"fillr", &[]);
        let filler_key = extract_template_key_inline(&filler, &lib_lookup, None, &test_cb_hasher());

        assert!(
            supp_key < filler_key,
            "sec/supp must sort at its template coordinate (pos1=1000 via tc), before \
             the pos1=1200 template; per-record keying misplaces it at the mate \
             coordinate (1400), after the 1200 template"
        );
    }

    /// Without a `tc` tag, a secondary/supplementary must fall back to the per-record
    /// (own-position / mate) computation — the prior behavior — so the tc branch
    /// never changes keying for reads zipper did not stamp.
    #[rstest]
    #[case::secondary(flags::SECONDARY)]
    #[case::supplementary(flags::SUPPLEMENTARY)]
    fn test_extract_template_key_secondary_supplementary_without_tc_falls_back(
        #[case] sec_supp_flag: u16,
    ) {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        // Same geometry as the tc test, but no `tc` tag.
        let supp = build_secondary_or_supplementary_bam(9000, 1400, sec_supp_flag, b"tmpl", &[]);
        let supp_key = extract_template_key_inline(&supp, &lib_lookup, None, &test_cb_hasher());

        // Fallback keys off the mate coordinate (1400), so it sorts AFTER the
        // pos1=1200 template — the unchanged, pre-tc behavior.
        let filler = build_mapped_bam(0, 1200, b"fillr", &[]);
        let filler_key = extract_template_key_inline(&filler, &lib_lookup, None, &test_cb_hasher());

        assert!(
            supp_key > filler_key,
            "without tc, a sec/supp keeps the per-record (mate=1400) placement"
        );
    }

    /// Tie-case: with `(tid1, pos1, tid2, pos2)` held identical between a primary
    /// pair and its tc-stamped sec/supp reads, the tc branch must place the
    /// sec/supp read at the *exact same* `TemplateKey` as its primary — not merely
    /// nearby. This is the strong property the reviewer asked for: it pins down the
    /// full key (position lanes, strand lanes, CB, library, MI, name, **and** the
    /// `is_upper`/`is_read2` tiebreak), verifying tc ordering matches the mapped
    /// path rather than only asserting a coarse `<`/`>` placement.
    ///
    /// Geometry is a standard FR pair on `chr0`: R1 forward at unclipped-5' 1000
    /// (the lower/lane-1 end), R2 reverse at 1400 (the upper/lane-2 end). The two
    /// primaries key through the per-record path; the two sec/supp reads key
    /// through the tc branch. `is_read2` mirrors the primary `is_upper` here (R1 ->
    /// false, R2 -> true), so each sec/supp key is byte-identical to its primary's.
    #[rstest]
    #[case::secondary(flags::SECONDARY)]
    #[case::supplementary(flags::SUPPLEMENTARY)]
    fn test_extract_template_key_secondary_supplementary_matches_primary_via_tc(
        #[case] sec_supp_flag: u16,
    ) {
        const FIRST_SEGMENT: u16 = 0x40;
        const SECOND_SEGMENT: u16 = 0x80;
        const REVERSE: u16 = 0x10;
        const MATE_REVERSE: u16 = 0x20;

        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let key =
            |bam: &[u8]| extract_template_key_inline(bam, &lib_lookup, None, &test_cb_hasher());

        // Primaries through the per-record path. R1 (forward, lower) is lane 1;
        // R2 (reverse, upper) is lane 2. Both canonicalize to the same coordinate.
        let r1_primary =
            build_primary_bam(1000, 1400, flags::PAIRED | FIRST_SEGMENT | MATE_REVERSE, b"tmpl");
        let r2_primary =
            build_primary_bam(1400, 1000, flags::PAIRED | SECOND_SEGMENT | REVERSE, b"tmpl");
        let r1_primary_key = key(&r1_primary);
        let r2_primary_key = key(&r2_primary);

        // The two primaries co-locate but differ only in the is_upper tiebreak.
        assert_eq!(
            r1_primary_key.core_cmp(&r2_primary_key),
            std::cmp::Ordering::Equal,
            "the pair's two primaries share (tid1, pos1, tid2, pos2)"
        );
        assert_ne!(
            r1_primary_key, r2_primary_key,
            "primaries differ only in the is_upper tiebreak (R1 lower, R2 upper)"
        );

        // Sec/supp reads through the tc branch, stamped with the pair's exact
        // canonical template coordinate, but aligned far away (own pos 9000).
        let tc = tc_aux(&[0, 1000, 0, 0, 1400, 1]);
        let r1_ss = build_secondary_or_supplementary_bam(
            9000,
            1400,
            sec_supp_flag | FIRST_SEGMENT,
            b"tmpl",
            &tc,
        );
        let r2_ss = build_secondary_or_supplementary_bam(
            9000,
            1000,
            sec_supp_flag | SECOND_SEGMENT,
            b"tmpl",
            &tc,
        );

        assert_eq!(
            key(&r1_ss),
            r1_primary_key,
            "R1 sec/supp keyed via tc must be byte-identical to the R1 primary key"
        );
        assert_eq!(
            key(&r2_ss),
            r2_primary_key,
            "R2 sec/supp keyed via tc must be byte-identical to the R2 primary key"
        );
    }

    #[test]
    fn test_extract_template_key_unmapped_with_cb() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");

        // Build unmapped read (flag = PAIRED | UNMAPPED | MATE_UNMAPPED = 0x0D)
        let mut bam = vec![0u8; 32];
        bam[8] = 6; // l_read_name = 6 ("read1" + null)
        bam[14..16].copy_from_slice(&0x000Du16.to_le_bytes()); // flags
        // ref_id = -1 (unmapped)
        bam[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        bam[4..8].copy_from_slice(&(-1i32).to_le_bytes()); // pos = -1
        bam[20..24].copy_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id = -1
        bam[24..28].copy_from_slice(&(-1i32).to_le_bytes()); // mate_pos = -1
        bam.extend_from_slice(b"read1\0");
        bam.extend_from_slice(&aux);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "unmapped read with CB should have non-zero cb_hash");
        assert_eq!(key.primary, u64::MAX, "unmapped both-mates should have MAX primary");
    }

    // ========================================================================
    // Consolidation chunk naming tests
    // ========================================================================

    /// Count records in a BAM file by reading with the raw BAM reader.
    fn count_bam_records(path: &std::path::Path) -> u64 {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader).count() as u64
    }

    /// Assert the output really is in `order`.
    ///
    /// Built on the crate's own [`crate::verify::verify_sort_order`] and the same
    /// key extractors the merge uses, so it covers every sort order rather than
    /// only coordinate. Record counts alone cannot distinguish "merged correctly"
    /// from "emitted every record in the wrong order", which is what a bad append
    /// would produce.
    fn assert_sorted_in(order: SortOrder, path: &std::path::Path) {
        use crate::keys::{RawQuerynameKey, RawQuerynameLexKey, SortContext};
        use crate::reader::RawBamRecordReader;
        use crate::verify::verify_sort_order;

        let (_, header) = create_raw_bam_reader(path, 1).expect("open bam for header");
        let new_reader = || {
            let file = std::fs::File::open(path).expect("open bam");
            let mut reader = RawBamRecordReader::new(file).expect("bam reader");
            reader.skip_header().expect("skip header");
            reader
        };

        let summary = match order {
            SortOrder::Coordinate => {
                let n_ref = u32::try_from(header.reference_sequences().len())
                    .expect("reference sequence count fits in u32");
                verify_sort_order(
                    new_reader(),
                    |bam| crate::inline::extract_coordinate_key_inline(bam, n_ref),
                    |key, prev| key < prev,
                )
            }
            SortOrder::Queryname(QuerynameComparator::Natural) => {
                let ctx = SortContext::from_header(&header);
                verify_sort_order(
                    new_reader(),
                    |bam| RawQuerynameKey::extract(bam, &ctx),
                    |key, prev| key < prev,
                )
            }
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                let ctx = SortContext::from_header(&header);
                verify_sort_order(
                    new_reader(),
                    |bam| RawQuerynameLexKey::extract(bam, &ctx),
                    |key, prev| key < prev,
                )
            }
            SortOrder::TemplateCoordinate => {
                let lib_lookup = LibraryLookup::from_header(&header);
                let hasher = cb_hasher();
                verify_sort_order(
                    new_reader(),
                    |bam| extract_template_key_inline(bam, &lib_lookup, None, &hasher),
                    |key, prev| key < prev,
                )
            }
        };
        let (_, violations, _) = summary.expect("verify should read the output");
        assert_eq!(violations, 0, "{order:?}: output is not sorted");
    }

    /// Verifies that sort with consolidation preserves all records.
    ///
    /// Uses a tiny memory limit and low `max_temp_files` to force many chunks and
    /// consolidation. Before the fix (chunk files named by `chunk_files.len()`),
    /// consolidation would drain entries from the vector, shrinking its length,
    /// causing new chunks to collide with existing non-consolidated chunk files.
    ///
    /// Parameterized over both spill codecs so the consolidation merge writer
    /// path is exercised for each format. zstd uses `temp_compression(1)` since
    /// zstd has no level-0 "stored" mode (the CLI rejects 0+zstd; here we make
    /// the same choice explicitly so the test matches user-facing behavior).
    #[rstest::rstest]
    #[case::coordinate_bgzf(SortOrder::Coordinate, false, crate::codec::SpillCodec::Bgzf, 0)]
    #[case::coordinate_with_index_bgzf(
        SortOrder::Coordinate,
        true,
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::queryname_bgzf(
        SortOrder::Queryname(QuerynameComparator::default()),
        false,
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::queryname_natural_bgzf(
        SortOrder::Queryname(QuerynameComparator::Natural),
        false,
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::template_coordinate_bgzf(
        SortOrder::TemplateCoordinate,
        false,
        crate::codec::SpillCodec::Bgzf,
        0
    )]
    #[case::coordinate_zstd(SortOrder::Coordinate, false, crate::codec::SpillCodec::Zstd, 1)]
    #[case::coordinate_with_index_zstd(
        SortOrder::Coordinate,
        true,
        crate::codec::SpillCodec::Zstd,
        1
    )]
    #[case::queryname_zstd(
        SortOrder::Queryname(QuerynameComparator::default()),
        false,
        crate::codec::SpillCodec::Zstd,
        1
    )]
    #[case::queryname_natural_zstd(
        SortOrder::Queryname(QuerynameComparator::Natural),
        false,
        crate::codec::SpillCodec::Zstd,
        1
    )]
    #[case::template_coordinate_zstd(
        SortOrder::TemplateCoordinate,
        false,
        crate::codec::SpillCodec::Zstd,
        1
    )]
    fn test_sort_with_consolidation_preserves_all_records(
        #[case] sort_order: SortOrder,
        #[case] write_index: bool,
        #[case] spill_codec: crate::codec::SpillCodec,
        #[case] temp_compression: u32,
    ) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 30;
        let mut builder = SamBuilder::new();
        // Descending coordinates: ascending input would fold into a single
        // run under natural run formation, and this test needs several
        // chunks to exercise consolidation.
        for i in 0..num_pairs {
            let descending = num_pairs - 1 - i;
            let _ = builder
                .add_pair()
                .name(&format!("read{descending:05}"))
                .start1(descending * 200 + 1)
                .start2(descending * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Tiny memory limit forces many chunks; low max_temp_files triggers consolidation
        let stats = RawExternalSorter::new(sort_order)
            .memory_limit(1024) // 1 KB — each chunk holds very few records
            .max_temp_files(4)
            .spill_codec(spill_codec)
            .temp_compression(temp_compression)
            .output_compression(0)
            .write_index(write_index)
            .sort(&input, &output)
            .expect("sort should succeed");

        assert!(
            stats.runs_written >= 5,
            "expected at least 5 chunks to exercise post-consolidation naming, got {}",
            stats.runs_written
        );

        // Count records in the output BAM to verify no data was lost
        let expected = (num_pairs * 2) as u64;
        let observed = count_bam_records(&output);
        assert_eq!(observed, expected, "chunk filename collision likely lost data");
    }

    /// Collect (name, pos) for every record in a BAM, in file order.
    fn collect_names_and_positions(path: &Path) -> Vec<(String, i32)> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader)
            .map(|rec| {
                let v = fgumi_raw_bam::RawRecordView::new(rec.as_ref());
                let name = String::from_utf8(v.read_name().to_vec()).expect("valid UTF-8 name");
                (name, v.pos())
            })
            .collect()
    }

    /// Audit C2: standalone `fgumi sort --order queryname` must preserve input
    /// order for records whose queryname keys are exactly equal (same name +
    /// same flag class), matching `samtools sort -n` (stable `ks_mergesort`),
    /// fgbio, and fgumi's own arena/runall path. The in-chunk sorts used
    /// `sort_unstable`, so exact ties could emerge reordered.
    ///
    /// Construction: several distinct names, each appearing many times at
    /// distinct, strictly increasing positions, added interleaved. A stable sort
    /// must group by name and keep every name's records in ingest (position)
    /// order; an unstable sort may reorder the equal-key runs. Exercised in both
    /// the single-threaded (`entries.sort_unstable_by`) and multi-threaded
    /// (`par_chunks_mut` + per-chunk `sort_unstable_by`) in-memory paths, and
    /// with a small memory limit to also cover the per-chunk spill sort.
    #[rstest::rstest]
    #[case::single_threaded_in_memory(1, 64 * 1024 * 1024)]
    #[case::multi_threaded_in_memory(4, 64 * 1024 * 1024)]
    #[case::single_threaded_spill(1, 4096)]
    #[case::multi_threaded_spill(4, 4096)]
    fn test_sort_queryname_preserves_exact_tie_input_order(
        #[case] threads: usize,
        #[case] memory_limit: usize,
    ) {
        use fgumi_sam::SamBuilder;

        const NAMES: usize = 8;
        const COPIES: usize = 16;

        let mut builder = SamBuilder::new();
        let mut pos = 1usize;
        for _copy in 0..COPIES {
            for n in 0..NAMES {
                let _ = builder.add_frag().name(&format!("read{n:02}")).start(pos).build();
                pos += 1;
            }
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .threads(threads)
            .memory_limit(memory_limit)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        // For each name, the positions must come out strictly increasing — i.e.
        // exactly the ingest order, since positions were assigned in ingest order.
        let out = collect_names_and_positions(&output);
        let mut last_pos_by_name: std::collections::HashMap<String, i32> =
            std::collections::HashMap::new();
        for (name, pos) in out {
            if let Some(&prev) = last_pos_by_name.get(&name) {
                assert!(
                    pos > prev,
                    "queryname sort reordered exact ties for {name}: saw pos {pos} after {prev} \
                     (stable sort must preserve ingest order for equal keys)"
                );
            }
            last_pos_by_name.insert(name, pos);
        }
    }

    /// A `memory_limit` the whole fixture fits under, so the sort stays in memory.
    const IN_MEMORY_MEMORY_LIMIT: usize = 64 * 1024 * 1024;

    /// A `memory_limit` the fixture exceeds several times over, forcing spills.
    const SPILLING_MEMORY_LIMIT: usize = 4096;

    /// Sorter-level output identity for both queryname comparators.
    ///
    /// `test_sort_queryname_preserves_exact_tie_input_order` pins the tie rule
    /// for the default comparator; this pins the *entire* emitted record
    /// sequence, for the natural comparator as well, against the sequence a
    /// stable sort of the input produces. Names are chosen so the two
    /// comparators disagree ("read2" < "read10" naturally, "read10" < "read2"
    /// lexicographically), and each name repeats at strictly increasing
    /// positions so exact-key ties are present in every run.
    ///
    /// The `threads` x `memory_limit` product covers the single-threaded and
    /// `par_chunks_mut` in-memory paths plus the spill path, where keys are
    /// serialized without the ingest position and rebuilt at merge time with
    /// position 0 — so identity there rests on the merge breaking cross-chunk
    /// ties by source (ingest) order. The body asserts which path each limit
    /// took, so the spill case cannot quietly become a third in-memory run if
    /// the per-record memory estimate or the fixture size changes.
    #[rstest::rstest]
    #[case::lexicographic(
        QuerynameComparator::Lexicographic,
        ["read1", "read10", "read2", "read20", "read3"]
    )]
    #[case::natural(
        QuerynameComparator::Natural,
        ["read1", "read2", "read3", "read10", "read20"]
    )]
    fn test_sort_queryname_output_identity_matches_stable_baseline(
        #[case] comparator: QuerynameComparator,
        #[case] expected_name_order: [&str; 5],
        #[values(1, 4)] threads: usize,
        #[values(IN_MEMORY_MEMORY_LIMIT, SPILLING_MEMORY_LIMIT)] memory_limit: usize,
    ) {
        use fgumi_sam::SamBuilder;

        const COPIES: usize = 16;
        const NAMES: [&str; 5] = ["read1", "read2", "read3", "read10", "read20"];

        // Ingest order: every name once per copy, each at the next position, so a
        // record's position is also its ingest rank. `start` is 1-based while the
        // BAM record — and so `collect_names_and_positions` — stores it 0-based.
        let mut builder = SamBuilder::new();
        let mut ingested: Vec<(String, i32)> = Vec::new();
        let mut pos = 0usize;
        for _copy in 0..COPIES {
            for name in NAMES {
                let _ = builder.add_frag().name(name).start(pos + 1).build();
                ingested.push((name.to_string(), i32::try_from(pos).expect("position fits i32")));
                pos += 1;
            }
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        let stats = RawExternalSorter::new(SortOrder::Queryname(comparator))
            .threads(threads)
            .memory_limit(memory_limit)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        // Pin which path each limit took. Identity through the spill path is the
        // interesting half of this test — that is where the ingest position is
        // dropped on serialization and the merge has to break ties by source
        // order instead — so a spilling limit that quietly stopped spilling
        // would leave that half uncovered with the test still green.
        assert_eq!(
            stats.runs_written > 0,
            memory_limit == SPILLING_MEMORY_LIMIT,
            "a {memory_limit}-byte limit wrote {} spill run(s), which is not the path this \
             case is meant to exercise",
            stats.runs_written
        );

        // The stable baseline: ingest order, stably sorted by the comparator's
        // name order. `sort_by_key` is stable, so records sharing a name keep
        // their ingest order — exactly what the sorter must reproduce.
        let mut expected = ingested;
        expected.sort_by_key(|(name, _)| {
            expected_name_order
                .iter()
                .position(|expected_name| expected_name == name)
                .expect("every ingested name appears in the expected order")
        });

        assert_eq!(
            collect_names_and_positions(&output),
            expected,
            "{comparator} sort with {threads} thread(s) and a {memory_limit}-byte limit did not \
             reproduce the stable baseline sequence"
        );
    }

    /// Sort `input` to `output` in `order`, spilling or not according to
    /// `memory_limit`.
    fn run_sort(
        order: SortOrder,
        input: &std::path::Path,
        output: &std::path::Path,
        memory_limit: usize,
    ) -> RawSortStats {
        RawExternalSorter::new(order)
            .memory_limit(memory_limit)
            .max_temp_files(0)
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(input, output)
            .expect("sort should succeed")
    }

    /// Memory limit large enough that nothing spills.
    const NO_SPILL_MEMORY: usize = 256 * 1024 * 1024;
    /// Memory limit small enough to force many chunks.
    const SPILL_MEMORY: usize = 8 * 1024;

    /// Produce a BAM genuinely in `order`, by sorting without spilling.
    ///
    /// Built by sorting rather than by construction: no order implies another.
    /// A coordinate-sorted file is not in template-coordinate order (template
    /// keys derive from the earlier mate), and neither is in queryname order.
    fn presorted_in(dir: &std::path::Path, order: SortOrder, num_pairs: usize) -> PathBuf {
        let raw = build_ordered_bam(dir, num_pairs);
        let sorted = dir.join(format!("presorted-{}.bam", order.header_so_tag()));
        let stats = run_sort(order, &raw, &sorted, NO_SPILL_MEMORY);
        assert_eq!(stats.runs_written, 0, "setup must not spill");
        assert_eq!(stats.total_records, stats.output_records);
        sorted
    }

    /// Every sort order: input already in the requested order must collapse to a
    /// single run, because run formation is what makes that case cheap.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_presorted_input_collapses_to_one_run(#[case] order: SortOrder) {
        let dir = tempfile::tempdir().expect("tempdir");
        let presorted = presorted_in(dir.path(), order, 300);
        let output = dir.path().join("output.bam");

        let stats = run_sort(order, &presorted, &output, SPILL_MEMORY);

        assert_eq!(
            stats.runs_written, 1,
            "{order:?}: input already in this order should spill exactly one run, got {}",
            stats.runs_written
        );
        assert_eq!(stats.total_records, stats.output_records, "{order:?}: lost records");
        assert_eq!(count_bam_records(&output), 600, "{order:?}: wrong record count");
        assert_sorted_in(order, &output);
    }

    /// Every sort order: input NOT in the requested order must not collapse.
    /// This is the guard against run formation appending where it has no right to.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_unsorted_input_never_collapses_to_one_run(#[case] order: SortOrder) {
        let dir = tempfile::tempdir().expect("tempdir");
        // Descending coordinates with names that ascend independently, so the
        // input is out of order for coordinate and template-coordinate, and its
        // name order does not track its coordinate order either.
        let input = build_shuffled_bam(dir.path(), 300);
        let output = dir.path().join("output.bam");

        let stats = run_sort(order, &input, &output, SPILL_MEMORY);

        assert!(
            stats.runs_written >= 2,
            "{order:?}: unsorted input must not collapse into one run, got {}",
            stats.runs_written
        );
        assert_eq!(stats.total_records, stats.output_records, "{order:?}: lost records");
        assert_eq!(count_bam_records(&output), 600, "{order:?}: wrong record count");
        assert_sorted_in(order, &output);
    }

    /// Every sort order: run formation must not change the output bytes.
    ///
    /// The appending path and the no-spill path must agree byte for byte. This is
    /// the strongest available check -- it subsumes ordering, and for template it
    /// is what confirms the buffer's lane-wise radix sort and the merge's
    /// `K::cmp` agree, which run formation depends on.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_run_formation_output_is_byte_identical(#[case] order: SortOrder) {
        let dir = tempfile::tempdir().expect("tempdir");
        let presorted = presorted_in(dir.path(), order, 300);

        let spilled = dir.path().join("spilled.bam");
        let in_memory = dir.path().join("in_memory.bam");

        let spilled_stats = run_sort(order, &presorted, &spilled, SPILL_MEMORY);
        assert_eq!(spilled_stats.runs_written, 1, "{order:?}: expected the appending path");

        let in_memory_stats = run_sort(order, &presorted, &in_memory, NO_SPILL_MEMORY);
        assert_eq!(in_memory_stats.runs_written, 0, "{order:?}: expected the no-spill path");

        assert_eq!(
            std::fs::read(&spilled).expect("read spilled"),
            std::fs::read(&in_memory).expect("read in-memory"),
            "{order:?}: appending runs changed the output bytes"
        );
    }

    /// A BAM whose coordinates descend while its names ascend, so it is out of
    /// order under every sort order this engine supports.
    fn build_shuffled_bam(dir: &std::path::Path, num_pairs: usize) -> PathBuf {
        use fgumi_sam::SamBuilder;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let descending = num_pairs - 1 - i;
            let _ = builder
                .add_pair()
                // Names descend too, so queryname order is also violated.
                .name(&format!("read{descending:05}"))
                .start1(descending * 200 + 1)
                .start2(descending * 200 + 101)
                .build();
        }
        let input = dir.join("shuffled.bam");
        builder.write_bam(&input).expect("write bam");
        input
    }

    /// Build a coordinate-ascending BAM: names and positions both increase.
    ///
    /// The starting point for every run-formation test. Tests that need input out
    /// of order sort this into the order under test first (`presorted_in`) or use
    /// [`build_shuffled_bam`].
    fn build_ordered_bam(dir: &std::path::Path, num_pairs: usize) -> PathBuf {
        use fgumi_sam::SamBuilder;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let n = i;
            let _ = builder
                .add_pair()
                .name(&format!("read{i:05}"))
                .start1(n * 200 + 1)
                .start2(n * 200 + 101)
                .build();
        }
        let input = dir.join("input.bam");
        builder.write_bam(&input).expect("write bam");
        input
    }

    /// Build a coordinate BAM that ascends, drops once, then ascends again.
    ///
    /// The descent is placed past the first spill boundary so it closes an open
    /// run rather than landing inside the first chunk, where sorting would
    /// absorb it.
    fn build_descent_bam(dir: &std::path::Path, num_pairs: usize) -> PathBuf {
        use fgumi_sam::SamBuilder;
        let mut builder = SamBuilder::new();
        let half = num_pairs / 2;
        for i in 0..num_pairs {
            // Second half restarts low, so exactly one chunk boundary sees a key
            // below the open run's last key.
            let n = if i < half { i } else { i - half };
            let _ = builder
                .add_pair()
                .name(&format!("read{i:05}"))
                .start1(n * 200 + 1)
                .start2(n * 200 + 101)
                .build();
        }
        let input = dir.join("descent.bam");
        builder.write_bam(&input).expect("write bam");
        input
    }

    /// The case that distinguishes run formation from both extremes: a descent
    /// closes the open run and opens exactly one more. Neither the all-sorted nor
    /// the reverse-sorted test reaches [`RunFormer::place`]'s fresh-path branch
    /// while a run is open -- one never leaves it, the other never enters it.
    #[test]
    fn test_descent_opens_exactly_one_more_run() {
        let dir = tempfile::tempdir().expect("tempdir");
        let input = build_descent_bam(dir.path(), 300);
        let output = dir.path().join("output.bam");

        let stats = run_sort(SortOrder::Coordinate, &input, &output, SPILL_MEMORY);

        assert!(
            stats.runs_written >= 2,
            "a mid-stream descent must close the open run, got {} run(s)",
            stats.runs_written
        );
        assert_eq!(stats.total_records, stats.output_records, "sort must not lose records");
        assert_eq!(count_bam_records(&output), 600);
        assert_sorted_in(SortOrder::Coordinate, &output);
    }

    /// Consolidation can merge the open run away, leaving `RunFormer` holding a
    /// path that no longer means what it did. Appending there would duplicate or
    /// drop records, so exercise it with consolidation actually enabled --
    /// every other run-formation test disables it.
    ///
    /// Parameterized over every order because each reaches
    /// [`RunFormer::place`]'s consolidation check with a different key type
    /// through the same generic [`chunk_bounds`].
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname_lex(SortOrder::Queryname(QuerynameComparator::Lexicographic))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_consolidation_while_a_run_is_open(#[case] order: SortOrder) {
        let dir = tempfile::tempdir().expect("tempdir");
        let input = build_descent_bam(dir.path(), 400);
        let output = dir.path().join("output.bam");

        let stats = RawExternalSorter::new(order)
            .memory_limit(8 * 1024)
            .max_temp_files(2) // small enough that consolidation can swallow the open run
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort with consolidation should succeed");

        assert_eq!(stats.total_records, stats.output_records, "{order:?}: lost records");
        assert_eq!(count_bam_records(&output), 800, "{order:?}: wrong record count");
        assert_sorted_in(order, &output);
    }

    /// Appending must work for both spill codecs. Zstd skips `ZSPILL_MAGIC` on
    /// append -- a second magic mid-file would be read as frame data -- and that
    /// branch is otherwise untested.
    #[rstest::rstest]
    #[case::bgzf(crate::codec::SpillCodec::Bgzf, 0)]
    #[case::zstd(crate::codec::SpillCodec::Zstd, 1)]
    fn test_append_across_spill_codecs(
        #[case] codec: crate::codec::SpillCodec,
        #[case] temp_compression: u32,
    ) {
        let dir = tempfile::tempdir().expect("tempdir");
        let input = build_ordered_bam(dir.path(), 300);
        let output = dir.path().join("output.bam");

        let stats = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024)
            .max_temp_files(0)
            .threads(1)
            .spill_codec(codec)
            .temp_compression(temp_compression)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        assert_eq!(stats.runs_written, 1, "sorted input should collapse to one run");
        assert_eq!(stats.total_records, stats.output_records, "sort must not lose records");
        assert_eq!(count_bam_records(&output), 600);
        assert_sorted_in(SortOrder::Coordinate, &output);
    }

    /// The indexed path got run formation too, so cover appending there: the
    /// merge sees one source instead of many while still emitting a BAI.
    #[test]
    fn test_appending_run_with_write_index() {
        let dir = tempfile::tempdir().expect("tempdir");
        let input = build_ordered_bam(dir.path(), 300);
        let output = dir.path().join("output.bam");

        let stats = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024)
            .max_temp_files(0)
            .threads(1)
            .write_index(true)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("indexed sort should succeed");

        assert_eq!(stats.runs_written, 1, "sorted input should collapse to one run");
        assert_eq!(stats.total_records, stats.output_records, "sort must not lose records");
        assert_eq!(count_bam_records(&output), 600);
        assert!(output.with_extension("bam.bai").exists(), "index should be written");
    }

    /// Verifies that sort with many chunks exercises the pool-integrated merge
    /// readers during the final merge phase (not just consolidation).
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_sort_many_chunks_with_semaphore(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 200;
        let mut builder = SamBuilder::new();
        // Emit descending coordinates so the input is NOT already in coordinate
        // order. Natural run formation appends a chunk to the open run whenever
        // its keys all sort after it, so ascending input collapses to a single
        // run -- and a single source exercises neither the k-way merge nor the
        // reader semaphore this test exists to cover.
        for i in 0..num_pairs {
            let descending = num_pairs - 1 - i;
            let _ = builder
                .add_pair()
                .name(&format!("read{descending:05}"))
                .start1(descending * 200 + 1)
                .start2(descending * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Small memory limit + many records = guaranteed spill to multiple chunks
        // (no consolidation, so the semaphore must cap concurrent readers).
        // 32 KiB is small enough to force several spills across 400 records but
        // large enough to avoid creating hundreds of tiny chunks that saturate
        // OS threads on the CI runner and cause timeouts.
        let stats = RawExternalSorter::new(sort_order)
            .memory_limit(32 * 1024)
            .max_temp_files(0) // disable consolidation
            .threads(2) // semaphore allows 2 concurrent readers
            .spill_codec(crate::codec::SpillCodec::Bgzf) // level 0 = stored mode (zstd has no level 0)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        assert!(
            stats.runs_written >= 2,
            "expected multiple chunks to exercise merge, got {}",
            stats.runs_written
        );

        let expected = (num_pairs * 2) as u64;
        let observed = count_bam_records(&output);
        assert_eq!(observed, expected, "semaphore-capped merge lost data");
    }

    // ========================================================================
    // Sub-array merge source tests
    // ========================================================================

    /// Verifies that multi-threaded sort produces the same output as single-threaded
    /// sort for each sort order.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_sort_sub_arrays_match_single_thread(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 50;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output_st = dir.path().join("output_1t.bam");
        let output_mt = dir.path().join("output_2t.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Sort single-threaded
        RawExternalSorter::new(sort_order)
            .memory_limit(16 * 1024) // force spill (50 pairs × ~300B ≈ 15KB) so merge path is exercised
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf) // level 0 = stored mode (zstd has no level 0)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_st)
            .expect("sort should succeed");

        // Sort multi-threaded
        RawExternalSorter::new(sort_order)
            .memory_limit(16 * 1024)
            .threads(2)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_mt)
            .expect("sort should succeed");

        let names_st = collect_read_names(&output_st);
        let names_mt = collect_read_names(&output_mt);

        let expected = num_pairs * 2;
        assert_eq!(names_st.len(), expected, "single-thread record count mismatch");
        assert_eq!(names_mt.len(), expected, "multi-thread record count mismatch");

        // Both outputs must have the same read names in the same order
        assert_eq!(names_st, names_mt, "multi-thread sort order differs from single-thread");
    }

    /// Verifies that the in-memory-only path (no spill to disk) works correctly.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::default()))]
    #[case::queryname_natural(SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_sort_sub_arrays_in_memory_only(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        let num_pairs = 20;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("failed to write BAM");

        // Large memory limit so everything stays in memory (no spill chunks).
        RawExternalSorter::new(sort_order)
            .memory_limit(10 * 1024 * 1024)
            .threads(2)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");

        let expected = (num_pairs * 2) as u64;
        let observed = count_bam_records(&output);
        assert_eq!(observed, expected, "sort lost data");
    }

    /// A record larger than one BGZF block must survive a spill followed by a
    /// pool-integrated merge.
    ///
    /// When records are spilled, the writer pre-flushes so a normal record
    /// never straddles a BGZF block boundary — only an oversized record
    /// (written via `write_chunked`) spans blocks. On merge that oversized
    /// record is the one case that cannot be borrowed in place from a single
    /// decompressed block, so it exercises the reassemble-into-scratch slow
    /// path that the zero-copy fast path falls back to. This guards that
    /// fast/slow split: the fast path handles the ~200 normal reads, the slow
    /// path the one oversized read, and every byte must round-trip.
    ///
    /// The split is asserted, not just described. It was described here (and in
    /// `parse_next_record`) while the keyed orders in fact reassembled *every*
    /// record: the borrow was implemented only for the `EMBEDDED_IN_RECORD`
    /// keys, so a template-coordinate merge paid a full-record copy 776,167,078
    /// times out of 776,167,078 on the one thread that touches every record.
    /// Prose cannot catch that regression; a counter can.
    ///
    /// The counters are process-global, so under `cargo ci-test` (nextest,
    /// process-per-test) these deltas are exactly this test's, while under a
    /// thread-parallel `cargo test` a concurrent test can only inflate them.
    /// Asserting lower bounds is therefore sound in both.
    #[rstest::rstest]
    #[case::coordinate(SortOrder::Coordinate)]
    #[case::template_coordinate(SortOrder::TemplateCoordinate)]
    fn test_oversized_record_survives_spill_and_pool_merge(#[case] sort_order: SortOrder) {
        use fgumi_sam::SamBuilder;

        // Well over one ~64 KB BGZF block; default refs are 200 MB so a mapped
        // read this long is still in-bounds. The bases are position-dependent
        // rather than a uniform run so they can serve as their own oracle: a
        // scratch reassembly that shuffled or duplicated bytes within the
        // record would preserve the length (and round-trip a run of `A`s) but
        // must not survive an exact base-by-base comparison below.
        let big_len = 70_000usize;
        let big_seq: String = {
            let mut state = 0x2545_f491_4f6c_dd1d_u64;
            (0..big_len)
                .map(|_| {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    ['A', 'C', 'G', 'T'][((state >> 33) & 0b11) as usize]
                })
                .collect()
        };

        let mut builder = SamBuilder::new();
        for i in 0..100 {
            let _ = builder.add_frag().name(&format!("n{i:04}")).contig(0).start(i + 1).build();
        }
        let _ = builder.add_frag().name("oversized").contig(0).start(50).bases(&big_seq).build();
        for i in 0..100 {
            let _ = builder.add_frag().name(&format!("m{i:04}")).contig(0).start(i + 1).build();
        }

        let dir = tempfile::tempdir().expect("tempdir");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("write_bam");

        let borrowed_before = RECORD_BORROWED.load(std::sync::atomic::Ordering::Relaxed);
        let reassembled_before = RECORD_REASSEMBLED.load(std::sync::atomic::Ordering::Relaxed);

        // Memory limit below the oversized record forces it into its own spill
        // chunk (so it is read back through the PoolDisk merge, not from memory);
        // two threads give a genuine multi-source merge.
        let stats = RawExternalSorter::new(sort_order)
            .memory_limit(64 * 1024)
            .threads(2)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("sort should succeed");
        assert!(
            stats.runs_written > 0,
            "test must spill to disk so the oversized record is read back through the \
             PoolDisk merge (slow path); got runs_written = 0"
        );

        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&output).expect("open output");
        let _header = reader.read_header().expect("read header");
        let mut count = 0u64;
        let mut oversized_bases: Option<Vec<u8>> = None;
        for result in reader.records() {
            let record = result.expect("read record");
            let name: &[u8] = record.name().expect("read name").as_ref();
            if name == b"oversized" {
                oversized_bases = Some(record.sequence().iter().collect());
            }
            count += 1;
        }
        assert_eq!(count, 201, "spill + pool merge lost records");

        // The fast/slow split the doc comment above describes.
        let borrowed = RECORD_BORROWED.load(std::sync::atomic::Ordering::Relaxed) - borrowed_before;
        let reassembled =
            RECORD_REASSEMBLED.load(std::sync::atomic::Ordering::Relaxed) - reassembled_before;
        assert!(
            reassembled >= 1,
            "the oversized record spans blocks and must reassemble into scratch, exercising the \
             slow path; got {reassembled} reassembled for {sort_order:?}"
        );
        // Measured: 100 of the 201 records come back through the PoolDisk merge
        // (the rest are served from the in-memory chunk, which never parses).
        // With the keyed borrow removed this is 0 for TemplateCoordinate and
        // unchanged at 100 for Coordinate, which is what makes this an
        // order-sensitive guard rather than a restatement of the record count.
        assert!(
            borrowed >= 90,
            "normal records are borrowed in place, so this merge must borrow ~100 of them -- 0 \
             means the zero-copy path never fires for {sort_order:?}, which is how the keyed \
             orders silently copied every record; got {borrowed} borrowed"
        );
        let oversized_bases = oversized_bases.expect("oversized read must survive the merge");
        assert_eq!(
            oversized_bases.len(),
            big_len,
            "oversized read lost bases across the spill + pool merge"
        );
        assert_eq!(
            oversized_bases,
            big_seq.as_bytes(),
            "oversized read's bases were corrupted across the spill + pool merge"
        );
    }

    // ========================================================================
    // merge_bams tests
    // ========================================================================

    /// Helper: create a `SamBuilder` with `num_pairs` pairs at non-overlapping positions,
    /// write an unsorted BAM, sort it with the given order, and return the sorted path.
    /// `start_offset` shifts all positions so different inputs have distinct read names/positions.
    fn create_sorted_bam(
        dir: &Path,
        prefix: &str,
        num_pairs: usize,
        start_offset: usize,
        sort_order: SortOrder,
    ) -> (PathBuf, Vec<String>) {
        let mut builder = SamBuilder::new();
        let mut names = Vec::with_capacity(num_pairs);
        for i in 0..num_pairs {
            let name = format!("{prefix}_read{i:04}");
            names.push(name.clone());
            let _ = builder
                .add_pair()
                .name(&name)
                .start1((start_offset + i * 200) + 1)
                .start2((start_offset + i * 200) + 101)
                .build();
        }
        let unsorted = dir.join(format!("{prefix}_unsorted.bam"));
        let sorted = dir.join(format!("{prefix}_sorted.bam"));
        builder.write_bam(&unsorted).expect("failed to write BAM");
        RawExternalSorter::new(sort_order)
            .output_compression(0)
            .sort(&unsorted, &sorted)
            .expect("sort should succeed");
        (sorted, names)
    }

    // ========================================================================
    // sort_records (in-process record stream) tests
    // ========================================================================

    /// Read a BAM into its header plus every record, for replaying through
    /// [`RawExternalSorter::sort_records`].
    fn read_bam_records(path: &Path) -> (Header, Vec<fgumi_raw_bam::RawRecord>) {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, header) =
            create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        (header, RawReadAheadReader::new(reader).collect())
    }

    /// Build an unsorted BAM of `num_pairs` pairs whose coordinates deliberately
    /// do not follow input order, so a sort has real work to do.
    fn build_unsorted_bam(dir: &Path, num_pairs: usize) -> PathBuf {
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            // Interleave high and low coordinates so input order != sort order.
            let start = if i % 2 == 0 { (num_pairs - i) * 200 + 1 } else { i * 200 + 1 };
            let _ = builder
                .add_pair()
                .name(&format!("read{i:04}"))
                .start1(start)
                .start2(start + 100)
                .build();
        }
        let path = dir.join("unsorted.bam");
        builder.write_bam(&path).expect("failed to write BAM");
        path
    }

    /// `sort_records` must produce byte-identical output to `sort` over the same
    /// records in the same order, for every sort order and on both the
    /// fits-in-memory and spill-and-merge paths.
    ///
    /// This is the invariant callers rely on when they replace an
    /// "write unsorted BAM, then sort it" round-trip with a direct stream:
    /// removing the intermediate file must not change a single output byte.
    #[rstest::rstest]
    fn test_sort_records_matches_sort_from_file(
        #[values(
            SortOrder::Coordinate,
            SortOrder::Queryname(QuerynameComparator::Lexicographic),
            SortOrder::Queryname(QuerynameComparator::Natural),
            SortOrder::TemplateCoordinate
        )]
        sort_order: SortOrder,
        // 8 MiB fits the whole 100-pair set in memory (no spill); 16 KiB forces
        // several spills and a k-way merge.
        #[values(8 * 1024 * 1024, 16 * 1024)] memory_limit: usize,
    ) {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = build_unsorted_bam(dir.path(), 100);
        let (header, records) = read_bam_records(&input);

        let from_file = dir.path().join("from_file.bam");
        let from_stream = dir.path().join("from_stream.bam");

        let sorter = || {
            RawExternalSorter::new(sort_order)
                .memory_limit(memory_limit)
                .threads(2)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0)
        };

        let file_stats = sorter().sort(&input, &from_file).expect("file sort should succeed");
        let stream_stats = sorter()
            .sort_records(records.into_iter().map(Ok), &header, &from_stream)
            .expect("stream sort should succeed");

        assert_eq!(file_stats.total_records, stream_stats.total_records);
        assert_eq!(file_stats.output_records, stream_stats.output_records);
        assert_eq!(
            file_stats.runs_written, stream_stats.runs_written,
            "stream path should spill exactly like the file path"
        );
        assert_eq!(
            std::fs::read(&from_file).expect("read file output"),
            std::fs::read(&from_stream).expect("read stream output"),
            "sort_records output must be byte-identical to sort()"
        );
    }

    /// An empty stream produces a header-only BAM, matching `sort` on an empty
    /// input rather than erroring.
    #[test]
    fn test_sort_records_empty_stream_writes_header_only() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = build_unsorted_bam(dir.path(), 1);
        let (header, _) = read_bam_records(&input);
        let output = dir.path().join("empty.bam");

        let stats = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .sort_records(std::iter::empty(), &header, &output)
            .expect("empty stream should succeed");

        assert_eq!(stats.total_records, 0);
        assert_eq!(count_bam_records(&output), 0);
    }

    /// A producer error must abort the sort and propagate, leaving no output
    /// file behind — a truncated "successful" result would be worse than a
    /// failure, since callers would ship a silently short BAM.
    ///
    /// Covered on both the fits-in-memory path and after a spill, since the
    /// error is only observed once the stream is drained.
    #[rstest::rstest]
    #[case::before_spill(8 * 1024 * 1024)]
    #[case::after_spill(16 * 1024)]
    fn test_sort_records_producer_error_aborts(#[case] memory_limit: usize) {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = build_unsorted_bam(dir.path(), 100);
        let (header, records) = read_bam_records(&input);
        let output = dir.path().join("aborted.bam");

        // Emit every record, then fail instead of ending cleanly.
        let stream = records
            .into_iter()
            .map(Ok)
            .chain(std::iter::once(Err(anyhow::anyhow!("producer exploded"))));

        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(memory_limit)
            .threads(2)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort_records(stream, &header, &output)
            .expect_err("producer error should abort the sort");

        assert!(
            format!("{err:#}").contains("producer exploded"),
            "producer error should be propagated, got: {err:#}"
        );
        assert!(!output.exists(), "aborted sort must not leave an output BAM behind");
    }

    /// A producer that fails before yielding anything aborts too, rather than
    /// falling through the empty-input path and writing a header-only BAM.
    #[test]
    fn test_sort_records_error_on_first_record_aborts() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let input = build_unsorted_bam(dir.path(), 1);
        let (header, _) = read_bam_records(&input);
        let output = dir.path().join("aborted_first.bam");

        let stream = std::iter::once(Err(anyhow::anyhow!("failed before first record")));

        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .sort_records(stream, &header, &output)
            .expect_err("producer error should abort the sort");

        assert!(
            format!("{err:#}").contains("failed before first record"),
            "producer error should be propagated, got: {err:#}"
        );
        assert!(!output.exists(), "aborted sort must not leave an output BAM behind");
    }

    /// Collect all read names from a BAM file as strings.
    fn collect_read_names(path: &Path) -> Vec<String> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader)
            .map(|rec| {
                let name_bytes = fgumi_raw_bam::RawRecordView::new(rec.as_ref()).read_name();
                String::from_utf8(name_bytes.to_vec()).expect("read name should be valid UTF-8")
            })
            .collect()
    }

    /// Collect (`ref_id`, pos) tuples for every record in a BAM.
    fn collect_positions(path: &Path) -> Vec<(i32, i32)> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("failed to create raw BAM reader");
        RawReadAheadReader::new(reader)
            .map(|rec| {
                let bytes = rec.as_ref();
                {
                    let v = fgumi_raw_bam::RawRecordView::new(bytes);
                    (v.ref_id(), v.pos())
                }
            })
            .collect()
    }

    /// Helper to build a merge header from the `SamBuilder` default header.
    fn default_merge_header() -> Header {
        SamBuilder::new().header.clone()
    }

    #[test]
    fn test_merge_bams_coordinate_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, SortOrder::Coordinate);
        let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, SortOrder::Coordinate);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Coordinate)
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        // 10 pairs * 2 records * 2 inputs = 40
        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify coordinate sort order: (ref_id, pos) is non-decreasing
        let positions = collect_positions(&merged);
        for w in positions.windows(2) {
            assert!(w[0] <= w[1], "coordinate sort violated: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn test_merge_bams_template_coordinate_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, SortOrder::TemplateCoordinate);
        let (bam_b, _) =
            create_sorted_bam(dir.path(), "b", 10, 10_000, SortOrder::TemplateCoordinate);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);
    }

    #[test]
    fn test_merge_bams_queryname_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(
            dir.path(),
            "a",
            10,
            0,
            SortOrder::Queryname(QuerynameComparator::default()),
        );
        let (bam_b, _) = create_sorted_bam(
            dir.path(),
            "b",
            10,
            10_000,
            SortOrder::Queryname(QuerynameComparator::default()),
        );

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify queryname sort order: read names are non-decreasing
        let names = collect_read_names(&merged);
        for w in names.windows(2) {
            assert!(w[0] <= w[1], "queryname sort violated: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn test_merge_bams_single_input() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 15, 0, SortOrder::Coordinate);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Coordinate)
            .output_compression(0)
            .merge_bams(&[bam_a], &header, &merged)
            .expect("sort should succeed");

        // 15 pairs * 2 = 30
        assert_eq!(count, 30);
        assert_eq!(count_bam_records(&merged), 30);
    }

    #[test]
    fn test_merge_bams_preserves_all_records() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let (bam_a, names_a) = create_sorted_bam(
            dir.path(),
            "a",
            5,
            0,
            SortOrder::Queryname(QuerynameComparator::default()),
        );
        let (bam_b, names_b) = create_sorted_bam(
            dir.path(),
            "b",
            5,
            10_000,
            SortOrder::Queryname(QuerynameComparator::default()),
        );

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::default()))
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("sort should succeed");

        let merged_names: std::collections::HashSet<String> =
            collect_read_names(&merged).into_iter().collect();

        // Every expected name (from both inputs) must appear in the merged output
        for name in names_a.iter().chain(names_b.iter()) {
            assert!(merged_names.contains(name), "read name {name:?} missing from merged output");
        }
    }

    #[test]
    fn test_merge_bams_many_inputs() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let k = 8;
        let pairs_per_input = 5;
        let mut inputs = Vec::with_capacity(k);

        for i in 0..k {
            let (bam, _) = create_sorted_bam(
                dir.path(),
                &format!("in{i}"),
                pairs_per_input,
                i * 50_000,
                SortOrder::Coordinate,
            );
            inputs.push(bam);
        }

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(SortOrder::Coordinate)
            .output_compression(0)
            .merge_bams(&inputs, &header, &merged)
            .expect("sort should succeed");

        let expected = (k * pairs_per_input * 2) as u64; // 8 * 5 * 2 = 80
        assert_eq!(count, expected);
        assert_eq!(count_bam_records(&merged), expected);

        // Verify coordinate sort order
        let positions = collect_positions(&merged);
        for w in positions.windows(2) {
            assert!(w[0] <= w[1], "coordinate sort violated with k={k}: {:?} > {:?}", w[0], w[1]);
        }
    }

    #[test]
    fn test_merge_bams_queryname_natural_sort() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let nat = SortOrder::Queryname(QuerynameComparator::Natural);
        let (bam_a, _) = create_sorted_bam(dir.path(), "a", 10, 0, nat);
        let (bam_b, _) = create_sorted_bam(dir.path(), "b", 10, 10_000, nat);

        let merged = dir.path().join("merged.bam");
        let header = default_merge_header();
        let count = RawExternalSorter::new(nat)
            .output_compression(0)
            .merge_bams(&[bam_a, bam_b], &header, &merged)
            .expect("merge should succeed");

        assert_eq!(count, 40);
        assert_eq!(count_bam_records(&merged), 40);

        // Verify natural queryname order: names are non-decreasing
        let names = collect_read_names(&merged);
        for w in names.windows(2) {
            assert!(w[0] <= w[1], "natural queryname sort violated: {:?} > {:?}", w[0], w[1]);
        }
    }

    // ========================================================================
    // SortPhaseTimer tests
    // ========================================================================

    #[test]
    fn test_sort_phase_timer_all_methods() {
        let mut timer = SortPhaseTimer::new();
        assert!(timer.overall_start.is_some());
        assert!(timer.read_span_start.is_some());

        // end_read_span accumulates and clears the span
        let elapsed = timer.end_read_span();
        assert!(timer.read_secs >= 0.0);
        assert!(timer.read_span_start.is_none());
        let _ = elapsed;

        // end_read_span with no active span → zero, no panic
        let elapsed2 = timer.end_read_span();
        assert_eq!(elapsed2, std::time::Duration::ZERO);
        assert!(timer.read_secs >= 0.0); // unchanged

        // begin_read_span restores the span
        timer.begin_read_span();
        assert!(timer.read_span_start.is_some());

        // time_sort measures elapsed
        timer.time_sort(|| {});
        assert!(timer.sort_secs >= 0.0);

        // time_spill_write: counts spills and returns the closure result
        let result = timer.time_spill_write(|| Ok::<u32, anyhow::Error>(42));
        assert_eq!(result.unwrap(), 42);
        assert_eq!(timer.spill_count, 1);
        assert!(timer.spill_write_secs >= 0.0);

        // record_spill_size: adds file size to total
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("spill.bin");
        std::fs::write(&path, b"hello world").expect("write");
        timer.record_spill_size(&path);
        assert_eq!(timer.total_spill_bytes, 11);

        // record_spill_size with nonexistent path: no-op, no panic
        timer.record_spill_size(&dir.path().join("nonexistent.bin"));
        assert_eq!(timer.total_spill_bytes, 11);

        // time_consolidate: fast closure (<1ms) → consolidate_count stays 0
        timer.time_consolidate(|| Ok(())).expect("consolidate ok");
        assert_eq!(timer.consolidate_count, 0);

        // time_consolidate: slow closure (>1ms) → counted
        timer
            .time_consolidate(|| {
                std::thread::sleep(std::time::Duration::from_millis(10));
                Ok(())
            })
            .expect("consolidate ok");
        assert_eq!(timer.consolidate_count, 1);
        assert!(timer.consolidate_secs > 0.0);

        // time_merge
        timer.time_merge(|| Ok::<(), anyhow::Error>(())).expect("merge ok");
        assert!(timer.merge_secs >= 0.0);

        // time_write_output
        timer.time_write_output(|| Ok(())).expect("write ok");
        assert!(timer.write_output_secs >= 0.0);

        // log_summary must not panic (output goes to log sink). `consolidate_count`
        // is 1 here, so the consolidation branch is exercised too.
        timer.log_summary(4, 4, 64);
    }

    // ========================================================================
    // Chunk boundary alignment tests
    // ========================================================================

    /// Regression test: `split_off` from the tail must align with `par_chunks_mut`
    /// boundaries so that every emitted chunk is individually sorted.
    ///
    /// The old code split fixed-size chunks from the tail, which crossed sorted
    /// chunk boundaries when `n % chunk_size != 0`.
    #[test]
    #[allow(clippy::cast_possible_truncation)]
    fn test_split_off_aligns_with_par_chunks_mut() {
        use rayon::prelude::*;

        let pool =
            rayon::ThreadPoolBuilder::new().num_threads(4).build().expect("build rayon pool");

        // Test a variety of (n, threads) combinations that exercise the tail case
        for &(n, threads) in &[(10, 3), (11, 4), (17, 3), (100, 7), (1000, 8), (13, 5)] {
            let mut entries: Vec<(u64, Vec<u8>)> =
                (0..n).rev().map(|i| (i as u64, vec![i as u8])).collect();

            let chunk_size = entries.len().div_ceil(threads);

            pool.install(|| {
                entries.par_chunks_mut(chunk_size).for_each(|chunk| {
                    chunk.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                });
            });

            // Apply the same split logic as production code
            let mut remaining = std::mem::take(&mut entries);
            let num_chunks = remaining.len().div_ceil(chunk_size);
            let mut chunks: Vec<Vec<(u64, Vec<u8>)>> = Vec::with_capacity(num_chunks);
            let tail_len = remaining.len() % chunk_size;
            if tail_len != 0 {
                let split_at = remaining.len() - tail_len;
                chunks.push(remaining.split_off(split_at));
            }
            while !remaining.is_empty() {
                let split_at = remaining.len().saturating_sub(chunk_size);
                chunks.push(remaining.split_off(split_at));
            }
            chunks.reverse();

            // Every chunk must be individually sorted
            for (ci, chunk) in chunks.iter().enumerate() {
                for i in 1..chunk.len() {
                    assert!(
                        chunk[i - 1].0 <= chunk[i].0,
                        "n={n} threads={threads}: chunk {ci} not sorted at index {i} \
                         ({} > {})",
                        chunk[i - 1].0,
                        chunk[i].0,
                    );
                }
            }

            // Total record count must match
            let total: usize = chunks.iter().map(Vec::len).sum();
            assert_eq!(total, n, "n={n} threads={threads}: total mismatch");
        }
    }

    /// End-to-end sort across two temp directories. Forces multiple spill
    /// chunks and verifies the output is still correct. The round-robin
    /// distribution itself is covered by `TmpDirAllocator`'s unit tests; this
    /// test proves the plumbing from `temp_dirs(...)` through to the sort
    /// fns works and that multi-dir mode produces byte-identical output to
    /// single-dir mode.
    #[test]
    fn test_sort_with_two_temp_dirs_matches_single_dir() {
        use fgumi_sam::SamBuilder;

        let num_pairs = 200;
        let mut builder = SamBuilder::new();
        // Descending coordinates so run formation does not collapse the input
        // to one run; this test compares multi-chunk spill across temp dirs.
        for i in 0..num_pairs {
            let descending = num_pairs - 1 - i;
            let _ = builder
                .add_pair()
                .name(&format!("read{descending:05}"))
                .start1(descending * 200 + 1)
                .start2(descending * 200 + 101)
                .build();
        }

        let workdir = tempfile::tempdir().expect("workdir");
        let input = workdir.path().join("input.bam");
        let output_multi = workdir.path().join("output_multi.bam");
        let output_single = workdir.path().join("output_single.bam");
        builder.write_bam(&input).expect("write bam");

        let tmp_a = tempfile::tempdir().expect("tmp a");
        let tmp_b = tempfile::tempdir().expect("tmp b");

        // 8 KiB memory limit forces several spills across the dir rotation.
        let stats_multi = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024)
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf) // level 0 = stored mode (zstd has no level 0)
            .temp_compression(0)
            .output_compression(0)
            .temp_dirs(vec![tmp_a.path().to_path_buf(), tmp_b.path().to_path_buf()])
            .sort(&input, &output_multi)
            .expect("multi-dir sort should succeed");

        assert!(stats_multi.runs_written >= 2, "expected multiple spill runs");

        RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024)
            .threads(1)
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output_single)
            .expect("single-dir sort should succeed");

        let names_multi = collect_read_names(&output_multi);
        let names_single = collect_read_names(&output_single);
        assert_eq!(names_multi.len(), num_pairs * 2, "record count mismatch");
        assert_eq!(
            names_multi, names_single,
            "multi-dir and single-dir sort produced different record orders"
        );
    }

    /// The indexed coordinate sort derives the sidecar path in two independent
    /// places: the in-memory path (everything fits under the memory limit) and
    /// the spilling merge path. This covers the in-memory one; the spilling
    /// sibling below covers the other. Both must land on the samtools sidecar
    /// path even when the output is not named `*.bam`.
    #[rstest::rstest]
    #[case::dot_bam_output("output.bam")]
    #[case::non_bam_output("output.sorted")]
    fn test_sort_coordinate_with_index_in_memory_writes_sidecar(#[case] output_name: &str) {
        use fgumi_sam::SamBuilder;

        let mut builder = SamBuilder::new();
        for i in 0..20 {
            let _ = builder
                .add_pair()
                .name(&format!("read{i:05}"))
                .start1(i * 200 + 1)
                .start2(i * 200 + 101)
                .build();
        }

        let workdir = tempfile::tempdir().expect("workdir");
        let input = workdir.path().join("input.bam");
        let output = workdir.path().join(output_name);
        builder.write_bam(&input).expect("write bam");

        let stats = RawExternalSorter::new(SortOrder::Coordinate)
            // Generous limit so nothing spills and the in-memory branch runs.
            .memory_limit(64 * 1024 * 1024)
            .threads(1)
            .write_index(true)
            .output_compression(0)
            .sort(&input, &output)
            .expect("indexed coordinate sort should succeed");

        assert_eq!(stats.runs_written, 0, "expected no spill for the in-memory path");
        assert_eq!(collect_read_names(&output).len(), 40, "record count mismatch");

        let bai = fgumi_bam_io::bai_sidecar_path(&output);
        assert!(bai.exists(), "index should exist at the sidecar path {}", bai.display());

        let stale = output.with_extension("bam.bai");
        if stale != bai {
            assert!(!stale.exists(), "no index at the extension-replaced path {}", stale.display());
        }
    }

    /// `output_name` covers both the `.bam` case (where the old
    /// `with_extension("bam.bai")` happened to be right) and a name that does
    /// not end in `.bam` (where it silently wrote the index under a different
    /// basename entirely).
    #[rstest::rstest]
    #[case::dot_bam_output("output.bam")]
    #[case::non_bam_output("output.sorted")]
    #[case::no_extension_output("output")]
    fn test_sort_coordinate_with_index_spilled_preserves_records(#[case] output_name: &str) {
        use fgumi_sam::SamBuilder;

        // Enough pairs that an 8 KiB memory limit forces multiple spills,
        // so the index-emitting merge runs over Disk sources + the in-memory
        // chunk simultaneously (the mixed-source path the refactor must keep
        // byte-identical).
        // Descending coordinates: the point is a merge over several Disk sources
        // plus the in-memory chunk, and natural run formation would fold
        // already-ascending input into one run, leaving nothing to merge.
        let num_pairs = 300;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let descending = num_pairs - 1 - i;
            let _ = builder
                .add_pair()
                .name(&format!("read{descending:05}"))
                .start1(descending * 200 + 1)
                .start2(descending * 200 + 101)
                .build();
        }

        let workdir = tempfile::tempdir().expect("workdir");
        let input = workdir.path().join("input.bam");
        let output = workdir.path().join(output_name);
        builder.write_bam(&input).expect("write bam");

        let stats = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(8 * 1024) // tiny → forces several spills
            .threads(1)
            .write_index(true) // routes through sort_coordinate_with_index → merge_chunks_with_index
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("indexed coordinate sort should succeed");

        assert!(stats.runs_written >= 2, "expected multiple spill runs");

        // All records preserved.
        let names = collect_read_names(&output);
        assert_eq!(names.len(), num_pairs * 2, "record count mismatch");

        // The index must land at exactly the samtools sidecar path -- append
        // `.bai` to the full output path. Asserting the exact path (rather than
        // accepting either naming) is the point: the old
        // `with_extension("bam.bai")` put it under a different basename for any
        // output not ending in `.bam`, where no indexed reader would find it.
        let bai = fgumi_bam_io::bai_sidecar_path(&output);
        assert!(
            bai.exists(),
            "index should exist at the sidecar path {} for output {}",
            bai.display(),
            output.display(),
        );

        // And nowhere else: a stray index under the rewritten basename would
        // mean the old expression is still in play somewhere.
        let stale = output.with_extension("bam.bai");
        if stale != bai {
            assert!(
                !stale.exists(),
                "no index should be written to the extension-replaced path {}",
                stale.display(),
            );
        }
    }

    /// A record spanning multiple BGZF blocks must resolve correct virtual
    /// offsets through the pool-integrated indexed merge.
    ///
    /// Normal records never straddle a block boundary (the spill writer
    /// pre-flushes), so an oversized record is the only thing that drives
    /// `BaiBuilder`'s cross-block `compute_end_vpos` walk while
    /// `merge_chunks_with_index` builds the index — a wrong block stride there
    /// would fail the index build or emit a corrupt `.bai`. This is the
    /// non-samtools counterpart to the `#[ignore]` samtools cross-check: it
    /// forces the oversized record through a spilled indexed merge, then loads
    /// the produced `.bai` back and checks it is well-formed.
    #[test]
    fn test_write_index_oversized_record_builds_loadable_bai() {
        use fgumi_sam::SamBuilder;

        let big_seq = "A".repeat(70_000); // > one ~64 KB BGZF block
        let mut builder = SamBuilder::new();
        for i in 0..80 {
            let _ = builder.add_frag().name(&format!("n{i:04}")).contig(0).start(i + 1).build();
        }
        let _ = builder.add_frag().name("oversized").contig(0).start(40).bases(&big_seq).build();
        for i in 0..80 {
            let _ = builder.add_frag().name(&format!("m{i:04}")).contig(0).start(i + 1).build();
        }
        let expected_refs = builder.header.reference_sequences().len();

        let dir = tempfile::tempdir().expect("tempdir");
        let input = dir.path().join("input.bam");
        let output = dir.path().join("output.bam");
        builder.write_bam(&input).expect("write_bam");

        let stats = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(64 * 1024) // below the oversized record → forces it into a spill chunk
            .threads(2) // multi-source pool-integrated indexed merge
            .write_index(true) // routes through merge_chunks_with_index
            .spill_codec(crate::codec::SpillCodec::Bgzf)
            .temp_compression(0)
            .output_compression(0)
            .sort(&input, &output)
            .expect("indexed coordinate sort should succeed");
        assert!(stats.runs_written > 0, "test must spill so the indexed pool merge runs");

        // Every record (including the oversized one) survives the merge.
        assert_eq!(count_bam_records(&output), 161);

        // The .bai the indexed merge produced must load cleanly and cover every
        // reference: a corrupt cross-block virtual offset would have failed the
        // index build above or produced an unreadable sidecar here.
        let bai = fgumi_bam_io::bai_sidecar_path(&output);
        let index = noodles::bam::bai::fs::read(&bai).expect("BAI must be loadable");
        assert_eq!(
            index.reference_sequences().len(),
            expected_refs,
            "index should cover every reference sequence"
        );
    }

    /// Per-phase thread counts are a scheduling knob only: whatever the split,
    /// the sorted output must be byte-identical to a plain `--threads` run.
    /// Also exercises the asymmetric cases in both directions, since Phase 1 and
    /// Phase 2 size the shared pool differently (`max` of the two) and the cap
    /// is flipped at the ingest/merge boundary -- a capped worker that failed to
    /// drain held items would strand records and change the output.
    ///
    /// `write_index` is crossed in via `#[values]` so both the plain and the
    /// `--write-index` output paths are covered for every split -- the indexed
    /// path routes through separate writers (`sort_coordinate_with_index` and
    /// `merge_chunks_with_index`) that must honor `merge_threads` the same way.
    /// `memory_limit` is crossed in to hit both indexed writer sites: the small
    /// limit spills and runs `merge_chunks_with_index`; the large limit keeps
    /// everything in memory and runs the in-memory indexed writer. Output is
    /// byte-identical regardless of the writer's compression-thread count, so
    /// this guards the path against regression rather than asserting the count.
    #[rstest::rstest]
    #[case::narrow_sort_wide_merge(1, 4)]
    #[case::wide_sort_narrow_merge(4, 1)]
    #[case::both_narrow(1, 1)]
    #[case::both_wide(3, 3)]
    fn test_per_phase_threads_produce_identical_output(
        #[case] sort_threads: usize,
        #[case] merge_threads: usize,
        #[values(false, true)] write_index: bool,
        #[values(16 * 1024, 64 * 1024 * 1024)] memory_limit: usize,
    ) {
        use fgumi_sam::SamBuilder;

        // Enough records that, at the small memory limit, Phase 1 spills several
        // chunks and Phase 2 runs a real k-way merge across them; at the large
        // limit everything fits in memory and the in-memory write path runs.
        let num_pairs = 400;
        let mut builder = SamBuilder::new();
        for i in 0..num_pairs {
            let _ = builder
                .add_pair()
                .name(&format!("read{i:05}"))
                .start1((num_pairs - i) * 100 + 1)
                .start2((num_pairs - i) * 100 + 51)
                .build();
        }

        let workdir = tempfile::tempdir().expect("workdir");
        let input = workdir.path().join("input.bam");
        builder.write_bam(&input).expect("write bam");

        let run = |sorter: RawExternalSorter, out: &std::path::Path| {
            let stats = sorter.sort(&input, out).expect("sort should succeed");
            (stats.total_records, std::fs::read(out).expect("read output"))
        };

        let baseline_out = workdir.path().join("baseline.bam");
        let (baseline_records, baseline_bytes) = run(
            RawExternalSorter::new(SortOrder::Coordinate)
                .memory_limit(memory_limit)
                .threads(2)
                .write_index(write_index)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0),
            &baseline_out,
        );

        let split_out = workdir.path().join("split.bam");
        let (split_records, split_bytes) = run(
            RawExternalSorter::new(SortOrder::Coordinate)
                .memory_limit(memory_limit)
                .threads(2)
                .sort_threads(sort_threads)
                .merge_threads(merge_threads)
                .write_index(write_index)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0),
            &split_out,
        );

        assert_eq!(
            split_records, baseline_records,
            "per-phase threads must not change the record count",
        );
        assert_eq!(
            split_bytes, baseline_bytes,
            "per-phase threads ({sort_threads}/{merge_threads}, write_index={write_index}) \
             changed the output bytes; this is a scheduling knob and must be output-identical",
        );
    }

    /// The per-phase counts fall back to `--threads` independently, so an unset
    /// override must not silently pin a phase to 1.
    #[rstest::rstest]
    #[case::neither_set(None, None, 6, 6)]
    #[case::sort_only(Some(2), None, 2, 6)]
    #[case::merge_only(None, Some(3), 6, 3)]
    #[case::both_set(Some(2), Some(3), 2, 3)]
    #[case::zero_clamps_to_one(Some(0), Some(0), 1, 1)]
    fn test_phase_thread_defaults(
        #[case] sort_threads: Option<usize>,
        #[case] merge_threads: Option<usize>,
        #[case] expected_phase1: usize,
        #[case] expected_phase2: usize,
    ) {
        let mut sorter = RawExternalSorter::new(SortOrder::Coordinate).threads(6);
        if let Some(n) = sort_threads {
            sorter = sorter.sort_threads(n);
        }
        if let Some(n) = merge_threads {
            sorter = sorter.merge_threads(n);
        }
        assert_eq!(sorter.phase1_threads(), expected_phase1);
        assert_eq!(sorter.phase2_threads(), expected_phase2);
    }

    // ========================================================================
    // read_exact_or_eof tests
    // ========================================================================

    /// A `Read` impl that returns at most one byte per `read()` call. This is
    /// legal for `Read` (which may return any `0 < n <= buf.len()`) and is what
    /// motivates using `read_exact_or_eof` over a single `read()` call for
    /// fixed-size magic-byte detection.
    struct OneBytePerRead<'a> {
        bytes: &'a [u8],
        pos: usize,
    }
    impl std::io::Read for OneBytePerRead<'_> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if buf.is_empty() || self.pos >= self.bytes.len() {
                return Ok(0);
            }
            buf[0] = self.bytes[self.pos];
            self.pos += 1;
            Ok(1)
        }
    }

    #[test]
    fn test_read_exact_or_eof_clean_eof() {
        // Empty stream -> Ok(false), buffer untouched.
        let mut src: &[u8] = &[];
        let mut buf = [0u8; 4];
        let filled = read_exact_or_eof(&mut src, &mut buf).expect("clean eof");
        assert!(!filled, "clean EOF should return Ok(false)");
        assert_eq!(buf, [0u8; 4]);
    }

    #[test]
    fn test_read_exact_or_eof_full_read() {
        // Source has exactly buf.len() bytes -> Ok(true), buffer filled.
        let mut src: &[u8] = b"ZSP1";
        let mut buf = [0u8; 4];
        let filled = read_exact_or_eof(&mut src, &mut buf).expect("full read");
        assert!(filled, "full read should return Ok(true)");
        assert_eq!(&buf, b"ZSP1");
    }

    #[test]
    fn test_read_exact_or_eof_short_reads_fill_buffer() {
        // Regression: a `Read` impl that returns 1 byte per call must still
        // fill the 4-byte buffer. A naive single `read()` would only get 1
        // byte and `SpillCodec::from_magic` would return None, silently
        // routing a real ZSP1 spill through the uncompressed fallback path.
        let mut src = OneBytePerRead { bytes: b"ZSP1", pos: 0 };
        let mut buf = [0u8; 4];
        let filled = read_exact_or_eof(&mut src, &mut buf).expect("short-read source");
        assert!(filled, "partial reads should still fill the buffer");
        assert_eq!(&buf, b"ZSP1");
        // And `from_magic` correctly identifies the codec from the filled
        // buffer — proving the original bug is fixed end-to-end.
        assert_eq!(
            crate::codec::SpillCodec::from_magic(&buf),
            Some(crate::codec::SpillCodec::Zstd)
        );
    }

    #[test]
    fn test_read_exact_or_eof_truncated_is_error() {
        // Source has fewer bytes than requested but is non-empty:
        // must error with UnexpectedEof (not silently report short fill).
        let mut src: &[u8] = b"ZSP"; // only 3 of 4 bytes
        let mut buf = [0u8; 4];
        let err = read_exact_or_eof(&mut src, &mut buf).expect_err("truncated read should error");
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
    }

    // ========================================================================
    // verify_dropped_lanes / dropped_lane_error tests
    // ========================================================================

    fn full_key(cb: u64, tertiary: u64) -> TemplateKey {
        TemplateKey { primary: 0, secondary: 0, cb_hash: cb, tertiary, name_hash_upper: 0 }
    }

    const LITE: TemplateKeyVariant = TemplateKeyVariant { cb: false, tertiary: false };

    #[test]
    fn verify_detects_cb_appearing() {
        let v = verify_dropped_lanes(&full_key(0, 0), &full_key(123, 0), LITE);
        assert_eq!(v, Some(DroppedLaneViolation::Cb));
        assert_eq!(v.unwrap().key_types_token(), "cb");
    }

    #[test]
    fn verify_detects_mi_change() {
        let first = full_key(0, 0);
        let cur = full_key(0, 0b10); // mi_value bit set, library still 0
        assert_eq!(verify_dropped_lanes(&first, &cur, LITE), Some(DroppedLaneViolation::Mi));
    }

    #[test]
    fn verify_detects_library_change() {
        let first = full_key(0, 0);
        let cur = full_key(0, 1u64 << 48); // library ordinal 1
        assert_eq!(verify_dropped_lanes(&first, &cur, LITE), Some(DroppedLaneViolation::Library));
    }

    #[test]
    fn verify_passes_when_dropped_lanes_constant() {
        let variant = TemplateKeyVariant { cb: true, tertiary: false };
        assert_eq!(verify_dropped_lanes(&full_key(1, 5), &full_key(2, 5), variant), None);
    }

    #[test]
    fn verify_error_message_names_field_and_token() {
        let e = dropped_lane_error("read42", DroppedLaneViolation::Library);
        let msg = e.to_string();
        assert!(
            msg.contains("read42")
                && msg.contains("library")
                && msg.contains("--key-types library"),
            "{msg}"
        );
    }

    // ========================================================================
    // select_template_variant tests
    // ========================================================================

    fn key(cb: u64, tertiary: u64) -> TemplateKey {
        TemplateKey { primary: 1, secondary: 2, cb_hash: cb, tertiary, name_hash_upper: 3 }
    }

    #[test]
    fn auto_lite_when_no_cb_no_tertiary() {
        let v = select_template_variant(Some(&key(0, 0)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: false });
        assert_eq!(v.lanes(), 3);
    }

    #[test]
    fn auto_cb_only_when_first_record_has_cb() {
        let v = select_template_variant(Some(&key(42, 0)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: true, tertiary: false });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_tertiary_only_when_first_record_has_mi_or_library() {
        // 0xABCD is in the low-48 bits => MI present => keep tertiary, even though
        // the header library does not vary.
        let v = select_template_variant(Some(&key(0, 0xABCD)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_full_when_both_present() {
        let v = select_template_variant(Some(&key(7, 9)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: true, tertiary: true });
        assert_eq!(v.lanes(), 5);
    }

    // tertiary high bits only = a single library ordinal, no MI. With a non-varying
    // header, the constant library lane is dropped -> lite.
    #[test]
    fn auto_drops_constant_single_library_lane() {
        let k = key(0, 1u64 << 48); // library ordinal 1, no MI
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: false });
        assert_eq!(v.lanes(), 3);
    }

    #[test]
    fn auto_keeps_tertiary_when_header_library_varies() {
        let k = key(0, 1u64 << 48); // library ordinal 1, no MI
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, true); // >1 library
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_keeps_tertiary_for_mi_even_with_constant_library() {
        // MI present (low bits) -> keep tertiary even though header library is constant.
        let k = key(0, (1u64 << 48) | 0b10); // library ordinal 1 + MI value bit
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
    }

    #[test]
    fn none_forces_lite() {
        assert_eq!(
            select_template_variant(Some(&key(99, 99)), KeyTypesSpec::None, true).lanes(),
            3
        );
    }

    #[test]
    fn full_forces_all_lanes() {
        assert_eq!(select_template_variant(Some(&key(0, 0)), KeyTypesSpec::Full, false).lanes(), 5);
    }

    #[test]
    fn explicit_spec_passthrough() {
        assert_eq!(
            select_template_variant(
                None,
                KeyTypesSpec::Explicit { cb: true, tertiary: false },
                false
            )
            .lanes(),
            4
        );
    }

    #[test]
    fn empty_input_defaults_to_lite() {
        assert_eq!(select_template_variant(None, KeyTypesSpec::Auto, false).lanes(), 3);
    }

    // ========================================================================
    // Capacity estimator consistency tests (T8)
    // ========================================================================

    #[test]
    fn estimator_bytes_per_record_matches_ref_size() {
        use crate::inline::TemplateKey32;
        for (actual, expected) in [
            (std::mem::size_of::<TemplateRecordRef<TemplateKey24>>(), 40usize),
            (std::mem::size_of::<TemplateRecordRef<TemplateKey32>>(), 48),
            (std::mem::size_of::<TemplateRecordRef<TemplateKey40>>(), 56),
        ] {
            assert_eq!(actual, expected, "ref size must equal key + 16 B offset/len/pad");
        }
        let bpr_lite = (TEMPLATE_HEADER_SIZE + EST_BAM_BYTES_PER_TEMPLATE_RECORD)
            + std::mem::size_of::<TemplateRecordRef<TemplateKey24>>();
        assert_eq!(bpr_lite, 8 + 250 + 40);
    }

    // ========================================================================
    // T10: narrow-key byte-identity correctness proof
    // ========================================================================
    //
    // The whole point of `--key-types`: a narrow-key template-coordinate sort
    // must produce a record stream byte-identical to the full-key (`--key-types
    // full`) baseline on the same input. These tests build one unsorted BAM per
    // mode (with genuinely varying tags so the narrowed lane actually matters),
    // sort it under both the narrow spec and the full spec, and compare the
    // per-record BAM bodies in file order. They are non-vacuous: a sanity helper
    // asserts the constructed input produces the nonzero cb_hash / tertiary that
    // the mode claims, and the bonus `*_lite_sort_hard_errors` assertions prove a
    // mis-narrowed (lite) sort of the same input would refuse to run rather than
    // silently mis-sort.

    /// Number of distinct-position pairs per mode input — large enough that the
    /// spill-forced variant crosses a tiny memory limit and exercises the Phase-2
    /// k-way merge. Tied clusters (see `CLUSTER_SIZE` / `CLUSTER_COUNT`) are added
    /// on top of these.
    const T10_DISTINCT_PAIRS: usize = 300;

    /// Number of records in each coordinate-tied cluster. Within a cluster every
    /// record shares the same primary (tid, start1) AND secondary (start2, strands)
    /// and differs ONLY in the kept optional lane (`cb_hash` or tertiary), so that
    /// lane is the load-bearing tiebreaker above `name_hash_upper`. Made large so
    /// the chance that `name_hash` order coincidentally equals kept-lane order for
    /// every cluster (which would re-vacuate the test) is negligible.
    const CLUSTER_SIZE: usize = 8;

    /// Number of coordinate-tied clusters, each at a distinct shared coordinate.
    /// Several clusters further drive the coincidence probability to zero.
    const CLUSTER_COUNT: usize = 4;

    /// Total pairs (distinct + clustered) each mode input contains.
    const T10_PAIRS: usize = T10_DISTINCT_PAIRS + CLUSTER_SIZE * CLUSTER_COUNT;

    /// Collect each record's raw BAM bytes (excluding the header) in file order.
    ///
    /// Compares at the record-stream level rather than the file level: BGZF block
    /// batching is nondeterministic, so two semantically identical BAMs can differ
    /// byte-for-byte on disk while their decoded record bodies are identical.
    fn collect_record_bytes(path: &Path) -> Vec<Vec<u8>> {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, _) = create_raw_bam_reader(path, 1).expect("reader");
        RawReadAheadReader::new(reader).map(|rec| rec.as_ref().to_vec()).collect()
    }

    /// Sort `input` with the given key-types spec and return the record-byte stream.
    ///
    /// `cell_tag(SamTag::CB)` is REQUIRED: `RawExternalSorter::new` defaults the
    /// cell tag to `None`, so without this the CB lane is never populated and any
    /// single-cell assertion would pass vacuously (`cb_hash` stays 0). This mirrors
    /// the CLI's `parse_cell_tag`.
    fn sort_and_collect(input: &Path, dir: &Path, tag: &str, spec: KeyTypesSpec) -> Vec<Vec<u8>> {
        let out = dir.join(format!("{tag}.bam"));
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(spec)
            .sort(input, &out)
            .expect("sort");
        collect_record_bytes(&out)
    }

    /// Insert a read group carrying an `LB:` (library) field into a builder header.
    ///
    /// Used by the multi-library inputs to realize two distinct library ordinals
    /// (sorted alphabetically, then 1-based) so the tertiary lane genuinely varies.
    fn add_rg_with_library(header: &mut Header, rg_id: &str, library: &str) {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from(library))
            .build()
            .expect("valid read group");
        header.read_groups_mut().insert(BString::from(rg_id), rg);
    }

    /// Extract the full-width template key for the first record of a BAM, using a
    /// fresh `LibraryLookup` from its header. Lets a builder sanity-check that the
    /// tags it wrote actually produce the nonzero `cb_hash` / tertiary it claims.
    fn first_record_full_key(path: &Path) -> TemplateKey40 {
        use crate::read_ahead::RawReadAheadReader;
        let (reader, header) = create_raw_bam_reader(path, 1).expect("reader");
        let lib = LibraryLookup::from_header(&header);
        let first = RawReadAheadReader::new(reader).next().expect("at least one record");
        extract_template_key_inline(first.as_ref(), &lib, Some(SamTag::CB), &cb_hasher())
    }

    /// Scramble a cluster member's slot index into a kept-lane value index so the
    /// kept-lane order is a NON-identity permutation of the name order within the
    /// cluster. `7` is coprime to `CLUSTER_SIZE` (8), so this is a bijection that
    /// fixes no element to its name position (the `+3` shift removes the only fixed
    /// point). This anti-correlates the kept lane with `name_hash` so a wrong
    /// narrowing that falls back to `name_hash` cannot coincidentally reproduce the
    /// correct (kept-lane) order.
    fn scrambled_lane_index(member: usize) -> usize {
        (member * 7 + 3) % CLUSTER_SIZE
    }

    /// Shared start1 for cluster `c`. Far above any distinct-position record
    /// (`i * 200` for `i < T10_DISTINCT_PAIRS`, i.e. below ~`60_000`) so clusters
    /// never collide with the general path and each cluster sits at its own
    /// coordinate. All members of a cluster share this exact coordinate.
    fn cluster_start1(c: usize) -> usize {
        10_000_000 + c * 10_000
    }

    /// Build an unsorted BAM with `T10_DISTINCT_PAIRS` uniquely-positioned pairs
    /// PLUS `CLUSTER_COUNT` coordinate-tied clusters of `CLUSTER_SIZE` pairs each.
    ///
    /// `tagger` decorates a uniquely-positioned pair (CB/MI/RG). `cluster_tagger`
    /// decorates a cluster member, receiving the scrambled lane index so the caller
    /// can set the kept optional lane (CB / MI / RG) to a value whose order is
    /// anti-correlated with the member's name. The uniquely-positioned records
    /// exercise the general sort/spill path; the clusters make the kept lane the
    /// load-bearing tiebreaker (shared primary+secondary, differ only in kept lane).
    fn build_mode_bam(
        dir: &Path,
        name: &str,
        mut header: Header,
        tagger: impl Fn(PairBuilder<'_>, usize) -> PairBuilder<'_>,
        cluster_tagger: impl Fn(PairBuilder<'_>, usize) -> PairBuilder<'_>,
    ) -> PathBuf {
        let mut builder = SamBuilder::new();
        builder.header = std::mem::take(&mut header);

        // Coordinate-tied clusters first: each member shares start1/start2/strands
        // with its cluster siblings and differs only in the kept optional lane.
        for c in 0..CLUSTER_COUNT {
            let start1 = cluster_start1(c);
            let start2 = start1 + 100;
            for member in 0..CLUSTER_SIZE {
                let pair = builder
                    .add_pair()
                    .name(&format!("read_c{c}_{member:02}"))
                    .start1(start1)
                    .start2(start2);
                let _ = cluster_tagger(pair, scrambled_lane_index(member)).build();
            }
        }

        // Uniquely-positioned pairs: interleave positions (even i ascend, odd i
        // descend) so file order is far from template-coordinate order.
        for i in 0..T10_DISTINCT_PAIRS {
            let pos = if i % 2 == 0 { i * 200 } else { (T10_DISTINCT_PAIRS - i) * 200 };
            let pair =
                builder.add_pair().name(&format!("read{i:05}")).start1(pos + 1).start2(pos + 101);
            let _ = tagger(pair, i).build();
        }

        let path = dir.join(format!("{name}.bam"));
        builder.write_bam(&path).expect("write bam");
        path
    }

    /// Distinct CB barcodes, one per cluster slot, so a tied cluster's members get
    /// distinct `cb_hash` values (the load-bearing single-cell tiebreaker).
    const CLUSTER_BARCODES: [&str; CLUSTER_SIZE] = [
        "AAAACCCC", "GGGGTTTT", "ACGTACGT", "TTTTGGGG", "CCCCAAAA", "TTTTACGT", "GACTGACT",
        "TGCATGCA",
    ];

    /// bulk (lite): no CB, no MI, single default read group. No load-bearing opt
    /// lane, so the clusters carry no distinguishing tag — bulk byte-identity is
    /// inherently about the general path, and the cluster members tie down to
    /// `name_hash` in both baseline and narrow (lite drops both lanes anyway).
    fn build_bulk_bam(dir: &Path) -> PathBuf {
        build_mode_bam(
            dir,
            "bulk",
            SamBuilder::new().header.clone(),
            |pair, _| pair,
            |pair, _| pair,
        )
    }

    /// single-cell: every pair carries `CB:Z:<barcode>`; no MI. Cluster members get
    /// DISTINCT barcodes (scrambled vs name) so `cb_hash` is the load-bearing lane.
    fn build_single_cell_bam(dir: &Path) -> PathBuf {
        const BARCODES: [&str; 4] = ["AAAACCCC", "GGGGTTTT", "ACGTACGT", "TTTTGGGG"];
        build_mode_bam(
            dir,
            "single_cell",
            SamBuilder::new().header.clone(),
            |pair, i| pair.attr("CB", BARCODES[i % BARCODES.len()]),
            |pair, lane| pair.attr("CB", CLUSTER_BARCODES[lane]),
        )
    }

    /// post-group: every pair carries `MI:Z:<id>`; no CB; single library. Cluster
    /// members get DISTINCT MI ids (scrambled vs name), so the tertiary lane
    /// (MI bits) is the load-bearing tiebreaker. Ids are `1..=CLUSTER_SIZE` so
    /// tertiary != 0 for every member.
    fn build_post_group_bam(dir: &Path) -> PathBuf {
        build_mode_bam(
            dir,
            "post_group",
            SamBuilder::new().header.clone(),
            |pair, i| pair.attr("MI", format!("{}", (i % 5) + 1)),
            |pair, lane| pair.attr("MI", format!("{}", lane + 1)),
        )
    }

    /// multi-lib: header with TWO `@RG` lines whose `LB:` differ; pairs split across
    /// the two RGs; no CB, no MI. The two libraries realize ordinals 1 and 2. Within
    /// a cluster, members alternate read groups by the SCRAMBLED lane index so the
    /// tertiary library bits (high-16) vary load-bearingly across the cluster and
    /// the library order is anti-correlated with the name order.
    fn build_multi_lib_bam(dir: &Path) -> PathBuf {
        let mut header = SamBuilder::new().header.clone();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        add_rg_with_library(&mut header, "rgB", "LibBeta");
        build_mode_bam(
            dir,
            "multi_lib",
            header,
            |pair, i| pair.attr("RG", if i % 2 == 0 { "rgA" } else { "rgB" }),
            |pair, lane| pair.attr("RG", if lane % 2 == 0 { "rgA" } else { "rgB" }),
        )
    }

    /// full: CB + MI + multi-library together — all optional lanes vary. Cluster
    /// members get distinct CB and MI (scrambled vs name) plus alternating RG.
    /// Baseline (full) == narrow (full) so this mode is inherently trivial, but the
    /// clusters keep it meaningful (the kept lanes still vary load-bearingly).
    fn build_full_bam(dir: &Path) -> PathBuf {
        const BARCODES: [&str; 4] = ["AAAACCCC", "GGGGTTTT", "ACGTACGT", "TTTTGGGG"];
        let mut header = SamBuilder::new().header.clone();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        add_rg_with_library(&mut header, "rgB", "LibBeta");
        build_mode_bam(
            dir,
            "full",
            header,
            |pair, i| {
                pair.attr("CB", BARCODES[i % BARCODES.len()])
                    .attr("MI", format!("{}", i % 5))
                    .attr("RG", if i % 2 == 0 { "rgA" } else { "rgB" })
            },
            |pair, lane| {
                pair.attr("CB", CLUSTER_BARCODES[lane])
                    .attr("MI", format!("{}", lane + 1))
                    .attr("RG", if lane % 2 == 0 { "rgA" } else { "rgB" })
            },
        )
    }

    /// For each mode, the narrow-key sort must produce a record stream IDENTICAL to
    /// the `--key-types full` baseline, and the `Auto`-detected variant must match
    /// it too. The narrow spec is the minimal one for the mode; an over-narrow
    /// (lite) sort of single-cell / post-group / multi-lib inputs would instead
    /// hard-error (see `*_lite_sort_hard_errors`), proving these are not vacuous.
    #[rstest::rstest]
    #[case::bulk(build_bulk_bam as fn(&Path) -> PathBuf, KeyTypesSpec::None)]
    #[case::single_cell(
        build_single_cell_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: true, tertiary: false }
    )]
    #[case::post_group(
        build_post_group_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::multi_lib(
        build_multi_lib_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::full(build_full_bam as fn(&Path) -> PathBuf, KeyTypesSpec::Full)]
    fn narrow_key_sort_byte_identical_to_full(
        #[case] build: fn(&Path) -> PathBuf,
        #[case] narrow: KeyTypesSpec,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let input = build(dir.path());

        let baseline = sort_and_collect(&input, dir.path(), "full_baseline", KeyTypesSpec::Full);
        assert_eq!(baseline.len(), T10_PAIRS * 2, "baseline must retain every record");

        let narrowed = sort_and_collect(&input, dir.path(), "narrow", narrow);
        assert_eq!(narrowed, baseline, "narrow-key record stream must equal full-key stream");

        let auto = sort_and_collect(&input, dir.path(), "auto", KeyTypesSpec::Auto);
        assert_eq!(auto, baseline, "auto-detected variant must equal full-key stream");
    }

    /// single-library: one `@RG` carrying an `LB:`, every pair assigned to it; no CB,
    /// no MI. The lone library realizes ordinal 1, so the first record's tertiary is
    /// nonzero (high bits) yet the library lane is CONSTANT across every record — the
    /// WES/WGS scenario. `Auto` must therefore drop the tertiary lane to the 24-byte
    /// lite key while staying byte-identical to the full-key baseline.
    fn build_single_library_bam(dir: &Path) -> PathBuf {
        // Attach a library to the builder's single default read group ("A") so the
        // header realizes exactly ONE distinct library ordinal (1). Reads default to
        // RG "A", so every record maps to that lone library — the tertiary high bits
        // are a constant nonzero ordinal, with no MI and no CB.
        let mut header = SamBuilder::new().header.clone();
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibSolo"))
            .build()
            .expect("valid read group");
        // The builder defaults to a single read group "A"; replace it in place with a
        // library-bearing one so there is exactly one read group (one ordinal).
        header.read_groups_mut().clear();
        header.read_groups_mut().insert(BString::from("A"), rg);
        build_mode_bam(dir, "single_library", header, |pair, _| pair, |pair, _| pair)
    }

    /// WES/WGS regression: a single-library input (one constant library ordinal, no
    /// CB, no MI) must auto-narrow to the 3-lane lite key AND sort byte-identically to
    /// the full-key baseline. Proves the constant library lane is provisionally
    /// droppable under `Auto` and that dropping it does not change the record stream.
    #[test]
    fn auto_drops_constant_library_lane_byte_identical() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_single_library_bam(dir.path());

        // Production provisioning logic, reproduced: the first record carries a
        // nonzero (library) tertiary, but the header realizes a single ordinal, so
        // Auto must pick the 3-lane lite variant.
        let first = first_record_full_key(&input);
        let (_, header) = create_raw_bam_reader(&input, 1).expect("reader");
        let lib_lookup = LibraryLookup::from_header(&header);
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        assert!(!header_library_varies, "single-library header must not vary");
        assert_ne!(first.tertiary, 0, "first record carries a (constant) library ordinal");

        let variant =
            select_template_variant(Some(&first), KeyTypesSpec::Auto, header_library_varies);
        assert_eq!(variant.lanes(), 3, "auto must narrow a single-library input to the lite key");

        // Byte-identity: auto (lite) == full-key baseline.
        let baseline = sort_and_collect(&input, dir.path(), "full_baseline", KeyTypesSpec::Full);
        assert_eq!(baseline.len(), T10_PAIRS * 2, "baseline must retain every record");
        let auto = sort_and_collect(&input, dir.path(), "auto", KeyTypesSpec::Auto);
        assert_eq!(auto, baseline, "auto (lite) record stream must equal full-key stream");

        // Cross-check: an explicit lite (`--key-types none`) sort is also identical,
        // proving the lite key is genuinely valid for this input.
        let lite = sort_and_collect(&input, dir.path(), "lite", KeyTypesSpec::None);
        assert_eq!(lite, baseline, "explicit lite record stream must equal full-key stream");
    }

    /// Spill-forced byte-identity: a tiny memory limit forces Phase 1 to spill and
    /// Phase 2 to k-way merge `PooledChunkWriter::<K>` chunks. This is the path
    /// where the serialized narrow-key width must match the Phase-1 width; a
    /// mismatch would corrupt the merge. Asserted on the bulk (lite) input.
    #[test]
    fn narrow_key_sort_byte_identical_when_spilling() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_bulk_bam(dir.path());

        let sort_spilling = |tag: &str, spec: KeyTypesSpec| -> (Vec<Vec<u8>>, usize) {
            let out = dir.path().join(format!("{tag}.bam"));
            let stats = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(32 * 1024)
                .max_temp_files(0)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0)
                .cell_tag(SamTag::CB)
                .key_types(spec)
                .sort(&input, &out)
                .expect("sort");
            (collect_record_bytes(&out), stats.runs_written)
        };

        let (baseline, base_chunks) = sort_spilling("spill_full", KeyTypesSpec::Full);
        let (narrowed, narrow_chunks) = sort_spilling("spill_narrow", KeyTypesSpec::None);

        assert!(base_chunks >= 2, "baseline must spill multiple chunks, got {base_chunks}");
        assert!(narrow_chunks >= 2, "narrow must spill multiple chunks, got {narrow_chunks}");
        assert_eq!(narrowed, baseline, "spilled narrow stream must equal spilled full stream");
    }

    /// Spill-forced byte-identity for ALL key widths: a tiny memory limit forces Phase 1
    /// to spill and Phase 2 to k-way merge `PooledChunkWriter::<K>` chunks across all four
    /// key widths (24-byte lite, 32-byte `CbKey32`, 32-byte `TertKey32`, 40-byte `TemplateKey40`).
    ///
    /// The existing `narrow_key_sort_byte_identical_when_spilling` test only covers the
    /// 24-byte (bulk/lite) path. A width/serialization regression in the 32- or 40-byte
    /// Phase-2 merge path would not be caught by that test alone. This test closes the gap
    /// by asserting both that each width's narrow-key spilled stream is byte-identical to
    /// the Full baseline spilled stream, AND that both runs actually spilled (>=2 chunks),
    /// so the Phase-2 merge path at each key width is genuinely exercised.
    #[rstest::rstest]
    #[case::bulk(build_bulk_bam as fn(&Path) -> PathBuf, KeyTypesSpec::None)]
    #[case::single_cell(
        build_single_cell_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: true, tertiary: false }
    )]
    #[case::post_group(
        build_post_group_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::multi_lib(
        build_multi_lib_bam as fn(&Path) -> PathBuf,
        KeyTypesSpec::Explicit { cb: false, tertiary: true }
    )]
    #[case::full(build_full_bam as fn(&Path) -> PathBuf, KeyTypesSpec::Full)]
    fn narrow_key_spill_merge_byte_identical_all_widths(
        #[case] build: fn(&Path) -> PathBuf,
        #[case] narrow: KeyTypesSpec,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let input = build(dir.path());

        let sort_spilling = |tag: &str, spec: KeyTypesSpec| -> (Vec<Vec<u8>>, usize) {
            let out = dir.path().join(format!("{tag}.bam"));
            let stats = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(32 * 1024)
                .max_temp_files(0)
                .spill_codec(crate::codec::SpillCodec::Bgzf)
                .temp_compression(0)
                .output_compression(0)
                .cell_tag(SamTag::CB)
                .key_types(spec)
                .sort(&input, &out)
                .expect("sort");
            (collect_record_bytes(&out), stats.runs_written)
        };

        let (baseline, base_chunks) = sort_spilling("spill_full_baseline", KeyTypesSpec::Full);
        let (narrowed, narrow_chunks) = sort_spilling("spill_narrow", narrow);

        assert!(
            base_chunks >= 2,
            "Full-key baseline must spill multiple chunks (Phase-2 merge must run), \
             got {base_chunks}"
        );
        assert!(
            narrow_chunks >= 2,
            "Narrow-key sort must spill multiple chunks (Phase-2 merge must run at this \
             key width), got {narrow_chunks}"
        );
        assert_eq!(
            narrowed, baseline,
            "Narrow-key spilled record stream must be byte-identical to Full-key spilled \
             record stream — a mismatch indicates a width/serialization regression in the \
             Phase-2 merge path"
        );
    }

    /// Sanity (anti-vacuity): each mode's input must produce the nonzero `cb_hash` /
    /// tertiary it claims on its first record, otherwise the byte-identity tests
    /// would compare two identically-narrowed streams and prove nothing.
    #[test]
    fn mode_inputs_realize_claimed_lanes() {
        let dir = tempfile::tempdir().unwrap();

        let bulk = first_record_full_key(&build_bulk_bam(dir.path()));
        assert_eq!(bulk.cb_hash, 0, "bulk must have no CB");
        assert_eq!(bulk.tertiary, 0, "bulk must have no MI/library");

        let cell = first_record_full_key(&build_single_cell_bam(dir.path()));
        assert_ne!(cell.cb_hash, 0, "single-cell first record must have nonzero cb_hash");
        assert_eq!(cell.tertiary, 0, "single-cell must have no MI/library");

        let pg = first_record_full_key(&build_post_group_bam(dir.path()));
        assert_eq!(pg.cb_hash, 0, "post-group must have no CB");
        assert_ne!(pg.tertiary, 0, "post-group first record must have nonzero MI tertiary");
        assert_eq!(pg.tertiary >> 48, 0, "post-group must have no library (ordinal 0)");

        let multi = first_record_full_key(&build_multi_lib_bam(dir.path()));
        assert_eq!(multi.cb_hash, 0, "multi-lib must have no CB");
        assert_ne!(multi.tertiary >> 48, 0, "multi-lib first record must have nonzero library");

        let full = first_record_full_key(&build_full_bam(dir.path()));
        assert_ne!(full.cb_hash, 0, "full first record must have nonzero cb_hash");
        assert_ne!(full.tertiary, 0, "full first record must have nonzero tertiary");
        assert_ne!(full.tertiary >> 48, 0, "full first record must have nonzero library ordinal");
    }

    /// Guard the single-cell cluster against a silent `cb_hash` collision: every
    /// barcode in `CLUSTER_BARCODES` must hash to a DISTINCT `cb_hash` value under
    /// `cb_hasher()`, otherwise two cluster members would share the same sort key
    /// and the tiebreaker designed to make the narrow-key test non-vacuous would
    /// silently collapse.
    ///
    /// Uses `cb_hasher().hash_one(bc.as_bytes())` — exactly the expression in
    /// `extract_template_key_inline` (`cb_hasher.hash_one(cb_bytes)` where
    /// `cb_bytes` is the raw string value of the CB aux tag, i.e. ASCII barcode
    /// bytes without NUL terminator).
    #[test]
    fn cluster_barcodes_have_distinct_cb_hashes() {
        let hasher = cb_hasher();
        let mut hash_vals: Vec<u64> =
            CLUSTER_BARCODES.iter().map(|bc| hasher.hash_one(bc.as_bytes())).collect();
        let n = hash_vals.len();
        hash_vals.sort_unstable();
        hash_vals.dedup();
        assert_eq!(
            hash_vals.len(),
            n,
            "cluster barcodes must hash to distinct cb_hash values, else the single-cell \
             cluster is vacuous"
        );
    }

    /// A narrow (tertiary-only) sort output must re-verify correct at FULL width via
    /// `core_cmp`. This locks the consistency between the narrowed key used for
    /// sorting and the full key used by `fgumi sort --verify`: if narrowing dropped
    /// ordering information the merge depended on, the full-width re-check would
    /// find a violation.
    #[test]
    fn narrow_sort_output_passes_full_width_verify() {
        use std::cmp::Ordering;

        let dir = tempfile::tempdir().unwrap();
        let input = build_post_group_bam(dir.path());
        let out = dir.path().join("narrow.bam");
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::Explicit { cb: false, tertiary: true })
            .sort(&input, &out)
            .expect("sort");

        // Header only (for LibraryLookup); the auto reader is dropped immediately.
        let (_, header) = create_raw_bam_reader(&out, 1).expect("header");
        let lib = LibraryLookup::from_header(&header);
        let hasher = cb_hasher();

        // verify_sort_order wants RawBamRecordReader<File> — build it like execute_verify.
        let file = std::fs::File::open(&out).expect("open");
        let mut raw_reader = crate::reader::RawBamRecordReader::new(file).expect("reader");
        raw_reader.skip_header().expect("skip header");

        let (total, violations, first) = crate::verify::verify_sort_order(
            raw_reader,
            |bam| extract_template_key_inline(bam, &lib, Some(SamTag::CB), &hasher),
            |cur: &TemplateKey40, prev: &TemplateKey40| cur.core_cmp(prev) == Ordering::Less,
        )
        .expect("verify runs");

        assert_eq!(total, (T10_PAIRS * 2) as u64, "verify must see every record");
        assert_eq!(
            violations, 0,
            "narrow output must re-verify at full width (first={first:?}, total={total})"
        );
    }

    // ------------------------------------------------------------------------
    // T6 (deferred) hard-error integration tests
    //
    // Build a BAM whose FIRST pair lacks an optional field and a LATER pair
    // carries/changes it, then sort with the over-narrow spec that drops that
    // lane. The decode-time `verify_dropped_lanes` check must hard-error with a
    // message naming the `--key-types <token>` that re-includes the lane, rather
    // than silently mis-sorting. Each test forces the relevant lane to vary
    // across records while holding the dropped lanes constant; attribution
    // (cb / mi / library) follows `verify_dropped_lanes`.
    // ------------------------------------------------------------------------

    /// Build a two-pair BAM where the first pair omits a tag the second supplies.
    fn build_first_then_appears(
        dir: &Path,
        name: &str,
        header: Header,
        first: impl Fn(PairBuilder<'_>) -> PairBuilder<'_>,
        second: impl Fn(PairBuilder<'_>) -> PairBuilder<'_>,
    ) -> PathBuf {
        let mut builder = SamBuilder::new();
        builder.header = header;
        let _ = first(builder.add_pair().name("read00000").start1(1).start2(101)).build();
        let _ = second(builder.add_pair().name("read00001").start1(201).start2(301)).build();
        let path = dir.join(format!("{name}.bam"));
        builder.write_bam(&path).expect("write bam");
        path
    }

    /// CB appears after the first record under `--key-types none` → hard error
    /// instructing `--key-types cb`.
    #[test]
    fn lite_sort_hard_errors_when_cb_appears() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_first_then_appears(
            dir.path(),
            "cb_appears",
            SamBuilder::new().header.clone(),
            |pair| pair,
            |pair| pair.attr("CB", "AAAACCCC"),
        );
        let out = dir.path().join("out.bam");
        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::None)
            .sort(&input, &out)
            .expect_err("CB appearing after a CB-free first record must hard-error under lite");
        let msg = err.to_string();
        assert!(msg.contains("--key-types cb"), "expected --key-types cb guidance, got: {msg}");
    }

    /// MI appears/changes after the first record under `--key-types none` → hard
    /// error instructing `--key-types mi`.
    #[test]
    fn lite_sort_hard_errors_when_mi_appears() {
        let dir = tempfile::tempdir().unwrap();
        let input = build_first_then_appears(
            dir.path(),
            "mi_appears",
            SamBuilder::new().header.clone(),
            |pair| pair,
            |pair| pair.attr("MI", "7"),
        );
        let out = dir.path().join("out.bam");
        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::None)
            .sort(&input, &out)
            .expect_err("MI appearing after an MI-free first record must hard-error under lite");
        let msg = err.to_string();
        assert!(msg.contains("--key-types mi"), "expected --key-types mi guidance, got: {msg}");
    }

    /// Library differs after the first record under `--key-types none` → hard error
    /// instructing `--key-types library`. The first pair carries a read group whose
    /// `LB:` realizes ordinal 1; the second pair carries NO RG tag, realizing
    /// ordinal 0 — exercising the "header undercounts realized libraries" case
    /// where a later record drops to an unseen library ordinal.
    #[test]
    fn lite_sort_hard_errors_when_library_differs() {
        let dir = tempfile::tempdir().unwrap();
        // Header has exactly one read group (with LB) → ordinal 1; the default
        // SamBuilder "A" group is replaced so the only library is LibAlpha.
        let mut header = Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(200_000_000).expect("nonzero"),
                ),
            )
            .build();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");

        // First pair: RG=rgA (ordinal 1). Second pair: RG points at an id absent
        // from the header → realized ordinal 0, differing from the first.
        let input = build_first_then_appears(
            dir.path(),
            "library_differs",
            header,
            |pair| pair.attr("RG", "rgA"),
            |pair| pair.attr("RG", "rgUnknown"),
        );
        let out = dir.path().join("out.bam");
        let err = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::None)
            .sort(&input, &out)
            .expect_err("a differing library ordinal must hard-error under lite");
        let msg = err.to_string();
        assert!(
            msg.contains("--key-types library"),
            "expected --key-types library guidance, got: {msg}"
        );
    }

    /// Regression (#375): a single-library header with MULTIPLE read groups (all
    /// sharing one `LB:`) must sort under auto `--key-types` WITHOUT a false
    /// dropped-lane "library" violation.
    ///
    /// Every read group realizes the SAME library ordinal, so
    /// `distinct_header_ordinals() == 1`, the tertiary library lane is genuinely
    /// constant, and auto correctly drops it. Decode-time verify must therefore
    /// see a constant library ordinal across records even though they carry
    /// DIFFERENT (but same-library) RG tags. Mirrors production multi-lane WGS
    /// BAMs (e.g. 12 `@RG` lines, one `LB`) that regressed to a spurious
    /// "carries a library value absent from the input's first record" error.
    #[test]
    fn auto_single_library_multi_readgroup_sorts_without_false_violation() {
        let dir = tempfile::tempdir().unwrap();
        // Build the header from scratch (no SamBuilder default read group, which
        // would add a second library) so the ONLY libraries are the ones we add.
        let mut header = Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(200_000_000).expect("nonzero"),
                ),
            )
            .build();
        // Two read groups, SAME library -> exactly one distinct library ordinal.
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        add_rg_with_library(&mut header, "rgB", "LibAlpha");
        assert_eq!(
            LibraryLookup::from_header(&header).distinct_header_ordinals(),
            1,
            "two read groups sharing one library must realize a single ordinal"
        );

        // Pairs split across the two read groups; both resolve to LibAlpha.
        let input = build_mode_bam(
            dir.path(),
            "single_lib_multi_rg",
            header,
            |pair, i| pair.attr("RG", if i % 2 == 0 { "rgA" } else { "rgB" }),
            |pair, lane| pair.attr("RG", if lane % 2 == 0 { "rgA" } else { "rgB" }),
        );
        let out = dir.path().join("out.bam");
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::Auto)
            .sort(&input, &out)
            .expect(
                "single-library input with multiple read groups must sort without a \
                 false dropped-lane library violation",
            );
    }

    /// Regression (#375): a single-library input that contains a COMPLETELY
    /// UNMAPPED read (both mates unmapped) must sort under auto `--key-types`
    /// WITHOUT a spurious dropped-lane "library" violation.
    ///
    /// `extract_template_key_inline` packs the read's library ordinal into the
    /// full key for mapped reads (and for unmapped-with-mapped-mate reads), but
    /// `TemplateKey::unmapped` zeroed the whole `tertiary` lane — so a
    /// fully-unmapped read carrying a valid RG realized library ordinal 0 while
    /// its mapped, same-library peers realized ordinal 1. With a single library
    /// the auto path drops the (provably-constant) library lane, and decode-time
    /// verify then saw 1 (the mapped first record) vs 0 (the unmapped read) and
    /// hard-errored with "carries a library value absent from the input's first
    /// record". Reproduces the 1kg WGS regression (unmapped pairs at the tail).
    #[test]
    fn auto_single_library_unmapped_read_sorts_without_false_violation() {
        let dir = tempfile::tempdir().unwrap();
        // Single read group with a library -> exactly one library ordinal, so the
        // auto path drops the tertiary library lane and relies on decode verify.
        let mut header = Header::builder()
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(200_000_000).expect("nonzero"),
                ),
            )
            .build();
        add_rg_with_library(&mut header, "rgA", "LibAlpha");
        assert_eq!(LibraryLookup::from_header(&header).distinct_header_ordinals(), 1);

        let mut builder = SamBuilder::new();
        builder.header = header;
        // A normal mapped pair (library ordinal 1) ...
        let _ =
            builder.add_pair().name("mapped").start1(1_000).start2(1_100).attr("RG", "rgA").build();
        // ... and a COMPLETELY unmapped pair carrying the same RG. Its realized
        // library ordinal is also 1, but the old unmapped key zeroed it.
        let _ =
            builder.add_pair().name("unmapped").unmapped1().unmapped2().attr("RG", "rgA").build();
        let input = dir.path().join("in.bam");
        builder.write_bam(&input).expect("write bam");

        let out = dir.path().join("out.bam");
        RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .output_compression(0)
            .cell_tag(SamTag::CB)
            .key_types(KeyTypesSpec::Auto)
            .sort(&input, &out)
            .expect(
                "a completely-unmapped read in a single-library input must not trigger a \
                 false dropped-lane library violation",
            );
    }
}
