#![deny(unsafe_code)]
//! Shared CLI types and helpers for fgumi commands.
//!
//! Holds the pieces every command needs in common: the [`Command`] trait, the [`FgumiError`]
//! type, host-capacity detection, memory-budget parsing and resolution, and the clap argument
//! groups shared across subcommands.

// ─────────────────────────────────────────────────────────────────────────────
// Command trait
// ─────────────────────────────────────────────────────────────────────────────

/// Trait implemented by all fgumi CLI commands.
///
/// Each command provides an `execute` method that runs the command's main logic.
/// The `command_line` parameter contains the full command invocation for @PG records.
///
/// Commands are dispatched by a hand-written `match` in the binary rather than by
/// `#[enum_dispatch]`: that macro pairs a trait with its enum through a registry local to a
/// single proc-macro invocation, so it cannot link a trait defined here to an enum defined
/// in another crate.
pub trait Command {
    #[allow(clippy::missing_errors_doc)]
    fn execute(&self, command_line: &str) -> anyhow::Result<()>;
}

// ─────────────────────────────────────────────────────────────────────────────
// Error types
// ─────────────────────────────────────────────────────────────────────────────

use thiserror::Error;

/// Result type alias for fgumi operations (preferred in standalone crates).
pub type FgumiResult<T> = std::result::Result<T, FgumiError>;

/// Unqualified result alias used by the umbrella crate's error/validation modules.
///
/// Both `FgumiResult<T>` and `Result<T>` are the same type; the two names exist
/// so callers that shadow `std::result::Result` with `use crate::errors::Result`
/// (umbrella convention) resolve to the same `FgumiError`-based alias.
pub type Result<T> = std::result::Result<T, FgumiError>;

/// Error type for fgumi operations
#[derive(Error, Debug)]
pub enum FgumiError {
    /// Invalid parameter value provided
    #[error("Invalid parameter '{parameter}': {reason}")]
    InvalidParameter {
        /// The parameter name
        parameter: String,
        /// Explanation of why it's invalid
        reason: String,
    },

    /// Invalid frequency threshold
    #[error("Invalid frequency threshold: {value} (must be between {min} and {max})")]
    InvalidFrequency {
        /// The invalid frequency value
        value: f64,
        /// Minimum valid value
        min: f64,
        /// Maximum valid value
        max: f64,
    },

    /// Invalid quality threshold
    #[error("Invalid quality threshold: {value} (must be between 0 and {max})")]
    InvalidQuality {
        /// The invalid quality value
        value: u8,
        /// Maximum valid value (usually 93 for SAM/BAM)
        max: u8,
    },

    /// File format error
    #[error("Invalid {file_type} file '{path}': {reason}")]
    InvalidFileFormat {
        /// Type of file (e.g., "BAM", "FASTQ")
        file_type: String,
        /// Path to the file
        path: String,
        /// Explanation of the problem
        reason: String,
    },

    /// Required reference sequence not found
    #[error("Reference sequence '{ref_name}' not found in header")]
    ReferenceNotFound {
        /// The reference sequence name
        ref_name: String,
    },

    /// Invalid memory size string
    #[error("Invalid memory size: {reason}")]
    InvalidMemorySize {
        /// Explanation of why the value is invalid
        reason: String,
    },
}

// ─────────────────────────────────────────────────────────────────────────────
// System detection
// ─────────────────────────────────────────────────────────────────────────────

/// Returns the effective total memory available to this process in bytes.
///
/// Checks cgroup memory limits for container environments, falling back to
/// physical RAM on bare-metal or macOS. Always returns `min(cgroup, physical)`.
///
/// Call this once and reuse the result: it constructs a `sysinfo::System`, reads
/// `/proc/meminfo`, and probes the cgroup v1/v2 limit files, none of which is free.
/// [`resolve_memory_budget`] deliberately calls it exactly once.
#[must_use]
pub fn detect_total_memory() -> usize {
    let mut system = sysinfo::System::new();
    system.refresh_memory();
    let physical = system.total_memory();
    let bytes = system.cgroup_limits().map_or(physical, |c| c.total_memory.min(physical));
    // Saturate at usize::MAX / 2 rather than usize::MAX on 32-bit platforms so
    // the downstream `budget > total` overflow check in `resolve_memory_budget`
    // can fire correctly (no value can exceed usize::MAX, so using it as the
    // fallback renders the check dead).
    usize::try_from(bytes).unwrap_or(usize::MAX / 2)
}

/// Returns the number of logical CPUs available to this process.
///
/// Honors cgroup CPU quotas (e.g. `--cpus` in Docker or Kubernetes resource limits), so a
/// container gets its quota rather than the host's core count. Returns 1 when the platform
/// cannot report parallelism.
#[must_use]
pub fn detect_cpu_count() -> usize {
    std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
}

// ─────────────────────────────────────────────────────────────────────────────
// Formatting helpers
// ─────────────────────────────────────────────────────────────────────────────

// Pure formatting with no CLI dependency, so it lives in the `fgumi-fmt` leaf crate that
// lower layers can reach too. Re-exported here for callers already depending on this crate.
pub use fgumi_fmt::{format_count, format_duration, format_rate};

// ─────────────────────────────────────────────────────────────────────────────
// Operation timer
// ─────────────────────────────────────────────────────────────────────────────

/// Operation timing and summary helper.
///
/// Tracks operation timing and provides formatted summary output.
pub struct OperationTimer {
    operation: String,
    start_time: std::time::Instant,
}

impl OperationTimer {
    /// Creates a new operation timer and logs the start.
    #[must_use]
    pub fn new(operation: &str) -> Self {
        log::info!("{operation} ...");
        Self { operation: operation.to_string(), start_time: std::time::Instant::now() }
    }

    /// Logs the completion with item count and rate.
    pub fn log_completion(&self, count: u64) {
        let duration = self.start_time.elapsed();
        log::info!(
            "{} completed: {} in {} ({})",
            self.operation,
            format_count(count),
            format_duration(duration),
            format_rate(count, duration)
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation helpers
// ─────────────────────────────────────────────────────────────────────────────

use bytesize::ByteSize;
use std::path::Path;

/// Validate that a file exists.
///
/// # Errors
///
/// Returns [`FgumiError::InvalidFileFormat`] if the file does not exist.
pub fn validate_file_exists<P: AsRef<Path>>(path: P, description: &str) -> FgumiResult<()> {
    let path_ref = path.as_ref();
    if !path_ref.exists() {
        return Err(FgumiError::InvalidFileFormat {
            file_type: description.to_string(),
            path: path_ref.display().to_string(),
            reason: "File does not exist".to_string(),
        });
    }
    Ok(())
}

/// Returns `true` if `s` contains a genuine decimal scientific-notation
/// mantissa-exponent (e.g. `1e5`, `1.5e10`, `2E+3`).
///
/// Only an `e`/`E` that is preceded by a digit or `.` **and** followed by a
/// digit or a sign counts. This deliberately does *not* match the exabyte units
/// `EB`/`EiB` (where `E` is followed by `B`/`i`), which `bytesize` accepts — the
/// earlier `contains('e') || contains('E')` heuristic wrongly rejected them.
fn looks_like_scientific_notation(s: &str) -> bool {
    let bytes = s.as_bytes();
    bytes.iter().enumerate().any(|(i, &b)| {
        if b != b'e' && b != b'E' {
            return false;
        }
        let prev_is_mantissa = i > 0 && (bytes[i - 1].is_ascii_digit() || bytes[i - 1] == b'.');
        let next_is_exponent = matches!(
            bytes.get(i + 1),
            Some(&c) if c.is_ascii_digit() || c == b'+' || c == b'-'
        );
        prev_is_mantissa && next_is_exponent
    })
}

/// Parses a memory size string into bytes.
///
/// Accepts both plain numbers (interpreted as MiB) and human-readable formats like:
/// - "2GB", "2G" -> 2 gigabytes (decimal: 2,000,000,000)
/// - "1.5GB" -> 1.5 gigabytes
/// - "1024MB", "1024M" -> 1024 megabytes (decimal)
/// - "512MiB" -> 512 mebibytes (binary: 536,870,912)
/// - "1EB", "1EiB" -> 1 exabyte / exbibyte (decimal / binary; the largest units
///   accepted — values above `u64::MAX` saturate rather than wrapping)
/// - "768" -> 768 MiB (plain numbers are interpreted as mebibytes)
///
/// # Errors
///
/// Returns [`FgumiError::InvalidMemorySize`] if the string cannot be parsed.
pub fn parse_memory_size(size_str: &str) -> FgumiResult<u64> {
    let trimmed = size_str.trim();
    if trimmed.is_empty() {
        return Err(FgumiError::InvalidMemorySize {
            reason: "Memory size cannot be empty".to_string(),
        });
    }

    if trimmed.starts_with('-') {
        return Err(FgumiError::InvalidMemorySize {
            reason: format!("Memory size cannot be negative: '{trimmed}'"),
        });
    }

    if let Ok(mb_value) = trimmed.parse::<u64>() {
        if mb_value == 0 {
            return Err(FgumiError::InvalidMemorySize {
                reason: "Memory size cannot be zero".to_string(),
            });
        }
        if mb_value > 1_000_000 {
            return Err(FgumiError::InvalidMemorySize {
                reason: format!(
                    "Plain number memory size too large: {} MiB. Use human-readable format like '{}GiB' instead.",
                    mb_value,
                    mb_value / 1024
                ),
            });
        }
        return mb_value.checked_mul(1024 * 1024).ok_or_else(|| FgumiError::InvalidMemorySize {
            reason: format!("Memory size calculation overflow for {mb_value} MiB"),
        });
    }

    if looks_like_scientific_notation(trimmed) {
        return Err(FgumiError::InvalidMemorySize {
            reason: format!(
                "Scientific notation not supported: '{trimmed}'. Use integer values or human-readable formats like '2GB'."
            ),
        });
    }

    if trimmed.contains('.') && trimmed.chars().all(|c| c.is_ascii_digit() || c == '.') {
        return Err(FgumiError::InvalidMemorySize {
            reason: format!(
                "Plain decimal numbers not supported: '{trimmed}'. Use an integer for MiB (e.g. '768') or a human-readable format (e.g. '1.5GB')."
            ),
        });
    }

    match trimmed.parse::<ByteSize>() {
        Ok(size) => {
            if size.0 == 0 {
                return Err(FgumiError::InvalidMemorySize {
                    reason: format!("Memory size cannot be zero: '{trimmed}'"),
                });
            }
            Ok(size.0)
        }
        Err(_) => Err(FgumiError::InvalidMemorySize {
            reason: format!(
                "Invalid memory size '{trimmed}'. Valid formats:\n\
                 - Plain numbers (interpreted as MiB): '768', '4096'\n\
                 - Human-readable (decimal): '2GB', '1024MB'\n\
                 - Human-readable (binary): '1GiB', '512MiB'"
            ),
        }),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Memory/compression options
// ─────────────────────────────────────────────────────────────────────────────

/// A memory limit, either auto-detected from the host or a fixed byte count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryLimit {
    /// Detect the (cgroup-aware) host memory and subtract the reserve.
    Auto,
    /// Use a fixed memory limit in bytes.
    Fixed(usize),
}

impl Default for MemoryLimit {
    fn default() -> Self {
        Self::Fixed(768 * 1024 * 1024)
    }
}

impl std::fmt::Display for MemoryLimit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => f.write_str("auto"),
            Self::Fixed(bytes) => format_binary_bytes(*bytes, f),
        }
    }
}

/// How much memory to reserve for other processes when [`MemoryLimit::Auto`] is used.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MemoryReserve {
    /// Automatic: `min(10 GiB, 50% of host memory)`.
    /// Matches the clap `default_value = "auto"` on `SortOptions::memory_reserve`.
    #[default]
    Auto,
    /// Reserve a fixed number of bytes.
    Fixed(usize),
}

impl std::fmt::Display for MemoryReserve {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Auto => f.write_str("auto"),
            Self::Fixed(bytes) => format_binary_bytes(*bytes, f),
        }
    }
}

/// Format a byte count in the largest binary unit that divides cleanly.
///
/// Emits the compact, unspaced form the CLI accepts as input (`768MiB`), so a default shown
/// in `--help` can be pasted straight back on the command line. `bytesize::ByteSize`'s own
/// `Display` also round-trips through [`parse_memory`], but renders as `768.0 MiB` — a space
/// and a redundant `.0`, which reads poorly as an argument default.
fn format_binary_bytes(bytes: usize, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    const K: usize = 1024;
    const M: usize = K * 1024;
    const G: usize = M * 1024;
    if bytes >= G && bytes.is_multiple_of(G) {
        write!(f, "{}GiB", bytes / G)
    } else if bytes >= M && bytes.is_multiple_of(M) {
        write!(f, "{}MiB", bytes / M)
    } else if bytes >= K && bytes.is_multiple_of(K) {
        write!(f, "{}KiB", bytes / K)
    } else {
        write!(f, "{bytes}B")
    }
}

/// The minimum per-thread memory budget (256 MiB).
pub const MIN_MEMORY_PER_THREAD: usize = 256 * 1024 * 1024;

/// Default auto-reserve cap: 10 GiB.
///
/// Computed with `saturating_mul` so it does not const-overflow `usize` on
/// 32-bit targets (where `usize == u32` cannot hold 10 GiB). There it saturates
/// to `usize::MAX`, so `resolve_reserve`'s `AUTO_RESERVE_CAP.min(total / 2)`
/// simply yields `total / 2` — consistent with `detect_total_memory`, which also
/// deliberately supports 32-bit by saturating.
const AUTO_RESERVE_CAP: usize = 10usize.saturating_mul(1024 * 1024 * 1024);

/// Parse a memory-limit string (e.g. "512M", "1G", "768", "auto").
///
/// # Errors
///
/// Returns an error string if parsing fails.
pub fn parse_memory(s: &str) -> std::result::Result<MemoryLimit, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryLimit::Auto);
    }
    Ok(MemoryLimit::Fixed(parse_memory_bytes(s, "Memory size")?))
}

/// Parse a memory-reserve string (e.g. "10G", "auto").
///
/// # Errors
///
/// Returns an error string if parsing fails.
pub fn parse_memory_reserve(s: &str) -> std::result::Result<MemoryReserve, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryReserve::Auto);
    }
    Ok(MemoryReserve::Fixed(parse_memory_bytes(s, "Memory reserve")?))
}

/// Resolve a [`MemoryReserve`] to a concrete byte count given total host memory.
#[must_use]
pub fn resolve_reserve(reserve: MemoryReserve, total_memory: usize) -> usize {
    match reserve {
        MemoryReserve::Fixed(bytes) => bytes,
        MemoryReserve::Auto => AUTO_RESERVE_CAP.min(total_memory / 2),
    }
}

/// Resolve a memory budget to a concrete byte count.
///
/// # Errors
///
/// Returns an error if `threads` is 0 or the multiplication overflows.
pub fn resolve_memory_budget(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    per_thread: bool,
) -> anyhow::Result<usize> {
    resolve_memory_budget_with_total(limit, reserve, threads, per_thread, detect_total_memory())
}

/// Pure resolver behind [`resolve_memory_budget`], with `total` injected for testability.
fn resolve_memory_budget_with_total(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    per_thread: bool,
    total: usize,
) -> anyhow::Result<usize> {
    if threads == 0 {
        anyhow::bail!("--threads must be at least 1");
    }

    let budget = match limit {
        MemoryLimit::Fixed(bytes) => {
            if per_thread {
                bytes
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("memory limit × {threads} threads overflowed"))?
            } else {
                bytes
            }
        }
        MemoryLimit::Auto => {
            let margin = resolve_reserve(reserve, total);
            let available = total.saturating_sub(margin);
            let target = if per_thread {
                (available / threads)
                    .max(MIN_MEMORY_PER_THREAD)
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("auto memory budget overflowed"))?
            } else {
                available.max(MIN_MEMORY_PER_THREAD)
            };
            let budget = target.min(available);
            if budget < target {
                log::warn!(
                    "Auto memory: capping budget to host-available {} (minimum viable target {} \
                     exceeds it after reserve {}); throughput may drop but the run stays within memory",
                    bytesize::ByteSize(budget as u64),
                    bytesize::ByteSize(target as u64),
                    bytesize::ByteSize(margin as u64),
                );
            }
            // Only the per-thread arm allocates `budget` as `threads` independent slices;
            // otherwise it is one shared pool and `budget / threads` would misdescribe it.
            if per_thread {
                log::debug!(
                    "Auto memory: {} of {} ({}/thread × {} threads, reserve {})",
                    bytesize::ByteSize(budget as u64),
                    bytesize::ByteSize(total as u64),
                    bytesize::ByteSize((budget / threads) as u64),
                    threads,
                    bytesize::ByteSize(margin as u64),
                );
            } else {
                log::debug!(
                    "Auto memory: {} of {} (shared across {} threads, reserve {})",
                    bytesize::ByteSize(budget as u64),
                    bytesize::ByteSize(total as u64),
                    threads,
                    bytesize::ByteSize(margin as u64),
                );
            }
            budget
        }
    };

    if budget > total {
        log::warn!(
            "Memory budget {} exceeds total host memory {}; this may cause OOM (or, for sort, earlier spill-to-disk)",
            bytesize::ByteSize(budget as u64),
            bytesize::ByteSize(total as u64),
        );
    }

    Ok(budget)
}

/// Parse a memory size string into `usize` bytes (private helper).
fn parse_memory_bytes(s: &str, label: &str) -> std::result::Result<usize, String> {
    let bytes = parse_memory_size(s).map_err(|e| e.to_string())?;
    usize::try_from(bytes).map_err(|_| format!("{label} too large: {bytes}"))
}

/// Parses a boolean value from a string, accepting: true/false, yes/no, y/n, t/f
/// (case-insensitive). Matches sopt/fgbio behavior.
///
/// # Errors
///
/// Returns an error string if the input is not a recognized boolean.
pub fn parse_bool(s: &str) -> std::result::Result<bool, String> {
    match s.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => Ok(true),
        "false" | "f" | "no" | "n" => Ok(false),
        _ => Err(format!("Invalid boolean value '{s}'. Expected: true|false|yes|no|y|n|t|f")),
    }
}

/// Options for output compression.
///
/// Controls BGZF compression level for BAM output files.
#[derive(Debug, Clone, clap::Args)]
pub struct CompressionOptions {
    /// Compression level for output BAM (0-12).
    ///
    /// Level 0 disables compression (uncompressed BGZF blocks).
    /// Level 1 is fastest of the compressing levels with larger files;
    /// level 12 produces the smallest files but is slowest.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u32).range(0..=12))]
    pub compression_level: u32,
}

impl Default for CompressionOptions {
    /// Mirrors the clap `default_value_t = 1` so programmatic/default-constructed
    /// callers emit level-1 compression rather than the `u32` default of `0`
    /// (uncompressed BGZF).
    fn default() -> Self {
        Self { compression_level: 1 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::time::Duration;

    #[test]
    fn test_detect_total_memory_nonzero() {
        let total = detect_total_memory();
        assert!(total > 0, "expected non-zero total memory, got {total}");
    }

    #[test]
    fn test_detect_cpu_count_at_least_one() {
        assert!(detect_cpu_count() >= 1);
    }

    // `format_count`, `format_duration` and `format_rate` are re-exports from `fgumi-fmt`,
    // which owns their case tables and doctests. The `pub use` is compile-checked, so
    // re-asserting their behavior here would only duplicate coverage.

    #[test]
    fn test_parse_memory_size_plain() {
        assert_eq!(parse_memory_size("768").unwrap(), 768 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_size_human() {
        assert_eq!(parse_memory_size("2GB").unwrap(), 2 * 1000 * 1000 * 1000);
        assert_eq!(parse_memory_size("512MiB").unwrap(), 512 * 1024 * 1024);
    }

    #[test]
    fn test_parse_memory_size_errors() {
        assert!(parse_memory_size("").is_err());
        assert!(parse_memory_size("-1").is_err());
        assert!(parse_memory_size("0").is_err());
    }

    /// Audit C3: the scientific-notation guard rejected any `e`/`E`, which also
    /// rejected the exabyte units `EB`/`EiB` that `bytesize` accepts. It must
    /// only fire for a genuine mantissa-exponent (a digit/`.`, then `e`/`E`,
    /// then a digit or sign).
    #[rstest]
    #[case::exabyte("1EB", true)]
    #[case::exbibyte("1EiB", true)]
    #[case::exabyte_multi("2EB", true)]
    #[case::gibibyte("512MiB", true)]
    #[case::scientific_lower("1e5", false)]
    #[case::scientific_upper("1E5", false)]
    #[case::scientific_decimal("1.5e10", false)]
    #[case::scientific_signed("2E+3", false)]
    #[case::scientific_neg("3e-2", false)]
    fn test_parse_memory_size_scientific_notation_guard(#[case] input: &str, #[case] ok: bool) {
        assert_eq!(
            parse_memory_size(input).is_ok(),
            ok,
            "parse_memory_size({input:?}) expected ok={ok}, got {:?}",
            parse_memory_size(input)
        );
    }

    #[test]
    fn memory_limit_auto_displays_as_auto() {
        assert_eq!(MemoryLimit::Auto.to_string(), "auto");
    }

    /// A rendered limit must survive a round trip through `parse_memory`, so a default shown
    /// in `--help` can be pasted back on the command line unchanged.
    #[rstest]
    #[case::gibibytes(2 * 1024 * 1024 * 1024)]
    #[case::mebibytes(768 * 1024 * 1024)]
    #[case::kibibytes(4 * 1024)]
    #[case::bare_bytes(100)]
    #[case::not_a_clean_multiple(1536)]
    fn memory_display_round_trips_through_parse_memory(#[case] bytes: usize) {
        let rendered = MemoryLimit::Fixed(bytes).to_string();
        assert_eq!(
            parse_memory(&rendered).unwrap(),
            MemoryLimit::Fixed(bytes),
            "{rendered} did not round-trip"
        );
    }

    #[test]
    fn test_resolve_reserve_auto() {
        // Target-aware: a fixed 32 GiB `usize` literal cannot compile on 32-bit (`usize == u32`).
        // Derive `total` from the cap so the test builds on every target, and assert the exact
        // `resolve_reserve` formula (`AUTO_RESERVE_CAP.min(total / 2)`): on 64-bit the 10 GiB cap
        // strictly binds (total/2 = 15 GiB > cap); on 32-bit the cap saturates to `usize::MAX`, so
        // the expected value is `total / 2` — matching the documented 32-bit saturation behavior.
        let total = AUTO_RESERVE_CAP.saturating_mul(3);
        let r = resolve_reserve(MemoryReserve::Auto, total);
        assert_eq!(r, AUTO_RESERVE_CAP.min(total / 2));
    }

    #[test]
    fn test_resolve_memory_budget_fixed() {
        let budget = resolve_memory_budget(
            MemoryLimit::Fixed(512 * 1024 * 1024),
            MemoryReserve::Auto,
            4,
            true,
        )
        .unwrap();
        assert_eq!(budget, 4 * 512 * 1024 * 1024);
    }

    #[rstest]
    #[case::true_word("true", Some(true))]
    #[case::t_short("t", Some(true))]
    #[case::yes_word("yes", Some(true))]
    #[case::y_short("y", Some(true))]
    #[case::false_word("false", Some(false))]
    #[case::f_short("f", Some(false))]
    #[case::no_word("no", Some(false))]
    #[case::n_short("n", Some(false))]
    #[case::mixed_case("TrUe", Some(true))]
    #[case::garbage("maybe", None)]
    #[case::empty("", None)]
    fn parse_bool_accepts_the_documented_spellings(
        #[case] input: &str,
        #[case] expected: Option<bool>,
    ) {
        assert_eq!(parse_bool(input).ok(), expected);
    }

    #[derive(clap::Parser)]
    struct CompressionHarness {
        #[command(flatten)]
        compression: CompressionOptions,
    }

    // In-range values (the 0-12 bounds; 0 = uncompressed) parse.
    #[rstest]
    #[case::uncompressed(0_u32)]
    #[case::fastest(1_u32)]
    #[case::midpoint(6_u32)]
    #[case::maximum(12_u32)]
    fn test_compression_level_accepts_in_range(#[case] level: u32) {
        use clap::Parser;

        assert_eq!(
            CompressionHarness::try_parse_from(["prog", "--compression-level", &level.to_string()])
                .unwrap()
                .compression
                .compression_level,
            level
        );
    }

    // The programmatic `Default` must match the clap `default_value_t = 1`; a
    // derived `Default` would silently yield level 0 (uncompressed BGZF).
    #[test]
    fn test_compression_options_default_matches_cli_default() {
        assert_eq!(CompressionOptions::default().compression_level, 1);
    }

    #[test]
    fn test_compression_level_default_and_rejects_out_of_range() {
        use clap::Parser;

        // Default is 1 when the flag is omitted.
        assert_eq!(CompressionHarness::parse_from(["prog"]).compression.compression_level, 1);

        // Out-of-range values are rejected at parse time rather than silently accepted.
        assert!(CompressionHarness::try_parse_from(["prog", "--compression-level", "13"]).is_err());
        assert!(CompressionHarness::try_parse_from(["prog", "--compression-level", "99"]).is_err());
    }

    #[test]
    fn test_resolve_memory_budget_auto_low_available() {
        // When available memory per thread is below MIN_MEMORY_PER_THREAD, the
        // budget is floored to MIN_MEMORY_PER_THREAD × threads, then capped at
        // available (which is less), so the result equals available.
        let total = 2 * MIN_MEMORY_PER_THREAD; // very tight: only 2 × floor per 4 threads
        let reserve = 0;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Fixed(reserve),
            4,
            true,
            total,
        )
        .unwrap();
        // available = total - 0 = 2×MIN; per-thread = 2×MIN/4 < MIN → floored to MIN;
        // target = MIN × 4 = 4×MIN > available → capped at available.
        assert_eq!(budget, total);
    }

    #[test]
    fn test_resolve_memory_budget_auto_margin_exceeds_total() {
        // When the reserve margin >= total, saturating_sub → 0 available.
        // Budget is then floored to MIN_MEMORY_PER_THREAD, capped at 0
        // (available), so result == 0 (the cap wins).
        let total = 512 * 1024 * 1024_usize; // 512 MiB
        let margin = total + 1; // margin exceeds total
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Fixed(margin),
            1,
            false,
            total,
        )
        .unwrap();
        // available = total.saturating_sub(margin) = 0; target = max(0, MIN) = MIN;
        // budget = MIN.min(0) = 0.
        assert_eq!(budget, 0);
    }

    #[test]
    fn test_resolve_memory_budget_threads_zero_bails() {
        // `threads == 0` is an invalid input and must bail before any
        // per-thread arithmetic.
        let err = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(512 * 1024 * 1024),
            MemoryReserve::Auto,
            0,
            true,
            32 * 1024 * 1024 * 1024,
        )
        .expect_err("threads == 0 must be rejected");
        let msg = format!("{err:#}");
        assert!(msg.contains("--threads must be at least 1"), "got: {msg}");
    }

    #[test]
    fn test_resolve_memory_budget_fixed_per_thread_overflow() {
        // Fixed(usize::MAX) × 2 threads overflows the `checked_mul`, which must
        // surface an error rather than wrap.
        let err = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(usize::MAX),
            MemoryReserve::Auto,
            2,
            true,
            32 * 1024 * 1024 * 1024,
        )
        .expect_err("Fixed per-thread multiplication must overflow");
        let msg = format!("{err:#}");
        assert!(msg.contains("overflowed"), "got: {msg}");
    }

    #[test]
    fn test_resolve_memory_budget_auto_per_thread_overflow() {
        // Construct an Auto config whose per-thread floor × threads overflows:
        // with `total == usize::MAX` and a tiny reserve, `available` is enormous,
        // so `(available / threads).max(MIN) * threads` overflows `usize`.
        let err = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Fixed(0),
            usize::MAX,
            true,
            usize::MAX,
        )
        .expect_err("Auto per-thread multiplication must overflow");
        let msg = format!("{err:#}");
        assert!(msg.contains("auto memory budget overflowed"), "got: {msg}");
    }

    /// `format_rate` has three regimes: a sub-millisecond guard that avoids dividing by ~zero,
    /// a normal items/s path, and an items/min path once the rate drops below one per second.
    #[rstest]
    #[case::sub_millisecond_guard(5, Duration::from_micros(100), "5 items/s")]
    #[case::exactly_one_per_second(1, Duration::from_secs(1), "1 items/s")]
    #[case::thousands_per_second(1000, Duration::from_secs(1), "1,000 items/s")]
    #[case::below_one_per_second(30, Duration::from_secs(60), "30.0 items/min")]
    #[case::far_below_one_per_second(1, Duration::from_secs(120), "0.5 items/min")]
    fn format_rate_switches_units_below_one_per_second(
        #[case] count: u64,
        #[case] duration: Duration,
        #[case] expected: &str,
    ) {
        assert_eq!(format_rate(count, duration), expected);
    }

    #[test]
    fn operation_timer_reports_completion_without_panicking() {
        let timer = OperationTimer::new("test operation");
        timer.log_completion(1234);
        assert!(timer.start_time.elapsed() >= Duration::ZERO);
        assert_eq!(timer.operation, "test operation");
    }

    #[test]
    fn validate_file_exists_accepts_a_real_file() {
        let file = tempfile::NamedTempFile::new().expect("failed to create temp file");
        validate_file_exists(file.path(), "input BAM").expect("existing file must validate");
    }

    #[test]
    fn validate_file_exists_rejects_a_missing_path() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let missing = dir.path().join("definitely-not-here.bam");
        let err = validate_file_exists(&missing, "input BAM")
            .expect_err("a missing path must not validate");
        let msg = err.to_string();
        assert!(msg.contains("input BAM"), "error should name the description; got: {msg}");
        assert!(msg.contains("File does not exist"), "got: {msg}");
    }

    /// The error arms of `parse_memory_size` that the happy-path tests do not reach. Each input
    /// is rejected for a different documented reason, so the assertions pin the reason text.
    #[rstest]
    #[case::plain_number_too_large("2000000", "Plain number memory size too large")]
    #[case::plain_decimal("1.5", "Plain decimal numbers not supported")]
    #[case::human_readable_zero("0B", "Memory size cannot be zero")]
    #[case::unparseable("not-a-size", "Invalid memory size")]
    fn parse_memory_size_rejects_with_a_specific_reason(
        #[case] input: &str,
        #[case] expected_reason: &str,
    ) {
        let err = parse_memory_size(input).expect_err("input must be rejected");
        let msg = err.to_string();
        assert!(msg.contains(expected_reason), "expected {expected_reason:?}, got: {msg}");
    }

    #[test]
    fn memory_limit_default_is_768_mib() {
        assert_eq!(MemoryLimit::default(), MemoryLimit::Fixed(768 * 1024 * 1024));
    }

    /// `format_binary_bytes` picks the largest binary unit that divides cleanly, falling back to
    /// bytes. It backs the `Display` impls of both `MemoryLimit` and `MemoryReserve`.
    #[rstest]
    #[case::gibibytes(2 * 1024 * 1024 * 1024, "2GiB")]
    #[case::mebibytes(512 * 1024 * 1024, "512MiB")]
    #[case::kibibytes(4 * 1024, "4KiB")]
    #[case::bare_bytes(100, "100B")]
    #[case::not_a_clean_multiple(1536, "1536B")]
    fn memory_display_uses_the_largest_clean_binary_unit(
        #[case] bytes: usize,
        #[case] expected: &str,
    ) {
        assert_eq!(MemoryLimit::Fixed(bytes).to_string(), expected);
        assert_eq!(MemoryReserve::Fixed(bytes).to_string(), expected);
    }

    #[test]
    fn memory_reserve_auto_displays_as_auto() {
        assert_eq!(MemoryReserve::Auto.to_string(), "auto");
        assert_eq!(MemoryReserve::default(), MemoryReserve::Auto);
    }

    #[rstest]
    #[case::lowercase_auto("auto", MemoryLimit::Auto)]
    #[case::uppercase_auto("AUTO", MemoryLimit::Auto)]
    #[case::surrounding_whitespace("  auto  ", MemoryLimit::Auto)]
    #[case::human_readable("512MiB", MemoryLimit::Fixed(512 * 1024 * 1024))]
    #[case::plain_number_is_mib("768", MemoryLimit::Fixed(768 * 1024 * 1024))]
    fn parse_memory_accepts_auto_and_sizes(#[case] input: &str, #[case] expected: MemoryLimit) {
        assert_eq!(parse_memory(input).unwrap(), expected);
    }

    #[rstest]
    #[case::lowercase_auto("auto", MemoryReserve::Auto)]
    #[case::uppercase_auto("AUTO", MemoryReserve::Auto)]
    #[case::human_readable("10GiB", MemoryReserve::Fixed(10 * 1024 * 1024 * 1024))]
    fn parse_memory_reserve_accepts_auto_and_sizes(
        #[case] input: &str,
        #[case] expected: MemoryReserve,
    ) {
        assert_eq!(parse_memory_reserve(input).unwrap(), expected);
    }

    #[test]
    fn parse_memory_surfaces_the_underlying_parse_error() {
        let err = parse_memory("not-a-size").expect_err("garbage must not parse");
        assert!(err.contains("Invalid memory size"), "got: {err}");
        let err = parse_memory_reserve("-1").expect_err("negative must not parse");
        assert!(err.contains("cannot be negative"), "got: {err}");
    }

    #[test]
    fn resolve_memory_budget_fixed_without_per_thread_ignores_thread_count() {
        // The non-per-thread Fixed arm returns the limit verbatim, regardless of `threads`.
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(512 * 1024 * 1024),
            MemoryReserve::Auto,
            8,
            false,
            16 * 1024 * 1024 * 1024,
        )
        .unwrap();
        assert_eq!(budget, 512 * 1024 * 1024);
    }

    #[test]
    fn resolve_memory_budget_allows_a_limit_above_total_host_memory() {
        // Over-committing is a warning, not an error: the caller may know better (e.g. swap), and
        // for sort it just means spilling to disk earlier.
        let total = 1024 * 1024 * 1024;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(4 * 1024 * 1024 * 1024),
            MemoryReserve::Auto,
            1,
            false,
            total,
        )
        .unwrap();
        assert_eq!(budget, 4 * 1024 * 1024 * 1024);
        assert!(budget > total, "this case must exercise the over-commit path");
    }
}
