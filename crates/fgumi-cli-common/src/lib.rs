#![deny(unsafe_code)]
//! Shared CLI types and helpers for fgumi commands.

// ─────────────────────────────────────────────────────────────────────────────
// Command trait
// ─────────────────────────────────────────────────────────────────────────────

use enum_dispatch::enum_dispatch;

/// Trait implemented by all fgumi CLI commands.
///
/// Each command provides an `execute` method that runs the command's main logic.
/// The `command_line` parameter contains the full command invocation for @PG records.
#[enum_dispatch]
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
/// Uses `num_cpus::get()`, which attempts to honor cgroup CPU quotas (e.g.
/// `--cpus` in Docker or Kubernetes resource limits) but may fall back to the
/// physical core count — its cgroup v2 quota handling is incomplete (see
/// [num_cpus#122](https://github.com/seanmonstar/num_cpus/issues/122)). Returns
/// at least 1.
#[must_use]
pub fn detect_cpu_count() -> usize {
    num_cpus::get().max(1)
}

// ─────────────────────────────────────────────────────────────────────────────
// Formatting helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Format an integer with comma separators (e.g. 1234567 → "1,234,567").
#[must_use]
pub fn format_count(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    let offset = s.len() % 3;
    for (i, c) in s.chars().enumerate() {
        // Insert a comma before every group of three digits, but not before
        // the first character. `offset` is the length of the leading partial
        // group, so commas fall at indices `offset, offset+3, offset+6, …`
        // (skipping index 0 when `offset == 0`).
        if i >= offset && i > 0 && (i - offset).is_multiple_of(3) {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// Formats a duration in human-readable form.
///
/// # Examples
///
/// ```
/// use fgumi_cli_common::format_duration;
/// use std::time::Duration;
///
/// assert_eq!(format_duration(Duration::from_secs(45)), "45s");
/// assert_eq!(format_duration(Duration::from_secs(135)), "2m 15s");
/// assert_eq!(format_duration(Duration::from_secs(5400)), "1h 30m");
/// ```
#[must_use]
pub fn format_duration(duration: std::time::Duration) -> String {
    let secs = duration.as_secs();
    if secs < 60 {
        format!("{secs}s")
    } else if secs < 3600 {
        let mins = secs / 60;
        let remaining_secs = secs % 60;
        if remaining_secs == 0 { format!("{mins}m") } else { format!("{mins}m {remaining_secs}s") }
    } else {
        let hours = secs / 3600;
        let mins = (secs % 3600) / 60;
        if mins == 0 { format!("{hours}h") } else { format!("{hours}h {mins}m") }
    }
}

/// Formats a rate (items per second) with appropriate units.
///
/// # Examples
///
/// ```
/// use fgumi_cli_common::format_rate;
/// use std::time::Duration;
///
/// assert_eq!(format_rate(1000, Duration::from_secs(1)), "1,000 items/s");
/// assert_eq!(format_rate(600, Duration::from_secs(60)), "10 items/s");
/// ```
#[must_use]
#[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn format_rate(count: u64, duration: std::time::Duration) -> String {
    let secs = duration.as_secs_f64();
    if secs < 0.001 {
        return format!("{} items/s", format_count(count));
    }

    let rate = count as f64 / secs;
    if rate >= 1.0 {
        format!("{} items/s", format_count(rate as u64))
    } else {
        let items_per_min = count as f64 / (secs / 60.0);
        format!("{items_per_min:.1} items/min")
    }
}

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

/// Parses a memory size string into bytes.
///
/// Accepts both plain numbers (interpreted as MiB) and human-readable formats like:
/// - "2GB", "2G" -> 2 gigabytes (decimal: 2,000,000,000)
/// - "1.5GB" -> 1.5 gigabytes
/// - "1024MB", "1024M" -> 1024 megabytes (decimal)
/// - "512MiB" -> 512 mebibytes (binary: 536,870,912)
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
                    "Plain number memory size too large: {} MiB. Use human-readable format like '{}GB' instead.",
                    mb_value,
                    mb_value / 1000
                ),
            });
        }
        return mb_value.checked_mul(1024 * 1024).ok_or_else(|| FgumiError::InvalidMemorySize {
            reason: format!("Memory size calculation overflow for {mb_value} MiB"),
        });
    }

    if trimmed.contains('e') || trimmed.contains('E') {
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
const AUTO_RESERVE_CAP: usize = 10 * 1024 * 1024 * 1024;

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
            log::debug!(
                "Auto memory: {} of {} ({}/thread × {} threads, reserve {})",
                bytesize::ByteSize(budget as u64),
                bytesize::ByteSize(total as u64),
                bytesize::ByteSize((budget / threads) as u64),
                threads,
                bytesize::ByteSize(margin as u64),
            );
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

    #[test]
    fn test_detect_total_memory_nonzero() {
        let total = detect_total_memory();
        assert!(total > 0, "expected non-zero total memory, got {total}");
    }

    #[test]
    fn test_detect_cpu_count_at_least_one() {
        assert!(detect_cpu_count() >= 1);
    }

    #[test]
    fn test_format_count_no_commas() {
        assert_eq!(format_count(0), "0");
        assert_eq!(format_count(999), "999");
    }

    #[test]
    fn test_format_count_with_commas() {
        assert_eq!(format_count(1000), "1,000");
        assert_eq!(format_count(1_234_567), "1,234,567");
    }

    #[test]
    fn test_format_duration() {
        use std::time::Duration;
        assert_eq!(format_duration(Duration::from_secs(45)), "45s");
        assert_eq!(format_duration(Duration::from_secs(135)), "2m 15s");
        assert_eq!(format_duration(Duration::from_secs(5400)), "1h 30m");
    }

    #[test]
    fn test_format_rate() {
        use std::time::Duration;
        assert_eq!(format_rate(1000, Duration::from_secs(1)), "1,000 items/s");
    }

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

    #[test]
    fn test_memory_limit_display() {
        assert_eq!(MemoryLimit::Auto.to_string(), "auto");
        assert_eq!(MemoryLimit::Fixed(768 * 1024 * 1024).to_string(), "768MiB");
        assert_eq!(MemoryLimit::Fixed(1024 * 1024 * 1024).to_string(), "1GiB");
    }

    #[test]
    fn test_resolve_reserve_auto() {
        let total = 32 * 1024 * 1024 * 1024_usize; // 32 GiB
        let r = resolve_reserve(MemoryReserve::Auto, total);
        assert_eq!(r, AUTO_RESERVE_CAP); // capped at 10 GiB
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

    #[test]
    fn test_parse_bool() {
        assert!(parse_bool("true").unwrap());
        assert!(!parse_bool("false").unwrap());
        assert!(parse_bool("yes").unwrap());
        assert!(!parse_bool("no").unwrap());
        assert!(parse_bool("maybe").is_err());
    }

    #[derive(clap::Parser)]
    struct CompressionHarness {
        #[command(flatten)]
        compression: CompressionOptions,
    }

    // In-range values (the 0-12 bounds; 0 = uncompressed) parse.
    #[rstest]
    #[case(0_u32)]
    #[case(1_u32)]
    #[case(6_u32)]
    #[case(12_u32)]
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
}
