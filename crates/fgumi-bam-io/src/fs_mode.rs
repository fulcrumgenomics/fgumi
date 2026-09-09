//! Output file-permission helper for the atomic temp-file-then-rename pattern.
//!
//! The commands and sinks in the fgumi family write their outputs atomically:
//! write to a `tempfile::NamedTempFile`, then `persist` it onto the final path
//! so a reader never observes a partial file and a failed write leaves nothing
//! behind. `NamedTempFile` hard-codes an owner-only `0o600` mode (a security
//! default for scratch files), and `persist` is a rename that keeps that mode —
//! so, without intervention, an output written this way lands `0o600` instead of
//! the mode a plain `File::create` would have produced. That surprises group
//! readers and downstream tools running as another user.
//!
//! [`restamp_for_persist`] re-stamps the temp file to that `File::create` mode
//! before it is persisted. It derives the new-file mode by observing a one-time
//! probe `File::create` (the kernel applies the umask atomically at open time),
//! so this crate stays `#![deny(unsafe_code)]` — no libc / `unsafe` — and never
//! mutates the process-wide umask. It is the single home for this logic across
//! the mode-preserving callers: the `compare bam-roundtrip` command, the BAI
//! index writer, and `fgumi-sort`'s merge output all route through it.
//!
//! A sibling module (`fgumi-metrics`'s `write_metrics_atomic`) solves the same
//! "don't inherit `NamedTempFile`'s `0o600`" problem the other way — it stages
//! the temp via `File::create` so the kernel applies the umask at creation. That
//! is simpler but only covers the new-file case and couples correctness to *how*
//! each caller creates its temp. This helper is instead **self-contained**: it
//! is correct no matter how the temp was created, and it additionally preserves
//! an existing destination's mode on overwrite (which the merge output requires,
//! and the create-time approach cannot express). Keeping the temp owner-only
//! until this final re-stamp also narrows the window in which the not-yet-
//! finished output is group-readable.

use std::fs::File;
use std::io;
use std::path::Path;

/// Re-stamp `temp`'s permissions to the mode a plain `File::create(dest)` would
/// leave the output with, so an atomic `NamedTempFile` → `persist(dest)` does
/// not inherit `NamedTempFile`'s owner-only `0o600`.
///
/// Matches `File::create` semantics on Unix:
/// - `dest` already exists → keep its current mode (create-truncate never
///   changes an existing file's permissions; `std::fs::metadata` follows
///   symlinks, so an existing symlink target's mode is used);
/// - `dest` is new → `0o666 & !umask`.
///
/// Call this on the temp file **before** persisting it onto `dest`. No-op on
/// non-Unix targets (no umask model; `NamedTempFile` does not force `0o600`).
///
/// # Errors
///
/// Returns the I/O error from `set_permissions` if the mode cannot be applied.
#[cfg_attr(not(unix), allow(clippy::unnecessary_wraps))]
pub fn restamp_for_persist(temp: &File, dest: &Path) -> io::Result<()> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mode = target_file_mode(dest);
        temp.set_permissions(std::fs::Permissions::from_mode(mode))?;
    }
    #[cfg(not(unix))]
    {
        let _ = (temp, dest);
    }
    Ok(())
}

/// The `0o777` permission bits a `File::create(dest)` would leave `dest` with:
/// the existing file's permission bits when `dest` exists, else `0o666 & !umask`.
///
/// Only the low `0o777` bits are considered; setuid/setgid/sticky are not
/// re-applied (they are vanishingly rare on a BAM/index output, and this matches
/// the merge-output path this helper replaces).
#[cfg(unix)]
fn target_file_mode(dest: &Path) -> u32 {
    use std::os::unix::fs::PermissionsExt;
    match std::fs::metadata(dest) {
        // Overwriting: `File::create` never changes an existing file's mode.
        Ok(meta) => meta.permissions().mode() & 0o777,
        // New file: the mode `File::create` would leave, `0o666 & !umask`.
        Err(_) => new_file_mode(),
    }
}

/// The mode a fresh `File::create` produces (`0o666 & !umask`), observed once by
/// creating a throwaway probe file and reading its resulting permissions.
///
/// POSIX has no read-only `umask(2)`: the value can only be read by *setting* a
/// new mask and capturing the previous one. Doing that — even set-to-`0` then
/// restore — mutates the **process-wide** umask for the duration, so a
/// `File::create` running concurrently on another thread (fgumi-sort creates
/// spill files in parallel) could momentarily observe the cleared mask and land
/// a world-writable file. Creating a probe file instead lets the kernel apply
/// the umask atomically at `open(O_CREAT)` time, reading the effective new-file
/// mode without ever touching the shared umask. The result is cached in a
/// process-wide `OnceLock`, so the probe's handful of syscalls run at most once.
///
/// If the probe cannot be created (e.g. no writable temp dir), fall back to
/// `0o644` — the mode `File::create` yields under the ubiquitous `umask 022`,
/// and never world-writable.
#[cfg(unix)]
fn new_file_mode() -> u32 {
    use std::sync::OnceLock;
    static MODE: OnceLock<u32> = OnceLock::new();
    *MODE.get_or_init(|| {
        probe_new_file_mode().unwrap_or_else(|e| {
            log::warn!("could not probe the umask; new outputs fall back to mode 0o644: {e}");
            0o644
        })
    })
}

/// Create a throwaway file via `File::create` and return its `0o777` permission
/// bits — i.e. `0o666 & !umask`, with the umask applied by the kernel rather
/// than read from the process. The probe file and its enclosing temp dir are
/// removed when the `TempDir` drops at the end of this function.
#[cfg(unix)]
fn probe_new_file_mode() -> io::Result<u32> {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir()?;
    let probe = File::create(dir.path().join("umask-probe"))?;
    Ok(probe.metadata()?.permissions().mode() & 0o777)
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    fn mode_of(path: &Path) -> u32 {
        std::fs::metadata(path).expect("stat").permissions().mode() & 0o777
    }

    /// A persisted `NamedTempFile` re-stamped for a **new** destination lands at
    /// the same mode a plain `File::create` produces — i.e. no longer the
    /// owner-only `0o600` that `NamedTempFile` forces. Compared against a live
    /// `File::create` so the assertion tracks the runner's ambient umask rather
    /// than hard-coding `0o644`.
    #[test]
    fn new_dest_matches_file_create_mode() {
        let dir = tempfile::tempdir().expect("temp dir");
        let reference_mode = {
            let p = dir.path().join("reference");
            std::fs::File::create(&p).expect("create reference");
            mode_of(&p)
        };

        let tmp = tempfile::NamedTempFile::new_in(dir.path()).expect("temp");
        assert_eq!(mode_of(tmp.path()), 0o600, "NamedTempFile should start owner-only");

        let dest = dir.path().join("out"); // does not exist yet
        restamp_for_persist(tmp.as_file(), &dest).expect("restamp");
        tmp.persist(&dest).expect("persist");
        assert_eq!(mode_of(&dest), reference_mode);
    }

    /// Overwriting an **existing** destination preserves that file's mode, exactly
    /// as `File::create` (open-truncate) does — not the umask default.
    #[test]
    fn existing_dest_mode_is_preserved() {
        let dir = tempfile::tempdir().expect("temp dir");
        let dest = dir.path().join("out");
        std::fs::File::create(&dest).expect("create dest");
        std::fs::set_permissions(&dest, std::fs::Permissions::from_mode(0o640))
            .expect("chmod dest");

        let tmp = tempfile::NamedTempFile::new_in(dir.path()).expect("temp");
        restamp_for_persist(tmp.as_file(), &dest).expect("restamp");
        tmp.persist(&dest).expect("persist");

        assert_eq!(mode_of(&dest), 0o640, "overwriting must preserve the existing file's mode");
    }

    /// Triggering the internal new-file-mode probe must not perturb the mode a
    /// plain `File::create` produces. The probe observes the umask by creating a
    /// file, so — unlike a `umask(0)`/restore read — it never mutates the
    /// process-wide umask a concurrent `File::create` depends on. Assert a
    /// `File::create` lands the same mode before and after a restamp. nextest
    /// runs each test in its own process, so the `OnceLock` probe here is
    /// genuinely first-init.
    #[test]
    fn probe_does_not_disturb_file_create_mode() {
        let dir = tempfile::tempdir().expect("temp dir");

        std::fs::File::create(dir.path().join("before")).expect("create before");
        let mode_before = mode_of(&dir.path().join("before"));

        // Trigger the internal umask probe via a restamp for a new destination.
        let tmp = tempfile::NamedTempFile::new_in(dir.path()).expect("temp");
        restamp_for_persist(tmp.as_file(), &dir.path().join("new")).expect("restamp");

        std::fs::File::create(dir.path().join("after")).expect("create after");
        let mode_after = mode_of(&dir.path().join("after"));

        assert_eq!(mode_before, mode_after, "the umask probe must not change File::create's mode");
    }
}
