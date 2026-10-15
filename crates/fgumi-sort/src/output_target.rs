//! The atomic output of a merge: stage in a same-directory temp, persist only
//! on success.

use std::path::{Path, PathBuf};

use anyhow::{Context as _, Result};
use fgumi_bam_io::is_stdout_path;

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
/// direct write (`File::create`), [`persist`](MergeOutputTarget::persist) re-stamps
/// the temp to the mode the final file must carry — via the shared
/// [`fgumi_bam_io::restamp_for_persist`] — before the rename.
pub enum MergeOutputTarget {
    /// A regular-file output, staged in a same-directory temp. `dest` is the
    /// resolved final path the temp is atomically renamed onto (a symlinked
    /// output is followed to its real target, so the temp is staged next to —
    /// and renamed onto — the linked file rather than the link itself). The temp
    /// is re-stamped to `dest`'s `File::create` mode at persist time (see
    /// [`persist`](MergeOutputTarget::persist)).
    Temp {
        /// The staged temp.
        temp: tempfile::NamedTempFile,
        /// The resolved destination.
        dest: PathBuf,
    },
    /// Stdout, written directly.
    Stdout(PathBuf),
}

impl MergeOutputTarget {
    /// Stage `output` (see the type docs).
    ///
    /// # Errors
    ///
    /// An existing destination that is not a regular file, a symlink chain that
    /// does not resolve, or a temp that cannot be created.
    pub fn create(output: &Path) -> Result<Self> {
        Self::create_for(output, "merge")
    }

    /// [`Self::create`] for `command`'s output: `command` names the staged
    /// temp (`.fgumi-<command>-*.tmp`) and the error for a non-regular
    /// destination.
    ///
    /// # Errors
    ///
    /// As [`Self::create`].
    pub fn create_for(output: &Path, command: &str) -> Result<Self> {
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
                "{command} output must be a regular file or stdout: {}",
                dest.display()
            );
        }
        let dir = dest
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
        let temp = tempfile::Builder::new()
            .prefix(&format!(".fgumi-{command}-"))
            .suffix(".tmp")
            .tempfile_in(&dir)
            .with_context(|| {
                format!("failed to create a {command} temp file in {}", dir.display())
            })?;
        Ok(Self::Temp { temp, dest })
    }

    /// The path the merged BAM is written to.
    #[must_use]
    pub fn path(&self) -> &Path {
        match self {
            Self::Temp { temp, .. } => temp.path(),
            Self::Stdout(path) => path,
        }
    }

    /// Renames the finished temp onto its resolved destination (a no-op for
    /// stdout), first stamping the resolved mode so the output is not left
    /// temp-private (`0600`). Consumes `self`, disarming the RAII auto-remove.
    ///
    /// The merge wrote the temp through its own descriptor, which its writer
    /// already closed with the result checked (`OutputFile`); the temp's own
    /// handle is closed here, also checked, before the rename. The rename is
    /// atomic against a failure of this process, not a host crash: the data is
    /// not synced first (see `fgumi_bam_io::persist_after_close`).
    ///
    /// # Errors
    ///
    /// Re-stamping, closing or renaming the temp failed.
    pub fn persist(self) -> Result<()> {
        self.persist_with(|file| fgumi_bam_io::OutputFile::from(file).close())
    }

    /// [`persist`](Self::persist) with the temp handle's close supplied.
    /// Production always passes a checked `OutputFile` close; the parameter
    /// exists so tests can make the close fail.
    fn persist_with<F>(self, close: F) -> Result<()>
    where
        F: FnOnce(std::fs::File) -> std::io::Result<()>,
    {
        if let Self::Temp { temp, dest } = self {
            // Re-stamp the temp to the mode `File::create(dest)` would produce
            // before the rename: `NamedTempFile` is `0600`, and `persist` keeps
            // that mode, so merged output would otherwise be owner-only. The
            // shared helper derives that mode by probing `File::create` once per
            // process (cached) rather than mutating the process-wide umask, so it
            // is safe under concurrent output creation regardless of when it
            // first runs.
            fgumi_bam_io::restamp_for_persist(temp.as_file(), &dest).with_context(|| {
                format!("failed to set mode on merge temp for {}", dest.display())
            })?;
            fgumi_bam_io::persist_after_close(temp, &dest, close).with_context(|| {
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
///
/// # Errors
///
/// A link that cannot be read, or too many levels of links.
pub(crate) fn resolve_symlink_output(output: &Path) -> Result<PathBuf> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    // ========================================================================
    // Merge output permissions
    // ========================================================================

    /// `MergeOutputTarget::persist` must re-stamp the temp to the mode a plain
    /// `File::create(dest)` would produce, so merged output is never left with
    /// `NamedTempFile`'s owner-only `0o600`. The mode-derivation itself is unit-
    /// tested in `fgumi-bam-io`; this is the merge-level guard that the
    /// `restamp_for_persist` call is actually wired into `persist` — without it,
    /// dropping that call ships `0o600` and nothing in this crate fails.
    ///
    /// Two cases mirror `File::create` semantics: a new destination lands at
    /// `0o666 & !umask` (compared against a live `File::create` reference so the
    /// assertion tracks the runner's umask), and overwriting an existing file
    /// keeps that file's mode.
    #[cfg(unix)]
    #[rstest]
    #[case::new_dest(None)]
    #[case::overwrite_existing(Some(0o640))]
    fn merge_persist_restamps_output_to_file_create_mode(#[case] preexisting_mode: Option<u32>) {
        use std::os::unix::fs::PermissionsExt;

        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let dest = dir.path().join("merged.bam");

        // The mode `File::create(dest)` would leave: an existing file's own mode
        // on overwrite, else `0o666 & !umask` observed via a live reference.
        let expected_mode = if let Some(mode) = preexisting_mode {
            std::fs::File::create(&dest).expect("pre-create destination");
            std::fs::set_permissions(&dest, std::fs::Permissions::from_mode(mode))
                .expect("chmod destination");
            mode
        } else {
            let reference = dir.path().join("reference");
            std::fs::File::create(&reference).expect("create reference");
            std::fs::metadata(&reference).expect("stat reference").permissions().mode() & 0o777
        };

        let target = MergeOutputTarget::create(&dest).expect("create merge target");
        target.persist().expect("persist merge output");

        let mode =
            std::fs::metadata(&dest).expect("stat merged output").permissions().mode() & 0o777;
        // Guard against the temp file's owner-only 0o600 leaking through -- but only when
        // File::create's own mode differs from 0o600. Under a restrictive umask (e.g. 0o077)
        // File::create legitimately yields 0o600, so this guard would otherwise fire on a
        // correct result; the assert_eq below fully pins the mode in that case.
        if expected_mode != 0o600 {
            assert_ne!(mode, 0o600, "merged output must not inherit the temp's owner-only 0o600");
        }
        assert_eq!(mode, expected_mode, "persisted merge output must match File::create's mode");
    }

    /// The merge temp is closed before it is renamed, and a failed close fails
    /// the merge: nothing is renamed onto `dest` (an existing output is left
    /// as it was) and the temp is removed.
    #[test]
    fn merge_persist_close_error_leaves_dest_and_removes_temp() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let dest = dir.path().join("merged.bam");
        std::fs::write(&dest, b"previous").expect("pre-create destination");

        let target = MergeOutputTarget::create(&dest).expect("create merge target");
        let temp_path = target.path().to_path_buf();
        std::fs::write(&temp_path, b"merged").expect("write merge temp");
        let err = target
            .persist_with(|_| {
                assert_eq!(std::fs::read(&dest).unwrap(), b"previous", "closed before rename");
                Err(std::io::Error::other("input/output error at close"))
            })
            .expect_err("a failed close must fail the merge");

        assert!(format!("{err:#}").contains("error at close"), "{err:#}");
        assert_eq!(std::fs::read(&dest).expect("read dest"), b"previous");
        assert!(!temp_path.exists(), "the merge temp must be removed");
    }

    /// A clean persist renames the merged bytes onto `dest` and leaves no temp.
    #[test]
    fn merge_persist_renames_temp_onto_dest() {
        let dir = tempfile::tempdir().expect("failed to create temp directory");
        let dest = dir.path().join("merged.bam");
        let target = MergeOutputTarget::create(&dest).expect("create merge target");
        let temp_path = target.path().to_path_buf();
        std::fs::write(&temp_path, b"merged").expect("write merge temp");
        target.persist().expect("persist merge output");
        assert_eq!(std::fs::read(&dest).expect("read dest"), b"merged");
        assert!(!temp_path.exists(), "the temp must be renamed away");
    }

    // ========================================================================
    // resolve_symlink_output
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
}
