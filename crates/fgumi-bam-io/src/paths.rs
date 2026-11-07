//! Path utility functions for stdin and stdout detection.

use std::path::Path;

/// Returns `true` if the path refers to stdin (`-` or `/dev/stdin`).
///
/// # Examples
///
/// ```
/// use fgumi_bam_io::paths::is_stdin_path;
/// use std::path::Path;
///
/// assert!(is_stdin_path(Path::new("-")));
/// assert!(is_stdin_path(Path::new("/dev/stdin")));
/// assert!(!is_stdin_path(Path::new("input.bam")));
/// ```
pub fn is_stdin_path<P: AsRef<Path>>(path: P) -> bool {
    let path_str = path.as_ref().to_string_lossy();
    path_str == "-" || path_str == "/dev/stdin"
}

/// Returns `true` if the input at `path` arrives as a stream rather than from a
/// regular file: a FIFO, a pipe, a `/dev/fd/N` from process substitution, a
/// character device. For stdin (`-` or `/dev/stdin`) that is whatever fd 0 is,
/// so `-i - < in.bam` is a regular file and `cat in.bam | ... -i -` a stream. A
/// path that cannot be stat'ed, or a directory, is not a stream; opening it
/// reports the problem.
///
/// Only a stream has a writer that can block when its reader pauses, which is
/// what decides whether reading it through a
/// [`ReadBuffer`](crate::read_buffer::ReadBuffer) is worth its memory. `stat`
/// does not open the path, so a FIFO with no writer yet does not block here.
///
/// Two other `stat`-based checks answer different questions on purpose and are
/// not this one: `sam_input`'s `input_is_pipe_like` treats stdin as pipe-like
/// whatever fd 0 is and a character device as not (it only words the error for
/// an empty input, where `-` and a FIFO point back upstream), and
/// [`open_seekable_bgzf`](crate::pread::open_seekable_bgzf) reads only a
/// regular file positionally, so stdin, a directory and a stat failure are all
/// "not seekable" there (the stat failure an error) rather than "not a stream".
#[must_use]
pub fn is_stream_input<P: AsRef<Path>>(path: P) -> bool {
    let metadata = if is_stdin_path(&path) { stdin_metadata() } else { std::fs::metadata(path) };
    metadata.is_ok_and(|meta| !meta.is_file() && !meta.is_dir())
}

/// The metadata of whatever fd 0 is, from `fstat` on a duplicate of it.
#[cfg(unix)]
fn stdin_metadata() -> std::io::Result<std::fs::Metadata> {
    use std::os::fd::AsFd;
    let fd = std::io::stdin().as_fd().try_clone_to_owned()?;
    std::fs::File::from(fd).metadata()
}

/// Off Unix, stdin cannot be inspected; reported as unavailable.
#[cfg(not(unix))]
fn stdin_metadata() -> std::io::Result<std::fs::Metadata> {
    Err(std::io::Error::other("stdin metadata is not available on this platform"))
}

/// Returns `true` if the path refers to stdout (`-` or `/dev/stdout`).
///
/// # Examples
///
/// ```
/// use fgumi_bam_io::paths::is_stdout_path;
/// use std::path::Path;
///
/// assert!(is_stdout_path(Path::new("-")));
/// assert!(is_stdout_path(Path::new("/dev/stdout")));
/// assert!(!is_stdout_path(Path::new("output.bam")));
/// ```
pub fn is_stdout_path<P: AsRef<Path>>(path: P) -> bool {
    let path_str = path.as_ref().to_string_lossy();
    path_str == "-" || path_str == "/dev/stdout"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_stdin_path() {
        assert!(is_stdin_path("-"));
        assert!(is_stdin_path("/dev/stdin"));
        assert!(is_stdin_path(Path::new("-")));
        assert!(is_stdin_path(Path::new("/dev/stdin")));

        assert!(!is_stdin_path("input.bam"));
        assert!(!is_stdin_path("/path/to/file.bam"));
        assert!(!is_stdin_path(""));
        assert!(!is_stdin_path("/dev/null"));
    }

    #[test]
    fn test_is_stream_input() {
        let dir = tempfile::tempdir().expect("tempdir");
        let file = dir.path().join("in.bam");
        std::fs::write(&file, b"x").expect("write");
        // Stdin's answer depends on what the test runner attached to fd 0, so
        // it is covered by the CLI tests instead.
        assert!(is_stream_input("/dev/null"));
        assert!(!is_stream_input(&file));
        assert!(!is_stream_input(dir.path().join("missing.bam")));
        assert!(!is_stream_input(dir.path()));
    }

    #[test]
    fn test_is_stdout_path() {
        assert!(is_stdout_path("-"));
        assert!(is_stdout_path("/dev/stdout"));
        assert!(is_stdout_path(Path::new("-")));
        assert!(is_stdout_path(Path::new("/dev/stdout")));

        assert!(!is_stdout_path("output.bam"));
        assert!(!is_stdout_path("/path/to/file.bam"));
        assert!(!is_stdout_path(""));
        assert!(!is_stdout_path("/dev/null"));
    }
}
