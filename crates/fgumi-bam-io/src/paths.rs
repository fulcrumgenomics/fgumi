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
