//! Pipeline engine error types.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum PipelineError {
    #[error("stage type mismatch: expected {expected}, got {actual}")]
    TypeMismatch { expected: &'static str, actual: &'static str },

    #[error("pipeline was cancelled")]
    Cancelled,

    #[error("deadlock detected: no queue progress for {timeout_secs}s")]
    Deadlock { timeout_secs: u64 },

    #[error("stage '{name}' failed: {source}")]
    Stage {
        name: String,
        #[source]
        source: anyhow::Error,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_type_mismatch_message() {
        let err = PipelineError::TypeMismatch { expected: "u32", actual: "String" };
        assert_eq!(format!("{err}"), "stage type mismatch: expected u32, got String");
    }

    #[test]
    fn test_deadlock_message() {
        let err = PipelineError::Deadlock { timeout_secs: 10 };
        assert_eq!(format!("{err}"), "deadlock detected: no queue progress for 10s");
    }
}
