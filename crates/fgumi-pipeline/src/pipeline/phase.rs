//! Pipeline phase state for the two-phase architecture.
//!
//! The pipeline runs in two phases:
//! - **Phase A (Alignment):** Extract through `SortAccumulate`, with an external aligner.
//! - **Phase B (Processing):** `SortMerge` through Write (group, consensus, filter, compress).
//!
//! [`PhaseState`] provides an atomic flag that all workers read to determine which phase is active.

use std::sync::atomic::{AtomicU8, Ordering};

/// Pipeline phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum PipelinePhase {
    /// Phase A: alignment pipeline (Extract through `SortAccumulate`).
    Alignment = 0,
    /// Phase B: processing pipeline (`SortMerge` through Write).
    Processing = 1,
    /// Both phases complete.
    Done = 2,
}

/// Atomic phase flag read by all workers.
pub struct PhaseState {
    phase: AtomicU8,
}

impl PhaseState {
    /// Create a new phase state starting at the given phase.
    #[must_use]
    pub fn new(initial: PipelinePhase) -> Self {
        Self { phase: AtomicU8::new(initial as u8) }
    }

    /// Get the current phase.
    #[must_use]
    pub fn current(&self) -> PipelinePhase {
        match self.phase.load(Ordering::Acquire) {
            0 => PipelinePhase::Alignment,
            1 => PipelinePhase::Processing,
            _ => PipelinePhase::Done,
        }
    }

    /// Atomically transition from Alignment to Processing.
    /// Returns true if this call performed the transition.
    pub fn transition_to_processing(&self) -> bool {
        self.phase
            .compare_exchange(
                PipelinePhase::Alignment as u8,
                PipelinePhase::Processing as u8,
                Ordering::AcqRel,
                Ordering::Acquire,
            )
            .is_ok()
    }

    /// Transition to Done.
    pub fn transition_to_done(&self) {
        self.phase.store(PipelinePhase::Done as u8, Ordering::Release);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initial_phase() {
        let state = PhaseState::new(PipelinePhase::Alignment);
        assert_eq!(state.current(), PipelinePhase::Alignment);

        let state = PhaseState::new(PipelinePhase::Processing);
        assert_eq!(state.current(), PipelinePhase::Processing);

        let state = PhaseState::new(PipelinePhase::Done);
        assert_eq!(state.current(), PipelinePhase::Done);
    }

    #[test]
    fn test_transition_to_processing() {
        let state = PhaseState::new(PipelinePhase::Alignment);
        assert!(state.transition_to_processing());
        assert_eq!(state.current(), PipelinePhase::Processing);
    }

    #[test]
    fn test_transition_idempotent() {
        let state = PhaseState::new(PipelinePhase::Alignment);
        assert!(state.transition_to_processing());
        assert!(!state.transition_to_processing());
        assert_eq!(state.current(), PipelinePhase::Processing);
    }

    #[test]
    fn test_transition_from_wrong_phase() {
        let state = PhaseState::new(PipelinePhase::Processing);
        assert!(!state.transition_to_processing());
        assert_eq!(state.current(), PipelinePhase::Processing);
    }

    #[test]
    fn test_transition_to_done() {
        let state = PhaseState::new(PipelinePhase::Processing);
        state.transition_to_done();
        assert_eq!(state.current(), PipelinePhase::Done);
    }
}
