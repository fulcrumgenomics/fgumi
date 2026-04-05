//! Phase A types and step definitions for the alignment pipeline.
//!
//! Phase A covers Extract through `SortAccumulate`, with an external aligner in the middle.
//! [`PhaseAStep`] enumerates the work-stealing pool stages. Zipper and `SortAccumulate` are
//! dedicated threads (not pool stages). [`PhaseAScheduler`] implements drain-first priority
//! with sticky behavior for work-stealing scheduling.

/// Phase A pool stages (work-stealing).
///
/// Only Extract and `ToFastq` run in the pool. Zipper, stdout-reader, and
/// sort-accumulate are dedicated threads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhaseAStep {
    Extract = 0,
    Correct = 1,
    ToFastq = 2,
}

impl PhaseAStep {
    /// Number of Phase A pool steps.
    pub const NUM_STEPS: usize = 3;

    /// All steps in priority order (downstream-first for drain).
    pub const ALL: [PhaseAStep; 3] = [
        PhaseAStep::ToFastq, // drain q_extracted first
        PhaseAStep::Correct,
        PhaseAStep::Extract, // producer last
    ];

    /// Step name for logging.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Extract => "extract",
            Self::Correct => "correct",
            Self::ToFastq => "to_fastq",
        }
    }
}

/// Backpressure state for Phase A scheduling.
pub struct PhaseABackpressure {
    /// `q_extracted` is getting full — deprioritize Extract.
    pub extracted_high: bool,
    /// `ToFastq` input queue (`q_corrected` or `q_extracted`) is getting full — deprioritize Correct.
    pub tofastq_input_high: bool,
    /// Aligner batch gate — can `ToFastq` send?
    pub can_send_fastq: bool,
}

impl Default for PhaseABackpressure {
    fn default() -> Self {
        Self { extracted_high: false, tofastq_input_high: false, can_send_fastq: true }
    }
}

/// Phase A scheduler: drain-first priority with sticky behavior.
pub struct PhaseAScheduler {
    /// Current preferred step.
    current_step: PhaseAStep,
    /// Priority buffer.
    priority_buffer: [PhaseAStep; PhaseAStep::NUM_STEPS],
    /// Which steps are active (some may be skipped via `--start-from`).
    active_steps: [bool; PhaseAStep::NUM_STEPS],
}

impl PhaseAScheduler {
    /// Create a new scheduler with the given active steps.
    #[must_use]
    pub fn new(active_steps: [bool; PhaseAStep::NUM_STEPS]) -> Self {
        // Start on the highest-priority active step.
        let initial = PhaseAStep::ALL
            .iter()
            .find(|s| active_steps[**s as usize])
            .copied()
            .unwrap_or(PhaseAStep::ToFastq);

        Self {
            current_step: initial,
            priority_buffer: [PhaseAStep::Extract; PhaseAStep::NUM_STEPS],
            active_steps,
        }
    }

    /// Get priorities for this iteration.
    ///
    /// Returns a slice of steps in priority order. The current (sticky) step comes first,
    /// followed by remaining active steps in drain-first order. Steps gated by backpressure
    /// (`ToFastq` when batch gate is closed, Extract when extracted queue is high) are pushed
    /// to the end.
    pub fn get_priorities(&mut self, bp: &PhaseABackpressure) -> &[PhaseAStep] {
        let mut count = 0;

        // 1. Current step first (sticky on success).
        if self.active_steps[self.current_step as usize] {
            self.priority_buffer[count] = self.current_step;
            count += 1;
        }

        // 2. Drain-first order: ToFastq > Correct > Extract.
        //    Skip ToFastq if batch gate is closed; skip Extract if q_extracted is high.
        for &step in &PhaseAStep::ALL {
            if step == self.current_step {
                continue;
            }
            if !self.active_steps[step as usize] {
                continue;
            }
            if step == PhaseAStep::ToFastq && !bp.can_send_fastq {
                continue;
            }
            if step == PhaseAStep::Correct && bp.tofastq_input_high {
                continue;
            }
            if step == PhaseAStep::Extract && bp.extracted_high {
                continue;
            }

            self.priority_buffer[count] = step;
            count += 1;
        }

        // 3. Add backpressure-gated steps at the end (still available, just deprioritized).
        if self.active_steps[PhaseAStep::ToFastq as usize]
            && !bp.can_send_fastq
            && self.current_step != PhaseAStep::ToFastq
        {
            self.priority_buffer[count] = PhaseAStep::ToFastq;
            count += 1;
        }
        if self.active_steps[PhaseAStep::Correct as usize]
            && bp.tofastq_input_high
            && self.current_step != PhaseAStep::Correct
        {
            self.priority_buffer[count] = PhaseAStep::Correct;
            count += 1;
        }
        if self.active_steps[PhaseAStep::Extract as usize]
            && bp.extracted_high
            && self.current_step != PhaseAStep::Extract
        {
            self.priority_buffer[count] = PhaseAStep::Extract;
            count += 1;
        }

        &self.priority_buffer[..count]
    }

    /// Record the outcome of attempting a step.
    ///
    /// On success, the step becomes the new sticky step. On failure, the current sticky step
    /// is unchanged and the priority list will naturally move to the next option.
    pub fn record_outcome(&mut self, step: PhaseAStep, success: bool) {
        if success {
            self.current_step = step;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: all steps active.
    fn all_active() -> [bool; PhaseAStep::NUM_STEPS] {
        [true; PhaseAStep::NUM_STEPS]
    }

    #[test]
    fn test_default_priority_order() {
        let mut scheduler = PhaseAScheduler::new(all_active());
        let bp = PhaseABackpressure::default();
        let priorities = scheduler.get_priorities(&bp);

        // First should be ToFastq (highest priority, also initial sticky step).
        assert_eq!(priorities[0], PhaseAStep::ToFastq);
        // Should contain all 3 steps.
        assert_eq!(priorities.len(), PhaseAStep::NUM_STEPS);
        // Remaining should be in drain-first order.
        assert_eq!(priorities[1], PhaseAStep::Correct);
        assert_eq!(priorities[2], PhaseAStep::Extract);
    }

    #[test]
    fn test_sticky_on_success() {
        let mut scheduler = PhaseAScheduler::new(all_active());
        let bp = PhaseABackpressure::default();

        // Record success on Extract — it becomes the sticky step.
        scheduler.record_outcome(PhaseAStep::Extract, true);
        let priorities = scheduler.get_priorities(&bp);
        assert_eq!(priorities[0], PhaseAStep::Extract);
    }

    #[test]
    fn test_batch_gate_skips_to_fastq() {
        let mut scheduler = PhaseAScheduler::new(all_active());
        // Make Extract the sticky step so ToFastq is handled by the drain-first loop.
        scheduler.record_outcome(PhaseAStep::Extract, true);
        let bp = PhaseABackpressure {
            tofastq_input_high: false,
            can_send_fastq: false,
            ..Default::default()
        };
        let priorities = scheduler.get_priorities(&bp);

        // ToFastq should be last (deprioritized, not removed).
        assert_eq!(*priorities.last().unwrap(), PhaseAStep::ToFastq);

        // ToFastq should not appear in the middle.
        let middle = &priorities[1..priorities.len() - 1];
        assert!(!middle.contains(&PhaseAStep::ToFastq));
    }

    #[test]
    fn test_correct_backpressure_deprioritizes() {
        let mut scheduler = PhaseAScheduler::new(all_active());
        scheduler.record_outcome(PhaseAStep::Extract, true);
        let bp = PhaseABackpressure { tofastq_input_high: true, ..Default::default() };
        let priorities = scheduler.get_priorities(&bp);

        // Correct should be last (deprioritized because q_corrected is full).
        assert_eq!(*priorities.last().unwrap(), PhaseAStep::Correct);
    }

    #[test]
    fn test_inactive_steps_skipped() {
        let mut active = all_active();
        active[PhaseAStep::Correct as usize] = false;

        let mut scheduler = PhaseAScheduler::new(active);
        let bp = PhaseABackpressure::default();
        let priorities = scheduler.get_priorities(&bp);

        assert_eq!(priorities.len(), PhaseAStep::NUM_STEPS - 1);
        assert!(!priorities.contains(&PhaseAStep::Correct));
    }
}
