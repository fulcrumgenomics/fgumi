//! Deterministic, contention-free hashing for hot per-item maps.
//!
//! Re-exports the shared fixed-seed hasher so the seed table lives in exactly one place; see
//! [`fgumi_dna::hashing::deterministic_state`] for the full rationale.

pub(crate) use fgumi_dna::deterministic_state;
