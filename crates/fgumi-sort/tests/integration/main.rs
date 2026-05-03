//! Integration tests for the fgumi-sort crate.
//!
//! These tests validate end-to-end behavior that spans multiple modules
//! (sort engine, verifier, pool-integrated reader) and use shared SAM
//! fixtures, so they live under `tests/integration/` per the project
//! convention.

mod test_read_ahead;
mod test_verify;
