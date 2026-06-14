//! Compile-fail tests for the unified pipeline core. Each `.rs` file under
//! `tests/compile-fail/` is expected to fail to compile with the trait
//! bounds and type checks the framework promises.
//!
//! Per the Phase 0 design (line 530), trybuild covers type-mismatch
//! compile-fail invariants.

#[test]
fn pipeline_core_compile_fail() {
    let t = trybuild::TestCases::new();
    t.compile_fail("tests/compile-fail/*.rs");
}
