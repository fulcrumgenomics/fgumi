//! Compile-fail coverage for the build-time diagnostics `multi_options` documents.
//!
//! The macro reports every unsupported input by returning a spanned `syn::Error`
//! and emitting it via `into_compile_error()`, so these surface as ordinary
//! compilation errors — they cannot be exercised with `#[should_panic]`, because
//! the macro runs while the test crate is being compiled, not while it runs.
//! `trybuild` compiles each case and diffs the emitted error against its
//! committed `.stderr`, which pins both the message and the span it points at.

#[test]
fn documented_invalid_inputs_fail_to_compile() {
    trybuild::TestCases::new().compile_fail("tests/ui/*.rs");
}
