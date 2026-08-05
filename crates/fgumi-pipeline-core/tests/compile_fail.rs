//! Compile-fail tests for the unified pipeline core. Each `.rs` file under
//! `tests/compile-fail/` is expected to fail to compile with the trait
//! bounds and type checks the framework promises.
//!
//! trybuild covers the type-mismatch compile-fail invariants: that a chain
//! link's `Step::Input` must match its producer's item type
//! (`tests/compile-fail/chain_input_type_mismatch.rs`), and that
//! `OrderedBytesSingle<T>` requires both `HeapSize` and `Ordered`
//! (`ordered_bytes_single_requires_{heapsize,ordered}.rs`).
//!
//! ## The `.stderr` fixtures are whole-diagnostic, not assertions
//!
//! trybuild compares each fixture to rustc's normalized stderr for **equality**;
//! it has no wildcard or subset matching, so a `.stderr` file records the entire
//! diagnostic — including rustc's `the following other types implement trait
//! HeapSize` / `the trait Ordered is implemented for Sequenced<T>` help blocks.
//! Those blocks are an inventory of the crate's impls, not part of the bound
//! under test, so adding a `HeapSize` or `Ordered` impl anywhere rewrites them
//! and fails this test for a reason unrelated to `OrderedBytesSingle`.
//!
//! They cannot simply be deleted from the fixtures — that is a mismatch like any
//! other and fails immediately. So when a diff here is confined to those help
//! blocks (or to rustc's wording), it is a **re-bless**, not a contract
//! regression: check that the `error[E0277]` lines and the `required by a bound
//! in OrderedBytesSingle` notes still say what they should, then regenerate with
//! `TRYBUILD=overwrite cargo nextest run -p fgumi-pipeline-core --test compile_fail`.
//! A diff that touches an `E0277` line or a `required by a bound` note is the
//! real signal and must not be blessed away.
//!
//! That re-bless runs through a bare `cargo nextest run`, not the repo's
//! `cargo ci-test` alias: the alias hardcodes `--workspace`, which wins over the
//! `-p` filter and also selects `fgumi-cli-macros`' own `compile_fail` binary, so
//! blessing through it would rewrite that crate's fixtures in the same pass.

/// The number of fixtures `tests/compile-fail/` is expected to hold. Asserted
/// before handing the glob to trybuild because `compile_fail` **passes when the
/// glob matches nothing** — so renaming or moving the fixture directory would
/// silently retire every compile-time contract above while CI stayed green.
/// Raise this when adding a fixture.
const EXPECTED_FIXTURES: usize = 3;

#[test]
fn pipeline_core_compile_fail() {
    let fixtures = std::fs::read_dir("tests/compile-fail")
        .expect("tests/compile-fail must exist")
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "rs"))
        .count();
    assert!(
        fixtures >= EXPECTED_FIXTURES,
        "expected at least {EXPECTED_FIXTURES} compile-fail fixtures, found {fixtures}; \
         a zero-match glob makes `compile_fail` a silent no-op"
    );
    let t = trybuild::TestCases::new();
    t.compile_fail("tests/compile-fail/*.rs");
}
