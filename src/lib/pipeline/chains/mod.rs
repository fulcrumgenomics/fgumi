#![allow(rustdoc::broken_intra_doc_links)]

//! Unified chain-builder API for the typed-step pipeline framework.
//!
//! Single declarative entry point — [`build_for`] — that constructs a
//! ready-to-run [`BuiltPipeline`] from a [`ChainSpec`]. Both standalone
//! commands and `runall` route through here; the only behavioral
//! difference between them is how the spec is constructed.
//!
//! See `docs/design/refactor/2026-05-20-unified-chain-builder-design.md`
//! for the architectural rationale.

pub mod stage;
pub use stage::Stage;

pub mod source_spec;
pub use source_spec::SourceSpec;

pub mod sink_spec;
pub use sink_spec::SinkSpec;

pub mod finalize;
pub use finalize::{
    BuiltPipeline, FinalizeHook, PipelineStatsFinalizeHook, StageTimingFinalizeHook,
};

pub mod options_bag;
pub use options_bag::StageOptionsBag;

pub mod spec;
pub use spec::ChainSpec;

pub mod validate;

pub mod build;
pub use build::build_for;

pub mod builder;
pub use builder::ChainBuilder;

pub mod build_helpers;
pub use build_helpers::build_pipeline_config_for_chain;

pub mod commands;
