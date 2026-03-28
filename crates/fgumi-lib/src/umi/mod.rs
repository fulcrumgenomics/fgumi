//! UMI (Unique Molecular Identifier) utilities
//!
//! This module re-exports functionality from the `fgumi-umi` crate.

// Re-export the assigner module and top-level items from fgumi-umi
pub use fgumi_umi::{
    MoleculeId, TagInfo, TagSets, UmiValidation, assigner, extract_mi_base, validate_umi,
};

// Re-export commonly used items from the assigner module for convenience
pub use fgumi_umi::assigner::{
    AdjacencyUmiAssigner, IdentityUmiAssigner, PairedUmiAssigner, SimpleErrorUmiAssigner, Strategy,
    Umi, UmiAssigner,
};

pub mod parallel_assigner;
