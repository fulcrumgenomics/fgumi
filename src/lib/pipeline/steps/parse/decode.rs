//! `DecodeRecords` mid-step. `Parallel + ByItemOrdinal`. Takes a
//! `DecompressedBlock` (record-aligned bytes from `FindBamBoundaries`),
//! parses the records, and computes a `GroupKey` for each in a single
//! pass — emits a `DecodedRecordBatch`.
//!
//! Mirrors the legacy pipeline's combined "Decode" step
//! (`pipeline/bam.rs:329-403`, `try_step_decode`): parse +
//! key extraction in one pass, parallelizable, so the downstream
//! Serial `GroupBam` step only does fast hash-and-name comparisons
//! under its mutex.
//!
//! ## Why parse + decode in one step (was two)
//!
//! An earlier shape split parse and decode into two `Parallel` steps
//! (`ParseBamRecords` → `DecodeRecords`) under the assumption that
//! some non-grouping callers would want the parse without the key
//! extraction. In practice every grouping command (`clip`, `group`,
//! `dedup`, `simplex`, `duplex`, `correct`, `filter`) does want the
//! keys, so the split forced one extra `Sequenced<T>`-wrap +
//! `ReorderStage` + queue traversal per batch with no payoff. Profiling
//! at scale (twist-umi 16 M records, threads=8) showed ~16M extra
//! atomic ops + ordinal allocations from the redundant edge —
//! collapsing them halves the per-record framework dispatch on the
//! upstream half of the chain.
//!
//! `ParseBamRecords` is still available as a standalone step for callers
//! who genuinely don't need group keys (e.g., `roundtrip` no-op chains,
//! integration tests that only mutate raw bytes).

use std::io;
use std::sync::Arc;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::parse::bam::parse_records;
use crate::pipeline::steps::types::{DecodedRecordBatch, DecompressedBlock, RecordBatch};
use fgumi_bam_io::{DecodedRecord, GroupKeyConfig, compute_group_key_from_raw, name_hash_key};

/// Fail closed on a raw BAM record body too short for the unchecked field
/// accessors used during group-key extraction.
///
/// The parse steps that feed this one ([`parse_records`] and the
/// `RecordBatch` walked by [`DecodeFromRecords`]) accept any record whose
/// `block_size` prefix is internally consistent — they do **not** enforce a
/// per-record minimum body size (see `parse_records`' contract). Both
/// [`name_hash_key`] and [`compute_group_key_from_raw`] then call
/// `fgumi_raw_bam::read_name`, which reads `raw[8]` (`l_read_name`, including
/// the trailing NUL) and slices `raw[32..32 + l_read_name - 1]` with **no**
/// bounds check; the full-key path additionally reads the 32-byte fixed header
/// (`flags`, `ref_id`, mate fields, …). A truncated or malformed record would
/// therefore panic on out-of-bounds indexing. We reject it here with
/// `InvalidData` before any key is computed so the pipeline surfaces a clean
/// error instead of unwinding a worker thread.
///
/// ## Minimum-length reasoning
///
/// - `raw.len() >= MIN_BAM_RECORD_LEN` (32) covers every fixed-offset field,
///   including the `raw[8]` read of `l_read_name` itself.
/// - A `>= MIN_BAM_RECORD_LEN` check alone is **not** sufficient: with
///   `l_read_name` as large as 255 the name slice
///   `raw[32..32 + l_read_name - 1]` ends well past a 32-byte record. So when
///   `l_read_name > 1` we additionally require
///   `raw.len() >= 32 + l_read_name - 1` (the exclusive end of that slice).
///   When `l_read_name <= 1`, `read_name` returns `&[]` without slicing, so the
///   32-byte header is enough.
fn validate_record_for_decode(raw: &[u8]) -> io::Result<()> {
    if raw.len() < fgumi_raw_bam::MIN_BAM_RECORD_LEN {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "DecodeRecords: BAM record too short for group-key extraction \
                 ({} byte(s) < {}-byte fixed header)",
                raw.len(),
                fgumi_raw_bam::MIN_BAM_RECORD_LEN,
            ),
        ));
    }
    // `raw[8]` is in bounds because `raw.len() >= MIN_BAM_RECORD_LEN` (>= 9).
    // `read_name` only slices the name region when `l_read_name > 1`; that
    // slice's exclusive end index is `32 + l_read_name - 1`, which must not run
    // past the record end.
    let l_read_name = raw[8] as usize;
    if l_read_name > 1 {
        let name_end = 32 + l_read_name - 1;
        if name_end > raw.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "DecodeRecords: BAM record read-name region runs past record end \
                     (l_read_name={l_read_name} needs {name_end} byte(s), record is {} byte(s))",
                    raw.len(),
                ),
            ));
        }
    }
    Ok(())
}

/// Cache the UMI tag's value position on a `DecodedRecord` during decode.
///
/// Looks up `umi_tag` in the record's aux data and, when present and Z-typed,
/// stores the record-relative offset of its value bytes (excluding the trailing
/// NUL) via [`DecodedRecord::set_cached_umi`]. The UMI-assignment step can then
/// slice the value through [`DecodedRecord::cached_umi`] instead of re-scanning
/// aux data once per template. When the tag is absent or not Z-typed the cache
/// is left in the `UMI_OFFSET_UNCACHED` state and downstream code falls back to
/// scanning aux data. Issue #334.
fn cache_umi_position(decoded: &mut DecodedRecord, umi_tag: [u8; 2]) {
    let raw = decoded.raw_bytes();
    let Some(aux_offset) = fgumi_raw_bam::aux_data_offset_from_record(raw) else {
        return;
    };
    let aux = &raw[aux_offset..];
    let Some((value_off_in_aux, value_len)) = fgumi_raw_bam::find_string_tag_position(aux, umi_tag)
    else {
        return;
    };
    let Ok(aux_offset_u32) = u32::try_from(aux_offset) else {
        return;
    };
    let Some(value_offset) = aux_offset_u32.checked_add(value_off_in_aux) else {
        return;
    };
    decoded.set_cached_umi(value_offset, value_len);
}

/// `Parallel + ByItemOrdinal` parse + decode. `DecompressedBlock →
/// DecodedRecordBatch`.
pub struct DecodeRecords {
    /// Shared across all worker clones — `GroupKeyConfig` holds an
    /// `Arc<LibraryIndex>` internally, so cloning the config is cheap.
    key_config: GroupKeyConfig,
    held: HeldSlot<Unpushed<DecodedRecordBatch>>,
    output_byte_limit: u64,
}

impl DecodeRecords {
    #[must_use]
    pub fn new(key_config: GroupKeyConfig, output_byte_limit: u64) -> Self {
        Self { key_config, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for DecodeRecords {
    fn clone(&self) -> Self {
        Self {
            key_config: self.key_config.clone(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

impl Step for DecodeRecords {
    type Input = DecompressedBlock;
    type Outputs = OrderedBytesSingle<DecodedRecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "DecodeRecords",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(block) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let DecompressedBlock { batch_serial, bytes } = block;

        // Single pass: parse the record-aligned bytes, then compute the
        // group key per record. Both halves are parallelizable per record
        // and don't share state, so the framework's `Parallel` kind covers
        // the whole pass.
        let records = parse_records(&bytes)?;
        let library_index: &Arc<_> = &self.key_config.library_index;
        let cell_tag = self.key_config.cell_tag;
        let name_hash_only = self.key_config.name_hash_only;
        let umi_tag = self.key_config.umi_tag;
        let decoded: Vec<DecodedRecord> = records
            .into_iter()
            .map(|raw| -> io::Result<DecodedRecord> {
                // Fail closed before the unchecked raw-field accessors run: a
                // truncated/malformed body would otherwise panic in
                // `read_name` (see `validate_record_for_decode`).
                validate_record_for_decode(raw.as_ref())?;
                // Queryname-grouping stages (e.g. correct) read only
                // `key.name_hash`; skip the CIGAR position walk + aux-tag pass.
                let key = if name_hash_only {
                    name_hash_key(raw.as_ref())
                } else {
                    compute_group_key_from_raw(raw.as_ref(), library_index, cell_tag)
                };
                let mut decoded = DecodedRecord::from_raw_bytes(raw, key);
                // Cache the UMI value position so the Group step's assignment
                // pass can slice it without re-scanning aux data (#334).
                if let Some(umi_tag) = umi_tag {
                    cache_umi_position(&mut decoded, umi_tag);
                }
                Ok(decoded)
            })
            .collect::<io::Result<Vec<_>>>()?;

        let out = DecodedRecordBatch::new(batch_serial, decoded);
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

/// `Parallel + ByItemOrdinal` group-key extractor. `RecordBatch →
/// DecodedRecordBatch`. Sibling of [`DecodeRecords`] that skips the
/// parse phase — useful when an earlier step (e.g. the `Sort` typed
/// step's output) already produced `RecordBatch`es and we need to
/// re-attach group keys without re-serializing through framed bytes.
///
/// Used by the runall fusion path:
///   `... → ParseBamRecords → Sort → DecodeFromRecords → GroupByPosition`
/// instead of the (no-sort) shorter
///   `... → DecodeRecords → GroupByPosition`.
pub struct DecodeFromRecords {
    key_config: GroupKeyConfig,
    held: HeldSlot<Unpushed<DecodedRecordBatch>>,
    output_byte_limit: u64,
}

impl DecodeFromRecords {
    #[must_use]
    pub fn new(key_config: GroupKeyConfig, output_byte_limit: u64) -> Self {
        Self { key_config, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for DecodeFromRecords {
    fn clone(&self) -> Self {
        Self {
            key_config: self.key_config.clone(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

impl Step for DecodeFromRecords {
    type Input = RecordBatch;
    type Outputs = OrderedBytesSingle<DecodedRecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "DecodeFromRecords",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(batch) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        let batch_serial = batch.batch_serial();

        let library_index: &Arc<_> = &self.key_config.library_index;
        let cell_tag = self.key_config.cell_tag;
        let name_hash_only = self.key_config.name_hash_only;
        let umi_tag = self.key_config.umi_tag;
        // `DecodedRecord` owns its bytes (via `RawRecord`), so we materialize
        // a heap-allocated copy here. The `RecordBatch`'s shared backing
        // buffer is dropped once this batch is consumed.
        let decoded: Vec<DecodedRecord> = batch
            .iter_record_bytes()
            .map(|bytes| -> io::Result<DecodedRecord> {
                // Fail closed before the unchecked raw-field accessors run: a
                // truncated/malformed body would otherwise panic in
                // `read_name` (see `validate_record_for_decode`).
                validate_record_for_decode(bytes)?;
                // Mirror `DecodeRecords::try_run`: queryname-grouping stages
                // (e.g. correct) read only `key.name_hash`, so skip the CIGAR
                // position walk + aux-tag pass when `name_hash_only` is set.
                let key = if name_hash_only {
                    name_hash_key(bytes)
                } else {
                    compute_group_key_from_raw(bytes, library_index, cell_tag)
                };
                let mut decoded = DecodedRecord::from_raw_bytes(
                    fgumi_raw_bam::RawRecord::from(bytes.to_vec()),
                    key,
                );
                // Cache the UMI value position so the Group step's assignment
                // pass can slice it without re-scanning aux data (#334).
                if let Some(umi_tag) = umi_tag {
                    cache_umi_position(&mut decoded, umi_tag);
                }
                Ok(decoded)
            })
            .collect::<io::Result<Vec<_>>>()?;

        let out = DecodedRecordBatch::new(batch_serial, decoded);
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = DecodeRecords::new(GroupKeyConfig::default(), 1024);
        let p = s.profile();
        assert_eq!(p.name, "DecodeRecords");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    #[test]
    fn from_records_profile_advertises_parallel_byordinal() {
        let s = DecodeFromRecords::new(GroupKeyConfig::default(), 1024);
        let p = s.profile();
        assert_eq!(p.name, "DecodeFromRecords");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }

    /// The `name_hash_only` decode path (used by queryname-grouping stages
    /// like `correct`) must produce the *same* `name_hash` the full
    /// `compute_group_key_from_raw` would, with every other `GroupKey` field
    /// left at its default. This is what makes the optimization output-safe:
    /// `GroupByQueryname` reads only `key.name_hash`.
    #[test]
    fn name_hash_only_key_matches_full_key_name_hash_and_zeros_rest() {
        use fgumi_bam_io::{GroupKey, LibraryIndex};
        use fgumi_raw_bam::SamBuilder;
        use fgumi_raw_bam::flags::{PAIRED, SECONDARY, UNMAPPED};
        use noodles::sam::alignment::record::data::field::Tag;

        let lib = LibraryIndex::default();
        let cb = Tag::from([b'C', b'B']);

        let make = |name: &[u8], flags: u16, mapped: bool| -> fgumi_raw_bam::RawRecord {
            let mut b = SamBuilder::new();
            b.read_name(name).flags(flags).sequence(b"ACGT").qualities(b"IIII");
            if mapped {
                b.ref_id(0).pos(100).cigar_ops(&[4u32 << 4]); // 4M
            }
            b.build()
        };

        let records = [
            make(b"q-mapped-paired", PAIRED, true),
            make(b"q-unmapped", UNMAPPED, false),
            make(b"q-secondary", SECONDARY, true),
            make(b"q", 0, true), // short name
        ];

        for raw in &records {
            let full = compute_group_key_from_raw(raw.as_ref(), &lib, Some(cb));
            let name_only = name_hash_key(raw.as_ref());
            assert_eq!(
                name_only.name_hash, full.name_hash,
                "name_hash_only must reproduce the full key's name_hash"
            );
            assert_eq!(
                name_only,
                GroupKey { name_hash: full.name_hash, ..GroupKey::default() },
                "name_hash_only must leave all non-name fields at default"
            );
        }
    }

    // ========================================================================
    // Fail-closed validation of undersized / malformed records
    // ========================================================================

    /// A record body shorter than the 32-byte BAM fixed header must be rejected
    /// with `InvalidData` rather than panicking in the unchecked field
    /// accessors. Mirrors the 8-byte records built by
    /// `parse_records_decodes_block_size_prefix` in `bam.rs` — exactly the
    /// shape `parse_records` accepts but the key extractors cannot decode.
    #[test]
    fn validate_record_for_decode_rejects_record_shorter_than_fixed_header() {
        let too_short = [0u8; 8];
        let err = validate_record_for_decode(&too_short).expect_err("must fail closed");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// A record long enough for the fixed header but whose declared
    /// `l_read_name` runs past the record end must also be rejected. This is
    /// the case a bare `len >= MIN_BAM_RECORD_LEN` check would miss: without
    /// the read-name bound, `read_name` would slice `raw[32..32 + 50 - 1]` on a
    /// 32-byte record and panic.
    #[test]
    fn validate_record_for_decode_rejects_truncated_read_name() {
        let mut rec = [0u8; fgumi_raw_bam::MIN_BAM_RECORD_LEN];
        rec[8] = 50; // l_read_name claims 49 name bytes + NUL, none of which exist
        let err = validate_record_for_decode(&rec).expect_err("must fail closed");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    /// The boundary cases that must be accepted: a 32-byte record with an empty
    /// name (`l_read_name <= 1`, where `read_name` never slices) and a real
    /// well-formed record built by `SamBuilder`. Proves the guard does not
    /// reject valid input.
    #[test]
    fn validate_record_for_decode_accepts_minimal_and_valid_records() {
        use fgumi_raw_bam::SamBuilder;

        // 32-byte record, l_read_name = 1 (just the NUL): read_name short-circuits.
        let mut minimal = [0u8; fgumi_raw_bam::MIN_BAM_RECORD_LEN];
        minimal[8] = 1;
        validate_record_for_decode(&minimal).expect("minimal empty-name record is valid");

        // A real record exercising read_name's slice path.
        let mut b = SamBuilder::new();
        b.read_name(b"read1").sequence(b"ACGT").qualities(&[30u8; 4]);
        let raw = b.build();
        validate_record_for_decode(raw.as_ref()).expect("well-formed record is valid");
    }

    /// End-to-end fail-closed proof: the same undersized record fed through the
    /// `DecodeRecords` map logic (the production `name_hash_only` and full-key
    /// branches) yields an `io::Error`, not a panic. Reproduces the map +
    /// `collect::<io::Result<Vec<_>>>()` shape used by `try_run`.
    #[test]
    fn decode_map_fails_closed_on_undersized_record_for_both_key_branches() {
        use fgumi_bam_io::{GroupKey, LibraryIndex};
        use noodles::sam::alignment::record::data::field::Tag;

        let lib = LibraryIndex::default();
        let cell_tag = Some(Tag::from([b'C', b'B']));
        // Below the fixed-header minimum — would panic in `read_name` if it
        // reached the key extractors.
        let undersized: &[u8] = &[0u8; 8];

        for name_hash_only in [true, false] {
            let result: io::Result<Vec<GroupKey>> = [undersized]
                .into_iter()
                .map(|raw| -> io::Result<GroupKey> {
                    validate_record_for_decode(raw)?;
                    let key = if name_hash_only {
                        name_hash_key(raw)
                    } else {
                        compute_group_key_from_raw(raw, &lib, cell_tag)
                    };
                    Ok(key)
                })
                .collect();
            let err = result.expect_err("undersized record must fail closed");
            assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        }
    }

    // ========================================================================
    // UMI position caching during Decode (issue #334)
    // ========================================================================

    #[test]
    fn cache_umi_position_records_value_offset_matching_fallback() {
        use crate::sam::SamTag;
        use fgumi_bam_io::GroupKey;
        use fgumi_raw_bam::{SamBuilder, flags};

        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT);
        b.add_string_tag(SamTag::RX, b"ACGTACGT");
        let raw = b.build();

        let mut decoded = DecodedRecord::from_raw_bytes(raw, GroupKey::default());
        cache_umi_position(&mut decoded, *SamTag::RX);

        // The cached slice must equal the canonical find_string_tag answer.
        let aux = fgumi_raw_bam::aux_data_slice(decoded.raw_bytes());
        let expected = fgumi_raw_bam::find_string_tag(aux, *SamTag::RX).expect("RX present");
        assert_eq!(decoded.cached_umi(), Some(expected));
        assert_eq!(decoded.cached_umi().unwrap(), b"ACGTACGT");
    }

    #[test]
    fn cache_umi_position_leaves_cache_unset_when_tag_missing() {
        use crate::sam::SamTag;
        use fgumi_bam_io::GroupKey;
        use fgumi_raw_bam::{SamBuilder, flags};

        // Record with no RX tag — caching must be a no-op.
        let mut b = SamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(flags::PAIRED | flags::FIRST_SEGMENT);
        let raw = b.build();

        let mut decoded = DecodedRecord::from_raw_bytes(raw, GroupKey::default());
        cache_umi_position(&mut decoded, *SamTag::RX);

        assert!(decoded.cached_umi().is_none());
        assert_eq!(decoded.cached_umi_position().0, DecodedRecord::UMI_OFFSET_UNCACHED);
    }
}
