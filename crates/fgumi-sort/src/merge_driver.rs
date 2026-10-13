//! The k-way merge driver over [`SortMergeSlot`]s and in-memory chunks: the
//! slot parser, [`MergeDriver`] and its type-erased dispatch.

use crate::inline::{ChunkRef, InMemoryChunk};
use crate::keys::RawSortKey;
use crate::merge_slots::SortMergeSlot;
use crate::spill_block_reader::read_exact_or_eof;
use anyhow::{Context as _, Result};
use fgumi_bam_io::ProgressTracker;
use log::info;
use std::sync::Arc;

/// Open `path` as a single [`SortMergeSlot`] with the given `file_id`.
///
/// The spill codec is auto-detected from the file magic: for zstd (`ZSP1`) the
/// frames start after the 4-byte magic; for BGZF the `1f 8b` is part of the
/// first block, so the frames start at byte 0. The slot carries the file as a
/// positional handle with that body range and allocates no read buffer;
/// `SortSpillDecompress` reads `slot.codec` to decompress accordingly.
///
/// `file_id` is the slot's stable merge-order identifier — `SortMerge` orders
/// sources by it so the `LoserTree` tie-break for equal sort keys matches the
/// legacy chunk-files order. The spill writer passes a run's **logical spill
/// index** here so the tie-break is independent of which worker happened to
/// compress which block.
///
/// # Errors
///
/// Returns an error if the file cannot be opened, or its magic or length
/// cannot be read.
pub fn open_spill_slot(path: &std::path::Path, file_id: u32) -> Result<Arc<SortMergeSlot>> {
    open_spill_slot_as(path, file_id, file_id)
}

/// Open a merge input — a BAM file or stream behind `source` — as a
/// body-keyed slot: its BGZF stream read from byte 0, the `header_len`
/// decompressed header bytes skipped, keys extracted from the record bodies.
/// `len` is the file's length, or `None` for a stream (a sequential slot,
/// whose EOF is its first short read).
///
/// # Errors
///
/// `InvalidInput` when `file_id` is in the consolidation slot-id range.
pub fn open_bam_slot(
    source: Arc<dyn fgumi_bam_io::pread::PositionalSource>,
    file_id: u32,
    header_len: u64,
    len: Option<u64>,
) -> Result<Arc<SortMergeSlot>> {
    if file_id >= crate::merge_slots::CONSOLIDATION_SLOT_IDS_START {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("merge input {file_id} is in the consolidation slot-id range"),
        )
        .into());
    }
    let codec = crate::codec::SpillCodec::Bgzf;
    let slot = match len {
        Some(len) => SortMergeSlot::new(file_id, source, 0, len, codec),
        None => SortMergeSlot::new_sequential(file_id, source, codec),
    };
    Ok(Arc::new(slot.with_skip_bytes(header_len).body_keyed()))
}

/// [`open_spill_slot`] under supply identity `slot_id`
/// ([`SortMergeSlot::slot_id`]): a consolidation input, opened under a fresh
/// id so it never shares a supply-map key with another slot of its `file_id`.
///
/// # Errors
///
/// As [`open_spill_slot`]; also `InvalidInput` when `file_id` is at or above
/// [`crate::merge_slots::CONSOLIDATION_SLOT_IDS_START`] (the minted range).
pub fn open_spill_slot_as(
    path: &std::path::Path,
    file_id: u32,
    slot_id: u32,
) -> Result<Arc<SortMergeSlot>> {
    if file_id >= crate::merge_slots::CONSOLIDATION_SLOT_IDS_START {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "spill run file id {file_id} is in the consolidation slot-id range (>= {})",
                crate::merge_slots::CONSOLIDATION_SLOT_IDS_START
            ),
        )
        .into());
    }
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("failed to open spill chunk {}", path.display()))?;
    let mut magic = [0u8; 4];
    let filled = read_exact_or_eof(&mut file, &mut magic)
        .with_context(|| format!("failed to read spill chunk magic {}", path.display()))?;
    let codec = if filled {
        crate::codec::SpillCodec::from_magic(&magic).unwrap_or(crate::codec::SpillCodec::Bgzf)
    } else {
        crate::codec::SpillCodec::Bgzf
    };
    let body_start = crate::spill_block::spill_magic(codec).len() as u64;
    let len = file
        .metadata()
        .with_context(|| format!("failed to stat spill chunk {}", path.display()))?
        .len();
    Ok(Arc::new(
        SortMergeSlot::new(file_id, Arc::new(file), body_start, len, codec).with_slot_id(slot_id),
    ))
}

// ============================================================================
// The merge driver (the chain's `SortMerge` and `ConsolidateRuns` steps)
// ============================================================================
//
// The slot-backed, *non-blocking* merge driver (`MergeDriver::from_slots`). It
// reads records out of `SortMergeSlot` queues, which the `SortSpillDecompress`
// step has *already* filled with decompressed `Vec<u8>` blocks — so this path
// is **codec-agnostic** (it never touches disk or the spill codec).

/// Result of a non-blocking source read. `WouldBlock` is only ever produced by
/// `Slot` sources whose decompressed queue is momentarily empty and not yet at
/// EOF; an in-memory chunk source always returns `Ready`.
pub(crate) enum TryRead<T> {
    /// A record was read (`Some`) or the source is at clean EOF (`None`).
    Ready(Option<T>),
    /// The source has no bytes available right now but is not at EOF; the
    /// caller should yield and retry on a later dispatch.
    WouldBlock,
}

/// Where a merge source's current record lives, so the driver can emit it by
/// reference instead of copying it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RecLoc {
    /// `[start, end)` of the source's `SlotParserState::current_buf`: the
    /// record lies wholly inside the current block.
    Block { start: usize, end: usize },
    /// Record `idx` of the source's shared residual chunk (`MemoryShared`).
    Chunk { idx: usize },
    /// The leaf's scratch buffer: a record that spanned blocks.
    Scratch,
}

/// Per-slot byte-stream parser state for the non-blocking merge reader.
///
/// As decompressed blocks are popped from a `SortMergeSlot`, their bytes are
/// stashed in `current_buf` and consumed left-to-right by the record parser,
/// with resumable cross-block framing state.
struct SlotParserState {
    /// Current decompressed block being consumed.
    current_buf: Vec<u8>,
    /// Read position within `current_buf`.
    current_pos: usize,
    /// Resumable framing state for [`slot_try_next_record`]. `Some` while a
    /// record spans a block boundary and the next block was not yet available,
    /// so the read returned `WouldBlock` mid-record; the next call resumes from
    /// here.
    pending: Option<PendingRecord>,
    /// Decompressed bytes still to skip before the first record (a BAM
    /// input's header; [`SortMergeSlot::skip_bytes`]).
    skip_remaining: u64,
    /// The slot's records carry no key prefix: keys come from their bodies
    /// ([`SortMergeSlot::is_body_keyed`]).
    body_keyed: bool,
    /// Applied to every record body before its key is extracted (the merge's
    /// read-tag rewrite), with the slot's `file_id`; the rewritten record is
    /// what the merge emits.
    rewrite: Option<RecordRewrite>,
}

/// A rewrite applied to each record of a merge source before its key is
/// extracted: `(file_id, record)`, editing the raw record in place.
pub type RecordRewrite = Arc<dyn Fn(u32, &mut Vec<u8>) + Send + Sync>;

/// Partially-read record, retained across a `WouldBlock` so the non-blocking
/// slot reader can resume. The wire layout is `[key:ksize][len:4][body:len]`
/// for non-embedded keys (`ksize = K::SERIALIZED_SIZE`) and `[len:4][body:len]`
/// for embedded keys (`ksize = 0`, key extracted from the body). We collect the
/// key prefix (if any), then the length prefix, then the body.
#[derive(Default)]
struct PendingRecord {
    /// Expected serialized key size (`0` for embedded keys).
    key_size: usize,
    /// Key-prefix bytes collected so far (length grows to `key_size`).
    key_buf: Vec<u8>,
    /// Length-prefix bytes collected so far (`0..=4`).
    len_buf: [u8; 4],
    /// Count of `len_buf` bytes filled.
    len_got: usize,
    /// Total body length, known once `len_got == 4`. The body accumulates
    /// directly into the caller's scratch buffer; progress is `out.len()`.
    body_len: Option<usize>,
}

impl SlotParserState {
    fn new() -> Self {
        Self {
            current_buf: Vec::new(),
            current_pos: 0,
            pending: None,
            skip_remaining: 0,
            body_keyed: false,
            rewrite: None,
        }
    }

    /// The parser for `slot`: its header skip and key layout.
    fn for_slot(slot: &crate::merge_slots::SortMergeSlot) -> Self {
        Self { skip_remaining: slot.skip_bytes(), body_keyed: slot.is_body_keyed(), ..Self::new() }
    }

    fn remaining(&self) -> usize {
        self.current_buf.len() - self.current_pos
    }
}

/// Outcome of a non-blocking attempt to pull the next decompressed block from a
/// slot into `parser.current_buf`.
enum BlockLoad {
    /// A block was moved into `parser.current_buf`.
    Loaded,
    /// The slot is fully delivered (queue empty AND `queue_eof`).
    Eof,
    /// The queue is empty but the producer is still feeding this slot.
    WouldBlock,
}

/// Pop one decompressed block (FIFO) into `parser`, surfacing a decompression
/// error in preference to EOF, reporting clean EOF, or `WouldBlock` while the
/// producer is still active. Never parks.
fn slot_try_load_block(
    slot: &Arc<crate::merge_slots::SortMergeSlot>,
    parser: &mut SlotParserState,
) -> Result<BlockLoad> {
    use std::sync::atomic::Ordering;

    let mut guard = slot.decompressed.lock().expect("slot.decompressed mutex poisoned");
    if let Some(data) = slot.pop_locked(&mut guard) {
        drop(guard);
        parser.current_buf = data;
        parser.current_pos = 0;
        return Ok(BlockLoad::Loaded);
    }
    // Queue empty. Both atomics are visible under the lock release-acquire
    // chain (producer stores them while holding this same mutex). Surface
    // error in preference to EOF.
    if slot.decomp_error.load(Ordering::Acquire) {
        drop(guard);
        anyhow::bail!("spill decompression error on slot {}", slot.file_id);
    }
    if slot.queue_eof.load(Ordering::Acquire) {
        return Ok(BlockLoad::Eof);
    }
    Ok(BlockLoad::WouldBlock)
}

/// Serialized key size on the wire for a slot source: `0` for embedded keys
/// (key extracted from the body), else `K::SERIALIZED_SIZE`. Variable-length
/// non-embedded keys are unsupported on the slot path (none exist in
/// production — every slot-path key is embedded).
pub(crate) fn slot_key_size<K: RawSortKey>() -> Result<usize> {
    if K::EMBEDDED_IN_RECORD {
        Ok(0)
    } else {
        K::SERIALIZED_SIZE.ok_or_else(|| {
            anyhow::anyhow!("non-embedded slot key must have a fixed serialized size")
        })
    }
}

/// Parse a sort key from `bytes`: extracted from the body for embedded keys, or
/// deserialized from the `key_size`-byte prefix for non-embedded keys.
#[cfg(test)]
pub(crate) fn slot_parse_key<K: RawSortKey + Default + 'static>(
    ctx: &K::MergeCtx,
    hint: &mut K::MergeHint,
    key_bytes: &[u8],
    body: &[u8],
) -> Result<K> {
    slot_parse_key_in::<K>(ctx, hint, key_bytes, body, false)
}

/// `slot_parse_key` for a slot whose records may be body-keyed: a
/// body-keyed slot's key is extracted from the body whatever `K`'s spill
/// layout ([`RawSortKey::extract_from_record_in`] with the merge's context).
pub(crate) fn slot_parse_key_in<K: RawSortKey + Default + 'static>(
    ctx: &K::MergeCtx,
    hint: &mut K::MergeHint,
    key_bytes: &[u8],
    body: &[u8],
    body_keyed: bool,
) -> Result<K> {
    if body_keyed {
        // A body-keyed slot's records come from a user's BAM, not from a spill
        // run fgumi wrote: hold them to the BAM layout before any extractor
        // indexes into them, as the sort's ingest does, so a short or
        // inconsistent record is an error naming the input (`MergeDriver`
        // adds the input and record), not a panic.
        check_bam_record_layout(body)?;
        return Ok(K::extract_from_record_in(ctx, hint, body));
    }
    if K::EMBEDDED_IN_RECORD {
        // An embedded key is read out of the body, so an empty body cannot
        // carry one. Every extractor indexes straight into the slice
        // (`RawCoordinateKey` reads `ref_id` off the front), so without this an
        // empty body panics mid-merge instead of surfacing a corrupt-slot
        // error — and it is reachable: a corrupt length prefix of 0 makes the
        // body loop break immediately with `out` still empty.
        //
        // This is a floor, not a full bounds check. The minimum body an
        // extractor needs is key-specific (8 bytes here, 16 for the coordinate
        // key, more for queryname), and `RawSortKey` does not expose it, so a
        // corrupt *small but non-zero* length can still panic inside an
        // extractor. Closing that needs a per-key minimum on the trait.
        anyhow::ensure!(
            !body.is_empty(),
            "record body is empty and cannot carry an embedded sort key (corrupt slot)",
        );
        Ok(K::extract_from_record_in(ctx, hint, body))
    } else {
        let mut r = key_bytes;
        K::read_from(&mut r).map_err(|e| anyhow::anyhow!("malformed slot key: {e}"))
    }
}

/// A user BAM record's layout holds: at least the 32-byte fixed core, and the
/// read name, CIGAR, sequence and qualities its core declares fit inside it.
fn check_bam_record_layout(body: &[u8]) -> Result<()> {
    anyhow::ensure!(
        body.len() >= fgumi_raw_bam::MIN_BAM_RECORD_LEN,
        "block_size {} is shorter than a BAM record ({} bytes)",
        body.len(),
        fgumi_raw_bam::MIN_BAM_RECORD_LEN,
    );
    let needed = fgumi_raw_bam::aux_data_offset_from_record(body).unwrap_or(usize::MAX);
    anyhow::ensure!(
        needed <= body.len(),
        "block_size {} is shorter than the {needed} bytes its read name, CIGAR, sequence and \
         qualities need",
        body.len(),
    );
    Ok(())
}

/// Non-blocking, resumable read of the next record from a slot.
///
/// Returns `Ready(Some((key, loc)))` with `loc` naming where the record body
/// lives: [`RecLoc::Block`] for a record wholly inside the current block (left
/// in place, not copied), or [`RecLoc::Scratch`] for one assembled in `out`
/// because it spans blocks. Returns `Ready(None)` at clean source EOF, or
/// `WouldBlock` when a needed block is not yet available (progress is saved in
/// `parser.pending` so the next call resumes).
///
/// # Errors
///
/// Truncation (EOF mid-record), a malformed key, or a producer-side
/// decompression error.
fn slot_try_next_record<K: RawSortKey + Default + 'static>(
    slot: &Arc<crate::merge_slots::SortMergeSlot>,
    parser: &mut SlotParserState,
    out: &mut Vec<u8>,
    ctx: &K::MergeCtx,
    hint: &mut K::MergeHint,
) -> Result<TryRead<(K, RecLoc)>> {
    let key_size = if parser.body_keyed { 0 } else { slot_key_size::<K>()? };

    if parser.pending.is_none() {
        // A BAM input's header precedes its first record: skip its
        // decompressed bytes, across as many blocks as it takes.
        while parser.skip_remaining > 0 {
            if parser.remaining() == 0 {
                match slot_try_load_block(slot, parser)? {
                    BlockLoad::Loaded => {}
                    BlockLoad::Eof => {
                        return Err(truncated(format!(
                            "slot {} ends inside its {}-byte header",
                            slot.file_id,
                            slot.skip_bytes()
                        )));
                    }
                    BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
                }
            }
            let take = usize::try_from(parser.skip_remaining)
                .unwrap_or(usize::MAX)
                .min(parser.remaining());
            parser.current_pos += take;
            parser.skip_remaining -= take as u64;
        }
        // Ensure at least one byte is available, or detect clean EOF at a
        // record boundary.
        if parser.remaining() == 0 {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => return Ok(TryRead::Ready(None)),
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        // Fast path: the entire record lives in the current block → parse in
        // place and leave the body there; the driver emits it by reference.
        let header = key_size + 4;
        if parser.remaining() >= header {
            let p = parser.current_pos;
            let len = u32::from_le_bytes(
                parser.current_buf[p + key_size..p + header].try_into().expect("4-byte slice"),
            ) as usize;
            if parser.remaining() >= header + len {
                let start = p + header;
                let end = start + len;
                if let Some(rewrite) = &parser.rewrite {
                    // A rewritten record may change length, so it is emitted
                    // from the scratch, not the block.
                    out.clear();
                    out.extend_from_slice(&parser.current_buf[start..end]);
                    parser.current_pos = end;
                    // The rewriter walks the aux data, so a user's record is
                    // held to the BAM layout before it, not only after.
                    if parser.body_keyed {
                        check_bam_record_layout(out)?;
                    }
                    rewrite(slot.file_id, out);
                    let key = slot_parse_key_in::<K>(ctx, hint, &[], out, parser.body_keyed)?;
                    return Ok(TryRead::Ready(Some((key, RecLoc::Scratch))));
                }
                let key = slot_parse_key_in::<K>(
                    ctx,
                    hint,
                    &parser.current_buf[p..p + key_size],
                    &parser.current_buf[start..end],
                    parser.body_keyed,
                )?;
                parser.current_pos = end;
                return Ok(TryRead::Ready(Some((key, RecLoc::Block { start, end }))));
            }
        }
        // Slow path: record spans a block boundary (or fewer than `header`
        // bytes remain). Begin a resumable pending record.
        parser.pending = Some(PendingRecord { key_size, ..PendingRecord::default() });
    }

    Ok(match slot_collect_pending::<K>(slot, parser, out, ctx, hint)? {
        TryRead::Ready(key) => TryRead::Ready(key.map(|k| (k, RecLoc::Scratch))),
        TryRead::WouldBlock => TryRead::WouldBlock,
    })
}

/// A slot that ends inside a record: an `UnexpectedEof` I/O error carrying
/// `msg`, so a caller can tell truncation from other merge failures.
fn truncated(msg: String) -> anyhow::Error {
    std::io::Error::new(std::io::ErrorKind::UnexpectedEof, msg).into()
}

/// Resumable continuation of [`slot_try_next_record`]'s slow path. Collects the
/// key prefix (if any), the 4-byte length prefix, then the body across as many
/// blocks as needed, returning `WouldBlock` (state retained in
/// `parser.pending`) whenever the next block isn't ready yet.
#[allow(clippy::too_many_lines)]
fn slot_collect_pending<K: RawSortKey + Default + 'static>(
    slot: &Arc<crate::merge_slots::SortMergeSlot>,
    parser: &mut SlotParserState,
    out: &mut Vec<u8>,
    ctx: &K::MergeCtx,
    hint: &mut K::MergeHint,
) -> Result<TryRead<K>> {
    // Stage 0: key prefix (non-embedded keys only; `key_size == 0` is a no-op).
    loop {
        let (key_size, key_have) = {
            let pending = parser.pending.as_ref().expect("pending set");
            (pending.key_size, pending.key_buf.len())
        };
        let need = key_size - key_have;
        if need == 0 {
            break;
        }
        if parser.current_pos == parser.current_buf.len() {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => {
                    return Err(truncated(format!(
                        "truncated record key in slot {} ({key_have} of {key_size} key bytes at \
                         EOF)",
                        slot.file_id,
                    )));
                }
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        let avail = parser.current_buf.len() - parser.current_pos;
        let take = need.min(avail);
        let src = &parser.current_buf[parser.current_pos..parser.current_pos + take];
        parser.pending.as_mut().expect("pending set").key_buf.extend_from_slice(src);
        parser.current_pos += take;
    }

    // Stage 1: length prefix.
    loop {
        let len_got = parser.pending.as_ref().expect("pending set").len_got;
        if len_got == 4 {
            break;
        }
        if parser.current_pos == parser.current_buf.len() {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => {
                    return Err(truncated(format!(
                        "truncated record length in slot {} ({len_got} of 4 length bytes at EOF)",
                        slot.file_id,
                    )));
                }
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        let avail = parser.current_buf.len() - parser.current_pos;
        let take = (4 - len_got).min(avail);
        let src = &parser.current_buf[parser.current_pos..parser.current_pos + take];
        let pending = parser.pending.as_mut().expect("pending set");
        pending.len_buf[len_got..len_got + take].copy_from_slice(src);
        pending.len_got += take;
        parser.current_pos += take;
    }

    // Compute body length once and reset `out` for accumulation.
    {
        let pending = parser.pending.as_mut().expect("pending set");
        if pending.body_len.is_none() {
            let len = u32::from_le_bytes(pending.len_buf) as usize;
            // The prefix is untrusted — it comes straight off a decompressed
            // spill block — and the next statement hands it to `reserve`, which
            // answers an allocation failure by aborting the process rather than
            // returning an error. Bound it first, so a corrupt slot surfaces
            // `Err` like the `decomp_error` and truncation arms do. The ceiling
            // is the one a single record already cannot exceed: a record must
            // fit within one sort segment (see `SORT_SEGMENT_SIZE`).
            anyhow::ensure!(
                len <= crate::inline::SORT_SEGMENT_SIZE,
                "implausible record length {len} in slot {} (exceeds the \
                 {}-byte segment a single record must fit in)",
                slot.file_id,
                crate::inline::SORT_SEGMENT_SIZE,
            );
            pending.body_len = Some(len);
            out.clear();
            out.reserve(len);
        }
    }

    // Stage 2: body.
    loop {
        let body_len = parser.pending.as_ref().expect("pending set").body_len.expect("len known");
        if out.len() >= body_len {
            break;
        }
        if parser.current_pos == parser.current_buf.len() {
            match slot_try_load_block(slot, parser)? {
                BlockLoad::Loaded => {}
                BlockLoad::Eof => {
                    return Err(truncated(format!(
                        "truncated record body in slot {} ({} of {body_len} bytes at EOF)",
                        slot.file_id,
                        out.len(),
                    )));
                }
                BlockLoad::WouldBlock => return Ok(TryRead::WouldBlock),
            }
        }
        let avail = parser.current_buf.len() - parser.current_pos;
        let take = (body_len - out.len()).min(avail);
        out.extend_from_slice(&parser.current_buf[parser.current_pos..parser.current_pos + take]);
        parser.current_pos += take;
    }

    if let Some(rewrite) = &parser.rewrite {
        // As on the fast path: check the layout before the rewriter reads it.
        if parser.body_keyed {
            check_bam_record_layout(out)?;
        }
        rewrite(slot.file_id, out);
    }
    let key = {
        let pending = parser.pending.as_ref().expect("pending set");
        slot_parse_key_in::<K>(ctx, hint, &pending.key_buf, out, parser.body_keyed)?
    };
    parser.pending = None;
    Ok(TryRead::Ready(Some(key)))
}

/// Merge source for the slot-backed `MergeDriver`: either a decompressing slot
/// (read non-blocking) or an already-sorted in-memory residual chunk.
enum SlotMergeSource<K: RawSortKey + Default + 'static, R: ChunkRef<K> = K> {
    /// Pipeline-integrated source — `SortSpillDecompress` fills
    /// `slot.decompressed` with decompressed blocks; this driver reads them
    /// in-order via the per-source [`SlotParserState`].
    /// `hint` is the key's per-source merge hint, carried from record to
    /// record (see [`RawSortKey::MergeHint`]).
    Slot {
        slot: Arc<crate::merge_slots::SortMergeSlot>,
        parser: SlotParserState,
        hint: K::MergeHint,
    },
    /// In-memory sorted residual from an inline-buffer (coordinate / template)
    /// sort, sharing an `Arc<SegmentedBuf>` — zero per-record allocation. A
    /// read names the record's index; the driver emits its bytes straight from
    /// the shared buffer.
    MemoryShared { chunk: InMemoryChunk<R>, idx: usize },
}

impl<K: RawSortKey + Default + 'static, R: ChunkRef<K>> SlotMergeSource<K, R> {
    /// Non-blocking read of the next record, returning its sort key and where
    /// its bytes live ([`RecLoc`]). Only a record that spans blocks is placed
    /// in `scratch`. `Slot` sources may return
    /// `WouldBlock`; memory sources always return `Ready`.
    fn try_next_record(
        &mut self,
        scratch: &mut Vec<u8>,
        ctx: &K::MergeCtx,
    ) -> Result<TryRead<(K, RecLoc)>> {
        match self {
            SlotMergeSource::Slot { slot, parser, hint } => {
                slot_try_next_record::<K>(slot, parser, scratch, ctx, hint)
            }
            SlotMergeSource::MemoryShared { chunk, idx } => {
                if *idx < chunk.len() {
                    // The bytes stay in the shared `Arc<SegmentedBuf>`; the
                    // driver emits them from there by index.
                    let key = chunk.key_at(*idx).merge_key(chunk.record_bytes(*idx));
                    let loc = RecLoc::Chunk { idx: *idx };
                    *idx += 1;
                    Ok(TryRead::Ready(Some((key, loc))))
                } else {
                    Ok(TryRead::Ready(None))
                }
            }
        }
    }
}

/// One step of a non-blocking k-way merge.
pub enum MergeStep<'a> {
    /// The next merged record. The bytes borrow the winning source's current
    /// block, its shared residual chunk, or the driver's scratch buffer (a
    /// record that spanned blocks), and are valid only until the next
    /// `try_step` call — the caller must copy them out before stepping again.
    /// The borrow checker enforces this, which is what makes the
    /// deferred-refill protocol sound.
    Produced(&'a [u8]),
    /// A source's decompressed queue is momentarily empty (and not at EOF), or
    /// a source hasn't been primed yet. No record is available right now; the
    /// caller yields and retries on a later dispatch. The slot named by
    /// [`MergeDriverDyn::stalled_slot`] is the one to await: its producer wakes
    /// the merge when the awaited block, EOF or failure lands
    /// ([`crate::MergeDemand`]).
    Stalled,
    /// The merge is exhausted.
    Done,
}

impl std::fmt::Debug for MergeStep<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Produced(bytes) => write!(f, "Produced({} bytes)", bytes.len()),
            Self::Stalled => f.write_str("Stalled"),
            Self::Done => f.write_str("Done"),
        }
    }
}

/// One step of a non-blocking k-way merge that frames the merged record into a
/// caller's buffer ([`MergeDriverDyn::try_step_framed`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MergeStepFramed {
    /// The next merged record was appended to the buffer.
    Produced,
    /// A source is momentarily empty (see [`MergeStep::Stalled`]).
    Stalled,
    /// The merge is exhausted.
    Done,
}

/// Object-safe view over `MergeDriver<K>` so the pipeline `SortMerge` step can
/// hold a single concrete type regardless of sort order.
pub trait MergeDriverDyn: Send {
    /// Produce the next merged record **without blocking**. Any step that needs
    /// a not-yet-ready slot block returns [`MergeStep::Stalled`].
    ///
    /// # Errors
    ///
    /// Propagates spill-file corruption or decompression errors from the
    /// underlying slot sources.
    fn try_step(&mut self) -> Result<MergeStep<'_>>;
    /// Like [`Self::try_step`], but append the winner to `out` in the spill
    /// framing — `[key][u32 len][body]` ([`crate::frame_keyed_record_into`]),
    /// or `[u32 len][body]` for a key embedded in the record — read through the
    /// borrowed winner, so a consolidated run is framed exactly as a run the
    /// spill writer wrote.
    ///
    /// # Errors
    ///
    /// As [`Self::try_step`], plus a framing error (a record too long for its
    /// length prefix).
    fn try_step_framed(&mut self, out: &mut Vec<u8>) -> Result<MergeStepFramed>;
    /// The spill key kind the driver was built for ([`merge_driver_for_kind`];
    /// tests).
    #[cfg(test)]
    fn key_kind(&self) -> Option<crate::SpillKeyKind>;
    /// Total records emitted since the driver was constructed.
    fn records_merged(&self) -> u64;
    /// The slot whose empty FIFO caused the last [`MergeStep::Stalled`], or
    /// `None` after a `Produced` (a memory source never stalls). The merge
    /// awaits this slot so the delivery that refills it can wake the merge.
    fn stalled_slot(&self) -> Option<&Arc<crate::merge_slots::SortMergeSlot>>;
    /// The slot whose read failed the last `try_step`, if that read was a
    /// slot's (the caller names the failed input).
    fn failed_slot(&self) -> Option<&Arc<crate::merge_slots::SortMergeSlot>>;
    /// `(borrowed, assembled)`: records emitted by reference to a source's
    /// block or shared residual chunk, and records emitted from the driver's
    /// scratch buffer.
    fn record_counts(&self) -> (u64, u64);
    /// [`MergeDriver::with_input_order_check`] on a built driver.
    fn set_input_order_check(&mut self, inputs: Vec<String>, order: &'static str);
    /// [`MergeDriver::with_record_rewrite`] on a built driver.
    fn set_record_rewrite(&mut self, rewrite: &RecordRewrite);
}

/// Internal phase of the [`MergeDriver`] non-blocking state machine.
enum MergePhase<K> {
    /// Lazily priming the loser tree: each source's first record is pulled one
    /// at a time (non-blocking). `next` is the index of the next source to
    /// prime; `keys` accumulates the primed initial keys; `rec` is the scratch
    /// for the source currently being primed, which holds a first record that
    /// spans blocks (retained across a `Stalled` so it resumes correctly) and
    /// becomes that leaf's scratch.
    Priming { next: usize, keys: Vec<K>, rec: Vec<u8> },
    /// Active k-way merge. `pending_refill` is set after a record is emitted
    /// (`Produced`) and means the winner's source must be refilled before the
    /// next winner is produced.
    Merging { tree: crate::loser_tree::LoserTree<K>, pending_refill: bool },
    /// Exhausted.
    Done,
}

/// Pausable state machine for a non-blocking k-way merge over slot + in-memory
/// sources. Constructed via [`MergeDriver::from_slots`]; driven by the pipeline
/// `SortMerge` step via [`MergeDriverDyn::try_step`].
pub struct MergeDriver<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K> = K> {
    sources: Vec<SlotMergeSource<K, R>>,
    /// Per-active-source scratch buffer, indexed by loser-tree leaf (parallel
    /// to `source_map`): holds the leaf's current record only when `locs` says
    /// [`RecLoc::Scratch`]. Grows during `Priming`. While a deferred refill is
    /// `WouldBlock`ed on a record that spans blocks, the scratch holds that
    /// partial *next* record and `locs` still names the already-emitted one;
    /// no `Produced` reads it until the refill completes and updates `locs`.
    records: Vec<Vec<u8>>,
    /// Per loser-tree leaf, where its current record lives.
    locs: Vec<RecLoc>,
    /// Records emitted from a scratch buffer. The rest of `records_merged`
    /// were emitted by reference (see [`MergeDriverDyn::record_counts`]).
    assembled: u64,
    /// Loser-tree leaf → `sources` index. Built during `Priming`.
    source_map: Vec<usize>,
    /// Non-blocking merge state machine.
    phase: MergePhase<K>,
    records_merged: u64,
    progress: ProgressTracker,
    /// `sources` index of the slot that caused the last `Stalled`; cleared on
    /// `Produced`.
    stalled_src: Option<usize>,
    /// `sources` index of the source whose read returned the last error.
    failed_src: Option<usize>,
    /// The sort's merge-wide demand; the driver publishes predictions and the
    /// frontier through it. `None` publishes nothing.
    demand: Option<Arc<crate::MergeDemand>>,
    /// Per `sources` index, the slot id of the spill file it reads (`None` for
    /// a memory source, which has no file to read ahead): what predictions and
    /// the frontier publish.
    source_slot_id: Vec<Option<u32>>,
    /// Per `sources` index, whether a slot source has reached its end.
    slot_drained: Vec<bool>,
    /// `sources` index of the last record's source (a change is a switch).
    last_src: Option<usize>,
    /// The file the last switch predicted, scored on the next switch.
    last_predicted: Option<u32>,
    /// The key's merge-wide context for rebuilding slot keys (see
    /// [`RawSortKey::MergeCtx`]).
    ctx: K::MergeCtx,
    /// The kind [`merge_driver_for_kind`] built this driver for (tests).
    #[cfg(test)]
    kind_for_test: Option<crate::SpillKeyKind>,
    /// Check each source is sorted ([`Self::with_input_order_check`]): the
    /// inputs' names (by `file_id`) and the order's CLI name.
    order_check: Option<(Vec<String>, &'static str)>,
    /// Per `sources` index, records read so far (for the order check's
    /// 1-based record number).
    pulled: Vec<u64>,
}

/// The error of a merge input that is not sorted in the merge order:
/// `record_no` (1-based within `input`) sorts before the record read before it.
#[must_use]
pub fn input_order_violation(input: &str, order: &str, record_no: u64) -> String {
    format!(
        "Input '{input}' is not sorted in {order} order: record {record_no} sorts \
         before a preceding record from the same input, so the k-way merge would \
         corrupt the output. Sort that input in {order} order first \
         (`fgumi sort -i <input> -o sorted.bam --order {order}`), or locate the \
         violation with `fgumi sort -i <input> --verify --order {order}`."
    )
}

impl<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K>> MergeDriver<K, R> {
    /// Construct a driver from `SortMergeSlot`s and already-sorted in-memory
    /// chunks (arena-backed residuals sharing an `Arc<SegmentedBuf>`: zero
    /// per-record allocation, each record's stored ref `R` yielding its merge
    /// key through [`ChunkRef`]). Used by the cooperative pipeline `SortMerge`
    /// step.
    ///
    /// Each slot becomes a slot source read **non-blocking** by
    /// [`MergeDriverDyn::try_step`]. Priming is **lazy**: the loser tree is
    /// built on the first `try_step` call(s) as each source yields its first
    /// record, so construction never blocks on a not-yet-decompressed slot.
    /// Construction is infallible (read errors surface from `try_step`), and
    /// the empty-input case is reported by the first `try_step` returning
    /// [`MergeStep::Done`]. `ctx` is the key's merge-wide context
    /// ([`RawSortKey::MergeCtx`]; `()` for keys that need none).
    #[must_use]
    pub fn from_slots(
        slots: Vec<Arc<crate::merge_slots::SortMergeSlot>>,
        memory: Vec<InMemoryChunk<R>>,
        total_records: u64,
        ctx: K::MergeCtx,
    ) -> Self {
        let num_slots = slots.len();
        let num_memory = memory.iter().filter(|c| !c.is_empty()).count();
        let num_sources = num_slots + num_memory;

        if num_slots > 0 {
            info!(
                "Pipeline-integrated merge: {num_slots} slot sources + {num_memory} memory sources"
            );
        }

        let mut sources: Vec<SlotMergeSource<K, R>> = Vec::with_capacity(num_sources);
        // Order by `file_id`, do not trust the caller's `Vec` order. `file_id`
        // is documented on `SortMergeSlot` as the stable merge-order identity,
        // and the `LoserTree` breaks ties between equal sort keys by leaf index
        // — so if the leaf order came from push order instead, two records with
        // the same key could be emitted in a different order than the legacy
        // chunk-files path produces. That is an output-identity divergence in a
        // tool whose sort output is checked byte-for-byte, and it would be
        // invisible except on inputs that happen to contain ties. Sorting a
        // handful of slots once per merge costs nothing next to that.
        let mut slots = slots;
        slots.sort_by_key(|s| s.file_id);
        for slot in slots {
            sources.push(SlotMergeSource::Slot {
                parser: SlotParserState::for_slot(&slot),
                slot,
                hint: K::MergeHint::default(),
            });
        }
        for chunk in memory {
            if !chunk.is_empty() {
                sources.push(SlotMergeSource::MemoryShared { chunk, idx: 0 });
            }
        }

        let progress = ProgressTracker::new("Merged records")
            .with_interval(1_000_000)
            .with_total(total_records);
        let source_slot_id: Vec<Option<u32>> = sources
            .iter()
            .map(|s| match s {
                SlotMergeSource::Slot { slot, .. } => Some(slot.slot_id()),
                SlotMergeSource::MemoryShared { .. } => None,
            })
            .collect();
        let slot_drained = vec![false; sources.len()];

        Self {
            sources,
            records: Vec::with_capacity(num_sources),
            locs: Vec::with_capacity(num_sources),
            assembled: 0,
            source_map: Vec::with_capacity(num_sources),
            // Lazy priming: tree built on first try_step as sources yield.
            phase: MergePhase::Priming {
                next: 0,
                keys: Vec::with_capacity(num_sources),
                rec: Vec::new(),
            },
            records_merged: 0,
            progress,
            stalled_src: None,
            failed_src: None,
            demand: None,
            source_slot_id,
            slot_drained,
            last_src: None,
            last_predicted: None,
            ctx,
            #[cfg(test)]
            kind_for_test: None,
            order_check: None,
            pulled: vec![0; num_sources],
        }
    }

    /// Fail when a source is not sorted: a record whose key sorts before the
    /// previous record of its source ([`RawSortKey::input_order_cmp`]) is an
    /// [`input_order_violation`]`(inputs[file_id], order, record_no)` error.
    /// For a merge of independently sorted inputs, whose order the merge
    /// assumes; a spill run is sorted by construction.
    #[must_use]
    pub fn with_input_order_check(mut self, inputs: Vec<String>, order: &'static str) -> Self {
        MergeDriverDyn::set_input_order_check(&mut self, inputs, order);
        self
    }

    /// Apply `rewrite` to every slot source's records before their keys are
    /// extracted (see [`RecordRewrite`]).
    #[must_use]
    pub fn with_record_rewrite(mut self, rewrite: &RecordRewrite) -> Self {
        MergeDriverDyn::set_record_rewrite(&mut self, rewrite);
        self
    }

    /// Record `src` as the failed source and, for a merge of named inputs
    /// ([`Self::with_input_order_check`]), name the input and the record being
    /// read in `e`.
    fn source_error(&mut self, src: usize, e: anyhow::Error) -> anyhow::Error {
        self.failed_src = Some(src);
        let Some((inputs, _)) = &self.order_check else { return e };
        let SlotMergeSource::Slot { slot, .. } = &self.sources[src] else { return e };
        let name = inputs.get(slot.file_id as usize).map_or("?", String::as_str);
        e.context(format!("Input '{name}': record {}", self.pulled[src] + 1))
    }

    /// The order-check error for `src`'s latest record, if it sorts before the
    /// record `prev` its source emitted last under
    /// [`RawSortKey::input_order_cmp`].
    fn order_violation(&self, src: usize, key: &K, prev: &K) -> Option<anyhow::Error> {
        let (inputs, order) = self.order_check.as_ref()?;
        if key.input_order_cmp(prev) != std::cmp::Ordering::Less {
            return None;
        }
        let SlotMergeSource::Slot { slot, .. } = &self.sources[src] else { return None };
        let name = inputs.get(slot.file_id as usize).map_or("?", String::as_str);
        Some(anyhow::anyhow!(input_order_violation(name, order, self.pulled[src])))
    }

    /// The sort's merge-wide demand; the driver publishes predictions and the
    /// frontier through it. Publishes the initial frontier (the lowest slot's
    /// file) at once, so a reader never sees "every file drained" before the
    /// merge has read anything.
    #[must_use]
    pub fn with_demand(mut self, demand: Arc<crate::MergeDemand>) -> Self {
        self.demand = Some(demand);
        self.publish_frontier();
        self
    }

    /// Mark `src` as having reached its end; a slot source's drain moves the
    /// frontier.
    fn on_source_end(&mut self, src: usize) {
        if self.source_slot_id[src].is_some() {
            self.slot_drained[src] = true;
            self.publish_frontier();
        }
    }

    /// On a source switch publish the runner-up's file (a memory source is
    /// never predicted — it has no file to read — so it publishes `None`) and
    /// score the previous switch's prediction: a hit when the new winner is
    /// the source it named, a memory source included. The first record's
    /// switch has no earlier prediction to score. This runs on every switch,
    /// which in an interleaved merge is most records: the caller's
    /// `runner_up` walk (O(log k) key compares), one or two relaxed stat
    /// increments, and a store only when the prediction changes.
    fn on_source_switch(&mut self, runner_up_src: Option<usize>, src: usize) {
        let first = self.last_src.is_none();
        self.last_src = Some(src);
        let Some(d) = &self.demand else { return };
        let now = self.source_slot_id[src];
        if !first {
            d.stats().record_prediction(now == self.last_predicted);
        }
        let next = runner_up_src.and_then(|s| self.source_slot_id[s]);
        self.last_predicted = next;
        d.set_predicted(next);
    }

    /// The current record bytes of loser-tree `leaf`. Valid until that leaf's
    /// source is next read, which only `try_step`'s deferred refill or priming
    /// does, after the caller has consumed the `Produced` borrow (`Produced`
    /// borrows `&mut self`). A source's block is replaced only by its own read
    /// (`slot_try_load_block`); prediction and frontier publishing never touch
    /// a source's parser.
    fn record_bytes(&self, leaf: usize) -> &[u8] {
        match self.locs[leaf] {
            RecLoc::Scratch => &self.records[leaf],
            RecLoc::Block { start, end } => match &self.sources[self.source_map[leaf]] {
                SlotMergeSource::Slot { parser, .. } => &parser.current_buf[start..end],
                SlotMergeSource::MemoryShared { .. } => {
                    unreachable!("a block location on a non-slot source")
                }
            },
            RecLoc::Chunk { idx } => match &self.sources[self.source_map[leaf]] {
                SlotMergeSource::MemoryShared { chunk, .. } => chunk.record_bytes(idx),
                SlotMergeSource::Slot { .. } => {
                    unreachable!("a chunk location on a non-chunk source")
                }
            },
        }
    }

    /// Publish the lowest undrained slot source's slot id (slot sources come
    /// first, in `file_id` order). O(k), run only on a drain.
    fn publish_frontier(&self) {
        if let Some(d) = &self.demand {
            let frontier = (0..self.sources.len())
                .find(|&i| self.source_slot_id[i].is_some() && !self.slot_drained[i])
                .and_then(|i| self.source_slot_id[i]);
            d.set_frontier(frontier);
        }
    }
}

impl<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K>> MergeDriverDyn
    for MergeDriver<K, R>
{
    fn try_step(&mut self) -> Result<MergeStep<'_>> {
        // Move the phase out so the body works with owned `tree`/`keys`/`rec`
        // and can borrow `self.{sources,records,source_map}` freely without
        // aliasing the phase; the phase is written back before every return.
        // `Done` is the placeholder; any early `?` leaves the driver `Done`,
        // which is safe because the whole merge aborts on error.
        loop {
            match std::mem::replace(&mut self.phase, MergePhase::Done) {
                MergePhase::Priming { mut next, mut keys, mut rec } => {
                    let mut stalled = false;
                    while next < self.sources.len() {
                        let read = self.sources[next].try_next_record(&mut rec, &self.ctx);
                        let read = read.map_err(|e| self.source_error(next, e))?;
                        match read {
                            TryRead::WouldBlock => {
                                stalled = true;
                                self.stalled_src = Some(next);
                                break;
                            }
                            TryRead::Ready(Some((key, loc))) => {
                                self.pulled[next] += 1;
                                keys.push(key);
                                // Only a scratch-held record hands its buffer to
                                // the leaf; otherwise `rec` (empty) serves the
                                // next source.
                                let scratch = if loc == RecLoc::Scratch {
                                    std::mem::take(&mut rec)
                                } else {
                                    Vec::new()
                                };
                                self.records.push(scratch);
                                self.locs.push(loc);
                                self.source_map.push(next);
                                next += 1;
                            }
                            TryRead::Ready(None) => {
                                rec.clear();
                                self.on_source_end(next);
                                next += 1;
                            }
                        }
                    }
                    if stalled {
                        self.phase = MergePhase::Priming { next, keys, rec };
                        return Ok(MergeStep::Stalled);
                    }
                    if keys.is_empty() {
                        return Ok(MergeStep::Done);
                    }
                    info!("Merging from {} sources...", keys.len());
                    let tree = crate::loser_tree::LoserTree::new(keys);
                    self.phase = MergePhase::Merging { tree, pending_refill: false };
                    // Loop to emit the first winner.
                }
                MergePhase::Merging { mut tree, pending_refill } => {
                    // Resolve the refill deferred from the previous `Produced`
                    // (the read that would overwrite the just-emitted winner's
                    // buffer). On `WouldBlock`, yield with the refill still
                    // pending.
                    if pending_refill && tree.winner_is_active() {
                        let winner = tree.winner();
                        let src = self.source_map[winner];
                        let read =
                            self.sources[src].try_next_record(&mut self.records[winner], &self.ctx);
                        let read = read.map_err(|e| self.source_error(src, e))?;
                        match read {
                            TryRead::WouldBlock => {
                                self.stalled_src = Some(src);
                                self.phase = MergePhase::Merging { tree, pending_refill: true };
                                return Ok(MergeStep::Stalled);
                            }
                            TryRead::Ready(Some((key, loc))) => {
                                self.pulled[src] += 1;
                                if let Some(e) = self.order_violation(src, &key, tree.winner_key())
                                {
                                    return Err(e);
                                }
                                self.locs[winner] = loc;
                                tree.replace_winner(key);
                            }
                            TryRead::Ready(None) => {
                                tree.remove_winner();
                                self.on_source_end(src);
                            }
                        }
                    }
                    if !tree.winner_is_active() {
                        return Ok(MergeStep::Done);
                    }
                    let winner = tree.winner();
                    let src = self.source_map[winner];
                    if self.last_src != Some(src) {
                        let runner_up = tree.runner_up().map(|leaf| self.source_map[leaf]);
                        self.on_source_switch(runner_up, src);
                    }
                    self.records_merged += 1;
                    self.progress.log_if_needed(1);
                    // Defer this winner's refill to the next call so the
                    // returned borrow stays valid until the caller copies it.
                    self.phase = MergePhase::Merging { tree, pending_refill: true };
                    self.stalled_src = None;
                    if self.locs[winner] == RecLoc::Scratch {
                        self.assembled += 1;
                    }
                    return Ok(MergeStep::Produced(self.record_bytes(winner)));
                }
                MergePhase::Done => return Ok(MergeStep::Done),
            }
        }
    }

    fn try_step_framed(&mut self, out: &mut Vec<u8>) -> Result<MergeStepFramed> {
        match self.try_step()? {
            MergeStep::Produced(_) => {}
            MergeStep::Stalled => return Ok(MergeStepFramed::Stalled),
            MergeStep::Done => return Ok(MergeStepFramed::Done),
        }
        // `Produced` defers the winner's refill to the next step, so the tree's
        // winner is still the record just produced, with its key.
        let MergePhase::Merging { tree, .. } = &self.phase else {
            unreachable!("a produced record leaves the driver merging")
        };
        let winner = tree.winner();
        let body = self.record_bytes(winner);
        if K::EMBEDDED_IN_RECORD {
            crate::spill_block::frame_embedded_record_into(out, body)?;
        } else {
            crate::spill_block::frame_keyed_record_into(out, tree.winner_key(), body)?;
        }
        Ok(MergeStepFramed::Produced)
    }

    #[cfg(test)]
    fn key_kind(&self) -> Option<crate::SpillKeyKind> {
        self.kind_for_test
    }

    fn records_merged(&self) -> u64 {
        self.records_merged
    }

    fn record_counts(&self) -> (u64, u64) {
        // Every `Produced` counts in `records_merged`; only a scratch one also
        // counts in `assembled`, so the borrowed share is the difference.
        (self.records_merged - self.assembled, self.assembled)
    }

    fn set_input_order_check(&mut self, inputs: Vec<String>, order: &'static str) {
        self.order_check = Some((inputs, order));
    }

    fn set_record_rewrite(&mut self, rewrite: &RecordRewrite) {
        for s in &mut self.sources {
            if let SlotMergeSource::Slot { parser, .. } = s {
                parser.rewrite = Some(Arc::clone(rewrite));
            }
        }
    }

    fn stalled_slot(&self) -> Option<&Arc<crate::merge_slots::SortMergeSlot>> {
        match self.sources.get(self.stalled_src?)? {
            SlotMergeSource::Slot { slot, .. } => Some(slot),
            SlotMergeSource::MemoryShared { .. } => None,
        }
    }

    fn failed_slot(&self) -> Option<&Arc<crate::merge_slots::SortMergeSlot>> {
        match self.sources.get(self.failed_src?)? {
            SlotMergeSource::Slot { slot, .. } => Some(slot),
            SlotMergeSource::MemoryShared { .. } => None,
        }
    }
}

/// The merge-wide context of a kind-dispatched driver
/// ([`merge_driver_for_kind`]): the class table for the queryname kinds, none
/// for the rest.
#[derive(Clone)]
pub enum MergeCtxErased {
    /// No context (coordinate and template kinds; a queryname kind then uses an
    /// empty class table, which is correct, only slower).
    None,
    /// The sort's queryname class table.
    Queryname(Arc<crate::QnClassTable>),
    /// The full template key's context: a merge of independently sorted
    /// template-coordinate inputs (`fgumi merge`, body-keyed slots).
    Template(Arc<crate::TemplateKeyContext>),
}

/// A merge driver over `slots` (spill runs, in run order) for runs written with
/// `kind`, publishing its demand through `demand`: the dispatch from a run's
/// recorded key kind to `MergeDriver<K>`. Queryname kinds merge on the class
/// window ([`crate::QnMergeKey`]) with `ctx`'s table; the framing they write
/// back is the embedded `[u32 len][body]` either way.
#[must_use]
pub fn merge_driver_for_kind(
    kind: crate::SpillKeyKind,
    slots: Vec<Arc<SortMergeSlot>>,
    ctx: MergeCtxErased,
    demand: Arc<crate::MergeDemand>,
) -> Box<dyn MergeDriverDyn + Send> {
    use crate::SpillKeyKind as Kind;
    use crate::{QnKey, QnLex, QnMergeKey, QnNatural};
    fn build<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K>>(
        slots: Vec<Arc<SortMergeSlot>>,
        ctx: K::MergeCtx,
        demand: Arc<crate::MergeDemand>,
    ) -> Box<dyn MergeDriverDyn + Send> {
        #[cfg_attr(not(test), allow(unused_mut))]
        let mut d = MergeDriver::<K, R>::from_slots(slots, Vec::new(), 0, ctx).with_demand(demand);
        #[cfg(test)]
        {
            d.kind_for_test = Some(spill_kind_of::<K>());
        }
        Box::new(d)
    }
    let table = |ctx: MergeCtxErased| match ctx {
        MergeCtxErased::Queryname(t) => t,
        MergeCtxErased::None | MergeCtxErased::Template(_) => crate::QnClassTable::new(),
    };
    let template = |ctx: MergeCtxErased| match ctx {
        MergeCtxErased::Template(c) => Some(c),
        MergeCtxErased::None | MergeCtxErased::Queryname(_) => None,
    };
    match kind {
        Kind::Coordinate => {
            build::<crate::keys::RawCoordinateKey, crate::keys::RawCoordinateKey>(slots, (), demand)
        }
        Kind::QuerynameLex => build::<QnMergeKey<QnLex>, QnKey>(slots, table(ctx), demand),
        Kind::QuerynameNatural => build::<QnMergeKey<QnNatural>, QnKey>(slots, table(ctx), demand),
        Kind::TemplateK24 => build::<crate::TemplateKey24, crate::TemplateKey24>(slots, (), demand),
        Kind::TemplateCb32 => build::<crate::CbKey32, crate::CbKey32>(slots, (), demand),
        Kind::TemplateTert32 => build::<crate::TertKey32, crate::TertKey32>(slots, (), demand),
        Kind::TemplateK40 => {
            build::<crate::TemplateKey, crate::TemplateKey>(slots, template(ctx), demand)
        }
    }
}

/// The spill kind whose runs `K` merges, derived from the key type itself, so a
/// dispatch arm instantiating the wrong key type is caught (tests).
#[cfg(test)]
fn spill_kind_of<K: 'static>() -> crate::SpillKeyKind {
    use crate::SpillKeyKind as Kind;
    use crate::{QnLex, QnMergeKey, QnNatural};
    use std::any::TypeId;
    let t = TypeId::of::<K>();
    [
        (TypeId::of::<crate::keys::RawCoordinateKey>(), Kind::Coordinate),
        (TypeId::of::<QnMergeKey<QnLex>>(), Kind::QuerynameLex),
        (TypeId::of::<QnMergeKey<QnNatural>>(), Kind::QuerynameNatural),
        (TypeId::of::<crate::TemplateKey24>(), Kind::TemplateK24),
        (TypeId::of::<crate::CbKey32>(), Kind::TemplateCb32),
        (TypeId::of::<crate::TertKey32>(), Kind::TemplateTert32),
        (TypeId::of::<crate::TemplateKey>(), Kind::TemplateK40),
    ]
    .into_iter()
    .find_map(|(id, kind)| (id == t).then_some(kind))
    .unwrap_or_else(|| panic!("no spill kind merges {}", std::any::type_name::<K>()))
}

#[cfg(test)]
impl<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K>> MergeDriver<K, R> {
    /// Whether `ptr` lies inside loser-tree `leaf`'s slot source's current block.
    fn block_contains_for_test(&self, leaf: usize, ptr: *const u8) -> bool {
        match &self.sources[self.source_map[leaf]] {
            SlotMergeSource::Slot { parser, .. } => {
                parser.current_buf.as_ptr_range().contains(&ptr)
            }
            SlotMergeSource::MemoryShared { .. } => false,
        }
    }

    /// Capacity of loser-tree `leaf`'s scratch buffer.
    fn scratch_capacity_for_test(&self, leaf: usize) -> usize {
        self.records[leaf].capacity()
    }
}

impl<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K>> Drop for MergeDriver<K, R> {
    fn drop(&mut self) {
        self.progress.log_final();
    }
}

#[cfg(test)]
mod tests;
