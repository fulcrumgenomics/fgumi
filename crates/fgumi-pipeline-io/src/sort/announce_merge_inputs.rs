//! `AnnounceMergeInputs` — the source of `fgumi merge`'s chain: announces each
//! input BAM as a merge slot, so the rest of the chain is the sort's phase 2.
//!
//! Each input becomes a *body-keyed* slot ([`fgumi_sort::open_bam_slot`]): its
//! BGZF stream read from byte 0 through the spill read planner, its decompressed
//! header skipped, its records keyed from their bodies. A file input is a
//! positional slot of the file's length; a stream input (pipe, FIFO, process
//! substitution) is a sequential slot over a
//! [`SequentialSource`](fgumi_bam_io::pread::SequentialSource), whose header
//! bytes the command's header parse already consumed (the source replays them).
//! The step emits one `SpillReady` per input, in input order (the merge's
//! tie-break order), then `AllAnnounced`.

use std::collections::VecDeque;
use std::io;
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_bam_io::pread::PositionalSource;
use fgumi_pipeline_core::{
    HeldRetry, Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

use super::protocol::SortPhase1Event;

/// One input of a merge.
pub struct MergeInput {
    /// The input's name (errors, `SpillReady.path`).
    pub path: PathBuf,
    /// The input's bytes: a file, or a stream's sequential source.
    pub source: Arc<dyn PositionalSource>,
    /// The file's length; `None` for a stream.
    pub len: Option<u64>,
    /// For a stream, the bytes its header parse consumed: the header is found
    /// within them, and reading past them would move the stream's cursor
    /// before the slot reads it. `None` for a file.
    pub prefix_len: Option<u64>,
}

/// Bytes read per step while locating the end of a header.
const HEADER_READ_BYTES: usize = 64 * 1024;

/// The decompressed length of the BAM header at the front of `source` (magic,
/// text and references: the bytes a body-keyed slot skips), reading complete
/// BGZF blocks from offset 0 and never past `limit`.
///
/// # Errors
///
/// A read error, a malformed block, or a header that does not end within the
/// bytes read.
pub fn bam_header_len_of(
    source: &dyn PositionalSource,
    path: &std::path::Path,
    limit: Option<u64>,
) -> io::Result<u64> {
    let mut compressed: Vec<u8> = Vec::new();
    let mut decompressed: Vec<u8> = Vec::new();
    let mut decompressor = fgumi_bgzf::Decompressor::new();
    let mut at = 0usize; // next undecoded byte of `compressed`
    let mut offset = 0u64;
    let mut eof = false;
    let corrupt = |reason: String| {
        io::Error::new(io::ErrorKind::InvalidData, format!("{}: {reason}", path.display()))
    };
    loop {
        // Decode every complete block buffered so far, framed by the shared
        // BGZF header check and `BSIZE` floor.
        while compressed.len() - at >= fgumi_bgzf::BGZF_HEADER_SIZE {
            let block = &compressed[at..];
            fgumi_bgzf::header::validate(block).map_err(|e| corrupt(e.to_string()))?;
            let block_size = fgumi_bgzf::block_size_checked(block).ok_or_else(|| {
                corrupt(format!(
                    "a BGZF block's BSIZE is below the {} byte minimum",
                    fgumi_bgzf::MIN_BLOCK_SIZE
                ))
            })?;
            if block.len() < block_size {
                break;
            }
            fgumi_bgzf::decompress_block_slice_into(
                &block[..block_size],
                &mut decompressor,
                &mut decompressed,
            )?;
            at += block_size;
            if let Some(n) = crate::boundaries::bam_header_len(&decompressed)? {
                return Ok(n as u64);
            }
        }
        if eof {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("{}: the input ends inside its BAM header", path.display()),
            ));
        }
        let want = match limit {
            Some(l) => usize::try_from(l.saturating_sub(offset)).unwrap_or(usize::MAX),
            None => HEADER_READ_BYTES,
        }
        .min(HEADER_READ_BYTES);
        if want == 0 {
            eof = true;
            continue;
        }
        let mut buf = vec![0u8; want];
        let n = source.read_at(&mut buf, offset)?;
        if n == 0 {
            eof = true;
        }
        compressed.extend_from_slice(&buf[..n]);
        offset += n as u64;
    }
}

/// `Exclusive` source announcing the merge inputs (see the module docs).
pub struct AnnounceMergeInputs {
    inputs: Vec<MergeInput>,
    outbox: VecDeque<SortPhase1Event>,
    held: HeldSlot<Unpushed<SortPhase1Event>>,
    announced: bool,
    output_byte_limit: u64,
}

impl AnnounceMergeInputs {
    /// A source announcing `inputs`, in order.
    #[must_use]
    pub fn new(inputs: Vec<MergeInput>, output_byte_limit: u64) -> Self {
        Self {
            inputs,
            outbox: VecDeque::new(),
            held: HeldSlot::new(),
            announced: false,
            output_byte_limit,
        }
    }

    /// Open every input as a slot and queue the events.
    fn announce(&mut self) -> io::Result<()> {
        let n = u32::try_from(self.inputs.len())
            .map_err(|_| io::Error::other("more merge inputs than slot ids"))?;
        for (i, input) in std::mem::take(&mut self.inputs).into_iter().enumerate() {
            let MergeInput { path, source, len, prefix_len } = input;
            let header_len = bam_header_len_of(&*source, &path, prefix_len)?;
            let file_id = u32::try_from(i).expect("bounded by n above");
            let slot = fgumi_sort::open_bam_slot(source, file_id, header_len, len)
                .map_err(|e| io::Error::other(format!("{}: {e:#}", path.display())))?;
            self.outbox.push_back(SortPhase1Event::SpillReady {
                slot,
                path,
                records_ingested_so_far: 0,
            });
        }
        self.outbox.push_back(SortPhase1Event::AllAnnounced {
            slot_count: n,
            memory_chunk_count: 0,
            total_records: 0,
        });
        Ok(())
    }
}

impl Step for AnnounceMergeInputs {
    type Input = ();
    type Outputs = Single<SortPhase1Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "AnnounceMergeInputs",
            kind: StepKind::Exclusive,
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if matches!(ctx.outputs.retry_held(&mut self.held), HeldRetry::StillHeld) {
            return Ok(StepOutcome::Contention);
        }
        if !self.announced {
            self.announced = true;
            self.announce()?;
        }
        if let Some(ev) = self.outbox.pop_front() {
            if let Err(u) = ctx.outputs.push(ev) {
                self.held.put(u);
            }
            return Ok(StepOutcome::Progress);
        }
        Ok(StepOutcome::Finished)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A positional source over `bytes`.
    fn source(bytes: &[u8]) -> std::fs::File {
        let mut f = tempfile::tempfile().unwrap();
        io::Write::write_all(&mut f, bytes).unwrap();
        f
    }

    /// A BGZF block holding `bytes`.
    fn bgzf_block(bytes: &[u8]) -> Vec<u8> {
        let mut w = noodles::bgzf::io::Writer::new(Vec::new());
        io::Write::write_all(&mut w, bytes).unwrap();
        w.finish().unwrap()
    }

    /// A header split over two blocks is found; a block whose `BSIZE` is below
    /// the BGZF minimum, or a stream that is not BGZF, is an `InvalidData`
    /// error naming the input, not a slice past the buffer. Catches: framing by
    /// a raw `BSIZE` with no floor.
    #[test]
    fn header_length_is_read_through_checked_bgzf_framing() {
        let header = noodles::sam::Header::default();
        let mut raw = Vec::new();
        fgumi_bam_io::write_bam_header(&mut raw, &header).unwrap();
        let (a, b) = raw.split_at(raw.len() / 2);
        let mut file = bgzf_block(a);
        file.truncate(file.len() - fgumi_bgzf::BGZF_EOF.len());
        file.extend(bgzf_block(b));
        let path = std::path::Path::new("in.bam");
        assert_eq!(bam_header_len_of(&source(&file), path, None).unwrap(), raw.len() as u64);

        let mut bad = file.clone();
        bad[16..18].copy_from_slice(&3u16.to_le_bytes());
        let err = bam_header_len_of(&source(&bad), path, None).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("in.bam") && err.to_string().contains("BSIZE"), "{err}");

        let err =
            bam_header_len_of(&source(b"@HD\tVN:1.6\n@SQ\tSN:x\tLN:9\n"), path, None).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{err}");
    }
}
