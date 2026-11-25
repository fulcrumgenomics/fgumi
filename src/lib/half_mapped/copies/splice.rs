//! Puts half-mapped templates into a BAM that has been written without them.
//!
//! Each batch of position groups is written as its own BGZF blocks, so those of the batches that
//! get no template are copied as they are.

use std::fs::{self, File};
use std::io::{self, BufReader, Read, Seek, SeekFrom, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result, ensure};
use fgumi_bgzf::{
    BGZF_EOF, Decompressor, InlineBgzfCompressor, decompress_block_into, read_raw_blocks,
};
use fgumi_pipeline_io::sink::write_bgzf::BlockEnds;
use log::info;
use rayon::prelude::*;

use super::{BatchCopies, CopyStore, FamilySizeChange, insert_framed};
use crate::pipeline::chains::FinalizeHook;

/// Batches rewritten at a time, for each thread.
const BATCHES_PER_THREAD: usize = 32;

/// Bytes copied at a time.
const COPY_BUFFER_BYTES: usize = 1 << 20;

/// BGZF blocks read at a time.
const BLOCKS_PER_READ: usize = 16;

/// Rewrites the output once the pipeline has written it, if there is anything to put into it.
pub(crate) struct SpliceFinalizeHook {
    pub(crate) store: Arc<CopyStore>,
    pub(crate) ends: Arc<BlockEnds>,
    pub(crate) path: PathBuf,
    pub(crate) compression_level: u32,
    pub(crate) threads: usize,
}

/// Whether `path` is, or is to become, a file that can be rewritten.
pub(crate) fn can_rewrite(path: &Path) -> bool {
    !fgumi_bam_io::is_stdout_path(path) && fs::metadata(path).map_or(true, |file| file.is_file())
}

impl FinalizeHook for SpliceFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        if self.store.is_empty() || !self.ends.is_complete() {
            return Ok(());
        }
        info!("Putting {} half-mapped templates into their molecules", self.store.len());
        self.splice().with_context(|| format!("failed to rewrite {}", self.path.display()))
    }
}

impl SpliceFinalizeHook {
    fn splice(&self) -> Result<()> {
        let path = fs::canonicalize(&self.path)?;
        let ends = self.ends.take();
        let batches = self.store.by_batch();
        ensure!(
            batches.last().is_none_or(|(last, _)| last + 1 < ends.len()),
            "fewer batches were written than were grouped"
        );

        let directory = path.parent().unwrap_or(Path::new("."));
        let staged = tempfile::Builder::new().prefix(".fgumi-group-").tempfile_in(directory)?;
        let mut input = File::open(&path)?;
        let mut output = staged.as_file();
        let mut buffer = vec![0; COPY_BUFFER_BYTES];
        let pool = rayon::ThreadPoolBuilder::new().num_threads(self.threads).build()?;

        let mut copied = 0;
        for window in batches.chunks(BATCHES_PER_THREAD * self.threads.max(1)) {
            let mut inputs = Vec::with_capacity(window.len());
            for (batch, refs) in window {
                let (start, end) = (ends[*batch], ends[batch + 1]);
                input.seek(SeekFrom::Start(start))?;
                let mut blocks = vec![0; usize::try_from(end - start)?];
                input.read_exact(&mut blocks)?;
                inputs.push((start, end, blocks, self.store.read(refs)?));
            }
            let rewritten: Vec<_> = pool.install(|| {
                let rewrite = |(start, end, blocks, copies): (_, _, Vec<u8>, _)| {
                    let (blocks, changes) = rewrite(&blocks, copies, self.compression_level)?;
                    Ok((start, end, blocks, changes))
                };
                inputs.into_par_iter().map(rewrite).collect::<io::Result<_>>()
            })?;
            for (start, end, blocks, changes) in rewritten {
                copy_range(&input, output, copied..start, &mut buffer)?;
                output.write_all(&blocks)?;
                copied = end;
                self.store.add_size_changes(&changes);
            }
        }
        let len = input.metadata()?.len();
        copy_range(&input, output, copied..len, &mut buffer)?;
        output.write_all(&BGZF_EOF)?;

        staged.as_file().set_permissions(fs::metadata(&path)?.permissions())?;
        staged.persist(&path)?;
        Ok(())
    }
}

fn copy_range(
    mut input: &File,
    mut output: &File,
    range: Range<u64>,
    buffer: &mut [u8],
) -> io::Result<()> {
    input.seek(SeekFrom::Start(range.start))?;
    let len = range.end - range.start;
    if let Some(buffer) = usize::try_from(len).ok().and_then(|len| buffer.get_mut(..len)) {
        input.read_exact(buffer)?;
        return output.write_all(buffer);
    }
    let copied =
        io::copy(&mut BufReader::with_capacity(buffer.len(), input.take(len)), &mut output)?;
    if copied == len { Ok(()) } else { Err(io::ErrorKind::UnexpectedEof.into()) }
}

/// Puts `copies` into the records of some BGZF blocks.
pub(super) fn rewrite(
    blocks: &[u8],
    copies: BatchCopies,
    compression_level: u32,
) -> io::Result<(Vec<u8>, Vec<FamilySizeChange>)> {
    let (bytes, changes) = insert_framed(&decompress(blocks)?, copies)?;
    Ok((compress(&bytes, compression_level)?, changes))
}

pub(super) fn decompress(mut blocks: &[u8]) -> io::Result<Vec<u8>> {
    let mut decompressor = Decompressor::new();
    let mut bytes = Vec::new();
    loop {
        let read = read_raw_blocks(&mut blocks, BLOCKS_PER_READ)?;
        if read.is_empty() {
            return Ok(bytes);
        }
        for block in read {
            decompress_block_into(&block, &mut decompressor, &mut bytes)?;
        }
    }
}

pub(super) fn compress(bytes: &[u8], compression_level: u32) -> io::Result<Vec<u8>> {
    let mut compressor = InlineBgzfCompressor::new(compression_level);
    compressor.write_all(bytes)?;
    compressor.flush()?;
    let mut blocks = Vec::new();
    compressor.write_blocks_to(&mut blocks)?;
    Ok(blocks)
}
