//! `AlignAndMerge`: a [`SpecialStage`] that wraps an external aligner
//! subprocess (e.g. `bwa mem`) with integrated FASTQ generation and
//! zipper merge.
//!
//! ## Input
//!
//! [`SerializedBatch`] — length-prefixed unmapped BAM record bytes
//! (typically emitted by `ExtractStage` or `CorrectStage`).
//!
//! ## Output
//!
//! [`SerializedBatch`] — length-prefixed merged BAM record bytes with
//! alignment information from the aligner plus the UMI/auxiliary tags
//! preserved from the unmapped input.
//!
//! ## Ordering guarantees
//!
//! Best-effort: thread A forwards the unmapped bytes of each input batch
//! to the merge thread in input-ordinal order, and the aligner emits
//! SAM in the order it consumed FASTQ. However, `bwa mem -t N` with
//! `N > 1` reorders its output when worker threads finish alignments
//! out of submission order, so the merged output is **not**
//! deterministically ordered for `-t > 1`. Downstream stages that need
//! canonical order must call a reorder/sort stage. The emitted
//! `ordinal` is a fresh monotonic counter stamped at merge time and is
//! therefore meaningful only for in-run batch identity, not for
//! reconstructing input order.
//!
//! ## Memory model
//!
//! Three cooperating threads share bounded crossbeam channels
//! (`unmapped_tx` capped at 50 000 batches of raw unmapped bytes,
//! `mapped_tx` capped at 10 000 mapped templates). Each input batch's
//! primary bytes are cloned once for the merge thread; FASTQ is
//! rebuilt per batch into a reusable buffer. Peak memory is dominated
//! by the in-flight aligner pipe + these two channels.
//!
//! ## Determinism
//!
//! Not bit-identical across runs when the aligner uses more than one
//! thread (`bwa mem -t > 1`): template arrival order at the merge
//! thread depends on aligner-internal scheduling. At aligner `-t 1`
//! plus pool `-t 1`, runs are byte-identical.

use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::thread;

use anyhow::{Result, anyhow};
use fgumi_raw_bam::RawRecord;
use fgumi_umi::TagInfo;
use noodles::sam::Header;

use crate::aligner::AlignerProcess;
use crate::commands::zipper::{build_output_header, merge_raw as merge};
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::output_types::{RawBytes, SerializedBatch};
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::special_stage::SpecialStage;
use crate::runall::engine::stage::SequencedItem;
use crate::runall::engine::stages::tofastq::{DEFAULT_EXCLUDE_FLAGS, bam_records_to_fastq};

/// Configuration for the `AlignAndMerge` special stage.
pub struct AlignAndMerge {
    /// Shell command for the aligner (e.g. `bwa mem -p -K 10000000 -t 4 ref.fa /dev/stdin`).
    aligner_command: String,
    /// Reference FASTA path (retained for future use, e.g. bisulfite restoration).
    #[allow(dead_code)]
    reference: PathBuf,
    /// Path to the reference dictionary file.
    dict_path: PathBuf,
    /// Unmapped BAM header (from the upstream extract/correct stage).
    unmapped_header: Header,
    /// Tag info for zipper merge (tags to remove, reverse, revcomp).
    tag_info: TagInfo,
    /// Skip adding PA (primary alignment) tags during merge.
    skip_pa_tags: bool,
    /// Stderr ring buffer size for error reporting.
    stderr_ring_size: usize,
    /// Whether to suppress /1 /2 suffixes in FASTQ output.
    no_read_suffix: bool,
    /// Shared slot populated at runtime with the merged output header
    /// (built from unmapped + aligner-produced mapped header + dict). When
    /// the plan wires a [`BamSink`] via `with_deferred_header`, that sink
    /// blocks on this slot and uses its value as the BAM writer's header.
    /// When `None`, the merged header is still computed (for early
    /// validation) but discarded.
    ///
    /// [`BamSink`]: crate::runall::engine::sink::BamSink
    deferred_header: Option<crate::runall::engine::sink::DeferredHeader>,
}

impl AlignAndMerge {
    /// Construct a new `AlignAndMerge` stage.
    #[must_use]
    pub fn new(
        aligner_command: String,
        reference: PathBuf,
        dict_path: PathBuf,
        unmapped_header: Header,
        tag_info: TagInfo,
        skip_pa_tags: bool,
        no_read_suffix: bool,
    ) -> Self {
        Self {
            aligner_command,
            reference,
            dict_path,
            unmapped_header,
            tag_info,
            skip_pa_tags,
            stderr_ring_size: 50,
            no_read_suffix,
            deferred_header: None,
        }
    }

    /// Attach the downstream `BamSink`'s deferred-header slot. When the plan
    /// includes both this stage and a [`BamSink`], the planner calls this so
    /// the sink can use the aligner-aware output header instead of the
    /// plan-time approximation produced by `build_zipper_sink_header`.
    ///
    /// [`BamSink`]: crate::runall::engine::sink::BamSink
    #[must_use]
    pub fn with_deferred_header(
        mut self,
        deferred_header: crate::runall::engine::sink::DeferredHeader,
    ) -> Self {
        self.deferred_header = Some(deferred_header);
        self
    }
}

impl SpecialStage for AlignAndMerge {
    type Input = SerializedBatch;
    type Output = SerializedBatch;

    #[allow(
        clippy::too_many_lines,
        reason = "thread spawning + join logic is cohesive and readable as one block"
    )]
    #[tracing::instrument(name = "align_and_merge", skip_all)]
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<SerializedBatch>>,
        output: Box<dyn OutputQueue<SerializedBatch>>,
        cancel: CancelToken,
    ) -> Result<()> {
        // 1. Spawn aligner subprocess.
        let mut aligner = AlignerProcess::spawn(&self.aligner_command, self.stderr_ring_size)?;
        let mut stdin = aligner.take_stdin().expect("aligner stdin was piped");
        let stdout = aligner.take_stdout().expect("aligner stdout was piped");

        // 2. Channel for the mapped header: thread B -> thread C (single-shot).
        let (header_tx, header_rx) = std::sync::mpsc::sync_channel::<Header>(1);

        // 3. Thread A -> Thread C: unmapped BAM record bytes, one chunk per input batch.
        //    Capped at 256 in-flight chunks (~10 MB at typical batch sizes);
        //    backpressure flows into the upstream `StageQueue` which is itself
        //    `MemoryTracker`-bounded, so this channel doesn't need to be.
        let (unmapped_tx, unmapped_rx) = crossbeam_channel::bounded::<Vec<u8>>(256);
        // 4. Thread B -> Thread C: mapped templates from aligner stdout.
        let (mapped_tx, mapped_rx) =
            crossbeam_channel::bounded::<Result<crate::template::Template>>(1024);

        let cancel_a = cancel.clone();
        let cancel_b = cancel.clone();
        let cancel_c = cancel.clone();

        let no_read_suffix = self.no_read_suffix;

        // Thread A: stdin writer + unmapped record forwarder.
        let thread_a = thread::Builder::new().name("align_and_merge_stdin".into()).spawn(
            move || -> Result<()> {
                let mut fastq_buf: Vec<u8> = Vec::with_capacity(256 * 1024);

                loop {
                    if cancel_a.is_cancelled() {
                        break;
                    }
                    let Some(item) = input.pop() else {
                        if input.is_drained() {
                            break;
                        }
                        thread::yield_now();
                        continue;
                    };

                    let batch = item.item;

                    // Convert to FASTQ and write to aligner stdin (borrows the
                    // batch bytes). Then move those bytes to thread C for merge
                    // — no clone, one allocation, one memcpy-equivalent pass.
                    fastq_buf.clear();
                    bam_records_to_fastq(
                        &batch.primary.data,
                        &mut fastq_buf,
                        no_read_suffix,
                        DEFAULT_EXCLUDE_FLAGS,
                        0,
                    )?;
                    if !fastq_buf.is_empty()
                        && let Err(e) = stdin.write_all(&fastq_buf)
                    {
                        // Broken pipe = aligner died; the real error surfaces via wait().
                        tracing::debug!("stdin write failed (aligner likely exited): {e}");
                        break;
                    }

                    if unmapped_tx.send(batch.primary.data).is_err() {
                        break; // thread C dropped
                    }
                }

                drop(stdin); // close aligner stdin -> aligner sees EOF
                drop(unmapped_tx); // signal thread C: no more unmapped records
                Ok(())
            },
        )?;

        // Thread B: stdout reader (SAM parsing -> templates).
        let thread_b = thread::Builder::new().name("align_and_merge_stdout".into()).spawn(
            move || -> Result<()> {
                let buf_reader = BufReader::with_capacity(256 * 1024, stdout);
                let mut sam_reader = noodles::sam::io::Reader::new(buf_reader);
                let mapped_header = sam_reader.read_header()?;

                header_tx
                    .send(mapped_header.clone())
                    .map_err(|_| anyhow!("AlignAndMerge: merge thread dropped header channel"))?;

                let mh = mapped_header;
                let mut encoder = fgumi_raw_bam::RecordBufEncoder::new(&mh);
                let record_iter = sam_reader.record_bufs(&mh).map(move |r| {
                    r.map_err(anyhow::Error::from).and_then(|rb| {
                        encoder
                            .encode(&rb)
                            .map_err(|e| anyhow!("AlignAndMerge: encode failed: {e}"))
                    })
                });
                let iter = crate::template::TemplateIterFromRecords::new(record_iter);
                for template in iter {
                    if cancel_b.is_cancelled() {
                        break;
                    }
                    if mapped_tx.send(template).is_err() {
                        break; // thread C dropped
                    }
                }
                drop(mapped_tx); // signal thread C: no more mapped templates
                Ok(())
            },
        )?;

        // Thread C: zipper merge + output.
        let unmapped_header = self.unmapped_header;
        let dict_path = self.dict_path;
        let tag_info = self.tag_info;
        let skip_pa_tags = self.skip_pa_tags;
        let deferred_header = self.deferred_header.clone();

        let thread_c = thread::Builder::new().name("align_and_merge_merge".into()).spawn(
            move || -> Result<()> {
                // Ensure the output queue is closed on every exit path so downstream
                // stages observe EOF even if the merge aborts with an error.
                let result = run_merge(
                    header_rx,
                    unmapped_rx,
                    mapped_rx,
                    &*output,
                    &cancel_c,
                    &unmapped_header,
                    &dict_path,
                    &tag_info,
                    skip_pa_tags,
                    deferred_header.as_ref(),
                );
                output.close();
                result
            },
        )?;

        // Join all threads, recording the first error. Cancel the rest on failure.
        let mut first_error: Option<anyhow::Error> = None;

        match thread_a.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                cancel.cancel();
                first_error.get_or_insert(e);
            }
            Err(p) => {
                cancel.cancel();
                first_error.get_or_insert_with(|| anyhow!("stdin thread panicked: {p:?}"));
            }
        }
        match thread_b.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                cancel.cancel();
                first_error.get_or_insert(e);
            }
            Err(p) => {
                cancel.cancel();
                first_error.get_or_insert_with(|| anyhow!("stdout thread panicked: {p:?}"));
            }
        }
        match thread_c.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => {
                cancel.cancel();
                first_error.get_or_insert(e);
            }
            Err(p) => {
                cancel.cancel();
                first_error.get_or_insert_with(|| anyhow!("merge thread panicked: {p:?}"));
            }
        }

        // Wait for aligner process and check exit status.
        if let Err(e) = aligner.wait()
            && first_error.is_none()
        {
            first_error = Some(e);
        }

        match first_error {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    fn name(&self) -> &'static str {
        "AlignAndMerge"
    }
}

/// Run the zipper merge, pulling unmapped record bytes and mapped templates from
/// two channels and pushing merged [`SerializedBatch`] items to `output`. The
/// caller is responsible for closing `output` once this function returns.
#[expect(clippy::too_many_arguments, reason = "internal helper threaded through closure")]
#[expect(
    clippy::needless_pass_by_value,
    reason = "receivers are owned by this function and dropped on return"
)]
fn run_merge(
    header_rx: std::sync::mpsc::Receiver<Header>,
    unmapped_rx: crossbeam_channel::Receiver<Vec<u8>>,
    mapped_rx: crossbeam_channel::Receiver<Result<crate::template::Template>>,
    output: &dyn OutputQueue<SerializedBatch>,
    cancel: &CancelToken,
    unmapped_header: &Header,
    dict_path: &Path,
    tag_info: &TagInfo,
    skip_pa_tags: bool,
    deferred_header: Option<&crate::runall::engine::sink::DeferredHeader>,
) -> Result<()> {
    let mapped_header = header_rx.recv().map_err(|_| {
        anyhow!("AlignAndMerge: failed to receive mapped header from stdout reader")
    })?;

    // Build the coherent output header from unmapped + aligner-produced mapped
    // + reference dict. This mirrors `fgumi zipper`'s merged header. When the
    // planner wired a BamSink `with_deferred_header`, populate its slot so the
    // sink opens the BAM writer with the aligner-aware header. When no sink is
    // waiting, the call still surfaces header incompatibility errors early.
    let output_header = build_output_header(unmapped_header, &mapped_header, dict_path)?;
    if let Some(slot) = deferred_header {
        // `set` returns Err if already populated; that's benign on retries so
        // we silently ignore the second attempt.
        let _ = slot.set(output_header);
    }

    // Parse raw unmapped BAM chunks into RecordBuf and group into templates.
    let unmapped_record_iter = RawBamChunkIterator::new(unmapped_rx, unmapped_header);
    let unmapped_templates = crate::template::TemplateIterFromRecords::new(unmapped_record_iter);

    let mut mapped_peek: Option<crate::template::Template> = None;
    let mut ordinal: u64 = 0;

    for unmapped_result in unmapped_templates {
        if cancel.is_cancelled() {
            break;
        }
        let unmapped_template = unmapped_result?;

        if mapped_peek.is_none() {
            mapped_peek = match mapped_rx.recv() {
                Ok(Ok(t)) => Some(t),
                Ok(Err(e)) => return Err(e),
                Err(_) => None, // channel closed
            };
        }

        let Some(ref mut mapped_template) = mapped_peek else {
            break; // no more mapped templates; remaining unmapped are orphans
        };

        if mapped_template.name != unmapped_template.name {
            // Orphan unmapped read: skip (matches zipper exclude_missing_reads=true).
            continue;
        }

        merge(&unmapped_template, mapped_template, tag_info, skip_pa_tags)?;

        // Two to ~eight records per template (primary + supplementaries); at
        // ~200 bytes each, 2 KiB covers the typical case without realloc.
        let mut data: Vec<u8> = Vec::with_capacity(2 * 1024);
        let mut record_count: u64 = 0;
        for rec in mapped_template.records() {
            let raw: &[u8] = rec.as_ref();
            let block_size = u32::try_from(raw.len())
                .map_err(|_| anyhow!("AlignAndMerge: encoded record exceeds u32::MAX bytes"))?;
            data.extend_from_slice(&block_size.to_le_bytes());
            data.extend_from_slice(raw);
            record_count += 1;
        }

        let mem_est = data.capacity();
        let batch =
            SerializedBatch { primary: RawBytes { data, record_count }, secondary: None, ordinal };
        let seq_item = SequencedItem::new(ordinal, batch, mem_est);
        if output.push_until_cancelled(seq_item, cancel).is_err() {
            break;
        }
        ordinal += 1;
        mapped_peek = None;
    }

    // Sanity check: any leftover mapped templates indicate a name mismatch.
    let leftover = mapped_peek.map(Ok).or_else(|| mapped_rx.recv().ok()).transpose()?;
    if let Some(remaining) = leftover {
        return Err(anyhow!(
            "AlignAndMerge: processed all unmapped reads but mapped reads remain. \
             Found template '{}'. Ensure unmapped and mapped inputs have matching read names.",
            String::from_utf8_lossy(&remaining.name)
        ));
    }

    Ok(())
}

/// Iterator that consumes raw BAM byte chunks from a channel and yields
/// individual [`RawRecord`] values (zero-copy from the chunk bytes).
struct RawBamChunkIterator {
    rx: crossbeam_channel::Receiver<Vec<u8>>,
    current_chunk: Vec<u8>,
    cursor: usize,
}

impl RawBamChunkIterator {
    #[allow(clippy::needless_pass_by_value)]
    fn new(rx: crossbeam_channel::Receiver<Vec<u8>>, _header: &Header) -> Self {
        Self { rx, current_chunk: Vec::new(), cursor: 0 }
    }

    /// Advance to the next non-empty chunk when the current one is exhausted.
    fn ensure_data(&mut self) -> bool {
        while self.cursor >= self.current_chunk.len() {
            match self.rx.recv() {
                Ok(chunk) => {
                    self.current_chunk = chunk;
                    self.cursor = 0;
                }
                Err(_) => return false, // channel closed
            }
        }
        true
    }
}

impl Iterator for RawBamChunkIterator {
    type Item = Result<RawRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.ensure_data() {
            return None;
        }

        let data = &self.current_chunk;
        if data.len() < self.cursor + 4 {
            return Some(Err(anyhow!(
                "RawBamChunkIterator: truncated block_size at offset {}",
                self.cursor
            )));
        }
        let block_size =
            u32::from_le_bytes(data[self.cursor..self.cursor + 4].try_into().unwrap()) as usize;
        let rec_start = self.cursor + 4;
        let rec_end = rec_start + block_size;
        if data.len() < rec_end {
            return Some(Err(anyhow!(
                "RawBamChunkIterator: truncated record at offset {}",
                self.cursor
            )));
        }

        let raw = data[rec_start..rec_end].to_vec();
        self.cursor = rec_end;

        Some(Ok(RawRecord::from(raw)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_and_merge_name() {
        let stage = AlignAndMerge::new(
            "echo".to_string(),
            PathBuf::from("/ref.fa"),
            PathBuf::from("/ref.dict"),
            Header::default(),
            TagInfo::new(vec![], vec![], vec![]),
            false,
            false,
        );
        assert_eq!(stage.name(), "AlignAndMerge");
    }
}
