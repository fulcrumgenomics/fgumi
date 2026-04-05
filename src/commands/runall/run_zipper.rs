//! Pipeline entry point: run from zipper merge (unmapped BAM + mapped SAM/BAM).

use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use anyhow::{Context, Result, anyhow};
use fgumi_lib::bam_io::create_bam_reader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::template::{Template, TemplateIterator};
use fgumi_lib::vendored::encode_record_buf;
use fgumi_umi::TagInfo;
use log::info;
use noodles::sam::Header;

use super::Runall;
use crate::commands::zipper::{build_output_header, merge};

impl Runall {
    /// Run the pipeline starting from zipper (unmapped BAM + mapped SAM/BAM).
    ///
    /// Merges unmapped BAM (with UMI tags) and mapped SAM/BAM by queryname, then feeds
    /// merged records through sort -> group -> consensus -> filter. Both inputs must be
    /// queryname sorted or grouped.
    pub(super) fn run_from_zipper(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let timer = OperationTimer::new("Pipeline (from zipper)");
        let input = self.single_input()?;

        let unmapped_path = self.unmapped.as_ref().expect("--unmapped validated in validate()");
        let (_reference, dict_path) = self.require_reference()?;

        info!("Starting pipeline (from zipper)");
        info!("Mapped input:   {}", input.display());
        info!("Unmapped input: {}", unmapped_path.display());
        info!("Output: {}", self.output.display());

        // Open unmapped BAM in a reader thread (avoids lifetime issues with record_bufs)
        let (mut unmapped_reader, unmapped_header) = create_bam_reader(unmapped_path, 1)?;
        let unmapped_header_clone = unmapped_header.clone();
        let (unmapped_tx, unmapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(256);
        std::thread::spawn(move || {
            let iter = TemplateIterator::new(
                unmapped_reader
                    .record_bufs(&unmapped_header_clone)
                    .map(|r| r.map_err(anyhow::Error::from)),
            );
            for template in iter {
                if unmapped_tx.send(template).is_err() {
                    break;
                }
            }
        });
        let unmapped_iter = std::iter::from_fn(move || unmapped_rx.recv().ok());

        // Open mapped SAM/BAM in a reader thread
        let mapped_header: Header;
        let mapped_rx: std::sync::mpsc::Receiver<Result<Template>>;

        let ext = input.extension().and_then(|e| e.to_str()).unwrap_or("");
        if ext == "sam" {
            let file = std::fs::File::open(input).context("Failed to open mapped SAM")?;
            let buf_reader = BufReader::with_capacity(256 * 1024, file);
            let mut sam_reader = noodles::sam::io::Reader::new(buf_reader);
            mapped_header = sam_reader.read_header()?;
            let mh = mapped_header.clone();
            let (tx, rx) = std::sync::mpsc::sync_channel::<Result<Template>>(256);
            std::thread::spawn(move || {
                let iter = TemplateIterator::new(
                    sam_reader.record_bufs(&mh).map(|r| r.map_err(anyhow::Error::from)),
                );
                for template in iter {
                    if tx.send(template).is_err() {
                        break;
                    }
                }
            });
            mapped_rx = rx;
        } else {
            let (mut reader, header) = create_bam_reader(input, 1)?;
            mapped_header = header.clone();
            let (tx, rx) = std::sync::mpsc::sync_channel::<Result<Template>>(256);
            std::thread::spawn(move || {
                let iter = TemplateIterator::new(
                    reader.record_bufs(&header).map(|r| r.map_err(anyhow::Error::from)),
                );
                for template in iter {
                    if tx.send(template).is_err() {
                        break;
                    }
                }
            });
            mapped_rx = rx;
        }
        let mapped_iter: Box<dyn Iterator<Item = Result<Template>> + Send> =
            Box::new(std::iter::from_fn(move || mapped_rx.recv().ok()));

        // Build output header from unmapped + mapped + reference dict
        let output_header = build_output_header(&unmapped_header, &mapped_header, &dict_path)?;
        let output_header = crate::commands::common::add_pg_record(output_header, command_line)?;

        // Build tag info for zipper merge
        let tag_info = TagInfo::new(
            self.zipper_opts.zipper_tags_to_remove.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_reverse.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_revcomp.clone().unwrap_or_default(),
        );

        // Stream merged records through a bounded channel to avoid materializing all
        // records in memory. The producer thread runs the zipper merge; the consumer
        // feeds the sort/consensus pipeline.
        let (merge_tx, merge_rx) = crossbeam_channel::bounded::<Vec<u8>>(10_000);
        let merge_header = output_header.clone();
        let merge_tag_info = tag_info.clone();
        let merge_skip_pa = self.zipper_opts.zipper_skip_pa_tags;

        let producer = std::thread::Builder::new()
            .name("zipper-merge".into())
            .spawn(move || -> Result<u64> {
                zipper_merge_to_channel(
                    unmapped_iter,
                    mapped_iter,
                    &merge_header,
                    &merge_tag_info,
                    merge_skip_pa,
                    merge_tx,
                )
            })
            .context("Failed to spawn zipper merge thread")?;

        // Sort the streamed records and run group → consensus → filter.
        self.sort_and_process_records(
            merge_rx.iter(),
            &output_header,
            command_line,
            tmp_output,
            cancel,
            &timer,
            "from zipper",
        )?;

        // Wait for the producer to finish and propagate any errors.
        let record_count =
            producer.join().map_err(|_| anyhow!("Zipper merge thread panicked"))??;
        info!("Zipper merge complete, {record_count} raw records streamed");

        Ok(())
    }
}

// ============================================================================
// Zipper merge: merge unmapped + mapped templates into raw BAM bytes
// ============================================================================

/// Stream merged records from unmapped + mapped template iterators into a channel.
///
/// Iterates both inputs in queryname order (like `process_singlethreaded` in zipper.rs),
/// calls `merge()` on each matched pair, and serializes the resulting `RecordBuf`s to
/// raw BAM bytes, sending each one through the provided channel sender.
///
/// This avoids materializing all merged records in memory at once.
pub(super) fn zipper_merge_to_channel(
    unmapped_iter: impl Iterator<Item = Result<Template>>,
    mut mapped_iter: Box<dyn Iterator<Item = Result<Template>> + Send>,
    output_header: &Header,
    tag_info: &TagInfo,
    skip_pa_tags: bool,
    tx: crossbeam_channel::Sender<Vec<u8>>,
) -> Result<u64> {
    let progress = ProgressTracker::new("Zipper merged records").with_interval(1_000_000);
    let mut mapped_peek: Option<Template> = None;
    let mut record_count: u64 = 0;

    for unmapped_result in unmapped_iter {
        let unmapped_template = unmapped_result?;

        if mapped_peek.is_none() {
            mapped_peek = mapped_iter.next().transpose()?;
        }

        if let Some(ref mut mapped_template) = mapped_peek {
            if mapped_template.name == unmapped_template.name {
                // Merge unmapped metadata into the mapped template
                merge(&unmapped_template, mapped_template, tag_info, skip_pa_tags)?;

                // Serialize all records from the merged template to raw BAM bytes
                for rec in mapped_template.all_reads() {
                    let mut buf = Vec::new();
                    encode_record_buf(&mut buf, output_header, rec).map_err(|e| {
                        anyhow!("Failed to encode merged record to raw BAM bytes: {e}")
                    })?;
                    if tx.send(buf).is_err() {
                        // Receiver dropped — consumer encountered an error or was cancelled
                        anyhow::bail!("Zipper merge: consumer dropped the channel");
                    }
                    record_count += 1;
                    progress.log_if_needed(1);
                }

                mapped_peek = None;
            }
            // If names don't match, the unmapped read has no corresponding mapped read.
            // We skip orphan unmapped reads (like zipper with exclude_missing_reads=true).
        }
        // If no more mapped reads, remaining unmapped reads are skipped.
    }

    // Check for leftover mapped reads
    let remaining = mapped_peek.or_else(|| mapped_iter.next().transpose().ok().flatten());
    if let Some(remaining) = remaining {
        anyhow::bail!(
            "Zipper merge error: processed all unmapped reads but mapped reads remain. \
             Found template '{}'. Ensure unmapped and mapped inputs have matching read names.",
            String::from_utf8_lossy(&remaining.name)
        );
    }

    progress.log_final();
    Ok(record_count)
}
