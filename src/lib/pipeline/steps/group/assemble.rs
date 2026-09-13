//! `AssembleTemplates` — the parallel queryname grouper.
//!
//! When the source cutter runs in `BatchCut::Queryname`
//! ([`FindBamBoundaries`](crate::pipeline::steps::boundaries::bam::FindBamBoundaries)
//! / [`ReadSamChunks`](crate::pipeline::steps::source::read_sam_chunks::ReadSamChunks)),
//! every [`DecodedRecordBatch`] is *closed under queryname*: no queryname run
//! straddles a batch boundary. Grouping records into [`Template`]s then needs no
//! cross-batch state — it is a pure per-batch function — so this step is
//! `StepKind::Parallel` and replaces the serial
//! [`GroupByQueryname`](super::queryname::GroupByQueryname) on the fast path.
//!
//! Contract: the input batch **must** be marked
//! [`closed_under_queryname`](DecodedRecordBatch::closed_under_queryname). The
//! flag is a producer-side promise (set by the decode step from the same
//! builder decision that configured the cutter) that this step cannot itself
//! verify — a batch cannot see its neighbours. So it fails closed on the first
//! batch whose flag is clear rather than silently splitting a straddling
//! template.

use std::io;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{BamTemplateBatch, DecodedRecordBatch};
use crate::template::Template;
use fgumi_bam_io::DecodedRecord;

/// Group one queryname-closed batch of decoded records into templates.
///
/// Run-length groups consecutive records by `key.name_hash` (with a read-name
/// byte compare on hash-equal, matching `GroupByQueryname::process_record`'s
/// collision guard), builds each run into a [`Template`] via
/// [`Template::from_records`], and returns the templates plus the summed
/// per-template [`Template::heap_size`] so the caller can build a
/// [`BamTemplateBatch`] without a second walk.
///
/// The batch is assumed closed under queryname: the first and last runs are
/// complete, so every run here becomes a whole template. This is checked by the
/// caller ([`AssembleTemplates::try_run`]) before this runs.
///
/// # Errors
///
/// Propagates [`Template::from_records`] validation errors (truncated record,
/// QNAME mismatch, multiple primary R1/R2).
pub(crate) fn assemble_closed(records: Vec<DecodedRecord>) -> io::Result<(Vec<Template>, usize)> {
    let mut templates: Vec<Template> = Vec::new();
    let mut total_heap: usize = 0;

    // Current run state.
    let mut run: Vec<fgumi_raw_bam::RawRecord> = Vec::new();
    let mut run_hash: u64 = 0;
    let mut run_name: Vec<u8> = Vec::new();

    let flush = |run: &mut Vec<fgumi_raw_bam::RawRecord>,
                 templates: &mut Vec<Template>,
                 total_heap: &mut usize|
     -> io::Result<()> {
        if run.is_empty() {
            return Ok(());
        }
        let t = Template::from_records(std::mem::take(run))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        *total_heap += t.heap_size();
        templates.push(t);
        Ok(())
    };

    for decoded in records {
        let name_hash = decoded.key.name_hash;
        let raw = decoded.into_raw_bytes();
        let read_name = fgumi_raw_bam::read_name(raw.as_ref());
        let same_run = !run.is_empty() && run_hash == name_hash && run_name == read_name;
        if same_run {
            run.push(raw);
        } else {
            flush(&mut run, &mut templates, &mut total_heap)?;
            run_hash = name_hash;
            run_name.clear();
            run_name.extend_from_slice(read_name);
            run.push(raw);
        }
    }
    flush(&mut run, &mut templates, &mut total_heap)?;

    Ok((templates, total_heap))
}

/// `Parallel + ByItemOrdinal` queryname grouper over
/// queryname-closed [`DecodedRecordBatch`]es. One output batch per input batch,
/// same `batch_serial`. No cross-batch state — each worker clone holds only its
/// held slot.
pub struct AssembleTemplates {
    held: HeldSlot<Unpushed<BamTemplateBatch>>,
    output_byte_limit: u64,
}

impl AssembleTemplates {
    #[must_use]
    pub fn new(output_byte_limit: u64) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for AssembleTemplates {
    fn clone(&self) -> Self {
        Self { held: HeldSlot::new(), output_byte_limit: self.output_byte_limit }
    }
}

impl Step for AssembleTemplates {
    type Input = DecodedRecordBatch;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "AssembleTemplates",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Drain the held slot first.
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
            // No input this call. On drain, this clone will never push again —
            // report Finished (the last Parallel clone closes the shared output,
            // gated by the driver's StepDrainCounter).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        // Fail closed: this step trusts the closure invariant it cannot verify.
        // A batch reaching here with the flag clear is a builder mis-wiring
        // (a queryname grouper fed by a non-queryname cutter); splitting a
        // straddling template would corrupt output silently, so error instead.
        if !batch.closed_under_queryname() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "AssembleTemplates received a batch not marked closed under queryname; \
                 it must be fed by a queryname-cut source (BatchCut::Queryname)",
            ));
        }

        let serial = batch.batch_serial();
        let (templates, total_heap) = assemble_closed(batch.into_records())?;
        let out = BamTemplateBatch::from_parts(serial, templates, total_heap);
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
    use fgumi_bam_io::library::LibraryIndex;
    use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags};

    const R1: u16 = flags::PAIRED | flags::FIRST_SEGMENT;
    const R2: u16 = flags::PAIRED | flags::LAST_SEGMENT;

    fn raw_named(name: &[u8], flag: u16) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name).sequence(b"ACGT").qualities(&[30; 4]).flags(flag);
        b.build()
    }

    fn decoded(name: &[u8], flag: u16) -> DecodedRecord {
        let raw = raw_named(name, flag);
        let key = fgumi_bam_io::GroupKey {
            name_hash: LibraryIndex::hash_name(Some(name)),
            ..fgumi_bam_io::GroupKey::default()
        };
        DecodedRecord::from_raw_bytes(raw, key)
    }

    /// `assemble_closed` groups consecutive same-name records into one template
    /// each and sums their heap size.
    #[test]
    fn assemble_closed_groups_runs() {
        let records = vec![
            decoded(b"read1", R1),
            decoded(b"read1", R2),
            decoded(b"read2", R1),
            decoded(b"read2", R2),
        ];
        let (templates, total) = assemble_closed(records).expect("assemble ok");
        let names: Vec<&[u8]> = templates.iter().map(|t| t.name.as_slice()).collect();
        assert_eq!(names, vec![b"read1".as_slice(), b"read2".as_slice()]);
        assert!(templates.iter().all(|t| t.records().len() == 2));
        let expect: usize = templates.iter().map(Template::heap_size).sum();
        assert_eq!(total, expect);
    }

    /// A single-record template (single-end) is its own run.
    #[test]
    fn assemble_closed_single_record_templates() {
        let records = vec![decoded(b"a", 0), decoded(b"b", 0), decoded(b"c", 0)];
        let (templates, _) = assemble_closed(records).expect("assemble ok");
        assert_eq!(templates.len(), 3);
    }

    /// Distinct non-adjacent recurrences of a name are separate templates (the
    /// closed batch never merges non-adjacent runs).
    #[test]
    fn assemble_closed_nonadjacent_names_are_separate() {
        let records = vec![decoded(b"x", R1), decoded(b"y", R1), decoded(b"x", R1)];
        let (templates, _) = assemble_closed(records).expect("assemble ok");
        let names: Vec<&[u8]> = templates.iter().map(|t| t.name.as_slice()).collect();
        assert_eq!(names, vec![b"x".as_slice(), b"y".as_slice(), b"x".as_slice()]);
    }

    /// An empty batch assembles to no templates.
    #[test]
    fn assemble_closed_empty() {
        let (templates, total) = assemble_closed(Vec::new()).expect("assemble ok");
        assert!(templates.is_empty());
        assert_eq!(total, 0);
    }

    #[test]
    fn profile_is_parallel_byitemordinal() {
        let s = AssembleTemplates::new(1 << 20);
        let p = s.profile();
        assert_eq!(p.name, "AssembleTemplates");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    }
}
