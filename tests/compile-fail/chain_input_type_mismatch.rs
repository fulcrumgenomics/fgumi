//! `Chain<Single<u32>>::chain<S>` requires `S::Input = u32`. Chaining a
//! step whose `Input = u64` must fail to compile, not at runtime.

use std::io;

use fgumi_lib::pipeline::core::PipelineBuilder;
use fgumi_lib::pipeline::core::outputs::Single;
use fgumi_lib::pipeline::core::queues::QueueSpec;
use fgumi_lib::pipeline::core::reorder::BranchOrdering;
use fgumi_lib::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

#[derive(Clone)]
struct U32Source;
impl Step for U32Source {
    type Input = ();
    type Outputs = Single<u32>;
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "U32Source",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: 1 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }
    fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        Ok(StepOutcome::Finished)
    }
}

#[derive(Clone)]
struct U64Sink;
impl Step for U64Sink {
    type Input = u64; // mismatch with U32Source's Single<u32> output
    type Outputs = ();
    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "U64Sink",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }
    fn try_run(&mut self, _ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        Ok(StepOutcome::NoProgress)
    }
}

fn main() {
    let builder = PipelineBuilder::new();
    // Type error: cannot chain U64Sink (Input = u64) onto Chain<Single<u32>>.
    let _ = builder.chain(U32Source).chain(U64Sink);
}
