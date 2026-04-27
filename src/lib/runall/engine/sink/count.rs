//! Drops items. Used when downstream of `CountStage` (which produces `()`).

use anyhow::Result;

use crate::runall::engine::backoff::Backoff;
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::sink::{InputQueue, Sink};

pub struct CountSink;

impl Sink for CountSink {
    type Input = ();

    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let mut backoff = Backoff::new();
        loop {
            if cancel.is_cancelled() {
                return Ok(());
            }
            if input.pop().is_some() {
                backoff.reset();
            } else {
                if input.is_drained() {
                    return Ok(());
                }
                backoff.snooze();
            }
        }
    }

    fn name(&self) -> &'static str {
        "CountSink"
    }
}
