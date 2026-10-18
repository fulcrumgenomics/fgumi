//! Per-worker output sizing for the consensus steps.
//!
//! A consensus batch body serializes every consensus record of a batch into one
//! output buffer. Starting that buffer empty makes it regrow by doubling — a
//! realloc and a copy of everything written so far, several times per batch.
//! [`OutputSizer`] reserves the buffer once instead, from the batch's input
//! size and the output/input ratio this worker has observed so far. It lives in
//! each worker's own state, so it needs no sharing and is correct at any worker
//! count.

use crate::pipeline::steps::types::fit_for_queue;

/// Largest output/input ratio the estimate may hold, in thousandths (4.0).
///
/// Consensus output is normally a fraction of its input. A batch that runs far
/// above it (tiny framed input, large tagged output) would otherwise push the
/// estimate, and the next batch's reservation, far past anything the step needs.
const MAX_RATIO_MILLI: u64 = 4_000;

/// Per-worker estimate of a consensus step's output/input byte ratio, used to
/// reserve the output buffer once per batch.
///
/// The output buffer becomes a `DecompressedBlock`, which its byte-bounded
/// queue charges by **capacity**, so the reservation is kept within config and
/// close to what the batch will hold:
///
/// - **No guess before evidence.** Until the worker has seen a batch there is
///   no estimate, and nothing is reserved: the first batch grows its buffer as
///   before. A seeded guess would over-reserve every worker's first batch.
/// - **Estimate.** After that, an exponential moving average (alpha = 1/2) of
///   the observed ratios, in thousandths, clamped to [`MAX_RATIO_MILLI`]. The
///   first observation sets it directly.
/// - **Reservation.** 5/4 of the estimated output plus 4 KiB of headroom, so a
///   batch whose ratio matches recent ones fits without regrowing, capped at
///   `ceiling` — the step's per-queue byte budget, i.e. config.
/// - **Shrink before queueing.** [`Self::finish`] observes the batch and fits
///   the buffer for its queue ([`fit_for_queue`]: capacity at most
///   `2 * len + 64 KiB`), so an estimate that overshoots — say, after one very
///   deep family — costs one shrink, never an inflated queue charge.
#[derive(Debug, Clone, Copy)]
pub(crate) struct OutputSizer {
    /// Estimated output bytes per 1000 input bytes; `None` until a batch with
    /// input has been observed.
    ratio_milli: Option<u64>,
    /// Largest reservation, in bytes.
    ceiling: usize,
    /// Capacity of the buffer after the most recent [`Self::reserve`].
    #[cfg(test)]
    pub(crate) last_reserved: usize,
}

impl OutputSizer {
    /// A sizer with no estimate yet whose reservations never exceed `ceiling`.
    pub(crate) fn new(ceiling: usize) -> Self {
        Self {
            ratio_milli: None,
            ceiling,
            #[cfg(test)]
            last_reserved: 0,
        }
    }

    /// Reserve `buf` for a batch of `input_len` framed input bytes. Reserves
    /// nothing until a batch has been observed.
    pub(crate) fn reserve(&mut self, buf: &mut Vec<u8>, input_len: usize) {
        if let Some(ratio_milli) = self.ratio_milli {
            let estimate =
                u64::try_from(input_len).unwrap_or(u64::MAX).saturating_mul(ratio_milli) / 1000;
            let want = usize::try_from(estimate.saturating_mul(5) / 4)
                .unwrap_or(usize::MAX)
                .saturating_add(4096)
                .min(self.ceiling);
            buf.reserve_exact(want);
        }
        #[cfg(test)]
        {
            self.last_reserved = buf.capacity();
        }
    }

    /// Fold a finished batch's `input_len` and `buf.len()` into the estimate,
    /// then fit `buf` for its queue.
    pub(crate) fn finish(&mut self, buf: &mut Vec<u8>, input_len: usize) {
        self.observe(input_len, buf.len());
        fit_for_queue(buf);
    }

    /// Fold one batch's observed `output_len / input_len` into the estimate.
    fn observe(&mut self, input_len: usize, output_len: usize) {
        if input_len == 0 {
            return;
        }
        let ratio = (u64::try_from(output_len).unwrap_or(u64::MAX).saturating_mul(1000)
            / u64::try_from(input_len).unwrap_or(u64::MAX))
        .min(MAX_RATIO_MILLI);
        self.ratio_milli = Some(match self.ratio_milli {
            None => ratio,
            Some(previous) => previous / 2 + ratio / 2,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::{MAX_RATIO_MILLI, OutputSizer};

    /// Nothing is reserved before the first observation.
    #[test]
    fn output_sizer_reserves_nothing_before_it_has_seen_a_batch() {
        let mut sizer = OutputSizer::new(usize::MAX);
        let mut buf = Vec::new();
        sizer.reserve(&mut buf, 1 << 20);
        assert_eq!(buf.capacity(), 0);
    }

    /// The first observation sets the estimate; later ones are averaged with
    /// alpha = 1/2: from 300 toward 100 the error after 20 steps is
    /// 200 / 2^20 < 1 (integer halving may leave 1).
    #[test]
    fn output_sizer_converges_with_alpha_one_half() {
        let mut sizer = OutputSizer::new(usize::MAX);
        sizer.observe(10_000, 3_000);
        assert_eq!(sizer.ratio_milli, Some(300));
        for _ in 0..20 {
            sizer.observe(10_000, 1_000);
        }
        let ratio = sizer.ratio_milli.expect("observed");
        assert!(ratio.abs_diff(100) <= 1, "{ratio}");
    }

    /// An empty batch carries no ratio and leaves the estimate unchanged.
    #[test]
    fn output_sizer_ignores_empty_input() {
        let mut sizer = OutputSizer::new(usize::MAX);
        sizer.observe(0, 5_000);
        assert_eq!(sizer.ratio_milli, None);
    }

    /// A batch whose output dwarfs its input is clamped to the ratio ceiling.
    #[test]
    fn output_sizer_clamps_the_ratio() {
        let mut sizer = OutputSizer::new(usize::MAX);
        sizer.observe(10, 1_000_000);
        assert_eq!(sizer.ratio_milli, Some(MAX_RATIO_MILLI));
    }

    /// A reservation never exceeds the ceiling, whatever the estimate.
    #[test]
    fn output_sizer_caps_the_reservation_at_the_ceiling() {
        let mut sizer = OutputSizer::new(64 * 1024);
        sizer.observe(1_000, 1_000);
        let mut buf = Vec::new();
        sizer.reserve(&mut buf, 100 << 20);
        assert_eq!(buf.capacity(), 64 * 1024);
    }

    /// `finish` shrinks a buffer the estimate over-reserved, so its queue is
    /// charged about what it holds.
    #[test]
    fn output_sizer_finish_fits_an_overshot_buffer() {
        let mut sizer = OutputSizer::new(usize::MAX);
        sizer.observe(1_000, 1_000);
        let mut buf = Vec::new();
        sizer.reserve(&mut buf, 1 << 20);
        buf.extend_from_slice(&[0; 1000]);
        sizer.finish(&mut buf, 1 << 20);
        assert!(buf.capacity() <= 2 * buf.len() + 64 * 1024, "{}", buf.capacity());
    }
}
