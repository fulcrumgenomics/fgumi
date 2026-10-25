//! Re-alignment of consensus reads against their molecule's raw-read window.
//!
//! [`window_group`] computes each molecule's raw-read window when grouping
//! (`fgumi group --realign-window`) and stamps it on the raw reads as the `rw`
//! tag, which the consensus callers copy onto the consensus reads. After the
//! consensus reads are re-aligned, [`select`] prefers, per end, the alignment
//! that falls inside that window: bwa's primary if it does, else the first `XA`
//! alternative that does, else an in-window split (`SA`) piece of the same read,
//! which swaps roles with the primary. The two ends of a pair bwa flagged
//! proper move together or not at all: neither swap is made for one end unless
//! its mate is also outside its own window with an in-window `XA` hit or split
//! piece of its own (the proper-pair guard).

pub(crate) mod select;
pub(crate) mod window_group;
