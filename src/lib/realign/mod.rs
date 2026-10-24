//! Re-alignment of consensus reads against their molecule's raw-read window.
//!
//! [`window_group`] computes each molecule's raw-read window when grouping
//! (`fgumi group --realign-window`) and stamps it on the raw reads as the `rw`
//! tag, so a later re-alignment of the consensus reads can prefer the
//! alignment inside it.

pub(crate) mod window_group;
