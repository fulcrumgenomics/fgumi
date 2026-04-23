#![deny(unsafe_code)]

//! Validated SAM/BAM aux-tag identifier — see [`SamTag`].
//!
//! This crate exists to break a dependency cycle between `fgumi-sam` and
//! `fgumi-raw-bam`. Both need access to `SamTag`; putting it here keeps
//! either crate from depending on the other.

mod tag;
pub use tag::SamTag;

/// Convert a value to a reference to a two-byte SAM/BAM tag identifier.
///
/// Implemented for [`SamTag`], `[u8; 2]`, and `&[u8; 2]`. Use as a bound
/// (`tag: impl AsTagBytes`) on tag-parameter functions so that both
/// `SamTag::RX` and `b"RX"` literals can be passed without conversion.
pub trait AsTagBytes {
    /// Return a reference to the underlying two-byte tag.
    fn as_tag_bytes(&self) -> &[u8; 2];
}

impl AsTagBytes for [u8; 2] {
    #[inline]
    fn as_tag_bytes(&self) -> &[u8; 2] {
        self
    }
}

impl AsTagBytes for SamTag {
    #[inline]
    fn as_tag_bytes(&self) -> &[u8; 2] {
        self.as_ref()
    }
}

impl<T: AsTagBytes + ?Sized> AsTagBytes for &T {
    #[inline]
    fn as_tag_bytes(&self) -> &[u8; 2] {
        T::as_tag_bytes(self)
    }
}
