//! Bitmask type for SIMD-classified 64-byte blocks.

/// Bitmask output from lexing a 64-byte block.
///
/// Each bit position corresponds to a byte in the block.
/// Bit `i` is set if the byte at position `i` matches the target character.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FastqBitmask {
    /// Bit `i` = 1 if byte `i` is `\n` (0x0A).
    pub newlines: u64,
    /// Bit `i` = 1 if byte `i` is ACGT (case-insensitive).
    pub is_acgt: u64,
    /// 2-bit packed encoding of each byte position, using fgumi's `BitEnc` mapping:
    /// A/a=0b00, C/c=0b01, G/g=0b10, T/t=0b11.
    ///
    /// Bits `[2i, 2i+1]` encode byte position `i`. Only meaningful where `is_acgt`
    /// is set; non-ACGT positions contain garbage.
    pub two_bits: u128,
}
