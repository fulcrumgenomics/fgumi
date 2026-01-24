//! Encoder and Decoder implementing `std::io::Read`
use crate::*;
use std::io;

/// Streaming compression for input streams implementing `std::io::Read`.
///
/// Notes
/// -----
/// One should consider using `crate::compress` or `crate::compress_into` if possible.
/// In that context, we do not need to hold and maintain intermediate buffers for reading and writing.
///
/// Example
/// -------
/// ```
/// use std::{io, io::Read};
/// use isal::{read::Encoder, CompressionLevel, decompress, Codec};
/// let data = b"Hello, World!".to_vec();
///
/// let mut encoder = Encoder::new(data.as_slice(), CompressionLevel::Three, Codec::Gzip);
/// let mut compressed = vec![];
///
/// // Number of compressed bytes written to `output`
/// let n = io::copy(&mut encoder, &mut compressed).unwrap();
/// assert_eq!(n as usize, compressed.len());
///
/// let decompressed = decompress(compressed.as_slice(), Codec::Gzip).unwrap();
/// assert_eq!(decompressed.as_slice(), data);
/// ```
pub struct Encoder<R: io::Read> {
    inner: R,
    stream: ZStream,
    in_buf: [u8; BUF_SIZE],
}

impl<R: io::Read> Encoder<R> {
    /// Create a new `Encoder` which implements the `std::io::Read` trait.
    pub fn new(reader: R, level: CompressionLevel, codec: Codec) -> Encoder<R> {
        let in_buf = [0_u8; BUF_SIZE];

        let mut zstream = ZStream::new(level, ZStreamKind::Stateful);
        zstream.stream.end_of_stream = 0;
        zstream.stream.flush = FlushFlags::NoFlush as _;
        zstream.stream.gzip_flag = codec as _;

        Self {
            inner: reader,
            stream: zstream,
            in_buf,
        }
    }

    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner
    }
}

impl<R: io::Read> io::Read for Encoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.stream.stream.avail_out = buf.len() as _;
        self.stream.stream.next_out = buf.as_mut_ptr();

        // If we have an intermediate compressed format (ICF) buffer full, we just need
        // to provide more output and return
        if self.stream.stream.internal_state.state
            == isal::isal_zstate_state_ZSTATE_TMP_FLUSH_ICF_BUFFER
        {
            self.stream.deflate()?;
            return Ok(buf.len() - self.stream.stream.avail_out as usize);
        }

        let mut nbytes = 0;
        while self.stream.stream.avail_out > 0 {
            // Read out next buf len worth to compress; filling intermediate out_buf
            // Check if current next_in ptr is valid; pointing to somewhere in the in buffer
            // it may not be after a reset. If not, just read into in_buf and update
            // otherwise we can shift it to the front and read more w/o reallocation
            let next_in_ptr = self.stream.stream.next_in;
            let base_ptr = self.in_buf.as_mut_ptr();
            if next_in_ptr.is_null()
                || next_in_ptr < base_ptr
                || next_in_ptr > unsafe { base_ptr.add(self.in_buf.len()) }
            {
                self.stream.stream.avail_in = self.inner.read(&mut self.in_buf)? as _;
                self.stream.stream.next_in = self.in_buf.as_mut_ptr();
                self.stream.stream.end_of_stream = (self.stream.stream.avail_in == 0) as _;
            } else {
                let idx = unsafe { next_in_ptr.offset_from(base_ptr) };

                // Extra sanity check - this probably / should never happen after ptr check above.
                if idx < 0
                    || idx + self.stream.stream.avail_in as isize > self.in_buf.len() as isize
                {
                    let msg = format!(
                            "Compression error refilling in buffer, index is probably wrong: {}, avail_in: {}, avail_out: {}",
                            idx, self.stream.stream.avail_in, self.stream.stream.avail_out
                        );
                    let err = Error::Other((None, msg));
                    return Err(io::Error::new(io::ErrorKind::Other, err));
                }

                // Can safely shift next input to the front and read more into tail of in buf
                let idx = idx as usize;
                self.in_buf
                    .copy_within(idx..idx + self.stream.stream.avail_in as usize, 0);
                let avail_in = self
                    .inner
                    .read(&mut self.in_buf[self.stream.stream.avail_in as usize..])?
                    as usize;
                self.stream.stream.avail_in += avail_in as u32;
                self.stream.stream.end_of_stream = (avail_in == 0) as _;
                self.stream.stream.next_in = self.in_buf[0..].as_mut_ptr();
            }

            self.stream.deflate()?;
            nbytes = buf.len() - self.stream.stream.avail_out as usize;

            if self.stream.stream.internal_state.state == isal::isal_zstate_state_ZSTATE_END {
                break;
            }
        }

        Ok(nbytes)
    }
}

/// Streaming compression for input streams implementing `std::io::Read`.
///
/// Notes
/// -----
/// One should consider using `crate::decompress` or `crate::decompress_into` if possible.
/// In that context, we do not need to hold and maintain intermediate buffers for reading and writing.
///
/// Example
/// -------
/// ```
/// use std::{io, io::Read};
/// use isal::{read::Decoder, CompressionLevel, compress, Codec};
/// let data = b"Hello, World!".to_vec();
///
/// let compressed = compress(data.as_slice(), CompressionLevel::Three, Codec::Gzip).unwrap();
/// let mut decoder = Decoder::new(compressed.as_slice(), Codec::Gzip);
/// let mut decompressed = vec![];
///
/// // Numbeer of compressed bytes written to `output`
/// let n = io::copy(&mut decoder, &mut decompressed).unwrap();
/// assert_eq!(n as usize, data.len());
/// assert_eq!(decompressed.as_slice(), data);
/// ```
pub struct Decoder<R: io::Read> {
    inner: R,
    zst: InflateState,
    in_buf: [u8; BUF_SIZE],
    #[allow(dead_code)]
    codec: Codec,
}

impl<R: io::Read> Decoder<R> {
    pub fn new(reader: R, codec: Codec) -> Decoder<R> {
        let zst = InflateState::new(codec);

        Self {
            inner: reader,
            zst,
            in_buf: [0u8; BUF_SIZE],
            codec,
        }
    }

    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner
    }
}

impl<R: io::Read> io::Read for Decoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.zst.state.next_out = buf.as_mut_ptr();
        self.zst.state.avail_out = buf.len() as _;

        // keep writing as much as possible to the output buf
        let mut n_bytes = 0;
        while self.zst.state.avail_out > 0 {
            // Keep in_buf full for avail_in
            if self.zst.state.avail_in < self.in_buf.len() as _ {
                // if null, it's our first assignment of compressed data in
                if self.zst.state.next_in.is_null() {
                    self.zst.state.avail_in = self.inner.read(&mut self.in_buf)? as u32;
                    self.zst.state.next_in = self.in_buf.as_mut_ptr();

                // Otherwise shift current available in section to the front and read more
                } else {
                    let idx = unsafe {
                        self.zst.state.next_in.offset_from(self.in_buf.as_ptr()) as usize
                    };
                    self.in_buf
                        .copy_within(idx..idx + self.zst.state.avail_in as usize, 0);
                    self.zst.state.next_in = self.in_buf.as_mut_ptr();

                    // read some more and update avail_in
                    let n = self
                        .inner
                        .read(&mut self.in_buf[self.zst.state.avail_in as usize..])?;
                    self.zst.state.avail_in += n as u32;
                };
            }
            // Block finished, reset if we have more compressed data, otherwise exit
            if self.zst.block_state() == isal::isal_block_state_ISAL_BLOCK_FINISH {
                if self.zst.state.avail_in == 0 {
                    break;
                }
                self.zst.reset();
            }

            self.zst.step_inflate()?;
            n_bytes = buf.len() - self.zst.state.avail_out as usize;
        }

        Ok(n_bytes)
    }
}

/// Deflate decompression
/// Basically a wrapper to `Encoder` which sets the codec for you.
pub struct DeflateEncoder<R: io::Read> {
    inner: Encoder<R>,
}

impl<R: io::Read> DeflateEncoder<R> {
    pub fn new(reader: R, level: CompressionLevel) -> Self {
        Self {
            inner: Encoder::new(reader, level, Codec::Deflate),
        }
    }
    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner.inner
    }
}

impl<R: io::Read> io::Read for DeflateEncoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

/// Deflate decompression
/// Basically a wrapper to `Decoder` which sets the codec for you.
pub struct DeflateDecoder<R: io::Read> {
    inner: Decoder<R>,
}

impl<R: io::Read> DeflateDecoder<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: Decoder::new(reader, Codec::Deflate),
        }
    }
    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner.inner
    }
}

impl<R: io::Read> io::Read for DeflateDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}
/// Zlib decompression
/// Basically a wrapper to `Encoder` which sets the codec for you.
pub struct ZlibEncoder<R: io::Read> {
    inner: Encoder<R>,
}

impl<R: io::Read> ZlibEncoder<R> {
    pub fn new(reader: R, level: CompressionLevel) -> Self {
        Self {
            inner: Encoder::new(reader, level, Codec::Zlib),
        }
    }
    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner.inner
    }
}

impl<R: io::Read> io::Read for ZlibEncoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

/// Zlib decompression
/// Basically a wrapper to `Decoder` which sets the codec for you.
pub struct ZlibDecoder<R: io::Read> {
    inner: Decoder<R>,
}

impl<R: io::Read> ZlibDecoder<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: Decoder::new(reader, Codec::Zlib),
        }
    }
    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner.inner
    }
}

impl<R: io::Read> io::Read for ZlibDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}
/// Gzip decompression
/// Basically a wrapper to `Encoder` which sets the codec for you.
pub struct GzipEncoder<R: io::Read> {
    inner: Encoder<R>,
}

impl<R: io::Read> GzipEncoder<R> {
    pub fn new(reader: R, level: CompressionLevel) -> Self {
        Self {
            inner: Encoder::new(reader, level, Codec::Gzip),
        }
    }
    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner.inner
    }
}

impl<R: io::Read> io::Read for GzipEncoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

/// Gzip decompression
/// Basically a wrapper to `Decoder` which sets the codec for you.
pub struct GzipDecoder<R: io::Read> {
    inner: Decoder<R>,
}

impl<R: io::Read> GzipDecoder<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: Decoder::new(reader, Codec::Gzip),
        }
    }
    /// Mutable reference to underlying reader, not advisable to modify during reading.
    pub fn get_ref_mut(&mut self) -> &mut R {
        &mut self.inner.inner
    }

    // Reference to underlying reader
    pub fn get_ref(&self) -> &R {
        &self.inner.inner
    }
}

impl<R: io::Read> io::Read for GzipDecoder<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}
