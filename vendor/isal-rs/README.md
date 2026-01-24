# isal-rs

[![CI](https://github.com/milesgranger/isal-rs/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/milesgranger/isal-rs/actions/workflows/CI.yml)
[![Latest version](https://img.shields.io/crates/v/isal-rs.svg)](https://crates.io/crates/isal-rs)
[![Documentation](https://docs.rs/isal-rs/badge.svg)](https://docs.rs/isal-rs)
![License](https://img.shields.io/crates/l/isal-rs.svg)

Rust bindings for [isa-l](https://github.com/intel/isa-l)

---

Supports the following codecs using the ISA-L library under the hood:

- GZIP 
  - `isal::read::GzipEncoder/GzipDecoder`
  - `isal::write::GzipEncoder/GzipDecoder`
- DEFLATE
  - `isal::read::DeflateEncoder/DeflateDecoder`
  - `isal::write::DeflateEncoder/DeflateDecoder`
- ZLIB
  - `isal::read::ZlibEncoder/ZlibDecoder`
  - `isal::write::ZlibEncoder/ZlibDecoder`

Or can use functions of `de/compress` and `de/compress_into`

---

Building requires some system tools like `autotools`, `nasm`, `make`, and anything the official ISA-L repo suggests. 
On Windows the build is invoked with `nmake`, other systems use the `./autogen.sh` and `./configure` setups.

---

### Examples:

#### Functions like `compress_into` and `decompress`
(Similar functionality with `compress` and `decompress_into`)
```rust
use isal::{CompressionLevel, Codec, compress_into, decompress};

let mut compressed = vec![0u8; 100];
let nbytes = compress_into(b"foobar", &mut compressed, CompressionLevel::Three, Codec::Gzip).unwrap();

let decompressed = decompress(&compressed[..nbytes], Codec::Gzip).unwrap();
assert_eq!(decompressed.as_slice(), b"foobar");
```

#### Compress using the Encoder/Decoder structs implementing `io::Read`

```rust
use std::{io, io::Read};
use isal::{read::{Encoder, GzipEncoder}, CompressionLevel, decompress, Codec};

let data = b"Hello, World!".to_vec();

// Note these two encoders are equivelent...
let mut encoder = GzipEncoder::new(data.as_slice(), CompressionLevel::Three);
let mut encoder = Encoder::new(data.as_slice(), CompressionLevel::Three, Codec::Gzip);

// Number of compressed bytes written to `output`
let mut compressed = vec![];
let n = io::copy(&mut encoder, &mut compressed).unwrap();
assert_eq!(n as usize, compressed.len());

let decompressed = decompress(compressed.as_slice(), Codec::Gzip).unwrap();
assert_eq!(decompressed.as_slice(), data);
```

#### Decompress using the Encoder/Decoder structs implementing `io::Write`

```rust
use std::{io, io::Write};
use isal::{write::Decoder, CompressionLevel, compress, Codec};

let data = b"Hello, World!".to_vec();
let compressed = compress(data.as_slice(), CompressionLevel::Three, Codec::Gzip).unwrap();

let mut decompressed = vec![];
let mut decoder = Decoder::new(&mut decompressed, Codec::Gzip);

// Number of compressed bytes written to `output`
let n = io::copy(&mut io::Cursor::new(&compressed), &mut decoder).unwrap();
assert_eq!(n as usize, compressed.len());

assert_eq!(decompressed.as_slice(), data);
```

---

### Benchmarks

TL/DR: It's roughly 5-10x faster than flate2 with the default features, and
~2-3x faster using flate2 with zlib-ng backend.

Checkout the [README](./benches/README.md) in the benches directory.
Criterion benchmark report available here: https://milesgranger.github.io/isal-rs/benches/criterion/report/

_NOTE:_ ISA-L supports compression levels 0, 1, 3. These benchmarks compare gzip using compression level 3
for both flate2 and isal-rs.

---

### Versioning: 
Versions are specified in normal SemVer format, and a trailing "`+<< commit hash >>`" to indicate
which commit in [isa-l](https://github.com/intel/isa-l) the crate is built against. ie: `0.1.0+62519d9`
