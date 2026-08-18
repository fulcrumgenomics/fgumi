//! Compares the two `read_raw_record` body-read strategies head-to-head over
//! the SAME in-memory record stream, in one bench binary:
//!
//! * `zerofill` -- the OLD strategy: `resize(block_size, 0)` then
//!   `read_exact` (the memset PR1 removed).
//! * `spare_capacity` -- the NEW strategy: `clear()` + `reserve(block_size)`
//!   then a spare-capacity `read_exact` + `set_len` (what `read_raw_record`
//!   does today).
//!
//! Both arms read from a `Cursor` over the identical byte buffer, so
//! criterion isolates the single-line difference with no cross-build /
//! page-cache confound. The `zerofill` arm reconstructs the pre-change body
//! against a plain, locally-owned `Vec<u8>` -- the library's `RawRecord`
//! inner `Vec` is private outside the crate -- but calls the same public
//! `read_block_size` helper and does the identical `resize`/`read_exact`
//! pair `read_raw_record` used to do, over the same bytes. The
//! `spare_capacity` arm calls the real, current `read_raw_record` directly,
//! which is the most faithful "new" measurement.
//!
//! Run with: `cargo bench -p fgumi-raw-bam --bench record_read`
//! Quick sanity: `cargo bench -p fgumi-raw-bam --bench record_read -- --quick`

// Benches are a separate crate from the library; the library itself stays
// `#![deny(unsafe_code)]`. Nothing below currently needs `unsafe`, but this
// mirrors the brief's isolation of any bench-local unsafe spare-capacity
// experiments to this file only, never the library.
#![allow(unsafe_code)]

use std::hint::black_box;
use std::io::{self, Cursor, Read};

use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_raw_bam::{RawRecord, SamTag, UnmappedSamBuilder, read_block_size, read_raw_record};

/// Number of records in the synthetic stream. Large enough to amortize
/// per-iteration fixed costs (cursor/buffer setup) and to exercise many
/// resize/reserve cycles on the reused buffer.
const NUM_RECORDS: usize = 200_000;

/// Realistic short-read body sizes to cycle through, so the reused buffer's
/// resize/reserve path sees both growing and shrinking records rather than a
/// single fixed size.
const SEQ_LENS: [usize; 3] = [100, 150, 250];

/// Build a stream of `num_records` `block_size`-prefixed BAM records back to
/// back, each with a packed sequence + quality of a cycling length plus two
/// aux tags (`RG:Z`, `NM:i`) -- a realistic short-read record body of a few
/// hundred bytes.
fn make_record_stream(num_records: usize) -> Vec<u8> {
    let mut builder = UnmappedSamBuilder::new();
    let mut buf = Vec::new();
    for i in 0..num_records {
        let len = SEQ_LENS[i % SEQ_LENS.len()];
        let bases: Vec<u8> = (0..len).map(|j| b"ACGT"[j % 4]).collect();
        let quals = vec![30u8; len];
        let name = format!("read:{i}");
        builder.build_record(name.as_bytes(), 0, &bases, &quals);
        builder.append_string_tag(SamTag::RG, b"sample1");
        builder.append_int_tag(SamTag::NM, 0);
        builder.write_with_block_size(&mut buf);
    }
    buf
}

/// The OLD `read_raw_record` body-read strategy: `resize(block_size, 0)`
/// (zero-fill) then `read_exact` overwrites it. See the module doc for why
/// this operates on a plain `Vec<u8>` rather than a `RawRecord`.
fn read_zerofill<R: Read>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize> {
    let block_size = match read_block_size(reader)? {
        0 => return Ok(0),
        n => n,
    };
    buf.resize(block_size, 0);
    reader.read_exact(buf)?;
    Ok(block_size)
}

fn bench_record_read(c: &mut Criterion) {
    let stream = make_record_stream(NUM_RECORDS);

    let mut group = c.benchmark_group("record_read");

    group.bench_function("zerofill", |b| {
        b.iter(|| {
            let mut cursor = Cursor::new(stream.as_slice());
            let mut buf = Vec::new();
            let mut count = 0u64;
            loop {
                let n = read_zerofill(&mut cursor, &mut buf).expect("read should succeed");
                if n == 0 {
                    break;
                }
                count += 1;
            }
            black_box(count)
        });
    });

    group.bench_function("spare_capacity", |b| {
        b.iter(|| {
            let mut cursor = Cursor::new(stream.as_slice());
            let mut record = RawRecord::new();
            let mut count = 0u64;
            loop {
                let n = read_raw_record(&mut cursor, &mut record).expect("read should succeed");
                if n == 0 {
                    break;
                }
                count += 1;
            }
            black_box(count)
        });
    });

    group.finish();
}

criterion_group!(benches, bench_record_read);
criterion_main!(benches);
