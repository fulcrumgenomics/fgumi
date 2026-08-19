//! Per-record cost of the template-coordinate ingest path.
//!
//! Phase 1 of a spill-heavy sort is bound by its serial main thread: on
//! `1kg-wgs-HG00096` at 16 threads that thread is busy 219.5s of a 240.7s
//! phase, and `perf` puts 39.5% of it in `extract_template_key_inline` with a
//! further 12.3% in the slice-iterator `next()` the aux-tag scan walks bytes
//! with. This bench sizes that function directly, so a change to it is measured
//! before a seven-minute whole-genome run is spent on it.
//!
//! The records are built to the aux layout the measured sample actually carries
//! (`PG:Z AS:i XS:i MD:Z NM:i RG:Z MQ:i MC:Z`, 118 bytes) rather than a minimal
//! one: the scan's cost is per aux byte, so a record with two tags would report
//! a number that has nothing to do with the workload. `with_xa` adds the long
//! `XA:Z` alt-hit tag that a fraction of records carry, since it roughly doubles
//! the aux data and is where the scan's tail lives.

use bstr::BString;
use criterion::{Criterion, criterion_group, criterion_main};
use fgumi_raw_bam::SamTag;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use std::hint::black_box;

const READ_NAME: &[u8] = b"A00132:53:HFHJKDSXX:1:1646:26467:33332\0";
const RG_ID: &str = "HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L001";
const SEQ_LEN: usize = 151;

/// Append a `Z`-typed aux tag (two tag bytes, `Z`, value, NUL).
fn push_z(aux: &mut Vec<u8>, tag: SamTag, value: &[u8]) {
    aux.extend_from_slice(&*tag);
    aux.push(b'Z');
    aux.extend_from_slice(value);
    aux.push(0);
}

/// Append an `i`-typed (32-bit signed) aux tag.
fn push_i(aux: &mut Vec<u8>, tag: SamTag, value: i32) {
    aux.extend_from_slice(&*tag);
    aux.push(b'i');
    aux.extend_from_slice(&value.to_le_bytes());
}

/// One BAM record body (from `ref_id`, i.e. without the `block_size` prefix),
/// carrying the aux layout of the measured 1000 Genomes WGS sample.
fn record(pos: i32, mate_pos: i32, with_xa: bool) -> Vec<u8> {
    let mut aux = Vec::with_capacity(256);
    push_z(&mut aux, SamTag::PG, b"MarkDuplicates");
    push_i(&mut aux, SamTag::AS, 64);
    push_i(&mut aux, SamTag::XS, 61);
    push_z(&mut aux, SamTag::MD, b"0N0N0N0N0N2A61");
    push_i(&mut aux, SamTag::NM, 6);
    push_z(&mut aux, SamTag::RG, RG_ID.as_bytes());
    push_i(&mut aux, SamTag::MQ, 0);
    push_z(&mut aux, SamTag::MC, b"81S69M");
    if with_xa {
        push_z(
            &mut aux,
            SamTag::new(b'X', b'A'),
            b"chr3,+198173832,34M116S,0;chr12,-108091,117S33M,0;chr1,-180805,118S32M,0;\
              chr1,-10052,118S32M,0;chr3,-10519,118S32M,0;",
        );
    }

    let n_cigar_op: u16 = 2;
    let mut rec = Vec::with_capacity(32 + READ_NAME.len() + 8 + SEQ_LEN * 2 + aux.len());
    rec.extend_from_slice(&0i32.to_le_bytes()); // ref_id
    rec.extend_from_slice(&pos.to_le_bytes()); // pos
    rec.push(u8::try_from(READ_NAME.len()).expect("read name fits a u8"));
    rec.push(60); // mapq
    rec.extend_from_slice(&0u16.to_le_bytes()); // bin
    rec.extend_from_slice(&n_cigar_op.to_le_bytes());
    // PAIRED | PROPER_PAIR | REVERSE | FIRST_SEGMENT. The mate is forward, so
    // the mate lane resolves through `unclipped_other_start` (leading clips).
    rec.extend_from_slice(&0x0053u16.to_le_bytes());
    rec.extend_from_slice(&i32::try_from(SEQ_LEN).expect("seq len fits").to_le_bytes());
    rec.extend_from_slice(&0i32.to_le_bytes()); // next_ref_id
    rec.extend_from_slice(&mate_pos.to_le_bytes());
    rec.extend_from_slice(&0i32.to_le_bytes()); // tlen
    rec.extend_from_slice(READ_NAME);
    // CIGAR 81S69M: (len << 4) | op, S = 4, M = 0.
    rec.extend_from_slice(&((81u32 << 4) | 4).to_le_bytes());
    rec.extend_from_slice(&(69u32 << 4).to_le_bytes());
    rec.resize(rec.len() + SEQ_LEN.div_ceil(2), 0x11); // packed seq
    rec.resize(rec.len() + SEQ_LEN, 30); // qual
    rec.extend_from_slice(&aux);
    rec
}

/// A header declaring the sample's single read group, so the RG lookup resolves
/// rather than missing (a miss and a hit take different paths through the map).
fn header() -> Header {
    let read_group = Map::<ReadGroup>::builder()
        .insert(rg_tag::LIBRARY, String::from("lib1"))
        .build()
        .expect("a read group with only LB is valid");
    Header::builder().add_read_group(BString::from(RG_ID), read_group).build()
}

fn bench_template_key(c: &mut Criterion) {
    let lookup = fgumi_sort::LibraryLookup::from_header(&header());
    let hasher = fgumi_sort::cb_hasher();

    for (label, with_xa) in [("aux-122b", false), ("aux-with-xa", true)] {
        // A batch rather than one record: the scan is memory-bound over the aux
        // bytes, and re-reading one cache-hot record would flatter it.
        let records: Vec<Vec<u8>> =
            (0..1024).map(|i| record(10_000 + i, 133_000_000 + i, with_xa)).collect();
        let aux_bytes: usize =
            records.iter().map(|r| fgumi_raw_bam::aux_data_slice(r).len()).sum::<usize>()
                / records.len();

        let mut group = c.benchmark_group("extract_template_key_inline");
        group.throughput(criterion::Throughput::Elements(records.len() as u64));
        group.bench_function(format!("{label}-mean-aux-{aux_bytes}"), |b| {
            b.iter(|| {
                for rec in &records {
                    let key = fgumi_sort::extract_template_key_inline(
                        black_box(rec.as_slice()),
                        &lookup,
                        None,
                        &hasher,
                    );
                    black_box(&key);
                }
            });
        });
        group.finish();
    }
}

criterion_group!(benches, bench_template_key);
criterion_main!(benches);
