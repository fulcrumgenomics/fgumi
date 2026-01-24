use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use isal as igzip;
use std::io::{self, Cursor};

static DATA: [(&'static str, &'static [u8]); 25] = [
    ("COPYING", include_bytes!("data/COPYING")),
    ("alice29.txt", include_bytes!("data/alice29.txt")),
    ("asyoulik.txt", include_bytes!("data/asyoulik.txt")),
    ("fireworks.jpeg", include_bytes!("data/fireworks.jpeg")),
    ("geo.protodata", include_bytes!("data/geo.protodata")),
    ("html", include_bytes!("data/html")),
    ("html_x_4", include_bytes!("data/html_x_4")),
    ("kppkn.gtb", include_bytes!("data/kppkn.gtb")),
    ("lcet10.txt", include_bytes!("data/lcet10.txt")),
    (
        "Mark.Twain-Tom.Sawyer.txt",
        include_bytes!("data/Mark.Twain-Tom.Sawyer.txt"),
    ),
    ("mr", include_bytes!("data/mr")),
    ("ooffice", include_bytes!("data/ooffice")),
    ("osdb", include_bytes!("data/osdb")),
    ("paper-100k.pdf", include_bytes!("data/paper-100k.pdf")),
    ("plrabn12.txt", include_bytes!("data/plrabn12.txt")),
    ("sao", include_bytes!("data/sao")),
    ("urls.10K", include_bytes!("data/urls.10K")),
    ("xml", include_bytes!("data/xml")),
    ("x-ray", include_bytes!("data/x-ray")),
    ("reymont", include_bytes!("data/reymont")),
    ("dickens", include_bytes!("data/dickens")),
    ("samba", include_bytes!("data/samba")),
    ("nci", include_bytes!("data/nci")),
    ("mozilla", include_bytes!("data/mozilla")),
    ("webster", include_bytes!("data/webster")),
];

// Get a subset of the data, using all of it takes a very long time
// so can mix-match when in doubt or debugging. As of now, it selects
// a decent mix I think, last one being a larger file 'webster' @ ~40M
static SAMPLE: [(&'static str, &'static [u8]); 6] =
    [DATA[1], DATA[5], DATA[10], DATA[15], DATA[20], DATA[24]];

fn roundtrip(c: &mut Criterion) {
    for (name, data) in SAMPLE {
        let mut group = c.benchmark_group("roundtrip");

        group.bench_with_input(BenchmarkId::new("isal", name), data, |b, data| {
            b.iter(|| {
                let compressed =
                    isal::compress(data, isal::CompressionLevel::Three, isal::Codec::Gzip).unwrap();
                let decompressed =
                    isal::decompress(compressed.as_slice(), isal::Codec::Gzip).unwrap();
                assert_eq!(data.len(), decompressed.len());
            })
        });
        group.bench_with_input(BenchmarkId::new("flate2", name), data, |b, data| {
            b.iter(|| {
                let mut compressed = vec![];
                let mut encoder = flate2::read::GzEncoder::new(data, flate2::Compression::new(3));
                let _ = io::copy(&mut encoder, &mut compressed).unwrap();

                let mut decompressed = vec![];
                let mut decoder = flate2::read::GzDecoder::new(compressed.as_slice());
                let _ = io::copy(&mut decoder, &mut decompressed).unwrap();

                assert_eq!(data.len(), decompressed.len());
            })
        });
        group.finish()
    }
}

fn read_encoder(c: &mut Criterion) {
    for (name, data) in SAMPLE {
        let mut group = c.benchmark_group("io::Read Gzip Encoder");

        group.bench_with_input(BenchmarkId::new("isal", name), data, |b, data| {
            b.iter(|| {
                let mut output = vec![];
                let mut encoder =
                    igzip::read::GzipEncoder::new(data, igzip::CompressionLevel::Three);
                let _ = io::copy(&mut encoder, &mut output).unwrap();
            })
        });
        group.bench_with_input(BenchmarkId::new("flate2", name), data, |b, data| {
            b.iter(|| {
                let mut output = vec![];
                let mut encoder = flate2::read::GzEncoder::new(data, flate2::Compression::new(3));
                let _ = io::copy(&mut encoder, &mut output).unwrap();
            })
        });
        group.finish()
    }
}

fn write_encoder(c: &mut Criterion) {
    for (name, data) in SAMPLE {
        let mut group = c.benchmark_group("io::Write Grzip Encoder");

        group.bench_with_input(BenchmarkId::new("isal", name), data, |b, data| {
            b.iter(|| {
                let mut output = vec![];
                let mut encoder =
                    igzip::write::GzipEncoder::new(&mut output, igzip::CompressionLevel::Three);
                let _ = io::copy(&mut Cursor::new(&data), &mut encoder).unwrap();
            })
        });
        group.bench_with_input(BenchmarkId::new("flate2", name), data, |b, data| {
            b.iter(|| {
                let mut output = vec![];
                let mut encoder =
                    flate2::write::GzEncoder::new(&mut output, flate2::Compression::new(3));
                let _ = io::copy(&mut Cursor::new(&data), &mut encoder).unwrap();
            })
        });
        group.finish()
    }
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(50);
    targets =
        roundtrip,
        read_encoder,
        write_encoder,
}
criterion_main!(benches);
