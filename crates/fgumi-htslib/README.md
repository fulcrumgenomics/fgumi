# fgumi-htslib

CRAM reading and writing for the [fgumi](https://github.com/fulcrumgenomics/fgumi) family of crates, through a vendored [htslib](https://github.com/samtools/htslib) whose CPU-heavy work runs on the caller's threads.

htslib parallelizes CRAM on its own thread pool. This crate exposes the same three-stage shape, so the stages can run on fgumi's pipeline workers:

- **Read:** `CramReader::next_slice` (serial) yields `SliceJob`s, and `SliceJob::decode` (any thread) returns a slice's records as framed BAM.
- **Write:** `CramWriter::write_raw` (serial) packs records into containers. `EncodeJob::encode` (any thread) compresses one container, and `EncodedContainer::write` (serial, in order) writes it. `CramWriter::seal` and `CramFinisher::finish` close the file.

The vendored htslib (1.23, with htscodecs 1.6.5) carries one change, `vendor/container-executor.patch`. It adds `cram_set_container_executor`, which hands full containers to the caller instead of to htslib's thread pool. The same change has been proposed upstream as [samtools/htslib#2157](https://github.com/samtools/htslib/pull/2157); once it is released, the patch can be dropped. The read side needs no patch: `csrc/fgumi_hts.c` drives htslib's slice decoding directly.

This is fgumi's only crate with `unsafe` FFI. Its invariants are documented in the repository's `CLAUDE.md` ("Approved FFI crate"). Building it compiles htslib's C sources, so a C toolchain is required. fgumi enables it with the opt-in `cram` feature.

CRAM input older than version 2.1 is refused: those versions have no EOF container, so a truncated file cannot be told from a complete one. Re-encode such a file with `samtools view -C` (or convert it to BAM) first.

CRAM output is not byte-for-byte reproducible from run to run, even with identical input and settings, because htslib chooses each block's codec adaptively from trial compressions. The decoded records are the same every time, so compare CRAM outputs by their records, not by checksum.
