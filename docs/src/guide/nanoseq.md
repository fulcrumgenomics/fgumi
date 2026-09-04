# NanoSeq (Duplex-Seq) Guide

This guide describes how to process **NanoSeq** duplex sequencing data through fgumi's
consensus pipeline. NanoSeq (and the closely related BotSeqS/Duplex-Seq protocols) tag both
strands of each source molecule, so it fits fgumi's `paired` (duplex) grouping and `duplex`
consensus calling without any NanoSeq-specific code — the only protocol-specific choice is the
[read structure](read-structures.md).

## Background

NanoSeq is a whole-genome **duplex sequencing** protocol (Abascal et al., *Nature* 2021) designed
to call somatic mutations at single-molecule accuracy. Unlike hybridization-capture UMI panels,
where duplex recovery is rare, NanoSeq is engineered so that **most molecules are observed on both
strands** — which is exactly the regime `fgumi duplex` exists for.

Two design details matter for UMI processing:

- **A short inline barcode.** Each read begins with a 3&nbsp;bp molecular barcode, followed by a
  few constant bases from the protocol's restriction-based fragmentation that carry no molecular
  information and are skipped.
- **A read-barcode / mate-barcode pair.** The molecule's identity is the *pair* of the two reads'
  3&nbsp;bp barcodes together with the fragment's genomic coordinates. The two physical strands of
  one molecule are sequenced as separate templates whose barcode pairs are the **same two
  barcodes in swapped order**.

That swapped-pair structure is precisely what fgumi's `paired` grouping strategy canonicalizes: it
combines the two reads' UMIs into one order-independent molecule key and assigns the two strands
the `/A` and `/B` suffixes on the `MI` tag. No NanoSeq-aware logic is required — extracting the
barcodes into the standard `RX` tag and grouping with `--strategy paired` reproduces NanoSeq's
read-bundle model.

## Read Structure

NanoSeq's canonical extraction (the `cancerit/NanoSeq` `extract_tags.py -m 3 -s 4`) trims a
3&nbsp;bp barcode and skips 4&nbsp;bp of constant sequence on **each** read. In fgumi read-structure
terms that is:

```text
3M4S+T 3M4S+T
```

- `3M` — the 3&nbsp;bp molecular barcode (UMI), stored in `RX`.
- `4S` — 4&nbsp;bp of restriction-site/constant sequence, skipped.
- `+T` — the remaining bases are template (insert).

Both reads use the same structure because both carry a barcode. After extraction each read's `RX`
holds its own barcode plus its mate's, e.g. `RX:Z:GGG-GCA`; the opposite strand of the same
molecule carries the swapped pair `RX:Z:GCA-GGG`.

> If your NanoSeq library was prepared with a different barcode or skip length, adjust the `3M` /
> `4S` counts to match; check the `extract_tags` parameters used for your data.

## Pipeline

The NanoSeq pipeline is the [standard duplex pipeline](best-practices.md) with the NanoSeq read
structure at extraction. Alignment should be **genome-wide** (NanoSeq is WGS): the molecular key
depends on true fragment coordinates, so mapping to a subset reference and letting off-target reads
force-map would merge unrelated molecules.

```text
FASTQ → extract → fastq | aligner | zipper → sort → group → duplex → duplex-metrics
```

### Step 1: UMI Extraction

```bash
fgumi extract \
    --inputs R1.fastq.gz R2.fastq.gz \
    --output extracted.bam \
    --read-structures "3M4S+T" "3M4S+T" \
    --sample my_sample --library lib1
```

### Step 2: Alignment

Convert back to FASTQ, align genome-wide, and merge the UMI tags back onto the aligned reads with
`zipper`. `zipper` walks the extracted and aligned BAMs in lockstep, so the aligner must
**preserve read names exactly and emit records in the same queryname grouping and order as the
extracted BAM** — both inputs must carry the same set of read names in the same order, with
same-named records consecutive. Aligning the interleaved FASTQ with `bwa mem -p` (below) satisfies
this; use the same reference for `zipper --reference`.

```bash
samtools fastq extracted.bam \
    | bwa mem -t 16 -p -K 150000000 -Y hg38.fa - \
    | fgumi zipper \
        --unmapped extracted.bam \
        --reference hg38.fa \
        --output mapped.bam \
        --threads 16
```

### Step 3: Sort

Group requires template-coordinate order:

```bash
fgumi sort --input mapped.bam --output sorted.bam --order template-coordinate --threads 16
```

### Step 4: UMI Grouping (paired / duplex)

Use `--strategy paired` so the two strands of each molecule are grouped together and assigned
`/A` and `/B`. Pass `--edits 0` for exact barcode matching: NanoSeq's 3&nbsp;bp-per-read barcodes
give only 4096 possible pairs, so two distinct molecules that share a fragment's exact start/stop
(a routine coordinate duplicate) can fall within the default `--edits 1` of each other and be
merged. Canonical NanoSeq matches barcodes exactly — dropping a barcode-error read rather than
risking a merge — which preserves single-molecule specificity:

```bash
fgumi group --input sorted.bam --output grouped.bam --strategy paired --edits 0 --threads 16
```

### Step 5: Duplex Consensus

Call duplex consensus reads from the grouped BAM. See [Duplex Consensus
Calling](duplex-consensus-calling.md) for the `--min-reads` semantics (one, two, or three values,
high to low):

```bash
fgumi duplex --input grouped.bam --output consensus.bam --min-reads 1,1,0
```

`--min-reads 1,1,0` keeps every molecule (including single-strand ones); require `1,1,1` (or
higher on each strand) if you want **only** true duplex molecules.

### Step 6: QC Metrics

```bash
fgumi duplex-metrics --input grouped.bam --output my_sample
```

## Interpreting the Duplex Metrics

`duplex-metrics` writes `<prefix>.duplex_yield_metrics.txt`, whose last row (100% of the data)
summarizes duplex recovery. The families are counted at three levels of specificity:

- **CS** — coordinate + strand (physical fragments, pre-UMI).
- **SS** — single-strand tag families (coordinate + strand + UMI).
- **DS** — double-strand families (coordinate + UMI, collapsing the two strands).

The headline number is **`ds_fraction_duplexes`**: the fraction of double-strand families for which
**both** strands were actually observed — i.e. the true duplex rate. For a healthy NanoSeq library
this is high (a well-covered library commonly lands above 0.5); a value near zero means the data is
effectively simplex and duplex calling will yield little. `ds_fraction_duplexes_ideal` is the
asymptotic value at infinite sampling, so the gap between the two tells you how much duplex yield
more sequencing would recover.

The `<prefix>.duplex_family_sizes.txt` file breaks families down by `ab_size` / `ba_size` — the
number of reads observed on each strand — which is the distribution the `duplex --min-reads`
thresholds act on.

## Notes

- **Map genome-wide.** The molecule key is coordinate-anchored, so a subset reference distorts
  grouping. To restrict analysis to a region, align to the whole genome first and *then* subset by
  the aligned position, rather than aligning to a partial reference.
- **`paired` is required for duplex.** Only `--strategy paired` resolves the `/A` and `/B` strands;
  `adjacency`/`edit`/`identity` produce single-strand families and cannot feed `fgumi duplex`.
- **No UMI correction step; match barcodes exactly (`--edits 0`).** NanoSeq barcodes are random
  3-mers, not a fixed whitelist, so there is no `fgumi correct` step. Group with `--edits 0` rather
  than the default `--edits 1`, which tolerates one mismatch across the *entire* 6&nbsp;bp
  concatenated barcode pair (see [UMI Grouping](umi-grouping.md)) — enough that two unrelated
  molecules at a coordinate-duplicate position can be merged. Exact matching instead drops a
  barcode-error read from its family, the conservative failure mode NanoSeq's design assumes.
