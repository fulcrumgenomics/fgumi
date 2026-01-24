# fgumi simulate CLI Reference

Generate synthetic sequencing data for testing and benchmarking the fgumi pipeline.

**Requires:** `cargo build --features simulate`

## Commands Overview

| Command | Output | Pipeline Stage | Reference FASTA |
|---------|--------|----------------|-----------------|
| `fgumi simulate fastq-reads` | R1/R2 FASTQ.gz | Input for `extract` | Not used |
| `fgumi simulate mapped-reads` | Template-coord sorted BAM | Input for `group` | Optional |
| `fgumi simulate grouped-reads` | Template-coord sorted BAM with MI tags | Input for `simplex`/`duplex`/`codec` | Optional |
| `fgumi simulate consensus-reads` | Unmapped BAM with consensus tags | Input for `filter` | Not used |
| `fgumi simulate correct-reads` | Unmapped BAM + includelist | Input for `correct` | Not used |

**Note:** Commands that support `--reference` will sample positions from real chromosomes and extract actual genomic sequences. Without a reference, synthetic random sequences are generated.

---

## fgumi simulate fastq-reads

Generate paired-end FASTQ files with UMI sequences for input to `fgumi extract`.

### Usage

```bash
fgumi simulate fastq-reads \
    --r1 output_R1.fastq.gz \
    --r2 output_R2.fastq.gz \
    [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `-1, --r1` | PATH | Output R1 FASTQ.gz file |
| `-2, --r2` | PATH | Output R2 FASTQ.gz file |
| `--truth` | PATH | Output truth TSV file (for validation) |

### Simulation Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --num-molecules` | INT | 1000 | Number of unique molecules to simulate |
| `-l, --read-length` | INT | 150 | Length of each read in bases |
| `-u, --umi-length` | INT | 8 | Length of UMI sequence in bases |
| `--read-structure-r1` | STRING | `8M+T` | Read structure for R1 (fgbio notation) |
| `--read-structure-r2` | STRING | `+T` | Read structure for R2 (fgbio notation) |
| `--seed` | INT | (random) | Random seed for reproducibility |

### Quality Model Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--warmup-bases` | INT | 10 | Number of bases before peak quality is reached |
| `--warmup-quality` | INT | 25 | Starting quality score during warmup phase |
| `--peak-quality` | INT | 37 | Peak quality score (Phred) |
| `--decay-start` | INT | 100 | Position where quality decay begins |
| `--decay-rate` | FLOAT | 0.08 | Quality drop per base after decay starts |
| `--quality-noise` | FLOAT | 2.0 | Standard deviation of quality noise |
| `--r2-quality-offset` | INT | -2 | Quality offset for R2 reads (typically negative) |

### Family Size Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--family-size-dist` | STRING | `lognormal` | Distribution: `lognormal`, `negbin`, or path to histogram |
| `--family-size-mean` | FLOAT | 3.0 | Mean family size (for lognormal) |
| `--family-size-stddev` | FLOAT | 2.0 | Family size standard deviation (for lognormal) |
| `--family-size-r` | FLOAT | 2.0 | r parameter (for negative binomial) |
| `--family-size-p` | FLOAT | 0.5 | p parameter (for negative binomial) |
| `--min-family-size` | INT | 1 | Minimum reads per family |

### Insert Size Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--insert-size-mean` | FLOAT | 300.0 | Mean insert size |
| `--insert-size-stddev` | FLOAT | 50.0 | Insert size standard deviation |
| `--insert-size-min` | INT | 50 | Minimum insert size |
| `--insert-size-max` | INT | 800 | Maximum insert size |

### Truth File Format

The truth TSV file contains ground truth for validation:

| Column | Description |
|--------|-------------|
| `read_name` | Read name (matches FASTQ header) |
| `true_umi` | The true UMI sequence (before any errors) |
| `molecule_id` | Unique molecule identifier |
| `family_id` | Family within the molecule |
| `strand` | Strand (A or B for duplex) |

### Example

```bash
# Generate 10,000 molecules with 8bp UMIs
fgumi simulate fastq-reads \
    --r1 sim_R1.fastq.gz \
    --r2 sim_R2.fastq.gz \
    --truth sim_truth.tsv \
    --num-molecules 10000 \
    --umi-length 8 \
    --read-structure-r1 "8M142T" \
    --read-structure-r2 "150T" \
    --seed 42
```

---

## fgumi simulate mapped-reads

Generate template-coordinate sorted BAM with paired alignments for input to `fgumi group`.

### Usage

```bash
fgumi simulate mapped-reads \
    --output output.bam \
    [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `-o, --output` | PATH | Output BAM file (template-coordinate sorted) |
| `--truth` | PATH | Output truth TSV file (for validation) |

### Simulation Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --num-molecules` | INT | 1000 | Number of unique molecules to simulate |
| `-l, --read-length` | INT | 150 | Length of each read in bases |
| `-u, --umi-length` | INT | 8 | Length of UMI sequence in bases |
| `--seed` | INT | (random) | Random seed for reproducibility |
| `-t, --threads` | INT | 1 | Number of writer threads |

### Reference Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-r, --reference` | PATH | (none) | Reference FASTA file (sequences sampled from here) |
| `--ref-name` | STRING | `chr1` | Synthetic reference name (only used if no --reference) |
| `--ref-length` | INT | 10000000 | Synthetic reference length (only used if no --reference) |

When `--reference` is provided:
- Positions are sampled from real chromosomes (weighted by length)
- Read sequences are extracted from the reference
- BAM header contains actual contig names and lengths
- Reads will map correctly if used with an aligner

When `--reference` is not provided:
- A synthetic single-contig reference is used
- Random sequences are generated
- Useful for quick testing without a reference file

### Alignment Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--mapq` | INT | 60 | Mapping quality for aligned reads |
| `--unmapped-fraction` | FLOAT | 0.0 | Fraction of reads to leave unmapped |

### Position Distribution Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--num-positions` | INT | (derived) | Number of genomic positions to use (default: num-molecules) |
| `--umis-per-position` | INT | 1 | Number of unique UMIs per position |

By default, each molecule gets a unique position. For high-depth benchmarking (testing MIH optimization in `group`), use fewer positions with many UMIs per position.

### Quality Model Options

(Same as fastq-reads)

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--warmup-bases` | INT | 10 | Number of bases before peak quality |
| `--warmup-quality` | INT | 25 | Starting quality during warmup |
| `--peak-quality` | INT | 37 | Peak quality score |
| `--decay-start` | INT | 100 | Position where decay begins |
| `--decay-rate` | FLOAT | 0.08 | Quality drop per base |
| `--quality-noise` | FLOAT | 2.0 | Quality noise std dev |
| `--r2-quality-offset` | INT | -2 | R2 quality offset |

### Family Size Options

(Same as fastq-reads)

### Insert Size Options

(Same as fastq-reads)

### Output Tags

| Tag | Type | Description |
|-----|------|-------------|
| `RX` | String | Raw UMI sequence |
| `RG` | String | Read group (default: "A") |

### Truth File Format

| Column | Description |
|--------|-------------|
| `read_name` | Read name (matches BAM QNAME) |
| `true_umi` | The true UMI sequence (before any errors) |
| `molecule_id` | Unique molecule identifier |
| `chrom` | Chromosome/contig name |
| `position` | 1-based genomic position |
| `strand` | Strand (A or B for duplex) |

### Example

```bash
# Generate mapped reads from a reference FASTA
fgumi simulate mapped-reads \
    --output sim_mapped.bam \
    --truth sim_truth.tsv \
    --reference hg38.fa \
    --num-molecules 5000 \
    --seed 42 \
    --threads 4

# Generate mapped reads with synthetic reference (no FASTA needed)
fgumi simulate mapped-reads \
    --output sim_mapped.bam \
    --truth sim_truth.tsv \
    --ref-name chr1 \
    --ref-length 50000000 \
    --num-molecules 5000 \
    --seed 42

# High-depth mode: many UMIs at few positions (for MIH/group benchmarking)
fgumi simulate mapped-reads \
    --output sim_high_depth.bam \
    --truth sim_truth.tsv \
    --num-molecules 50000 \
    --num-positions 100 \
    --umis-per-position 500 \
    --seed 42
```

---

## fgumi simulate grouped-reads

Generate template-coordinate sorted BAM with MI (molecule ID) tags for input to consensus callers (`simplex`, `duplex`, `codec`).

### Usage

```bash
fgumi simulate grouped-reads \
    --output output.bam \
    [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `-o, --output` | PATH | Output BAM file (template-coordinate sorted) |
| `--truth` | PATH | Output truth TSV file (for validation) |

### Simulation Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --num-molecules` | INT | 1000 | Number of unique molecules to simulate |
| `-l, --read-length` | INT | 150 | Length of each read in bases |
| `-u, --umi-length` | INT | 8 | Length of UMI sequence in bases |
| `--seed` | INT | (random) | Random seed for reproducibility |
| `-t, --threads` | INT | 1 | Number of writer threads |

### Duplex Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--duplex` | FLAG | false | Generate duplex-style MI tags (e.g., "1/A", "1/B") |
| `--strand-alpha` | FLOAT | 5.0 | Beta distribution alpha for A/B strand ratio |
| `--strand-beta` | FLOAT | 5.0 | Beta distribution beta for A/B strand ratio |

### Reference Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-r, --reference` | PATH | (none) | Reference FASTA file (sequences sampled from here) |
| `--ref-name` | STRING | `chr1` | Synthetic reference name (only used if no --reference) |
| `--ref-length` | INT | 10000000 | Synthetic reference length (only used if no --reference) |

See `mapped-reads` for details on reference vs. synthetic mode.

### Quality Model Options

(Same as fastq-reads)

### Family Size Options

(Same as fastq-reads)

### Insert Size Options

(Same as fastq-reads)

### Output Tags

| Tag | Type | Description |
|-----|------|-------------|
| `RX` | String | Raw UMI sequence |
| `MI` | String | Molecule ID (integer for simplex, "N/A" or "N/B" for duplex) |
| `RG` | String | Read group (default: "A") |

### Truth File Format

| Column | Description |
|--------|-------------|
| `read_name` | Read name (matches BAM QNAME) |
| `true_umi` | The true UMI sequence (before any errors) |
| `molecule_id` | Unique molecule identifier |
| `expected_mi` | Expected MI tag value after grouping |
| `chrom` | Chromosome/contig name |
| `position` | 1-based genomic position |
| `strand` | Strand (A or B for duplex) |

### Example

```bash
# Generate simplex grouped reads with reference
fgumi simulate grouped-reads \
    --output sim_grouped.bam \
    --truth sim_truth.tsv \
    --reference hg38.fa \
    --num-molecules 5000 \
    --seed 42

# Generate duplex grouped reads with strand bias (synthetic reference)
fgumi simulate grouped-reads \
    --output sim_duplex_grouped.bam \
    --truth sim_truth.tsv \
    --num-molecules 5000 \
    --duplex \
    --strand-alpha 5.0 \
    --strand-beta 5.0 \
    --seed 42
```

---

## fgumi simulate consensus-reads

Generate unmapped BAM with consensus tags (cD, cM, cE, etc.) for input to `fgumi filter`.

### Usage

```bash
fgumi simulate consensus-reads \
    --output output.bam \
    [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `-o, --output` | PATH | Output BAM file (unmapped) |

### Simulation Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --num-reads` | INT | 1000 | Number of consensus read pairs to generate |
| `-l, --read-length` | INT | 150 | Length of each read in bases |
| `--seed` | INT | (random) | Random seed for reproducibility |
| `-t, --threads` | INT | 1 | Number of writer threads |

### Consensus Tag Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--min-depth` | INT | 1 | Minimum consensus depth (cM tag) |
| `--max-depth` | INT | 10 | Maximum consensus depth (cD tag) |
| `--depth-mean` | FLOAT | 5.0 | Mean depth for sampling |
| `--depth-stddev` | FLOAT | 2.0 | Depth standard deviation |
| `--error-rate-mean` | FLOAT | 0.01 | Mean error rate (cE tag) |
| `--error-rate-stddev` | FLOAT | 0.005 | Error rate standard deviation |

### Duplex Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--duplex` | FLAG | false | Generate duplex consensus tags (aD, bD, aM, bM, aE, bE) |
| `--strand-alpha` | FLOAT | 5.0 | Beta distribution alpha for A/B depth ratio |
| `--strand-beta` | FLOAT | 5.0 | Beta distribution beta for A/B depth ratio |

### Quality Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--consensus-quality` | INT | 40 | Base quality for consensus reads |

### Output Tags (Simplex)

| Tag | Type | Description |
|-----|------|-------------|
| `cD` | Int | Maximum per-base depth |
| `cM` | Int | Minimum per-base depth |
| `cE` | Float | Consensus error rate |
| `cd` | IntArray | Per-base depth array |
| `ce` | IntArray | Per-base error count array |

### Output Tags (Duplex, in addition to above)

| Tag | Type | Description |
|-----|------|-------------|
| `aD` | Int | A-strand maximum depth |
| `aM` | Int | A-strand minimum depth |
| `aE` | Float | A-strand error rate |
| `bD` | Int | B-strand maximum depth |
| `bM` | Int | B-strand minimum depth |
| `bE` | Float | B-strand error rate |
| `ad` | IntArray | A-strand per-base depth |
| `bd` | IntArray | B-strand per-base depth |

### Example

```bash
# Generate simplex consensus reads
fgumi simulate consensus-reads \
    --output sim_consensus.bam \
    --num-reads 10000 \
    --min-depth 2 \
    --max-depth 20 \
    --seed 42

# Generate duplex consensus reads
fgumi simulate consensus-reads \
    --output sim_duplex_consensus.bam \
    --num-reads 10000 \
    --duplex \
    --seed 42
```

---

## fgumi simulate correct-reads

Generate unmapped BAM and UMI includelist for input to `fgumi correct`.

### Usage

```bash
fgumi simulate correct-reads \
    --output output.bam \
    --includelist umis.txt \
    [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `-o, --output` | PATH | Output BAM file (unmapped) |
| `-i, --includelist` | PATH | Output UMI includelist file (one UMI per line) |
| `--truth` | PATH | Output truth TSV file (for validation) |

### Simulation Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --num-reads` | INT | 10000 | Number of reads to generate |
| `--num-umis` | INT | 1000 | Number of unique UMIs in includelist |
| `-u, --umi-length` | INT | 8 | Length of UMI sequence in bases |
| `--read-length` | INT | 100 | Length of template sequence |
| `--seed` | INT | (random) | Random seed for reproducibility |
| `-t, --threads` | INT | 1 | Number of writer threads |

### Error Distribution Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--exact-fraction` | FLOAT | 0.4 | Fraction with exact UMI match (0 edits) |
| `--edit1-fraction` | FLOAT | 0.3 | Fraction with 1 edit distance |
| `--edit2-fraction` | FLOAT | 0.2 | Fraction with 2 edit distance |
| `--multi-fraction` | FLOAT | 0.1 | Fraction with 3+ edits (should not correct) |

### Quality Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--quality` | INT | 30 | Base quality for all bases |

### Output Tags

| Tag | Type | Description |
|-----|------|-------------|
| `RX` | String | Observed (possibly erroneous) UMI sequence |

### Includelist Format

The includelist is a plain text file with one UMI per line, sorted alphabetically:

```
AAAACCCC
AAAACCCT
AAAACCGG
...
```

### Truth File Format

| Column | Description |
|--------|-------------|
| `read_name` | Read name (matches BAM QNAME) |
| `true_umi` | The correct UMI from the includelist |
| `observed_umi` | The (possibly erroneous) UMI in the RX tag |
| `expected_correction` | Expected UMI after correction |
| `edit_distance` | Edit distance between observed and true |
| `error_type` | Type: `exact`, `edit1`, `edit2`, or `multi` |

### Example

```bash
# Generate correction test data
fgumi simulate correct-reads \
    --output sim_correct.bam \
    --includelist sim_umis.txt \
    --truth sim_truth.tsv \
    --num-reads 100000 \
    --num-umis 5000 \
    --umi-length 12 \
    --exact-fraction 0.5 \
    --edit1-fraction 0.3 \
    --edit2-fraction 0.15 \
    --multi-fraction 0.05 \
    --seed 42
```

---

## Shared Concepts

### Read Structure Notation

Uses fgbio-style read structure notation:
- `M` = Molecular barcode (UMI)
- `T` = Template
- `S` = Skip (ignored bases)
- `B` = Sample barcode
- `+` = Variable length (consumes remaining bases)

Examples:
- `8M142T` = 8bp UMI followed by 142bp template
- `8M+T` = 8bp UMI followed by variable-length template
- `+T` = All template (no UMI)

### Family Size Distributions

Three distribution types are supported:

1. **Log-normal** (default): Natural for PCR amplification
   - Parameters: `--family-size-mean`, `--family-size-stddev`

2. **Negative binomial**: Alternative PCR model
   - Parameters: `--family-size-r`, `--family-size-p`

3. **Empirical**: Load from `fgumi group -f` histogram output
   - Parameter: `--family-size-dist /path/to/histogram.tsv`

### Quality Score Model

Quality scores follow a three-phase model:

1. **Warmup** (positions 0 to `warmup-bases`): Quality ramps from `warmup-quality` to `peak-quality`
2. **Peak** (positions `warmup-bases` to `decay-start`): Quality stays at `peak-quality`
3. **Decay** (positions after `decay-start`): Quality decreases by `decay-rate` per base

R2 reads have an additional offset (`r2-quality-offset`, typically -2) applied.

### Template-Coordinate Sorting

For `mapped-reads` and `grouped-reads`, output is sorted by template coordinate:
- Primary sort: 5' position of the leftmost read in the pair
- Secondary sort: Read name (for determinism)

BAM header includes: `SO:unsorted`, `GO:query`, `SS:template-coordinate`

### High-Depth Benchmarking

To benchmark the MIH (Multiple Identical Hits) optimization in `fgumi group`, use the position distribution options to create data with many UMIs at few positions:

```bash
# 100 positions, each with 500 unique UMIs, ~10 reads per UMI = 500,000 reads
fgumi simulate mapped-reads \
    --output high_depth.bam \
    --truth high_depth_truth.tsv \
    --num-molecules 500000 \
    --num-positions 100 \
    --umis-per-position 500 \
    --family-size-mean 10 \
    --seed 42
```

This creates the kind of position-clustered data that stresses the UMI assignment algorithm and tests the MIH optimization path in `group`.
