# Glossary

Definitions of the terms and tags used throughout this guide.

## Concepts

**UMI (Unique Molecular Identifier).** A short, often random, sequence attached to each original molecule before amplification. Reads that share a UMI and a mapping position are presumed to come from the same source molecule, which lets fgumi collapse PCR and sequencing duplicates into a consensus.

**Template.** The pair of reads (R1 and R2) sequenced from one molecule. fgumi groups and sorts by template, not by individual reads.

**Read structure.** A compact description of where UMI, template, sample-barcode, cell-barcode, and skip bases sit within each sequencing read, e.g. `8M+T` = an 8 bp UMI followed by all remaining template bases. See [Read Structures](read-structures.md).

**Family / family size.** A *family* is the set of reads grouped to one source molecule; *family size* is how many templates it contains (on paired-end data a read pair counts as one template), matching the `family_sizes.txt` histogram. The distribution of family sizes is a key QC signal — see [Working with Metrics](working-with-metrics.md).

**Simplex (single-strand) consensus.** A consensus built from reads of a single strand of the source molecule. Produced by `fgumi simplex`. See [Consensus Calling](consensus-calling.md).

**Duplex (double-strand) consensus.** A consensus that combines evidence from *both* strands of the source molecule, which suppresses errors that occur on only one strand. Produced by `fgumi duplex`. See [Duplex Consensus Calling](duplex-consensus-calling.md).

**CODEC.** A library protocol (Bae et al. 2023) in which a single read pair sequences both strands of the duplex molecule, so even one read pair can yield a duplex consensus. Produced by `fgumi codec`.

**Template-coordinate order.** A sort order that keeps all reads of a template together and orders templates by the outer 5′ coordinates of the pair. It is the required input order for `fgumi group`, `fgumi dedup`, and `fgumi downsample`, and is produced by `fgumi sort --order template-coordinate`. (Use `fgumi sort`, not `samtools sort`, ahead of these commands — see [Troubleshooting](troubleshooting.md).)

**Phred quality.** A base-quality score `Q = -10·log₁₀(P_error)`; higher is better (Q30 ≈ 1 error in 1000).

**Pre-UMI / post-UMI error rate.** The consensus model splits sequencing error into errors that occur *before* the UMI is integrated into the molecule (`--error-rate-pre-umi`) and *after* (`--error-rate-post-umi`). See [Consensus Calling](consensus-calling.md).

## Strand labels

**`/A` and `/B`.** Suffixes on the `MI` tag that `fgumi group --strategy paired` assigns to the two single-strand sub-families of one duplex molecule (e.g. `1/A` and `1/B`). The `/A` sub-family is the read pair whose read 1 5′ end comes at or before read 2's, ignoring soft-clipping and reference strand. See [Tracking Reads](tracking-reads.md).

**AB / BA strand.** The two strand orientations combined into a duplex consensus. Per-strand thresholds in `fgumi duplex` and `fgumi filter` take up to three values, `[duplex, AB, BA]`.

## BAM tags

| Tag | Set by | Meaning |
|-----|--------|---------|
| `RX` | `extract` | Raw UMI bases (input to `group`/`correct`). |
| `QX` | `extract` | UMI base qualities (optional). |
| `OX` | `correct` | Original UMI, before correction. |
| `MI` | `group` | Molecule ID assigning each read to a UMI family; carries `/A`·`/B` suffixes under the `paired` strategy. |
| `CB` | upstream | Cell barcode for single-cell data; reads are partitioned by `CB` before grouping, sorting, and deduplication. Not corrected by fgumi. |
| `tc` | `zipper` | Template-coordinate sort key added to secondary/supplementary reads so they sort and deduplicate correctly. |
| `pa` | `zipper` | Primary-alignment template sort-key coordinates on secondary/supplementary reads. |

Consensus callers add a family of per-read and per-base tags (`cD`, `cM`, `cE`, `cd`, `ce`, and the `a*`/`b*` single-strand variants for duplex). See [Tracking Reads](tracking-reads.md) for the full scheme.
