#!/usr/bin/env python3
"""Generate synthetic paired-end reads from a reference FASTA.

Unlike generate_pipeline_data.py which uses random sequences, this script
samples reads directly from a reference genome so they will map correctly.
This enables proper benchmarking of the full pipeline including mapping.

Features:
- Samples read positions from reference genome
- Generates UMI families with realistic size distributions
- Adds sequencing errors at configurable rates
- Outputs gzipped FASTQ pairs ready for the pipeline
"""

import argparse
import gzip
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Tuple

import numpy as np
import pysam


@dataclass
class ReferenceReadConfig:
    """Configuration for reference-based read generation."""

    n_read_pairs: int = 10_000_000
    read_length: int = 150
    umi_length: int = 10
    skip_length: int = 1

    # Fragment size distribution (for paired-end)
    fragment_mean: int = 300
    fragment_std: int = 50

    # UMI family size distribution (log-normal)
    family_size_mu: float = 2.5
    family_size_sigma: float = 1.0

    # Error rates
    umi_error_rate: float = 0.01
    base_error_rate: float = 0.001

    # Quality scores
    mean_quality: int = 30
    quality_std: float = 5.0

    seed: int = 42

    # Duplex mode: generate both A and B strand reads
    duplex: bool = True
    b_strand_fraction: float = 0.5


def load_reference(fasta_path: Path) -> dict[str, str]:
    """Load reference sequences from FASTA file."""
    ref = {}
    with pysam.FastaFile(str(fasta_path)) as fa:
        for name in fa.references:
            ref[name] = fa.fetch(name)
    return ref


def generate_random_umi(length: int, rng: random.Random) -> str:
    """Generate a random UMI sequence."""
    return "".join(rng.choices("ACGT", k=length))


def introduce_errors(seq: str, error_rate: float, rng: random.Random) -> str:
    """Introduce random substitution errors into a sequence."""
    if error_rate <= 0:
        return seq

    seq_list = list(seq)
    for i in range(len(seq_list)):
        if rng.random() < error_rate:
            current = seq_list[i]
            alternatives = [b for b in "ACGT" if b != current]
            seq_list[i] = rng.choice(alternatives)
    return "".join(seq_list)


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


def generate_quality_scores(length: int, mean_q: int, std_q: float, rng: random.Random) -> str:
    """Generate realistic quality scores (Phred+33 encoding)."""
    scores = []
    for _ in range(length):
        q = int(rng.gauss(mean_q, std_q))
        q = max(2, min(40, q))
        scores.append(chr(q + 33))
    return "".join(scores)


def sample_positions(
    ref_seqs: dict[str, str],
    n_positions: int,
    fragment_mean: int,
    read_length: int,
    rng: random.Random,
) -> List[Tuple[str, int]]:
    """Sample valid positions from reference for read generation."""
    positions = []

    # Get valid chromosomes and their lengths
    valid_chroms = [(name, len(seq)) for name, seq in ref_seqs.items()
                    if len(seq) > fragment_mean + read_length * 2]

    if not valid_chroms:
        raise ValueError("No chromosomes long enough for read generation")

    # Weight by chromosome length
    total_len = sum(length for _, length in valid_chroms)
    weights = [length / total_len for _, length in valid_chroms]

    for _ in range(n_positions):
        # Select chromosome
        chrom_idx = rng.choices(range(len(valid_chroms)), weights=weights)[0]
        chrom_name, chrom_len = valid_chroms[chrom_idx]

        # Select position (leaving room for fragment)
        max_pos = chrom_len - fragment_mean - read_length
        if max_pos < 0:
            continue
        pos = rng.randint(0, max_pos)
        positions.append((chrom_name, pos))

    return positions


def generate_reads_from_reference(
    config: ReferenceReadConfig,
    ref_seqs: dict[str, str],
    rng: random.Random,
    np_rng: np.random.Generator,
) -> Iterator[Tuple[str, str, str, str, str, str, str, int, str, str]]:
    """Generate paired reads from reference positions.

    For duplex mode:
    - Each family has two UMI halves: UMI_A and UMI_B
    - A-strand reads: R1 has UMI_A, R2 has UMI_B -> extract produces RX:Z:UMI_A-UMI_B
    - B-strand reads: R1 has UMI_B, R2 has UMI_A -> extract produces RX:Z:UMI_B-UMI_A
    - The paired grouper recognizes these as the same molecule on opposite strands

    Yields: (name, seq1, qual1, seq2, qual2, umi_display, chrom, pos, family_id, strand)
    """
    # Calculate number of positions needed
    # Estimate based on average family size
    avg_family_size = int(np.exp(config.family_size_mu + config.family_size_sigma**2 / 2))
    n_positions = max(1000, config.n_read_pairs // avg_family_size)

    print(f"Sampling {n_positions} genomic positions...", file=sys.stderr)
    positions = sample_positions(
        ref_seqs, n_positions, config.fragment_mean, config.read_length, rng
    )

    # Template length (after UMI + skip)
    template_len = config.read_length - config.umi_length - config.skip_length

    read_id = 0
    # Store UMI pairs (A, B) for each position/family
    position_umis: dict[Tuple[str, int], dict[int, Tuple[str, str]]] = {}

    # Generate reads with family structure
    pos_idx = 0
    while read_id < config.n_read_pairs:
        chrom, pos = positions[pos_idx % len(positions)]
        pos_idx += 1

        # Get or create UMI families for this position
        pos_key = (chrom, pos)
        if pos_key not in position_umis:
            position_umis[pos_key] = {}

        # Generate family sizes
        family_id = len(position_umis[pos_key])
        family_size = max(1, int(np_rng.lognormal(config.family_size_mu, config.family_size_sigma)))
        family_size = min(family_size, config.n_read_pairs - read_id)

        # Create UMI pair for this family
        umi_a = generate_random_umi(config.umi_length, rng)
        if config.duplex:
            umi_b = generate_random_umi(config.umi_length, rng)
        else:
            umi_b = umi_a
        position_umis[pos_key][family_id] = (umi_a, umi_b)

        # Sample fragment size
        frag_size = max(
            config.read_length,
            int(rng.gauss(config.fragment_mean, config.fragment_std))
        )

        # Get template sequences from reference
        ref_seq = ref_seqs[chrom]
        if pos + frag_size > len(ref_seq):
            continue

        # R1 forward, R2 reverse complement
        template1 = ref_seq[pos:pos + template_len]
        r2_start = pos + frag_size - template_len
        template2 = reverse_complement(ref_seq[r2_start:r2_start + template_len])

        # Skip if templates contain Ns
        if "N" in template1.upper() or "N" in template2.upper():
            continue

        for _ in range(family_size):
            if read_id >= config.n_read_pairs:
                break

            # Determine strand (A or B) for this read
            if config.duplex and rng.random() < config.b_strand_fraction:
                strand = "B"
                # B-strand: R1 gets UMI_B, R2 gets UMI_A
                umi_r1, umi_r2 = umi_b, umi_a
            else:
                strand = "A"
                # A-strand: R1 gets UMI_A, R2 gets UMI_B
                umi_r1, umi_r2 = umi_a, umi_b

            # Introduce UMI sequencing errors
            observed_umi_r1 = introduce_errors(umi_r1, config.umi_error_rate, rng)
            observed_umi_r2 = introduce_errors(umi_r2, config.umi_error_rate, rng)

            # Skip base (random)
            skip_base1 = rng.choice("ACGT")
            skip_base2 = rng.choice("ACGT")

            # Add template errors
            obs_template1 = introduce_errors(template1, config.base_error_rate, rng)
            obs_template2 = introduce_errors(template2, config.base_error_rate, rng)

            # Construct full read sequences: UMI + skip + template
            seq1 = observed_umi_r1 + skip_base1 + obs_template1
            seq2 = observed_umi_r2 + skip_base2 + obs_template2

            # Generate quality scores
            qual1 = generate_quality_scores(
                config.read_length, config.mean_quality, config.quality_std, rng
            )
            qual2 = generate_quality_scores(
                config.read_length, config.mean_quality, config.quality_std, rng
            )

            name = f"synth.{read_id}"
            # UMI display for truth file: canonical form (A-B)
            umi_display = f"{umi_a}-{umi_b}"
            yield (name, seq1, qual1, seq2, qual2, umi_display, chrom, pos, str(family_id), strand)

            read_id += 1

            if read_id % 1_000_000 == 0:
                print(f"  Generated {read_id:,} read pairs...", file=sys.stderr)


def write_fastq_pair(
    config: ReferenceReadConfig,
    ref_path: Path,
    output_r1: Path,
    output_r2: Path,
    output_truth: Path | None = None,
) -> None:
    """Write paired FASTQ files from reference-based reads."""
    rng = random.Random(config.seed)
    np_rng = np.random.default_rng(config.seed)

    print(f"Loading reference: {ref_path}", file=sys.stderr)
    ref_seqs = load_reference(ref_path)
    print(f"  Loaded {len(ref_seqs)} sequences", file=sys.stderr)

    print(f"Generating {config.n_read_pairs:,} read pairs...", file=sys.stderr)
    print(f"  Read length: {config.read_length}", file=sys.stderr)
    print(f"  UMI structure: {config.umi_length}M{config.skip_length}S+T", file=sys.stderr)
    print(f"  Duplex mode: {config.duplex}", file=sys.stderr)
    if config.duplex:
        print(f"  B-strand fraction: {config.b_strand_fraction}", file=sys.stderr)

    truth_file = None
    if output_truth:
        truth_file = open(output_truth, "w")
        truth_file.write("read_name\ttrue_umi\tchrom\tpos\tfamily\tstrand\n")

    a_strand_count = 0
    b_strand_count = 0

    with gzip.open(output_r1, "wt") as f1, gzip.open(output_r2, "wt") as f2:
        for name, seq1, qual1, seq2, qual2, umi, chrom, pos, fam, strand in generate_reads_from_reference(
            config, ref_seqs, rng, np_rng
        ):
            # Write R1
            f1.write(f"@{name} 1:N:0:0\n")
            f1.write(f"{seq1}\n")
            f1.write("+\n")
            f1.write(f"{qual1}\n")

            # Write R2
            f2.write(f"@{name} 2:N:0:0\n")
            f2.write(f"{seq2}\n")
            f2.write("+\n")
            f2.write(f"{qual2}\n")

            if strand == "A":
                a_strand_count += 1
            else:
                b_strand_count += 1

            if truth_file:
                truth_file.write(f"{name}\t{umi}\t{chrom}\t{pos}\t{fam}\t{strand}\n")

    if truth_file:
        truth_file.close()

    print(f"Wrote: {output_r1}", file=sys.stderr)
    print(f"Wrote: {output_r2}", file=sys.stderr)
    if output_truth:
        print(f"Wrote: {output_truth}", file=sys.stderr)
    if config.duplex:
        print(f"  A-strand reads: {a_strand_count:,}", file=sys.stderr)
        print(f"  B-strand reads: {b_strand_count:,}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic paired-end reads from a reference genome"
    )
    parser.add_argument(
        "-r", "--reference", type=Path, required=True,
        help="Reference FASTA file"
    )
    parser.add_argument(
        "-o", "--output-prefix", type=str, required=True,
        help="Output prefix (will create PREFIX_1.fastq.gz and PREFIX_2.fastq.gz)"
    )
    parser.add_argument(
        "-n", "--n-read-pairs", type=int, default=10_000_000,
        help="Number of read pairs to generate (default: 10M)"
    )
    parser.add_argument(
        "-l", "--read-length", type=int, default=150,
        help="Read length (default: 150)"
    )
    parser.add_argument(
        "-u", "--umi-length", type=int, default=10,
        help="UMI length (default: 10)"
    )
    parser.add_argument(
        "-k", "--skip-length", type=int, default=1,
        help="Skip bases after UMI (default: 1)"
    )
    parser.add_argument(
        "--fragment-mean", type=int, default=300,
        help="Mean fragment size (default: 300)"
    )
    parser.add_argument(
        "--fragment-std", type=int, default=50,
        help="Fragment size standard deviation (default: 50)"
    )
    parser.add_argument(
        "--family-mu", type=float, default=2.5,
        help="Log-normal mu for family size distribution (default: 2.5)"
    )
    parser.add_argument(
        "--family-sigma", type=float, default=1.0,
        help="Log-normal sigma for family size distribution (default: 1.0)"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--truth", type=Path, default=None,
        help="Optional output path for truth file"
    )
    parser.add_argument(
        "--duplex", action="store_true", default=True,
        help="Generate duplex data with A and B strand reads (default: True)"
    )
    parser.add_argument(
        "--no-duplex", action="store_true",
        help="Generate simplex data (same UMI for R1 and R2)"
    )
    parser.add_argument(
        "--b-strand-fraction", type=float, default=0.5,
        help="Fraction of reads that are B-strand in duplex mode (default: 0.5)"
    )

    args = parser.parse_args()

    # Handle duplex/no-duplex flags
    duplex = args.duplex and not args.no_duplex

    config = ReferenceReadConfig(
        n_read_pairs=args.n_read_pairs,
        read_length=args.read_length,
        umi_length=args.umi_length,
        skip_length=args.skip_length,
        fragment_mean=args.fragment_mean,
        fragment_std=args.fragment_std,
        family_size_mu=args.family_mu,
        family_size_sigma=args.family_sigma,
        seed=args.seed,
        duplex=duplex,
        b_strand_fraction=args.b_strand_fraction,
    )

    output_r1 = Path(f"{args.output_prefix}_1.fastq.gz")
    output_r2 = Path(f"{args.output_prefix}_2.fastq.gz")

    write_fastq_pair(config, args.reference, output_r1, output_r2, args.truth)


if __name__ == "__main__":
    main()
