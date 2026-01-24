#!/usr/bin/env python3
"""Generate synthetic paired-end FASTQ data for pipeline benchmarking.

This script creates realistic FASTQ pairs with:
- Configurable read counts to achieve target runtimes
- UMI family distributions (log-normal for realistic data)
- Proper read structure (UMI + skip + template)
- Gzipped output for realistic I/O

The data flows through: extract -> group -> simplex/duplex/codec -> filter
"""

import argparse
import gzip
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Tuple

import numpy as np


@dataclass
class PipelineConfig:
    """Configuration for synthetic pipeline data generation."""

    n_read_pairs: int = 10_000_000
    read_length: int = 150
    umi_length: int = 10
    skip_length: int = 1

    # Genomic distribution
    n_positions: int = 50_000

    # UMI family size distribution (log-normal)
    # mean ~20 reads per family with long tail
    family_size_mu: float = 2.5
    family_size_sigma: float = 1.0

    # Error rates
    umi_error_rate: float = 0.01  # 1% per base
    base_error_rate: float = 0.001  # 0.1% per base

    # Quality score parameters (Phred)
    mean_quality: int = 30
    quality_std: float = 5.0

    seed: int = 42

    # Duplex mode: generate both A and B strand reads
    # When True, each family has reads with UMI_A-UMI_B (A-strand) and UMI_B-UMI_A (B-strand)
    duplex: bool = True

    # Fraction of reads that are B-strand (only used when duplex=True)
    b_strand_fraction: float = 0.5


def generate_random_sequence(length: int, rng: random.Random) -> str:
    """Generate a random DNA sequence."""
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


def generate_quality_scores(length: int, mean_q: int, std_q: float, rng: random.Random) -> str:
    """Generate realistic quality scores (Phred+33 encoding)."""
    scores = []
    for _ in range(length):
        q = int(rng.gauss(mean_q, std_q))
        q = max(2, min(40, q))  # Clamp to valid range
        scores.append(chr(q + 33))
    return "".join(scores)


def generate_umi_families(
    n_total_reads: int,
    n_positions: int,
    family_mu: float,
    family_sigma: float,
    rng: random.Random,
    np_rng: np.random.Generator,
) -> Iterator[Tuple[int, int]]:
    """Generate UMI family assignments with log-normal distribution.

    Yields: (position_id, family_id_at_position)
    """
    # Distribute reads across positions
    reads_per_position = n_total_reads // n_positions
    remainder = n_total_reads % n_positions

    for pos_id in range(n_positions):
        n_reads_this_pos = reads_per_position + (1 if pos_id < remainder else 0)
        if n_reads_this_pos == 0:
            continue

        # Generate family sizes for this position using log-normal
        # Keep generating until we have enough reads
        family_id = 0
        reads_assigned = 0

        while reads_assigned < n_reads_this_pos:
            # Sample family size from log-normal
            family_size = int(np_rng.lognormal(family_mu, family_sigma))
            family_size = max(1, min(family_size, n_reads_this_pos - reads_assigned))

            for _ in range(family_size):
                yield (pos_id, family_id)
                reads_assigned += 1
                if reads_assigned >= n_reads_this_pos:
                    break

            family_id += 1


def generate_fastq_records(
    config: PipelineConfig,
    rng: random.Random,
    np_rng: np.random.Generator,
) -> Iterator[Tuple[str, str, str, str, str, str, str, str, str]]:
    """Generate paired FASTQ records.

    Yields: (name, seq1, qual1, seq2, qual2, umi_display, pos_id, family_id, strand)

    For duplex mode:
    - Each family has two UMI halves: UMI_A and UMI_B
    - A-strand reads: R1 has UMI_A, R2 has UMI_B -> extract produces RX:Z:UMI_A-UMI_B
    - B-strand reads: R1 has UMI_B, R2 has UMI_A -> extract produces RX:Z:UMI_B-UMI_A
    - The paired grouper recognizes these as the same molecule on opposite strands
    """
    # Pre-generate UMI pairs for each position/family
    # In duplex mode, each family has two UMI halves (A and B)
    position_umis: dict[int, dict[int, Tuple[str, str]]] = {}

    # Template length (after UMI + skip)
    template_len = config.read_length - config.umi_length - config.skip_length

    read_id = 0
    for pos_id, family_id in generate_umi_families(
        config.n_read_pairs,
        config.n_positions,
        config.family_size_mu,
        config.family_size_sigma,
        rng,
        np_rng,
    ):
        # Get or create UMI pair for this family
        if pos_id not in position_umis:
            position_umis[pos_id] = {}
        if family_id not in position_umis[pos_id]:
            umi_a = generate_random_sequence(config.umi_length, rng)
            if config.duplex:
                # Generate distinct UMI_B for duplex
                umi_b = generate_random_sequence(config.umi_length, rng)
            else:
                # For simplex, use same UMI for both
                umi_b = umi_a
            position_umis[pos_id][family_id] = (umi_a, umi_b)

        true_umi_a, true_umi_b = position_umis[pos_id][family_id]

        # Determine strand (A or B) for this read
        if config.duplex and rng.random() < config.b_strand_fraction:
            strand = "B"
            # B-strand: R1 gets UMI_B, R2 gets UMI_A
            umi_r1 = true_umi_b
            umi_r2 = true_umi_a
        else:
            strand = "A"
            # A-strand: R1 gets UMI_A, R2 gets UMI_B
            umi_r1 = true_umi_a
            umi_r2 = true_umi_b

        # Introduce UMI sequencing errors
        observed_umi_r1 = introduce_errors(umi_r1, config.umi_error_rate, rng)
        observed_umi_r2 = introduce_errors(umi_r2, config.umi_error_rate, rng)

        # Generate skip base (random)
        skip_base1 = rng.choice("ACGT")
        skip_base2 = rng.choice("ACGT")

        # Generate template sequences (R1 and R2)
        template1 = generate_random_sequence(template_len, rng)
        template2 = generate_random_sequence(template_len, rng)

        # Add errors to templates
        template1 = introduce_errors(template1, config.base_error_rate, rng)
        template2 = introduce_errors(template2, config.base_error_rate, rng)

        # Construct full read sequences: UMI + skip + template
        seq1 = observed_umi_r1 + skip_base1 + template1
        seq2 = observed_umi_r2 + skip_base2 + template2

        # Generate quality scores
        qual1 = generate_quality_scores(config.read_length, config.mean_quality, config.quality_std, rng)
        qual2 = generate_quality_scores(config.read_length, config.mean_quality, config.quality_std, rng)

        # Read name
        name = f"synth.{read_id}"

        # UMI display for truth file: show the canonical form (A-B)
        umi_display = f"{true_umi_a}-{true_umi_b}"

        yield (name, seq1, qual1, seq2, qual2, umi_display, str(pos_id), str(family_id), strand)
        read_id += 1

        # Progress reporting
        if read_id % 1_000_000 == 0:
            print(f"  Generated {read_id:,} read pairs...", file=sys.stderr)


def write_fastq_pair(
    config: PipelineConfig,
    output_r1: Path,
    output_r2: Path,
    output_truth: Path | None = None,
) -> None:
    """Write paired FASTQ files (gzipped).

    Args:
        config: Generation configuration
        output_r1: Path to R1 FASTQ output (gzipped)
        output_r2: Path to R2 FASTQ output (gzipped)
        output_truth: Optional path to truth file for validation
    """
    rng = random.Random(config.seed)
    np_rng = np.random.default_rng(config.seed)

    print(f"Generating synthetic FASTQ data:", file=sys.stderr)
    print(f"  Read pairs: {config.n_read_pairs:,}", file=sys.stderr)
    print(f"  Positions: {config.n_positions:,}", file=sys.stderr)
    print(f"  Read length: {config.read_length}", file=sys.stderr)
    print(f"  UMI length: {config.umi_length}", file=sys.stderr)
    print(f"  Read structure: {config.umi_length}M{config.skip_length}S+T", file=sys.stderr)
    print(f"  Duplex mode: {config.duplex}", file=sys.stderr)
    if config.duplex:
        print(f"  B-strand fraction: {config.b_strand_fraction}", file=sys.stderr)

    truth_file = None
    if output_truth:
        truth_file = open(output_truth, "w")
        truth_file.write("read_name\ttrue_umi\tposition\tfamily\tstrand\n")

    a_strand_count = 0
    b_strand_count = 0

    with gzip.open(output_r1, "wt") as f1, gzip.open(output_r2, "wt") as f2:
        for name, seq1, qual1, seq2, qual2, umi, pos, fam, strand in generate_fastq_records(config, rng, np_rng):
            # Write R1
            f1.write(f"@{name} 1/1\n")
            f1.write(f"{seq1}\n")
            f1.write("+\n")
            f1.write(f"{qual1}\n")

            # Write R2
            f2.write(f"@{name} 2/2\n")
            f2.write(f"{seq2}\n")
            f2.write("+\n")
            f2.write(f"{qual2}\n")

            if strand == "A":
                a_strand_count += 1
            else:
                b_strand_count += 1

            if truth_file:
                truth_file.write(f"{name}\t{umi}\t{pos}\t{fam}\t{strand}\n")

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
        description="Generate synthetic paired-end FASTQ data for pipeline benchmarking"
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
        "-p", "--n-positions", type=int, default=50_000,
        help="Number of genomic positions (default: 50000)"
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
        "--umi-error-rate", type=float, default=0.01,
        help="UMI sequencing error rate (default: 0.01)"
    )
    parser.add_argument(
        "--base-error-rate", type=float, default=0.001,
        help="Template base error rate (default: 0.001)"
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

    config = PipelineConfig(
        n_read_pairs=args.n_read_pairs,
        read_length=args.read_length,
        umi_length=args.umi_length,
        skip_length=args.skip_length,
        n_positions=args.n_positions,
        family_size_mu=args.family_mu,
        family_size_sigma=args.family_sigma,
        umi_error_rate=args.umi_error_rate,
        base_error_rate=args.base_error_rate,
        seed=args.seed,
        duplex=duplex,
        b_strand_fraction=args.b_strand_fraction,
    )

    output_r1 = Path(f"{args.output_prefix}_1.fastq.gz")
    output_r2 = Path(f"{args.output_prefix}_2.fastq.gz")

    write_fastq_pair(config, output_r1, output_r2, args.truth)


if __name__ == "__main__":
    main()
