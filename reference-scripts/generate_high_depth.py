#!/usr/bin/env python3
"""Generate synthetic high-depth BAM data for MIH benchmarking.

This script creates synthetic BAM files with controlled numbers of unique UMIs
per genomic position, enabling benchmarking of the MIH optimization which only
activates when there are many UMIs at a single position.

The generated data simulates amplicon sequencing or other high-depth scenarios
where a single position may have hundreds or thousands of unique UMIs.
"""

import argparse
import random
import sys
from pathlib import Path

import pysam


def generate_random_umi(length: int = 12) -> str:
    """Generate a random UMI sequence."""
    return "".join(random.choice("ACGT") for _ in range(length))


def generate_umis_with_errors(
    n_unique: int, n_total: int, umi_length: int = 12, error_rate: float = 0.1
) -> list[str]:
    """Generate UMIs with some sequencing errors.

    Creates n_unique "true" UMIs, then generates n_total UMIs by either:
    - Using the true UMI (90% of time)
    - Introducing 1 error (10% of time)

    This simulates realistic UMI data where some reads have sequencing errors.
    """
    # Generate unique true UMIs
    true_umis = [generate_random_umi(umi_length) for _ in range(n_unique)]

    # Generate observed UMIs with occasional errors
    result = []
    for _ in range(n_total):
        true_umi = random.choice(true_umis)
        if random.random() < error_rate:
            # Introduce 1 error
            pos = random.randint(0, len(true_umi) - 1)
            bases = [b for b in "ACGT" if b != true_umi[pos]]
            mutated = true_umi[:pos] + random.choice(bases) + true_umi[pos + 1 :]
            result.append(mutated)
        else:
            result.append(true_umi)

    return result


def create_synthetic_bam(
    output_path: Path,
    n_positions: int,
    unique_umis_per_position: int,
    reads_per_position: int,
    umi_length: int = 12,
    read_length: int = 150,
    reference_name: str = "chr1",
    reference_length: int = 1_000_000,
) -> None:
    """Create a synthetic BAM file with controlled UMI depth.

    Args:
        output_path: Path to write the BAM file
        n_positions: Number of genomic positions to generate
        unique_umis_per_position: Number of unique UMIs per position
        reads_per_position: Total reads per position
        umi_length: Length of UMI sequences
        read_length: Length of read sequences
        reference_name: Reference sequence name
        reference_length: Reference sequence length
    """
    # Create header
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"LN": reference_length, "SN": reference_name}],
            "RG": [{"ID": "A", "SM": "synthetic", "LB": "lib1", "PL": "ILLUMINA"}],
        }
    )

    # Generate positions spread across the reference
    position_spacing = reference_length // (n_positions + 1)
    positions = [position_spacing * (i + 1) for i in range(n_positions)]

    # Create reads
    reads = []
    read_id = 0

    for pos in positions:
        # Generate UMIs for this position
        umis = generate_umis_with_errors(
            n_unique=unique_umis_per_position,
            n_total=reads_per_position,
            umi_length=umi_length,
        )

        for umi in umis:
            cigar_str = f"{read_length}M"

            # Create R1 (forward read)
            r1 = pysam.AlignedSegment(header)
            r1.query_name = f"read_{read_id}"
            r1.query_sequence = "A" * read_length
            r1.query_qualities = pysam.qualitystring_to_array("I" * read_length)
            r1.reference_id = 0
            r1.reference_start = pos
            r1.mapping_quality = 60
            r1.cigar = [(0, read_length)]  # 150M
            r1.flag = 99  # paired, proper pair, mate reverse, first in pair
            r1.next_reference_id = 0
            r1.next_reference_start = pos + read_length + 100
            r1.template_length = read_length * 2 + 100
            r1.set_tag("RG", "A")
            r1.set_tag("RX", umi)  # UMI tag
            r1.set_tag("MC", cigar_str)  # Mate CIGAR for template-coordinate sort

            # Create R2 (reverse read)
            r2 = pysam.AlignedSegment(header)
            r2.query_name = f"read_{read_id}"
            r2.query_sequence = "T" * read_length
            r2.query_qualities = pysam.qualitystring_to_array("I" * read_length)
            r2.reference_id = 0
            r2.reference_start = pos + read_length + 100
            r2.mapping_quality = 60
            r2.cigar = [(0, read_length)]  # 150M
            r2.flag = 147  # paired, proper pair, reverse, second in pair
            r2.next_reference_id = 0
            r2.next_reference_start = pos
            r2.template_length = -(read_length * 2 + 100)
            r2.set_tag("RG", "A")
            r2.set_tag("RX", umi)  # Same UMI for mate
            r2.set_tag("MC", cigar_str)  # Mate CIGAR for template-coordinate sort

            reads.append(r1)
            reads.append(r2)
            read_id += 1

    # Sort by coordinate
    reads.sort(key=lambda r: (r.reference_start, r.query_name, r.is_read2))

    # Write BAM
    with pysam.AlignmentFile(str(output_path), "wb", header=header) as bam:
        for read in reads:
            bam.write(read)

    # Index the BAM
    pysam.index(str(output_path))


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic high-depth BAM data for MIH benchmarking"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Output BAM file path"
    )
    parser.add_argument(
        "-p",
        "--positions",
        type=int,
        default=100,
        help="Number of genomic positions (default: 100)",
    )
    parser.add_argument(
        "-u",
        "--unique-umis",
        type=int,
        default=50,
        help="Unique UMIs per position (default: 50)",
    )
    parser.add_argument(
        "-r",
        "--reads-per-position",
        type=int,
        default=500,
        help="Total reads per position (default: 500)",
    )
    parser.add_argument(
        "-l", "--umi-length", type=int, default=12, help="UMI length (default: 12)"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42, help="Random seed (default: 42)"
    )

    args = parser.parse_args()

    random.seed(args.seed)

    print(f"Generating synthetic BAM with:", file=sys.stderr)
    print(f"  Positions: {args.positions}", file=sys.stderr)
    print(f"  Unique UMIs per position: {args.unique_umis}", file=sys.stderr)
    print(f"  Reads per position: {args.reads_per_position}", file=sys.stderr)
    print(f"  Total reads: {args.positions * args.reads_per_position * 2}", file=sys.stderr)

    create_synthetic_bam(
        output_path=args.output,
        n_positions=args.positions,
        unique_umis_per_position=args.unique_umis,
        reads_per_position=args.reads_per_position,
        umi_length=args.umi_length,
    )

    print(f"Wrote: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
