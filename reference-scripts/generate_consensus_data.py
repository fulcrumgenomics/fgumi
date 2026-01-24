#!/usr/bin/env python3
"""Generate synthetic consensus BAM data for filter benchmarking.

This script creates consensus BAM files with the proper tags expected by
FilterConsensusReads, allowing filter benchmarks to be decoupled from the
full pipeline (which would require enormous input data).

Consensus tags generated:
- Simplex: cD, cE, cM, cd, ce
- Duplex: aD, bD, cD, aE, bE, cE, aM, bM, cM, ad, bd, cd, ae, be, ce
"""

import argparse
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pysam


@dataclass
class ConsensusConfig:
    """Configuration for synthetic consensus data generation."""

    n_consensus_reads: int = 2_000_000
    read_length: int = 289  # Match real data
    umi_length: int = 10

    # Depth distribution (log-normal to match real data)
    depth_mu: float = 2.5  # Mean depth ~12
    depth_sigma: float = 0.8

    # Error rate distribution
    mean_error_rate: float = 0.001
    error_rate_std: float = 0.0005

    # Duplex vs simplex
    consensus_type: str = "simplex"  # "simplex" or "duplex"

    # Quality score parameters
    mean_quality: int = 35
    quality_std: float = 3.0

    seed: int = 42


def generate_random_sequence(length: int, rng: random.Random) -> str:
    """Generate a random DNA sequence."""
    return "".join(rng.choices("ACGT", k=length))


def generate_quality_scores(length: int, mean_q: int, std_q: float, rng: random.Random) -> str:
    """Generate realistic quality scores (Phred+33 encoding)."""
    scores = []
    for _ in range(length):
        q = int(rng.gauss(mean_q, std_q))
        q = max(2, min(40, q))
        scores.append(chr(q + 33))
    return "".join(scores)


def generate_depth_array(length: int, base_depth: int, rng: random.Random) -> List[int]:
    """Generate per-base depth array with some variation."""
    depths = []
    for _ in range(length):
        # Small variation around base depth
        d = base_depth + rng.randint(-1, 1)
        d = max(1, d)
        depths.append(d)
    return depths


def generate_error_array(length: int, mean_error: float, rng: random.Random) -> List[int]:
    """Generate per-base error count array."""
    errors = []
    for _ in range(length):
        # Most bases have 0 errors, occasionally 1-2
        if rng.random() < mean_error * 100:  # Scale up for visibility
            errors.append(rng.randint(1, 2))
        else:
            errors.append(0)
    return errors


def generate_umi(umi_length: int, rng: random.Random, is_duplex: bool = False) -> str:
    """Generate UMI string (duplex has two UMIs separated by dash)."""
    umi1 = generate_random_sequence(umi_length, rng)
    if is_duplex:
        umi2 = generate_random_sequence(umi_length, rng)
        return f"{umi1}-{umi2}"
    return umi1


def create_simplex_read(
    read_id: int,
    config: ConsensusConfig,
    rng: random.Random,
    np_rng: np.random.Generator,
    header: pysam.AlignmentHeader,
) -> Tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
    """Create a simplex consensus read pair (R1 and R2)."""

    # Sample depth from log-normal
    depth = max(1, int(np_rng.lognormal(config.depth_mu, config.depth_sigma)))

    # Generate sequences
    seq1 = generate_random_sequence(config.read_length, rng)
    seq2 = generate_random_sequence(config.read_length, rng)

    # Generate quality scores
    qual1 = generate_quality_scores(config.read_length, config.mean_quality, config.quality_std, rng)
    qual2 = generate_quality_scores(config.read_length, config.mean_quality, config.quality_std, rng)

    # Generate UMI
    umi = generate_umi(config.umi_length, rng, is_duplex=False)

    # Generate per-base arrays
    cd = generate_depth_array(config.read_length, depth, rng)
    ce = generate_error_array(config.read_length, config.mean_error_rate, rng)

    # Calculate summary stats
    cD = depth
    cM = min(cd)
    cE = sum(ce) / (len(ce) * depth) if depth > 0 else 0.0

    # Create R1
    r1 = pysam.AlignedSegment(header)
    r1.query_name = f"library:{read_id}"
    r1.query_sequence = seq1
    r1.query_qualities = pysam.qualitystring_to_array(qual1)
    r1.flag = 77  # paired, unmapped, mate unmapped, first in pair
    r1.set_tag("RG", "A")
    r1.set_tag("MI", str(read_id))
    r1.set_tag("RX", umi)
    r1.set_tag("cD", cD)
    r1.set_tag("cM", cM)
    r1.set_tag("cE", float(cE))
    r1.set_tag("cd", cd)
    r1.set_tag("ce", ce)

    # Create R2
    r2 = pysam.AlignedSegment(header)
    r2.query_name = f"library:{read_id}"
    r2.query_sequence = seq2
    r2.query_qualities = pysam.qualitystring_to_array(qual2)
    r2.flag = 141  # paired, unmapped, mate unmapped, second in pair
    r2.set_tag("RG", "A")
    r2.set_tag("MI", str(read_id))
    r2.set_tag("RX", umi)
    r2.set_tag("cD", cD)
    r2.set_tag("cM", cM)
    r2.set_tag("cE", float(cE))
    r2.set_tag("cd", cd)
    r2.set_tag("ce", ce)

    return r1, r2


def create_duplex_read(
    read_id: int,
    config: ConsensusConfig,
    rng: random.Random,
    np_rng: np.random.Generator,
    header: pysam.AlignmentHeader,
) -> Tuple[pysam.AlignedSegment, pysam.AlignedSegment]:
    """Create a duplex consensus read pair (R1 and R2)."""

    # Sample depths for A-strand, B-strand
    aD = max(1, int(np_rng.lognormal(config.depth_mu, config.depth_sigma)))
    bD = max(1, int(np_rng.lognormal(config.depth_mu, config.depth_sigma)))
    cD = aD + bD

    # Generate sequences
    seq1 = generate_random_sequence(config.read_length, rng)
    seq2 = generate_random_sequence(config.read_length, rng)

    # Generate quality scores (higher for duplex)
    qual1 = generate_quality_scores(config.read_length, config.mean_quality + 5, config.quality_std, rng)
    qual2 = generate_quality_scores(config.read_length, config.mean_quality + 5, config.quality_std, rng)

    # Generate UMI (duplex has paired UMI)
    umi = generate_umi(config.umi_length, rng, is_duplex=True)

    # Generate per-base arrays for A, B, and combined
    ad = generate_depth_array(config.read_length, aD, rng)
    bd = generate_depth_array(config.read_length, bD, rng)
    cd = [a + b for a, b in zip(ad, bd)]

    ae = generate_error_array(config.read_length, config.mean_error_rate, rng)
    be = generate_error_array(config.read_length, config.mean_error_rate, rng)
    ce = [a + b for a, b in zip(ae, be)]

    # Calculate summary stats
    aM = min(ad)
    bM = min(bd)
    cM = min(cd)

    aE = sum(ae) / (len(ae) * aD) if aD > 0 else 0.0
    bE = sum(be) / (len(be) * bD) if bD > 0 else 0.0
    cE = sum(ce) / (len(ce) * cD) if cD > 0 else 0.0

    # Create R1
    r1 = pysam.AlignedSegment(header)
    r1.query_name = f"library:{read_id}"
    r1.query_sequence = seq1
    r1.query_qualities = pysam.qualitystring_to_array(qual1)
    r1.flag = 77
    r1.set_tag("RG", "A")
    r1.set_tag("MI", str(read_id))
    r1.set_tag("RX", umi)
    # A-strand tags
    r1.set_tag("aD", aD)
    r1.set_tag("aM", aM)
    r1.set_tag("aE", float(aE))
    # B-strand tags
    r1.set_tag("bD", bD)
    r1.set_tag("bM", bM)
    r1.set_tag("bE", float(bE))
    # Combined tags
    r1.set_tag("cD", cD)
    r1.set_tag("cM", cM)
    r1.set_tag("cE", float(cE))
    # Per-base arrays
    r1.set_tag("ad", ad)
    r1.set_tag("bd", bd)
    r1.set_tag("cd", cd)
    r1.set_tag("ae", ae)
    r1.set_tag("be", be)
    r1.set_tag("ce", ce)

    # Create R2
    r2 = pysam.AlignedSegment(header)
    r2.query_name = f"library:{read_id}"
    r2.query_sequence = seq2
    r2.query_qualities = pysam.qualitystring_to_array(qual2)
    r2.flag = 141
    r2.set_tag("RG", "A")
    r2.set_tag("MI", str(read_id))
    r2.set_tag("RX", umi)
    r2.set_tag("aD", aD)
    r2.set_tag("aM", aM)
    r2.set_tag("aE", float(aE))
    r2.set_tag("bD", bD)
    r2.set_tag("bM", bM)
    r2.set_tag("bE", float(bE))
    r2.set_tag("cD", cD)
    r2.set_tag("cM", cM)
    r2.set_tag("cE", float(cE))
    r2.set_tag("ad", ad)
    r2.set_tag("bd", bd)
    r2.set_tag("cd", cd)
    r2.set_tag("ae", ae)
    r2.set_tag("be", be)
    r2.set_tag("ce", ce)

    return r1, r2


def generate_consensus_bam(
    config: ConsensusConfig,
    output_path: Path,
) -> None:
    """Generate synthetic consensus BAM file.

    Args:
        config: Generation configuration
        output_path: Path to output BAM file
    """
    rng = random.Random(config.seed)
    np_rng = np.random.default_rng(config.seed)

    print(f"Generating synthetic {config.consensus_type} consensus BAM:", file=sys.stderr)
    print(f"  Consensus reads: {config.n_consensus_reads:,}", file=sys.stderr)
    print(f"  Read length: {config.read_length}", file=sys.stderr)
    print(f"  Depth mu/sigma: {config.depth_mu}/{config.depth_sigma}", file=sys.stderr)

    # Create header
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "RG": [{"ID": "A", "SM": "synthetic", "LB": "library"}],
    })

    create_read = create_simplex_read if config.consensus_type == "simplex" else create_duplex_read

    with pysam.AlignmentFile(str(output_path), "wb", header=header) as bam:
        for read_id in range(config.n_consensus_reads):
            r1, r2 = create_read(read_id, config, rng, np_rng, header)
            bam.write(r1)
            bam.write(r2)

            if (read_id + 1) % 500_000 == 0:
                print(f"  Generated {read_id + 1:,} consensus reads...", file=sys.stderr)

    print(f"Wrote: {output_path}", file=sys.stderr)
    print(f"Total reads in BAM: {config.n_consensus_reads * 2:,} (paired)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic consensus BAM data for filter benchmarking"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Output BAM file path"
    )
    parser.add_argument(
        "-n", "--n-consensus-reads", type=int, default=2_000_000,
        help="Number of consensus read pairs to generate (default: 2M)"
    )
    parser.add_argument(
        "-l", "--read-length", type=int, default=289,
        help="Read length (default: 289)"
    )
    parser.add_argument(
        "-t", "--consensus-type", choices=["simplex", "duplex"], default="simplex",
        help="Type of consensus reads (default: simplex)"
    )
    parser.add_argument(
        "--depth-mu", type=float, default=2.5,
        help="Log-normal mu for depth distribution (default: 2.5)"
    )
    parser.add_argument(
        "--depth-sigma", type=float, default=0.8,
        help="Log-normal sigma for depth distribution (default: 0.8)"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42,
        help="Random seed (default: 42)"
    )

    args = parser.parse_args()

    config = ConsensusConfig(
        n_consensus_reads=args.n_consensus_reads,
        read_length=args.read_length,
        consensus_type=args.consensus_type,
        depth_mu=args.depth_mu,
        depth_sigma=args.depth_sigma,
        seed=args.seed,
    )

    generate_consensus_bam(config, args.output)


if __name__ == "__main__":
    main()
