"""
Generate synthetic BAM data for CorrectUmis benchmarking.

This script creates:
1. A UMI whitelist with N random UMIs
2. A BAM file with reads containing UMIs at various edit distances
3. A truth file mapping each read to its expected correction

Error distribution (configurable):
- 40% exact matches (edit distance 0)
- 30% 1 error (edit distance 1, should correct)
- 20% 2 errors (edit distance 2, may or may not correct based on settings)
- 10% 3+ errors (should not correct)

This module can be used as a Snakemake script or imported as a library.
"""

import random
from dataclasses import dataclass
from pathlib import Path
from typing import List
from typing import Tuple

import pysam


@dataclass
class SyntheticConfig:
    """Configuration for synthetic data generation."""

    n_umis: int
    n_reads: int
    umi_length: int
    seed: int
    # Error distribution (must sum to 1.0)
    exact_frac: float = 0.40
    edit1_frac: float = 0.30
    edit2_frac: float = 0.20
    multi_frac: float = 0.10


def generate_random_umi(length: int, rng: random.Random) -> str:
    """Generate a random UMI sequence."""
    return "".join(rng.choices("ACGT", k=length))


def introduce_errors(umi: str, n_errors: int, rng: random.Random) -> str:
    """Introduce N substitution errors into a UMI."""
    if n_errors == 0:
        return umi

    umi_list = list(umi)
    positions = rng.sample(range(len(umi)), min(n_errors, len(umi)))

    for pos in positions:
        current = umi_list[pos]
        alternatives = [b for b in "ACGT" if b != current]
        umi_list[pos] = rng.choice(alternatives)

    return "".join(umi_list)


def generate_whitelist(config: SyntheticConfig, rng: random.Random) -> List[str]:
    """Generate a list of unique UMIs for the whitelist."""
    umis: set[str] = set()
    while len(umis) < config.n_umis:
        umis.add(generate_random_umi(config.umi_length, rng))
    return sorted(umis)


def generate_read(
    read_id: int,
    umi: str,
    error_type: str,
    rng: random.Random,
    _config: SyntheticConfig,  # noqa: ARG001
) -> Tuple[pysam.AlignedSegment, str, str]:
    """
    Generate a single synthetic read with UMI.

    Args:
        read_id: Unique read identifier
        umi: The true UMI from whitelist
        error_type: One of 'exact', 'edit1', 'edit2', 'multi'
        rng: Random number generator
        config: Synthetic data configuration

    Returns:
        Tuple of (read, observed_umi, expected_correction)
    """
    # Determine number of errors based on error type
    if error_type == "exact":
        n_errors = 0
    elif error_type == "edit1":
        n_errors = 1
    elif error_type == "edit2":
        n_errors = 2
    else:  # multi
        n_errors = rng.randint(3, 5)

    observed_umi = introduce_errors(umi, n_errors, rng)

    # Expected correction: original UMI if correctable, else observed
    # Typically CorrectUmis uses max-mismatches=1
    if n_errors <= 1:
        expected = umi
    else:
        expected = observed_umi  # Not correctable

    # Create read with random sequence
    read = pysam.AlignedSegment()
    read.query_name = f"read_{read_id:08d}"
    read.query_sequence = generate_random_umi(100, rng)  # Random 100bp sequence
    read.query_qualities = pysam.qualitystring_to_array("I" * 100)
    read.flag = 4  # Unmapped
    read.set_tag("BC", observed_umi)

    return read, observed_umi, expected


def generate_synthetic_data(
    config: SyntheticConfig,
    output_bam: Path,
    output_umis: Path,
    output_truth: Path,
) -> None:
    """
    Generate complete synthetic dataset.

    Args:
        config: Configuration for data generation
        output_bam: Path to output BAM file
        output_umis: Path to output UMI whitelist file
        output_truth: Path to output truth TSV file
    """
    rng = random.Random(config.seed)

    # Generate whitelist
    whitelist = generate_whitelist(config, rng)

    # Write whitelist
    with open(output_umis, "w") as f:
        for umi in whitelist:
            f.write(f"{umi}\n")

    # Determine error type for each read
    n_exact = int(config.n_reads * config.exact_frac)
    n_edit1 = int(config.n_reads * config.edit1_frac)
    n_edit2 = int(config.n_reads * config.edit2_frac)
    n_multi = config.n_reads - n_exact - n_edit1 - n_edit2

    error_types = (
        ["exact"] * n_exact + ["edit1"] * n_edit1 + ["edit2"] * n_edit2 + ["multi"] * n_multi
    )
    rng.shuffle(error_types)

    # Generate BAM
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "RG": [{"ID": "synthetic", "SM": "synthetic", "LB": "lib1"}],
    })

    truth_records: List[dict[str, str | None]] = []

    with pysam.AlignmentFile(str(output_bam), "wb", header=header) as bam:
        for read_id, error_type in enumerate(error_types):
            # Pick random UMI from whitelist
            true_umi = rng.choice(whitelist)

            read, observed, expected = generate_read(read_id, true_umi, error_type, rng, config)
            read.set_tag("RG", "synthetic")
            bam.write(read)

            truth_records.append({
                "read_id": read.query_name,
                "true_umi": true_umi,
                "observed_umi": observed,
                "expected_correction": expected,
                "error_type": error_type,
            })

    # Write truth file
    with open(output_truth, "w") as f:
        f.write("read_id\ttrue_umi\tobserved_umi\texpected_correction\terror_type\n")
        for rec in truth_records:
            f.write(
                f"{rec['read_id']}\t{rec['true_umi']}\t{rec['observed_umi']}\t"
                f"{rec['expected_correction']}\t{rec['error_type']}\n"
            )

    print(f"Generated {config.n_reads} reads with {config.n_umis} UMIs")
    print(f"  - Exact matches: {n_exact} ({config.exact_frac * 100:.0f}%)")
    print(f"  - Edit distance 1: {n_edit1} ({config.edit1_frac * 100:.0f}%)")
    print(f"  - Edit distance 2: {n_edit2} ({config.edit2_frac * 100:.0f}%)")
    print(f"  - Multi-error: {n_multi} ({config.multi_frac * 100:.0f}%)")


# Snakemake script entry point
if "snakemake" in dir():
    _snakemake_config = SyntheticConfig(
        n_umis=snakemake.params.n_umis,  # type: ignore[name-defined] # noqa: F821
        n_reads=snakemake.params.n_reads,  # type: ignore[name-defined] # noqa: F821
        umi_length=snakemake.params.umi_length,  # type: ignore[name-defined] # noqa: F821
        seed=snakemake.params.seed,  # type: ignore[name-defined] # noqa: F821
    )
    generate_synthetic_data(
        _snakemake_config,
        Path(snakemake.output.bam),  # type: ignore[name-defined] # noqa: F821
        Path(snakemake.output.umis),  # type: ignore[name-defined] # noqa: F821
        Path(snakemake.output.truth),  # type: ignore[name-defined] # noqa: F821
    )

elif __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate synthetic BAM data for CorrectUmis benchmarking"
    )
    parser.add_argument("--n-umis", type=int, default=96, help="Number of UMIs")
    parser.add_argument("--n-reads", type=int, default=10000, help="Number of reads")
    parser.add_argument("--umi-length", type=int, default=18, help="UMI length")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--output-bam", type=Path, required=True, help="Output BAM")
    parser.add_argument("--output-umis", type=Path, required=True, help="Output UMI whitelist")
    parser.add_argument("--output-truth", type=Path, required=True, help="Output truth TSV")

    args = parser.parse_args()

    config = SyntheticConfig(
        n_umis=args.n_umis,
        n_reads=args.n_reads,
        umi_length=args.umi_length,
        seed=args.seed,
    )

    generate_synthetic_data(
        config,
        args.output_bam,
        args.output_umis,
        args.output_truth,
    )
