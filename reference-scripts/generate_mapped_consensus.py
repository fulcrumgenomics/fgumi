#!/usr/bin/env python3
"""Generate synthetic mapped consensus BAM for filter benchmarking.

This script creates consensus BAMs with:
- Proper mapping coordinates (aligned to reference)
- Consensus tags (cD, cM, cE, cd, ce for simplex; adds aD, bD, etc. for duplex)
- Realistic quality scores and depth distributions

This allows benchmarking filter independently of the full pipeline.
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
class MappedConsensusConfig:
    """Configuration for mapped consensus generation."""

    n_consensus_reads: int = 2_000_000
    read_length: int = 150

    # Depth distribution (log-normal)
    depth_mu: float = 2.5
    depth_sigma: float = 0.8

    # Error rate distribution
    mean_error_rate: float = 0.001
    error_rate_std: float = 0.0005

    # Consensus type
    consensus_type: str = "simplex"  # "simplex" or "duplex"

    # Quality scores
    mean_quality: int = 35
    quality_std: float = 3.0

    # Fragment size for paired reads
    fragment_mean: int = 300
    fragment_std: int = 50

    # Output sort order: "coordinate" or "query-grouped"
    sort_order: str = "coordinate"

    seed: int = 42


def load_reference(fasta_path: Path) -> Tuple[dict[str, str], dict[str, int]]:
    """Load reference and return sequences and lengths."""
    seqs = {}
    lengths = {}
    with pysam.FastaFile(str(fasta_path)) as fa:
        for name in fa.references:
            seqs[name] = fa.fetch(name)
            lengths[name] = fa.get_reference_length(name)
    return seqs, lengths


def generate_quality_scores(length: int, mean_q: int, std_q: float, np_rng: np.random.Generator) -> List[int]:
    """Generate quality scores as integers."""
    scores = np_rng.normal(mean_q, std_q, length).astype(int)
    scores = np.clip(scores, 2, 40)
    return scores.tolist()


def generate_depth_array(length: int, base_depth: int, np_rng: np.random.Generator) -> List[int]:
    """Generate per-base depth array."""
    depths = base_depth + np_rng.integers(-1, 2, length)
    depths = np.maximum(depths, 1)
    return depths.tolist()


def generate_error_array(length: int, mean_error: float, np_rng: np.random.Generator) -> List[int]:
    """Generate per-base error count array."""
    mask = np_rng.random(length) < mean_error * 100
    errors = np.where(mask, np_rng.integers(1, 3, length), 0)
    return errors.tolist()


def sample_position(
    ref_lengths: dict[str, int],
    read_length: int,
    rng: random.Random,
) -> Tuple[str, int]:
    """Sample a valid position from reference."""
    # Get valid chromosomes
    valid_chroms = [(name, length) for name, length in ref_lengths.items()
                    if length > read_length * 2]

    if not valid_chroms:
        raise ValueError("No chromosomes long enough")

    # Weight by length
    total = sum(l for _, l in valid_chroms)
    weights = [l / total for _, l in valid_chroms]

    idx = rng.choices(range(len(valid_chroms)), weights=weights)[0]
    chrom, chrom_len = valid_chroms[idx]

    pos = rng.randint(0, chrom_len - read_length * 2)
    return chrom, pos


def create_simplex_read(
    read_id: int,
    chrom: str,
    pos: int,
    ref_seq: str,
    config: MappedConsensusConfig,
    rng: random.Random,
    np_rng: np.random.Generator,
    header: pysam.AlignmentHeader,
    is_r2: bool = False,
) -> pysam.AlignedSegment:
    """Create a simplex consensus read."""
    # Sample depth
    depth = max(1, int(np_rng.lognormal(config.depth_mu, config.depth_sigma)))

    # Get sequence from reference
    if is_r2:
        # R2 maps to opposite strand, downstream
        r2_pos = pos + config.fragment_mean
        seq = ref_seq[r2_pos:r2_pos + config.read_length]
        # Reverse complement for R2
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        seq = "".join(comp.get(b, "N") for b in reversed(seq))
        map_pos = r2_pos
        is_reverse = True
    else:
        seq = ref_seq[pos:pos + config.read_length]
        map_pos = pos
        is_reverse = False

    if len(seq) < config.read_length:
        # Pad with Ns if near end of chromosome
        seq = seq + "N" * (config.read_length - len(seq))

    # Generate quality scores
    quals = generate_quality_scores(config.read_length, config.mean_quality, config.quality_std, np_rng)

    # Generate per-base arrays
    cd = generate_depth_array(config.read_length, depth, np_rng)
    ce = generate_error_array(config.read_length, config.mean_error_rate, np_rng)

    # Summary stats
    cD = depth
    cM = min(cd)
    cE = sum(ce) / (len(ce) * depth) if depth > 0 else 0.0

    # Create read
    read = pysam.AlignedSegment(header)
    read.query_name = f"consensus:{read_id}"
    read.query_sequence = seq.upper()
    read.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))

    # Set flags
    flag = 0
    if is_r2:
        flag |= 0x80  # second in pair
        flag |= 0x10  # reverse strand
    else:
        flag |= 0x40  # first in pair

    flag |= 0x1  # paired
    flag |= 0x2  # proper pair
    read.flag = flag

    # Set mapping
    read.reference_id = header.get_tid(chrom)
    read.reference_start = map_pos
    read.mapping_quality = 60
    read.cigar = [(0, config.read_length)]  # All matches

    # Set tags
    read.set_tag("RG", "A")
    read.set_tag("MI", str(read_id))
    read.set_tag("cD", cD)
    read.set_tag("cM", cM)
    read.set_tag("cE", float(cE))
    read.set_tag("cd", cd)
    read.set_tag("ce", ce)

    return read


def create_duplex_read(
    read_id: int,
    chrom: str,
    pos: int,
    ref_seq: str,
    config: MappedConsensusConfig,
    rng: random.Random,
    np_rng: np.random.Generator,
    header: pysam.AlignmentHeader,
    is_r2: bool = False,
) -> pysam.AlignedSegment:
    """Create a duplex consensus read."""
    # Sample depths for A and B strands
    aD = max(1, int(np_rng.lognormal(config.depth_mu, config.depth_sigma)))
    bD = max(1, int(np_rng.lognormal(config.depth_mu, config.depth_sigma)))
    cD = aD + bD

    # Get sequence from reference
    if is_r2:
        r2_pos = pos + config.fragment_mean
        seq = ref_seq[r2_pos:r2_pos + config.read_length]
        comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        seq = "".join(comp.get(b, "N") for b in reversed(seq))
        map_pos = r2_pos
        is_reverse = True
    else:
        seq = ref_seq[pos:pos + config.read_length]
        map_pos = pos
        is_reverse = False

    if len(seq) < config.read_length:
        seq = seq + "N" * (config.read_length - len(seq))

    # Quality scores (higher for duplex)
    quals = generate_quality_scores(config.read_length, config.mean_quality + 5, config.quality_std, np_rng)

    # Per-base arrays
    ad = generate_depth_array(config.read_length, aD, np_rng)
    bd = generate_depth_array(config.read_length, bD, np_rng)
    cd = [a + b for a, b in zip(ad, bd)]

    ae = generate_error_array(config.read_length, config.mean_error_rate, np_rng)
    be = generate_error_array(config.read_length, config.mean_error_rate, np_rng)
    ce = [a + b for a, b in zip(ae, be)]

    # Summary stats
    aM = min(ad)
    bM = min(bd)
    cM = min(cd)

    aE = sum(ae) / (len(ae) * aD) if aD > 0 else 0.0
    bE = sum(be) / (len(be) * bD) if bD > 0 else 0.0
    cE = sum(ce) / (len(ce) * cD) if cD > 0 else 0.0

    # Create read
    read = pysam.AlignedSegment(header)
    read.query_name = f"consensus:{read_id}"
    read.query_sequence = seq.upper()
    read.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))

    flag = 0
    if is_r2:
        flag |= 0x80 | 0x10
    else:
        flag |= 0x40
    flag |= 0x1 | 0x2
    read.flag = flag

    read.reference_id = header.get_tid(chrom)
    read.reference_start = map_pos
    read.mapping_quality = 60
    read.cigar = [(0, config.read_length)]

    # Tags
    read.set_tag("RG", "A")
    read.set_tag("MI", str(read_id))
    # A-strand
    read.set_tag("aD", aD)
    read.set_tag("aM", aM)
    read.set_tag("aE", float(aE))
    # B-strand
    read.set_tag("bD", bD)
    read.set_tag("bM", bM)
    read.set_tag("bE", float(bE))
    # Combined
    read.set_tag("cD", cD)
    read.set_tag("cM", cM)
    read.set_tag("cE", float(cE))
    # Per-base arrays
    read.set_tag("ad", ad)
    read.set_tag("bd", bd)
    read.set_tag("cd", cd)
    read.set_tag("ae", ae)
    read.set_tag("be", be)
    read.set_tag("ce", ce)

    return read


def generate_mapped_consensus(
    config: MappedConsensusConfig,
    ref_path: Path,
    output_path: Path,
) -> None:
    """Generate mapped consensus BAM file."""
    rng = random.Random(config.seed)
    np_rng = np.random.default_rng(config.seed)

    print(f"Loading reference: {ref_path}", file=sys.stderr)
    ref_seqs, ref_lengths = load_reference(ref_path)
    print(f"  Loaded {len(ref_seqs)} sequences", file=sys.stderr)

    print(f"Generating {config.n_consensus_reads:,} {config.consensus_type} consensus reads...", file=sys.stderr)

    # Create header with reference info
    sort_order = config.sort_order
    if sort_order == "coordinate":
        hd_dict = {"VN": "1.6", "SO": "coordinate"}
    else:  # query-grouped
        hd_dict = {"VN": "1.6", "SO": "unsorted", "GO": "query"}

    header_dict = {
        "HD": hd_dict,
        "SQ": [{"SN": name, "LN": length} for name, length in ref_lengths.items()],
        "RG": [{"ID": "A", "SM": "synthetic", "LB": "library"}],
    }
    header = pysam.AlignmentHeader.from_dict(header_dict)

    create_read = create_simplex_read if config.consensus_type == "simplex" else create_duplex_read

    total_reads = 0
    reads = []
    progress_verb = "Generated" if sort_order == "coordinate" else "Written"

    with pysam.AlignmentFile(str(output_path), "wb", header=header) as bam:
        print(f"Writing {output_path} ...", file=sys.stderr)
        for read_id in range(config.n_consensus_reads):
            chrom, pos = sample_position(ref_lengths, config.read_length, rng)
            ref_seq = ref_seqs[chrom]

            # Create R1 and R2
            r1 = create_read(read_id, chrom, pos, ref_seq, config, rng, np_rng, header, is_r2=False)
            r2 = create_read(read_id, chrom, pos, ref_seq, config, rng, np_rng, header, is_r2=True)

            # Set mate info
            r1.next_reference_id = r2.reference_id
            r1.next_reference_start = r2.reference_start
            r2.next_reference_id = r1.reference_id
            r2.next_reference_start = r1.reference_start

            total_reads += 2

            if sort_order == "coordinate":
                # Collect for sorting later
                reads.append((r1.reference_id, r1.reference_start, r1))
                reads.append((r2.reference_id, r2.reference_start, r2))
            else:
                # Write immediately (query-grouped: R1 then R2)
                bam.write(r1)
                bam.write(r2)

            if (read_id + 1) % 500_000 == 0:
                print(f"  {progress_verb} {read_id + 1:,} consensus read pairs...", file=sys.stderr)

        if sort_order == "coordinate":
            print("Sorting reads by coordinate...", file=sys.stderr)
            reads.sort(key=lambda x: (x[0], x[1]))

            print(f"Writing sorted reads...", file=sys.stderr)
            for _, _, read in reads:
                bam.write(read)

    if sort_order == "coordinate":
        print("Indexing BAM...", file=sys.stderr)
        pysam.index(str(output_path))

    print(f"Done! Total reads: {total_reads:,}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic mapped consensus BAM for filter benchmarking"
    )
    parser.add_argument(
        "-r", "--reference", type=Path, required=True,
        help="Reference FASTA file"
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True,
        help="Output BAM file"
    )
    parser.add_argument(
        "-n", "--n-consensus-reads", type=int, default=2_000_000,
        help="Number of consensus read pairs (default: 2M)"
    )
    parser.add_argument(
        "-l", "--read-length", type=int, default=150,
        help="Read length (default: 150)"
    )
    parser.add_argument(
        "-t", "--consensus-type", choices=["simplex", "duplex"], default="simplex",
        help="Consensus type (default: simplex)"
    )
    parser.add_argument(
        "--depth-mu", type=float, default=2.5,
        help="Log-normal mu for depth (default: 2.5)"
    )
    parser.add_argument(
        "--depth-sigma", type=float, default=0.8,
        help="Log-normal sigma for depth (default: 0.8)"
    )
    parser.add_argument(
        "-s", "--seed", type=int, default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--sort-order", choices=["coordinate", "query-grouped"], default="coordinate",
        help="Output sort order: coordinate (sorted+indexed) or query-grouped (streamed, no index). Default: coordinate"
    )

    args = parser.parse_args()

    config = MappedConsensusConfig(
        n_consensus_reads=args.n_consensus_reads,
        read_length=args.read_length,
        consensus_type=args.consensus_type,
        depth_mu=args.depth_mu,
        depth_sigma=args.depth_sigma,
        sort_order=args.sort_order,
        seed=args.seed,
    )

    generate_mapped_consensus(config, args.reference, args.output)


if __name__ == "__main__":
    main()
