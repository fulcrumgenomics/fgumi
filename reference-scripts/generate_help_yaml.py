#!/usr/bin/env python3
"""Generate YAML help files for fgumi binaries.

This script parses the --help output from fgumi binaries and generates
structured YAML files mapping canonical option names to actual CLI options.

This enables version-aware argument selection in Snakemake rules, since
different fgumi versions may have different option names (e.g.,
--min-base-quality vs --min-input-base-quality).

Usage:
    python generate_help_yaml.py bin/fgumi-abc1234 bin/fgumi-abc1234.help.yaml
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import yaml


# Canonical option names mapped to known long-form variants (in order of preference)
# Only long-form options are used for clarity and consistency
CANONICAL_OPTIONS = {
    # Threading options
    "threads": ["--threads"],
    # Consensus calling options
    "input": ["--input"],
    "output": ["--output"],
    "min_reads": ["--min-reads"],
    "min_input_base_quality": ["--min-input-base-quality", "--min-base-quality"],
    "min_base_quality": ["--min-base-quality", "--min-input-base-quality"],
    "read_name_prefix": ["--read-name-prefix"],
    "read_group_id": ["--read-group-id"],
    "output_per_base_tags": ["--output-per-base-tags"],
    "stats": ["--stats"],
    "rejects": ["--rejects"],
    "trim": ["--trim"],
    "max_reads_per_strand": ["--max-reads-per-strand"],
    # Filter options
    "ref": ["--ref"],
    "require_single_strand_agreement": ["--require-single-strand-agreement"],
    "regenerate_tags": ["--regenerate-tags"],
    # Group options
    "strategy": ["--strategy"],
    "edits": ["--edits"],
    "min_map_q": ["--min-map-q"],
    "family_size_histogram": ["--family-size-histogram"],
    # Correct options
    "revcomp": ["--revcomp"],
    "umi_tag": ["--umi-tag"],
    "umi_files": ["--umi-files"],
    "max_mismatches": ["--max-mismatches"],
    "min_distance": ["--min-distance"],
    "metrics": ["--metrics"],
    # Extract options
    "inputs": ["--inputs"],
    "read_structures": ["--read-structures"],
    "sample": ["--sample"],
    "library": ["--library"],
    # Zipper options
    "unmapped": ["--unmapped"],
    "reference": ["--reference"],
    "tags_to_reverse": ["--tags-to-reverse"],
    "tags_to_revcomp": ["--tags-to-revcomp"],
    # Compare options
    "mode": ["--mode"],
    "ignore_order": ["--ignore-order"],
}


def parse_options_from_help(help_text: str) -> set[str]:
    """Extract option names from --help output."""
    options = set()
    # Match patterns like:
    # -i, --input <INPUT>
    # --threads <THREADS>
    # --output-per-base-tags
    pattern = r"^\s+(-[a-zA-Z],\s+)?(--[a-z][a-z0-9-]*)"
    for line in help_text.split("\n"):
        match = re.match(pattern, line)
        if match:
            # Add long option
            options.add(match.group(2))
            # Add short option if present
            if match.group(1):
                short = match.group(1).strip().rstrip(",")
                options.add(short)
    return options


def find_best_match(canonical_key: str, actual_options: set[str]) -> str | None:
    """Find the best matching option for a canonical key.

    Only matches exact options from the known variants list.
    This avoids false positives from fuzzy matching.
    """
    variants = CANONICAL_OPTIONS.get(canonical_key, [])

    # Only use exact matches from the known variants list
    for variant in variants:
        if variant in actual_options:
            return variant

    return None


def get_subcommands(binary: Path) -> list[str]:
    """Get list of subcommands from the binary."""
    result = subprocess.run(
        [str(binary), "--help"],
        capture_output=True,
        text=True,
    )

    # Parse subcommands from help output
    # They appear after "Commands:" line
    subcommands = []
    in_commands = False
    for line in result.stdout.split("\n"):
        if line.strip() == "Commands:":
            in_commands = True
            continue
        if in_commands:
            if line.strip() == "" or line.startswith("Options:"):
                break
            # Parse command name (first word)
            parts = line.strip().split()
            if parts and not parts[0].startswith("-"):
                cmd = parts[0]
                if cmd != "help":  # Skip the help subcommand
                    subcommands.append(cmd)

    return subcommands


def get_subcommand_help(binary: Path, subcommand: str) -> str:
    """Get help text for a specific subcommand."""
    result = subprocess.run(
        [str(binary), subcommand, "--help"],
        capture_output=True,
        text=True,
    )
    return result.stdout


def generate_help_yaml(binary: Path, output: Path) -> dict:
    """Generate YAML help file for a binary."""
    data = {"subcommands": {}}

    subcommands = get_subcommands(binary)

    for subcommand in subcommands:
        help_text = get_subcommand_help(binary, subcommand)
        actual_options = parse_options_from_help(help_text)

        # Build canonical -> actual option mapping
        options = {}
        for canonical_key in CANONICAL_OPTIONS:
            match = find_best_match(canonical_key, actual_options)
            if match:
                options[canonical_key] = match

        data["subcommands"][subcommand] = {"options": options}

    # Write YAML
    with open(output, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=True)

    return data


def main():
    parser = argparse.ArgumentParser(
        description="Generate YAML help file for fgumi binary"
    )
    parser.add_argument("binary", type=Path, help="Path to fgumi binary")
    parser.add_argument("output", type=Path, help="Output YAML file path")

    args = parser.parse_args()

    if not args.binary.exists():
        print(f"Error: Binary not found: {args.binary}", file=sys.stderr)
        sys.exit(1)

    if not args.binary.is_file():
        print(f"Error: Not a file: {args.binary}", file=sys.stderr)
        sys.exit(1)

    generate_help_yaml(args.binary, args.output)
    print(f"Generated: {args.output}")


if __name__ == "__main__":
    main()
