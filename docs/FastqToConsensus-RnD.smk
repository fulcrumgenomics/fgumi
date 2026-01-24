"""
fgumi FASTQ to Consensus Pipeline (R&D Version)

This Snakemake pipeline implements the best-practice workflow for processing
paired-end FASTQ files through to filtered consensus sequences using fgumi.

This is the "R&D" version that generates intermediate files for debugging and
parameter exploration. See best-practice-consensus-pipeline.md for details.

Usage:
    snakemake -s FastqToConsensus-RnD.smk --cores 16

Requirements:
    - fgumi (version 0.1+)
    - bwa mem (version 0.7.17+)
    - snakemake (version 7.0+)
"""

# =============================================================================
# Configuration
# =============================================================================

# Sample configuration
SAMPLE = "my_sample"

# Input files
R1_FASTQ = "r1.fq.gz"
R2_FASTQ = "r2.fq.gz"

# Reference genome
REFERENCE = "ref/genome.fa"

# Read structures (8bp UMI on R1 and R2 for duplex)
READ_STRUCTURE_R1 = "8M+T"
READ_STRUCTURE_R2 = "8M+T"

# UMI grouping strategy: "adjacency" for simplex, "paired" for duplex
GROUPING_STRATEGY = "paired"

# Consensus caller: "simplex" or "duplex"
CONSENSUS_CALLER = "duplex"

# Filtering parameters
MIN_READS = "3,2,2"  # For duplex: [duplex, AB, BA]. For simplex: single value
MAX_READ_ERROR_RATE = 0.025
MAX_BASE_ERROR_RATE = 0.1
MIN_BASE_QUALITY = 40
MAX_NO_CALL_FRACTION = 0.2

# Threading configuration
THREADS_EXTRACT = 4
THREADS_ALIGN = 16
THREADS_GROUP = 8
THREADS_CONSENSUS = 8
THREADS_FILTER = 8
THREADS_SORT = 8

# =============================================================================
# Target Rule
# =============================================================================

rule all:
    input:
        f"{SAMPLE}.consensus.filtered.bam",
        f"{SAMPLE}.consensus.filtered.bam.bai",
        f"{SAMPLE}.family_sizes.txt",
        f"{SAMPLE}.grouping_metrics.txt"


# =============================================================================
# Phase 1: FASTQ to Grouped BAM
# =============================================================================

rule extract_umis:
    """Extract UMIs from FASTQ and create unmapped BAM."""
    input:
        r1 = R1_FASTQ,
        r2 = R2_FASTQ
    output:
        bam = temp(f"{SAMPLE}.unmapped.bam")
    params:
        rs1 = READ_STRUCTURE_R1,
        rs2 = READ_STRUCTURE_R2,
        sample = SAMPLE,
        library = SAMPLE
    threads: THREADS_EXTRACT
    resources:
        mem_mb = 2000
    shell:
        """
        fgumi extract \
            --inputs {input.r1} {input.r2} \
            --read-structures {params.rs1} {params.rs2} \
            --sample {params.sample} \
            --library {params.library} \
            --output {output.bam} \
            --threads {threads} \
            --compression-level 1
        """


rule align_reads:
    """Align reads with bwa mem and restore tags with zipper."""
    input:
        unmapped = f"{SAMPLE}.unmapped.bam",
        ref = REFERENCE
    output:
        bam = temp(f"{SAMPLE}.aligned.bam")
    threads: THREADS_ALIGN
    resources:
        mem_mb = 16000
    shell:
        """
        fgumi fastq --input {input.unmapped} --threads 4 \
            | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - \
            | fgumi zipper \
                --unmapped {input.unmapped} \
                --reference {input.ref} \
                --output {output.bam} \
                --threads 4 \
                --compression-level 1
        """


rule group_reads:
    """Group reads by UMI using directed adjacency."""
    input:
        bam = f"{SAMPLE}.aligned.bam"
    output:
        bam = temp(f"{SAMPLE}.grouped.bam"),
        histogram = f"{SAMPLE}.family_sizes.txt",
        metrics = f"{SAMPLE}.grouping_metrics.txt"
    params:
        strategy = GROUPING_STRATEGY
    threads: THREADS_GROUP
    resources:
        mem_mb = 8000
    shell:
        """
        fgumi group \
            --input {input.bam} \
            --output {output.bam} \
            --strategy {params.strategy} \
            --edits 1 \
            --family-size-histogram {output.histogram} \
            --grouping-metrics {output.metrics} \
            --threads {threads} \
            --compression-level 1
        """


# =============================================================================
# Phase 2: Consensus Calling and Filtering
# =============================================================================

rule call_consensus:
    """Call consensus sequences from grouped reads."""
    input:
        bam = f"{SAMPLE}.grouped.bam"
    output:
        bam = temp(f"{SAMPLE}.consensus.unmapped.bam")
    params:
        caller = CONSENSUS_CALLER
    threads: THREADS_CONSENSUS
    resources:
        mem_mb = 8000
    run:
        if params.caller == "simplex":
            shell("""
                fgumi simplex \
                    --input {input.bam} \
                    --output {output.bam} \
                    --min-reads 1 \
                    --min-input-base-quality 20 \
                    --output-per-base-tags \
                    --threads {threads} \
                    --compression-level 1
            """)
        else:
            shell("""
                fgumi duplex \
                    --input {input.bam} \
                    --output {output.bam} \
                    --min-reads 1 \
                    --min-input-base-quality 20 \
                    --output-per-base-tags \
                    --threads {threads} \
                    --compression-level 1
            """)


rule align_consensus:
    """Re-align consensus reads."""
    input:
        unmapped = f"{SAMPLE}.consensus.unmapped.bam",
        ref = REFERENCE
    output:
        bam = temp(f"{SAMPLE}.consensus.aligned.bam")
    threads: THREADS_ALIGN
    resources:
        mem_mb = 8000
    shell:
        """
        fgumi fastq --input {input.unmapped} --threads 2 \
            | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - \
            | fgumi zipper \
                --unmapped {input.unmapped} \
                --reference {input.ref} \
                --output {output.bam} \
                --threads 2 \
                --compression-level 1
        """


rule filter_consensus:
    """Filter consensus reads and sort to coordinate order."""
    input:
        bam = f"{SAMPLE}.consensus.aligned.bam",
        ref = REFERENCE
    output:
        bam = f"{SAMPLE}.consensus.filtered.bam"
    params:
        min_reads = MIN_READS,
        max_read_error = MAX_READ_ERROR_RATE,
        max_base_error = MAX_BASE_ERROR_RATE,
        min_base_qual = MIN_BASE_QUALITY,
        max_no_call = MAX_NO_CALL_FRACTION,
        caller = CONSENSUS_CALLER
    threads: THREADS_FILTER
    resources:
        mem_mb = 8000
    run:
        extra_args = ""
        if params.caller == "duplex":
            extra_args = "--require-single-strand-agreement"

        shell("""
            fgumi filter \
                --input {input.bam} \
                --output {output.bam} \
                --ref {input.ref} \
                --min-reads {params.min_reads} \
                --max-read-error-rate {params.max_read_error} \
                --max-base-error-rate {params.max_base_error} \
                --min-base-quality {params.min_base_qual} \
                --max-no-call-fraction {params.max_no_call} \
                --reverse-per-base-tags \
                --sort-order coordinate \
                --threads {threads} \
                --compression-level 6 \
                {extra_args}
        """)


rule index_bam:
    """Create BAM index for final output."""
    input:
        bam = f"{SAMPLE}.consensus.filtered.bam"
    output:
        bai = f"{SAMPLE}.consensus.filtered.bam.bai"
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        fgumi sort \
            --input {input.bam} \
            --output {input.bam}.tmp \
            --order coordinate \
            --write-index \
            --threads 1
        mv {input.bam}.tmp {input.bam}
        mv {input.bam}.tmp.bai {output.bai}
        """


# =============================================================================
# Optional: Collect Duplex Metrics
# =============================================================================

rule duplex_metrics:
    """Collect duplex QC metrics (only for duplex workflow)."""
    input:
        bam = f"{SAMPLE}.grouped.bam"
    output:
        prefix = f"{SAMPLE}.duplex_metrics"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """
        fgumi duplex-metrics \
            --input {input.bam} \
            --output {output.prefix} \
            --description "{SAMPLE}"
        """
