"""
FastQ Screen multi-genome contamination check.  Both mates are screened in
one swarm line (independently, as FastQ Screen does for multiple files).

Single-cell mode: only R2 (cDNA) is screened; barcode/UMI reads in R1 do
not map to any genome.
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample, fastq_basename
from .base import Context, Job, biological_reads, compose, prepare_reads, q

TOOL = "FastQ Screen"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/fastq_screen"
    log     = f"{out_dir}/{sample.name}.fastq_screen.log"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx, biological_reads(sample, ctx))

    steps = prep + [
        f"module load {MODULES['fastq_screen']}",
        f"fastq_screen --force --aligner bowtie2 --subset 100000 "
        f"--threads {ctx.threads} --bowtie2 --local --outdir {q(out_dir)} "
        + " ".join(q(r) for r in reads),
    ]
    outputs = [f"{out_dir}/{fastq_basename(r)}_screen.txt" for r in reads]
    return [Job(TOOL, sample.name, log, outputs, compose(out_dir, log, steps, temp))]
