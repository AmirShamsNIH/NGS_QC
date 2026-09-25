"""
FastQC per-base / per-read quality metrics on all reads (R1 and R2 in one
swarm line, so both mates share one log and one subset step).
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample, fastq_basename
from .base import Context, Job, compose, prepare_reads, q

TOOL = "FastQC"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/fastqc"
    log     = f"{out_dir}/{sample.name}.fastqc.log"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx)

    steps = prep + [
        f"module load {MODULES['fastqc']}",
        f"fastqc -o {q(out_dir)} -f fastq --threads {ctx.threads} "
        + " ".join(q(r) for r in reads),
    ]
    outputs = [f"{out_dir}/{fastq_basename(r)}_fastqc.html" for r in reads]
    return [Job(TOOL, sample.name, log, outputs, compose(out_dir, log, steps, temp))]
