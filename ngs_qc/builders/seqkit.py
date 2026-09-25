"""
SeqKit stats: read counts, length distribution, GC and Q20/Q30 per file.
Runs on all reads; parsed natively by MultiQC.
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample
from .base import Context, Job, compose, prepare_reads, q

TOOL = "SeqKit"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/seqkit"
    log     = f"{out_dir}/{sample.name}.seqkit.log"
    tsv     = f"{out_dir}/{sample.name}.seqkit_stats.tsv"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx)

    steps = prep + [
        f"module load {MODULES['seqkit']}",
        f"seqkit stats --all --tabular -j {ctx.threads} "
        + " ".join(q(r) for r in reads) + f" -o {q(tsv)}",
    ]
    return [Job(TOOL, sample.name, log, [tsv], compose(out_dir, log, steps, temp))]
