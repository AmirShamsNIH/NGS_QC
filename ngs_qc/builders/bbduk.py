"""
BBDuk adapter and PhiX/spike-in content, **report only**: matched reads are
counted per reference sequence (``stats=``) and nothing is written out.
Uses the adapter and PhiX references bundled with BBTools.

BBDuk's own stderr goes to ``<sample>.bbduk_stderr.txt``: MultiQC only
recognises it if "Executing jgi.BBDuk" is in the first two lines, which the
``module load`` messages in the job log would push down.

Single-cell mode: only R2 (cDNA) is analysed; R1 holds barcodes + UMIs.
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample
from .base import Context, Job, biological_reads, compose, prepare_reads, q

TOOL = "BBDuk"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/bbduk"
    log     = f"{out_dir}/{sample.name}.bbduk.log"
    stats   = f"{out_dir}/{sample.name}.bbduk_stats.txt"
    stderr  = f"{out_dir}/{sample.name}.bbduk_stderr.txt"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx, biological_reads(sample, ctx))

    inputs = f"in1={q(reads[0])}" + (f" in2={q(reads[1])}" if len(reads) == 2 else "")
    steps = prep + [
        f"module load {MODULES['bbtools']}",
        f"bbtools bbduk {inputs} ref=adapters,phix k=23 hdist=1 "
        f"threads={ctx.threads} stats={q(stats)} 2> {q(stderr)}",
    ]
    return [Job(TOOL, sample.name, log, [stats, stderr], compose(out_dir, log, steps, temp))]
