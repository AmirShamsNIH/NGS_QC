"""
fastp quality statistics in **report-only** mode: adapter trimming, quality
filtering and length filtering are disabled and reads go to /dev/null, so
nothing is modified.  The JSON report is consumed by MultiQC.
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample
from .base import Context, Job, compose, prepare_reads, q

TOOL = "fastp"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/fastp"
    log     = f"{out_dir}/{sample.name}.fastp.log"
    json    = f"{out_dir}/{sample.name}.fastp.json"
    html    = f"{out_dir}/{sample.name}.fastp.html"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx)

    inputs = f"-i {q(reads[0])}" + (f" -I {q(reads[1])}" if len(reads) == 2 else "")
    steps = prep + [
        f"module load {MODULES['fastp']}",
        f"fastp --thread {ctx.threads} "
        f"--disable_adapter_trimming --disable_quality_filtering "
        f"--disable_length_filtering "
        f"--json {q(json)} --html {q(html)} {inputs} --stdout > /dev/null",
    ]
    return [Job(TOOL, sample.name, log, [json], compose(out_dir, log, steps, temp))]
