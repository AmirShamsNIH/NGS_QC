"""
AfterQC in ``--qc_only`` mode: quality/bias/overlap report without writing
any filtered reads.  Each sample gets its own report folder.

``--seq_len_req 1`` disables AfterQC's default 35 bp minimum: with short
reads (e.g. 28 bp single-cell R1) every read would otherwise be "filtered"
and AfterQC crashes plotting the empty post-filter statistics.

Even in QC-only mode AfterQC creates empty ``good/`` and ``bad/`` folders in
the working directory, so it is run from its own output folder and those
folders are removed afterwards.
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample
from .base import Context, Job, compose, prepare_reads, q

TOOL = "AfterQC"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/afterqc/{sample.name}"
    log     = f"{out_dir}/{sample.name}.afterqc.log"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx)

    inputs = f"-1 {q(reads[0])}" + (f" -2 {q(reads[1])}" if len(reads) == 2 else "")
    steps = prep + [
        f"module load {MODULES['afterqc']}",
        f"cd {q(out_dir)}",
        f"after.py --qc_only --seq_len_req 1 {inputs} --report_output_folder {q(out_dir)}",
    ]
    # AfterQC names reports after the input file; any HTML in the
    # per-sample folder counts.
    return [Job(TOOL, sample.name, log, [f"{out_dir}/*.html"],
                compose(out_dir, log, steps,
                        temp + [f"{out_dir}/good", f"{out_dir}/bad"]))]
