"""
SortMeRNA rRNA fraction detection on a 5 % BBTools sub-sample (or on the
first N reads in ``--subset`` mode).  The aligned-read log
``<sample>_rRNA.log`` reports the rRNA fraction and is parsed by MultiQC.

SortMeRNA 4.x refuses to start if its ``kvdb`` directory is not empty, so
the index/read databases from any previous run are cleared first, and the
sub-sampled reads are deleted afterwards.
"""

from __future__ import annotations

from ..config import MODULES, SORTMERNA_REFS, SORTMERNA_SUBSAMPLE_RATE
from ..fastq import Sample
from .base import Context, Job, compose, prepare_reads, q

TOOL = "SortMeRNA"


def build(sample: Sample, ctx: Context) -> list[Job]:
    work_dir = f"{ctx.study_dir}/sortmerna/{sample.name}"
    log      = f"{work_dir}/{sample.name}.sortmerna.log"
    aligned  = f"{work_dir}/{sample.name}_rRNA"
    stale    = [f"{work_dir}/{d}" for d in ("kvdb", "readb", "idx")]

    if ctx.subset_reads:
        prep, reads, temp = prepare_reads(sample, work_dir, ctx)
    else:
        reads = [f"{work_dir}/{sample.name}_subsample_R1.fastq.gz"]
        io = f"in={q(sample.fwd)} out={q(reads[0])}"
        if sample.paired:
            reads.append(f"{work_dir}/{sample.name}_subsample_R2.fastq.gz")
            io = (f"in1={q(sample.fwd)} in2={q(sample.rev)} "
                  f"out1={q(reads[0])} out2={q(reads[1])}")
        prep = [
            f"module load {MODULES['bbtools']}",
            f"bbtools reformat {io} samplerate={SORTMERNA_SUBSAMPLE_RATE} overwrite=t",
        ]
        temp = list(reads)

    ref_flags = " ".join(f"--ref {q(f'{ctx.sortmerna_db}/{r}')}" for r in SORTMERNA_REFS)
    steps = [f"rm -rf {' '.join(q(d) for d in stale)}"] + prep + [
        f"module load {MODULES['sortmerna']}",
        f"sortmerna --threads {ctx.threads} -m 40000 -v "
        f"--workdir {q(work_dir)} --aligned {q(aligned)} "
        + " ".join(f"--reads {q(r)}" for r in reads)
        + f" {ref_flags}",
    ]
    return [Job(TOOL, sample.name, log, [f"{aligned}.log"],
                compose(work_dir, log, steps, temp + stale))]
