"""
Insert-size distribution from read-pair overlap with BBMerge (``ihist=``).
Nothing is merged or written besides the histogram, which MultiQC parses
with its BBTools module.

Only paired-end bulk libraries: skipped for single-end samples, in
single-cell mode (R1 is barcode/UMI, so mates never overlap), and when
either mate's reads are shorter than ``INSERT_SIZE_MIN_READ_LEN``.  Pairs whose
insert is longer than the combined read length do not overlap and are not
counted, so the histogram under-represents long inserts.
"""

from __future__ import annotations

import logging

from ..config import INSERT_SIZE_MIN_READ_LEN, INSERT_SIZE_READS, MODULES
from ..fastq import Sample, first_read_length
from .base import Context, Job, compose, prepare_reads, q

TOOL = "Insert size"

logger = logging.getLogger(__name__)


def build(sample: Sample, ctx: Context) -> list[Job]:
    if not sample.paired or ctx.single_cell:
        return []
    lengths = [first_read_length(f) for f in sample.files]
    if any(n is not None and n < INSERT_SIZE_MIN_READ_LEN for n in lengths):
        logger.info("Insert size skipped for %s: read lengths %s, mates cannot "
                    "overlap (min %d bp)", sample.name, lengths,
                    INSERT_SIZE_MIN_READ_LEN)
        return []
    out_dir = f"{ctx.study_dir}/insert_size"
    log     = f"{out_dir}/{sample.name}.bbmerge.log"
    ihist   = f"{out_dir}/{sample.name}.ihist.txt"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx)

    limit = f" reads={INSERT_SIZE_READS}" if INSERT_SIZE_READS else ""
    steps = prep + [
        f"module load {MODULES['bbtools']}",
        f"bbtools bbmerge in1={q(reads[0])} in2={q(reads[1])} "
        f"ihist={q(ihist)} threads={ctx.threads}{limit}",
    ]
    return [Job(TOOL, sample.name, log, [ihist], compose(out_dir, log, steps, temp))]
