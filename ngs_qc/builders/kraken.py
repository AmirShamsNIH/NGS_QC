"""
Kraken2 taxonomic contamination screening, with Bracken chained after it.

Kraken2's per-read classification output is discarded (``--output
/dev/null``); only the summary report is kept.  Bracken then re-estimates
species-level abundance from that report.  Both run in one swarm line so
Bracken never starts before its Kraken2 report exists.

Single-cell mode: only R2 (cDNA) is classified.
"""

from __future__ import annotations

from ..config import BRACKEN_LEVEL, BRACKEN_THRESHOLD, MODULES
from ..fastq import Sample, is_gzipped
from .base import KRAKEN_SWARM, Context, Job, biological_reads, compose, prepare_reads, q

TOOL = "Kraken2"
BRACKEN_TOOL = "Bracken"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/kraken"
    br_dir  = f"{ctx.study_dir}/bracken"
    report  = f"{out_dir}/{sample.name}.kraken.txt"
    log     = f"{out_dir}/{sample.name}.kraken.log"

    reads_in = biological_reads(sample, ctx)
    prep, reads, temp = prepare_reads(sample, out_dir, ctx, reads_in)
    gz = " --gzip-compressed" if all(is_gzipped(r) for r in reads) else ""
    paired = " --paired" if len(reads) == 2 else ""

    steps = prep + [
        f"module load {MODULES['kraken2']}",
        f"kraken2 --threads {ctx.kraken_threads} --db {q(ctx.kraken_db)} "
        f"--report {q(report)} --output /dev/null{gz}{paired} "
        + " ".join(q(r) for r in reads),
    ]
    jobs = [Job(TOOL, sample.name, log, [report], None, KRAKEN_SWARM)]

    if ctx.bracken_read_length:
        br_log    = f"{br_dir}/{sample.name}.bracken.log"
        br_out    = f"{br_dir}/{sample.name}.bracken.txt"
        br_report = f"{br_dir}/{sample.name}.bracken_report.txt"
        # Bracken's log is written separately so its status is reported on
        # its own row; the Kraken2 log still captures the chain's exit.
        steps += [
            f"mkdir -p {q(br_dir)}",
            f"module load {MODULES['bracken']}",
            f"bracken -d {q(ctx.kraken_db)} -i {q(report)} -o {q(br_out)} "
            f"-w {q(br_report)} -r {ctx.bracken_read_length} "
            f"-l {BRACKEN_LEVEL} -t {BRACKEN_THRESHOLD} > {q(br_log)} 2>&1",
        ]
        jobs.append(Job(BRACKEN_TOOL, sample.name, br_log, [br_out, br_report],
                        None, KRAKEN_SWARM))

    jobs[0].cmd = compose(out_dir, log, steps, temp)
    return jobs
