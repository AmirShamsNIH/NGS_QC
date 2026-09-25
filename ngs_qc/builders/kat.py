"""
KAT k-mer frequency histogram.  Reveals contamination, GC bias, copy-number
artefacts and sequencing errors invisible to per-read tools.  KAT writes the
histogram to ``<prefix>`` (no extension), a plot to ``<prefix>.png`` and
``<prefix>.dist_analysis.json``, which MultiQC parses.

KAT 2.4.2 intermittently segfaults while binning k-mers (not reproducible
on re-run of the same input), so a failed run is retried once.  The retry
message is the marker ``collect_qc_logs.py`` uses to ignore the first
attempt's crash when the retry succeeds.

Single-cell mode: only R2 (cDNA) is analysed.
"""

from __future__ import annotations

from ..config import MODULES
from ..fastq import Sample
from .base import Context, Job, biological_reads, compose, prepare_reads, q

TOOL = "KAT"


def build(sample: Sample, ctx: Context) -> list[Job]:
    out_dir = f"{ctx.study_dir}/kat"
    prefix  = f"{out_dir}/{sample.name}.kat"
    log     = f"{out_dir}/{sample.name}.kat.log"
    prep, reads, temp = prepare_reads(sample, out_dir, ctx, biological_reads(sample, ctx))

    kat_cmd = (f"kat hist --threads {ctx.threads} --output_prefix {q(prefix)} "
               + " ".join(q(r) for r in reads))
    steps = prep + [
        f"module load {MODULES['kat']}",
        f"{{ {kat_cmd} || {{ echo 'KAT failed, retrying once'; {kat_cmd}; }}; }}",
    ]
    return [Job(TOOL, sample.name, log, [prefix], compose(out_dir, log, steps, temp))]
