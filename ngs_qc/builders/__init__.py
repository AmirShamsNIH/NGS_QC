"""
Script builders for each QC tool.

Each per-sample tool module exposes ``TOOL`` and ``build(sample, ctx)``;
``PER_SAMPLE_BUILDERS`` fixes the order tools appear in the swarm files
and the XLSX report.  Only tools whose output MultiQC parses are
included.  MultiQC and the log report are study-level scripts.
"""

from . import (
    afterqc,
    fastp,
    fastq_screen,
    fastqc,
    insert_size,
    kat,
    kraken,
    sortmerna,
)
from .base import KRAKEN_SWARM, QC_SWARM, Context, Job
from .log_report import build_log_report_script
from .multiqc import build_multiqc_config, build_multiqc_script

PER_SAMPLE_BUILDERS = [
    kraken,        # + Bracken, chained in the same line
    fastqc,
    fastp,
    afterqc,
    insert_size,
    fastq_screen,
    sortmerna,
    kat,
]

__all__ = [
    "PER_SAMPLE_BUILDERS",
    "KRAKEN_SWARM",
    "QC_SWARM",
    "Context",
    "Job",
    "build_multiqc_config",
    "build_multiqc_script",
    "build_log_report_script",
]
