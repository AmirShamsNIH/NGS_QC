"""
Shared types and shell helpers for the per-tool swarm builders.

Every tool module exposes ``build(sample, ctx) -> list[Job]``.  A job is one
swarm line plus the log and output files that prove it succeeded; the same
expectations are written to ``qc_manifest.json`` so ``collect_qc_logs.py``
checks exactly what was scheduled instead of guessing from file globs.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass, field

from ..fastq import Sample, subset_cmd

# which swarm file a job goes into
KRAKEN_SWARM = "kraken"
QC_SWARM = "qc"


@dataclass(frozen=True)
class Context:
    """Run-wide settings shared by all builders."""

    study_dir: str
    single_cell: bool = False
    subset_reads: int = 0
    kraken_db: str = ""
    sortmerna_db: str = ""
    bracken_read_length: int | None = None
    threads: int = 10
    kraken_threads: int = 20


@dataclass
class Job:
    """One sample × tool unit of work.

    *cmd* is None for a report-only entry whose work is done by the
    preceding job's command (e.g. Bracken chained after Kraken2).
    """

    tool: str
    sample: str
    log: str
    outputs: list[str]
    cmd: str | None
    swarm: str = QC_SWARM
    extra: dict = field(default_factory=dict)


def q(path: str) -> str:
    return shlex.quote(path)


def compose(
    workdir: str,
    log: str,
    steps: list[str],
    cleanup: list[str] | None = None,
) -> str:
    """Build one swarm line.

    Steps are chained with ``&&`` so a failed subset / module load stops the
    tool, all output goes to *log*, temporary files in *cleanup* are always
    removed, and the line's exit status is the tool's exit status.
    """
    body = " && ".join(s for s in steps if s)
    line = f"mkdir -p {q(workdir)} && {{ {body} ; }} > {q(log)} 2>&1"
    if cleanup:
        rm = " ".join(q(p) for p in cleanup)
        line = f"{line}; rc=$?; rm -rf {rm}; (exit $rc)"
    return line


def prepare_reads(
    sample: Sample,
    workdir: str,
    ctx: Context,
    reads: list[str] | None = None,
) -> tuple[list[str], list[str], list[str]]:
    """Return (prep_steps, read_paths, temp_files) for a tool.

    In subset mode the first N reads of each file are written into
    *workdir* (one copy per tool, so parallel swarm lines never share a
    temp file) and deleted when the tool finishes.
    """
    reads = reads if reads is not None else sample.files
    if not ctx.subset_reads:
        return [], reads, []

    steps, paths = [], []
    for src in reads:
        mate = "R2" if src == sample.rev else "R1"
        dest = f"{workdir}/{sample.name}_subset_{mate}.fastq.gz"
        steps.append(subset_cmd(q(src), q(dest), ctx.subset_reads))
        paths.append(dest)
    return steps, paths, list(paths)


def biological_reads(sample: Sample, ctx: Context) -> list[str]:
    """Reads carrying biological sequence.

    In single-cell mode R1 holds cell barcodes + UMIs, so only R2 (cDNA)
    is used by contamination / k-mer tools, fastp and BBDuk.
    """
    if ctx.single_cell and sample.paired:
        return [sample.rev]
    return sample.files
