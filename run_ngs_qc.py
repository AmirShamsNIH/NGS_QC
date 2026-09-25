#!/usr/bin/env python3
"""
run_ngs_qc.py – NGS Quality Control pipeline (v2)
==================================================

Scans a directory of FASTQ files, resolves paired-end vs single-end layout
per sample, and generates ready-to-submit HPC scripts for a comprehensive,
**read-only** QC run.  No FASTQ files are modified or filtered at any stage.

Files written per study:

  1. ``<output>/<study>_kraken.swarm`` and ``<output>/<study>_qc.swarm``
        One shell command per line; submitted as two swarm jobs so every
        sample × tool combination runs in parallel on the cluster.  Kraken2
        (+ Bracken) gets its own swarm because only it needs the
        large-memory allocation for the database.

  2. ``<output>/<study>_multiqc.sh`` + ``<output>/<study>/multiqc_config.yaml``
        sbatch script that aggregates all QC outputs with MultiQC.

  3. ``<output>/<study>_log_report.sh``
        sbatch script that runs ``collect_qc_logs.py`` after all QC jobs
        finish and writes a multi-tab XLSX report showing which samples
        passed or failed per tool.

  4. ``<output>/<study>/qc_manifest.json``
        Every scheduled sample × tool job with its log and expected
        outputs; read by ``collect_qc_logs.py``.

  5. ``<output>/<study>_execution.sh``
        Wrapper that submits both swarms, then chains MultiQC and the log
        report sbatch jobs via ``--dependency afterany:<kraken>:<qc>``.

QC tools
--------
  Kraken2 + Bracken  – taxonomic contamination screening / species abundance
  FastQC             – per-base / per-read quality metrics
  fastp              – quality statistics (report-only; no reads are modified)
  SeqKit stats       – read counts, lengths, GC, Q20/Q30
  AfterQC            – quality / bias report (--qc_only)
  BBDuk              – adapter and PhiX content (report-only)
  Insert size        – BBMerge overlap-based insert-size histogram (paired bulk)
  FastQ Screen       – multi-genome contamination check
  SortMeRNA          – rRNA fraction detection
  KAT                – k-mer histogram / contamination

Usage
-----
::

    python run_ngs_qc.py INPUT_DIR OUTPUT_DIR [OPTIONS]

    # Then execute the generated wrapper on the cluster:
    bash <output>/<study>_execution.sh

Positional arguments
--------------------
INPUT_DIR
    Directory containing raw FASTQ files (``*.fastq[.gz]`` / ``*.fq[.gz]``).
OUTPUT_DIR
    Root directory where all QC results and generated scripts are written.

Options
-------
--study-name NAME
    Label used for sub-directories and file prefixes.
    Defaults to the basename of INPUT_DIR.
--single-cell
    Contamination / k-mer tools (Kraken2, FastQ Screen, KAT) process
    only R2 (cDNA) reads, and insert size is skipped.  FastQC, fastp and
    the other per-read tools still run on all reads.
--kraken-db PATH
    Kraken2 / Bracken database directory (default: KRAKEN2_DB in config).
--sortmerna-db PATH
    SortMeRNA rRNA reference database directory (default: SORTMERNA_DB).
--bracken-read-length N
    Bracken k-mer distribution to use.  Defaults to the longest one
    available in the Kraken2 DB that does not exceed the read length.
--subset [N]
    Run all tools on the first N reads per sample (default when flag is
    present: 10 000).  Useful for rapid pipeline validation.
--log-level {DEBUG,INFO,WARNING,ERROR}
    Verbosity level (default: INFO).
--version
    Show version and exit.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import re
import shlex
import sys
from dataclasses import asdict
from pathlib import Path

from ngs_qc import __version__
from ngs_qc.builders import (
    KRAKEN_SWARM,
    PER_SAMPLE_BUILDERS,
    QC_SWARM,
    Context,
    Job,
    build_log_report_script,
    build_multiqc_config,
    build_multiqc_script,
)
from ngs_qc.config import (
    KRAKEN2_DB,
    KRAKEN_SWARM_MEMORY_GB,
    KRAKEN_SWARM_THREADS,
    MULTIQC_CPUS,
    MULTIQC_LSCRATCH_GB,
    MULTIQC_MEMORY,
    SORTMERNA_DB,
    SWARM_MEMORY_GB,
    SWARM_PARTITION,
    SWARM_THREADS,
    SWARM_TIME,
)
from ngs_qc.fastq import Sample, describe_layout, discover_samples, first_read_length

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s – %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("run_ngs_qc")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_execution_script(
    study: str,
    swarm_files: dict[str, str],
    multiqc_file: str,
    log_report_file: str,
    study_dir: str,
) -> str:
    """Return a bash wrapper that submits both swarms, then chains MultiQC
    and the log-report job (both depend on both swarm job IDs)."""
    logdir = shlex.quote(f"{study_dir}/swarm_logs/")

    def _swarm(var: str, mem_gb: int, threads: int, path: str) -> list[str]:
        return [
            f"{var}=$(swarm "
            f"-g {mem_gb} "
            f"-t {threads} "
            f"--time {SWARM_TIME} "
            f"--logdir {logdir} "
            f"--partition {SWARM_PARTITION} "
            f"--sbatch '--mail-type=FAIL' "
            f"-f {shlex.quote(path)})",
            f'[[ -n "${var}" ]] || {{ echo "swarm submission failed: {path}" >&2; exit 1; }}',
        ]

    deps = "--dependency afterany:$kraken_jobID:$qc_jobID"

    return "\n".join([
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        f"# Auto-generated execution wrapper for study: {study}",
        f"# Generated by run_ngs_qc.py v{__version__}",
        "",
        "# Step 1a – Kraken2 + Bracken swarm (large-memory: loads the Kraken2 DB)",
        *_swarm("kraken_jobID", KRAKEN_SWARM_MEMORY_GB, KRAKEN_SWARM_THREADS,
                swarm_files[KRAKEN_SWARM]),
        'echo "Kraken2 swarm submitted: $kraken_jobID"',
        "",
        "# Step 1b – all other per-sample QC tools",
        *_swarm("qc_jobID", SWARM_MEMORY_GB, SWARM_THREADS, swarm_files[QC_SWARM]),
        'echo "QC swarm submitted: $qc_jobID"',
        "",
        "# Step 2 – MultiQC after both swarms complete",
        f"sbatch "
        f"--mem={MULTIQC_MEMORY} "
        f"--cpus-per-task={MULTIQC_CPUS} "
        f"--partition={SWARM_PARTITION} "
        f"--gres=lscratch:{MULTIQC_LSCRATCH_GB} "
        f"--time={SWARM_TIME} "
        f"{deps} "
        f"{shlex.quote(multiqc_file)}",
        "",
        "# Step 3 – collect logs and generate XLSX execution report",
        f"sbatch "
        f"--mem=4G "
        f"--cpus-per-task=2 "
        f"--time=00:30:00 "
        f"--partition={SWARM_PARTITION} "
        f"{deps} "
        f"{shlex.quote(log_report_file)}",
        "",
    ])


def _pick_bracken_length(kraken_db: str, samples: list[Sample],
                         requested: int | None) -> int | None:
    """Choose the Bracken k-mer distribution for this run.

    Uses *requested* if given, else the longest distribution in the DB that
    does not exceed the observed read length.  Returns None (Bracken
    skipped) when no distribution files exist.
    """
    available = sorted(
        int(m.group(1))
        for f in glob.glob(os.path.join(kraken_db, "database*mers.kmer_distrib"))
        if (m := re.search(r"database(\d+)mers\.kmer_distrib$", f))
    )
    if not available:
        logger.warning("No Bracken kmer_distrib files in %s – Bracken skipped", kraken_db)
        return None
    if requested:
        if requested not in available:
            raise ValueError(
                f"--bracken-read-length {requested} not built in {kraken_db}; "
                f"available: {available}"
            )
        return requested

    # single-cell R1 is short barcode sequence; Kraken/Bracken see R2
    read_len = first_read_length(samples[0].rev or samples[0].fwd)
    if read_len is None:
        logger.warning("Could not read first FASTQ record; Bracken uses %d", available[0])
        return available[0]
    fitting = [n for n in available if n <= read_len]
    return fitting[-1] if fitting else available[0]


def _write(path: str, content: str) -> None:
    Path(path).write_text(content, encoding="utf-8")
    logger.info("  Written: %s", path)


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------

def run(
    input_dir: str,
    output_dir: str,
    study_name: str | None = None,
    single_cell: bool = False,
    subset_reads: int = 0,
    kraken_db: str | None = None,
    sortmerna_db: str | None = None,
    bracken_read_length: int | None = None,
) -> None:
    """Discover samples, build all QC scripts, and write them to *output_dir*.

    Parameters
    ----------
    input_dir:
        Directory containing raw FASTQ files.
    output_dir:
        Root directory for QC results and generated scripts.
    study_name:
        Label for sub-directories and file prefixes.
        Defaults to ``os.path.basename(input_dir)``.
    single_cell:
        When *True*, contamination / k-mer tools receive only R2 (cDNA
        reads) and insert size is skipped.
    subset_reads:
        When > 0, every tool only processes the first *subset_reads* reads
        per sample.  Useful for rapid pipeline validation.
    kraken_db, sortmerna_db:
        Override database paths from ``ngs_qc/config.py``.
    bracken_read_length:
        Force a Bracken k-mer distribution instead of auto-selecting.
    """
    input_dir  = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    study      = study_name or os.path.basename(input_dir.rstrip("/"))
    study_dir  = f"{output_dir}/{study}"

    kraken_db    = kraken_db    or KRAKEN2_DB
    sortmerna_db = sortmerna_db or SORTMERNA_DB

    samples = discover_samples(input_dir)
    layout  = describe_layout(samples)
    bracken_len = _pick_bracken_length(kraken_db, samples, bracken_read_length)

    logger.info("Study         : %s", study)
    logger.info("Input         : %s", input_dir)
    logger.info("Output        : %s", output_dir)
    logger.info("Samples       : %d (%s-end)", len(samples), layout)
    logger.info("Single-cell   : %s", single_cell)
    logger.info("Kraken2 DB    : %s", kraken_db)
    logger.info("Bracken       : %s", f"{bracken_len}-mer" if bracken_len else "skipped")
    logger.info("SortMeRNA DB  : %s", sortmerna_db)
    if subset_reads:
        logger.info("Subset mode   : first %d reads per sample", subset_reads)

    ctx = Context(
        study_dir=study_dir,
        single_cell=single_cell,
        subset_reads=subset_reads,
        kraken_db=kraken_db,
        sortmerna_db=sortmerna_db,
        bracken_read_length=bracken_len,
        threads=SWARM_THREADS,
        kraken_threads=KRAKEN_SWARM_THREADS,
    )

    # --- build every sample × tool job, grouped by tool ---
    header = [
        "#!/bin/bash",
        f"# NGS-QC swarm script  |  study: {study}  |  layout: {layout}",
        f"# Generated by run_ngs_qc.py v{__version__}",
        f"# Single-cell mode: {single_cell}",
        f"# Subset mode: {f'{subset_reads} reads' if subset_reads else 'off'}",
        "",
    ]
    swarm_lines: dict[str, list[str]] = {KRAKEN_SWARM: list(header), QC_SWARM: list(header)}
    all_jobs: list[Job] = []

    for builder in PER_SAMPLE_BUILDERS:
        jobs = [job for sample in samples for job in builder.build(sample, ctx)]
        if not jobs:
            logger.info("  %-12s skipped (not applicable)", builder.TOOL)
            continue
        lines = [j.cmd for j in jobs if j.cmd]
        swarm = jobs[0].swarm
        swarm_lines[swarm] += [f"# --- {builder.TOOL} ---", *lines, ""]
        all_jobs += jobs
        logger.info("  %-12s %d job(s)", builder.TOOL, len(lines))

    # --- output file paths ---
    os.makedirs(study_dir, exist_ok=True)
    swarm_files = {
        KRAKEN_SWARM: os.path.join(output_dir, f"{study}_kraken.swarm"),
        QC_SWARM:     os.path.join(output_dir, f"{study}_qc.swarm"),
    }
    multiqc_file    = os.path.join(output_dir, f"{study}_multiqc.sh")
    multiqc_config  = os.path.join(study_dir, "multiqc_config.yaml")
    log_report_file = os.path.join(output_dir, f"{study}_log_report.sh")
    manifest_file   = os.path.join(study_dir, "qc_manifest.json")
    exec_file       = os.path.join(output_dir, f"{study}_execution.sh")
    collect_py      = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "collect_qc_logs.py")

    manifest = {
        "study": study,
        "version": __version__,
        "single_cell": single_cell,
        "subset_reads": subset_reads,
        "tools": list(dict.fromkeys(j.tool for j in all_jobs)),
        "samples": [asdict(s) for s in samples],
        "jobs": [
            {"tool": j.tool, "sample": j.sample, "log": j.log, "outputs": j.outputs}
            for j in all_jobs
        ],
    }

    for swarm, path in swarm_files.items():
        _write(path, "\n".join(swarm_lines[swarm]))
    _write(multiqc_file, build_multiqc_script(study, study_dir, multiqc_config))
    _write(multiqc_config, build_multiqc_config(study))
    _write(log_report_file,
           build_log_report_script(study, output_dir, study_dir, collect_py))
    _write(manifest_file, json.dumps(manifest, indent=2) + "\n")
    _write(exec_file, _build_execution_script(
        study, swarm_files, multiqc_file, log_report_file, study_dir))

    logger.info(
        "Done. To submit, run on the cluster:\n\n"
        "    bash %s\n", exec_file
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="run_ngs_qc.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input_dir",
        metavar="INPUT_DIR",
        help="Directory containing raw FASTQ files.",
    )
    parser.add_argument(
        "output_dir",
        metavar="OUTPUT_DIR",
        help="Root directory for QC results and generated scripts.",
    )
    parser.add_argument(
        "--study-name",
        metavar="NAME",
        default=None,
        help=(
            "Study label used in sub-directory names and file prefixes. "
            "Defaults to the basename of INPUT_DIR."
        ),
    )
    parser.add_argument(
        "--single-cell",
        action="store_true",
        default=False,
        help=(
            "Single-cell mode: contamination / k-mer tools (Kraken2, FastQ "
            "Screen, KAT) process only R2 (cDNA reads) and insert size "
            "is skipped. Other tools still run on all reads."
        ),
    )
    parser.add_argument(
        "--kraken-db",
        metavar="PATH",
        default=None,
        help=f"Kraken2 / Bracken database directory. Default: {KRAKEN2_DB}",
    )
    parser.add_argument(
        "--sortmerna-db",
        metavar="PATH",
        default=None,
        help=f"SortMeRNA rRNA reference database directory. Default: {SORTMERNA_DB}",
    )
    parser.add_argument(
        "--bracken-read-length",
        metavar="N",
        type=int,
        default=None,
        help=(
            "Bracken k-mer distribution read length. Defaults to the longest "
            "one in the Kraken2 DB not exceeding the observed read length."
        ),
    )
    parser.add_argument(
        "--subset",
        nargs="?",
        const=10_000,
        default=0,
        type=int,
        metavar="N",
        help=(
            "Subset each sample to the first N reads before running any QC tool. "
            "Useful for rapid pipeline validation. "
            "Using --subset without a value defaults to 10 000 reads."
        ),
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging verbosity (default: INFO).",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    logging.getLogger().setLevel(getattr(logging, args.log_level))

    logger.info("NGS-QC pipeline v%s", __version__)

    try:
        run(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            study_name=args.study_name,
            single_cell=args.single_cell,
            subset_reads=args.subset,
            kraken_db=args.kraken_db,
            sortmerna_db=args.sortmerna_db,
            bracken_read_length=args.bracken_read_length,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("Fatal: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
