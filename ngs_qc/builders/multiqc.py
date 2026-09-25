"""
MultiQC aggregation step, submitted as a standalone ``sbatch`` job after
both swarms finish.

Report customisation lives in ``ngs_qc/multiqc_config.json`` (edit it to
change the report) and is copied into each study as ``multiqc_config.yaml``:

* sections that repeat another tool are removed: FastQC is the primary
  per-read QC, so fastp's quality / GC / N / duplication /
  overrepresented-sequence plots go, as does its filtering chart (always
  100 % in report-only mode); ``fastqc_sequence_counts`` is removed too
  (inaccurate for deduplicated / subsampled libraries);
* General Statistics shows each metric once: reads from SeqKit, GC / dups
  / length from FastQC, Q30 from fastp; the duplicate columns are hidden
  (still available via "Configure columns");
* the Kraken module runs twice so Kraken2 and Bracken reports appear as
  separate sections instead of overwriting each other's samples;
* FASTQ files (subset temporaries) are never scanned;
* sample names are normalised so every tool reports one row per sample:
  tool suffixes (``.kraken``, ``.bracken_report``, ``.bbduk_stats``,
  ``.kat``, ``_subsample``, ``_subset``) are stripped and R1 / R2 rows are
  merged under the sample (``table_sample_merge``);
* the data folder of a previous report is ignored: MultiQC >= 1.25
  re-imports ``multiqc.parquet`` from any folder it scans, which would
  bring back the old run's samples and names.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from ..config import MODULES
from .base import q

logger = logging.getLogger(__name__)

# Editable report customisation; see the README's MultiQC section.
CONFIG_TEMPLATE = Path(__file__).resolve().parent.parent / "multiqc_config.json"


def build_multiqc_config(study: str) -> str:
    """Return the MultiQC config for *study* from ``ngs_qc/multiqc_config.json``.

    ``{study}`` in the title is filled in; everything else is copied as is.
    JSON is valid YAML, so the result is written as ``multiqc_config.yaml``.
    """
    config = json.loads(CONFIG_TEMPLATE.read_text(encoding="utf-8"))
    config["title"] = config.get("title", "").format(study=study)
    return json.dumps(config, indent=2) + "\n"


def build_multiqc_script(study: str, study_dir: str, config_path: str) -> str:
    """Return a standalone bash script that runs MultiQC for *study*."""
    report_stem = f"{study}_multiqc_report"
    script_lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        "###############################################",
        "#MULTIQC",
        "",
        f"module load {MODULES['multiqc']}",
        f"cd {q(study_dir)}",
        f"multiqc --force --config {q(config_path)} "
        f"--ignore '*_multiqc_report_data' --filename {report_stem} .",
        "",
    ]
    logger.debug("Built MultiQC script for study '%s'", study)
    return "\n".join(script_lines)
