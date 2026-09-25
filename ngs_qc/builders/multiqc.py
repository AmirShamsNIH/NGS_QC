"""
MultiQC aggregation step, submitted as a standalone ``sbatch`` job after
both swarms finish.

Report customisation is done through a generated ``multiqc_config.yaml``
rather than by editing the HTML afterwards:

* the ``fastqc_sequence_counts`` section is removed (its count method is
  inaccurate for deduplicated / subsampled libraries);
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

from ..config import MODULES
from .base import q

logger = logging.getLogger(__name__)


def build_multiqc_config(study: str) -> str:
    """Return the MultiQC config (JSON is valid YAML) for *study*."""
    config = {
        "title": f"NGS-QC: {study}",
        "remove_sections": ["fastqc_sequence_counts"],
        "fn_ignore_files": ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"],
        "extra_fn_clean_exts": [
            ".kraken",
            ".bracken_report",
            ".bbduk_stats",
            ".kat",
            {"type": "remove", "pattern": "_subsample"},
            {"type": "remove", "pattern": "_subset"},
        ],
        "table_sample_merge": {
            "R1": [{"type": "regex", "pattern": "_R1(_\\d{3})?$"}],
            "R2": [{"type": "regex", "pattern": "_R2(_\\d{3})?$"}],
        },
        "module_order": [
            {"kraken": {"name": "Kraken2", "anchor": "kraken2",
                        "path_filters": ["kraken/*"]}},
            {"kraken": {"name": "Bracken", "anchor": "bracken",
                        "path_filters": ["bracken/*_report.txt"]}},
            "fastqc", "fastp", "seqkit", "afterqc", "bbduk", "bbmap",
            "fastq_screen", "sortmerna", "kat",
        ],
    }
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
