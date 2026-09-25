"""
sbatch script builder for the post-run log collection and XLSX report.

The generated script is submitted after both swarms finish
(``--dependency afterany:$kraken_jobID:$qc_jobID``).  It calls
``collect_qc_logs.py``, which checks every job listed in the study's
``qc_manifest.json`` and writes a multi-tab XLSX report showing which
samples passed or failed per tool.
"""

from __future__ import annotations

from ..config import MODULES, SWARM_PARTITION
from .base import q


def build_log_report_script(study: str, output_dir: str, study_dir: str,
                            collect_script: str) -> str:
    report_xlsx = f"{output_dir}/{study}_qc_execution_report.xlsx"
    log_out     = f"{output_dir}/{study}_log_report_job.log"

    return "\n".join([
        "#!/bin/bash",
        "#SBATCH --mem=4G",
        "#SBATCH --cpus-per-task=2",
        "#SBATCH --time=00:30:00",
        f"#SBATCH --partition={SWARM_PARTITION}",
        f"#SBATCH --job-name={study}_qc_report",
        f"#SBATCH --output={log_out}",
        "",
        f"module load {MODULES['python']}",
        "",
        f"python {q(collect_script)} {q(study_dir)} {q(study)} --output {q(report_xlsx)}",
        "",
    ])
