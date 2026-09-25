"""
ngs_qc – NGS Quality Control pipeline (v2)
===========================================

Generates HPC swarm scripts for read-only, per-sample QC (Kraken2 + Bracken,
FastQC, fastp, SeqKit, AfterQC, BBDuk, insert size, FastQ Screen,
SortMeRNA, KAT) plus MultiQC and XLSX log-report steps.

Typical usage
-------------
::

    python run_ngs_qc.py INPUT_DIR OUTPUT_DIR

See ``run_ngs_qc.py --help`` for the generated files and options.
"""

__version__ = "2.1.0"
__author__  = "Amir Shams"
__email__   = "amir.shams84@gmail.com"
