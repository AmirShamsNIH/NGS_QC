# NGS-QC v2

A modular Python pipeline that generates HPC **swarm** and **sbatch** scripts
for comprehensive, **read-only** quality control of NGS FASTQ datasets on
Slurm-based clusters (e.g. NIH Biowulf).

> **No reads are modified.**  Every tool runs in observation / reporting mode
> only.  fastp disables trimming and filtering, BBDuk and BBMerge only write
> statistics, and AfterQC runs with `--qc_only`.

---

## Table of Contents

1. [Overview](#overview)
2. [QC Tools](#qc-tools)
3. [Project Structure](#project-structure)
4. [Usage](#usage)
5. [Sample Discovery](#sample-discovery)
6. [Single-Cell Mode](#single-cell-mode)
7. [Subset Mode](#subset-mode)
8. [Generated Output Files](#generated-output-files)
9. [Execution on HPC (Swarm)](#execution-on-hpc-swarm)
10. [Configuration](#configuration)
11. [Author](#author)

---

## Overview

`run_ngs_qc.py` takes a FASTQ directory and an output directory and writes
everything needed to run QC on the cluster with one command:

```bash
python run_ngs_qc.py /path/to/fastq /path/to/output
bash /path/to/output/<study>_execution.sh
```

Supports **bulk** and **single-cell** data, paired-end, single-end, and
directories that mix both.

---

## QC Tools

| Tool | Module | Purpose | Swarm |
|------|--------|---------|-------|
| [Kraken2](https://ccb.jhu.edu/software/kraken2/) | kraken/2.1.2 | Taxonomic contamination screening | kraken |
| [Bracken](https://ccb.jhu.edu/software/bracken/) | bracken/2.8 | Species-level abundance re-estimation, chained after Kraken2 | kraken |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | fastqc/0.12.1 | Per-base / per-read quality metrics | qc |
| [fastp](https://github.com/OpenGene/fastp) | fastp/0.23.2 | Quality statistics — **report only** | qc |
| [SeqKit stats](https://bioinf.shenwei.me/seqkit/) | seqkit/2.13.0 | Read counts, lengths, GC, Q20/Q30 | qc |
| [AfterQC](https://github.com/OpenGene/AfterQC) | afterqc/0.9.7 | Quality / bias / overlap report (`--qc_only`) | qc |
| [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | bbtools/39.06 | Adapter and PhiX content — **report only** | qc |
| Insert size ([BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)) | bbtools/39.06 | Overlap-based insert-size histogram (paired bulk only) | qc |
| [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) | fastq_screen/0.15.3 | Multi-genome contamination check | qc |
| [SortMeRNA](https://github.com/sortmerna/sortmerna) | sortmerna/4.3.6 | rRNA fraction (on a 5 % sub-sample) | qc |
| [KAT](https://kat.readthedocs.io/) | kat/2.4.2 | K-mer histogram — contamination, GC bias, copy-number artefacts | qc |
| [MultiQC](https://multiqc.info/) | multiqc/1.34 | Aggregated report | sbatch |

Notes:

* **Bracken** uses the `databaseNNNmers.kmer_distrib` files in the Kraken2 DB.
  The read length is chosen automatically (longest distribution not
  exceeding the observed read length); override with `--bracken-read-length`.
  Bracken is skipped if the DB has no distribution files.
* **Insert size** is measured from read-pair overlap, so pairs whose insert
  is longer than the two reads combined are not counted.  It is skipped for
  single-end samples, in `--single-cell` mode, and when either mate is
  shorter than `INSERT_SIZE_MIN_READ_LEN` (50 bp, e.g. a 28 bp R1), since
  such mates almost never overlap.
* **KAT** distribution analysis is shown in MultiQC; histograms and PNGs
  are in `kat/`.

---

## Project Structure

```
NGS_QC_v2/
├── run_ngs_qc.py               # Main entry point (CLI)
├── collect_qc_logs.py          # Post-run log collector → multi-tab XLSX report
│
└── ngs_qc/
    ├── __init__.py
    ├── config.py               # Module versions, DB paths, swarm resources
    ├── fastq.py                # Sample discovery, pairing, subset helper
    └── builders/
        ├── __init__.py         # PER_SAMPLE_BUILDERS: tool order
        ├── base.py             # Context / Job types, shell helpers
        ├── kraken.py           # Kraken2 + Bracken
        ├── fastqc.py
        ├── fastp.py
        ├── seqkit.py
        ├── afterqc.py
        ├── bbduk.py
        ├── insert_size.py      # BBMerge ihist
        ├── fastq_screen.py
        ├── sortmerna.py
        ├── kat.py
        ├── multiqc.py          # MultiQC sbatch script + multiqc_config.yaml
        └── log_report.py       # Log collection sbatch script
```

Every tool's output is parsed by MultiQC; tools MultiQC cannot read are
not included.  To add a tool: create a module in `builders/` with `TOOL = "<name>"` and
`build(sample, ctx) -> list[Job]`, add it to `PER_SAMPLE_BUILDERS`, and add
its module to `MODULES` in `config.py`.  The XLSX report picks it up
automatically from the manifest.

**Python ≥ 3.10.**  Script generation uses only the standard library; the
log report needs `openpyxl` (included in Biowulf's `python` module).  All
tools are loaded as HPC modules inside the generated scripts.

---

## Usage

```bash
# Minimal
python run_ngs_qc.py /path/to/fastq /path/to/output

# Custom study name
python run_ngs_qc.py /path/to/fastq /path/to/output --study-name MyStudy

# Single-cell data
python run_ngs_qc.py /path/to/fastq /path/to/output --single-cell

# Quick validation on the first 10 000 reads (or --subset N)
python run_ngs_qc.py /path/to/fastq /path/to/output --subset

# Other options
python run_ngs_qc.py --help
```

---

## Sample Discovery

The input directory is scanned for `*.fastq`, `*.fq` (optionally `.gz`).

* A read tag `_R1` / `_R2`, optionally followed by a chunk number
  (`_R1_001`), directly before the extension marks a mate.  The sample name
  is everything before it: `LIB_1_S1_R1_001.fastq.gz` → `LIB_1_S1`.
* A sample whose R1 has a matching R2 is paired-end; files without a read
  tag, or R1/R2 without a mate, are single-end.  Layout is per sample, so
  mixed directories work.
* Index reads (`_I1`, `_I2`) are skipped.  `Undetermined_*` reads are QC'd
  as their own sample (e.g. `Undetermined_L001`), since their size and
  content reveal demultiplexing problems.
* Two files resolving to the same sample name is an error.
* Uncompressed FASTQ is supported.

---

## Single-Cell Mode

Single-cell libraries (e.g. 10x Chromium) put the cell barcode + UMI in R1
and the cDNA in R2.  Barcodes do not map to any genome, so with
`--single-cell`:

| Tool | Bulk (default) | `--single-cell` |
|------|---------------|-----------------|
| FastQC, fastp, SeqKit, AfterQC, BBDuk, SortMeRNA | R1 + R2 | R1 + R2 |
| Kraken2 (+ Bracken) | R1 + R2 paired | R2 only |
| FastQ Screen | R1 + R2 | R2 only |
| KAT | R1 + R2 | R2 only |
| Insert size | R1 + R2 | skipped |

---

## Subset Mode

`--subset` (10 000 reads) or `--subset N` restricts every tool to the first
N reads per sample.  Each tool writes its own temporary subset file in its
output directory and deletes it when it finishes, so parallel jobs never
share a temporary file.  For SortMeRNA the fixed-N subset replaces the
default 5 % BBTools sub-sample.

---

## Generated Output Files

```
<OUTPUT_DIR>/
├── <study>_kraken.swarm           ← Kraken2 + Bracken (large-memory swarm)
├── <study>_qc.swarm               ← all other per-sample tools
├── <study>_multiqc.sh             ← MultiQC sbatch script
├── <study>_log_report.sh          ← XLSX execution report sbatch script
├── <study>_execution.sh           ← submits everything with one command
├── <study>_qc_execution_report.xlsx   (after the run)
└── <study>/
    ├── qc_manifest.json           ← every scheduled job + expected outputs
    ├── multiqc_config.yaml
    ├── <study>_multiqc_report.html    (after the run)
    ├── kraken/  bracken/  fastqc/  fastp/  seqkit/  afterqc/  bbduk/
    ├── insert_size/  fastq_screen/  sortmerna/  kat/
    └── swarm_logs/
```

Each swarm line writes all of its output to one per-sample log
(`<tool dir>/<sample>.<tool>.log`).  Steps within a line are chained with
`&&`, so a failed subset or module load stops the tool and the swarm
subjob exits non-zero.

### XLSX execution report

`collect_qc_logs.py` reads `qc_manifest.json` and checks every scheduled
job.  A job **PASSes** only if all its expected outputs exist and are
non-empty and its log has no fatal error pattern.

| Sheet | Content |
|-------|---------|
| Summary | Sample × tool matrix: PASS (green), FAIL (red), MISSING (amber, never ran), N/A (grey, not applicable) |
| Failed Jobs | Every FAIL / MISSING with log path, missing outputs and fatal lines |
| One sheet per tool | Status, output files and log tail for every sample |

Run it manually at any time:

```bash
python collect_qc_logs.py /path/to/output/MyStudy MyStudy
```

### MultiQC report

Configured by `<study>/multiqc_config.yaml`: the FastQC sequence-count plot
is removed (inaccurate for deduplicated / subsampled libraries), Kraken2
and Bracken are shown as separate sections, and sample names are cleaned
so the General Statistics table has one row per sample with its R1 / R2
rows grouped underneath.  A previous report's `*_multiqc_report_data`
folder is ignored so re-runs never mix in old data.

---

## Execution on HPC (Swarm)

```bash
bash <OUTPUT_DIR>/<study>_execution.sh
```

The wrapper (with `set -euo pipefail`, stopping if a submission fails):

1. submits `<study>_kraken.swarm` (`-g 32 -t 20`) and `<study>_qc.swarm`
   (`-g 48 -t 10`);
2. submits MultiQC and the log report with
   `--dependency afterany:$kraken_jobID:$qc_jobID`, so they run once both
   swarms have finished, whether or not individual jobs failed.

---

## Configuration

Edit `ngs_qc/config.py`:

| Setting | Description |
|---------|-------------|
| `MODULES` | HPC module name/version for each tool |
| `KRAKEN2_DB` | Kraken2 / Bracken database directory |
| `BRACKEN_LEVEL`, `BRACKEN_THRESHOLD` | Bracken taxonomic level and read threshold |
| `SORTMERNA_DB` / `SORTMERNA_REFS` | SortMeRNA reference databases |
| `SORTMERNA_SUBSAMPLE_RATE` | Fraction of reads given to SortMeRNA (default 0.05) |
| `INSERT_SIZE_READS` | Read pairs examined by BBMerge (0 = all) |
| `INSERT_SIZE_MIN_READ_LEN` | Skip insert size when a mate is shorter than this |
| `KRAKEN_SWARM_MEMORY_GB`, `KRAKEN_SWARM_THREADS` | Kraken swarm resources |
| `SWARM_MEMORY_GB`, `SWARM_THREADS`, `SWARM_TIME`, `SWARM_PARTITION` | QC swarm resources |
| `MULTIQC_*` | MultiQC sbatch resources |

---

## Author

**Amir Shams**
amir.shams84@gmail.com
