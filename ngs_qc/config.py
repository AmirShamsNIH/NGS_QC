"""
Global configuration constants for the NGS-QC pipeline.

All tool module versions, database paths, and runtime defaults are
defined here so they can be updated in one place.
"""

# ---------------------------------------------------------------------------
# HPC module names  (used in `module load <name>` calls)
# ---------------------------------------------------------------------------
MODULES: dict[str, str] = {
    "kraken2":      "kraken/2.1.2",
    "bracken":      "bracken/2.8",
    "fastqc":       "fastqc/0.12.1",
    "fastq_screen": "fastq_screen/0.15.3",
    "fastp":        "fastp/0.23.2",
    "kat":          "kat/2.4.2",
    "sortmerna":    "sortmerna/4.3.6",
    "bbtools":      "bbtools/39.06",
    "afterqc":      "afterqc/0.9.7",
    "multiqc":      "multiqc/1.34",
    "python":       "python",
}

# ---------------------------------------------------------------------------
# Reference database paths
# ---------------------------------------------------------------------------
# Kraken2 PlusPFP-16 (hash.k2d is 16 GB); Bracken's databaseNNNmers.kmer_distrib
# files live in the same directory.
KRAKEN2_DB: str = "/data/RTB_GRS/references/kraken"

SORTMERNA_DB: str = (
    "/usr/local/apps/sortmeRNA/4.3.6/sortmerna-4.3.6/data/rRNA_databases"
)

SORTMERNA_REFS: list[str] = [
    "silva-euk-28s-id98.fasta",
    "rfam-5.8s-database-id98.fasta",
    "silva-arc-16s-id95.fasta",
    "silva-euk-18s-id95.fasta",
    "rfam-5s-database-id98.fasta",
    "silva-bac-23s-id98.fasta",
    "silva-bac-16s-id90.fasta",
    "silva-arc-23s-id98.fasta",
]

# ---------------------------------------------------------------------------
# Tool parameters
# ---------------------------------------------------------------------------
BRACKEN_LEVEL:     str = "S"   # species
BRACKEN_THRESHOLD: int = 10    # min reads for a taxon to be re-estimated

SORTMERNA_SUBSAMPLE_RATE: float = 0.05

# Read pairs examined by BBMerge for the insert-size histogram (0 = all).
INSERT_SIZE_READS: int = 2_000_000
# Skip insert size when either mate is shorter than this: BBMerge needs the
# mates to overlap, which a short R1 (e.g. 28 bp barcode/UMI read) almost
# never does, leaving a histogram built from a handful of pairs.
INSERT_SIZE_MIN_READ_LEN: int = 50

# ---------------------------------------------------------------------------
# Default swarm / sbatch runtime parameters
# ---------------------------------------------------------------------------
# Kraken2 swarm: the 16 GB DB is loaded into RAM once per sample
KRAKEN_SWARM_MEMORY_GB: int = 32
KRAKEN_SWARM_THREADS:   int = 20
# All other QC tools
SWARM_MEMORY_GB: int  = 48   # SortMeRNA is run with -m 40000 (MB)
SWARM_THREADS:   int  = 10
SWARM_TIME:      str  = "24:00:00"
SWARM_PARTITION: str  = "norm"

MULTIQC_MEMORY:  str  = "10G"
MULTIQC_CPUS:    int  = 10
MULTIQC_LSCRATCH_GB: int = 200
