"""
FASTQ file discovery and read-pair resolution.

Files are grouped into samples by their Illumina-style read tag
(``_R1`` / ``_R2``, optionally followed by a ``_NNN`` chunk number) which
must sit immediately before the FASTQ extension.  A library split into
several chunks (``_R1_001``, ``_R1_002`` ...) is QC'd one chunk per sample,
named ``<stem>_<chunk>``; a single chunk keeps the plain stem.  Layout is
resolved per sample: a sample whose R1 has a matching R2 on disk is paired-end,
everything else is single-end, so mixed directories are handled correctly.

Index reads (``_I1`` / ``_I2``) are skipped.  ``Undetermined_*`` reads are
QC'd as their own sample: their volume and content are the main evidence
of demultiplexing / sample-sheet problems.
"""

from __future__ import annotations

import gzip
import logging
import os
import re
from dataclasses import dataclass

logger = logging.getLogger(__name__)

_READS_PER_RECORD: int = 4  # FASTQ: ID + seq + '+' + qual

_FASTQ_EXT_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$")
# <stem>_R1[_001].fastq.gz ; the tag must be directly before the extension so
# that sample names such as "Lib_R10_x" are never mistaken for a read tag.
_READ_TAG_RE = re.compile(
    r"^(?P<stem>.+?)_R(?P<read>[12])(?P<chunk>_\d{3})?(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"
)
_SKIP_RE = re.compile(r"_I[12](?:_\d{3})?\.(?:fastq|fq)(?:\.gz)?$")


@dataclass(frozen=True)
class Sample:
    """One sequencing library: a single-end file or an R1/R2 pair."""

    name: str
    fwd: str
    rev: str | None = None

    @property
    def paired(self) -> bool:
        return self.rev is not None

    @property
    def files(self) -> list[str]:
        return [self.fwd, self.rev] if self.rev else [self.fwd]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def discover_samples(data_dir: str) -> list[Sample]:
    """Scan *data_dir* and return all samples, sorted by name.

    Raises
    ------
    FileNotFoundError
        If no usable FASTQ files are found.
    ValueError
        If two files resolve to the same sample name.
    """
    names = sorted(
        f for f in os.listdir(data_dir)
        if _FASTQ_EXT_RE.search(f) and os.path.isfile(os.path.join(data_dir, f))
    )

    skipped = [n for n in names if _SKIP_RE.search(n)]
    if skipped:
        logger.info("Skipping %d index-read file(s): %s",
                    len(skipped), ", ".join(skipped))
    names = [n for n in names if n not in skipped]

    if not names:
        raise FileNotFoundError(f"No FASTQ files found in: {data_dir}")

    present = set(names)

    # read-1 chunks per stem, to decide whether the chunk goes in the name
    chunks: dict[str, set[str]] = {}
    for name in names:
        if m := _READ_TAG_RE.match(name):
            chunks.setdefault(m.group("stem"), set()).add(m.group("chunk") or "")
    samples: dict[str, Sample] = {}

    def _add(sample: Sample) -> None:
        if sample.name in samples:
            raise ValueError(
                f"Duplicate sample name {sample.name!r}: "
                f"{samples[sample.name].fwd} and {sample.fwd}"
            )
        samples[sample.name] = sample

    for name in names:
        path = os.path.join(data_dir, name)
        m = _READ_TAG_RE.match(name)
        if not m:
            _add(Sample(_FASTQ_EXT_RE.sub("", name), path))
            continue

        stem, read, chunk, ext = m.group("stem", "read", "chunk", "ext")
        mate = f"{stem}_R{'2' if read == '1' else '1'}{chunk or ''}{ext}"
        if len(chunks[stem]) > 1:
            stem += chunk or ""

        if read == "2":
            if mate not in present:
                logger.warning("R2 without matching R1, treating as single-end: %s", name)
                _add(Sample(stem, path))
            continue

        rev = os.path.join(data_dir, mate) if mate in present else None
        _add(Sample(stem, path, rev))

    result = sorted(samples.values(), key=lambda s: s.name)
    logger.debug("Discovered %d sample(s) in %s", len(result), data_dir)
    return result


def describe_layout(samples: list[Sample]) -> str:
    """Return ``"paired"``, ``"single"`` or ``"mixed"`` for logging."""
    kinds = {s.paired for s in samples}
    if kinds == {True}:
        return "paired"
    if kinds == {False}:
        return "single"
    return "mixed"


def fastq_basename(path: str) -> str:
    """File name with the FASTQ extension stripped (FastQC / FastQ Screen naming)."""
    return _FASTQ_EXT_RE.sub("", os.path.basename(path))


def is_gzipped(path: str) -> bool:
    return path.endswith(".gz")


def first_read_length(path: str) -> int | None:
    """Length of the first read in *path*, or None if it cannot be read."""
    opener = gzip.open if is_gzipped(path) else open
    try:
        with opener(path, "rt") as fh:
            fh.readline()
            return len(fh.readline().strip()) or None
    except OSError:
        return None


def subset_cmd(src: str, dest: str, n_reads: int) -> str:
    """Bash command writing the first *n_reads* records of *src* to *dest* (gz).

    ``zcat -f`` passes uncompressed input through unchanged.
    """
    n_lines = n_reads * _READS_PER_RECORD
    return f"zcat -f {src} | head -n {n_lines} | gzip -c > {dest}"
