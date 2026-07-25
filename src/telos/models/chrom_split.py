"""
Chromosome naming and **train/validation split** helpers for autosome ranges.

Maps many contig styles (``chr1``, ``1``, ``NC_000001.11``) to primary autosome indices 1–22 so a
policy like ``chr1-10`` behaves consistently across references.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

class ChromSplitError(RuntimeError):
    """Error raised by chromosome split helpers."""


def primary_autosome_number(chrom: object) -> int | None:
    """
    Map a sequence name to a human-style autosome index 1-22 when unambiguous.

    Recognizes:
    - GENCODE/UCSC: ``chr1`` … ``chr22`` (any ``chr`` casing)
    - Ensembl-style: ``1`` … ``22`` (bare integer tokens only)
    - RefSeq replicons: ``NC_000001.11`` … ``NC_000022.*`` (human nuclear; accession
      number 1-22 after ``NC_``; optional version suffix), including values prefixed
      with ``chr`` (e.g. after :func:`telos.labels.site_labels.normalize_chrom_name`).

    Contigs and scaffolds (``NT_``, ``NW_``, ``NS_``, …, UCSC ``GL*`` / ``KI*`` / ``JH*``,
    ``chrUn_*``) always return ``None``.

    Returns ``None`` for X/Y/MT/scaffolds/alt contigs or unrecognized names so they
    fall into the validation split (not training) for ``chr1-10``-style policies.
    """
    if chrom is None:
        raise ChromSplitError("chrom is None: can not parse primary autosome number")
    s = str(chrom).strip()
    if not s:
        raise ChromSplitError("chrom is empty: can not parse primary autosome number")
    low = s.lower()
    if low.startswith("chr"):
        body = s[3:]
    else:
        body = s
    low_b = body.lower()

    # NCBI non-chromosomal / scaffold accessions (never primary autosomes for split policy).
    if low_b.startswith(("nt_", "nw_", "ns_", "nr_", "ng_", "nz_")):
        return None
    # UCSC hg* patch and random contigs (GL…, KI…, JH…), chrUn_…, etc.
    if low_b.startswith("un_") or re.match(r"^(?:gl|ki|jh)\d", low_b):
        return None

    m = re.fullmatch(r"NC_(0*\d+)(?:\.\d+)?", body, flags=re.IGNORECASE)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 22:
            return n
        return None

    m = re.fullmatch(r"(\d+)", low_b)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 22:
            return n

    return None


def parse_split_policy(policy: str) -> tuple[int, int]:
    """
    Parse policies like ``chr1-10`` or ``1-10`` into an inclusive autosome index range.

    Training rows are those for which :func:`primary_autosome_number` is in this range;
    validation rows are all others (including X/Y/MT and contigs).
    """
    raw = str(policy).strip().lower().replace(" ", "")
    m = re.fullmatch(r"(?:chr)?(\d+)-(?:chr)?(\d+)", raw)
    if not m:
        raise ValueError(
            f"Unsupported split_policy={policy!r}; expected form like 'chr1-10' or '1-10' or 'chr1-chr10'."
            f"Got {policy!r}."
        )
    lo, hi = int(m.group(1)), int(m.group(2))
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi


def split_train_val_masks(
    df: pd.DataFrame,
    autosome_train_range: tuple[int, int],
    *,
    chrom_col: str = "chrom",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Train = rows whose ``chrom`` maps to an autosome index inside ``autosome_train_range``.

    Uses :func:`primary_autosome_number` so RefSeq ``NC_*``, Ensembl ``1``…``22``, and
    ``chr1``…``chr22`` names all align with the same policy without user tuning.
    """
    lo, hi = autosome_train_range
    if lo > hi:
        lo, hi = hi, lo
    nums = df[chrom_col].map(primary_autosome_number)
    in_train = nums.notna() & (nums >= lo) & (nums <= hi)
    train_mask = np.asarray(in_train, dtype=bool)
    val_mask = ~train_mask
    return train_mask, val_mask


def seqnames_on_validation_split_from_gtf(
    gtf_path: Path,
    autosome_train_range: tuple[int, int],
) -> list[str]:
    """
    Unique GTF sequence names (column 1) that fall on the **validation** side of the
    autosome train range — same rule as :func:`split_train_val_masks` (train = primary
    autosome index in range; validation = all other rows' chromosomes).

    Used for transcript PR: build a chromosome list for ``gtfformat filter-chrom`` on assembly and
    reference before ``update-transcript-cov`` / ``gffcompare``.
    """
    lo, hi = autosome_train_range
    if lo > hi:
        lo, hi = hi, lo
    names: set[str] = set()
    with gtf_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 1:
                continue
            seq = cols[0].strip()
            if seq:
                names.add(seq)
    val: list[str] = []
    for seq in names:
        n = primary_autosome_number(seq)
        if n is None or not (lo <= n <= hi):
            val.append(seq)
    return sorted(val)


def write_chrom_split_debug_lists(
    df: pd.DataFrame,
    train_mask: np.ndarray,
    val_mask: np.ndarray,
    reports_dir: Path,
    *,
    chrom_col: str = "chrom",
) -> None:
    """
    Write ``validation_chromosomes.txt`` and ``train_chromosomes.txt`` under ``reports_dir``.

    Intended for debugging splits when training saves intermediates.
    """
    reports_dir.mkdir(parents=True, exist_ok=True)
    tr = df.loc[train_mask]
    va = df.loc[val_mask]
    (reports_dir / "validation_chromosomes.txt").write_text(
        "\n".join(str(x) for x in sorted(va[chrom_col].unique())) + "\n",
        encoding="utf-8",
    )
    (reports_dir / "train_chromosomes.txt").write_text(
        "\n".join(str(x) for x in sorted(tr[chrom_col].unique())) + "\n",
        encoding="utf-8",
    )
