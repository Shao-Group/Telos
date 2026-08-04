"""CLI handler for ``telos prepare-gtf``."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from telos.prepare.gtf import PrepareGtfError, prepare_assembly_gtf

logger = logging.getLogger(__name__)


@dataclass(frozen=True, kw_only=True)
class PrepareGtfIO:
    """Inputs for :func:`run_prepare_gtf`."""

    gtf: Path
    outdir: Path
    tpm: Path | None = None
    ref_gtf: Path | None = None
    make_tmap: bool = False
    gffcompare_bin: Path | None = None


def run_prepare_gtf(cfg: PrepareGtfIO) -> int:
    """
    Normalize an assembly GTF (optional TPM→``cov``) and optionally build a ``.tmap``.

    Returns:
        ``0`` on success, ``2`` on user/input error.
    """
    try:
        result = prepare_assembly_gtf(
            gtf=cfg.gtf,
            outdir=cfg.outdir,
            tpm=cfg.tpm,
            ref_gtf=cfg.ref_gtf,
            make_tmap=cfg.make_tmap,
            gffcompare_bin=cfg.gffcompare_bin,
        )
    except PrepareGtfError as exc:
        logger.error("%s", exc)
        return 2

    print("[telos] prepare-gtf complete")
    print(f"  gtf={result.out_gtf}")
    if result.out_tmap is not None:
        print(f"  tmap={result.out_tmap}")
    print(
        f"  transcripts={result.transcript_lines} "
        f"updated_cov={result.updated_cov} "
        f"missing_abundance={result.missing_abundance}"
    )
    return 0
