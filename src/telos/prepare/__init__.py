"""Input normalization helpers (assembler GTF → Telos-ready GTF / tmap)."""

from telos.prepare.gtf import (
    PrepareGtfError,
    PrepareGtfResult,
    inject_cov_from_tpm,
    load_abundance_table,
    normalize_transcript_cov_attrs,
    run_gffcompare_tmap,
)

__all__ = [
    "PrepareGtfError",
    "PrepareGtfResult",
    "inject_cov_from_tpm",
    "load_abundance_table",
    "normalize_transcript_cov_attrs",
    "run_gffcompare_tmap",
]
