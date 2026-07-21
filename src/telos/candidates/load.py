"""
Canonical loaders bridging GTF inputs to Stage I/II data structures.

Delegates candidate extraction to :mod:`telos.candidates.extract` and transcript tables to
:mod:`telos.backends.gtfformat`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from telos.backends.gtfformat import build_cov_dataframe
from telos.candidates.extract import CandidateSite, extract_candidate_sites_from_gtf


def load_candidates(gtf_path: Path) -> list[CandidateSite]:
    """
    Return all TSS/TES candidates derived from ``gtf_path`` transcript records.

    Thin wrapper around :func:`~telos.candidates.extract.extract_candidate_sites_from_gtf` — this
    is the function :mod:`telos.pipeline_core` calls so call sites stay stable.
    """
    return extract_candidate_sites_from_gtf(gtf_path)


def load_transcript_cov_dataframe(gtf_path: Path) -> pd.DataFrame:
    """
    Build the per-transcript coverage / exon-statistics dataframe for Stage II merges.

    Uses pure-Python :func:`~telos.backends.gtfformat.build_cov_dataframe` (no external gtfformat
    binary required). Column set matches what :mod:`telos.features.stage2` expects for joins.
    """
    return build_cov_dataframe(gtf_path)
