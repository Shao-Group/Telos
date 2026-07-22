"""
Stage I site labeling: build reference TSS/TES from a reference GTF and label candidate rows by proximity.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from telos.candidates.load import load_candidates

class SiteLabelError(RuntimeError):
    """Raised for invalid site labeling operations or missing inputs."""
    


def normalize_chrom_name(chrom: object) -> str:
    """
    Normalize contig names to ``chr*`` form for joins between BAM, GTF, and coverage tables.

    If the string already starts with ``chr`` (case-insensitive), lowercases that prefix to ``chr``
    and keeps the remainder; otherwise prepends ``chr``.
    """
    s = str(chrom).strip()
    if not s:
        return s
    if len(s) > 3 and s[:3].lower() == "chr":
        return "chr" + s[3:]
    return f"chr{s}"


def reference_sites_from_gtf(ref_gtf: Path) -> pd.DataFrame:
    """
    True TSS/TES positions from reference GTF transcript records
    (same geometry as extract_candidate_sites_from_gtf).
    """
    sites = load_candidates(ref_gtf)
    if not sites:
        return pd.DataFrame(columns=["site_type", "chrom", "position", "strand"])
    return pd.DataFrame(
        {
            "site_type": [s.site_type for s in sites],
            "chrom": [s.chrom for s in sites],
            "position": [s.position for s in sites],
            "strand": [s.strand for s in sites],
        }
    )


def novel_reference_sites_from_gtf(ref_gtf: Path, *, novel_prefix: str = "NOVEL_TX_") -> pd.DataFrame:
    """
    Reference TSS/TES for transcripts tagged as novel in the augmented reference GTF.

    Novel transcripts are identified by ``transcript_id`` starting with ``novel_prefix``.
    """
    out: list[dict[str, object]] = []
    sites = load_candidates(ref_gtf)
    for s in sites:
        if not str(s.transcript_id).startswith(novel_prefix):
            continue
        out.append(
            {
                "site_type": s.site_type,
                "chrom": s.chrom,
                "position": int(s.position),
                "strand": s.strand,
            }
        )
    if not out:
        raise SiteLabelError("No novel reference sites found in GTF: %s" % ref_gtf)
    return pd.DataFrame(out)


def label_sites_by_proximity(
    features_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    site_type: str,
    tolerance_bp: int,
) -> pd.Series:
    """
    Binary label per row: 1 if any reference site of the same site_type matches
    chrom (normalized), strand, and |Δpos| <= tolerance_bp.
    """
    st = site_type.upper()
    feat = features_df.copy()
    ref = ref_df[ref_df["site_type"].str.upper() == st].copy()
    feat["_chrom_n"] = feat["chrom"].map(normalize_chrom_name)
    ref["_chrom_n"] = ref["chrom"].map(normalize_chrom_name)

    labels = pd.Series(0, index=feat.index, dtype=int)
    if ref.empty:
        raise SiteLabelError("No reference sites found for site type: %s" % site_type)

    ref_grouped: dict[tuple[str, str], np.ndarray] = {}
    for (c, strand), g in ref.groupby(["_chrom_n", "strand"]):
        ref_grouped[(c, strand)] = g["position"].astype(int).values

    for (c, strand), g in feat.groupby(["_chrom_n", "strand"], sort=False):
        ref_pos = ref_grouped.get((c, strand))
        if ref_pos is None or len(ref_pos) == 0:
            raise SiteLabelError("No reference sites found for site type: %s on chromosome: %s and strand: %s" % (site_type, c, strand))

        pos = g["position"].astype(int).values
        diff = np.abs(pos[:, None] - ref_pos[None, :])
        hit = (diff <= int(tolerance_bp)).any(axis=1).astype(int)
        labels.loc[g.index] = hit

    return labels
