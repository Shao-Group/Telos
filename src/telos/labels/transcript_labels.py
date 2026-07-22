"""
Transcript-level labels from gffcompare **.tmap** files for Stage II training.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

class TmapLabelError(RuntimeError):
    """Raised for missing inputs, or invalid GFFcompare operations in this module."""

def load_tmap_labels_with_ref(tmap_path: Path) -> pd.DataFrame:
    """
    Read gffcompare ``.tmap`` with explicit reference id columns.

    Returns columns ``transcript_id``, ``label``, ``ref_id``, ``class_code`` where:
    - ``transcript_id`` = ``qry_id``
    - ``label`` = 1 when ``class_code == '='`` else 0
    """
    df = pd.read_csv(tmap_path, sep="\t", comment="#", header=0)
    need = {"qry_id", "class_code"}
    if not need.issubset(df.columns):
        raise TmapLabelError(f"tmap missing required columns {sorted(need)}: {tmap_path}")
    ref_col = "ref_id" 
    if ref_col not in df.columns:
        raise TmapLabelError(f"tmap missing required column {ref_col}: {tmap_path}")
    out = df[["qry_id", "class_code", ref_col]].rename(
        columns={"qry_id": "transcript_id"}
    )
    out["label"] = (out["class_code"] == "=").astype(int)
    return out[["transcript_id", "label", "ref_id", "class_code"]]


def load_tmap_labels(tmap_path: Path) -> pd.DataFrame:
    """
    Read a gffcompare ``.tmap`` and produce ``transcript_id`` + binary ``label``.

    Keeps columns ``qry_id`` (renamed to ``transcript_id``) and ``class_code``. Sets ``label=1`` when
    ``class_code == '='`` (exact match to reference transcript), ``0`` otherwise. Comment lines and
    standard CSV header are handled by :func:`pandas.read_csv`.
    """
    df = load_tmap_labels_with_ref(tmap_path)
    return df[["transcript_id", "label"]]
