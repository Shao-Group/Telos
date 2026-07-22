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

    Raises :class:`TmapLabelError` if ``qry_id`` / ``transcript_id`` is not unique (a well-formed
    gffcompare tmap has one row per query transcript).
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
    dup_mask = out["transcript_id"].duplicated(keep=False)
    if dup_mask.any():
        n_dup_rows = int(dup_mask.sum())
        n_dup_ids = int(out.loc[dup_mask, "transcript_id"].nunique())
        examples = (
            out.loc[dup_mask, "transcript_id"]
            .drop_duplicates()
            .head(5)
            .tolist()
        )
        raise TmapLabelError(
            f"tmap has non-unique transcript_id ({n_dup_ids} ids, {n_dup_rows} rows); "
            f"examples: {examples}: {tmap_path}"
        )
    out["label"] = (out["class_code"] == "=").astype(int)
    return out[["transcript_id", "label", "ref_id", "class_code"]]


def load_tmap_labels(tmap_path: Path) -> pd.DataFrame:
    """
    Read a gffcompare ``.tmap`` and produce ``transcript_id`` + binary ``label``.

    Sets ``label=1`` when ``class_code == '='`` (exact match to reference transcript), ``0``
    otherwise. Non-unique query ids raise :class:`TmapLabelError` (see
    :func:`load_tmap_labels_with_ref`).
    """
    df = load_tmap_labels_with_ref(tmap_path)
    return df[["transcript_id", "label"]]
