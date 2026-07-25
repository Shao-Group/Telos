"""
Apply saved Stage I bundles to a feature dataframe and write scored site tables.

Merges per-backend probabilities into ``p_site_rf`` / ``p_site_xgb`` columns aligned to each candidate row.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from telos.models import (
    SITE_PROB_COLUMN_RF,
    SITE_PROB_COLUMN_XGB,
    STAGE1_BACKEND_RF,
    STAGE1_BACKENDS,
    stage1_bundle_path,
)
from telos.models.stage1_train import build_stage1_feature_matrix, load_stage1_bundle, stage1_positive_probas


def _add_site_id(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``site_id`` string key ``chrom:position:strand:site_type`` for stable row identity."""
    out = df.copy()
    out["site_id"] = (
        out["chrom"].astype(str)
        + ":"
        + out["position"].astype(str)
        + ":"
        + out["strand"].astype(str)
        + ":"
        + out["site_type"].astype(str)
    )
    return out


_MERGE_KEYS = ["transcript_id", "site_type", "chrom", "position", "strand"]


def _score_stage1_backend(df_features: pd.DataFrame, model_dir: Path, backend: str) -> pd.DataFrame:
    """
    For one ``backend``, load TSS and TES bundles from ``model_dir``, predict proba on matching rows.

    Returns a skinny dataframe with merge keys + one probability column for that backend.
    """
    col = SITE_PROB_COLUMN_RF if backend == STAGE1_BACKEND_RF else SITE_PROB_COLUMN_XGB
    scored_parts: list[pd.DataFrame] = []
    for site_type in ("TSS", "TES"):
        path = model_dir / stage1_bundle_path(site_type, backend)
        if not path.is_file():
            raise FileNotFoundError(f"Stage I model not found: {path}")
        sub = df_features[df_features["site_type"].str.upper() == site_type].copy()
        if sub.empty:
            continue
        bundle = load_stage1_bundle(path)
        feats: list[str] = list(bundle["feature_names"])
        model = bundle["model"]
        X = build_stage1_feature_matrix(sub, feats)
        sub = sub.reset_index(drop=True)
        sub[col] = stage1_positive_probas(model, X)
        scored_parts.append(sub[_MERGE_KEYS + [col]])

    if not scored_parts:
        return pd.DataFrame(columns=_MERGE_KEYS + [col])

    return pd.concat(scored_parts, axis=0, ignore_index=True)


def score_stage1_dataframe(df_features: pd.DataFrame, model_dir: Path) -> pd.DataFrame:
    """
    Apply saved TSS/TES Stage I models (RF + XGB) to feature rows (must include ``site_type``).

    Returns columns ``site_id``, ``transcript_id``, ``site_type``, ``chrom``, ``position``, ``strand``,
    ``p_site_rf``, ``p_site_xgb``. Raises if any row is missing a backend probability.
    """
    base_cols = [
        "site_id",
        "transcript_id",
        "site_type",
        "chrom",
        "position",
        "strand",
        SITE_PROB_COLUMN_RF,
        SITE_PROB_COLUMN_XGB,
    ]
    if df_features.empty:
        return pd.DataFrame(columns=base_cols)

    out = df_features.copy().reset_index(drop=True)
    for backend in STAGE1_BACKENDS:
        col = SITE_PROB_COLUMN_RF if backend == STAGE1_BACKEND_RF else SITE_PROB_COLUMN_XGB
        side = _score_stage1_backend(df_features, model_dir, backend)
        out = out.merge(side, on=_MERGE_KEYS, how="left")
        n_missing = int(out[col].isna().sum())
        if n_missing:
            raise ValueError(
                f"Stage I scoring incomplete for backend={backend!r}: "
                f"{n_missing} row(s) missing {col} after merge on {_MERGE_KEYS}."
            )

    out = _add_site_id(out)
    return out[
        [
            "site_id",
            "transcript_id",
            "site_type",
            "chrom",
            "position",
            "strand",
            SITE_PROB_COLUMN_RF,
            SITE_PROB_COLUMN_XGB,
        ]
    ].copy()


def write_sites_scored_tsv(df: pd.DataFrame, out_tsv: Path) -> None:
    """Write ``df`` as tab-separated UTF-8 without index (standard ``sites.scored.tsv`` artifact)."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
