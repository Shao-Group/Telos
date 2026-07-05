"""
Stage II **tabular** merges: cov dataframe + TSS/TES site probabilities + derived numeric features.

Operates on :class:`pandas.DataFrame` objects produced from :mod:`telos.backends.gtfformat` and
Stage I scoring; no BAM I/O here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

# Join keys for coordinate-based cov × site merges (dropped from the result after each join).
STAGE2_SITE_JOIN_KEYS: tuple[str, ...] = ("chrom", "position", "strand")


def suffix_site_table_columns(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
    """Rename non-join site columns with ``suffix`` (e.g. ``probability`` → ``probability_tss``)."""
    join = set(STAGE2_SITE_JOIN_KEYS)
    rename = {col: f"{col}{suffix}" for col in df.columns if col not in join}
    return df.rename(columns=rename)


def read_cov_tsv(path: Path) -> pd.DataFrame:
    """
    Load a tab-separated cov table (e.g. from ``gtfformat get-cov`` or Python :func:`build_cov_dataframe` dump).

    Forces string dtypes for id/chrom columns for stable merges.
    """
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "transcript_id": str,
            "strand": str,
            "tss_chrom": str,
            "tes_chrom": str,
        },
    )


def site_table_for_stage2_merge(
    df_all: pd.DataFrame,
    site_type: str,
    df_scored: pd.DataFrame,
    *,
    prob_col: str = "p_site_rf",
) -> pd.DataFrame:
    """
    Reduce Stage I rows for one ``site_type`` to **one row per genomic site** with a probability column.

    Left-joins feature rows ``df_all`` with scored probabilities from ``df_scored`` on
    ``transcript_id, chrom, position, site_type``. Missing probs become ``0.5``.

    For Stage II coordinate-based merging we must avoid transcript-id deduplication (a transcript_id can
    appear with ``strand="."`` expanded into both ``+`` and ``-`` conventions). Instead we dedupe by
    ``(chrom, position, strand, site_type)`` keeping the max-prob row.

    Returns the full Stage I feature row per site (all columns from ``df_all`` plus ``probability``).
    """
    st = site_type.upper()
    sub = df_all[df_all["site_type"].str.upper() == st].copy()
    if "transcript_id" not in sub.columns:
        raise KeyError("Stage I feature table missing transcript_id")

    if prob_col not in df_scored.columns:
        raise KeyError(f"Stage I scored table missing column {prob_col!r} (expected dual Stage I outputs).")
    prob = df_scored[df_scored["site_type"].str.upper() == st][
        ["transcript_id", "chrom", "position", "site_type", prob_col]
    ].copy()
    prob = prob.rename(columns={prob_col: "probability"})
    out = sub.merge(prob, on=["transcript_id", "chrom", "position", "site_type"], how="left")
    if out["probability"].isna().any():
        print(f"WARNING: {out['probability'].isna().sum()} rows have missing probability from stage 1 scoring")
        out["probability"] = out["probability"].fillna(0.5)
    n_rows_before = int(len(out))
    n_unique_sites_before = int(
        out.drop_duplicates(["chrom", "position", "strand", "site_type"]).shape[0]
    )
    # Keep one representative row per genomic site (coordinate+strand) for stable Stage II merges.
    out = out.sort_values(
        ["probability", "chrom", "position"],
        ascending=[False, True, True],
    ).drop_duplicates(["chrom", "position", "strand", "site_type"], keep="first")
    n_rows_after = int(len(out))
    duplicate_rows_removed = n_rows_before - n_rows_after
    print(
        f"[telos] Stage II site merge ({st}, {prob_col}): "
        f"rows={n_rows_before}, unique_sites={n_unique_sites_before}, "
        f"dedup_rows_removed={duplicate_rows_removed}, "
        f"columns={len(out.columns)}"
    )
    return out.reset_index(drop=True)


def merge_cov_tss_tes(
    df_cov: pd.DataFrame,
    df_tss: pd.DataFrame,
    df_tes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Inner-join cov with full TSS/TES site tables on coordinates and strand.

    All non-join site columns are suffixed ``_tss`` / ``_tes`` **before** merging so every Stage I
    feature from both ends is retained for training while cov (transcript-level) columns stay bare.
    """
    tss = suffix_site_table_columns(
        df_tss.drop(columns=["transcript_id"], errors="ignore").copy(),
        "_tss",
    )
    tes = suffix_site_table_columns(
        df_tes.drop(columns=["transcript_id"], errors="ignore").copy(),
        "_tes",
    )
    df = df_cov.merge(
        tss,
        left_on=["tss_chrom", "tss_pos", "strand"],
        right_on=list(STAGE2_SITE_JOIN_KEYS),
        how="inner",
    )
    df = df.drop(columns=["chrom", "position"], errors="ignore")
    df = df.merge(
        tes,
        left_on=["tes_chrom", "tes_pos", "strand"],
        right_on=list(STAGE2_SITE_JOIN_KEYS),
        how="inner",
    )
    return df.drop(columns=["chrom", "position"], errors="ignore")


def add_stage2_derived_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add length, confidence, exon-shape, log, and interaction columns used by Stage II LightGBM.

    Expects merged cov+TSS+TES columns (e.g. ``probability_tss`` / ``probability_tes``, exon fields).
    Uses :func:`numpy.log1p` for skew transforms. Returns a **copy**; does not mutate input.
    """
    df = df.copy()
    df["transcript_length"] = np.abs(df["tes_pos"] - df["tss_pos"])
    df["log_transcript_length"] = np.log1p(df["transcript_length"])
    df["tss_confidence"] = df.get("probability_tss", 0.5)
    df["tes_confidence"] = df.get("probability_tes", 0.5)
    df["min_confidence"] = np.minimum(df["tss_confidence"], df["tes_confidence"])
    df["confidence_product"] = df["tss_confidence"] * df["tes_confidence"]

    if "exon_count" in df.columns:
        df["exon_density"] = df["exon_count"] / df["transcript_length"].clip(lower=1)
        df["confidence_exon_interaction"] = df["confidence_product"] * df["exon_count"]
        df["coverage_per_exon"] = df["coverage"] / df["exon_count"].clip(lower=1)

    if "total_exon_length" in df.columns and "exon_count" in df.columns:
        df["avg_exon_length"] = df["total_exon_length"] / df["exon_count"].clip(lower=1)

    if "max_exon_length" in df.columns and "min_exon_length" in df.columns:
        df["exon_length_ratio"] = df["max_exon_length"] / df["min_exon_length"].clip(lower=1)

    df["coverage_length_ratio"] = df["coverage"] / df["transcript_length"].clip(lower=1)
    df["confidence_coverage_interaction"] = df["confidence_product"] * df["coverage"]

    df["confidence_sum"] = df["tss_confidence"] + df["tes_confidence"]
    df["confidence_diff"] = abs(df["tss_confidence"] - df["tes_confidence"])

    if "first_exon_length" in df.columns and "last_exon_length" in df.columns:
        df["terminal_exon_ratio"] = df["first_exon_length"] / df["last_exon_length"].clip(lower=1)

    if "mean_exon_length" in df.columns and "transcript_length" in df.columns:
        df["exon_efficiency"] = df["mean_exon_length"] / df["transcript_length"].clip(lower=1)

    skewed_features = ["coverage", "total_reads", "transcript_length"]
    for feature in skewed_features:
        if feature in df.columns:
            df[f"log_{feature}"] = np.log1p(df[feature])

    if all(col in df.columns for col in ["max_exon_length", "min_exon_length", "mean_exon_length"]):
        exon_length_cv = (df["max_exon_length"] - df["min_exon_length"]) / df["mean_exon_length"].clip(
            lower=1
        )
        df["exon_length_entropy"] = np.log1p(exon_length_cv)

    if "exon_count" in df.columns and "total_exon_length" in df.columns:
        df["exon_length_std"] = (
            np.sqrt(df["exon_length_variance"]) if "exon_length_variance" in df.columns else 0
        )
        df["exon_length_skewness"] = (df["max_exon_length"] - df["mean_exon_length"]) / df[
            "exon_length_std"
        ].clip(lower=1)

    return df


STAGE2_DROP_COLS = [
    "chrom",
    "position",
    "chrom_tss",
    "position_tss",
    "chrom_tes",
    "position_tes",
    "tss_chrom",
    "tss_pos",
    "tes_chrom",
    "tes_pos",
    "site_type",
    "site_type_tss",
    "site_type_tes",
    "strand",
    "strand_tss",
    "strand_tes",
    "ref_id",
    "chrom_num",
    "transcript_id",
    "label",
    "label_tss",
    "label_tes",
]


def list_stage2_model_feature_names(df: pd.DataFrame) -> list[str]:
    """
    All ``df`` columns suitable as model inputs: exclude ids, coordinates, strand, and label columns.

    Only drops entries from :data:`STAGE2_DROP_COLS` that actually exist, so optional columns do not break.
    """
    return [c for c in df.columns if c not in STAGE2_DROP_COLS]
