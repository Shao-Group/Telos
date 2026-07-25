"""
Build Stage II training/inference frames (cov × TSS × TES), train LightGBM, save model + ranked TSV.

Training merges gffcompare ``tmap`` labels; prediction omits labels and uses the same feature pipeline without the label column.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    classification_report,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline

from telos.features.stage2 import (
    add_stage2_derived_features,
    list_stage2_model_feature_names,
    merge_cov_tss_tes,
    read_cov_tsv,
    site_table_for_stage2_merge,
)
from telos.labels.transcript_labels import load_tmap_labels
from telos.models import (
    stage2_feature_names_json_for_backend,
    stage2_model_joblib_for_backend,
    transcripts_ranked_tsv_for_backend,
)
from telos.models.chrom_split import split_train_val_masks, write_chrom_split_debug_lists


def _log_stage2_row_counts(
    *,
    n_cov: int,
    n_after_sites: int,
    n_after_labels: int | None,
    context: str,
) -> None:
    """Print merge cardinality and warnings when many transcripts drop at cov×site or tmap joins."""
    extra = (
        f" -> after tmap/labels={n_after_labels}"
        if n_after_labels is not None
        else ""
    )
    print(
        f"[telos] Stage II merge ({context}): cov transcripts={n_cov} "
        f"-> after TSS/TES join={n_after_sites}{extra}"
    )
    if n_cov > 20 and n_after_sites < n_cov:
        print(
            f"[telos] warning: {n_cov - n_after_sites} transcript(s) dropped in cov×TSS×TES join "
            "(missing Stage I site rows or transcript_id mismatch)."
        )
    if (
        n_after_labels is not None
        and n_after_sites > 20
        and n_after_labels < n_after_sites
    ):
        print(
            f"[telos] warning: {n_after_sites - n_after_labels} transcript(s) dropped merging tmap "
            "(check .tmap qry_id vs assembly transcript_id)."
        )


def _normalize_stage2_chroms(
    df_cov: pd.DataFrame, df_tss: pd.DataFrame, df_tes: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply :func:`~telos.labels.site_labels.normalize_chrom_name` to cov and site ``chrom`` columns (copies)."""
    from telos.labels.site_labels import normalize_chrom_name

    dc = df_cov.copy()
    dc["tss_chrom"] = dc["tss_chrom"].map(normalize_chrom_name)
    dc["tes_chrom"] = dc["tes_chrom"].map(normalize_chrom_name)
    ts = df_tss.copy()
    te = df_tes.copy()
    ts["chrom"] = ts["chrom"].map(normalize_chrom_name)
    te["chrom"] = te["chrom"].map(normalize_chrom_name)
    return dc, ts, te


def build_stage2_training_frame(
    cov_source: Path | pd.DataFrame,
    df_all_stage1: pd.DataFrame,
    df_sites_scored: pd.DataFrame,
    tmap_path: Path,
    *,
    site_prob_column: str = "p_site_rf",
) -> pd.DataFrame:
    """
    Training path: cov → TSS/TES site probs → inner join tmap labels → derived features.

    ``cov_source`` is a path to cov TSV or an in-memory dataframe. ``site_prob_column`` selects which
    Stage I backend’s site probabilities feed the merge (``p_site_rf`` or ``p_site_xgb``).

    Returns:
        Labeled dataframe with ``label`` and Stage II feature columns; ``chrom`` is taken from
        ``tes_chrom`` (see :func:`~telos.features.stage2.merge_cov_tss_tes`).
    """
    df_cov = read_cov_tsv(cov_source) if isinstance(cov_source, Path) else cov_source.copy()
    n_cov = len(df_cov)
    df_tss = site_table_for_stage2_merge(
        df_all_stage1, "TSS", df_sites_scored, prob_col=site_prob_column
    )
    df_tes = site_table_for_stage2_merge(
        df_all_stage1, "TES", df_sites_scored, prob_col=site_prob_column
    )
    df_cov, df_tss, df_tes = _normalize_stage2_chroms(df_cov, df_tss, df_tes)
    df = merge_cov_tss_tes(df_cov, df_tss, df_tes)
    n_sites = len(df)
    df_label = load_tmap_labels(tmap_path)
    df = df.merge(df_label, on="transcript_id", how="inner")
    _log_stage2_row_counts(
        n_cov=n_cov,
        n_after_sites=n_sites,
        n_after_labels=len(df),
        context="train",
    )
    return add_stage2_derived_features(df)


def build_stage2_inference_frame(
    cov_source: Path | pd.DataFrame,
    df_all_stage1: pd.DataFrame,
    df_sites_scored: pd.DataFrame,
    *,
    site_prob_column: str = "p_site_rf",
) -> pd.DataFrame:
    """
    Predict path: same cov × TSS × TES merge as training but **no** ``tmap`` / ``label`` join.

    Used by :mod:`telos.commands.predict` after Stage I scoring.
    """
    df_cov = read_cov_tsv(cov_source) if isinstance(cov_source, Path) else cov_source.copy()
    n_cov = len(df_cov)
    df_tss = site_table_for_stage2_merge(
        df_all_stage1, "TSS", df_sites_scored, prob_col=site_prob_column
    )
    df_tes = site_table_for_stage2_merge(
        df_all_stage1, "TES", df_sites_scored, prob_col=site_prob_column
    )
    df_cov, df_tss, df_tes = _normalize_stage2_chroms(df_cov, df_tss, df_tes)
    df = merge_cov_tss_tes(df_cov, df_tss, df_tes)
    _log_stage2_row_counts(
        n_cov=n_cov,
        n_after_sites=len(df),
        n_after_labels=None,
        context="predict",
    )
    return add_stage2_derived_features(df)


def build_stage2_lightgbm_pipeline(*, n_jobs: int = -1) -> Pipeline:
    """
    Default :class:`lightgbm.LGBMClassifier` wrapped in a one-step sklearn ``Pipeline``.

    ``n_jobs`` controls tree fitting parallelism (from YAML ``lightgbm.n_jobs`` in training).
    """
    try:
        import lightgbm as lgb  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Stage II training requires lightgbm. Install with: pip install lightgbm"
        ) from exc

    clf_lgb = lgb.LGBMClassifier(
        n_estimators=300,
        max_depth=6,
        learning_rate=0.03,
        num_leaves=15,
        feature_fraction=0.7,
        bagging_fraction=0.7,
        bagging_freq=3,
        reg_alpha=0.2,
        reg_lambda=0.2,
        min_child_samples=20,
        min_split_gain=0.01,
        subsample_for_bin=200000,
        random_state=42,
        n_jobs=int(n_jobs),
        verbose=-1,
        force_col_wise=True,
    )
    return Pipeline([("clf", clf_lgb)])


def stage2_proba_positive_binary(clf: Pipeline, X: pd.DataFrame) -> np.ndarray:
    """Column1 of ``predict_proba`` when binary; constant vector if only one class was trained."""
    proba = clf.predict_proba(X)
    if proba.shape[1] >= 2:
        return proba[:, 1]
    only = int(clf.classes_[0])
    return np.full(proba.shape[0], 1.0 if only == 1 else 0.0)


def train_and_save_stage2(
    df: pd.DataFrame,
    models_dir: Path,
    predictions_dir: Path,
    reports_dir: Path | None,
    *,
    autosome_train_range: tuple[int, int],
    stage1_backend_tag: str,
    save_intermediates: bool = False,
    lgbm_n_jobs: int = -1,
) -> dict[str, Any]:
    """
    Fit Stage II on ``df`` (must include ``label`` and ``chrom``), save artifacts, score all rows.

    Steps: autosome train/val split → fit pipeline → validation metrics → dump joblib + feature JSON →
    optional feature importance TSV → write full ranked TSV with true labels for training diagnostics.

    Args:
        df: Merged Stage II frame from :func:`build_stage2_training_frame`.
        models_dir: Directory for ``stage2_model_{tag}.joblib`` and feature JSON.
        predictions_dir: Where ``transcripts.ranked.{tag}.tsv`` is written.
        reports_dir: Optional; used for split lists and importance when ``save_intermediates``.
        autosome_train_range: Inclusive autosome index range (same as Stage I).
        stage1_backend_tag: ``rf`` or ``xgb`` (filename suffix).
        save_intermediates: If true and ``reports_dir`` set, write chrom lists and importance.
        lgbm_n_jobs: Thread/process count for LightGBM.

    Returns:
        Metrics dict printed and optionally logged by caller.
    """
    train_mask, val_mask = split_train_val_masks(df, autosome_train_range, chrom_col="chrom")
    if save_intermediates and reports_dir is not None:
        write_chrom_split_debug_lists(df, train_mask, val_mask, reports_dir, chrom_col="chrom")
    if not train_mask.any() or not val_mask.any():
        raise ValueError("Stage II split produced empty train or validation set (check chrom column).")

    features = list_stage2_model_feature_names(df)
    if not features:
        raise ValueError("No Stage II feature columns after merge (unexpected empty schema).")

    X_train = df.loc[train_mask, features]
    X_val = df.loc[val_mask, features]
    y_train = df.loc[train_mask, "label"].astype(int)
    y_val = df.loc[val_mask, "label"].astype(int)

    clf = build_stage2_lightgbm_pipeline(n_jobs=lgbm_n_jobs)
    clf.fit(X_train, y_train)

    y_pred_val = clf.predict(X_val)
    y_prob_val = stage2_proba_positive_binary(clf, X_val)

    print(
        f"[telos] Stage II: {len(features)} features, "
        f"train n={train_mask.sum()} val n={val_mask.sum()}, "
        f"train labels={np.bincount(y_train.to_numpy()).tolist()} "
        f"val labels={np.bincount(y_val.to_numpy()).tolist()}"
    )

    metrics: dict[str, Any] = {
        "stage1_backend": stage1_backend_tag,
        "n_features": len(features),
        "n_train": int(train_mask.sum()),
        "n_val": int(val_mask.sum()),
        "accuracy": float(accuracy_score(y_val, y_pred_val)),
        "classification_report_val": classification_report(y_val, y_pred_val, digits=4),
    }
    try:
        metrics["roc_auc"] = float(roc_auc_score(y_val, y_prob_val, average="macro"))
        metrics["aupr"] = float(average_precision_score(y_val, y_prob_val, average="macro"))
    except ValueError:
        metrics["roc_auc"] = None
        metrics["aupr"] = None

    print(
        f"[telos] Stage II val: accuracy={metrics['accuracy']:.4f}"
        + (f" roc_auc={metrics['roc_auc']:.4f}" if metrics.get("roc_auc") is not None else "")
        + (f" aupr={metrics['aupr']:.4f}" if metrics.get("aupr") is not None else "")
    )

    models_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir.mkdir(parents=True, exist_ok=True)
    if reports_dir is not None:
        reports_dir.mkdir(parents=True, exist_ok=True)

    model_name = stage2_model_joblib_for_backend(stage1_backend_tag)
    feat_name = stage2_feature_names_json_for_backend(stage1_backend_tag)
    model_path = models_dir / model_name
    joblib.dump(clf, model_path)
    (models_dir / feat_name).write_text(json.dumps(features, indent=2), encoding="utf-8")

    clf_step = clf.named_steps["clf"]
    if hasattr(clf_step, "feature_importances_"):
        fi = pd.DataFrame(
            {"feature": features, "importance": clf_step.feature_importances_}
        ).sort_values("importance", ascending=False)
        metrics["feature_importance_top"] = fi.head(40).to_dict(orient="records")
        if save_intermediates and reports_dir is not None:
            fi.to_csv(
                reports_dir / f"stage2_feature_importance_{stage1_backend_tag}.tsv",
                sep="\t",
                index=False,
            )

    X_full = df[features]
    y_prob_all = stage2_proba_positive_binary(clf, X_full)
    y_pred_all = clf.predict(X_full)
    ranked = pd.DataFrame(
        {
            "transcript_id": df["transcript_id"].values,
            "pred_prob": y_prob_all,
            "pred_label": y_pred_all,
            "label": df["label"].values,
        }
    ).sort_values("pred_prob", ascending=False)
    ranked_name = transcripts_ranked_tsv_for_backend(stage1_backend_tag)
    ranked.to_csv(predictions_dir / ranked_name, sep="\t", index=False)

    return metrics
