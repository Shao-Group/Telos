"""
Train and persist **Stage I** site classifiers (RandomForest or XGBoost) per TSS/TES.

Splits labeled feature tables by autosome policy, fits sklearn estimators, reports validation metrics,
and saves joblib bundles consumed by :mod:`telos.models.stage1_predict`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline

from telos.models import STAGE1_BACKEND_RF, STAGE1_BACKEND_XGB, STAGE1_BACKENDS


def _final_estimator(model: Any) -> Any:
    """Return the underlying classifier from a sklearn ``Pipeline`` or ``model`` itself."""
    if isinstance(model, Pipeline) and "clf" in model.named_steps:
        return model.named_steps["clf"]
    return model


def _binary_positive_proba(clf: Any, proba: np.ndarray) -> np.ndarray:
    """
    Extract probability of class ``1`` from ``predict_proba`` output.

    If only one class was fit, return all ``0.0`` or ``1.0`` according to that class label.
    """
    est = _final_estimator(clf)
    if proba.shape[1] >= 2:
        return proba[:, 1]
    only = int(est.classes_[0])
    return np.full(proba.shape[0], 1.0 if only == 1 else 0.0)


def stage1_positive_probas(model: Any, X: pd.DataFrame) -> np.ndarray:
    """Vector of P(class=1) for each row of ``X`` using ``model.predict_proba``."""
    return _binary_positive_proba(model, model.predict_proba(X))


META_KEYS = frozenset(
    {
        "transcript_id",
        "site_type",
        "chrom",
        "position",
        "strand",
        "label",
    }
)


def build_stage1_feature_matrix(df: pd.DataFrame, feature_names: list[str]) -> pd.DataFrame:
    """Reindex ``df`` to ``feature_names``, coerce numeric, fill NaN with ``0.0``."""
    X = df.reindex(columns=feature_names, fill_value=0.0)
    return X.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def infer_stage1_feature_columns(df: pd.DataFrame) -> list[str]:
    """All numeric columns except :data:`META_KEYS` (transcript_id, site_type, chrom, position, strand, label)."""
    cols: list[str] = []
    for c in df.columns:
        if c in META_KEYS:
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            cols.append(c)
    return cols


def build_stage1_xgb_pipeline(xgb_config: dict[str, Any], random_state: int) -> Pipeline:
    """
    Single-step sklearn ``Pipeline`` wrapping :class:`xgboost.XGBClassifier` with YAML-tunable hyperparameters.

    Raises:
        ImportError: If xgboost is not installed.
    """
    try:
        from xgboost import XGBClassifier  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "Stage I XGBoost training requires xgboost. Install with: pip install xgboost"
        ) from exc

    cfg = dict(xgb_config) if xgb_config else {}
    n_jobs = int(cfg.pop("n_jobs", 8))
    xgb_model = XGBClassifier(
        n_estimators=int(cfg.get("n_estimators", 300)),
        max_depth=int(cfg.get("max_depth", 6)),
        learning_rate=float(cfg.get("learning_rate", 0.1)),
        subsample=float(cfg.get("subsample", 0.8)),
        colsample_bytree=float(cfg.get("colsample_bytree", 0.8)),
        reg_lambda=float(cfg.get("reg_lambda", 3)),
        reg_alpha=float(cfg.get("reg_alpha", 0.5)),
        objective="binary:logistic",
        random_state=random_state,
        n_jobs=n_jobs,
        eval_metric="logloss",
    )
    return Pipeline([("clf", xgb_model)])


def train_stage1_site_classifier(
    df_labeled: pd.DataFrame,
    site_type: str,
    autosome_train_range: tuple[int, int],
    *,
    backend: str,
    rf_config: dict[str, Any],
    xgb_config: dict[str, Any],
    random_state: int = 42,
) -> tuple[dict[str, Any], Any, list[str]]:
    """
    Fit one Stage I classifier for ``site_type`` (``TSS``/``TES``) and ``backend`` (``rf``/``xgb``).

    Splits ``df_labeled`` with :func:`~telos.models.chrom_split.split_train_val_masks`, infers feature
    columns, trains, evaluates on validation (accuracy, PR-AUC, ROC-AUC when both classes present).

    Returns:
        ``(metrics_dict, fitted_model, feature_names)``.
    """
    from telos.models.chrom_split import split_train_val_masks

    if df_labeled.empty:
        raise ValueError(f"No rows to train Stage I {site_type}.")
    b = backend.strip().lower()
    if b not in STAGE1_BACKENDS:
        raise ValueError(f"backend must be one of {STAGE1_BACKENDS}, got {backend!r}")
    train_mask, val_mask = split_train_val_masks(df_labeled, autosome_train_range)
    if not train_mask.any():
        raise ValueError(
            "Train split is empty; check split_policy and that chrom names map to autosomes in range "
            "(e.g. chr1-22, Ensembl 1-22, RefSeq NC_000001-NC_000022)."
        )
    if not val_mask.any():
        raise ValueError(
            "Validation split is empty; check split_policy (validation requires rows outside train range "
            "or on non-autosomal / unmapped chrom names)."
        )

    feature_names = infer_stage1_feature_columns(df_labeled)
    if not feature_names:
        raise ValueError("No numeric feature columns found for Stage I training.")

    X_tr = build_stage1_feature_matrix(df_labeled.loc[train_mask], feature_names)
    y_tr = df_labeled.loc[train_mask, "label"].astype(int)
    X_va = build_stage1_feature_matrix(df_labeled.loc[val_mask], feature_names)
    y_va = df_labeled.loc[val_mask, "label"].astype(int)

    if b == STAGE1_BACKEND_RF:
        clf = RandomForestClassifier(
            n_estimators=int(rf_config.get("n_estimators", 300)),
            max_depth=rf_config.get("max_depth"),
            max_features=rf_config.get("max_features", "sqrt"),
            random_state=random_state,
            n_jobs=int(rf_config.get("n_jobs", -1)),
            class_weight=rf_config.get("class_weight", "balanced_subsample"),
        )
    elif b == STAGE1_BACKEND_XGB:
        clf = build_stage1_xgb_pipeline(xgb_config, random_state)
    else:
        raise ValueError(f"Unknown Stage I backend {backend!r}")

    clf.fit(X_tr, y_tr)

    proba_va = clf.predict_proba(X_va)
    y_prob = _binary_positive_proba(clf, proba_va)
    y_pred = (y_prob >= 0.5).astype(int)

    metrics: dict[str, Any] = {
        "site_type": site_type.upper(),
        "backend": b,
        "n_train": int(len(y_tr)),
        "n_val": int(len(y_va)),
        "train_pos_rate": float(y_tr.mean()),
        "val_pos_rate": float(y_va.mean()),
        "accuracy": float(accuracy_score(y_va, y_pred)),
        "precision": float(precision_score(y_va, y_pred, zero_division=0)),
        "recall": float(recall_score(y_va, y_pred, zero_division=0)),
        "f1": float(f1_score(y_va, y_pred, zero_division=0)),
        "confusion_matrix": confusion_matrix(y_va, y_pred).tolist(),
    }
    if y_va.nunique() < 2:
        metrics["aupr"] = None
        metrics["roc_auc"] = None
    else:
        metrics["aupr"] = float(average_precision_score(y_va, y_prob, average="macro"))
        metrics["roc_auc"] = float(roc_auc_score(y_va, y_prob, average="macro"))

    return metrics, clf, feature_names


def save_stage1_bundle(
    out_path: Path,
    site_type: str,
    backend: str,
    model: Any,
    feature_names: list[str],
) -> None:
    """Serialize ``model`` + ``feature_names`` + metadata with :func:`joblib.dump` to ``out_path``."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bundle = {
        "schema_version": 2,
        "site_type": site_type.upper(),
        "backend": backend,
        "feature_names": feature_names,
        "model": model,
    }
    joblib.dump(bundle, out_path)


def load_stage1_bundle(path: Path) -> dict[str, Any]:
    """Load joblib bundle; validate keys ``model`` and ``feature_names``."""
    obj = joblib.load(path)
    if not isinstance(obj, dict) or "model" not in obj or "feature_names" not in obj:
        raise ValueError(f"Invalid Stage I bundle: {path}")
    return obj


def write_train_metrics(metrics_by_site: dict[str, dict[str, Any]], out_json: Path) -> None:
    """Pretty-print nested metrics dict to JSON on disk (optional training diagnostic)."""
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(metrics_by_site, indent=2, sort_keys=True), encoding="utf-8")
