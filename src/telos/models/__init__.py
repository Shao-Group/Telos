"""
Model artifact layout and constants for Stage I (RF + XGB) and Stage II (LightGBM ×2).

Filenames here must stay stable so ``train``, ``predict``, and preflight agree on disk layout.
"""

from __future__ import annotations

# Stage I: one RandomForest and one XGBoost classifier per site type (TSS/TES).
STAGE1_BACKEND_RF = "rf"
STAGE1_BACKEND_XGB = "xgb"
STAGE1_BACKENDS: tuple[str, ...] = (STAGE1_BACKEND_RF, STAGE1_BACKEND_XGB)

STAGE1_TSS_RF_MODEL_NAME = "stage1_tss_rf_model.joblib"
STAGE1_TES_RF_MODEL_NAME = "stage1_tes_rf_model.joblib"
STAGE1_TSS_XGB_MODEL_NAME = "stage1_tss_xgb_model.joblib"
STAGE1_TES_XGB_MODEL_NAME = "stage1_tes_xgb_model.joblib"

# Site-level probability columns after dual Stage I scoring
SITE_PROB_COLUMN_RF = "p_site_rf"
SITE_PROB_COLUMN_XGB = "p_site_xgb"

# Stage II: one LightGBM model per Stage I backend (same feature names; different training signal).
STAGE2_MODEL_RF_JOBLIB = "stage2_model_rf.joblib"
STAGE2_MODEL_XGB_JOBLIB = "stage2_model_xgb.joblib"
STAGE2_FEATURE_NAMES_RF_JSON = "stage2_feature_names_rf.json"
STAGE2_FEATURE_NAMES_XGB_JSON = "stage2_feature_names_xgb.json"

# Ranked transcript tables (tab-separated; one per Stage I backend)
TRANSCRIPTS_RANKED_RF_TSV = "transcripts.ranked.rf.tsv"
TRANSCRIPTS_RANKED_XGB_TSV = "transcripts.ranked.xgb.tsv"


def stage1_bundle_path(site_type: str, backend: str) -> str:
    """Return the joblib filename for one Stage I bundle (``TSS``/``TES`` × ``rf``/``xgb``)."""
    st = site_type.strip().upper()
    b = backend.strip().lower()
    if b not in STAGE1_BACKENDS:
        raise ValueError(f"Unknown stage1 backend {backend!r}; expected one of {STAGE1_BACKENDS}")
    if st == "TSS":
        return STAGE1_TSS_RF_MODEL_NAME if b == STAGE1_BACKEND_RF else STAGE1_TSS_XGB_MODEL_NAME
    if st == "TES":
        return STAGE1_TES_RF_MODEL_NAME if b == STAGE1_BACKEND_RF else STAGE1_TES_XGB_MODEL_NAME
    raise ValueError(f"Unknown site_type {site_type!r}")


def stage2_model_joblib_for_backend(backend: str) -> str:
    """Return the Stage II LightGBM model filename for the given Stage I backend tag."""
    b = backend.strip().lower()
    if b == STAGE1_BACKEND_RF:
        return STAGE2_MODEL_RF_JOBLIB
    if b == STAGE1_BACKEND_XGB:
        return STAGE2_MODEL_XGB_JOBLIB
    raise ValueError(f"Unknown stage1 backend {backend!r}")


def stage2_feature_names_json_for_backend(backend: str) -> str:
    """Return the JSON filename listing Stage II feature column order for ``backend``."""
    b = backend.strip().lower()
    if b == STAGE1_BACKEND_RF:
        return STAGE2_FEATURE_NAMES_RF_JSON
    if b == STAGE1_BACKEND_XGB:
        return STAGE2_FEATURE_NAMES_XGB_JSON
    raise ValueError(f"Unknown stage1 backend {backend!r}")


def transcripts_ranked_tsv_for_backend(backend: str) -> str:
    """Return the ranked transcript TSV basename associated with ``backend``."""
    b = backend.strip().lower()
    if b == STAGE1_BACKEND_RF:
        return TRANSCRIPTS_RANKED_RF_TSV
    if b == STAGE1_BACKEND_XGB:
        return TRANSCRIPTS_RANKED_XGB_TSV
    raise ValueError(f"Unknown stage1 backend {backend!r}")
