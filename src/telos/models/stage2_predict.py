"""
Run Stage II **inference**: load LightGBM pipeline + feature list, score merged transcript dataframe.
"""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import pandas as pd

from telos.models import (
    stage2_feature_names_json_for_backend,
    stage2_model_joblib_for_backend,
    transcripts_ranked_tsv_for_backend,
)
from telos.models.stage2_train import stage2_proba_positive_binary


def run_stage2_predict(
    df_inference: pd.DataFrame,
    model_dir: Path,
    predictions_dir: Path,
    *,
    stage1_backend_tag: str,
) -> Path:
    """
    Predict transcript scores for every row in ``df_inference`` using the model for ``stage1_backend_tag``.

    Loads ``stage2_*_{tag}.joblib`` and matching JSON feature list from ``model_dir``, aligns columns,
    writes sorted ``transcripts.ranked.{rf,xgb}.tsv`` under ``predictions_dir``, returns output path.

    Does not require a ``label`` column (predict-only path).
    """
    model_name = stage2_model_joblib_for_backend(stage1_backend_tag)
    feat_name = stage2_feature_names_json_for_backend(stage1_backend_tag)
    clf = joblib.load(model_dir / model_name)
    features = json.loads((model_dir / feat_name).read_text(encoding="utf-8"))
    X = df_inference.reindex(columns=list(features), fill_value=0.0)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y_prob = stage2_proba_positive_binary(clf, X)
    y_pred = clf.predict(X)
    out = pd.DataFrame(
        {
            "transcript_id": df_inference["transcript_id"].values,
            "pred_prob": y_prob,
            "pred_label": y_pred,
        }
    ).sort_values("pred_prob", ascending=False)
    predictions_dir.mkdir(parents=True, exist_ok=True)
    path = predictions_dir / transcripts_ranked_tsv_for_backend(stage1_backend_tag)
    out.to_csv(path, sep="\t", index=False)
    return path
