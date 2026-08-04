"""
Training command: end-to-end Stage I + Stage II for one BAM/GTF/ref/tmap bundle.

Loads YAML, validates shape, runs preflight, builds Stage I tables, trains four Stage I bundles
(TSS/TES × RF/XGB), scores all candidates, trains two LightGBM Stage II models (RF- and XGB-driven
site probabilities), writes ranked transcript TSVs for both.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
from pandas.errors import MergeError

from telos.config_loader import get_nested, load_mapping_config
from telos.config_validation import validate_stage_config
from telos.config_models import TrainIO
from telos.labels.site_labels import label_sites_by_proximity, reference_sites_from_gtf
from telos.models import (
    SITE_PROB_COLUMN_RF,
    SITE_PROB_COLUMN_XGB,
    STAGE1_BACKEND_RF,
    STAGE1_BACKEND_XGB,
    STAGE1_BACKENDS,
    stage1_bundle_path,
    stage2_model_joblib_for_backend,
    transcripts_ranked_tsv_for_backend,
)
from telos.models.chrom_split import parse_split_policy
from telos.models.stage1_predict import score_stage1_dataframe, write_sites_scored_tsv
from telos.models.stage1_train import (
    save_stage1_bundle,
    train_stage1_site_classifier,
)
from telos.models.stage2_train import build_stage2_training_frame, train_and_save_stage2
from telos.pipeline_core import (
    build_stage1_inputs,
    build_stage1_inputs_multi_gtf,
    build_stage1_runtime_config,
)
from telos.validation.preflight import (
    PreflightError,
    ensure_run_layout,
    run_preflight_train,
)

logger = logging.getLogger(__name__)

_TRAIN_METRICS_COLUMNS = [
    "stage",
    "site_type",
    "backend",
    "n_train",
    "n_val",
    "train_pos_rate",
    "val_pos_rate",
    "accuracy",
    "precision",
    "recall",
    "f1",
    "aupr",
    "roc_auc",
    "tn",
    "fp",
    "fn",
    "tp",
]


def _confusion_parts(value: Any) -> tuple[int | None, int | None, int | None, int | None]:
    if (
        isinstance(value, list)
        and len(value) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in value)
    ):
        return int(value[0][0]), int(value[0][1]), int(value[1][0]), int(value[1][1])
    return None, None, None, None


def _write_train_metrics(metrics_payload: dict[str, Any], out_csv: Path) -> None:
    rows: list[dict[str, Any]] = []
    for site_type in ("tss", "tes"):
        for backend in STAGE1_BACKENDS:
            metrics = metrics_payload[f"{site_type}_{backend}"]
            tn, fp, fn, tp = _confusion_parts(metrics.get("confusion_matrix"))
            rows.append(
                {
                    "stage": "stage1",
                    "site_type": metrics.get("site_type", site_type.upper()),
                    "backend": metrics.get("backend", backend),
                    "n_train": metrics.get("n_train"),
                    "n_val": metrics.get("n_val"),
                    "train_pos_rate": metrics.get("train_pos_rate"),
                    "val_pos_rate": metrics.get("val_pos_rate"),
                    "accuracy": metrics.get("accuracy"),
                    "precision": metrics.get("precision"),
                    "recall": metrics.get("recall"),
                    "f1": metrics.get("f1"),
                    "aupr": metrics.get("aupr"),
                    "roc_auc": metrics.get("roc_auc"),
                    "tn": tn,
                    "fp": fp,
                    "fn": fn,
                    "tp": tp,
                }
            )
    for backend in STAGE1_BACKENDS:
        metrics = metrics_payload[f"stage2_{backend}"]
        rows.append(
            {
                "stage": "stage2",
                "backend": metrics.get("stage1_backend", backend),
                "n_train": metrics.get("n_train"),
                "n_val": metrics.get("n_val"),
                "accuracy": metrics.get("accuracy"),
                "aupr": metrics.get("aupr"),
                "roc_auc": metrics.get("roc_auc"),
            }
        )
    pd.DataFrame(rows, columns=_TRAIN_METRICS_COLUMNS).to_csv(out_csv, index=False)


def run_train(cfg: TrainIO) -> int:
    """
    Run the full Telos training pipeline for one dataset.

    **Steps (high level)**

    1. Load and validate Stage I config; preflight BAM/GTF/tmap.
    2. ``ensure_run_layout`` under ``cfg.outdir`` (models + predictions dirs).
    3. Build Stage I runtime + ``df_cov`` / ``df_all`` via :mod:`telos.pipeline_core`.
    4. For each site type (TSS, TES) and each Stage I backend, label rows from ``ref_gtf``, train,
       and save a joblib bundle under ``models/``.
    5. Score every candidate row with both backends; write ``predictions/sites.scored.tsv``.
    6. For each backend, merge cov + sites + ``tmap`` labels into Stage II frame; train LightGBM;
       write ``models/stage2_*.joblib`` and ranked ``predictions/transcripts.ranked.*.tsv``.

    Returns:
        ``0`` on success; ``2`` on config/preflight/stage failure; ``3`` is reserved upstream for
        empty candidates (some code paths print and return ``2`` here for merge failures—see messages).
    """
    try:
        cfg_map = load_mapping_config(cfg.config_file)
        validate_stage_config(cfg_map)
    except ValueError as exc:
        logger.error("config error: %s", exc)
        return 2
    try:
        run_preflight_train(cfg.bam, cfg.gtf, cfg.ref_gtf, cfg.tmap, cfg_map)
    except PreflightError as exc:
        logger.error("preflight failed: %s", exc)
        return 2

    layout = ensure_run_layout(cfg.outdir, create_aux_dirs=True)
    runtime_cfg = build_stage1_runtime_config(
        cfg_map,
        cli_no_parallel=cfg.stage1_no_parallel,
        cli_n_workers=cfg.stage1_n_workers,
    )
    gtf_train_pool = [cfg.gtf, *(list(cfg.gtf_pool) if cfg.gtf_pool else [])]
    tmap_train_pool = [cfg.tmap, *(list(cfg.tmap_pool) if cfg.tmap_pool else [])]
    if len(tmap_train_pool) != len(gtf_train_pool):
        logger.error(
            "Stage II pooled supervision requires one tmap per training gtf "
            "(got gtfs=%s tmaps=%s).",
            len(gtf_train_pool),
            len(tmap_train_pool),
        )
        return 2
    try:
        if len(gtf_train_pool) == 1:
            df_cov, df_all = build_stage1_inputs(
                bam=cfg.bam,
                gtf=cfg.gtf,
                runtime_cfg=runtime_cfg,
            )
        else:
            logger.info("pooled training GTFs: %s", len(gtf_train_pool))
            df_cov, df_all = build_stage1_inputs_multi_gtf(
                bam=cfg.bam,
                gtfs=gtf_train_pool,
                runtime_cfg=runtime_cfg,
            )
    except ValueError as exc:
        logger.error("stage1 input prep failed: %s", exc)
        return 2

    split_policy = str(
        cfg.split_policy
        if cfg.split_policy is not None
        else get_nested(cfg_map, ["stage1", "training", "split_policy"], "chr1-10")
    )
    try:
        autosome_train_range = parse_split_policy(split_policy)
    except ValueError as exc:
        logger.error("invalid split_policy: %s", exc)
        return 2
    logger.info(
        "train/val split_policy=%s (train autosomes %s-%s; all other contigs = validation)",
        split_policy,
        autosome_train_range[0],
        autosome_train_range[1],
    )

    tol = int(get_nested(cfg_map, ["stage1", "training", "site_label_tolerance_bp"], 50))
    rf_cfg = dict(get_nested(cfg_map, ["stage1", "training", "random_forest"], {}) or {})
    xgb_cfg = dict(get_nested(cfg_map, ["stage1", "training", "xgboost"], {}) or {})
    seed = int(get_nested(cfg_map, ["stage1", "training", "random_state"], 42))
    lgbm_n_jobs = int(get_nested(cfg_map, ["stage1", "training", "lightgbm", "n_jobs"], -1))
    if cfg.n_jobs is not None:
        rf_cfg["n_jobs"] = int(cfg.n_jobs)
        xgb_cfg["n_jobs"] = int(cfg.n_jobs)
        lgbm_n_jobs = int(cfg.n_jobs)
        logger.info("model n_jobs override=%s (RF, XGBoost, LightGBM)", cfg.n_jobs)
    ref_df = reference_sites_from_gtf(cfg.ref_gtf)
    metrics_payload: dict[str, Any] = {}

    for st in ("TSS", "TES"):
        labeled = df_all[df_all["site_type"].str.upper() == st].copy()
        if labeled.empty:
            logger.error("no %s feature rows; cannot train Stage I.", st)
            return 2
        labeled["label"] = label_sites_by_proximity(labeled, ref_df, st, tol)
        for backend in STAGE1_BACKENDS:
            try:
                m, clf, feats = train_stage1_site_classifier(
                    labeled,
                    st,
                    autosome_train_range,
                    backend=backend,
                    rf_config=rf_cfg,
                    xgb_config=xgb_cfg,
                    random_state=seed,
                )
            except ImportError as exc:
                logger.error("Stage I training failed (%s, %s): %s", st, backend, exc)
                return 2
            except ValueError as exc:
                logger.error("Stage I training failed (%s, %s): %s", st, backend, exc)
                return 2
            metrics_payload[f"{st.lower()}_{backend}"] = m
            fname = stage1_bundle_path(st, backend)
            save_stage1_bundle(layout.models_dir / fname, st, backend, clf, feats)

    sites_scored = score_stage1_dataframe(df_all, layout.models_dir)
    sites_path = layout.predictions_dir / "sites.scored.tsv"
    write_sites_scored_tsv(sites_scored, sites_path)

    prob_cols = {STAGE1_BACKEND_RF: SITE_PROB_COLUMN_RF, STAGE1_BACKEND_XGB: SITE_PROB_COLUMN_XGB}
    for backend in STAGE1_BACKENDS:
        prob_col = prob_cols[backend]
        try:
            stage2_parts = [
                build_stage2_training_frame(df_cov, df_all, sites_scored, tm, site_prob_column=prob_col)
                for tm in tmap_train_pool
            ]
            df_stage2 = (
                stage2_parts[0]
                if len(stage2_parts) == 1
                else pd.concat(stage2_parts, axis=0, ignore_index=True)
            )
        except (ValueError, KeyError, TypeError, OSError, MergeError) as exc:
            logger.error("Stage II feature merge failed (%s): %s", backend, exc)
            return 2
        if df_stage2.empty:
            logger.error(
                "Stage II merged table is empty after inner-joins "
                "(cov + site scores + tmap labels, backend=%s). "
                "Check transcript_id overlap between assembly GTF, Stage I candidates, "
                "and bundle .tmap qry_id.",
                backend,
            )
            return 2

        try:
            metrics_payload[f"stage2_{backend}"] = train_and_save_stage2(
                df_stage2,
                layout.models_dir,
                layout.predictions_dir,
                autosome_train_range=autosome_train_range,
                stage1_backend_tag=backend,
                lgbm_n_jobs=lgbm_n_jobs,
            )
        except ImportError as exc:
            logger.error("Stage II requires lightgbm: %s", exc)
            return 2
        except (ValueError, KeyError, OSError, RuntimeError) as exc:
            logger.error("Stage II training failed (%s): %s", backend, exc)
            return 2

    metrics_path = layout.reports_dir / "train_metrics.csv"
    _write_train_metrics(metrics_payload, metrics_path)

    ranked_rf = layout.predictions_dir / transcripts_ranked_tsv_for_backend(STAGE1_BACKEND_RF)
    ranked_xgb = layout.predictions_dir / transcripts_ranked_tsv_for_backend(STAGE1_BACKEND_XGB)
    s1_list = ", ".join(
        str(layout.models_dir / stage1_bundle_path(st, b))
        for st in ("TSS", "TES")
        for b in STAGE1_BACKENDS
    )
    s2_list = ", ".join(
        str(layout.models_dir / stage2_model_joblib_for_backend(b)) for b in STAGE1_BACKENDS
    )
    summary_lines = [
        "[telos] train complete",
        f"  bam={cfg.bam}",
        f"  gtf={cfg.gtf}",
        f"  ref_gtf={cfg.ref_gtf}",
        f"  tmap={cfg.tmap}",
        f"  split_policy={split_policy}",
        f"  outdir={layout.root}",
        f"  stage1_models={s1_list}",
        f"  sites_scored={sites_path}",
        f"  stage2_models={s2_list}",
        f"  transcripts_ranked_rf={ranked_rf}",
        f"  transcripts_ranked_xgb={ranked_xgb}",
        f"  train_metrics={metrics_path}",
    ]
    print("\n".join(summary_lines))
    return 0
