"""
Prediction command: Stage I scoring + Stage II ranking using saved models.

Does not read ``tmap``; Stage II uses only coverage and site probabilities from the new assembly.
"""

from __future__ import annotations

import logging
from pathlib import Path

from pandas.errors import MergeError

from telos.backends.gtfformat import GtfformatError, run_update_transcript_cov
from telos.config_loader import load_mapping_config
from telos.config_validation import validate_stage_config
from telos.config_models import PredictIO
from telos.models import (
    SITE_PROB_COLUMN_RF,
    SITE_PROB_COLUMN_XGB,
    STAGE1_BACKEND_RF,
    STAGE1_BACKEND_XGB,
    STAGE1_BACKENDS,
    stage2_model_joblib_for_backend,
    transcripts_filtered_gtf_for_backend,
    transcripts_scored_gtf_for_backend,
)
from telos.models.stage1_predict import score_stage1_dataframe, write_sites_scored_tsv
from telos.models.stage2_predict import run_stage2_predict
from telos.models.stage2_train import build_stage2_inference_frame
from telos.pipeline_core import build_stage1_inputs, build_stage1_runtime_config
from telos.postprocess.filter_gtf import filter_gtf_by_transcript_scores
from telos.validation.preflight import (
    PreflightError,
    ensure_run_layout,
    run_preflight_predict,
)

logger = logging.getLogger(__name__)


def run_predict(cfg: PredictIO) -> int:
    """
    Score one BAM/GTF pair with existing Stage I/II artifacts.

    1. Load/validate config; preflight inputs and ``model_dir`` contents.
    2. Rebuild ``df_cov`` / ``df_all`` (same feature pipeline as train).
    3. Load Stage I bundles from ``model_dir``; write ``sites.scored.tsv``.
    4. For RF and XGB, build Stage II inference frame and run LightGBM predict; write both ranked TSVs.

    Returns:
        ``0`` if both backends produce non-empty ranked outputs; ``2`` on any handled failure.
    """
    try:
        cfg_map = load_mapping_config(cfg.config_file)
        validate_stage_config(cfg_map)
    except ValueError as exc:
        logger.error("config error: %s", exc)
        return 2
    try:
        run_preflight_predict(cfg.bam, cfg.gtf, cfg.model_dir, cfg_map)
    except PreflightError as exc:
        logger.error("preflight failed: %s", exc)
        return 2

    layout = ensure_run_layout(cfg.outdir, create_aux_dirs=True)
    runtime_cfg = build_stage1_runtime_config(
        cfg_map,
        cli_no_parallel=cfg.stage1_no_parallel,
        cli_n_workers=cfg.stage1_n_workers,
    )
    try:
        df_cov, df_all = build_stage1_inputs(
            bam=cfg.bam,
            gtf=cfg.gtf,
            runtime_cfg=runtime_cfg,
        )
    except ValueError as exc:
        logger.error("stage1 input prep failed: %s", exc)
        return 2
    try:
        sites_scored = score_stage1_dataframe(df_all, cfg.model_dir)
    except FileNotFoundError as exc:
        logger.error("predict: %s", exc)
        return 2
    sites_path = layout.predictions_dir / "sites.scored.tsv"
    write_sites_scored_tsv(sites_scored, sites_path)

    prob_cols = {STAGE1_BACKEND_RF: SITE_PROB_COLUMN_RF, STAGE1_BACKEND_XGB: SITE_PROB_COLUMN_XGB}
    ranked_paths: dict[str, Path] = {}
    for backend in STAGE1_BACKENDS:
        prob_col = prob_cols[backend]
        try:
            df_stage2 = build_stage2_inference_frame(
                df_cov, df_all, sites_scored, site_prob_column=prob_col
            )
        except (ValueError, KeyError, TypeError, OSError, MergeError) as exc:
            logger.error("Stage II merge failed (%s): %s", backend, exc)
            return 2
        if df_stage2.empty:
            logger.error(
                "Stage II merged table is empty (cov × sites, backend=%s).", backend
            )
            return 2
        try:
            ranked_paths[backend] = run_stage2_predict(
                df_stage2, cfg.model_dir, layout.predictions_dir, stage1_backend_tag=backend
            )
        except (FileNotFoundError, ValueError, KeyError, OSError) as exc:
            logger.error("Stage II predict failed (%s): %s", backend, exc)
            return 2

    target_backends = (
        STAGE1_BACKENDS if cfg.backend == "both" else (cfg.backend,)
    )
    scored_gtfs: dict[str, Path] = {}
    filtered_gtfs: dict[str, Path] = {}
    try:
        for backend in target_backends:
            ranked = ranked_paths[backend]
            scored_path = (
                layout.predictions_dir / transcripts_scored_gtf_for_backend(backend)
            )
            run_update_transcript_cov(cfg.gtf, ranked, scored_path)
            scored_gtfs[backend] = scored_path
            if cfg.min_score is not None:
                filtered_path = (
                    layout.predictions_dir
                    / transcripts_filtered_gtf_for_backend(backend)
                )
                filter_gtf_by_transcript_scores(
                    cfg.gtf,
                    ranked,
                    filtered_path,
                    cfg.min_score,
                )
                filtered_gtfs[backend] = filtered_path
    except (GtfformatError, FileNotFoundError, OSError, TypeError, ValueError) as exc:
        logger.error("GTF output failed: %s", exc)
        return 2


    def _count_rows(path: Path) -> int | None:
        """Return data row count (excluding header) for logging; ``None`` on read errors."""
        try:
            with path.open(encoding="utf-8") as fh:
                return max(0, sum(1 for _ in fh) - 1)
        except OSError:
            return None

    n_rf = _count_rows(ranked_paths[STAGE1_BACKEND_RF])
    n_xgb = _count_rows(ranked_paths[STAGE1_BACKEND_XGB])
    n_ranked_note = (
        f"rf={n_rf}, xgb={n_xgb}" if n_rf is not None and n_xgb is not None else None
    )

    summary_path = layout.reports_dir / "summary.txt"
    summary_lines = [
        "[telos] predict complete (Stage I + Stage II)",
        f"  bam={cfg.bam}",
        f"  gtf={cfg.gtf}",
        f"  model_dir={cfg.model_dir}",
        f"  outdir={layout.root}",
        f"  sites_scored={sites_path}",
    ]
    for b in STAGE1_BACKENDS:
        summary_lines.append(
            f"  stage2_model_{b}={cfg.model_dir / stage2_model_joblib_for_backend(b)}"
        )
        summary_lines.append(f"  transcripts_ranked_{b}={ranked_paths[b]}")
    for b in target_backends:
        summary_lines.append(f"  transcripts_scored_{b}={scored_gtfs[b]}")
        if b in filtered_gtfs:
            summary_lines.append(f"  transcripts_filtered_{b}={filtered_gtfs[b]}")
    if n_ranked_note is not None:
        summary_lines.append(f"  ranked_rows={n_ranked_note}")
    summary_lines.extend(
        [
            f"  summary={summary_path}",
        ]
    )
    summary_path.write_text(
        "\n".join(summary_lines).rstrip() + "\n",
        encoding="utf-8",
    )
    print("\n".join(summary_lines))
    return 0
