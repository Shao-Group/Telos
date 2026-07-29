"""
Prediction command: Stage I scoring + Stage II ranking using saved models.

Does not read ``tmap``; Stage II uses only coverage and site probabilities from the new assembly.
"""

from __future__ import annotations

from pathlib import Path

from pandas.errors import MergeError

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
)
from telos.models.stage1_predict import score_stage1_dataframe, write_sites_scored_tsv
from telos.models.stage2_predict import run_stage2_predict
from telos.models.stage2_train import build_stage2_inference_frame
from telos.pipeline_core import build_stage1_inputs, build_stage1_runtime_config
from telos.validation.preflight import (
    PreflightError,
    ensure_run_layout,
    run_preflight_predict,
)


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
        print(f"[telos] config error: {exc}")
        return 2
    try:
        run_preflight_predict(cfg.bam, cfg.gtf, cfg.model_dir, cfg_map)
    except PreflightError as exc:
        print(f"[telos] preflight failed: {exc}")
        return 2

    layout = ensure_run_layout(cfg.outdir, save_intermediates=False, create_aux_dirs=False)
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
        print(f"[telos] stage1 input prep failed: {exc}")
        return 2
    try:
        sites_scored = score_stage1_dataframe(df_all, cfg.model_dir)
    except FileNotFoundError as exc:
        print(f"[telos] predict: {exc}")
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
            print(f"[telos] Stage II merge failed ({backend}): {exc}")
            return 2
        if df_stage2.empty:
            print(f"[telos] Stage II merged table is empty (cov × sites, backend={backend}).")
            return 2
        try:
            ranked_paths[backend] = run_stage2_predict(
                df_stage2, cfg.model_dir, layout.predictions_dir, stage1_backend_tag=backend
            )
        except (FileNotFoundError, ValueError, KeyError, OSError) as exc:
            print(f"[telos] Stage II predict failed ({backend}): {exc}")
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

    print("[telos] predict complete (Stage I + Stage II)")
    print(f"  bam={cfg.bam}")
    print(f"  gtf={cfg.gtf}")
    print(f"  model_dir={cfg.model_dir}")
    print(f"  outdir={cfg.outdir}")
    print(f"  sites_scored={sites_path}")
    for b in STAGE1_BACKENDS:
        print(f"  stage2_model_{b}={cfg.model_dir / stage2_model_joblib_for_backend(b)}")
        print(f"  transcripts_ranked_{b}={ranked_paths[b]}")
    if n_ranked_note is not None:
        print(f"  ranked_rows={n_ranked_note}")
    return 0
