"""
Shared Stage I data preparation used by both training and prediction.

Turns a BAM + assembly GTF + parsed YAML fragment into (1) a per-transcript coverage/features
dataframe and (2) a per-candidate-site feature dataframe ready for Stage I models.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

from telos.backends.gtfformat import GtfformatError
from telos.candidates.load import load_candidates, load_transcript_cov_dataframe
from telos.config_loader import get_nested
from telos.features.stage1 import Stage1FeatureConfig, compute_stage1_features


@dataclass(frozen=True)
class Stage1RuntimeConfig:
    """
    Concrete parameters for :func:`build_stage1_inputs` after merging YAML and CLI overrides.

    ``feature_cfg`` holds numeric windows and feature toggles. ``parallel`` is the effective flag
    after CLI can disable YAML parallelism. ``parallel_min_sites`` gates pool use for small jobs.
    ``n_workers`` is ``None`` to mean “default pool sizing” in the feature layer.
    """

    feature_cfg: Stage1FeatureConfig
    parallel: bool
    parallel_min_sites: int
    n_workers: int | None
    cache_dir: Path | None


def _file_fingerprint(path: Path) -> dict[str, Any]:
    """Cheap change detector for cache keys: absolute path + size + mtime_ns."""
    p = path.resolve()
    st = p.stat()
    return {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def _stage1_cache_key(*, bam: Path, gtf: Path, runtime_cfg: Stage1RuntimeConfig) -> str:
    """Stable hash for Stage I feature/cov cache entries."""
    payload = {
        "bam": _file_fingerprint(bam),
        "bam_index": _file_fingerprint(Path(str(bam) + ".bai")) if Path(str(bam) + ".bai").is_file() else None,
        "gtf": _file_fingerprint(gtf),
        "feature_cfg": {
            "window_size": runtime_cfg.feature_cfg.window_size,
            "density_window": runtime_cfg.feature_cfg.density_window,
            "coverage_window": runtime_cfg.feature_cfg.coverage_window,
            "soft_clip_window": runtime_cfg.feature_cfg.soft_clip_window,
            "min_mapq": runtime_cfg.feature_cfg.min_mapq,
            "splice_site_window": runtime_cfg.feature_cfg.splice_site_window,
            "gradient_analysis_range": runtime_cfg.feature_cfg.gradient_analysis_range,
            "extended_window": runtime_cfg.feature_cfg.extended_window,
            "k_values": list(runtime_cfg.feature_cfg.k_values),
        },
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:20]


def build_stage1_runtime_config(
    cfg_map: dict[str, Any], *, cli_no_parallel: bool, cli_n_workers: int | None
) -> Stage1RuntimeConfig:
    """
    Build :class:`Stage1RuntimeConfig` from the ``stage1.feature_extraction`` and ``stage1.features`` sections.

    Reads window sizes, mapq cutoff, k-mer list, YAML ``parallel`` / ``parallel_min_sites`` /
    ``n_workers``. If ``cli_no_parallel`` is true, forces ``parallel=False``. If ``cli_n_workers`` is
    not ``None``, it overrides YAML ``n_workers``.

    Returns:
        Frozen config passed to :func:`compute_stage1_features`.
    """
    feat_cfg = Stage1FeatureConfig(
        window_size=int(get_nested(cfg_map, ["stage1", "feature_extraction", "window_size"], 100)),
        density_window=int(get_nested(cfg_map, ["stage1", "feature_extraction", "density_window"], 100)),
        coverage_window=int(get_nested(cfg_map, ["stage1", "feature_extraction", "coverage_window"], 100)),
        soft_clip_window=int(get_nested(cfg_map, ["stage1", "feature_extraction", "soft_clip_window"], 10)),
        min_mapq=int(get_nested(cfg_map, ["stage1", "feature_extraction", "min_mapq"], 10)),
        splice_site_window=int(get_nested(cfg_map, ["stage1", "feature_extraction", "splice_site_window"], 300)),
        gradient_analysis_range=int(
            get_nested(cfg_map, ["stage1", "feature_extraction", "gradient_analysis_range"], 100)
        ),
        extended_window=int(get_nested(cfg_map, ["stage1", "feature_extraction", "extended_window"], 1000)),
        k_values=tuple(get_nested(cfg_map, ["stage1", "features", "k_values"], [3])),
    )
    stage1_parallel_yaml = bool(get_nested(cfg_map, ["stage1", "feature_extraction", "parallel"], True))
    stage1_parallel = stage1_parallel_yaml and not cli_no_parallel
    stage1_parallel_min = int(get_nested(cfg_map, ["stage1", "feature_extraction", "parallel_min_sites"], 50))
    stage1_n_workers_yaml = get_nested(cfg_map, ["stage1", "feature_extraction", "n_workers"], None)
    stage1_n_workers = cli_n_workers if cli_n_workers is not None else stage1_n_workers_yaml
    cache_env = os.environ.get("TELOS_STAGE1_CACHE_DIR", "").strip()
    cache_cfg_raw = get_nested(cfg_map, ["stage1", "feature_extraction", "cache_dir"], None)
    cache_raw = cache_env if cache_env else cache_cfg_raw
    cache_dir = Path(str(cache_raw)).expanduser().resolve() if cache_raw else None
    return Stage1RuntimeConfig(
        feature_cfg=feat_cfg,
        parallel=stage1_parallel,
        parallel_min_sites=stage1_parallel_min,
        n_workers=stage1_n_workers,
        cache_dir=cache_dir,
    )


def stage1_cache_dir_for(
    *,
    bam: Path,
    gtf: Path,
    runtime_cfg: Stage1RuntimeConfig,
    cache_root: Path,
) -> Path:
    """Resolved cache directory for ``(bam, gtf, runtime_cfg)`` under ``cache_root``."""
    return (cache_root / _stage1_cache_key(bam=bam, gtf=gtf, runtime_cfg=runtime_cfg)).resolve()


def build_stage1_inputs(
    *,
    bam: Path,
    gtf: Path,
    runtime_cfg: Stage1RuntimeConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Produce Stage I inputs: transcript-level cov table and candidate-site feature table.

    1. ``load_candidates(gtf)`` — if empty, raises ``ValueError`` (no TSS/TES candidates).
    2. ``load_transcript_cov_dataframe`` — wraps ``gtfformat`` coverage extraction; ``GtfformatError``
       is turned into ``ValueError`` with context.
    3. ``compute_stage1_features`` — BAM scans at each candidate using ``runtime_cfg``; if no rows,
       raises ``ValueError``.

    Returns:
        ``(df_cov, df_all)`` where ``df_cov`` is one row per transcript (Stage II join key) and
        ``df_all`` is one row per candidate site with engineered features.
    """
    cache_key = _stage1_cache_key(bam=bam, gtf=gtf, runtime_cfg=runtime_cfg)
    cache_root = runtime_cfg.cache_dir
    cache_hit = False
    if cache_root is not None:
        cache_dir = (cache_root / cache_key).resolve()
        cov_pkl = cache_dir / "df_cov.pkl"
        all_pkl = cache_dir / "df_all.pkl"
        if cov_pkl.is_file() and all_pkl.is_file():
            try:
                df_cov = pd.read_pickle(cov_pkl)
                df_all = pd.read_pickle(all_pkl)
                cache_hit = True
                print(f"[telos] Stage I cache hit: {cache_dir}")
                return df_cov, df_all
            except Exception:
                pass

    candidates = load_candidates(gtf)
    if not candidates:
        raise ValueError(
            "No candidate sites from input GTF; check transcript features and coordinates."
        )
    try:
        df_cov = load_transcript_cov_dataframe(gtf)
    except GtfformatError as exc:
        raise ValueError(f"transcript coverage extraction failed: {exc}") from exc
    rows = compute_stage1_features(
        bam,
        candidates,
        runtime_cfg.feature_cfg,
        parallel=runtime_cfg.parallel,
        parallel_min_sites=runtime_cfg.parallel_min_sites,
        n_workers=runtime_cfg.n_workers,
        progress=False,
    )
    if not rows:
        raise ValueError("Stage I feature extraction returned no rows.")
    df_all = pd.DataFrame(rows)
    if cache_root is not None and not cache_hit:
        cache_dir = (cache_root / cache_key).resolve()
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            df_cov.to_pickle(cache_dir / "df_cov.pkl")
            df_all.to_pickle(cache_dir / "df_all.pkl")
            (cache_dir / "meta.json").write_text(
                json.dumps(
                    {
                        "bam": str(bam.resolve()),
                        "gtf": str(gtf.resolve()),
                        "cache_key": cache_key,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            print(f"[telos] Stage I cache write: {cache_dir}")
        except OSError:
            pass
    return df_cov, df_all


def build_stage1_inputs_multi_gtf(
    *,
    bam: Path,
    gtfs: list[Path],
    runtime_cfg: Stage1RuntimeConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build pooled Stage I/II inputs by concatenating per-GTF inputs from multiple assemblies.

    Reuses :func:`build_stage1_inputs` for each GTF so cache semantics stay identical.
    """
    if not gtfs:
        raise ValueError("gtfs list is empty for pooled training.")
    cov_parts: list[pd.DataFrame] = []
    feat_parts: list[pd.DataFrame] = []
    for gtf in gtfs:
        df_cov, df_all = build_stage1_inputs(bam=bam, gtf=gtf, runtime_cfg=runtime_cfg)
        cov_parts.append(df_cov)
        feat_parts.append(df_all)
    cov = pd.concat(cov_parts, axis=0, ignore_index=True)
    feats = pd.concat(feat_parts, axis=0, ignore_index=True)
    feats = feats.drop_duplicates(
        subset=["transcript_id", "site_type", "chrom", "position", "strand"],
        keep="first",
    ).reset_index(drop=True)
    return cov, feats
