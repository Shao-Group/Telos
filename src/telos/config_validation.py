"""
Strict, shallow validation of Stage I and benchmark YAML shapes.

Unknown keys are rejected so typos fail fast. This is **not** a full semantic validator (paths are
not checked here; preflight handles existence).
"""

from __future__ import annotations

from typing import Any


def _expect_dict(v: Any, where: str) -> dict[str, Any]:
    """Return ``v`` if it is a ``dict``; else raise ``ValueError`` naming ``where``."""
    if not isinstance(v, dict):
        raise ValueError(f"{where} must be a mapping/object")
    return v


def _expect_bool(v: Any, where: str) -> None:
    """Require strict Python ``bool`` (not ``0``/``1``)."""
    if not isinstance(v, bool):
        raise ValueError(f"{where} must be a boolean")


def _expect_number(v: Any, where: str) -> None:
    """Require ``int`` or ``float`` excluding ``bool`` subclasses."""
    if not isinstance(v, (int, float)) or isinstance(v, bool):
        raise ValueError(f"{where} must be a number")


def _expect_int_or_none(v: Any, where: str) -> None:
    """Allow ``None`` or a non-boolean ``int``."""
    if v is None:
        return
    if not isinstance(v, int) or isinstance(v, bool):
        raise ValueError(f"{where} must be an integer or null")


def _expect_list_int(v: Any, where: str) -> None:
    """Non-empty list of non-boolean ints (e.g. ``k_values``)."""
    if not isinstance(v, list) or not v:
        raise ValueError(f"{where} must be a non-empty list")
    for i, x in enumerate(v):
        if not isinstance(x, int) or isinstance(x, bool):
            raise ValueError(f"{where}[{i}] must be an integer")


def _reject_unknown_keys(obj: dict[str, Any], allowed: set[str], where: str) -> None:
    """If ``obj`` contains keys not in ``allowed``, raise with sorted unknown key list."""
    unknown = sorted(set(obj.keys()) - allowed)
    if unknown:
        raise ValueError(f"Unknown config keys at {where}: {unknown}")


def validate_stage_config(cfg_map: dict[str, Any]) -> None:
    """
    Validate the ``stage1.*`` subtree allowed keys and scalar types for train/predict configs.

    Empty ``cfg_map`` is a no-op (allows callers to pass ``{}`` before defaults merge). Raises
    :class:`ValueError` with a path-like message on the first violation.
    """
    if not cfg_map:
        return
    _reject_unknown_keys(cfg_map, {"stage1"}, "root")
    stage1 = _expect_dict(cfg_map.get("stage1", {}), "stage1")
    _reject_unknown_keys(stage1, {"training", "feature_extraction", "features"}, "stage1")

    training = _expect_dict(stage1.get("training", {}), "stage1.training")
    _reject_unknown_keys(
        training,
        {
            "site_label_tolerance_bp",
            "split_policy",
            "random_state",
            "random_forest",
            "xgboost",
            "lightgbm",
        },
        "stage1.training",
    )
    if "site_label_tolerance_bp" in training:
        _expect_number(training["site_label_tolerance_bp"], "stage1.training.site_label_tolerance_bp")
    if "split_policy" in training and not isinstance(training["split_policy"], str):
        raise ValueError("stage1.training.split_policy must be a string")
    if "random_state" in training:
        _expect_number(training["random_state"], "stage1.training.random_state")
    if "random_forest" in training:
        _expect_dict(training["random_forest"], "stage1.training.random_forest")
    if "xgboost" in training:
        _expect_dict(training["xgboost"], "stage1.training.xgboost")
    if "lightgbm" in training:
        lg = _expect_dict(training["lightgbm"], "stage1.training.lightgbm")
        _reject_unknown_keys(lg, {"n_jobs"}, "stage1.training.lightgbm")
        if "n_jobs" in lg:
            v = lg["n_jobs"]
            if not isinstance(v, int) or isinstance(v, bool):
                raise ValueError("stage1.training.lightgbm.n_jobs must be an integer (e.g. -1 or 8)")

    feat = _expect_dict(stage1.get("feature_extraction", {}), "stage1.feature_extraction")
    _reject_unknown_keys(
        feat,
        {
            "window_size",
            "density_window",
            "coverage_window",
            "soft_clip_window",
            "min_mapq",
            "splice_site_window",
            "gradient_analysis_range",
            "extended_window",
            "parallel",
            "parallel_min_sites",
            "n_workers",
            "cache_dir",
        },
        "stage1.feature_extraction",
    )
    for k in (
        "window_size",
        "density_window",
        "coverage_window",
        "soft_clip_window",
        "min_mapq",
        "splice_site_window",
        "gradient_analysis_range",
        "extended_window",
        "parallel_min_sites",
    ):
        if k in feat:
            _expect_number(feat[k], f"stage1.feature_extraction.{k}")
    if "parallel" in feat:
        _expect_bool(feat["parallel"], "stage1.feature_extraction.parallel")
    if "n_workers" in feat:
        _expect_int_or_none(feat["n_workers"], "stage1.feature_extraction.n_workers")

    features = _expect_dict(stage1.get("features", {}), "stage1.features")
    _reject_unknown_keys(features, {"k_values"}, "stage1.features")
    if "k_values" in features:
        _expect_list_int(features["k_values"], "stage1.features.k_values")


def validate_benchmark_config(cfg_map: dict[str, Any]) -> None:
    """
    Validate top-level benchmark keys ``train``, ``tests``, optional ``execution`` and ``analysis``.

    Optional root keys ``generated_by`` and ``schema_version`` (written by ``benchmark-matrix``) are
    accepted and type-checked.

    Ensures ``train.mode`` is ``run`` or ``skip`` with required fields for each mode, ``tests`` is a
    non-empty list of mappings with allowed keys only, and nested analysis flags have correct types.
    """
    _expect_dict(cfg_map, "benchmark config")
    _reject_unknown_keys(
        cfg_map,
        {"train", "tests", "execution", "analysis", "generated_by", "schema_version"},
        "root",
    )
    if "generated_by" in cfg_map and not isinstance(cfg_map["generated_by"], str):
        raise ValueError("generated_by must be a string")
    if "schema_version" in cfg_map:
        sv = cfg_map["schema_version"]
        if not isinstance(sv, int) or isinstance(sv, bool):
            raise ValueError("schema_version must be an integer")
    train = _expect_dict(cfg_map.get("train"), "train")
    _reject_unknown_keys(
        train,
        {"mode", "bam", "gtf", "gtf_pool", "ref_gtf", "tmap", "tmap_pool", "config", "outdir", "model_dir"},
        "train",
    )
    train_mode = str(train.get("mode", "run")).strip().lower()
    if train_mode == "skip":
        md = train.get("model_dir")
        if md is None or str(md).strip() == "":
            raise ValueError("train.mode=skip requires train.model_dir")
    elif train_mode == "run":
        for key in ("bam", "gtf", "ref_gtf", "tmap"):
            v = train.get(key)
            if v is None or str(v).strip() == "":
                raise ValueError(f"train.mode=run requires train.{key}")
    else:
        raise ValueError("train.mode must be 'run' or 'skip'")
    tests = cfg_map.get("tests")
    if not isinstance(tests, list) or not tests:
        raise ValueError("tests must be a non-empty list")
    for i, t in enumerate(tests, start=1):
        td = _expect_dict(t, f"tests[{i}]")
        _reject_unknown_keys(
            td,
            {"id", "assembler_id", "bam", "gtf", "ref_gtf", "tmap", "config", "outdir", "model_dir"},
            f"tests[{i}]",
        )
    if "execution" in cfg_map:
        ex = _expect_dict(cfg_map["execution"], "execution")
        _reject_unknown_keys(ex, {"stop_on_error", "max_parallel_tests"}, "execution")
        if "stop_on_error" in ex:
            _expect_bool(ex["stop_on_error"], "execution.stop_on_error")
        if "max_parallel_tests" in ex:
            mpt = ex["max_parallel_tests"]
            if not isinstance(mpt, int) or isinstance(mpt, bool) or mpt < 1:
                raise ValueError("execution.max_parallel_tests must be a positive integer")
    if "analysis" in cfg_map:
        an = _expect_dict(cfg_map["analysis"], "analysis")
        _reject_unknown_keys(
            an, {"enabled", "benchmark_mode", "debug", "pr_vs_baseline"}, "analysis"
        )
        if "enabled" in an:
            _expect_bool(an["enabled"], "analysis.enabled")
        if "benchmark_mode" in an:
            bm = str(an["benchmark_mode"]).strip().lower()
            if bm not in ("minimal", "full"):
                raise ValueError("analysis.benchmark_mode must be 'minimal' or 'full'")
        if an.get("debug") is not None:
            dbg = _expect_dict(an["debug"], "analysis.debug")
            _reject_unknown_keys(dbg, {"keep_pr_work"}, "analysis.debug")
            if "keep_pr_work" in dbg:
                _expect_bool(dbg["keep_pr_work"], "analysis.debug.keep_pr_work")
        if "pr_vs_baseline" in an:
            pr = _expect_dict(an["pr_vs_baseline"], "analysis.pr_vs_baseline")
            _reject_unknown_keys(
                pr,
                {
                    "enabled",
                    "measure",
                    "plot",
                    "filter_validation_chroms",
                    "chromosomes_file",
                    "save_pr_tables",
                    "gffcompare_bin",
                },
                "analysis.pr_vs_baseline",
            )
