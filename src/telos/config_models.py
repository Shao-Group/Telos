"""
Frozen dataclasses describing inputs for train, predict, and benchmark runs.

These types are the boundary between :mod:`telos.cli` (argument parsing) and
:mod:`telos.commands` / :mod:`telos.benchmark` (execution). All paths are
:class:`~pathlib.Path` instances; callers are responsible for resolving them if needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, kw_only=True)
class RunIO:
    """Shared fields for any run that scores sites and writes under ``outdir``."""

    bam: Path
    gtf: Path
    outdir: Path
    config_file: Path | None = None
    save_intermediates: bool = False


@dataclass(frozen=True, kw_only=True)
class TrainIO(RunIO):
    """
    Inputs for :func:`~telos.commands.train.run_train`.

    ``ref_gtf`` supplies reference TSS/TES for Stage I labels. ``tmap`` is the gffcompare
    transcript mapping for Stage II supervision. ``stage1_*`` flags override YAML parallelism
    for BAM feature extraction only.
    """

    ref_gtf: Path
    tmap: Path | None = None
    gtf_pool: tuple[Path, ...] | None = None
    tmap_pool: tuple[Path, ...] | None = None
    stage1_no_parallel: bool = False
    stage1_n_workers: int | None = None


@dataclass(frozen=True, kw_only=True)
class PredictIO(RunIO):
    """
    Inputs for :func:`~telos.commands.predict.run_predict`.

    ``model_dir`` must contain Stage I bundles and Stage II joblib+json artifacts for both backends.
    """

    model_dir: Path
    backend: str = "xgb"
    min_score: float | None = None
    stage1_no_parallel: bool = False
    stage1_n_workers: int | None = None


@dataclass(frozen=True)
class BenchmarkIO:
    """
    Inputs for :func:`~telos.benchmark.orchestrator.run_benchmark`.

    ``config`` is the benchmark YAML/JSON path. ``outdir`` is the root for ``train/``, ``tests/``,
    and ``reports/`` unless overridden inside the benchmark file.
    """

    config: Path
    outdir: Path
