"""
Fail-fast validation of filesystem inputs before expensive BAM/GTF work.

Training and prediction call these checks after config load so missing indices or models surface
with clear :class:`PreflightError` messages. :class:`RunLayout` standardizes output subdirectory names.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


class PreflightError(RuntimeError):
    """Raised when a required path, index, or model artifact is missing or invalid."""


@dataclass(frozen=True)
class RunLayout:
    """Absolute paths for a standard Telos run directory tree under ``root``."""

    root: Path
    models_dir: Path
    predictions_dir: Path
    reports_dir: Path


def ensure_run_layout(
    outdir: Path,
    create_aux_dirs: bool = True,
    *,
    create_models_dir: bool = True,
) -> RunLayout:
    """
    Create output subdirectories and return their paths.

    Always creates ``predictions/``. Creates ``models/`` when ``create_models_dir`` is true
    (training writes models there; prediction reads models from ``--model-dir`` and does not
    need a local ``models/``). If ``create_aux_dirs``, also creates ``reports/``.

    Returns:
        :class:`RunLayout` with resolved absolute paths.
    """
    root = outdir.resolve()
    models_dir = root / "models"
    predictions_dir = root / "predictions"
    reports_dir = root / "reports"

    predictions_dir.mkdir(parents=True, exist_ok=True)
    if create_models_dir:
        models_dir.mkdir(parents=True, exist_ok=True)
    if create_aux_dirs:
        reports_dir.mkdir(parents=True, exist_ok=True)
    return RunLayout(
        root=root,
        models_dir=models_dir,
        predictions_dir=predictions_dir,
        reports_dir=reports_dir,
    )


def _validate_file_exists(path: Path, label: str) -> None:
    """Require ``path`` to exist and be a regular file; else raise :class:`PreflightError`."""
    if not path.exists():
        raise PreflightError(f"{label} not found: {path}")
    if not path.is_file():
        raise PreflightError(f"{label} must be a file: {path}")


def validate_input_file(path: Path, label: str) -> None:
    """Public wrapper to assert an arbitrary input file exists (same rules as ``_validate_file_exists``)."""
    _validate_file_exists(path, label)


def validate_input_gtf(gtf_path: Path) -> None:
    """Validate ``gtf_path`` using :func:`_validate_gtf_format` (existence + minimal GTF structure)."""
    _validate_gtf_format(gtf_path)


def _validate_gtf_format(gtf_path: Path) -> None:
    """
    Lightweight GTF sanity check without parsing attributes:

    - Path exists and is a file.
    - At least one non-empty, non-``#`` line exists.
    - The first such line splits to at least nine tab columns (GTF9 convention).

    Raises:
        PreflightError: On any violation (empty file, malformed first record, etc.).
    """
    _validate_file_exists(gtf_path, "GTF")
    non_comment_seen = False
    with gtf_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            non_comment_seen = True
            cols = line.split("\t")
            if len(cols) < 9:
                raise PreflightError(
                    f"GTF appears malformed (expected >=9 columns): {gtf_path}"
                )
            break
    if not non_comment_seen:
        raise PreflightError(f"GTF has no feature records: {gtf_path}")


def _candidate_bam_indices(bam_path: Path) -> Iterable[Path]:
    """
    Yield likely BAM index paths (``.bai`` / ``.csi``) for both ``name.bam.bai`` and ``name.bai`` styles.
    """
    yield bam_path.with_suffix(bam_path.suffix + ".bai")  # sample.bam.bai
    yield bam_path.with_suffix(".bai")  # sample.bai
    yield bam_path.with_suffix(bam_path.suffix + ".csi")  # sample.bam.csi
    yield bam_path.with_suffix(".csi")  # sample.csi


def _validate_bam_and_index(bam_path: Path) -> None:
    """
    Assert BAM exists, suffix is ``.bam``, a sidecar index exists, and pysam reports coordinate sort.

    Opens the BAM with pysam and reads header ``HD:SO``; must be ``coordinate``. Raises
    :class:`PreflightError` if pysam is missing, the file is not a BAM, or sort order is wrong.
    """
    _validate_file_exists(bam_path, "BAM")
    if bam_path.suffix.lower() != ".bam":
        raise PreflightError(f"Expected a .bam input: {bam_path}")

    has_index = any(idx.exists() and idx.is_file() for idx in _candidate_bam_indices(bam_path))
    if not has_index:
        raise PreflightError(
            f"BAM index not found for {bam_path}. Expected one of .bai/.csi variants."
        )

    # Try to verify coordinate sort order via pysam header.
    try:
        import pysam  # type: ignore
    except Exception as exc:  # pragma: no cover - env dependent
        raise PreflightError(
            "pysam is required for BAM validation but is not available."
        ) from exc

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            hd = dict(bam.header).get("HD", {})
            so = hd.get("SO", "")
            if str(so).lower() != "coordinate":
                raise PreflightError(
                    f"BAM must be coordinate-sorted (HD:SO=coordinate), got SO={so!r} for {bam_path}"
                )
    except PreflightError:
        raise
    except Exception as exc:
        raise PreflightError(f"Unable to open/validate BAM: {bam_path}") from exc


def run_preflight_train(
    bam: Path,
    gtf: Path,
    ref_gtf: Path,
    tmap: Path | None = None,
    cfg_map: dict | None = None,
) -> None:
    """
    Preflight for training: BAM+index, assembly GTF, reference GTF, non-optional tmap file.

    ``cfg_map`` is accepted for forward compatibility; currently unused. Raises :class:`PreflightError`
    on first failure.
    """
    _validate_bam_and_index(bam)
    _validate_gtf_format(gtf)
    _validate_gtf_format(ref_gtf)
    if tmap is None:
        raise PreflightError(
            "--tmap is required for training (Stage II transcript labels from gffcompare .tmap)."
        )
    _validate_file_exists(tmap, "TMAP")
    _ = cfg_map


def run_preflight_predict(
    bam: Path,
    gtf: Path,
    model_dir: Path,
    cfg_map: dict | None = None,
) -> None:
    """
    Preflight for prediction: BAM+index, GTF, and every Stage I/II artifact file under ``model_dir``.

    Checks both RF and XGB Stage I bundles for TSS/TES and both Stage II ``.joblib`` plus feature-json
    sidecars. ``cfg_map`` reserved for future use.
    """
    from telos.models import (
        STAGE1_BACKENDS,
        stage1_bundle_path,
        stage2_feature_names_json_for_backend,
        stage2_model_joblib_for_backend,
    )

    _validate_bam_and_index(bam)
    _validate_gtf_format(gtf)
    if not model_dir.exists() or not model_dir.is_dir():
        raise PreflightError(f"Model directory not found: {model_dir}")
    for site_type in ("TSS", "TES"):
        for backend in STAGE1_BACKENDS:
            name = stage1_bundle_path(site_type, backend)
            p = model_dir / name
            if not p.is_file():
                raise PreflightError(f"Missing Stage I model artifact: {p}")
    for backend in STAGE1_BACKENDS:
        m = model_dir / stage2_model_joblib_for_backend(backend)
        if not m.is_file():
            raise PreflightError(f"Missing Stage II model artifact: {m}")
        j = model_dir / stage2_feature_names_json_for_backend(backend)
        if not j.is_file():
            raise PreflightError(f"Missing Stage II feature list: {j}")
    _ = cfg_map


def run_preflight_benchmark(config: Path) -> None:
    """
    Minimal benchmark config check: path exists, non-empty; if content looks like JSON, parse it.

    Full structural validation happens later via :func:`telos.config_validation.validate_benchmark_config`.
    """
    _validate_file_exists(config, "Benchmark config")
    txt = config.read_text(encoding="utf-8", errors="replace").strip()
    if not txt:
        raise PreflightError(f"Benchmark config is empty: {config}")
    if txt.startswith("{"):
        try:
            json.loads(txt)
        except Exception as exc:
            raise PreflightError(
                f"Benchmark config appears to be JSON but is invalid: {config}"
            ) from exc
