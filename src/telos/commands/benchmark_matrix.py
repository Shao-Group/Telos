"""Legacy ``benchmark-matrix`` CLI redirect."""

from __future__ import annotations

from pathlib import Path


def run_benchmark_matrix(
    *,
    data_type: str,
    train_annotation: str,
    test_annotation: str,
    outdir: Path,
    bundles_root: Path | None = None,
    stage1_config: Path | None = None,
) -> int:
    """Redirect callers to the separate reproduction repository."""
    _ = (
        data_type,
        train_annotation,
        test_annotation,
        outdir,
        bundles_root,
        stage1_config,
    )
    print(
        "[telos] Benchmarking lives in the Telos-test repo: "
        "https://github.com/Shao-Group/Telos-test. See its README."
    )
    return 2
