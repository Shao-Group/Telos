"""Legacy ``benchmark`` CLI redirect."""

from __future__ import annotations

from telos.config_models import BenchmarkIO


def run_benchmark(cfg: BenchmarkIO) -> int:
    """Redirect callers to the separate reproduction repository."""
    _ = cfg
    print(
        "[telos] Benchmarking lives in the Telos-test repo: "
        "https://github.com/Shao-Group/Telos-test. See its README."
    )
    return 2
