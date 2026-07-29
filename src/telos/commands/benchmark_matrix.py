"""
``benchmark-matrix`` CLI handler: re-exports :func:`~telos.benchmark.matrix.run_benchmark_matrix`.

Logic stays in :mod:`telos.benchmark.matrix` so tests import the generator without pulling in CLI.
"""

from __future__ import annotations

from telos.benchmark.matrix import run_benchmark_matrix

__all__ = ["run_benchmark_matrix"]
