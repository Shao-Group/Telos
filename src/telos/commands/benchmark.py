"""
``benchmark`` CLI handler: re-exports :func:`~telos.benchmark.orchestrator.run_benchmark`.
"""

from __future__ import annotations

from telos.benchmark.orchestrator import run_benchmark

__all__ = ["run_benchmark"]
