"""Logging configuration for the Telos command-line interface."""

from __future__ import annotations

import logging
from argparse import Namespace


def configure_logging(args: Namespace) -> None:
    """Configure one concise stderr handler for the ``telos`` logger."""
    verbosity = int(getattr(args, "verbose", 0) or 0)
    quiet = bool(getattr(args, "quiet", False))
    level = (
        logging.ERROR
        if quiet
        else logging.DEBUG
        if verbosity >= 2
        else logging.INFO
        if verbosity == 1
        else logging.WARNING
    )

    telos_logger = logging.getLogger("telos")
    telos_logger.setLevel(level)
    telos_logger.propagate = False

    for handler in list(telos_logger.handlers):
        if getattr(handler, "_telos_cli_handler", False):
            telos_logger.removeHandler(handler)

    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("[telos] %(message)s"))
    handler._telos_cli_handler = True  # type: ignore[attr-defined]
    telos_logger.addHandler(handler)
