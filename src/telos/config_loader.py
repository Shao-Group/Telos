"""
Load and navigate Telos YAML/JSON configuration files.

Paths are resolved relative to the caller; the default Stage I file is package data under
``telos/configs/``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def default_stage1_config_path() -> Path:
    """
    Return the absolute path to packaged ``telos/configs/stage1.defaults.yaml``.

    Used when CLI ``--config`` is omitted.
    """
    return Path(__file__).resolve().parent / "configs" / "stage1.defaults.yaml"


def load_mapping_config(path: Path | None) -> dict[str, Any]:
    """
    Read a mapping (dict) from a YAML or JSON config file.

    Behavior:

    - ``path is None`` → return ``{}`` (predict may omit config in some flows; caller validates).
    - ``.json`` → parse with stdlib ``json``; top level must be an object.
    - ``.yaml`` / ``.yml`` → parse with PyYAML ``safe_load``; top level must be a mapping.
    - Other suffixes → :class:`ValueError`.

    Raises:
        ValueError: Missing file, parse error, wrong top-level type, or missing PyYAML for YAML.
    """
    if path is None:
        return {}
    if not path.exists() or not path.is_file():
        raise ValueError(f"Config file not found: {path}")

    text = path.read_text(encoding="utf-8", errors="replace")
    suffix = path.suffix.lower()
    if suffix == ".json":
        import json

        data = json.loads(text)
        if not isinstance(data, dict):
            raise ValueError("Top-level config must be a mapping/object.")
        return data

    if suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore
        except ModuleNotFoundError as exc:
            raise ValueError(
                "YAML config requested but PyYAML is not installed. "
                "Install pyyaml or use a JSON config."
            ) from exc
        try:
            data = yaml.safe_load(text) or {}
        except Exception as exc:
            raise ValueError(f"Invalid YAML in {path}: {exc}") from exc
        if not isinstance(data, dict):
            raise ValueError("Top-level config must be a mapping/object.")
        return data

    raise ValueError(f"Unsupported config extension: {path.suffix}")


def get_nested(cfg: dict[str, Any], keys: list[str], default: Any) -> Any:
    """
    Walk a nested dict using a key list; return ``default`` if any step is missing or not a dict.

    Example: ``get_nested(cfg, ["stage1", "training", "split_policy"], "chr1-10")`` returns the
    string policy or ``"chr1-10"`` when the chain does not exist. Does **not** copy; returns the
    actual nested value reference when present.
    """
    cur: Any = cfg
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur
