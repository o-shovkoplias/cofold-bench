"""Configuration loading with repo-relative path resolution."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config.yaml"


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load ``config.yaml`` (default: the one at the repository root)."""
    cfg_path = Path(path) if path else DEFAULT_CONFIG
    with open(cfg_path, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    cfg["_root"] = str(cfg_path.resolve().parent)
    return cfg


def resolve(cfg: dict[str, Any], key: str) -> Path:
    """Return ``cfg['paths'][key]`` as an absolute path under the repo root."""
    root = Path(cfg.get("_root", REPO_ROOT))
    p = Path(cfg["paths"][key])
    return p if p.is_absolute() else root / p
