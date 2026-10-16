"""Load and validate the TOML pipeline configuration."""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

REQUIRED_SECTIONS = (
    "output",
    "thresholds",
    "qc",
    "assembly",
    "annotation",
    "pangenome",
    "rnaseq",
    "structure",
)


class ConfigError(ValueError):
    pass


def load_config(path: str | Path, root: str | Path | None = None) -> dict[str, Any]:
    """Read the config and resolve every ``*_dir``/path-like string against ``root``.

    ``root`` defaults to the parent of the ``config/`` directory (the repository root).
    Missing input files are reported together, before any analysis starts.
    """
    path = Path(path)
    cfg = tomllib.loads(path.read_text(encoding="utf-8"))
    missing_sections = [s for s in REQUIRED_SECTIONS if s not in cfg]
    if missing_sections:
        raise ConfigError(f"{path}: missing sections {missing_sections}")
    base = Path(root) if root is not None else path.resolve().parent.parent
    cfg["_root"] = base

    missing: list[str] = []

    def resolve(value: Any) -> Any:
        if isinstance(value, str) and value.startswith("data/"):
            p = base / value
            if not p.exists():
                missing.append(value)
            return p
        if isinstance(value, list):
            return [resolve(v) for v in value]
        if isinstance(value, dict):
            return {k: resolve(v) for k, v in value.items()}
        return value

    for key in REQUIRED_SECTIONS:
        if key == "output":
            continue
        cfg[key] = resolve(cfg[key])
    cfg["output"] = {k: base / v for k, v in cfg["output"].items()}
    if missing:
        raise ConfigError("missing input files:\n  " + "\n  ".join(sorted(missing)))
    if not isinstance(cfg["rnaseq"], list) or not cfg["rnaseq"]:
        raise ConfigError("[[rnaseq]] must define at least one dataset")
    return cfg
