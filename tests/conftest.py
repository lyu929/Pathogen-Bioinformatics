from __future__ import annotations

import gzip
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


@pytest.fixture
def data_dir() -> Path:
    return DATA


@pytest.fixture
def write_text(tmp_path: Path):
    """Write ``text`` to ``tmp_path/name`` (gzip-compressed if the name ends in .gz)."""

    def _write(name: str, text: str) -> Path:
        path = tmp_path / name
        if name.endswith(".gz"):
            with gzip.open(path, "wt") as fh:
                fh.write(text)
        else:
            path.write_text(text)
        return path

    return _write
