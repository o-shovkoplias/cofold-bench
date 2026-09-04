"""Shared fixtures: repository paths and the first available reference structure."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def targets() -> pd.DataFrame:
    csv = ROOT / "data" / "targets" / "targets.csv"
    if not csv.exists():
        pytest.skip("targets.csv not present (run select stage)")
    return pd.read_csv(csv, dtype=str)


@pytest.fixture(scope="session")
def reference_cifs(targets: pd.DataFrame) -> list[Path]:
    cifs = [ROOT / "data" / "targets" / "cif" / f"{pid}.cif" for pid in targets["pdb_id"]]
    cifs = [c for c in cifs if c.exists()]
    if not cifs:
        pytest.skip("no reference mmCIF files (run download stage)")
    return cifs
