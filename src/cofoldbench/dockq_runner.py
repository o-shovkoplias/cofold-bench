"""Thin wrapper around the ``DockQ`` package (Mirabello & Wallner 2024, v2).

DockQ is imported as a module::

    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

``run_on_all_native_interfaces(model, native, chain_map={native: model})`` returns
``({(nat_chain1, nat_chain2): {...}}, total_dockq)``.  For a heterodimer there is
one interface; we try both chain assignments (A->a,B->b and A->b,B->a) and keep
the best, which is the "automatic chain mapping" of the DockQ CLI restricted to
two chains.  Note that the module API expects the *native:model* mapping, while
the CLI ``--mapping`` flag is *model:native*.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from DockQ.DockQ import load_PDB, run_on_all_native_interfaces


def dockq_two_chain(
    model_path: str | Path,
    native_path: str | Path,
    native_chains: tuple[str, str],
    model_chains: tuple[str, str] = ("A", "B"),
    try_swapped: bool = True,
) -> dict[str, Any]:
    """DockQ of a two-chain ``model`` against a two-chain ``native``.

    Returns the interface dictionary of the best mapping plus ``mapping`` (as
    ``native:model``), ``DockQ`` and the CAPRI class.
    """
    model = load_PDB(str(model_path))
    native = load_PDB(str(native_path))
    n1, n2 = native_chains
    m1, m2 = model_chains
    mappings = [{n1: m1, n2: m2}]
    if try_swapped:
        mappings.append({n1: m2, n2: m1})
    best: dict[str, Any] | None = None
    for cm in mappings:
        try:
            res, total = run_on_all_native_interfaces(model, native, chain_map=cm)
        except Exception as exc:  # noqa: BLE001 - alignment failures for a wrong mapping
            res, total = {}, -1.0
            err = f"{type(exc).__name__}: {exc}"
        else:
            err = ""
        if not res:
            cand = {"DockQ": float("nan"), "error": err or "no_interface", "mapping": cm}
        else:
            key = next(iter(res))
            cand = dict(res[key])
            cand["mapping"] = cm
            cand["error"] = ""
        if best is None or (cand.get("DockQ") or -1) > (best.get("DockQ") or -1):
            best = cand
    assert best is not None
    best["capri"] = capri_class(best.get("DockQ", float("nan")))
    return best


def capri_class(dockq: float, acceptable: float = 0.23, medium: float = 0.49, high: float = 0.80) -> str:
    """CAPRI-style quality class from DockQ (Basu & Wallner 2016 thresholds)."""
    if dockq != dockq:  # NaN
        return "n/a"
    if dockq >= high:
        return "high"
    if dockq >= medium:
        return "medium"
    if dockq >= acceptable:
        return "acceptable"
    return "incorrect"
