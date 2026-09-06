"""Thin wrapper around the ``DockQ`` package (Mirabello & Wallner 2024, v2).

DockQ is imported as a module::

    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

``run_on_all_native_interfaces(model, native, chain_map={native: model})`` returns
``({(nat_chain1, nat_chain2): {...}}, total_dockq)``.  For a heterodimer there is
one interface and two possible chain assignments (A->a, B->b or A->b, B->a).

Chain mapping.  The DockQ CLI groups model and native chains by sequence and only
permutes chains within a group, so for a heterodimer of two unrelated proteins it
considers exactly one mapping.  We do the same: the mapping is chosen by sequence
identity (Biotite local alignment of the observed chain sequences); only if the
two assignments are indistinguishable by sequence (pseudo-homodimers, which the
target selection excludes) are both scored and the higher DockQ kept.  Note that
the module API expects the *native:model* mapping, while the CLI ``--mapping``
flag is *model:native*.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import biotite.sequence as bseq
import biotite.sequence.align as balign
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

from . import structure as st

_MATRIX = balign.SubstitutionMatrix.std_protein_matrix()


def sequence_identity(seq_a: str, seq_b: str) -> float:
    """Identity of the best local alignment, normalised by the shorter sequence (0-1)."""
    if not seq_a or not seq_b:
        return 0.0
    ali = balign.align_optimal(
        bseq.ProteinSequence(seq_a), bseq.ProteinSequence(seq_b), _MATRIX,
        gap_penalty=(-10, -1), local=True, max_number=1,
    )[0]
    return float(balign.get_sequence_identity(ali, mode="shortest"))


def mapping_by_sequence(
    model_path: str | Path,
    native_path: str | Path,
    native_chains: tuple[str, str],
    model_chains: tuple[str, str],
    min_margin: float = 0.05,
) -> tuple[dict[str, str] | None, dict[str, float]]:
    """Choose the native:model chain mapping by sequence identity.

    Returns ``(mapping, identities)``; ``mapping`` is ``None`` when the straight and
    swapped assignments differ by less than ``min_margin`` in summed identity.
    """
    mod = st.read_structure(model_path)
    nat = st.read_structure(native_path)
    n1, n2 = native_chains
    m1, m2 = model_chains
    ns = {c: st.sequence_of_chain(nat, c) for c in native_chains}
    ms = {c: st.sequence_of_chain(mod, c) for c in model_chains}
    straight = sequence_identity(ns[n1], ms[m1]) + sequence_identity(ns[n2], ms[m2])
    swapped = sequence_identity(ns[n1], ms[m2]) + sequence_identity(ns[n2], ms[m1])
    ids = {"straight": straight / 2, "swapped": swapped / 2}
    if abs(straight - swapped) < min_margin:
        return None, ids
    return ({n1: m1, n2: m2} if straight > swapped else {n1: m2, n2: m1}), ids


def dockq_two_chain(
    model_path: str | Path,
    native_path: str | Path,
    native_chains: tuple[str, str],
    model_chains: tuple[str, str] = ("A", "B"),
    mapping: str = "sequence",
) -> dict[str, Any]:
    """DockQ of a two-chain ``model`` against a two-chain ``native``.

    ``mapping``: ``"sequence"`` (default) picks the chain assignment by sequence
    identity and falls back to ``"best"`` when indecisive; ``"best"`` scores both
    assignments and keeps the higher DockQ; ``"fixed"`` uses ``native_chains[i] ->
    model_chains[i]`` as given.  Returns the interface dictionary of the chosen
    mapping plus ``mapping`` (as ``native:model``), ``DockQ`` and the CAPRI class.
    """
    n1, n2 = native_chains
    m1, m2 = model_chains
    straight, swapped = {n1: m1, n2: m2}, {n1: m2, n2: m1}
    if mapping == "fixed":
        mappings = [straight]
    elif mapping == "best":
        mappings = [straight, swapped]
    elif mapping == "sequence":
        chosen, _ = mapping_by_sequence(model_path, native_path, native_chains, model_chains)
        mappings = [chosen] if chosen is not None else [straight, swapped]
    else:
        raise ValueError(f"unknown mapping mode {mapping!r}")

    model = load_PDB(str(model_path))
    native = load_PDB(str(native_path))
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
