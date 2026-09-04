"""pDockQ - predicted DockQ from pLDDT and interface contacts (Bryant et al. 2022).

Reference: Bryant P., Pozzati G., Elofsson A. "Improved prediction of
protein-protein interactions using AlphaFold2", Nat. Commun. 13, 1265 (2022).

    pDockQ = L / (1 + exp(-k (x - x0))) + b
    x      = <pLDDT_interface> * log10(N_contacts)

with L = 0.724, x0 = 152.611, k = 0.052, b = 0.018 and interface contacts
defined as CB-CB (CA for Gly) pairs between the two chains closer than 8 A.
If there are no contacts, pDockQ = b (= 0.018), as in the original script.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import structure as st


@dataclass(frozen=True)
class PDockQParams:
    """Sigmoid parameters of pDockQ (defaults from Bryant et al. 2022)."""

    L: float = 0.724
    x0: float = 152.611
    k: float = 0.052
    b: float = 0.018
    cb_cutoff: float = 8.0


DEFAULT = PDockQParams()


def pdockq_from_x(x: float, p: PDockQParams = DEFAULT) -> float:
    """Evaluate the sigmoid at ``x = <pLDDT_if> * log10(N_contacts)``."""
    z = np.clip(-p.k * (x - p.x0), -700.0, 700.0)  # avoid overflow warnings far from x0
    return float(p.L / (1.0 + np.exp(z)) + p.b)


def interface_stats(
    coords_a: np.ndarray, plddt_a: np.ndarray, coords_b: np.ndarray, plddt_b: np.ndarray, cutoff: float = 8.0
) -> tuple[int, float]:
    """Return ``(n_contacts, mean_interface_plddt)`` for two chains.

    ``n_contacts`` counts residue pairs (i in A, j in B) with |CB_i - CB_j| < cutoff.
    The interface pLDDT is averaged over the *residues* of both chains that take
    part in at least one contact (the definition of the original pDockQ script).
    """
    if len(coords_a) == 0 or len(coords_b) == 0:
        return 0, float("nan")
    d = np.linalg.norm(coords_a[:, None, :] - coords_b[None, :, :], axis=-1)
    mask = d < cutoff
    n = int(mask.sum())
    if n == 0:
        return 0, float("nan")
    ia = np.where(mask.any(axis=1))[0]
    ib = np.where(mask.any(axis=0))[0]
    if_plddt = float(np.concatenate([plddt_a[ia], plddt_b[ib]]).mean())
    return n, if_plddt


def pdockq_from_arrays(
    coords_a: np.ndarray, plddt_a: np.ndarray, coords_b: np.ndarray, plddt_b: np.ndarray, p: PDockQParams = DEFAULT
) -> dict[str, float]:
    """pDockQ from per-residue CB coordinates and pLDDT (0-100 scale)."""
    n, if_plddt = interface_stats(coords_a, plddt_a, coords_b, plddt_b, p.cb_cutoff)
    if n == 0:
        return {"pdockq": p.b, "n_if_contacts": 0, "if_plddt": float("nan"), "x": float("nan")}
    x = if_plddt * np.log10(n)
    return {"pdockq": pdockq_from_x(x, p), "n_if_contacts": n, "if_plddt": if_plddt, "x": float(x)}


def plddt_scale_factor(bfactors: np.ndarray) -> float:
    """Return the factor that brings B-factor pLDDT onto the 0-100 scale.

    Boltz-2 (>= 2.x, ``boltz/data/write/mmcif.py``: ``plddt * 100``) and ColabFold
    both write pLDDT as 0-100; some tools write 0-1.  If every value is <= 1 we
    assume the 0-1 convention and multiply by 100.
    """
    b = bfactors[np.isfinite(bfactors)]
    if b.size and np.nanmax(b) <= 1.0:
        return 100.0
    return 1.0


def pdockq_from_structure(
    atoms, chain_a: str, chain_b: str, p: PDockQParams = DEFAULT, plddt_scale: float | str = "auto"
) -> dict[str, float]:
    """pDockQ for a predicted model whose B-factor column holds pLDDT.

    ``plddt_scale`` rescales the B-factor to the 0-100 range; ``"auto"`` (default)
    uses :func:`plddt_scale_factor`.
    """
    ca, _, ba = st.cb_coordinates(atoms, chain_a)
    cb, _, bb = st.cb_coordinates(atoms, chain_b)
    scale = plddt_scale_factor(np.concatenate([ba, bb])) if plddt_scale == "auto" else float(plddt_scale)
    out = pdockq_from_arrays(ca, ba * scale, cb, bb * scale, p)
    out["plddt_scale"] = scale
    return out
