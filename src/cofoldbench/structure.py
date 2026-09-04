"""Biotite-based structure utilities shared by selection and scoring.

All structure I/O goes through :mod:`biotite.structure.io.pdbx` (mmCIF) or
:mod:`biotite.structure.io.pdb` (PDB, for ColabFold outputs).
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

# B-factor is needed because Boltz-2 and ColabFold store per-residue pLDDT there.
EXTRA = ["b_factor"]


def read_structure(path: str | Path, assembly: str | None = None, model: int = 1) -> struc.AtomArray:
    """Read a PDB/mmCIF/BinaryCIF file into a Biotite :class:`AtomArray`.

    Parameters
    ----------
    path
        ``.cif`` / ``.mmcif`` / ``.bcif`` / ``.pdb``.
    assembly
        If given (e.g. ``"1"``), build that biological assembly from the mmCIF
        ``pdbx_struct_assembly`` records; otherwise return the asymmetric unit.
    model
        Model number (1-based) for multi-model files.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        f = pdbx.CIFFile.read(str(path))
        if assembly is not None:
            return pdbx.get_assembly(f, assembly_id=assembly, model=model, use_author_fields=False, extra_fields=EXTRA)
        return pdbx.get_structure(f, model=model, use_author_fields=False, extra_fields=EXTRA)
    if suffix == ".bcif":
        f = pdbx.BinaryCIFFile.read(str(path))
        if assembly is not None:
            return pdbx.get_assembly(f, assembly_id=assembly, model=model, use_author_fields=False, extra_fields=EXTRA)
        return pdbx.get_structure(f, model=model, use_author_fields=False, extra_fields=EXTRA)
    if suffix in {".pdb", ".ent"}:
        return pdb.PDBFile.read(str(path)).get_structure(model=model, extra_fields=EXTRA)
    raise ValueError(f"Unsupported structure format: {path}")


def protein_atoms(atoms: struc.AtomArray) -> struc.AtomArray:
    """Return only amino-acid residues, heavy atoms, first altloc (Biotite default)."""
    mask = struc.filter_amino_acids(atoms) & (atoms.element != "H")
    return atoms[mask]


def protein_chain_ids(atoms: struc.AtomArray) -> list[str]:
    """Chain IDs (in order of appearance) that contain at least one amino acid."""
    prot = protein_atoms(atoms)
    seen: dict[str, None] = {}
    for c in prot.chain_id:
        seen.setdefault(str(c), None)
    return list(seen)


def residue_count(atoms: struc.AtomArray, chain_id: str) -> int:
    """Number of amino-acid residues observed in ``chain_id``."""
    prot = protein_atoms(atoms)
    sub = prot[prot.chain_id == chain_id]
    return int(struc.get_residue_count(sub)) if sub.array_length() else 0


def cb_coordinates(atoms: struc.AtomArray, chain_id: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-residue representative coordinates for pDockQ: CB, or CA for glycine.

    Returns ``(coords[N,3], res_ids[N], b_factors[N])``.  ``b_factors`` holds the
    per-residue pLDDT for predicted models (Boltz-2 and ColabFold both write
    pLDDT into the B-factor column).
    """
    prot = protein_atoms(atoms)
    sub = prot[prot.chain_id == chain_id]
    coords, res_ids, bfac = [], [], []
    has_b = "b_factor" in sub.get_annotation_categories()
    for res in struc.residue_iter(sub):
        cb = res[res.atom_name == "CB"]
        if cb.array_length() == 0:
            cb = res[res.atom_name == "CA"]
        if cb.array_length() == 0:
            continue
        coords.append(cb.coord[0])
        res_ids.append(int(cb.res_id[0]))
        bfac.append(float(cb.b_factor[0]) if has_b else np.nan)
    return np.asarray(coords, dtype=float).reshape(-1, 3), np.asarray(res_ids), np.asarray(bfac)


def interface_residue_pairs(
    atoms: struc.AtomArray, chain_a: str, chain_b: str, cutoff: float = 5.0
) -> int:
    """Count residue pairs (a in chain_a, b in chain_b) with any heavy-atom distance < cutoff.

    Used as a sanity check that the two chains of a reference actually form an
    interface (a heterodimer, not two unrelated molecules in one asymmetric unit).
    """
    prot = protein_atoms(atoms)
    a = prot[prot.chain_id == chain_a]
    b = prot[prot.chain_id == chain_b]
    if a.array_length() == 0 or b.array_length() == 0:
        return 0
    cell = struc.CellList(b, cell_size=cutoff)
    near = cell.get_atoms(a.coord, radius=cutoff)  # (n_a, k) indices into b, -1 padded
    pairs: set[tuple[int, int]] = set()
    for i in range(near.shape[0]):
        js = near[i][near[i] >= 0]
        if js.size == 0:
            continue
        ra = int(a.res_id[i])
        for j in js:
            pairs.add((ra, int(b.res_id[j])))
    return len(pairs)


def sequence_of_chain(atoms: struc.AtomArray, chain_id: str) -> str:
    """One-letter sequence of the *observed* residues of a protein chain."""
    prot = protein_atoms(atoms)
    sub = prot[prot.chain_id == chain_id]
    if sub.array_length() == 0:
        return ""
    seq = struc.to_sequence(sub, allow_hetero=True)[0]
    return "".join(str(s) for s in seq)


def chains_summary(atoms: struc.AtomArray) -> dict[str, int]:
    """``{chain_id: n_residues}`` for protein chains."""
    return {c: residue_count(atoms, c) for c in protein_chain_ids(atoms)}


def relabel_chains(atoms: struc.AtomArray, mapping: dict[str, str]) -> struc.AtomArray:
    """Return a copy with ``chain_id`` renamed through ``mapping`` (unknown ids kept)."""
    out = atoms.copy()
    new = np.array([mapping.get(str(c), str(c)) for c in out.chain_id])
    out.chain_id = new
    return out


def write_pdb(atoms: struc.AtomArray, path: str | Path) -> None:
    """Write an AtomArray as PDB (DockQ reads PDB and mmCIF)."""
    f = pdb.PDBFile()
    f.set_structure(atoms)
    f.write(str(path))


def write_cif(atoms: struc.AtomArray, path: str | Path) -> None:
    """Write an AtomArray as mmCIF (with ``occupancy``/``b_factor`` columns, which DockQ's parser requires)."""
    out = atoms.copy()
    cats = out.get_annotation_categories()
    if "occupancy" not in cats:
        out.set_annotation("occupancy", np.ones(out.array_length(), dtype=float))
    if "b_factor" not in cats:
        out.set_annotation("b_factor", np.zeros(out.array_length(), dtype=float))
    f = pdbx.CIFFile()
    pdbx.set_structure(f, out)
    f.write(str(path))


def iter_protein_chains(atoms: struc.AtomArray) -> Iterable[tuple[str, struc.AtomArray]]:
    """Yield ``(chain_id, atoms_of_that_chain)`` for protein chains."""
    prot = protein_atoms(atoms)
    for c in protein_chain_ids(prot):
        yield c, prot[prot.chain_id == c]
