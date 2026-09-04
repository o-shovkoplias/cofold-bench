"""DockQ of every reference against itself must be 1.0 (chain-mapping and I/O sanity)."""

import math

import pytest

from cofoldbench import structure as st
from cofoldbench.dockq_runner import dockq_two_chain
from cofoldbench.score import reference_to_pdb


@pytest.mark.parametrize("idx", [0, 1, 2])
def test_dockq_reference_vs_itself(idx, reference_cifs, targets, tmp_path):
    if idx >= len(reference_cifs):
        pytest.skip("fewer references")
    cif = reference_cifs[idx]
    row = targets[targets["pdb_id"] == cif.stem].iloc[0]
    native = reference_to_pdb(cif, tmp_path / f"{cif.stem}_native.pdb")
    # "model" = the same coordinates, chains renamed to A/B like a prediction
    asm = st.read_structure(cif, assembly="1")
    prot = st.protein_atoms(asm)
    model = tmp_path / f"{cif.stem}_model.pdb"
    st.write_pdb(st.relabel_chains(prot, {row["asm_chain_a"]: "A", row["asm_chain_b"]: "B"}), model)
    res = dockq_two_chain(model, native, native_chains=(row["asm_chain_a"], row["asm_chain_b"]))
    assert res["error"] == ""
    assert math.isclose(res["DockQ"], 1.0, abs_tol=1e-3), res
    assert res["fnat"] > 0.999 and res["iRMSD"] < 1e-3 and res["LRMSD"] < 1e-3
    assert res["capri"] == "high"


def test_dockq_swapped_model_chains_is_recovered(reference_cifs, targets, tmp_path):
    """If the model has chains A/B swapped relative to the native, mapping search must still give 1.0."""
    cif = reference_cifs[0]
    row = targets[targets["pdb_id"] == cif.stem].iloc[0]
    native = reference_to_pdb(cif, tmp_path / "native.pdb")
    asm = st.read_structure(cif, assembly="1")
    prot = st.protein_atoms(asm)
    model = tmp_path / "model_swapped.pdb"
    st.write_pdb(st.relabel_chains(prot, {row["asm_chain_a"]: "B", row["asm_chain_b"]: "A"}), model)
    res = dockq_two_chain(model, native, native_chains=(row["asm_chain_a"], row["asm_chain_b"]))
    assert math.isclose(res["DockQ"], 1.0, abs_tol=1e-3)
    assert res["mapping"] == {row["asm_chain_a"]: "B", row["asm_chain_b"]: "A"}
