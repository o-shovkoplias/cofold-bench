"""Biotite structure utilities tested on the downloaded references."""

import numpy as np

from cofoldbench import structure as st
from cofoldbench.pdockq import pdockq_from_structure


def test_reference_assembly_is_heterodimer(reference_cifs, targets):
    for cif in reference_cifs[:4]:
        row = targets[targets["pdb_id"] == cif.stem].iloc[0]
        asm = st.read_structure(cif, assembly="1")
        chains = st.protein_chain_ids(asm)
        assert len(chains) == 2, (cif.stem, chains)
        assert set(chains) == {row["asm_chain_a"], row["asm_chain_b"]}
        n = st.interface_residue_pairs(asm, chains[0], chains[1], cutoff=5.0)
        assert n >= 20
        # observed sequence is a substring-compatible subset of the canonical one
        seq = st.sequence_of_chain(asm, chains[0])
        assert 0 < len(seq) <= len(row["seq_a"]) + 5


def test_cb_coordinates_shapes(reference_cifs):
    asm = st.read_structure(reference_cifs[0], assembly="1")
    a, b = st.protein_chain_ids(asm)
    coords, res_ids, bfac = st.cb_coordinates(asm, a)
    assert coords.shape == (len(res_ids), 3)
    assert coords.shape[0] == st.residue_count(asm, a)
    assert np.isfinite(coords).all()


def test_pdockq_on_reference_with_fake_plddt(reference_cifs):
    """Experimental B-factors are not pLDDT; the call must still run and return a contact count."""
    asm = st.read_structure(reference_cifs[0], assembly="1")
    a, b = st.protein_chain_ids(asm)
    out = pdockq_from_structure(asm, a, b)
    assert out["n_if_contacts"] > 0
    assert 0.0 <= out["pdockq"] <= 0.75


def test_relabel_and_write_roundtrip(tmp_path, reference_cifs):
    asm = st.read_structure(reference_cifs[0], assembly="1")
    a, b = st.protein_chain_ids(asm)
    prot = st.protein_atoms(asm)
    renamed = st.relabel_chains(prot, {a: "X", b: "Y"})
    assert set(st.protein_chain_ids(renamed)) == {"X", "Y"}
    out = tmp_path / "ref.pdb"
    st.write_pdb(renamed, out)
    back = st.read_structure(out)
    assert set(st.protein_chain_ids(back)) == {"X", "Y"}
    assert st.residue_count(back, "X") == st.residue_count(renamed, "X")
