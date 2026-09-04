"""End-to-end scoring on a synthetic prediction built from a reference structure.

A fake Boltz-2 output (model = the reference itself with pLDDT-like B-factors and a
confidence JSON) and a fake ColabFold output (same coordinates, PDB, scores JSON)
are written into a temporary results tree; ``cofoldbench.score`` must produce
results.csv with DockQ = 1.0 for both and ``cofoldbench.plot`` must write figures.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from cofoldbench import plot, score
from cofoldbench import structure as st


def _fake_outputs(cif: Path, row: pd.Series, boltz_root: Path, cf_root: Path) -> None:
    pid = cif.stem
    asm = st.read_structure(cif, assembly="1")
    prot = st.relabel_chains(st.protein_atoms(asm), {row["asm_chain_a"]: "A", row["asm_chain_b"]: "B"})
    # Boltz-2 (2.2.1 writer): mmCIF with pLDDT*100 in the B-factor column
    b = prot.copy()
    b.b_factor = np.full(b.array_length(), 90.0)
    d = boltz_root / pid / f"boltz_results_{pid}" / "predictions" / pid
    d.mkdir(parents=True)
    st.write_cif(b, d / f"{pid}_model_0.cif")
    (d / f"confidence_{pid}_model_0.json").write_text(json.dumps({
        "confidence_score": 0.88, "ptm": 0.9, "iptm": 0.85, "protein_iptm": 0.85, "complex_plddt": 0.9,
        "complex_iplddt": 0.88, "complex_pde": 1.0, "complex_ipde": 1.5,
        "chains_ptm": {"0": 0.9, "1": 0.9}, "pair_chains_iptm": {"0": {"0": 0.9, "1": 0.85}, "1": {"0": 0.85, "1": 0.9}},
    }))
    # ColabFold: PDB with pLDDT 0-100 in B-factors
    c = prot.copy()
    c.b_factor = np.full(c.array_length(), 90.0)
    d2 = cf_root / pid
    d2.mkdir(parents=True)
    tag = "rank_001_alphafold2_multimer_v3_model_3_seed_000"
    st.write_pdb(c, d2 / f"{pid}_unrelaxed_{tag}.pdb")
    n = st.residue_count(c, "A") + st.residue_count(c, "B")
    (d2 / f"{pid}_scores_{tag}.json").write_text(json.dumps({"plddt": [90.0] * n, "max_pae": 5.0, "pae": [], "ptm": 0.9, "iptm": 0.8}))
    (d2 / f"{pid}.done.txt").write_text("")


def test_score_and_plot_on_synthetic_predictions(tmp_path, root, targets, reference_cifs):
    cfg = yaml.safe_load((root / "config.yaml").read_text())
    boltz_root, cf_root = tmp_path / "boltz", tmp_path / "colabfold"
    ids = [c.stem for c in reference_cifs[:3]]
    for cif in reference_cifs[:3]:
        _fake_outputs(cif, targets[targets["pdb_id"] == cif.stem].iloc[0], boltz_root, cf_root)
    cfg["paths"].update({
        "boltz_out": str(boltz_root), "colabfold_out": str(cf_root),
        "results_dir": str(tmp_path / "results"), "figures_dir": str(tmp_path / "figures"),
        "targets_dir": str(root / "data/targets"), "cif_dir": str(root / "data/targets/cif"),
    })
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    assert score.main(["--config", str(cfg_path), "--ids", *ids]) == 0
    res = pd.read_csv(tmp_path / "results" / "results.csv")
    assert len(res) == 2 * len(ids)
    assert set(res["predictor"]) == {"boltz2", "af2m"}
    assert (res["dockq_error"].fillna("") == "").all()
    assert all(math.isclose(v, 1.0, abs_tol=1e-3) for v in res["dockq"])
    assert res["capri"].eq("high").all()
    # B-factor scale handling: both predictors must yield the same interface pLDDT (90) and contacts
    piv = res.pivot(index="pdb_id", columns="predictor", values=["pdockq_if_plddt", "pdockq_n_if_contacts", "pdockq"])
    assert np.allclose(piv["pdockq_if_plddt"]["boltz2"], 90.0) and np.allclose(piv["pdockq_if_plddt"]["af2m"], 90.0)
    assert (piv["pdockq_n_if_contacts"]["boltz2"] == piv["pdockq_n_if_contacts"]["af2m"]).all()
    assert (res["pdockq"] > 0.2).all()  # real interfaces with pLDDT 90 -> at least "acceptable"-range pDockQ
    assert (res.loc[res.predictor == "boltz2", "iptm"] == 0.85).all()
    assert (tmp_path / "results" / "summary.csv").exists()

    assert plot.main(["--config", str(cfg_path)]) == 0
    for name in ("dockq_scatter.png", "dockq_bars.png", "iptm_vs_dockq.png"):
        assert (tmp_path / "figures" / name).stat().st_size > 1000
