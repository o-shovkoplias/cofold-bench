"""Confidence parsers on synthetic Boltz-2 / ColabFold output directories."""

import json

from cofoldbench.confidence import boltz_models, colabfold_models


def test_boltz_parser(tmp_path):
    d = tmp_path / "predictions" / "9ASS"
    d.mkdir(parents=True)
    for k, cs in enumerate([0.7, 0.9]):
        (d / f"confidence_9ASS_model_{k}.json").write_text(json.dumps({
            "confidence_score": cs, "ptm": 0.8, "iptm": 0.6 + k / 10, "protein_iptm": 0.6 + k / 10,
            "complex_plddt": 0.85, "complex_iplddt": 0.8, "complex_pde": 1.0, "complex_ipde": 2.0,
            "chains_ptm": {"0": 0.9, "1": 0.8}, "pair_chains_iptm": {"0": {"0": 0.9, "1": 0.6}, "1": {"0": 0.6, "1": 0.8}},
        }))
        (d / f"9ASS_model_{k}.cif").write_text("data_x\n")
    models = boltz_models(d)
    assert [m.rank for m in models] == [1, 2]
    assert models[0].model_file.endswith("9ASS_model_1.cif")  # higher confidence_score first
    assert models[0].iptm == 0.7 and models[0].plddt == 85.0


def test_colabfold_parser(tmp_path):
    d = tmp_path
    for rank, (iptm, ptm) in enumerate([(0.8, 0.7), (0.5, 0.6)], start=1):
        tag = f"rank_{rank:03d}_alphafold2_multimer_v3_model_{rank}_seed_000"
        (d / f"9ASS_scores_{tag}.json").write_text(json.dumps({
            "plddt": [90.0, 80.0], "max_pae": 10.0, "pae": [[0, 1], [1, 0]], "ptm": ptm, "iptm": iptm,
        }))
        (d / f"9ASS_unrelaxed_{tag}.pdb").write_text("END\n")
    models = colabfold_models(d, jobname="9ASS")
    assert [m.rank for m in models] == [1, 2]
    assert models[0].iptm == 0.8 and models[0].plddt == 85.0
    assert abs(models[0].ranking_score - (0.8 * 0.8 + 0.2 * 0.7)) < 1e-9
    assert models[0].model_file.endswith("model_1_seed_000.pdb")
