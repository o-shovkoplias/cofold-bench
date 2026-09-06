"""Score predictions against the experimental references and aggregate results.

For each target and predictor the *top-ranked* model (and, optionally, every
model) is scored with

* ipTM / pTM / mean pLDDT from the predictor's confidence JSON (:mod:`confidence`),
* pDockQ (:mod:`pdockq`) from the model's CB coordinates and B-factor pLDDT
  (0-100 for both predictors; a 0-1 scale is auto-detected and rescaled),
* DockQ / iRMSD / LRMSD / fnat (:mod:`dockq_runner`) against biological
  assembly 1 of the RCSB mmCIF, restricted to protein heavy atoms and written
  as PDB with the original chain ids.

Output: ``results/results.csv`` (one row per target x predictor x model) and
``results/summary.csv`` (top-1 per target x predictor, wide format).

Usage
-----
    python -m cofoldbench.score [--config config.yaml] [--all-models] [--ids ...]
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from . import structure as st
from .confidence import Confidence, boltz_models, colabfold_models
from .config import load_config, resolve
from .dockq_runner import dockq_two_chain
from .download import load_targets
from .pdockq import PDockQParams, pdockq_from_structure

log = logging.getLogger("cofoldbench.score")


def reference_to_pdb(cif: str | Path, out: str | Path, assembly: str = "1") -> Path:
    """Write protein heavy atoms of biological assembly ``assembly`` as PDB (for DockQ)."""
    asm = st.read_structure(cif, assembly=assembly)
    prot = st.protein_atoms(asm)
    st.write_pdb(prot, out)
    return Path(out)


def find_models(cfg: dict, predictor: str, pid: str) -> list[Confidence]:
    """Locate the models of ``predictor`` for target ``pid`` under results/raw."""
    if predictor == "boltz2":
        d = resolve(cfg, "boltz_out") / pid / "boltz_results_" / "predictions" / pid
        if not d.exists():  # boltz writes out_dir/boltz_results_<input_stem>/predictions/<name>
            cands = list((resolve(cfg, "boltz_out") / pid).glob("boltz_results_*/predictions/*"))
            d = cands[0] if cands else d
        return boltz_models(d) if d.exists() else []
    if predictor == "af2m":
        d = resolve(cfg, "colabfold_out") / pid
        return colabfold_models(d, jobname=pid) if d.exists() else []
    raise ValueError(predictor)


def score_model(
    model: Confidence, native_pdb: Path, native_chains: tuple[str, str], pdq: PDockQParams
) -> dict:
    """pDockQ + DockQ for a single model file (chains A/B)."""
    atoms = st.read_structure(model.model_file)
    chains = st.protein_chain_ids(atoms)
    if len(chains) != 2:
        raise ValueError(f"{model.model_file}: expected 2 protein chains, found {chains}")
    # both Boltz-2 (>= 2.x) and ColabFold write pLDDT as 0-100 in the B-factor column; auto-detected anyway
    pdq_out = pdockq_from_structure(atoms, chains[0], chains[1], pdq, plddt_scale="auto")
    dq = dockq_two_chain(model.model_file, native_pdb, native_chains, model_chains=(chains[0], chains[1]))
    row = model.as_row()
    row.update({f"pdockq_{k}" if k != "pdockq" else "pdockq": v for k, v in pdq_out.items()})
    row.update({
        "dockq": dq.get("DockQ"), "irmsd": dq.get("iRMSD"), "lrmsd": dq.get("LRMSD"),
        "fnat": dq.get("fnat"), "fnonnat": dq.get("fnonnat"), "capri": dq.get("capri"),
        "chain_mapping": str(dq.get("mapping")), "dockq_error": dq.get("error", ""),
    })
    return row


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--ids", nargs="*")
    ap.add_argument("--all-models", action="store_true", help="score every model, not only rank 1")
    ap.add_argument("--predictors", nargs="*", default=["boltz2", "af2m"])
    ap.add_argument("--boltz-out", help="override paths.boltz_out (e.g. results/raw/boltz_default for the MSA-depth side comparison)")
    ap.add_argument("--tag", help="write results_<tag>.csv / summary_<tag>.csv instead of results.csv / summary.csv")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    cfg = load_config(args.config)

    if args.boltz_out:

        cfg["paths"]["boltz_out"] = args.boltz_out
    pdq = PDockQParams(**cfg["scoring"]["pdockq"])
    df = load_targets(cfg)
    if args.ids:
        df = df[df["pdb_id"].isin([x.upper() for x in args.ids])]
    cif_dir = resolve(cfg, "cif_dir")
    res_dir = resolve(cfg, "results_dir")
    native_dir = res_dir / "raw" / "native_pdb"
    native_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    n_missing = 0
    for _, t in df.iterrows():
        pid = t["pdb_id"]
        cif = cif_dir / f"{pid}.cif"
        if not cif.exists():
            log.warning("%s: reference mmCIF missing", pid)
            continue
        native_pdb = native_dir / f"{pid}.pdb"
        if not native_pdb.exists():
            reference_to_pdb(cif, native_pdb)
        native_chains = (t["asm_chain_a"], t["asm_chain_b"])
        for pred in args.predictors:
            models = find_models(cfg, pred, pid)
            if not models:
                n_missing += 1
                log.info("%s/%s: no predictions found (pending phase 2)", pid, pred)
                continue
            for m in models if args.all_models else models[:1]:
                try:
                    row = score_model(m, native_pdb, native_chains, pdq)
                except Exception as exc:  # noqa: BLE001
                    log.error("%s/%s rank %d failed: %s", pid, pred, m.rank, exc)
                    row = m.as_row() | {"dockq_error": f"{type(exc).__name__}: {exc}"}
                row.update({"pdb_id": pid, "role": t["role"], "subset": t["subset"], "total_length": t["total_length"]})
                rows.append(row)
                log.info("%s %-6s rank%d ipTM=%.2f pDockQ=%s DockQ=%s", pid, pred, m.rank, m.iptm,
                         f"{row.get('pdockq', float('nan')):.3f}", f"{row.get('dockq', float('nan')):.3f}")

    res_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.tag}" if args.tag else ""
    out = res_dir / f"results{suffix}.csv"
    if not rows:
        log.warning("no predictions scored; %d target/predictor pairs pending. Nothing written.", n_missing)
        return 0
    res = pd.DataFrame(rows)
    lead = ["pdb_id", "role", "subset", "total_length", "predictor", "rank", "iptm", "ptm", "plddt", "ranking_score",
            "pdockq", "pdockq_n_if_contacts", "pdockq_if_plddt", "dockq", "irmsd", "lrmsd", "fnat", "fnonnat", "capri"]
    res = res[[c for c in lead if c in res.columns] + [c for c in res.columns if c not in lead]]
    res.to_csv(out, index=False)
    top = res[res["rank"] == 1].pivot_table(index="pdb_id", columns="predictor", values=["iptm", "pdockq", "dockq"])
    top.columns = [f"{a}_{b}" for a, b in top.columns]
    top.to_csv(res_dir / f"summary{suffix}.csv")
    log.info("wrote %s (%d rows) and summary.csv; %d target/predictor pairs still pending", out, len(res), n_missing)
    return 0


if __name__ == "__main__":
    sys.exit(main())
