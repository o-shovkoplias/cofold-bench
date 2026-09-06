"""Summarise ``results/results.csv`` into ``results/report.md`` (tables quoted in the README).

Everything printed here is recomputed from the CSVs; the README copies these tables verbatim.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .config import load_config, resolve

LABEL = {"boltz2": "Boltz-2", "af2m": "AF2-Multimer v3 (ColabFold)"}
CAPRI = [("incorrect", 0.0), ("acceptable", 0.23), ("medium", 0.49), ("high", 0.80)]


def capri_counts(x: pd.Series) -> dict:
    return {
        "acceptable+ (>=0.23)": int((x >= 0.23).sum()),
        "medium+ (>=0.49)": int((x >= 0.49).sum()),
        "high (>=0.80)": int((x >= 0.80).sum()),
        "incorrect (<0.23)": int((x < 0.23).sum()),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--driver-logs", nargs="*", default=[],
                    help="optional: driver logs (stdout of scripts/run_*.sh, not committed) with [run ]/[done] lines; adds wall-time sections")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    res_dir = resolve(cfg, "results_dir")
    d = pd.read_csv(res_dir / "results.csv")
    t = pd.read_csv(resolve(cfg, "targets_dir") / "targets.csv")
    top = d[d["rank"] == 1] if "rank" in d else d
    w = top.pivot(index="pdb_id", columns="predictor", values="dockq")
    w = w.join(top.pivot(index="pdb_id", columns="predictor", values="iptm").add_prefix("iptm_"))
    w = w.join(t.set_index("pdb_id")[["total_length", "subset", "title"]])
    lines = ["# cofold-bench report", "", f"Targets: {len(w)} heterodimers; predictors: " + ", ".join(LABEL[p] for p in top.predictor.unique()), ""]

    # 1. per-predictor summary
    lines += ["## Per-predictor summary (top-ranked model per target)", "", "| predictor | n | median DockQ | mean DockQ | acceptable+ | medium+ | high | incorrect | Spearman(ipTM, DockQ) | Spearman(pDockQ, DockQ) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for p in ["boltz2", "af2m"]:
        s = top[top.predictor == p]
        cc = capri_counts(s.dockq)
        rho_i = stats.spearmanr(s.iptm, s.dockq).statistic
        rho_p = stats.spearmanr(s.pdockq, s.dockq).statistic
        lines.append(f"| {LABEL[p]} | {len(s)} | {s.dockq.median():.3f} | {s.dockq.mean():.3f} | {cc['acceptable+ (>=0.23)']} | {cc['medium+ (>=0.49)']} | {cc['high (>=0.80)']} | {cc['incorrect (<0.23)']} | {rho_i:.2f} | {rho_p:.2f} |")
    # 2. paired comparison
    both = w.dropna(subset=["boltz2", "af2m"])
    diff = both["boltz2"] - both["af2m"]
    wil = stats.wilcoxon(both["boltz2"], both["af2m"]) if len(both) >= 6 else None
    lines += ["", "## Paired comparison (same targets)", "",
              f"- n = {len(both)}; Boltz-2 higher DockQ on {(diff > 0).sum()}, AF2-Multimer higher on {(diff < 0).sum()}, ties {(diff == 0).sum()}.",
              f"- median dDockQ (Boltz-2 - AF2-Multimer) = {diff.median():+.3f}; mean = {diff.mean():+.3f}; IQR [{diff.quantile(.25):+.3f}, {diff.quantile(.75):+.3f}].",
              (f"- Wilcoxon signed-rank p = {wil.pvalue:.3g}." if wil else "- too few pairs for a test."),
              f"- Both acceptable (>=0.23): {((both.boltz2>=0.23)&(both.af2m>=0.23)).sum()}; only Boltz-2: {((both.boltz2>=0.23)&(both.af2m<0.23)).sum()}; only AF2-Multimer: {((both.boltz2<0.23)&(both.af2m>=0.23)).sum()}; neither: {((both.boltz2<0.23)&(both.af2m<0.23)).sum()}.", ""]
    # 3. per-target table
    lines += ["## Per-target results (top-ranked model)", "", "| PDB | subset | length | DockQ Boltz-2 | DockQ AF2-M | ipTM Boltz-2 | ipTM AF2-M | title |", "|---|---|---|---|---|---|---|---|"]
    for pid, r in w.sort_values("boltz2", ascending=False).iterrows():
        lines.append(f"| {pid} | {r['subset']} | {int(r['total_length'])} | {r['boltz2']:.2f} | {r['af2m']:.2f} | {r['iptm_boltz2']:.2f} | {r['iptm_af2m']:.2f} | {str(r['title'])[:60]} |")
    # 4. confident failures
    lines += ["", "## Confident failures (ipTM >= 0.7 but DockQ < 0.23)", ""]
    cf = top[(top.iptm >= 0.7) & (top.dockq < 0.23)]
    if len(cf):
        for _, r in cf.iterrows():
            lines.append(f"- {r.pdb_id} / {LABEL[r.predictor]}: ipTM {r.iptm:.2f}, pDockQ {r.pdockq:.2f}, DockQ {r.dockq:.2f}")
    else:
        lines.append("- none")
    # 5. MSA-depth side comparison
    side = res_dir / "results_msa_default.csv"
    if side.exists():
        b = pd.read_csv(side).set_index("pdb_id")
        a = top[top.predictor == "boltz2"].set_index("pdb_id").loc[b.index]
        dd = (a.dockq - b.dockq)
        lines += ["", "## MSA-depth side comparison (Boltz-2, 7 targets that fit in 8 GB with default MSA)", "",
                  "| PDB | DockQ 512/128 | DockQ 8192/1024 | dDockQ | ipTM 512/128 | ipTM 8192/1024 |", "|---|---|---|---|---|---|"]
        for pid in b.index:
            lines.append(f"| {pid} | {a.loc[pid,'dockq']:.3f} | {b.loc[pid,'dockq']:.3f} | {dd[pid]:+.3f} | {a.loc[pid,'iptm']:.3f} | {b.loc[pid,'iptm']:.3f} |")
        lines.append(f"\nMean dDockQ (512/128 minus 8192/1024) = {dd.mean():+.3f}; max |dDockQ| = {dd.abs().max():.3f}.")
    # 6. timing from driver logs
    for log in args.driver_logs:
        p = Path(log)
        if not p.exists():
            continue
        txt = p.read_text(errors="ignore")
        runs = re.findall(r"^\[run \] (\S+)\s+(\S+)", txt, flags=re.M)
        dones = re.findall(r"^\[done\] (\S+)", txt, flags=re.M)
        times = {}
        stamps = [(pid, pd.Timestamp(ts)) for pid, ts in runs]
        for i, (pid, ts) in enumerate(stamps):
            nxt = stamps[i + 1][1] if i + 1 < len(stamps) else None
            if pid in dones and nxt is not None:
                times[pid] = (nxt - ts).total_seconds() / 60
        if times:
            s = pd.Series(times)
            lines += ["", f"## Wall time per target from `{p.name}` (RTX 4060 Laptop 8 GB; includes MSA server wait)", "",
                      f"- n = {len(s)} completed targets timed; median {s.median():.1f} min, min {s.min():.1f}, max {s.max():.1f} (last target of a batch is not timed)."]
    out = res_dir / "report.md"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")
    print("\n".join(lines[:14]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
