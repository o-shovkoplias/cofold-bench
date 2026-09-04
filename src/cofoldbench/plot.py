"""Figures from ``results/results.csv`` (top-ranked model per target and predictor).

1. ``figures/dockq_scatter.png``  DockQ Boltz-2 vs AF2-Multimer, CAPRI bands.
2. ``figures/dockq_bars.png``     per-target DockQ bars, both predictors.
3. ``figures/iptm_vs_dockq.png``  ipTM (and pDockQ) vs DockQ with Spearman rho.

Usage
-----
    python -m cofoldbench.plot [--config config.yaml]
"""

from __future__ import annotations

import argparse
import logging
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from .config import load_config, resolve  # noqa: E402

log = logging.getLogger("cofoldbench.plot")
LABEL = {"boltz2": "Boltz-2", "af2m": "AF2-Multimer v3 (ColabFold)"}
COLOR = {"boltz2": "#1b6ca8", "af2m": "#d1495b"}
CAPRI = [(0.23, "acceptable"), (0.49, "medium"), (0.80, "high")]


def top1(res: pd.DataFrame) -> pd.DataFrame:
    return res[res["rank"] == 1].copy()


def scatter(top: pd.DataFrame, out) -> None:
    w = top.pivot(index="pdb_id", columns="predictor", values="dockq").dropna()
    if w.empty or not {"boltz2", "af2m"} <= set(w.columns):
        log.warning("scatter skipped: need both predictors")
        return
    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    for y, name in CAPRI:
        ax.axhline(y, color="0.85", lw=0.8, zorder=0)
        ax.axvline(y, color="0.85", lw=0.8, zorder=0)
    ax.plot([0, 1], [0, 1], "--", color="0.6", lw=1)
    ax.scatter(w["af2m"], w["boltz2"], s=36, color="0.2", zorder=3)
    for pid, r in w.iterrows():
        ax.annotate(pid, (r["af2m"], r["boltz2"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"DockQ  {LABEL['af2m']}")
    ax.set_ylabel(f"DockQ  {LABEL['boltz2']}")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(f"n = {len(w)} heterodimers; Boltz-2 better in {(w['boltz2'] > w['af2m']).sum()}")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def bars(top: pd.DataFrame, out) -> None:
    w = top.pivot(index="pdb_id", columns="predictor", values="dockq")
    if w.empty:
        return
    w = w.sort_values(by=w.columns[0], ascending=False)
    x = np.arange(len(w))
    preds = [p for p in ("boltz2", "af2m") if p in w.columns]
    width = 0.8 / max(1, len(preds))
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(w)), 3.6))
    for i, p in enumerate(preds):
        ax.bar(x + (i - (len(preds) - 1) / 2) * width, w[p].fillna(0), width, label=LABEL[p], color=COLOR[p])
    for y, name in CAPRI:
        ax.axhline(y, color="0.7", lw=0.7, ls=":")
        ax.text(len(w) - 0.5, y, name, fontsize=6, va="bottom", ha="right", color="0.4")
    ax.set_xticks(x)
    ax.set_xticklabels(w.index, rotation=90, fontsize=7)
    ax.set_ylabel("DockQ (top-1 model)")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def iptm_vs_dockq(top: pd.DataFrame, out) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.9), sharey=True)
    for ax, xcol in zip(axes, ("iptm", "pdockq")):
        for p, g in top.groupby("predictor"):
            g = g.dropna(subset=[xcol, "dockq"])
            if len(g) < 3:
                continue
            rho, pval = spearmanr(g[xcol], g["dockq"])
            ax.scatter(g[xcol], g["dockq"], s=28, color=COLOR.get(p, "k"), label=f"{LABEL.get(p, p)}  rho={rho:.2f} (n={len(g)})")
        ax.set_xlabel(xcol if xcol != "iptm" else "ipTM")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=6.5, frameon=False, loc="upper left")
    axes[0].set_ylabel("DockQ (top-1 model)")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    cfg = load_config(args.config)
    res_csv = resolve(cfg, "results_dir") / "results.csv"
    fig_dir = resolve(cfg, "figures_dir")
    fig_dir.mkdir(parents=True, exist_ok=True)
    if not res_csv.exists():
        log.warning("%s not found - nothing to plot (predictions pending)", res_csv)
        return 0
    top = top1(pd.read_csv(res_csv))
    scatter(top, fig_dir / "dockq_scatter.png")
    bars(top, fig_dir / "dockq_bars.png")
    iptm_vs_dockq(top, fig_dir / "iptm_vs_dockq.png")
    log.info("figures written to %s", fig_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
