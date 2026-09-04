"""Download reference mmCIF files and write canonical FASTA files for all targets.

Reads ``data/targets/targets.csv`` (written by :mod:`cofoldbench.select_targets`),
fetches ``{id}.cif`` from files.rcsb.org (skipped when present) and writes
``data/targets/fasta/{id}.fasta`` with two records, chain ``A`` (entity 1) and
chain ``B`` (entity 2), using the canonical entity sequences from the Data API.

Usage
-----
    python -m cofoldbench.download [--config config.yaml] [--ids 9ASS 8RO8 ...]
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from . import rcsb
from .config import load_config, resolve

log = logging.getLogger("cofoldbench.download")


def load_targets(cfg: dict) -> pd.DataFrame:
    """Read ``targets.csv``."""
    return pd.read_csv(resolve(cfg, "targets_dir") / "targets.csv", dtype=str)


def write_fasta(row: pd.Series, fasta_dir: Path) -> Path:
    """Two-record FASTA: ``>{id}_A|entity_1|desc`` / ``>{id}_B|entity_2|desc``."""
    pid = row["pdb_id"]
    out = fasta_dir / f"{pid}.fasta"
    lines = []
    for tag, ent in (("A", "a"), ("B", "b")):
        desc = str(row[f"desc_{ent}"]).replace("\n", " ")
        lines.append(f">{pid}_{tag}|entity_{row[f'entity_{ent}']}|ref_chain_{row[f'label_asym_{ent}']}|{desc}")
        lines.append(str(row[f"seq_{ent}"]))
    out.write_text("\n".join(lines) + "\n")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--ids", nargs="*", help="subset of PDB ids (default: all targets)")
    ap.add_argument("--force", action="store_true", help="re-download existing mmCIF files")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    cfg = load_config(args.config)
    df = load_targets(cfg)
    if args.ids:
        df = df[df["pdb_id"].isin([x.upper() for x in args.ids])]
    cif_dir, fasta_dir = resolve(cfg, "cif_dir"), resolve(cfg, "fasta_dir")
    cif_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    n_dl = 0
    for _, row in df.iterrows():
        pid = row["pdb_id"]
        cif = cif_dir / f"{pid}.cif"
        if args.force or not cif.exists():
            rcsb.download_cif(pid, str(cif), base=cfg["selection"]["rcsb_files_url"])
            n_dl += 1
            log.info("downloaded %s", cif.name)
        write_fasta(row, fasta_dir)
    log.info("%d targets: %d mmCIF downloaded (%d already present), %d FASTA written",
             len(df), n_dl, len(df) - n_dl, len(df))
    return 0


if __name__ == "__main__":
    sys.exit(main())
