"""Parse predictor confidence files.

Boltz-2 (``predictions/<name>/confidence_<name>_model_<k>.json``, documented in
``docs/prediction.md`` of https://github.com/jwohlwend/boltz): keys
``confidence_score, ptm, iptm, ligand_iptm, protein_iptm, complex_plddt,
complex_iplddt, complex_pde, complex_ipde, chains_ptm, pair_chains_iptm``.
Per-token pLDDT is also written into the B-factor column of the mmCIF (0-1 scale).

ColabFold (``<jobname>_scores_rank_00k_alphafold2_multimer_v3_model_m_seed_00s.json``,
written in ``colabfold/batch.py``): keys ``plddt`` (per-residue list, 0-100),
``max_pae``, ``pae`` (matrix), ``ptm``, ``iptm`` and, for complexes in recent
versions, ipSAE-family scores (``ipsae``, ``pdockq2`` ...).  Ranking is by
``0.8*iptm + 0.2*ptm`` ("multimer" metric) for multimer models; the model file is
``<jobname>_unrelaxed_rank_00k_..._model_m_seed_00s.pdb`` with pLDDT in B-factors.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, asdict
from pathlib import Path


@dataclass
class Confidence:
    """Predictor-agnostic confidence summary of one model."""

    predictor: str
    model_file: str
    rank: int
    iptm: float
    ptm: float
    plddt: float
    ranking_score: float
    extra: dict

    def as_row(self) -> dict:
        d = asdict(self)
        d.pop("extra")
        return d


def _f(x, default=float("nan")) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------- Boltz-2
_BOLTZ_CONF = re.compile(r"^confidence_(?P<name>.+)_model_(?P<k>\d+)\.json$")


def boltz_models(pred_dir: str | Path) -> list[Confidence]:
    """All Boltz-2 models in ``predictions/<name>/``, sorted by confidence_score (desc)."""
    pred_dir = Path(pred_dir)
    out: list[Confidence] = []
    for js in sorted(pred_dir.glob("confidence_*_model_*.json")):
        m = _BOLTZ_CONF.match(js.name)
        if not m:
            continue
        k = int(m.group("k"))
        d = json.loads(js.read_text())
        model = pred_dir / f"{m.group('name')}_model_{k}.cif"
        if not model.exists():
            model = pred_dir / f"{m.group('name')}_model_{k}.pdb"
        out.append(
            Confidence(
                predictor="boltz2",
                model_file=str(model),
                rank=-1,
                iptm=_f(d.get("protein_iptm", d.get("iptm"))),
                ptm=_f(d.get("ptm")),
                plddt=_f(d.get("complex_plddt")) * 100.0,
                ranking_score=_f(d.get("confidence_score")),
                extra={k2: d.get(k2) for k2 in ("iptm", "protein_iptm", "complex_iplddt", "complex_pde", "complex_ipde", "pair_chains_iptm", "chains_ptm")},
            )
        )
    out.sort(key=lambda c: -c.ranking_score)
    for i, c in enumerate(out):
        c.rank = i + 1
    return out


# ---------------------------------------------------------------- ColabFold
_CF_SCORES = re.compile(
    r"^(?P<job>.+)_scores_rank_(?P<rank>\d{3})_(?P<model>alphafold2_[a-z0-9_]+_model_\d+_seed_\d+)\.json$"
)


def colabfold_models(result_dir: str | Path, jobname: str | None = None) -> list[Confidence]:
    """All ColabFold models in ``result_dir`` (one job), sorted by rank."""
    result_dir = Path(result_dir)
    out: list[Confidence] = []
    for js in sorted(result_dir.glob("*_scores_rank_*.json")):
        m = _CF_SCORES.match(js.name)
        if not m or (jobname and m.group("job") != jobname):
            continue
        d = json.loads(js.read_text())
        plddt = d.get("plddt", [])
        mean_plddt = sum(plddt) / len(plddt) if plddt else float("nan")
        iptm, ptm = _f(d.get("iptm")), _f(d.get("ptm"))
        rank = int(m.group("rank"))
        stem = f"{m.group('job')}_unrelaxed_rank_{rank:03d}_{m.group('model')}"
        model = result_dir / f"{stem}.pdb"
        relaxed = result_dir / f"{m.group('job')}_relaxed_rank_{rank:03d}_{m.group('model')}.pdb"
        out.append(
            Confidence(
                predictor="af2m",
                model_file=str(relaxed if relaxed.exists() else model),
                rank=rank,
                iptm=iptm,
                ptm=ptm,
                plddt=mean_plddt,
                ranking_score=0.8 * iptm + 0.2 * ptm,
                extra={k: d.get(k) for k in ("max_pae", "ipsae", "pdockq2", "actifptm") if k in d},
            )
        )
    out.sort(key=lambda c: c.rank)
    return out
