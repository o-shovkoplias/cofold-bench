"""Select post-cutoff protein-protein heterodimer targets from the RCSB PDB.

Pipeline
--------
1. RCSB Search API: entries with exactly 2 protein polymer entities, 1 copy each
   in the asymmetric unit, no nucleic acids, resolution <= 2.8 A, X-ray/cryo-EM,
   released after the configured date, deposited length in the 150-550 window.
2. RCSB Data API (``/core/entry``, ``/core/polymer_entity``): canonical sequences,
   lengths, descriptions, organisms, chain ids.  Filters: both chains >= 40 aa,
   canonical total length in window, antibody keywords -> labelled subset,
   one entry per (description pair).
3. Biotite verification on the downloaded mmCIF: biological assembly 1 has exactly
   two protein chains and they share >= ``min_interface_contacts`` residue pairs
   within 5 A (i.e. the pair really is a heterodimer with an interface).

Outputs (all under ``data/targets/``): ``candidates.csv`` (every entry examined,
with pass/fail reason), ``targets.csv`` (verified, ranked) and ``QUERY.md``
(the exact JSON query and the numbers at each stage).

Usage
-----
    python -m cofoldbench.select_targets [--config config.yaml] [--seed 2024] [--max-fetch 250]
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import re
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

from . import rcsb
from .config import load_config, resolve
from .structure import chains_summary, interface_residue_pairs, read_structure


def sequence_identity(seq_a: str, seq_b: str) -> float:
    """Global-alignment identity (BLOSUM62, Biotite) relative to the shorter sequence."""
    import biotite.sequence as bseq
    import biotite.sequence.align as balign

    a, b = bseq.ProteinSequence(seq_a), bseq.ProteinSequence(seq_b)
    matrix = balign.SubstitutionMatrix.std_protein_matrix()
    aln = balign.align_optimal(a, b, matrix, gap_penalty=(-10, -1), local=False, max_number=1)[0]
    ident = balign.get_sequence_identity(aln, mode="shortest")
    return float(ident)

log = logging.getLogger("cofoldbench.select")


def _first(x: Any, default: Any = None) -> Any:
    if isinstance(x, list):
        return x[0] if x else default
    return x if x is not None else default


def fetch_metadata(pdb_id: str, cache_dir: Path, cfg_sel: dict[str, Any]) -> dict[str, Any] | None:
    """Fetch entry + both polymer entities; cache the raw JSON under ``cache_dir``."""
    cache = cache_dir / f"{pdb_id}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    base = cfg_sel["rcsb_data_url"]
    e = rcsb.entry(pdb_id, base=base)
    if e is None:
        return None
    ent_ids = e.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    entities = []
    for eid in ent_ids:
        pe = rcsb.polymer_entity(pdb_id, eid, base=base)
        if pe is not None:
            entities.append(pe)
        time.sleep(0.05)
    meta = {"entry": e, "entities": entities}
    cache.write_text(json.dumps(meta))
    return meta


def summarise(meta: dict[str, Any]) -> dict[str, Any]:
    """Flatten the Data API JSON into the columns we keep."""
    e = meta["entry"]
    info = e.get("rcsb_entry_info", {})
    res = info.get("resolution_combined")
    row: dict[str, Any] = {
        "pdb_id": e["rcsb_id"],
        "title": e.get("struct", {}).get("title", ""),
        "method": _first([m.get("method") for m in e.get("exptl", [])], ""),
        "resolution": float(_first(res)) if res else None,
        "release_date": e.get("rcsb_accession_info", {}).get("initial_release_date", "")[:10],
        "n_entities": len(meta["entities"]),
    }
    for idx, pe in enumerate(meta["entities"]):
        tag = "ab"[idx] if idx < 2 else str(idx)
        poly = pe.get("entity_poly", {})
        ids = pe.get("rcsb_polymer_entity_container_identifiers", {})
        orgs = pe.get("rcsb_entity_source_organism", []) or []
        row[f"entity_{tag}"] = ids.get("entity_id")
        row[f"desc_{tag}"] = pe.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
        row[f"organism_{tag}"] = "; ".join(sorted({o.get("ncbi_scientific_name", "") for o in orgs if o.get("ncbi_scientific_name")}))
        row[f"type_{tag}"] = poly.get("rcsb_entity_polymer_type", "")
        row[f"seq_{tag}"] = (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
        row[f"len_{tag}"] = int(poly.get("rcsb_sample_sequence_length") or len(row[f"seq_{tag}"]))
        row[f"label_asym_{tag}"] = _first(ids.get("asym_ids"), "")
        row[f"auth_asym_{tag}"] = _first(ids.get("auth_asym_ids"), "")
    row["total_length"] = int(row.get("len_a", 0)) + int(row.get("len_b", 0))
    return row


def _has_keyword(text: str, keywords: list[str]) -> bool:
    text = text.lower()
    return any(re.search(rf"\b{re.escape(k)}\b", text) for k in keywords)


# Immunoglobulin variable-domain signatures (framework 2 / J segment).  Heavy (VH/VHH):
# "W[VIF]RQ...G" in FR2 and "WG.G...VTVSS" at the J end; light (VL): "WYQQ" FR2 and "FG.GT" J end.
_IG_HEAVY = re.compile(r"W[VIFL]R[QLKH].{2,6}G.{55,95}WG.G")
_IG_LIGHT = re.compile(r"W[YFL][QLH][QKH].{45,95}FG.G[TS]")


def is_antibody_chain(seq: str) -> bool:
    """Sequence heuristic: does this chain look like an immunoglobulin variable domain?"""
    return bool(_IG_HEAVY.search(seq) or _IG_LIGHT.search(seq))


def is_antibody_like(row: dict[str, Any], keywords: list[str]) -> bool:
    """True if the title/descriptions mention an antibody keyword or a chain looks like an Ig V-domain."""
    if _has_keyword(" ".join(str(row.get(k, "")) for k in ("title", "desc_a", "desc_b")), keywords):
        return True
    return is_antibody_chain(str(row.get("seq_a", ""))) or is_antibody_chain(str(row.get("seq_b", "")))


def is_antibody_only(row: dict[str, Any], keywords: list[str]) -> bool:
    """True if BOTH chains look like Ig variable domains (Fab/Fv heavy + light, no antigen)."""
    return is_antibody_chain(str(row.get("seq_a", ""))) and is_antibody_chain(str(row.get("seq_b", "")))


def metadata_filter(row: dict[str, Any], sel: dict[str, Any]) -> str:
    """Return '' if the entry passes the metadata filters, else the reason."""
    if row["n_entities"] != 2:
        return "not_two_entities"
    if row["type_a"] != "Protein" or row["type_b"] != "Protein":
        return "non_protein_entity"
    if min(row["len_a"], row["len_b"]) < sel["chain_length_min"]:
        return f"chain_shorter_than_{sel['chain_length_min']}"
    if not (sel["total_length_min"] <= row["total_length"] <= sel["total_length_max"]):
        return "canonical_total_length_out_of_window"
    if not row["seq_a"] or not row["seq_b"]:
        return "missing_canonical_sequence"
    if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", row["seq_a"] + row["seq_b"]):
        return "non_standard_residue_in_canonical_sequence"
    if row["desc_a"].strip().lower() == row["desc_b"].strip().lower():
        return "identical_entity_descriptions_homodimer_like"
    if sequence_identity(row["seq_a"], row["seq_b"]) >= sel.get("max_chain_identity", 0.9):
        return "chains_too_similar_homodimer_like"
    if is_antibody_only(row, sel["antibody_keywords"]):
        return "antibody_without_antigen"
    if row["release_date"] <= sel["min_release_date"]:
        return "release_date_not_after_cutoff"
    if row["resolution"] is None or row["resolution"] > sel["max_resolution"]:
        return "resolution"
    return ""


def biotite_verify(cif_path: Path, sel: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    """Check assembly 1 is a two-chain protein heterodimer with a real interface."""
    try:
        asm = read_structure(cif_path, assembly="1")
    except Exception as exc:  # noqa: BLE001 - report any parsing problem as a reason
        return f"biotite_assembly_error:{type(exc).__name__}", {}
    chains = chains_summary(asm)
    if len(chains) != 2:
        return f"assembly1_has_{len(chains)}_protein_chains", {"asm_chains": json.dumps(chains)}
    (ca, na), (cb, nb) = chains.items()
    contacts = interface_residue_pairs(asm, ca, cb, cutoff=5.0)
    extra = {
        "asm_chain_a": ca,
        "asm_chain_b": cb,
        "obs_len_a": na,
        "obs_len_b": nb,
        "interface_pairs_5A": contacts,
    }
    if contacts < sel["min_interface_contacts"]:
        return "too_few_interface_contacts", extra
    if min(na, nb) < sel["chain_length_min"]:
        return "observed_chain_too_short", extra
    return "", extra


def write_query_md(path: Path, query: dict[str, Any], stats: dict[str, Any], sel: dict[str, Any]) -> None:
    md = f"""# Target selection query (RCSB PDB)

Generated by `python -m cofoldbench.select_targets` on {stats['date']} (seed {stats['seed']}).

## Rationale

* **Heterodimer**: exactly 2 polymer entities, both protein, 2 deposited polymer
  instances (one copy of each entity in the asymmetric unit); no DNA/RNA/hybrid entities.
* **Quality**: X-ray or cryo-EM, resolution <= {sel['max_resolution']} A.
* **Size**: deposited polymer monomer count {sel['total_length_min']}-{sel['total_length_max']}
  (re-checked on canonical entity sequences) so both Boltz-2 and AF2-Multimer fit an 8 GB GPU;
  each chain >= {sel['chain_length_min']} aa (excludes peptide complexes).
* **Post-cutoff**: initial release date strictly after {sel['min_release_date']}.
  AlphaFold-Multimer v3 (AlphaFold v2.3.0) was trained on PDB entries released before
  2021-09-30 (DeepMind technical note v2.3.0); Boltz-2 on PDB entries released before
  2023-06-01 (Passaro et al. 2025, bioRxiv 10.1101/2025.06.14.659707, Data section).
* **Antibodies**: entries whose title/entity descriptions contain
  {sel['antibody_keywords']} are moved to a labelled `antibody` subset (max {sel['antibody_subset_max']}).
* **Diversity**: one entry per unordered pair of entity descriptions.
* **Biotite verification**: biological assembly 1 (`pdbx.get_assembly`) has exactly two
  protein chains sharing >= {sel['min_interface_contacts']} residue pairs within 5 A.

## Search API request (`POST {sel['rcsb_search_url']}`)

```json
{json.dumps(query, indent=2)}
```

## Data API

For every candidate: `GET {sel['rcsb_data_url']}/entry/{{id}}` and
`GET {sel['rcsb_data_url']}/polymer_entity/{{id}}/{{entity_id}}` for both entities
(canonical sequence `entity_poly.pdbx_seq_one_letter_code_can`, length
`entity_poly.rcsb_sample_sequence_length`, description, source organism,
`asym_ids`/`auth_asym_ids`).  mmCIF references from `{sel['rcsb_files_url']}/{{id}}.cif`.

## Numbers

| stage | count |
|---|---|
| Search API hits | {stats['n_hits']} |
| entries fetched from Data API (deterministic shuffle, seed {stats['seed']}) | {stats['n_fetched']} |
| pass metadata filters (non-antibody) | {stats['n_meta_pass']} |
| pass metadata filters (antibody-like) | {stats['n_meta_antibody']} |
| verified with Biotite (non-antibody) | {stats['n_verified']} |
| verified with Biotite (antibody subset) | {stats['n_verified_antibody']} |
| written to `targets.csv` | {stats['n_targets']} |

Per-entry reasons are in `candidates.csv` (column `status`).
"""
    path.write_text(md, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--seed", type=int, default=2024, help="shuffle seed for the candidate order")
    ap.add_argument("--max-fetch", type=int, default=250, help="max entries to pull from the Data API")
    ap.add_argument("--n-verify", type=int, default=None, help="stop after this many verified non-antibody targets (default: selection.n_candidates)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(levelname)s %(message)s")

    cfg = load_config(args.config)
    sel = cfg["selection"]
    targets_dir = resolve(cfg, "targets_dir")
    cif_dir = resolve(cfg, "cif_dir")
    raw_dir = Path(cfg["_root"]) / "results" / "raw" / "rcsb_meta"
    for d in (targets_dir, cif_dir, raw_dir):
        d.mkdir(parents=True, exist_ok=True)
    n_want = args.n_verify or sel["n_candidates"]

    query = rcsb.heterodimer_query(
        sel["min_release_date"], sel["max_resolution"], sel["methods"],
        sel["total_length_min"], sel["total_length_max"],
    )
    hits = rcsb.search_all(query, url=sel["rcsb_search_url"])
    log.info("Search API: %d hits", len(hits))
    (raw_dir / "search_hits.json").write_text(json.dumps(hits))

    order = list(hits)
    random.Random(args.seed).shuffle(order)

    rows: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str]] = set()
    n_ok = n_ab = 0
    for i, pid in enumerate(order[: args.max_fetch]):
        if n_ok >= n_want and n_ab >= sel["antibody_subset_max"]:
            break
        meta = fetch_metadata(pid, raw_dir, sel)
        if meta is None:
            rows.append({"pdb_id": pid, "status": "data_api_404"})
            continue
        row = summarise(meta)
        row["antibody_like"] = is_antibody_like(row, sel["antibody_keywords"])
        reason = metadata_filter(row, sel)
        if not reason:
            pair = tuple(sorted((row["desc_a"].lower(), row["desc_b"].lower())))
            if pair in seen_pairs:
                reason = "duplicate_description_pair"
            else:
                seen_pairs.add(pair)
        if not reason and row["antibody_like"] and n_ab >= sel["antibody_subset_max"]:
            reason = "antibody_subset_full"
        if not reason and not row["antibody_like"] and n_ok >= n_want:
            reason = "enough_targets"
        if reason:
            row["status"] = reason
            rows.append(row)
            log.debug("%s skip: %s", pid, reason)
            continue
        # Biotite verification (download the reference)
        cif_path = cif_dir / f"{pid}.cif"
        if not cif_path.exists():
            try:
                rcsb.download_cif(pid, str(cif_path), base=sel["rcsb_files_url"])
            except Exception as exc:  # noqa: BLE001
                row["status"] = f"cif_download_error:{type(exc).__name__}"
                rows.append(row)
                continue
        reason, extra = biotite_verify(cif_path, sel)
        row.update(extra)
        if reason:
            row["status"] = reason
            cif_path.unlink(missing_ok=True)
            rows.append(row)
            log.info("%s rejected by Biotite: %s", pid, reason)
            continue
        row["status"] = "ok"
        if row["antibody_like"]:
            n_ab += 1
        else:
            n_ok += 1
        rows.append(row)
        log.info("%s OK  %s | %s  (%d+%d aa, %.2f A, %s, %d iface pairs)", pid, row["desc_a"][:40],
                 row["desc_b"][:40], row["len_a"], row["len_b"], row["resolution"], row["release_date"],
                 row["interface_pairs_5A"])
        time.sleep(0.05)

    cand = pd.DataFrame(rows)
    cand.to_csv(targets_dir / "candidates.csv", index=False)

    ok = cand[cand["status"] == "ok"].copy()
    ok["subset"] = ok["antibody_like"].map({True: "antibody", False: "main"})
    main_df = ok[ok["subset"] == "main"].head(sel["n_candidates"]).copy()
    ab_df = ok[ok["subset"] == "antibody"].head(sel["antibody_subset_max"]).copy()
    main_df["role"] = ["primary" if i < sel["n_targets"] else "backup" for i in range(len(main_df))]
    ab_df["role"] = "antibody_subset"
    targets = pd.concat([main_df, ab_df], ignore_index=True)
    targets.insert(0, "rank", range(1, len(targets) + 1))
    cols = [
        "rank", "pdb_id", "role", "subset", "title", "method", "resolution", "release_date",
        "entity_a", "desc_a", "organism_a", "len_a", "label_asym_a", "auth_asym_a",
        "entity_b", "desc_b", "organism_b", "len_b", "label_asym_b", "auth_asym_b",
        "total_length", "asm_chain_a", "asm_chain_b", "obs_len_a", "obs_len_b", "interface_pairs_5A",
        "seq_a", "seq_b",
    ]
    targets[cols].to_csv(targets_dir / "targets.csv", index=False)

    # keep only per-target files that belong to selected targets
    keep = set(targets["pdb_id"])
    for d, pat in ((cif_dir, "*.cif"), (resolve(cfg, "fasta_dir"), "*.fasta"),
                   (resolve(cfg, "boltz_inputs"), "*.yaml"), (resolve(cfg, "colabfold_inputs"), "*.fasta")):
        if d.exists():
            for f in d.glob(pat):
                if f.stem not in keep:
                    f.unlink()

    stats = {
        "date": time.strftime("%Y-%m-%d"), "seed": args.seed, "n_hits": len(hits),
        "n_fetched": int(cand["pdb_id"].nunique()),
        "n_meta_pass": int(((cand["status"] == "ok") | cand["status"].str.contains("interface|assembly|observed|biotite|cif_", regex=True)).sum() - int((cand.get("antibody_like", pd.Series(dtype=bool)) == True).sum())) if len(cand) else 0,
        "n_meta_antibody": int((cand.get("antibody_like", pd.Series(dtype=bool)) == True).sum()) if len(cand) else 0,
        "n_verified": len(main_df), "n_verified_antibody": len(ab_df), "n_targets": len(targets),
    }
    write_query_md(targets_dir / "QUERY.md", query, stats, sel)
    log.info("targets.csv: %d rows (%d primary, %d backup, %d antibody)", len(targets),
             int((targets["role"] == "primary").sum()), int((targets["role"] == "backup").sum()), len(ab_df))
    return 0


if __name__ == "__main__":
    sys.exit(main())
