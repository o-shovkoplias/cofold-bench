"""Minimal clients for the RCSB Search, Data and Files APIs (no third-party SDK)."""

from __future__ import annotations

import json
import time
from typing import Any, Iterable

import requests

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/rest/v1/core"
FILES_URL = "https://files.rcsb.org/download"

_HEADERS = {"User-Agent": "cofold-bench/0.1 (https://github.com/oshovkoplias/cofold-bench)"}


def _terminal(attribute: str, operator: str, value: Any) -> dict[str, Any]:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def heterodimer_query(
    min_release_date: str,
    max_resolution: float,
    methods: Iterable[str],
    total_length_min: int,
    total_length_max: int,
    rows: int = 1000,
    start: int = 0,
) -> dict[str, Any]:
    """Build the RCSB Search API JSON for post-cutoff protein-protein heterodimers.

    The query asks for entries with exactly two polymer entities, both protein,
    two deposited polymer instances in total (i.e. one copy of each entity in
    the asymmetric unit), no DNA/RNA/hybrid entities, resolution <= max_resolution,
    X-ray or cryo-EM, released strictly after ``min_release_date`` and with a
    deposited polymer monomer count inside the length window.  The length window
    is re-checked later from canonical entity sequences.
    """
    nodes = [
        _terminal("rcsb_entry_info.polymer_entity_count", "equals", 2),
        _terminal("rcsb_entry_info.polymer_entity_count_protein", "equals", 2),
        _terminal("rcsb_entry_info.deposited_polymer_entity_instance_count", "equals", 2),
        _terminal("rcsb_entry_info.polymer_entity_count_DNA", "equals", 0),
        _terminal("rcsb_entry_info.polymer_entity_count_RNA", "equals", 0),
        _terminal("rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid", "equals", 0),
        _terminal("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution),
        _terminal("exptl.method", "in", list(methods)),
        _terminal(
            "rcsb_accession_info.initial_release_date",
            "greater",
            f"{min_release_date}T00:00:00Z",
        ),
        _terminal(
            "rcsb_entry_info.deposited_polymer_monomer_count",
            "range",
            {
                "from": total_length_min,
                "to": total_length_max,
                "include_lower": True,
                "include_upper": True,
            },
        ),
    ]
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "sort": [{"sort_by": "rcsb_accession_info.initial_release_date", "direction": "asc"}],
            "results_content_type": ["experimental"],
        },
    }


def search(query: dict[str, Any], url: str = SEARCH_URL, timeout: int = 120) -> tuple[int, list[str]]:
    """POST a Search API query; return ``(total_count, [pdb_id, ...])``."""
    r = requests.post(url, json=query, headers=_HEADERS, timeout=timeout)
    if r.status_code == 204:  # no results
        return 0, []
    r.raise_for_status()
    d = r.json()
    return int(d.get("total_count", 0)), [x["identifier"] for x in d.get("result_set", [])]


def search_all(query: dict[str, Any], url: str = SEARCH_URL, page: int = 1000) -> list[str]:
    """Paginate through all hits of a Search API query."""
    ids: list[str] = []
    start = 0
    while True:
        q = json.loads(json.dumps(query))
        q["request_options"]["paginate"] = {"start": start, "rows": page}
        total, chunk = search(q, url=url)
        ids.extend(chunk)
        start += page
        if not chunk or start >= total:
            break
        time.sleep(0.2)
    return ids


def _get(url: str, timeout: int = 60, retries: int = 3) -> dict[str, Any] | None:
    for attempt in range(retries):
        r = requests.get(url, headers=_HEADERS, timeout=timeout)
        if r.status_code == 404:
            return None
        if r.status_code == 429 or r.status_code >= 500:
            time.sleep(1.5 * (attempt + 1))
            continue
        r.raise_for_status()
        return r.json()
    r.raise_for_status()
    return None


def entry(pdb_id: str, base: str = DATA_URL) -> dict[str, Any] | None:
    """``/core/entry/{id}``."""
    return _get(f"{base}/entry/{pdb_id.upper()}")


def polymer_entity(pdb_id: str, entity_id: str | int, base: str = DATA_URL) -> dict[str, Any] | None:
    """``/core/polymer_entity/{id}/{entity_id}``."""
    return _get(f"{base}/polymer_entity/{pdb_id.upper()}/{entity_id}")


def assembly(pdb_id: str, assembly_id: str | int = 1, base: str = DATA_URL) -> dict[str, Any] | None:
    """``/core/assembly/{id}/{assembly_id}``."""
    return _get(f"{base}/assembly/{pdb_id.upper()}/{assembly_id}")


def download_cif(pdb_id: str, dest: str, base: str = FILES_URL, timeout: int = 120) -> str:
    """Download ``{id}.cif`` from files.rcsb.org to ``dest``; return the path."""
    r = requests.get(f"{base}/{pdb_id.upper()}.cif", headers=_HEADERS, timeout=timeout)
    r.raise_for_status()
    with open(dest, "wb") as fh:
        fh.write(r.content)
    return dest
