"""Mouse-GEM (COBRA SBML) metabolite crosswalk.

Given accepted KEGG/HMDB/ChEBI ids, map to GEM MAM species (base id across
compartments) via the model's xref annotations. Reuses the verified A4 logic. The
model is large (~35 MB SBML); it is loaded once and the xref index is cached per path.
"""

from __future__ import annotations

import re
from collections import defaultdict
from functools import lru_cache

from ..config import GEM_MODEL, GEM_XREF_KEYS


def _anno_ids(m, key: str) -> list[str]:
    v = m.annotation.get(key)
    if v is None:
        return []
    if isinstance(v, str):
        v = [v]
    return [str(x) for x in v]


@lru_cache(maxsize=4)
def _index(model_path: str) -> dict:
    """Build {kegg,hmdb,chebi} -> set(base MAM id) plus base->name/compartments."""
    import cobra

    model = cobra.io.read_sbml_model(model_path)
    kegg: dict[str, set] = defaultdict(set)
    hmdb: dict[str, set] = defaultdict(set)
    chebi: dict[str, set] = defaultdict(set)
    meta: dict[str, dict] = defaultdict(lambda: {"name": "", "mams": set()})
    for m in model.metabolites:
        base = re.sub(r"[a-z]$", "", m.id)  # MAM00001c -> MAM00001
        meta[base]["name"] = m.name
        meta[base]["mams"].add(m.id)
        for k in _anno_ids(m, GEM_XREF_KEYS["kegg"]):
            kegg[k].add(base)
        for h in _anno_ids(m, GEM_XREF_KEYS["hmdb"]):
            hmdb[h].add(base)
        for c in _anno_ids(m, GEM_XREF_KEYS["chebi"]):
            c2 = c if c.upper().startswith("CHEBI:") else f"CHEBI:{c}"
            chebi[c2].add(base)
            chebi[c2.replace("CHEBI:", "")].add(base)
    return {
        "kegg": kegg, "hmdb": hmdb, "chebi": chebi, "meta": meta,
        "n_met": len(model.metabolites), "n_rxn": len(model.reactions),
        "n_base": len(meta),
    }


def crosswalk(kegg: str | None = None, hmdb: str | None = None,
              chebi: str | None = None, model_path: str | None = None) -> dict:
    """Map one metabolite (by any xref) to GEM MAM base ids.

    Priority KEGG > HMDB > ChEBI (KEGG is the most complete xref in Mouse-GEM).
    Returns {mapped, via, mam_bases[], mam_species[], gem_name}.
    """
    idx = _index(model_path or GEM_MODEL)
    for via, val, table in (("kegg", kegg, idx["kegg"]),
                            ("hmdb", hmdb, idx["hmdb"]),
                            ("chebi", chebi, idx["chebi"])):
        if not val:
            continue
        bases = table.get(str(val)) or table.get(str(val).replace("CHEBI:", ""))
        if bases:
            b = sorted(bases)
            species = sorted({s for base in b for s in idx["meta"][base]["mams"]})
            return {"mapped": True, "via": via, "mam_bases": b,
                    "mam_species": species, "gem_name": idx["meta"][b[0]]["name"]}
    return {"mapped": False, "via": None, "mam_bases": [], "mam_species": [],
            "gem_name": None}


def model_stats(model_path: str | None = None) -> dict:
    idx = _index(model_path or GEM_MODEL)
    return {k: idx[k] for k in ("n_met", "n_rxn", "n_base")}
