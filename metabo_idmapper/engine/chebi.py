"""ChEBI lookup via EBI OLS4 REST (no local ChEBI dump needed).

Verified approach from the skill: lowercase query + prefer a label-exact hit. Returns
candidate ChEBI ids with labels as evidence; the LLM decides which (if any) is right.
"""

from __future__ import annotations

import requests

_OLS = "https://www.ebi.ac.uk/ols4/api/search"
_HEADERS = {"User-Agent": "metabo-idmapper/0.1 (research)"}


def search(term: str, rows: int = 10) -> list[dict]:
    """Return [{chebi, label, exact}] for a name query against ChEBI."""
    try:
        r = requests.get(
            _OLS,
            params={"q": term.lower(), "ontology": "chebi", "rows": rows,
                    "fieldList": "obo_id,label"},
            headers=_HEADERS, timeout=30,
        )
        if r.status_code != 200:
            return []
        docs = r.json().get("response", {}).get("docs", [])
    except Exception:
        return []
    out = []
    tl = term.strip().lower()
    for d in docs:
        obo = d.get("obo_id")
        label = d.get("label")
        if not obo:
            continue
        out.append({"chebi": obo, "label": label,
                    "exact": bool(label and label.strip().lower() == tl)})
    # exact-label hits first
    out.sort(key=lambda x: (not x["exact"]))
    return out
