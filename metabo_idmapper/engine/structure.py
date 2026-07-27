"""PubChem structure lookup: name -> CID / InChIKey / formula / monoisotopic mass.

Uses the public PUG-REST endpoint via requests (pubchempy is a fallback). Returns
evidence only; the LLM decides whether the returned compound is the intended one.
"""

from __future__ import annotations

import time
import urllib.parse

import requests

from ..config import PUBCHEM_DELAY

_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_HEADERS = {"User-Agent": "metabo-idmapper/0.1 (research)"}


def _get(url: str, timeout: int = 30) -> dict | None:
    try:
        r = requests.get(url, headers=_HEADERS, timeout=timeout)
        time.sleep(PUBCHEM_DELAY)
        if r.status_code == 200:
            return r.json()
    except Exception:
        return None
    return None


def lookup(name: str) -> dict:
    """Return {found, cid, inchikey, formula, mono_mass, synonyms_sample, source}."""
    q = urllib.parse.quote(name)
    props = _get(
        f"{_BASE}/compound/name/{q}/property/"
        f"MolecularFormula,InChIKey,MonoisotopicMass,IUPACName/JSON"
    )
    out = {
        "found": False, "cid": None, "inchikey": None, "formula": None,
        "mono_mass": None, "iupac": None, "synonyms_sample": [], "source": "pubchem-rest",
    }
    if not props:
        return out
    try:
        p = props["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError):
        return out
    out.update(
        found=True,
        cid=p.get("CID"),
        inchikey=p.get("InChIKey"),
        formula=p.get("MolecularFormula"),
        mono_mass=float(p["MonoisotopicMass"]) if p.get("MonoisotopicMass") else None,
        iupac=p.get("IUPACName"),
    )
    # a few synonyms often carry KEGG/HMDB/ChEBI ids the LLM can hand to bridge_xref
    syn = _get(f"{_BASE}/compound/cid/{out['cid']}/synonyms/JSON")
    try:
        names = syn["InformationList"]["Information"][0]["Synonym"]
        out["synonyms_sample"] = names[:25]
    except (KeyError, IndexError, TypeError):
        pass
    return out
