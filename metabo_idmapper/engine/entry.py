"""ID -> its OWN database record (the back-check).

Every other evidence route in this package goes NAME -> id. That direction cannot detect a
plausible-but-wrong id, because the wrong id arrives attached to the right name (a bad xref
link, a PubChem entry that merges two compounds as synonyms, a class id handed to a species
query). The only evidence that settles it is the record the id itself resolves to: KEGG
C05472's own NAME field says "Urocortisol", which no amount of name->id searching reveals.

So this module fetches, for a given id, what that id IS: the DB's own name list, formula and
mass. Callers hand the returned names to `isomer.compare` to get a verdict. Evidence only.

KEGG goes through KEGGREST (batched, 10 ids per request — the KEGG API limit); ChEBI through
OLS4; PubChem through PUG-REST. HMDB has no open record API: back-check an HMDB accession by
bridging it to PubChem/InChIKey and checking THAT record (see tools.id_name_check).
"""

from __future__ import annotations

import urllib.parse

import requests

from . import chebi as _chebi
from . import structure as _struct
from .rcall import run_r

_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_HEADERS = {"User-Agent": "metabo-idmapper/0.2 (research)"}


def kegg_entries(ids: list[str]) -> dict[str, dict]:
    """{kegg_id: {found, names[], formula, exact_mass, chebi, pubchem, lipidmaps}} via keggGet."""
    ids = [str(i).strip() for i in ids if i]
    if not ids:
        return {}
    r = run_r("kegg_entry.R", {"ids": ids})
    if not r["ok"]:
        return {i: {"found": False, "error": r.get("error"), "stderr": r.get("stderr"),
                    "names": []} for i in ids}
    out: dict[str, dict] = {}
    for row in r["result"]:
        names = row.get("names") or []
        if isinstance(names, str):
            names = [names]
        out[str(row.get("kegg"))] = {
            "found": bool(row.get("found")), "names": [str(n) for n in names],
            "formula": row.get("formula"), "exact_mass": row.get("exact_mass"),
            "mol_weight": row.get("mol_weight"), "chebi": row.get("chebi"),
            "pubchem": row.get("pubchem"), "lipidmaps": row.get("lipidmaps"),
            "source": "keggrest-keggGet"}
    for i in ids:
        out.setdefault(i, {"found": False, "names": [], "source": "keggrest-keggGet"})
    return out


# OLS4 keeps the structural facts under `annotation`, and the formula key is NOT "formula" —
# reading the obvious name returned None for every ChEBI id, so the formula cross-check on that
# route quietly compared nothing at all.
_CHEBI_FORMULA_KEYS = ("generalized_empirical_formula", "formula", "Formula")


def _anno(term: dict, *keys):
    anno = term.get("annotation") or {}
    for k in keys:
        v = anno.get(k)
        if isinstance(v, list):
            v = v[0] if v else None
        if v:
            return str(v)
    return None


def chebi_entry(chebi_id: str) -> dict:
    """ChEBI id -> {found, names[], formula, mono_mass, inchikey} via OLS4 term lookup."""
    cid = str(chebi_id).strip()
    if not cid.upper().startswith("CHEBI:"):
        cid = f"CHEBI:{cid}"
    out = {"found": False, "names": [], "formula": None, "mono_mass": None,
           "inchikey": None, "source": "chebi-ols4-term"}
    try:
        r = requests.get("https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms",
                         params={"obo_id": cid}, headers=_HEADERS, timeout=30)
        if r.status_code != 200:
            return out
        terms = r.json().get("_embedded", {}).get("terms", [])
    except Exception:
        return out
    if not terms:
        return out
    t = terms[0]
    names = [t.get("label")] + list(t.get("synonyms") or [])
    mass = _anno(t, "monoisotopic_mass", "mass")
    out.update(found=True, names=[n for n in names if n],
               formula=_anno(t, *_CHEBI_FORMULA_KEYS),
               mono_mass=float(mass) if mass else None,
               inchikey=_anno(t, "inchi_key_string", "inchikey"))
    return out


def pubchem_entry(cid: str) -> dict:
    """PubChem CID -> {found, names[], formula, mono_mass, inchikey}.

    The synonym list is the important part: when one PubChem entry carries the names of TWO
    different compounds (the verified Lc3Cer / nLc4Cer case), those names sit side by side in
    this list — which is how a shared-id collision becomes visible.
    """
    cid = str(cid).strip()
    out = {"found": False, "names": [], "formula": None, "mono_mass": None,
           "inchikey": None, "source": "pubchem-rest-cid"}
    props = _struct.get_json(f"{_BASE}/compound/cid/{urllib.parse.quote(cid)}/property/"
                         f"MolecularFormula,InChIKey,MonoisotopicMass/JSON")
    if not props:
        return out
    try:
        p = props["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError):
        return out
    out.update(found=True, formula=p.get("MolecularFormula"), inchikey=p.get("InChIKey"),
               mono_mass=float(p["MonoisotopicMass"]) if p.get("MonoisotopicMass") else None)
    syn = _struct.get_json(f"{_BASE}/compound/cid/{cid}/synonyms/JSON")
    try:
        out["names"] = syn["InformationList"]["Information"][0]["Synonym"][:25]
    except (KeyError, IndexError, TypeError):
        pass
    return out


def fetch(db: str, ids: list[str]) -> dict[str, dict]:
    """Uniform entry fetch for one db over many ids: {id: record}."""
    db = db.lower()
    if db == "kegg":
        return kegg_entries(ids)
    if db == "chebi":
        return {str(i): chebi_entry(i) for i in ids if i}
    if db == "pubchem":
        return {str(i): pubchem_entry(i) for i in ids if i}
    if db == "hmdb":
        return {str(i): {"found": False, "names": [], "source": "none",
                         "note": "HMDB has no open record API — bridge the accession to "
                                 "PubChem/InChIKey (bridge_xref) and back-check that record"}
                for i in ids if i}
    return {str(i): {"found": False, "names": [], "error": f"unsupported db {db}"} for i in ids}


def search_chebi(term: str, rows: int = 6) -> list[dict]:
    """Thin re-export so callers need only this module for entry-level evidence."""
    return _chebi.search(term, rows=rows)
