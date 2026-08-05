"""Genome-scale model (COBRA SBML/JSON) metabolite crosswalk + direct model search.

Two routes, because the xref route alone is not enough. In a verified Recon3D run only 8 of
27 compounds crosswalked by KEGG/HMDB/ChEBI annotation — GEMs annotate lipid species very
sparsely — while the model actually CONTAINED 21 of them under names no xref would reach
("Chenodeoxyglycocholate" for GCDCA, "Sialyl-3-paragloboside" for the a2,3-sialyl GSL,
"Lactoneotetraosylceramide" for nLc4Cer). So:

  crosswalk()  — annotation xrefs (KEGG > HMDB > ChEBI), exact and cheap.
  search()     — the model's own NAME / FORMULA / MASS / id text, which is how the remaining
                 species are found. Returns same-formula groups so a formula-identical
                 isomer pair (pcholar_hs vs pcholn204_hs for LPC 20:4) is surfaced as an
                 ambiguity to resolve rather than silently collapsed to the first hit.

Both emit evidence only. `pool_species` marks an R-group placeholder (crm_hs, sphmyln_hs):
present in the model but unable to represent chain length, so it is a CLASS-level input.

The model is large (35 MB SBML); the index is built once and cached per path.
"""

from __future__ import annotations

import re
from collections import defaultdict
from functools import lru_cache

from ..config import GEM_MODEL, GEM_XREF_KEYS

_WORD_RE = re.compile(r"[a-z0-9]+")


def _anno_ids(m, key: str) -> list[str]:
    v = m.annotation.get(key)
    if v is None:
        return []
    if isinstance(v, str):
        v = [v]
    return [str(x) for x in v]


def base_id(met) -> str:
    """Strip the compartment suffix from a metabolite id.

    Handles both id conventions: Mouse-GEM 'MAM01234c' -> 'MAM01234' and Recon3D
    'tdchola_c' / 'met__L_c' -> 'tdchola' / 'met__L'. The naive "drop one trailing lowercase
    letter" rule left a dangling underscore on Recon3D ids ('tdchola_'), which is not a
    usable accession — a model input built from it silently fails to resolve.
    """
    mid = str(met.id)
    comp = str(getattr(met, "compartment", "") or "")
    if comp and mid.endswith(comp) and len(mid) > len(comp):
        return mid[: -len(comp)].rstrip("_") or mid
    m = re.match(r"^(.+?)_([a-z]{1,2})$", mid)
    if m:
        return m.group(1)
    return re.sub(r"[a-z]$", "", mid)


def _load(model_path: str):
    import cobra

    p = str(model_path)
    if p.endswith(".json"):
        return cobra.io.load_json_model(p)
    if p.endswith((".mat",)):
        return cobra.io.load_matlab_model(p)
    if p.endswith((".yml", ".yaml")):
        return cobra.io.load_yaml_model(p)
    return cobra.io.read_sbml_model(p)


def _tokens(s: str) -> set[str]:
    return set(_WORD_RE.findall((s or "").lower()))


def _grams(s: str, n: int = 3) -> set[str]:
    """Character n-grams of the alphanumeric-only lowercase string.

    Chemical names are agglutinative — "Glycochenodeoxycholic acid" and
    "Chenodeoxyglycocholate" are the same compound but share NO whole token. Character
    n-grams do overlap, which is what makes the model's own trivial names reachable.
    """
    flat = re.sub(r"[^a-z0-9]", "", (s or "").lower())
    return {flat[i:i + n] for i in range(max(len(flat) - n + 1, 0))}


def _dice(a: set[str], b: set[str]) -> float:
    return 0.0 if not a or not b else 2 * len(a & b) / (len(a) + len(b))


def _canon_formula(f: str | None) -> str | None:
    """Canonical element-count string, so C28H51NO7P and c28h51no7p compare equal.
    R-group / placeholder formulas are returned uppercased as-is (they cannot be parsed)."""
    if not f:
        return None
    raw = str(f).strip()
    try:
        from molmass import Formula
        return str(Formula(raw).formula)
    except Exception:
        return raw.upper()


@lru_cache(maxsize=4)
def _index(model_path: str) -> dict:
    """Build the xref, name, formula and mass indexes for one model (cached per path)."""
    model = _load(model_path)
    kegg: dict[str, set] = defaultdict(set)
    hmdb: dict[str, set] = defaultdict(set)
    chebi: dict[str, set] = defaultdict(set)
    by_formula: dict[str, set] = defaultdict(set)
    meta: dict[str, dict] = {}
    for m in model.metabolites:
        base = base_id(m)
        rec = meta.setdefault(base, {"name": "", "names": set(), "mams": set(),
                                     "compartments": set(), "formula": None,
                                     "charge": None, "kegg": set(), "hmdb": set(),
                                     "chebi": set()})
        rec["name"] = rec["name"] or (m.name or "")
        if m.name:
            rec["names"].add(m.name)
        rec["mams"].add(m.id)
        if getattr(m, "compartment", None):
            rec["compartments"].add(str(m.compartment))
        rec["formula"] = rec["formula"] or getattr(m, "formula", None)
        if rec["charge"] is None:
            rec["charge"] = getattr(m, "charge", None)
        for k in _anno_ids(m, GEM_XREF_KEYS["kegg"]):
            kegg[k].add(base)
            rec["kegg"].add(k)
        for h in _anno_ids(m, GEM_XREF_KEYS["hmdb"]):
            hmdb[h].add(base)
            rec["hmdb"].add(h)
        for c in _anno_ids(m, GEM_XREF_KEYS["chebi"]):
            c2 = c if c.upper().startswith("CHEBI:") else f"CHEBI:{c}"
            chebi[c2].add(base)
            chebi[c2.replace("CHEBI:", "")].add(base)
            rec["chebi"].add(c2)

    for base, rec in meta.items():
        rec["mams"] = sorted(rec["mams"])
        rec["compartments"] = sorted(rec["compartments"])
        rec["names"] = sorted(rec["names"])
        for k in ("kegg", "hmdb", "chebi"):
            rec[k] = sorted(rec[k])
        rec["canon_formula"] = _canon_formula(rec["formula"])
        rec["pool_species"] = bool(rec["formula"] and re.search(r"[RX]", str(rec["formula"])))
        rec["mono_mass"] = None
        if rec["canon_formula"] and not rec["pool_species"]:
            try:
                from molmass import Formula
                rec["mono_mass"] = float(Formula(rec["canon_formula"]).monoisotopic_mass)
            except Exception:
                rec["mono_mass"] = None
        rec["tokens"] = _tokens(" ".join(rec["names"]))
        # grams of the name AND of the name with trailing qualifiers stripped: a model name
        # like "1-Linoleoylglycerophosphocholine (Delta 9,12)" is the same species name plus a
        # parenthetical, and the qualifier alone must not push it below a generic class hit.
        rec["grams"] = {nm: _grams(nm) for nm in rec["names"]}
        rec["grams_bare"] = {nm: _grams(re.sub(r"\([^)]*\)", " ", nm)) for nm in rec["names"]}
        if rec["canon_formula"]:
            by_formula[rec["canon_formula"]].add(base)

    return {"kegg": kegg, "hmdb": hmdb, "chebi": chebi, "meta": meta,
            "by_formula": by_formula, "path": model_path,
            "n_met": len(model.metabolites), "n_rxn": len(model.reactions),
            "n_base": len(meta)}


def crosswalk(kegg: str | None = None, hmdb: str | None = None,
              chebi: str | None = None, model_path: str | None = None) -> dict:
    """Map one metabolite (by any xref) to model base ids.

    Priority KEGG > HMDB > ChEBI. Returns {mapped, via, mam_bases[], mam_species[], gem_name}.
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


def _hit(base: str, rec: dict, score: float, kind: str) -> dict:
    return {"mam_base": base, "gem_name": rec["name"], "all_names": rec["names"],
            "formula": rec["formula"], "canon_formula": rec["canon_formula"],
            "mono_mass": rec["mono_mass"], "charge": rec["charge"],
            "compartments": rec["compartments"], "mam_species": rec["mams"],
            "pool_species": rec["pool_species"], "kegg": rec["kegg"], "hmdb": rec["hmdb"],
            "chebi": rec["chebi"], "score": round(score, 1), "match": kind}


def search(query: str | None = None, queries: list[str] | None = None,
           formula: str | None = None,
           mass: float | None = None, mass_tol_ppm: float = 25.0,
           kegg: str | None = None, hmdb: str | None = None, chebi: str | None = None,
           model_path: str | None = None, max_hits: int = 25,
           min_score: float = 30.0) -> dict:
    """Search a model's own metabolite records by name / id text / formula / mass / xref.

    Pass several name spellings at once with `queries=` (e.g. the reported name plus
    `shorthand.variants(name)`); each hit records which query matched it.

    Name scoring (deterministic): exact name 100, exact id 96, whole-query substring 70-85 by
    length ratio, token-overlap Jaccard x 60, character-trigram Dice x 58 (the route that
    reaches agglutinative trivial names). Returns hits sorted by score plus
    `same_formula_groups` — formula-identical hit sets, i.e. positional/stereo isomer
    candidates that no formula or mass check can separate.
    """
    idx = _index(model_path or GEM_MODEL)
    meta = idx["meta"]
    hits: dict[str, dict] = {}

    def put(base: str, score: float, kind: str, matched: str | None = None) -> None:
        cur = hits.get(base)
        if cur is None or score > cur["score"]:
            h = _hit(base, meta[base], score, kind)
            if matched:
                h["matched_query"] = matched
            hits[base] = h

    for via, val in (("kegg", kegg), ("hmdb", hmdb), ("chebi", chebi)):
        if not val:
            continue
        for base in sorted(idx[via].get(str(val), set())
                           | idx[via].get(str(val).replace("CHEBI:", ""), set())):
            put(base, 100.0, f"xref:{via}")

    if formula:
        canon = _canon_formula(formula)
        for base in sorted(idx["by_formula"].get(canon, set())):
            put(base, 95.0, "formula")

    if mass is not None:
        win = float(mass) * mass_tol_ppm / 1e6
        for base, rec in meta.items():
            mm = rec["mono_mass"]
            if mm is not None and abs(mm - float(mass)) <= win:
                put(base, 85.0 - abs(mm - float(mass)) / max(win, 1e-9) * 5.0, "mass")

    qlist = [q.strip() for q in ([query] if query else []) + list(queries or []) if q and q.strip()]
    for q in qlist:
        ql = q.lower()
        qtok, qgram = _tokens(q), _grams(q)
        for base, rec in meta.items():
            best, kind = 0.0, ""
            if base.lower() == ql or any(s.lower() == ql for s in rec["mams"]):
                best, kind = 96.0, "id-exact"
            elif ql and ql in base.lower():
                best, kind = 60.0 + 20.0 * len(ql) / max(len(base), 1), "id-substring"
            for nm in rec["names"]:
                nl = nm.lower()
                if nl == ql:
                    best, kind = max(best, 100.0), "name-exact"
                    continue
                if ql and (ql in nl or nl in ql):
                    # length-ratio weighted: a SHORT generic name sitting inside a long
                    # species query (…"sn-glycero-3-phosphocholine") is a weak, class-level
                    # match and must not outrank a near-identical species name.
                    ratio = min(len(ql), len(nl)) / max(len(ql), len(nl), 1)
                    s = 50.0 + 35.0 * ratio
                    if s > best:
                        best, kind = s, "name-substring"
                ntok = _tokens(nm)
                if qtok and ntok:
                    j = len(qtok & ntok) / len(qtok | ntok)
                    s = 60.0 * j + (10.0 if qtok <= ntok else 0.0)
                    if s > best:
                        best, kind = s, "name-tokens"
                d = max(_dice(qgram, rec["grams"].get(nm, set())),
                        _dice(qgram, rec["grams_bare"].get(nm, set())))
                # a trigram Dice above 0.8 means the two strings are near-identical spellings
                # of one name ("1-palmitoylglycero-3-phosphocholine" vs
                # "1-Palmitoylglycerophosphocholine") — the strongest non-exact evidence there is
                s = 100.0 * d if d > 0.8 else 60.0 * d
                if s > best:
                    best, kind = s, "name-trigram"
            if best >= min_score:
                put(base, best, kind, matched=q)

    ranked = sorted(hits.values(), key=lambda h: (-h["score"], h["mam_base"]))[:max_hits]
    groups: dict[str, list[str]] = defaultdict(list)
    for h in ranked:
        if h["canon_formula"]:
            groups[h["canon_formula"]].append(h["mam_base"])
    same_formula = {f: b for f, b in groups.items() if len(b) > 1}
    return {"model": idx["path"], "n_hits": len(ranked), "hits": ranked,
            "same_formula_groups": same_formula,
            "note": "evidence only — pick the species yourself and record it with gem_assign. "
                    "same_formula_groups are formula-identical isomer candidates: a formula or "
                    "mass check CANNOT separate them, so disambiguate by name/biology."}


def species_detail(mam: str, model_path: str | None = None) -> dict:
    """Look up one base id (or full species id) in the model. {found, ...record}."""
    idx = _index(model_path or GEM_MODEL)
    meta = idx["meta"]
    if mam in meta:
        return {"found": True, **_hit(mam, meta[mam], 100.0, "base-id")}
    for base, rec in meta.items():
        if mam in rec["mams"]:
            return {"found": True, **_hit(base, rec, 100.0, "species-id")}
    return {"found": False, "mam_base": mam, "model": idx["path"]}


def model_stats(model_path: str | None = None) -> dict:
    idx = _index(model_path or GEM_MODEL)
    return {k: idx[k] for k in ("n_met", "n_rxn", "n_base")}


def model_label(model_path: str | None = None) -> str:
    """Short label for a model path, used to key per-model results in the ledger."""
    from pathlib import Path
    return Path(model_path or GEM_MODEL).name.split(".")[0]
