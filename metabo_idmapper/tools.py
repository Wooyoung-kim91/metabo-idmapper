"""Tool registry — deterministic operations that EMIT evidence.

Every tool takes a `workdir` (the session) unless it is a pure calculator. Tools attach
candidates / verification / crosswalks to the ledger and return a small JSON summary.
The identity CALL is made by the LLM via `record_decision`, which is the only tool that
sets an entry's final_class/confidence — and it refuses IDs no other tool produced.
"""

from __future__ import annotations

from pathlib import Path

from . import guidance as _guidance
from .config import BRIDGE_DB, GEM_MODEL
from .engine import chebi as _chebi
from .engine import gem as _gem
from .engine import lexicon as _lex
from .engine import normalize as _norm
from .engine import structure as _struct
from .engine import verify as _verify
from .engine.rcall import run_r
from .state import ALL_CLASSES, PRIMARY_CLASSES, Entry, Session

# BridgeDb system names for the ids we handle.
_SRC = {"kegg": "KEGG Compound", "hmdb": "HMDB", "chebi": "ChEBI",
        "pubchem": "PubChem-compound", "inchikey": "InChIKey"}


# ----------------------------------------------------------------------- guidance/state
def midmap_guidance() -> dict:
    """Return the canonical workflow + operating rules. Call once at the start."""
    return {"version": "0.0.1", "guidance": _guidance.CANONICAL}


def detect_state(workdir: str) -> dict:
    """Report what is mapped vs pending and suggest the next tool(s)."""
    s = Session(workdir)
    counts = s.class_counts()
    n = len(s.entries)
    pend = s.pending_ids()
    if n == 0:
        nxt = ["ingest_names"]
    elif not any(e.get("candidates") for e in s.entries.values()):
        nxt = ["exact_match"]
    elif pend:
        nxt = ["structure_lookup", "search_synonym", "bridge_xref",
               "verify_candidate", "record_decision"]
    else:
        nxt = ["gem_crosswalk", "coverage_summary"]
    return {"workdir": str(s.workdir), "n_entries": n, "counts": counts,
            "pending": len(pend), "pending_sample": pend[:15],
            "suggested_next_tools": nxt}


def ingest_names(workdir: str, names: list[str] | None = None,
                 xlsx: str | None = None, sheet: str | int = 0,
                 column: str | None = None) -> dict:
    """Load raw names (list OR an xlsx column), normalize, and seed the ledger.

    Flags entries the LLM should look at: parenthetical abbreviations, combined names,
    and isomer-sensitive digit locants (never auto-substitute those).
    """
    s = Session(workdir)
    if xlsx:
        import pandas as pd
        df = pd.read_excel(xlsx, sheet_name=sheet)
        col = column or df.columns[0]
        names = [str(x) for x in df[col].dropna().tolist()]
        s.data["source"] = {"xlsx": str(xlsx), "sheet": sheet, "column": col}
    if not names:
        return {"error": "no names provided (pass names=[...] or xlsx=...&column=...)"}

    flagged = []
    xeno = []
    start = len(s.entries)
    for i, raw in enumerate(names):
        fid = f"M{start + i:04d}"
        norm = _norm.normalize(raw)
        e = Entry(feature_id=fid, original_name=raw, normalized=norm)
        tags = _lex.xenobiotic_hits(raw)
        if tags:
            e.decisions.append({"action": "screen", "rationale": "lexicon xenobiotic cue",
                                "tags": tags})
            xeno.append({"feature_id": fid, "name": raw, "tags": tags})
        s.upsert(e)
        if tags:
            s.entries[fid].setdefault("_flags", []).append("possible_xenobiotic")
        if norm["flags"]:
            flagged.append({"feature_id": fid, "name": raw,
                            "parenthetical": norm["parenthetical"],
                            "flags": norm["flags"]})
    s.save()
    return {"workdir": str(s.workdir), "ingested": len(names),
            "total_entries": len(s.entries), "flagged_for_review": flagged[:40],
            "n_flagged": len(flagged), "possible_xenobiotic": xeno[:40],
            "n_possible_xenobiotic": len(xeno),
            "note": "flags/tags are cues, not decisions — verify before accepting. "
                    "possible_xenobiotic entries are contaminant-class name matches."}


# ----------------------------------------------------------------------- evidence tools
def exact_match(workdir: str, feature_ids: list[str] | None = None,
                auto_accept_exact: bool = True) -> dict:
    """MetaboAnalyst batch name->ID. Exact DB-name matches are a lookup, not a judgement:
    with auto_accept_exact they are committed at confidence M1 (KEGG/HMDB/structure by
    what resolves); everything non-exact stays pending for the LLM."""
    s = Session(workdir)
    ids = feature_ids or list(s.entries)
    todo = [(fid, s.entries[fid]) for fid in ids if fid in s.entries]
    # map normalized display name (fallback original)
    q_by_name: dict[str, list[str]] = {}
    for fid, e in todo:
        nm = e["normalized"].get("normalized") or e["original_name"]
        q_by_name.setdefault(nm, []).append(fid)
    work = str((Path(s.workdir) / "_metaboanalyst_work").resolve())
    r = run_r("metaboanalyst_map.R", {"names": list(q_by_name), "workdir": work})
    if not r["ok"]:
        return {"error": "MetaboAnalyst failed", **r}
    accepted = matched = 0
    flagged: list[dict] = []
    for row in r["result"]:
        fids = q_by_name.get(row.get("query"), [])
        cand = {"source": "metaboanalyst", "query": row.get("query"),
                "kegg": row.get("kegg"), "hmdb": row.get("hmdb"),
                "chebi": row.get("chebi"), "pubchem": row.get("pubchem"),
                "note": row.get("match")}
        cand = {k: v for k, v in cand.items() if v not in (None, "null", "NA", "")}
        cand["source"] = "metaboanalyst"
        is_exact = row.get("status") == "matched"
        for fid in fids:
            s.add_candidate(fid, dict(cand))
            if is_exact:
                matched += 1
                if auto_accept_exact and not s.entries[fid].get("final_class"):
                    _commit(s, fid, accepted={k: cand[k] for k in
                            ("kegg", "hmdb", "chebi", "pubchem") if k in cand},
                            confidence="M1", rationale="exact MetaboAnalyst DB name match",
                            auto=True)
                    accepted += 1
                    # Flag risky auto-accepts: messy original name (abbrev/combined/trade-
                    # name-like single token) or a Match that diverges from the query.
                    e = s.entries[fid]
                    nm = e["normalized"]
                    risky = bool(set(nm.get("flags", [])) & {"abbrev_parenthetical", "combined_name"})
                    tok = (nm.get("normalized") or "").split()
                    if len(tok) == 1 and tok and tok[0][:1].isupper() and tok[0].isalpha():
                        risky = True  # single capitalized token → possible trade/brand name
                    if risky:
                        flagged.append({"feature_id": fid, "name": e["original_name"],
                                        "matched_as": row.get("match"),
                                        "kegg": cand.get("kegg"), "hmdb": cand.get("hmdb")})
                        e.setdefault("_flags", []).append("auto_accept_review")
    s.save()
    return {"queried": len(q_by_name), "matched": matched,
            "auto_accepted_M1": accepted, "flagged_auto_accepts": flagged,
            "n_flagged_auto_accepts": len(flagged), "counts": s.class_counts(),
            "note": "non-exact entries remain pending; gather evidence then record_decision. "
                    "flagged_auto_accepts have messy/ambiguous names — send them to the reviewer."}


def structure_lookup(workdir: str, feature_id: str, name: str | None = None) -> dict:
    """PubChem: name -> CID / InChIKey / formula / monoisotopic mass. Attaches evidence."""
    s = Session(workdir)
    e = s.get(feature_id)
    if e is None:
        return {"error": f"unknown feature_id {feature_id}"}
    query = name or e["normalized"].get("normalized") or e["original_name"]
    res = _struct.lookup(query)
    if res["found"]:
        s.add_candidate(feature_id, {
            "source": "pubchem", "query": query, "pubchem": str(res["cid"]),
            "inchikey": res["inchikey"], "formula": res["formula"],
            "mass": res["mono_mass"]})
        s.save()
        return {"feature_id": feature_id, "query": query, **res}
    # Not found: hand the driver actionable cues instead of a silent miss.
    alts = e["normalized"].get("alternatives", [])
    flags = e["normalized"].get("flags", [])
    hint = ("no PubChem hit for this string — likely a typo, abbreviation, or spacing "
            "variant. Retry structure_lookup with a corrected name=, or use search_synonym "
            "with your proposed expansion. Try the alternatives below first.")
    return {"feature_id": feature_id, "query": query, **res,
            "hint": hint, "alternatives": alts, "flags": flags}


def search_synonym(workdir: str, feature_id: str, queries: list[str],
                   dbs: list[str] | None = None, max_per_term: int = 6) -> dict:
    """Search KEGG / PubChem / ChEBI for LLM-proposed query strings (abbrev expansions,
    typo fixes, synonyms). Returns candidates; the LLM picks + verifies. Never auto-accepts."""
    s = Session(workdir)
    e = s.get(feature_id)
    if e is None:
        return {"error": f"unknown feature_id {feature_id}"}
    dbs = dbs or ["kegg", "pubchem", "chebi"]
    found = {"kegg": [], "pubchem": [], "chebi": []}
    if "kegg" in dbs:
        r = run_r("kegg_search.R", {"terms": queries, "max_per_term": max_per_term})
        if r["ok"]:
            for row in r["result"]:
                found["kegg"].append(row)
                s.add_candidate(feature_id, {"source": "kegg-search",
                    "query": row["term"], "kegg": row["kegg"], "note": row["kegg_name"]})
    if "pubchem" in dbs:
        for q in queries:
            hit = _struct.lookup(q)
            if hit["found"]:
                found["pubchem"].append({"query": q, "cid": hit["cid"],
                    "inchikey": hit["inchikey"], "formula": hit["formula"]})
                s.add_candidate(feature_id, {"source": "pubchem-search", "query": q,
                    "pubchem": str(hit["cid"]), "inchikey": hit["inchikey"],
                    "formula": hit["formula"], "mass": hit["mono_mass"]})
    if "chebi" in dbs:
        for q in queries:
            for hit in _chebi.search(q, rows=max_per_term):
                found["chebi"].append({"query": q, **hit})
                if hit["exact"]:
                    s.add_candidate(feature_id, {"source": "chebi-ols4", "query": q,
                        "chebi": hit["chebi"], "note": hit["label"]})
    s.save()
    return {"feature_id": feature_id, "queries": queries, "candidates": found,
            "note": "verify_candidate (formula/mass) before record_decision."}


def bridge_xref(workdir: str | None, queries: list[dict], db: str | None = None) -> dict:
    """Promote ids across systems via BridgeDb (batch; DB loads once).

    Each query: {id, source, targets:[...], feature_id?}. source/target in
    {kegg,hmdb,chebi,pubchem,inchikey}. If feature_id + workdir are given, bridged ids
    are attached to that entry as a 'bridgedb' candidate.
    """
    payload_q = []
    for q in queries:
        payload_q.append({
            "id": q["id"], "source": _SRC.get(q["source"], q["source"]),
            "targets": [_SRC.get(t, t) for t in q["targets"]]})
    r = run_r("bridge_xref.R", {"db": db or BRIDGE_DB, "queries": payload_q})
    if not r["ok"]:
        return {"error": "BridgeDb failed", **r}
    # normalize back to short keys
    inv = {v: k for k, v in _SRC.items()}
    out = []
    s = Session(workdir) if workdir else None
    for q, res in zip(queries, r["result"]):
        mp = {inv.get(t, t): v for t, v in res.get("mappings", {}).items()}
        out.append({"feature_id": q.get("feature_id"), "id": q["id"],
                    "source": q["source"], "mappings": mp})
        if s and q.get("feature_id") in (s.entries if s else {}):
            cand = {"source": "bridgedb", "query": f"{q['source']}:{q['id']}"}
            for k in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
                if mp.get(k):
                    v = mp[k][0] if isinstance(mp[k], list) else mp[k]
                    cand[k] = _lex.pad_hmdb(v) if k == "hmdb" else v
            s.add_candidate(q["feature_id"], cand)
    if s:
        s.save()
    return {"results": out, "db": db or BRIDGE_DB}


def verify_candidate(proposed_formula: str | None = None,
                     proposed_mass: float | None = None,
                     observed_formula: str | None = None,
                     observed_mass: float | None = None,
                     mass_tol_ppm: float = 15.0,
                     workdir: str | None = None, feature_id: str | None = None,
                     rationale: str | None = None) -> dict:
    """Deterministic formula/mass consistency gate (molmass). Passing this is required
    before accepting any fuzzy/synonym/mass-only candidate. Optionally logs to the ledger."""
    rep = _verify.verify(proposed_formula, proposed_mass, observed_formula,
                         observed_mass, mass_tol_ppm)
    if workdir and feature_id:
        s = Session(workdir)
        if s.get(feature_id) is not None:
            s.add_decision(feature_id, "verify", rationale or "formula/mass check", **rep)
            s.save()
    return {"feature_id": feature_id, **rep}


def mass_match_candidates(mz: float, adducts: list[str] | None = None,
                          tol_ppm: float = 15.0) -> dict:
    """m/z + adduct -> neutral monoisotopic mass windows (weak evidence). Hand the window
    to search_synonym / a DB search; never accept a mass-only hit as primary."""
    return {"mz": mz, "windows": _verify.mass_candidates(mz, adducts, tol_ppm)}


# ----------------------------------------------------------------------- the CALL
def _commit(s: Session, fid: str, accepted: dict, confidence: str,
            rationale: str, final_class: str | None = None, auto: bool = False) -> str:
    """Set an entry's accepted ids + final_class + confidence. Derives class from ids
    if not given. Enforces provenance: every accepted id must appear in a candidate."""
    e = s.entries[fid]
    accepted = {k: str(v) for k, v in accepted.items() if v}
    if accepted.get("hmdb"):
        accepted["hmdb"] = _lex.pad_hmdb(accepted["hmdb"])  # stable 7-digit form
    # anti-fabrication: accepted ids must be backed by a candidate
    seen = set()
    for c in e["candidates"]:
        for k in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
            if c.get(k):
                v = _lex.pad_hmdb(c[k]) if k == "hmdb" else str(c[k])
                seen.add((k, v))
    unbacked = [(k, v) for k, v in accepted.items() if (k, v) not in seen]
    if unbacked and not auto:
        raise ValueError(f"accepted ids not produced by any tool for {fid}: {unbacked} "
                         "(run structure_lookup/search_synonym/bridge_xref first)")
    if final_class is None:
        if accepted.get("kegg"):
            final_class = "KEGG-mapped"
        elif accepted.get("hmdb"):
            final_class = "HMDB-mapped"
        elif accepted.get("chebi") or accepted.get("pubchem") or accepted.get("inchikey"):
            final_class = "structure-only"
        else:
            final_class = "unmapped"
    e["accepted"] = accepted
    e["final_class"] = final_class
    e["confidence"] = confidence
    s.add_decision(fid, "auto_accept" if auto else "accept", rationale,
                   accepted=accepted, final_class=final_class, confidence=confidence)
    return final_class


def record_decision(workdir: str, feature_id: str, rationale: str,
                    accepted: dict | None = None, final_class: str | None = None,
                    confidence: str | None = None) -> dict:
    """Commit the LLM's identity CALL for one entry — the ONLY tool that sets
    final_class/confidence. For an exclusion, pass final_class='exogenous-excluded'
    (or 'unmapped') with accepted={} and a rationale."""
    s = Session(workdir)
    e = s.get(feature_id)
    if e is None:
        return {"error": f"unknown feature_id {feature_id}"}
    if final_class and final_class not in ALL_CLASSES:
        return {"error": f"final_class must be one of {sorted(ALL_CLASSES)}"}
    accepted = accepted or {}
    if final_class in {"exogenous-excluded", "unmapped"} and not accepted:
        e["final_class"] = final_class
        e["confidence"] = confidence or ("X" if final_class == "exogenous-excluded" else "U")
        s.add_decision(feature_id, "exclude", rationale, final_class=final_class)
        s.save()
        return {"feature_id": feature_id, "final_class": final_class,
                "counts": s.class_counts()}
    try:
        fc = _commit(s, feature_id, accepted, confidence or "M2", rationale, final_class)
    except ValueError as ex:
        return {"error": str(ex)}
    # Non-blocking warnings for the reasoning layer to reconsider.
    warnings = []
    if fc in PRIMARY_CLASSES or fc == "structure-only":
        tags = _lex.xenobiotic_hits(e["original_name"])
        if tags:
            warnings.append(f"name matches contaminant class {tags} but kept endogenous — "
                            "confirm it is a real metabolite, not a xenobiotic/additive.")
        if "has_digit_locant" in e["normalized"].get("flags", []) and (confidence or "M2") != "M1":
            warnings.append("isomer/locant-sensitive name accepted via a non-exact route — "
                            "verify the accepted structure has the SAME locant/anomer.")
    s.save()
    out = {"feature_id": feature_id, "final_class": fc, "accepted": e["accepted"],
           "confidence": e["confidence"], "counts": s.class_counts()}
    if warnings:
        out["warnings"] = warnings
    return out


def screen_exogenous(workdir: str) -> dict:
    """Deterministic contaminant-class screen (LC-MS additives / surfactants / plasticizers)
    over all entries. EMITS evidence — the endogenous-vs-xenobiotic CALL stays with the LLM.
    Surfaces entries currently kept endogenous whose name matches a contaminant class."""
    s = Session(workdir)
    hits, kept_endogenous = [], []
    for fid, e in s.entries.items():
        tags = _lex.xenobiotic_hits(e["original_name"])
        if not tags:
            continue
        row = {"feature_id": fid, "name": e["original_name"], "tags": tags,
               "final_class": e.get("final_class")}
        hits.append(row)
        if e.get("final_class") in PRIMARY_CLASSES or e.get("final_class") == "structure-only":
            kept_endogenous.append(row)
    return {"n_contaminant_class": len(hits), "contaminant_hits": hits,
            "kept_endogenous_review": kept_endogenous,
            "note": "kept_endogenous_review entries look like xenobiotics but are mapped as "
                    "metabolites — reconsider exogenous-excluded. Drugs/dietary compounds are "
                    "NOT in this lexicon; judge those yourself."}


def harness_audit(workdir: str) -> dict:
    """GOVERNANCE / process-completeness auditor (read-only, no identity judgement): reads the
    session ledger + workdir artifacts and checks that the reasoning layer actually HONORED
    metabo-idmapper's OWN contract — no fabricated ID, no mass-only (W) candidate used as
    primary, final_class↔confidence coherent, fuzzy (M2/M3) accepts carry formula/mass
    verification, locant/anomer-sensitive names re-checked, exogenous-exclusion applied
    consistently, flagged auto-accepts reviewed, id-gap KEGG re-tried, every CALL has a
    rationale, no entry left pending, gem_crosswalk ran, and the stage-7 always-emit artifacts
    exist. Emits a per-check pass/warn/fail scorecard so 'rules defined but not followed' is
    caught. Changes nothing. Run LAST (after coverage_summary)."""
    from . import harness as _harness
    return _harness.audit(workdir)


# ----------------------------------------------------------------------- GEM + export
def gem_crosswalk(workdir: str, feature_ids: list[str] | None = None,
                  model_path: str | None = None) -> dict:
    """Map accepted KEGG/HMDB/ChEBI ids to Mouse-GEM MAM species (flux input).
    Sets gem_mam + gem_cause ('id-gap' when it has ids but the GEM lacks the xref)."""
    s = Session(workdir)
    ids = feature_ids or [fid for fid, e in s.entries.items()
                          if e.get("final_class") in PRIMARY_CLASSES]
    mapped = 0
    id_gap_with_kegg: list[dict] = []
    stats = _gem.model_stats(model_path)
    for fid in ids:
        e = s.get(fid)
        if e is None:
            continue
        acc = e.get("accepted", {})
        res = _gem.crosswalk(kegg=acc.get("kegg"), hmdb=acc.get("hmdb"),
                             chebi=acc.get("chebi"), model_path=model_path)
        if res["mapped"]:
            e["gem_mam"] = res["mam_bases"]
            e["gem_cause"] = None
            mapped += 1
        else:
            e["gem_mam"] = []
            e["gem_cause"] = "id-gap"  # has ids but GEM has no xref for it
            # A KEGG-bearing id-gap is often an anomer/stereo-specific KEGG the GEM lacks
            # (e.g. C03251 vs generic C00092). Surface it so the driver retries generic KEGG.
            if acc.get("kegg"):
                id_gap_with_kegg.append({"feature_id": fid, "name": e["original_name"],
                                         "kegg": acc["kegg"]})
                e.setdefault("_flags", []).append("id_gap_try_generic_kegg")
    s.save()
    return {"model": model_path or GEM_MODEL, "model_stats": stats,
            "considered": len(ids), "gem_mapped": mapped,
            "id_gap_with_kegg": id_gap_with_kegg,
            "note": "unmapped gem_cause='id-gap'. For id_gap_with_kegg, the KEGG may be an "
                    "anomer/stereo-specific id the GEM lacks — search the GENERIC KEGG and "
                    "re-decide. Genuinely out-of-model compounds → reclassify 'model-scope-absent'."}


def backfill_hmdb(workdir: str, db: str | None = None) -> dict:
    """Bridge missing HMDB for every non-exogenous entry lacking an accepted HMDB, using its
    available ids (InChIKey/KEGG/ChEBI/PubChem) via BridgeDb in ONE batched call, then accept.

    Reflects the run-2 finding that HMDB rides along with KEGG matching and is under-counted;
    this maximizes HMDB coverage. Does not touch KEGG assignments or final_class of KEGG-mapped
    entries (HMDB is added as an extra xref); a structure-only entry that gains HMDB → HMDB-mapped.
    """
    s = Session(workdir)
    targets = [(fid, e) for fid, e in s.entries.items()
               if e.get("final_class") in ("KEGG-mapped", "structure-only")
               and not e.get("accepted", {}).get("hmdb")]
    queries = []
    for fid, e in targets:
        seen_src = {}
        for src in ("inchikey", "kegg", "chebi", "pubchem"):
            v = e.get("accepted", {}).get(src) or next(
                (c[src] for c in e.get("candidates", []) if c.get(src)), None)
            if v:
                seen_src[src] = v
        for src, idv in seen_src.items():
            queries.append({"feature_id": fid, "id": idv, "source": src, "targets": ["hmdb"]})
    if queries:
        bridge_xref(workdir, queries, db=db)  # attaches padded HMDB candidates
    s = Session(workdir)
    gained = []
    for fid, _ in targets:
        e = s.get(fid)
        hmdbs = [c["hmdb"] for c in e["candidates"] if c.get("hmdb")]
        if not hmdbs:
            continue
        acc = dict(e.get("accepted", {})); acc["hmdb"] = hmdbs[0]
        record_decision(workdir, fid, rationale="HMDB backfill (BridgeDb)",
                        accepted=acc, confidence=e.get("confidence") or "M2")
        gained.append({"feature_id": fid, "name": e["original_name"], "hmdb": hmdbs[0]})
        s = Session(workdir)
    return {"targets_missing_hmdb": len(targets), "hmdb_gained": len(gained),
            "gained": gained[:40], "counts": Session(workdir).class_counts(),
            "note": "run before coverage_summary to normalize HMDB coverage; structure-only "
                    "with no id-bridgeable HMDB remain (need name→HMDB, a separate step)."}


_ROUTE = {
    "kegg-search": "KEGG name/synonym search (keggFind) + verify",
    "pubchem-search": "PubChem name search",
    "chebi-ols4": "ChEBI OLS4 name search",
    "bridgedb": "xref bridge (BridgeDb)",
    "pubchem": "PubChem structure lookup",
    "metaboanalyst": "MetaboAnalyst exact",
}


def mapping_provenance(workdir: str, export: bool = True) -> dict:
    """Recovery provenance for the report: what was mapped to KEGG / HMDB BEYOND the
    MetaboAnalyst 1st pass and by which logic, plus the unmapped (harmonized-name-tried,
    reason). Writes kegg_recovered.tsv, hmdb_recovered.tsv, unmapped_harmonization.tsv."""
    import csv as _csv

    s = Session(workdir)

    def ma(e):
        return next((c for c in e.get("candidates", []) if c.get("source") == "metaboanalyst"), {})

    def cand_for(e, db, val):
        return next((c for c in e.get("candidates", []) if str(c.get(db)) == str(val)
                     and c.get("source") != "metaboanalyst"), None)

    def rationales(e):
        return [d.get("rationale", "") for d in e.get("decisions", [])
                if d.get("action") in ("accept", "auto_accept")]

    def harmon_name(e, c):
        # the query that matched, unless it is a bridge-style "src:id" -> use the compound name
        q = (c or {}).get("query", "")
        if not q or (":" in q and q.split(":", 1)[0] in
                     ("kegg", "hmdb", "chebi", "pubchem", "inchikey")):
            return e["normalized"].get("normalized", "")
        return q

    def exo_category(reason):
        r = (reason or "").lower()
        if any(w in r for w in ("antibiotic", "drug", "pharmac", "sulfonamide", "appetite",
                                "antihistamine", "alzheimer")):
            return "drug"
        if any(w in r for w in ("surfactant", "detergent", "diethanolamide", "sulfonic",
                                "amidopropyl", "alkyl sulfate")):
            return "surfactant"
        if any(w in r for w in ("plasticizer", "phthalate", "phosphate ester", "flame")):
            return "plasticizer"
        if any(w in r for w in ("additive", "trifluoroacetic", "mobile phase", "ion-pair", "ion pairing")):
            return "LC-MS additive"
        if any(w in r for w in ("antioxidant", "nitrophenol", "tert-butyl")):
            return "industrial antioxidant"
        if any(w in r for w in ("reagent", "building block", "industrial", "thiocyanate", "synthetic")):
            return "industrial/reagent"
        return "xenobiotic"

    kegg_rows, hmdb_rows, unmapped_rows, exo_rows = [], [], [], []
    for fid, e in s.entries.items():
        a = e.get("accepted", {})
        m = ma(e)
        if e.get("final_class") == "exogenous-excluded":
            tags = sorted({t for d in e.get("decisions", []) if d.get("action") == "screen"
                           for t in (d.get("tags") or [])})
            excl = [d.get("rationale", "") for d in e.get("decisions", []) if d.get("action") == "exclude"]
            if not excl:
                excl = [d.get("rationale", "") for d in e.get("decisions", [])
                        if "exclud" in d.get("rationale", "").lower() or "contaminant" in d.get("rationale", "").lower()]
            reason = excl[-1] if excl else (e.get("decisions", [{}])[-1].get("rationale", ""))
            ids = ";".join(f"{k}:{a[k]}" for k in ("kegg", "hmdb", "chebi", "pubchem") if a.get(k))
            category = ";".join(tags) if tags else exo_category(reason)
            exo_rows.append([fid, e["original_name"], category, ids, reason])
        if a.get("kegg") and not m.get("kegg"):
            c = cand_for(e, "kegg", a["kegg"])
            route = _ROUTE.get(c["source"], c["source"]) if c else "?"
            harmon = harmon_name(e, c)
            logic = next((r for r in rationales(e) if "re-bridge" not in r and "backfill" not in r), "")
            kegg_rows.append([fid, e["original_name"], harmon, a["kegg"],
                              (c or {}).get("note", ""), route, logic])
        if a.get("hmdb") and not m.get("hmdb"):
            c = cand_for(e, "hmdb", a["hmdb"])
            route = _ROUTE.get(c["source"], c["source"]) if c else "?"
            harmon = harmon_name(e, c)
            logic = next((r for r in rationales(e) if "hmdb" in r.lower() or "re-bridge" in r
                          or "backfill" in r.lower()), rationales(e)[0] if rationales(e) else "")
            hmdb_rows.append([fid, e["original_name"], harmon, a["hmdb"], route, logic])
        if e.get("final_class") in ("structure-only", "unmapped"):
            tried = sorted({c.get("query") for c in e.get("candidates", []) if c.get("query")})
            sid = ";".join(f"{k}:{a[k]}" for k in ("chebi", "pubchem", "inchikey") if a.get(k))
            reason = next((r for r in rationales(e)), "") or "no KEGG/HMDB returned by any source"
            unmapped_rows.append([fid, e["original_name"], " | ".join(tried),
                                  e.get("final_class"), sid, reason])

    tables = {
        "kegg_recovered.tsv": (["feature_id", "original_name", "harmonized_name", "kegg",
                                "kegg_name", "route", "logic"], kegg_rows),
        "hmdb_recovered.tsv": (["feature_id", "original_name", "harmonized_name", "hmdb",
                                "route", "logic"], hmdb_rows),
        "unmapped_harmonization.tsv": (["feature_id", "original_name", "harmonized_names_tried",
                                        "final_class", "structure_ids", "reason"], unmapped_rows),
        "exogenous_excluded.tsv": (["feature_id", "original_name", "category", "ids", "reason"],
                                   exo_rows),
    }
    paths = []
    if export:
        for fname, (cols, rows) in tables.items():
            p = Session(workdir).workdir / fname
            with open(p, "w", newline="") as f:
                w = _csv.writer(f, delimiter="\t")
                w.writerow(cols)
                w.writerows(rows)
            paths.append(str(p))
    return {"kegg_recovered": len(kegg_rows), "hmdb_recovered": len(hmdb_rows),
            "unmapped": len(unmapped_rows), "exogenous": len(exo_rows), "exports": paths,
            "note": "recovered = mapped BEYOND MetaboAnalyst 1st pass; route/logic per row."}


def annotate_source(workdir: str, source: str | None = None, sheet: "str | int" = 0,
                    header: int = 0, name_column: str | None = None,
                    out: str | None = None) -> dict:
    """Append the final ID-result columns to the ORIGINAL data file — the intensity matrix
    keeps all its columns and gains ID_final_class / ID_confidence / ID_kegg / ID_hmdb /
    ID_chebi / ID_pubchem / ID_inchikey / ID_gem_mam per compound (matched by name).

    Uses the source recorded by ingest_names if `source` is omitted. `header` is the 0-based
    header row (some vendor sheets have a preamble row, e.g. header=1). Writes <source>_annotated
    .xlsx + .tsv into the workdir."""
    import pandas as pd

    s = Session(workdir)
    rec = s.data.get("source", {})
    source = source or rec.get("xlsx")
    if sheet == 0 and rec.get("sheet") is not None:
        sheet = rec["sheet"]
    name_column = name_column or rec.get("column")
    if not source:
        return {"error": "no source file (pass source=... ; ingest_names records it when xlsx= used)"}
    df = pd.read_excel(source, sheet_name=sheet, header=header)
    col = name_column
    if col is None or col not in df.columns:
        # pick the column whose values best overlap the ledger names
        names = {e["original_name"].strip().lower() for e in s.entries.values()}
        best, bestn = df.columns[0], -1
        for c in df.columns:
            hit = sum(str(v).strip().lower() in names for v in df[c].dropna())
            if hit > bestn:
                best, bestn = c, hit
        col = best
    idmap = {e["original_name"].strip().lower(): e for e in s.entries.values()}
    add = ["final_class", "confidence", "kegg", "hmdb", "chebi", "pubchem", "inchikey", "gem_mam"]
    for a in add:
        df["ID_" + a] = ""
    matched = 0
    for i, v in df[col].items():
        e = idmap.get(str(v).strip().lower())
        if not e:
            continue
        matched += 1
        acc = e.get("accepted", {})
        df.at[i, "ID_final_class"] = e.get("final_class", "")
        df.at[i, "ID_confidence"] = e.get("confidence", "")
        for k in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
            df.at[i, "ID_" + k] = acc.get(k, "")
        df.at[i, "ID_gem_mam"] = ";".join(e.get("gem_mam", []))
    stem = str(Path(source).name).rsplit(".", 1)[0]
    out = out or str(Session(workdir).workdir / f"{stem}_annotated.xlsx")
    df.to_excel(out, index=False)
    tsv = out.rsplit(".", 1)[0] + ".tsv"
    df.to_csv(tsv, sep="\t", index=False)
    return {"source": source, "name_column": col, "rows": len(df), "matched": matched,
            "columns_added": ["ID_" + a for a in add], "annotated": [out, tsv],
            "note": "intensity matrix + final ID columns; unmatched rows have empty ID_* fields."}


def export_code(workdir: str) -> dict:
    """Write standalone reproduction code that reproduces the run using the ORIGINAL library
    APIs (MetaboAnalystR/BridgeDbR/KEGGREST via Rscript, PubChem PUG-REST, molmass, COBRApy,
    matplotlib) — NOT the tool wrappers. Emits both code/reproduce_mapping.py (flow-based,
    reads the saved ledger) and code/reproduce_mapping.ipynb (same flow unrolled into linear
    cells, no def, with detailed input/output/reuse comments per cell)."""
    from . import scriptgen
    try:
        script = scriptgen.generate(workdir)
        notebook = scriptgen.generate_notebook(workdir)
    except Exception as e:
        return {"error": f"scriptgen failed: {type(e).__name__}: {e}"}
    return {"script": script, "notebook": notebook,
            "note": "raw-API reproduction; no metabo_idmapper import."}


def export_report_ppt(workdir: str, out: str | None = None) -> dict:
    """Build a PPTX report from the run's outputs + figures (Title · Coverage KPI · Methods ·
    Pipeline · UpSet · Improvement · Recovery cause→fix · KEGG/HMDB recovered · Unmapped ·
    Exogenous · Outputs). Reads coverage_summary.tsv + provenance tsvs + figures; run
    coverage_summary first. Requires python-pptx + pillow."""
    from . import report_ppt
    try:
        path = report_ppt.generate(workdir, out=out)
    except Exception as e:
        return {"error": f"ppt export failed: {type(e).__name__}: {e}"}
    return {"pptx": path, "note": "values read from run artifacts; no fabrication."}


def plot_coverage(workdir: str) -> dict:
    """Always-on DB-matching figures: figures/db_matching_upset.png (5-DB UpSet + enriched_xref.tsv)
    and figures/db_matching_improvement.png (MetaboAnalyst baseline vs current logic)."""
    from .engine import plots as _plots
    try:
        up = _plots.upset(workdir)
        imp = _plots.improvement(workdir)
    except Exception as e:  # never let a plotting error break the run
        return {"error": f"plot failed: {type(e).__name__}: {e}"}
    return {"upset": up, "improvement": imp}


def coverage_summary(workdir: str, export: bool = True, figures: bool = True) -> dict:
    """Compute class/confidence coverage and (optionally) write the master ledger tsv,
    coverage tsv, and provenance tsv. With figures=True (default) ALWAYS emits the
    db_matching_upset + db_matching_improvement figures + enriched_xref.tsv."""
    import csv

    s = Session(workdir)
    counts = s.class_counts()
    n = len(s.entries)
    primary = sum(counts[c] for c in PRIMARY_CLASSES)
    has_kegg = sum(1 for e in s.entries.values() if e.get("accepted", {}).get("kegg"))
    has_hmdb = sum(1 for e in s.entries.values() if e.get("accepted", {}).get("hmdb"))
    gem_mapped = sum(1 for e in s.entries.values() if e.get("gem_mam"))
    summary = {"total": n, "counts": counts,
               "primary_usable": primary, "has_kegg": has_kegg, "has_hmdb": has_hmdb,
               "gem_mapped": gem_mapped}
    if export:
        led = s.workdir / "master_ledger.tsv"
        cols = ["feature_id", "original_name", "normalized", "final_class", "confidence",
                "kegg", "hmdb", "chebi", "pubchem", "inchikey", "gem_mam", "gem_cause",
                "n_candidates", "n_decisions"]
        with open(led, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(cols)
            for fid, e in s.entries.items():
                acc = e.get("accepted", {})
                w.writerow([fid, e["original_name"],
                            e["normalized"].get("normalized", ""),
                            e.get("final_class", ""), e.get("confidence", ""),
                            acc.get("kegg", ""), acc.get("hmdb", ""), acc.get("chebi", ""),
                            acc.get("pubchem", ""), acc.get("inchikey", ""),
                            ";".join(e.get("gem_mam", [])), e.get("gem_cause") or "",
                            len(e.get("candidates", [])), len(e.get("decisions", []))])
        cov = s.workdir / "coverage_summary.tsv"
        with open(cov, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["metric", "value", "denom", "pct"])
            for c in sorted(ALL_CLASSES):
                w.writerow([f"class:{c}", counts[c], n,
                            round(100 * counts[c] / n, 1) if n else 0])
            for name, val in (("primary_usable", primary), ("has_kegg", has_kegg),
                              ("has_hmdb", has_hmdb), ("gem_mapped", gem_mapped)):
                w.writerow([name, val, n, round(100 * val / n, 1) if n else 0])
        summary["exports"] = [str(led), str(cov)]
    if figures:
        summary["figures"] = plot_coverage(workdir)
    if export:
        summary["provenance"] = mapping_provenance(workdir)
        summary["code"] = export_code(workdir)
        if Session(workdir).data.get("source", {}).get("xlsx"):
            try:  # append final ID columns to the original file (best-effort)
                summary["annotated_source"] = annotate_source(workdir)
            except Exception as e:
                summary["annotated_source"] = {"error": f"{type(e).__name__}: {e}"}
    return summary


# registry consumed by mcp_server
REGISTRY = [
    midmap_guidance, detect_state, ingest_names, exact_match, structure_lookup,
    search_synonym, bridge_xref, verify_candidate, mass_match_candidates,
    screen_exogenous, record_decision, backfill_hmdb, gem_crosswalk,
    mapping_provenance, annotate_source, plot_coverage, export_code,
    export_report_ppt, coverage_summary, harness_audit,
]
