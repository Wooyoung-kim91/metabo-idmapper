"""Tool registry — deterministic operations that EMIT evidence.

Every tool takes a `workdir` (the session) unless it is a pure calculator. Tools attach
candidates / verification / crosswalks to the ledger and return a small JSON summary.
The identity CALL is made by the LLM via `record_decision`, which is the only tool that
sets an entry's final_class/confidence — and it refuses IDs no other tool produced.
"""

from __future__ import annotations

from pathlib import Path

from . import collide as _collide
from . import guidance as _guidance
from .config import BRIDGE_DB, GEM_MODEL
from .engine import chebi as _chebi
from .engine import entry as _entry
from .engine import gem as _gem
from .engine import isomer as _isomer
from .engine import lexicon as _lex
from .engine import normalize as _norm
from .engine import shorthand as _short
from .engine import structure as _struct
from .engine import verify as _verify
from .engine.rcall import run_r
from .state import (ALL_CLASSES, CLASS_FOR_ORIGIN, EXCLUDED_CLASSES, GEM_MAPPED_RELATIONS,
                    GEM_RELATIONS, ID_KEYS, LEGACY_CLASSES, PRIMARY_CLASSES, VALID_ORIGINS,
                    Entry, Session)

# BridgeDb system names for the ids we handle.
_SRC = {"kegg": "KEGG Compound", "hmdb": "HMDB", "chebi": "ChEBI",
        "pubchem": "PubChem-compound", "inchikey": "InChIKey"}


# ----------------------------------------------------------------------- guidance/state
def midmap_guidance() -> dict:
    """Return the canonical workflow + operating rules. Call once at the start."""
    return {"version": "0.2.0", "guidance": _guidance.CANONICAL}


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
               "verify_candidate", "id_name_check", "record_decision"]
    else:
        nxt = ["id_name_check", "collision_check", "gem_crosswalk", "gem_search",
               "gem_assign", "coverage_summary"]
    # entry-level flags the driver must clear before finalizing
    open_flags: dict[str, list[str]] = {}
    for fid, e in s.entries.items():
        for f in e.get("_flags", []) or []:
            open_flags.setdefault(f, []).append(fid)
    unchecked = [fid for fid, e in s.entries.items()
                 if (e.get("accepted") or {}).get("kegg")
                 and not _backchecked(e, "kegg", str(e["accepted"]["kegg"]))]
    gem_open = [fid for fid, e in s.entries.items()
                if e.get("final_class") in (PRIMARY_CLASSES | {"exogenous"})
                and (e.get("gem_relation") in (None, "id-gap"))]
    return {"workdir": str(s.workdir), "n_entries": n, "counts": counts,
            "pending": len(pend), "pending_sample": pend[:15],
            "flags": {k: len(v) for k, v in sorted(open_flags.items())},
            "flagged_sample": {k: v[:8] for k, v in sorted(open_flags.items())},
            "kegg_not_backchecked": len(unchecked),
            "gem_unresolved": len(gem_open), "gem_unresolved_sample": gem_open[:10],
            "suggested_next_tools": nxt,
            "note": "flags are open review items (id_name_conflict / shared_id_collision / "
                    "isomer_token_conflict / class_level_id / hmdb_backfill_conflict / "
                    "auto_accept_review / id_gap_try_generic_kegg) — resolve each before "
                    "finalizing. gem_unresolved entries still need gem_search + gem_assign "
                    "(including an explicit 'model-scope-absent' where that is the answer)."}


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
    isomer_flagged: list[dict] = []
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
                    # The DB's Match name is evidence in its own right: when it disagrees with
                    # the query on a discriminant axis (ether linkage, glycan series, sialyl
                    # linkage, chains) the "exact" match is an isomer swap or a class id
                    # answering a species query — neither is an M1 identity.
                    rep = _isomer.compare(nm.get("normalized") or e["original_name"],
                                          row.get("match") or "")
                    if rep["verdict"] in ("conflict", "class-level"):
                        risky = True
                        isomer_flagged.append({
                            "feature_id": fid, "name": e["original_name"],
                            "matched_as": row.get("match"), "verdict": rep["verdict"],
                            "conflicts": rep["conflicts"], "kegg": cand.get("kegg"),
                            "hmdb": cand.get("hmdb")})
                        flag = ("isomer_token_conflict" if rep["verdict"] == "conflict"
                                else "class_level_id")
                        if flag not in e.setdefault("_flags", []):
                            e["_flags"].append(flag)
                    if risky:
                        flagged.append({"feature_id": fid, "name": e["original_name"],
                                        "matched_as": row.get("match"),
                                        "kegg": cand.get("kegg"), "hmdb": cand.get("hmdb")})
                        if "auto_accept_review" not in e.setdefault("_flags", []):
                            e["_flags"].append("auto_accept_review")
    s.save()
    return {"queried": len(q_by_name), "matched": matched,
            "auto_accepted_M1": accepted, "flagged_auto_accepts": flagged,
            "n_flagged_auto_accepts": len(flagged),
            "isomer_or_class_mismatch": isomer_flagged,
            "n_isomer_or_class_mismatch": len(isomer_flagged), "counts": s.class_counts(),
            "note": "non-exact entries remain pending; gather evidence then record_decision. "
                    "flagged_auto_accepts have messy/ambiguous names — send them to the "
                    "reviewer. isomer_or_class_mismatch entries were auto-accepted but the DB's "
                    "own Match name disagrees on a discriminant axis (or is a class entry for a "
                    "species name): re-decide those, and run id_name_check on their KEGG."}


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


def search_synonym(workdir: str, feature_id: str, queries: list[str] | None = None,
                   dbs: list[str] | None = None, max_per_term: int = 6,
                   expand_shorthand: bool = False) -> dict:
    """Search KEGG / PubChem / ChEBI for LLM-proposed query strings (abbrev expansions,
    typo fixes, synonyms). Returns candidates; the LLM picks + verifies. Never auto-accepts.

    Every hit that carries a name is screened with `isomer.compare` against the entry's name,
    so a keggFind hit that is the WRONG isomer arrives labelled as such — the verified case
    being nLc4Cer, where C04910's own name reads "…1,3-beta-D-Galactosyl…; Lc4Cer", i.e. the
    beta1-3 LACTO isomer of a neolacto query. Hits are returned safest-verdict first.

    With `expand_shorthand` the deterministic lipid/glycolipid shorthand variants of the entry
    name are searched too (LysoPC(16:0) -> 1-palmitoyl-sn-glycero-3-phosphocholine, ...).
    """
    s = Session(workdir)
    e = s.get(feature_id)
    if e is None:
        return {"error": f"unknown feature_id {feature_id}"}
    queries = list(queries or [])
    if expand_shorthand:
        base = _entry_query_name(e)
        queries += [base] + _short.variants(base) + _short.variants(e["original_name"])
    queries = list(dict.fromkeys(q for q in queries if q and q.strip()))
    if not queries:
        return {"error": "pass queries=[...] (or expand_shorthand=True to derive them)"}
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
    qname = _entry_query_name(e)
    screened = {
        "kegg": _isomer.screen(qname, found["kegg"], name_key="kegg_name"),
        "pubchem": found["pubchem"],
        "chebi": _isomer.screen(qname, found["chebi"], name_key="label"),
    }
    n_conflict = sum(1 for db in ("kegg", "chebi")
                     for h in screened[db] if h.get("isomer_verdict") == "conflict")
    return {"feature_id": feature_id, "queries": queries, "candidates": screened,
            "n_isomer_conflicts": n_conflict,
            "note": "candidates are sorted safest-verdict first and carry isomer_verdict: "
                    "'conflict' = the hit's own DB name disagrees with this name on a "
                    "discriminant axis (wrong isomer/compound) — do not accept it; "
                    "'class-level' = generic class entry for a species name. Then "
                    "verify_candidate (formula/mass) and id_name_check before record_decision."}


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


# --------------------------------------------------------- back-check + isomer + collisions
def _entry_query_name(e: dict) -> str:
    return e["normalized"].get("normalized") or e["original_name"]


def _best_verdict(query: str, db_names: list[str]) -> dict:
    """Verdict of a query name against ALL synonyms of one DB entry. A CONFLICT WINS.

    A DB entry lists many synonyms and they are not always consistent with each other. KEGG
    C04910 lists both "Lc4Cer" (the beta1-3 LACTO isomer) and "Paragloboside" (which IS
    neolacto) on the same entry: a rule that let a matching synonym win would report a
    neolacto query as `ok` against the lacto entry — the exact silent isomer swap this tool
    exists to prevent. So any conflicting synonym decides the verdict, and both sides are
    returned (`conflicting_names` / `matching_names`) for the driver to read.
    """
    reps = [_isomer.compare(query, nm) for nm in db_names if nm]
    if not reps:
        return {"verdict": "unknown", "matched_name": None, "conflicts": [],
                "note": "DB entry returned no name to compare"}
    bad = [r for r in reps if r["verdict"] == "conflict"]
    good = [r for r in reps if r["verdict"] == "ok"]
    lead = bad[0] if bad else (good[0] if good else None)
    if lead is None:
        for want in ("class-level", "ambiguous"):
            lead = next((r for r in reps if r["verdict"] == want), None)
            if lead:
                break
    lead = lead or reps[0]
    out = {"verdict": lead["verdict"], "matched_name": lead["candidate"],
           "conflicts": lead["conflicts"], "ambiguities": lead.get("ambiguities", []),
           "candidate_resolution": lead["candidate_resolution"],
           "conflicting_names": [r["candidate"] for r in bad],
           "matching_names": [r["candidate"] for r in good]}
    if bad and good:
        out["note"] = ("this DB entry's OWN synonyms disagree: some name this compound, others "
                       "name a different isomer. Decide from the systematic name/linkage, not "
                       "from the synonym list.")
    return out


def id_name_check(workdir: str | None = None, feature_ids: list[str] | None = None,
                  dbs: list[str] | None = None, ids: list[dict] | None = None,
                  only_accepted: bool = True) -> dict:
    """BACK-CHECK: fetch what each accepted ID *itself* resolves to and compare it to the name.

    Every other route here runs name -> id and therefore cannot detect a plausible-but-WRONG
    id, because the wrong id arrives attached to the right name (a bad BridgeDb link, a
    PubChem entry merging two compounds as synonyms, a class id answering a species query).
    This tool runs the other direction: id -> the DB's own NAME/FORMULA record (KEGG keggGet,
    batched; ChEBI OLS4; PubChem CID) -> `isomer.compare` against the entry's name.

    In a verified run this is what caught a taurochenodeoxycholate entry carrying C05472,
    whose own KEGG name is "Urocortisol / Tetrahydrocortisol" — a cortisol metabolite, not a
    bile acid. Mass/formula checks cannot catch that class of error and isomers are formula-
    identical, so run this on every KEGG accepted via a bridge or a search route.

    dbs defaults to ["kegg"] (batched, cheap). Pass ["kegg","chebi","pubchem"] for a full
    pass. HMDB has no open record API — bridge it to PubChem/InChIKey and check that instead.
    Records a `backcheck` decision per entry and flags conflicts; makes no CALL of its own.
    """
    dbs = [d.lower() for d in (dbs or ["kegg"])]

    # --- ad-hoc mode: check bare {db, id, name} triples without a session ---
    if ids:
        out = []
        by_db: dict[str, list[str]] = {}
        for q in ids:
            by_db.setdefault(q["db"].lower(), []).append(str(q["id"]))
        recs = {db: _entry.fetch(db, vals) for db, vals in by_db.items()}
        for q in ids:
            rec = recs.get(q["db"].lower(), {}).get(str(q["id"]), {})
            row = {"db": q["db"], "id": q["id"], "found": rec.get("found"),
                   "db_names": rec.get("names", [])[:8], "formula": rec.get("formula")}
            if q.get("name"):
                row.update(_best_verdict(q["name"], rec.get("names", [])))
            out.append(row)
        return {"checked": len(out), "results": out,
                "note": "ad-hoc back-check (no ledger written). verdict 'conflict' means the "
                        "id resolves to a different compound/isomer than the name."}

    if not workdir:
        return {"error": "pass workdir (+ optional feature_ids) or ids=[{db,id,name}]"}
    s = Session(workdir)
    targets = [(fid, s.entries[fid]) for fid in (feature_ids or list(s.entries))
               if fid in s.entries]

    # collect ids per db, then fetch each db ONCE (keggGet is batched 10 per request)
    want: dict[str, set[str]] = {db: set() for db in dbs}
    per_entry: dict[str, list[tuple[str, str]]] = {}
    for fid, e in targets:
        pairs: list[tuple[str, str]] = [(db, str(v)) for db, v in e.get("accepted", {}).items()
                                        if db in dbs and v]
        if not only_accepted:
            # check EVERY candidate id too, not one arbitrary representative per db
            pairs += [(db, str(c[db])) for c in e.get("candidates", []) for db in dbs
                      if c.get(db)]
        pairs = list(dict.fromkeys(pairs))
        if pairs:
            per_entry[fid] = pairs
            for db, val in pairs:
                want[db].add(val)
    recs = {db: _entry.fetch(db, sorted(vals)) for db, vals in want.items() if vals}

    rows, conflicts, fetch_errors = [], [], []
    for fid, pairs in per_entry.items():
        e = s.entries[fid]
        qname = _entry_query_name(e)
        for db, val in pairs:
            rec = recs.get(db, {}).get(val, {})
            if not rec.get("found"):
                fetch_errors.append({"feature_id": fid, "db": db, "id": val,
                                     "reason": rec.get("error") or rec.get("note")
                                     or "not found in DB"})
                continue
            rep = _best_verdict(qname, rec.get("names", []))
            fm = None
            obs = next((c.get("formula") for c in e.get("candidates", []) if c.get("formula")),
                       None)
            if rec.get("formula") and obs:
                fm = _verify.verify(proposed_formula=rec["formula"], observed_formula=obs)
            row = {"feature_id": fid, "name": e["original_name"], "db": db, "id": val,
                   "db_names": rec.get("names", [])[:6], "db_formula": rec.get("formula"),
                   "verdict": rep["verdict"], "matched_name": rep.get("matched_name"),
                   "conflicts": rep.get("conflicts", []),
                   "db_resolution": rep.get("candidate_resolution"),
                   "formula_check": fm}
            rows.append(row)
            s.add_decision(fid, "backcheck",
                           f"{db}:{val} resolves to {rec.get('names', [''])[:2]} — "
                           f"isomer verdict {rep['verdict']}",
                           db=db, id=val, db_names=rec.get("names", [])[:6],
                           verdict=rep["verdict"], conflicts=rep.get("conflicts", []),
                           db_formula=rec.get("formula"))
            flags = e.setdefault("_flags", [])
            if rep["verdict"] == "conflict":
                conflicts.append(row)
                if "id_name_conflict" not in flags:
                    flags.append("id_name_conflict")
            elif rep["verdict"] == "class-level" and "class_level_id" not in flags:
                flags.append("class_level_id")
    s.save()
    return {"checked": len(rows), "n_conflicts": len(conflicts), "conflicts": conflicts,
            "results": rows, "fetch_errors": fetch_errors,
            "by_verdict": {v: sum(1 for r in rows if r["verdict"] == v)
                           for v in ("ok", "class-level", "ambiguous", "conflict", "unknown")},
            "note": "CONFLICT = the id's own DB record is a different compound/isomer than the "
                    "name: re-decide it (search the right id, or drop the id). class-level = "
                    "the id is a generic class entry for a species-level name: keep it only if "
                    "you report it AS class-level. ok = the DB entry states this compound."}


def isomer_guard(query: str, candidate_names: list[str], workdir: str | None = None,
                 feature_id: str | None = None) -> dict:
    """Compare a metabolite name against candidate DB/model entry NAMES on the discriminant
    axes — skeleton, ether linkage (plasmenyl P- / plasmanyl O- / acyl), glycan series
    (neolacto vs lacto) and length, sialyl linkage (a2,3 vs a2,6), chains, double-bond
    positions, omega, oxidation vs peroxidation, hydroxy/acetyl count, polyamine backbone —
    and report which axes DISAGREE. Pure string logic: no network, no lookups.

    Use it whenever a candidate's name is available but formula/mass cannot decide, which is
    the normal case for lipids and glycolipids (isomers are formula-identical). Also reports
    `resolution` per candidate: species-level vs class-level (generic acyl/alkyl/R-group).
    """
    ranked = _isomer.screen(query, [{"name": n} for n in candidate_names])
    if workdir and feature_id:
        s = Session(workdir)
        if s.get(feature_id) is not None:
            s.add_decision(feature_id, "isomer_guard",
                           f"name-token screen of {len(candidate_names)} candidate name(s)",
                           query=query,
                           verdicts={c["name"]: c["isomer_verdict"] for c in ranked})
            s.save()
    return {"query": query, "candidates": ranked,
            "n_conflict": sum(1 for c in ranked if c["isomer_verdict"] == "conflict"),
            "query_features": _isomer.features(query),
            "note": "conflict = different compound/isomer, do NOT accept as the same identity "
                    "(record it as a surrogate if you use it anyway). class-level = generic "
                    "class entry for a species query. ambiguous = names cannot decide; get "
                    "another evidence source."}


def collision_check(workdir: str, model_label: str | None = None) -> dict:
    """Find identifiers (and model species) shared by MORE THAN ONE entry.

    A shared id is silent and destructive: in a verified run Lc3Cer and nLc4Cer — the
    numerator and denominator of one marker ratio — both ended up on PubChem CID 131770449 /
    HMDB0062485 because PubChem lists the two names as synonyms of one entry, which would
    have collapsed that ratio to exactly 1. Nothing in an id-to-id comparison reveals it.

    Each shared id is classified from the run's own evidence (see `collide`): a real
    `collision` between different compounds, a legitimate `class-id-shared` (the DB has only a
    class-level entry, e.g. LysoPC 16:0/18:2/20:4 all on KEGG C04230 — keep it, but report it
    as class-level), or a harmless `duplicate-name`. Back-checking an id with `id_name_check` is
    what turns a suspected collision into a documented class-level share.

    Shared GEM species are reported too, but a shared CLASS-PROXY pool species (crm_hs for two
    ceramides) is expected and reported as `class_proxy_shared`, not an error. Read-only apart
    from raising review flags.
    """
    s = Session(workdir)
    gem_hits, proxy_hits = [], []
    id_hits = _collide.scan(s.entries, s.accepted_index())
    for row in id_hits:
        for f in row["feature_ids"]:
            flags = s.entries[f].setdefault("_flags", [])
            if row["severity"] == "collision" and "shared_id_collision" not in flags:
                flags.append("shared_id_collision")
    label = model_label
    for mam, fids in sorted(s.gem_index(label).items()):
        if len(fids) < 2:
            continue
        rels = {f: (s.entries[f].get("gem_relation") or "unspecified") for f in fids}
        row = {"mam_base": mam, "feature_ids": fids,
               "names": [s.entries[f]["original_name"] for f in fids], "relations": rels}
        (proxy_hits if set(rels.values()) <= {"class-proxy"} else gem_hits).append(row)
    s.save()
    by_sev: dict[str, list[dict]] = {}
    for r in id_hits:
        by_sev.setdefault(r["severity"], []).append(r)
    return {"n_shared_ids": len(id_hits),
            "n_collisions": len(by_sev.get("collision", [])),
            "collisions": by_sev.get("collision", []),
            "class_id_shared": by_sev.get("class-id-shared", []),
            "duplicate_names": by_sev.get("duplicate-name", []),
            "gem_species_shared": gem_hits, "class_proxy_shared": proxy_hits,
            "note": "collision = one id on two DIFFERENT compounds — split them (find the "
                    "species-specific id for each) before finalizing; if either is the "
                    "numerator/denominator of a ratio the ratio is destroyed. class_id_shared "
                    "= the DB only has a class-level entry, so the share is legitimate but the "
                    "entries must be reported as class-level. duplicate-name = the same "
                    "compound listed twice. class_proxy_shared = several species legitimately "
                    "sharing one model class/pool species. Run id_name_check first: a "
                    "back-checked class-level id downgrades a suspected collision."}


# ----------------------------------------------------------------------- the CALL
# candidate sources that are a SEARCH or a BRIDGE rather than a curated exact name match:
# an id from one of these has never been confirmed against its own DB record.
_INDIRECT_SOURCES = {"bridgedb", "kegg-search", "pubchem-search", "chebi-ols4", "pubchem"}


def _candidates_for(e: dict, db: str, val: str) -> list[dict]:
    v = _lex.pad_hmdb(str(val)) if db == "hmdb" else str(val)
    out = []
    for c in e.get("candidates", []):
        cv = _lex.pad_hmdb(str(c[db])) if (db == "hmdb" and c.get(db)) else (
            str(c[db]) if c.get(db) else None)
        if cv == v:
            out.append(c)
    return out


def _backchecked(e: dict, db: str, val: str) -> bool:
    return any(d.get("action") == "backcheck" and d.get("db") == db
               and str(d.get("id")) == str(val) for d in e.get("decisions", []))


def _accept_warnings(s: Session, fid: str, e: dict) -> list[str]:
    """Warn about the three silent-error routes an id-level check cannot see: an id shared
    with another entry, an id whose own DB name disagrees with this name, and a searched or
    bridged primary id that was never back-checked against its DB record."""
    warns: list[str] = []
    accepted = e.get("accepted") or {}
    index = s.accepted_index()
    for db, val in accepted.items():
        others = [f for f in index.get((db, str(val)), []) if f != fid]
        if others:
            names = [s.entries[o]["original_name"] for o in others]
            verdicts = {n: _isomer.compare(e["original_name"], n)["verdict"] for n in names}
            distinct = [n for n, v in verdicts.items() if v == "conflict"] or [
                n for n in names if n.strip().lower() != e["original_name"].strip().lower()]
            if distinct:
                warns.append(f"{db}:{val} is ALSO accepted for {distinct} — one id on two "
                             "different compounds destroys any ratio between them; find the "
                             "species-specific id for each (run collision_check).")
        for c in _candidates_for(e, db, str(val)):
            nm = c.get("note") or ""
            if not nm:
                continue
            rep = _isomer.compare(_entry_query_name(e), nm)
            if rep["verdict"] == "conflict":
                warns.append(f"{db}:{val} was returned as '{nm}', which conflicts with this "
                             f"name on {[x['axis'] for x in rep['conflicts']]} — wrong isomer/"
                             "compound unless you can explain the difference.")
            elif rep["verdict"] == "class-level":
                warns.append(f"{db}:{val} ('{nm}') is a CLASS-level entry for a species-level "
                             "name — keep it only if you report it as class-level.")
    # KEGG only: it is the id that drives pathway + model input, and the one DB here with a
    # record API. A bridged HMDB is checked indirectly (collision_check / bridge to PubChem).
    val = accepted.get("kegg")
    if val and not _backchecked(e, "kegg", str(val)):
        srcs = {c.get("source") for c in _candidates_for(e, "kegg", str(val))}
        if srcs & _INDIRECT_SOURCES:
            warns.append(f"kegg:{val} came from {sorted(srcs & _INDIRECT_SOURCES)} and its own "
                         "KEGG record was never checked — run id_name_check (this is how a "
                         "bile-acid entry ended up carrying a cortisol-metabolite KEGG).")
    return warns


def _commit(s: Session, fid: str, accepted: dict, confidence: str,
            rationale: str, final_class: str | None = None, auto: bool = False,
            origin: str | None = None) -> str:
    """Set an entry's accepted ids + final_class + confidence (+ origin). Derives class from
    ids if not given. Enforces provenance: every accepted id must appear in a candidate."""
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
    # a KEPT endogenous mapping with no explicit origin defaults to 'endogenous'
    if origin is None and final_class in (PRIMARY_CLASSES | {"structure-only"}):
        origin = "endogenous"
    e["accepted"] = accepted
    e["final_class"] = final_class
    e["confidence"] = confidence
    e["origin"] = origin
    s.add_decision(fid, "auto_accept" if auto else "accept", rationale,
                   accepted=accepted, final_class=final_class, confidence=confidence,
                   origin=origin)
    return final_class


def record_decision(workdir: str, feature_id: str, rationale: str,
                    accepted: dict | None = None, final_class: str | None = None,
                    confidence: str | None = None, origin: str | None = None) -> dict:
    """Commit the LLM's identity CALL for one entry — the ONLY tool that sets
    final_class/confidence/origin. `origin` is the compound's PROVENANCE and, for non-
    endogenous compounds, decides the class:

      - BIOLOGICAL but from outside the host -> origin in {diet, drug, microbial, plant};
        final_class='exogenous' (KEPT + tagged; it is a real signal, may carry KEGG/HMDB IDs).
      - NON-biological / technical -> origin in {contaminant, industrial, additive, surfactant,
        plasticizer, reagent}; final_class='xenobiotic-excluded' (dropped from analysis).
      - host-produced -> origin='endogenous' (default for KEGG/HMDB/structure-only mappings).

    For an exclusion pass final_class='xenobiotic-excluded' (or 'unmapped') with accepted={};
    origin is auto-suggested from the contaminant lexicon when omitted."""
    s = Session(workdir)
    e = s.get(feature_id)
    if e is None:
        return {"error": f"unknown feature_id {feature_id}"}
    final_class = LEGACY_CLASSES.get(final_class, final_class)  # accept old 'exogenous-excluded'
    if final_class and final_class not in ALL_CLASSES:
        return {"error": f"final_class must be one of {sorted(ALL_CLASSES)}"}
    accepted = accepted or {}
    warnings = []

    # --- origin: auto-suggest + coherence with the class ---
    if origin is None and final_class == "xenobiotic-excluded":
        origin = _lex.xenobiotic_category(e["original_name"]) or "contaminant"
    if origin is not None and origin not in VALID_ORIGINS:
        warnings.append(f"origin '{origin}' is not a recognized tag {sorted(VALID_ORIGINS)}.")
    implied = CLASS_FOR_ORIGIN.get(origin)
    if implied and final_class and implied != final_class:
        warnings.append(f"origin '{origin}' implies final_class '{implied}', not "
                        f"'{final_class}' — exogenous(diet/drug/microbial/plant) is KEPT; "
                        "xenobiotic(contaminant/industrial/additive/...) is EXCLUDED.")

    # --- terminal classes that carry no accepted IDs (excluded / unmapped) ---
    if final_class in (EXCLUDED_CLASSES | {"unmapped"}) and not accepted:
        e["final_class"] = final_class
        e["confidence"] = confidence or ("X" if final_class in EXCLUDED_CLASSES else "U")
        e["origin"] = origin
        s.add_decision(feature_id, "exclude", rationale, final_class=final_class, origin=origin)
        s.save()
        out = {"feature_id": feature_id, "final_class": final_class, "origin": origin,
               "counts": s.class_counts()}
        if warnings:
            out["warnings"] = warnings
        return out

    # --- exogenous kept WITHOUT a resolved DB id (real signal, structure not DB-mapped) ---
    if final_class == "exogenous" and not accepted:
        if origin is None:
            warnings.append("exogenous entry has no origin — tag diet/drug/microbial/plant.")
        e["final_class"] = "exogenous"
        e["confidence"] = confidence or "M4"
        e["origin"] = origin
        s.add_decision(feature_id, "accept", rationale, final_class="exogenous", origin=origin)
        s.save()
        out = {"feature_id": feature_id, "final_class": "exogenous", "origin": origin,
               "counts": s.class_counts()}
        if warnings:
            out["warnings"] = warnings
        return out

    # --- accept WITH ids (endogenous primary/structure-only, or exogenous carrying an id) ---
    try:
        fc = _commit(s, feature_id, accepted, confidence or "M2", rationale, final_class, origin=origin)
    except ValueError as ex:
        return {"error": str(ex)}
    if fc in (PRIMARY_CLASSES | {"structure-only", "exogenous"}):
        tags = _lex.xenobiotic_hits(e["original_name"])
        if tags and fc != "xenobiotic-excluded":
            warnings.append(f"name matches NON-biological class {tags} but kept as '{fc}' — "
                            "reconsider final_class='xenobiotic-excluded' (contaminant), not a "
                            "real/exogenous metabolite.")
        if fc == "exogenous" and e.get("origin") is None:
            warnings.append("exogenous entry has no origin — tag diet/drug/microbial/plant.")
        if "has_digit_locant" in e["normalized"].get("flags", []) and (confidence or "M2") != "M1":
            warnings.append("isomer/locant-sensitive name accepted via a non-exact route — "
                            "verify the accepted structure has the SAME locant/anomer.")
        warnings += _accept_warnings(s, feature_id, e)
    s.save()
    out = {"feature_id": feature_id, "final_class": fc, "accepted": e["accepted"],
           "confidence": e["confidence"], "origin": e.get("origin"), "counts": s.class_counts()}
    if warnings:
        out["warnings"] = warnings
    return out


def screen_exogenous(workdir: str) -> dict:
    """Deterministic NON-biological (xenobiotic) screen: LC-MS additives / surfactants /
    plasticizers / industrial reagents. EMITS evidence only. This lexicon detects the
    *xenobiotic → xenobiotic-excluded* classes; it deliberately does NOT flag biologically
    exogenous compounds (diet/drug/microbial/plant → the KEPT 'exogenous' class) because those
    require reasoning-layer judgement. Surfaces xenobiotic-class names that are still kept in
    the analysable set so the driver can reclassify them as xenobiotic-excluded."""
    s = Session(workdir)
    kept = PRIMARY_CLASSES | {"structure-only", "exogenous"}
    hits, kept_analysable = [], []
    for fid, e in s.entries.items():
        tags = _lex.xenobiotic_hits(e["original_name"])
        if not tags:
            continue
        row = {"feature_id": fid, "name": e["original_name"], "tags": tags,
               "suggested_origin": _lex.xenobiotic_category(e["original_name"]),
               "final_class": e.get("final_class"), "origin": e.get("origin")}
        hits.append(row)
        if e.get("final_class") in kept:
            kept_analysable.append(row)
    return {"n_xenobiotic_class": len(hits), "xenobiotic_hits": hits,
            "kept_analysable_review": kept_analysable,
            "note": "kept_analysable_review entries match a NON-biological class but are still "
                    "in the analysable set — reclassify final_class='xenobiotic-excluded' with "
                    "the suggested_origin. Biologically exogenous compounds (drugs, dietary/gut- "
                    "microbial/plant metabolites) are NOT in this lexicon — judge those yourself "
                    "and record them as final_class='exogenous' with origin diet/drug/microbial/plant."}


def harness_audit(workdir: str) -> dict:
    """GOVERNANCE / process-completeness auditor (read-only, no identity judgement): reads the
    session ledger + workdir artifacts and checks that the reasoning layer actually HONORED
    metabo-idmapper's OWN contract — no fabricated ID, no mass-only (W) candidate used as
    primary, final_class↔confidence coherent, fuzzy (M2/M3) accepts carry formula/mass
    verification, locant/anomer-sensitive names re-checked, xenobiotic-exclusion applied
    consistently, origin↔class coherent (exogenous kept vs xenobiotic excluded), flagged
    auto-accepts reviewed, id-gap KEGG re-tried, every CALL has a
    rationale, no entry left pending, gem_crosswalk ran, and the stage-7 always-emit artifacts
    exist. Emits a per-check pass/warn/fail scorecard so 'rules defined but not followed' is
    caught. Changes nothing. Run LAST (after coverage_summary)."""
    from . import harness as _harness
    return _harness.audit(workdir)


# ----------------------------------------------------------------------- GEM + export
def _gem_write(e: dict, label: str, model_path: str, mam: list[str], relation: str,
               rationale: str, names: list[str], curated: bool) -> None:
    """Store one model's result on an entry, per model label, and mirror the active view."""
    e.setdefault("gem", {})[label] = {"mam": mam, "relation": relation, "gem_names": names,
                                      "rationale": rationale, "model_path": model_path,
                                      "curated": curated}
    e["gem_mam"] = mam
    e["gem_relation"] = relation
    e["gem_model"] = label
    e["gem_cause"] = None if relation in GEM_MAPPED_RELATIONS else relation
    if curated:
        e["gem_curated"] = True


def _gem_name_suggest(e: dict, model_path: str | None, max_hits: int) -> list[dict]:
    """Search the model by the entry's own name + deterministic shorthand variants."""
    qname = _entry_query_name(e)
    alts = _short.variants(e["original_name"]) + _short.variants(qname)
    res = _gem.search(query=qname, queries=sorted(set(alts)), model_path=model_path,
                      max_hits=max_hits)
    hits = _isomer.screen(qname, res["hits"], name_key="gem_name")
    for h in hits:
        h.pop("all_names", None)
    return hits


def gem_crosswalk(workdir: str, feature_ids: list[str] | None = None,
                  model_path: str | None = None, model_label: str | None = None,
                  suggest_by_name: bool = True, max_suggestions: int = 6,
                  overwrite_curated: bool = False) -> dict:
    """Map accepted KEGG/HMDB/ChEBI ids to a genome-scale model's species (flux input).

    The xref route alone under-reports badly: GEMs annotate lipid species sparsely, so a
    verified Recon3D run crosswalked only 8 of 27 compounds by xref while the model actually
    contained 21 of them under other names. So every xref MISS is followed by a search of the
    model's own names (`suggest_by_name`, on by default) and the candidates come back with an
    isomer verdict per hit — then you commit the pick with `gem_assign`.

    Results are stored PER MODEL LABEL, so crosswalking against a second model (Recon3D as
    well as Mouse-GEM) does not overwrite the first. A curated `gem_assign` result is never
    silently overwritten by a re-run (pass overwrite_curated=True to force it).

    Sets gem_mam + gem_relation ('exact', or 'class-proxy' when the model only carries the
    generic R-group pool species) and gem_cause='id-gap' on a miss.
    """
    s = Session(workdir)
    path = model_path or GEM_MODEL
    label = model_label or _gem.model_label(path)
    s.data.setdefault("models", {})[label] = str(path)
    # endogenous primary mappings + kept exogenous compounds (both may resolve to a species);
    # xenobiotic-excluded / unmapped / structure-only-without-id are skipped.
    ids = feature_ids or [fid for fid, e in s.entries.items()
                          if e.get("final_class") in (PRIMARY_CLASSES | {"exogenous"})]
    mapped = curated_kept = 0
    by_relation: dict[str, int] = {}
    id_gap_with_kegg: list[dict] = []
    suggestions: list[dict] = []
    stats = _gem.model_stats(path)
    for fid in ids:
        e = s.get(fid)
        if e is None:
            continue
        if e.get("gem_curated") and e.get("gem_model") == label and not overwrite_curated:
            curated_kept += 1
            rel = e.get("gem_relation") or "exact"
            by_relation[rel] = by_relation.get(rel, 0) + 1
            mapped += bool(e.get("gem_mam"))
            continue
        acc = e.get("accepted", {})
        res = _gem.crosswalk(kegg=acc.get("kegg"), hmdb=acc.get("hmdb"),
                             chebi=acc.get("chebi"), model_path=path)
        if res["mapped"]:
            det = [_gem.species_detail(b, model_path=path) for b in res["mam_bases"]]
            pool = any(d.get("pool_species") for d in det)
            relation = "class-proxy" if pool else "exact"
            _gem_write(e, label, str(path), res["mam_bases"], relation,
                       f"xref {res['via']} -> model species", [d.get("gem_name") or ""
                                                               for d in det], curated=False)
            mapped += 1
            by_relation[relation] = by_relation.get(relation, 0) + 1
        else:
            _gem_write(e, label, str(path), [], "id-gap",
                       "no KEGG/HMDB/ChEBI xref in this model", [], curated=False)
            by_relation["id-gap"] = by_relation.get("id-gap", 0) + 1
            # A KEGG-bearing id-gap is often an anomer/stereo-specific KEGG the model lacks
            # (e.g. C03251 vs generic C00092). Surface it so the driver retries generic KEGG.
            if acc.get("kegg"):
                id_gap_with_kegg.append({"feature_id": fid, "name": e["original_name"],
                                         "kegg": acc["kegg"]})
                if "id_gap_try_generic_kegg" not in e.setdefault("_flags", []):
                    e["_flags"].append("id_gap_try_generic_kegg")
            if suggest_by_name:
                hits = _gem_name_suggest(e, path, max_suggestions)
                e["gem_candidates"] = hits
                if hits:
                    suggestions.append({"feature_id": fid, "name": e["original_name"],
                                        "candidates": hits})
    s.save()
    return {"model": str(path), "model_label": label, "model_stats": stats,
            "considered": len(ids), "gem_mapped": mapped, "by_relation": by_relation,
            "curated_kept": curated_kept, "id_gap_with_kegg": id_gap_with_kegg,
            "name_suggestions": suggestions,
            "note": "xref misses come back as name_suggestions from the model's OWN names, "
                    "each with an isomer verdict: 'ok'/'class-level' candidates are usually "
                    "the right species (commit with gem_assign), 'conflict' ones are a "
                    "different isomer — if you use one anyway record it as relation="
                    "'isomer-surrogate', and if nothing fits record 'model-scope-absent'. "
                    "class-proxy = the model only has the generic R-group pool species."}


def gem_search(query: str | None = None, formula: str | None = None,
               mass: float | None = None, mass_tol_ppm: float = 25.0,
               kegg: str | None = None, hmdb: str | None = None, chebi: str | None = None,
               model_path: str | None = None, workdir: str | None = None,
               feature_id: str | None = None, expand_shorthand: bool = True,
               max_hits: int = 15) -> dict:
    """Search a genome-scale model's OWN metabolite records — name, id text, formula, mass,
    xref — instead of relying on its xref annotations.

    This is the route that finds the species an xref crosswalk misses, because GEMs store
    lipids under names no identifier reaches: GCDCA lives in Recon3D as `dgchol`
    "Chenodeoxyglycocholate", the a2,3-sialyl GSL as "Sialyl-3-paragloboside", nLc4Cer as
    "Lactoneotetraosylceramide", and its glycolipids as sugar compositions
    ("(Gal)1 (Glc)1 (GlcNAc)1 (Cer)1"). With `expand_shorthand` (default) the query is also
    tried as its deterministic shorthand variants (LysoPC(16:0) ->
    "1-palmitoylglycero-3-phosphocholine" etc.).

    Every hit carries: the model name, formula, compartments, species ids, existing xrefs,
    `pool_species` (an R-group placeholder = CLASS-level, cannot represent chain length), and
    an `isomer_verdict` against the query. `same_formula_groups` lists formula-identical hits —
    positional-isomer candidates (pcholar_hs vs pcholn204_hs for LPC 20:4) that no formula or
    mass check can separate, so they must be decided on name/biology.

    Evidence only. With workdir+feature_id the hits are attached to the entry so `gem_assign`
    can verify the species you commit came from a search.
    """
    s = Session(workdir) if workdir else None
    e = s.get(feature_id) if (s and feature_id) else None
    if query is None and e is not None:
        query = _entry_query_name(e)
    if query is None and not (formula or mass or kegg or hmdb or chebi):
        return {"error": "pass query=, formula=, mass= or an xref (or workdir+feature_id)"}
    alts: list[str] = []
    if expand_shorthand and query:
        alts = _short.variants(query)
        if e is not None:
            alts += _short.variants(e["original_name"])
    res = _gem.search(query=query, queries=sorted(set(alts)), formula=formula, mass=mass,
                      mass_tol_ppm=mass_tol_ppm, kegg=kegg, hmdb=hmdb, chebi=chebi,
                      model_path=model_path, max_hits=max_hits)
    hits = _isomer.screen(query, res["hits"], name_key="gem_name") if query else res["hits"]
    if s and e is not None:
        e["gem_candidates"] = [{k: v for k, v in h.items() if k != "all_names"} for h in hits]
        s.add_decision(feature_id, "gem_search",
                       f"model search '{query}' (+{len(alts)} shorthand variants) -> "
                       f"{len(hits)} candidate species",
                       model=res["model"], candidates=[h["mam_base"] for h in hits])
        s.save()
    return {"query": query, "tried_variants": alts, "model": res["model"],
            "n_hits": len(hits), "hits": hits,
            "same_formula_groups": res["same_formula_groups"],
            "note": res["note"] + " Commit your pick with gem_assign (relation exact / "
                    "class-proxy / isomer-surrogate), or record model-scope-absent."}


def gem_assign(workdir: str, feature_id: str, rationale: str,
               mam: "list[str] | str | None" = None, relation: str = "exact",
               model_path: str | None = None, model_label: str | None = None) -> dict:
    """Commit the GEM CALL for one entry: which model species it maps to, and HOW.

    This is the model-side counterpart of `record_decision` — the only tool that sets
    gem_relation — and it exists because a curated mapping used to be un-recordable: the
    verified Recon3D run found 13 species by hand that `gem_crosswalk` could not store, so
    they lived only in a report while the ledger still said 'id-gap'.

    relation (choose deliberately — this is what a flux result gets interpreted against):
      exact              the species IS this compound.
      class-proxy        the model only carries the generic class/pool (R-group) species:
                         usable, but chain length/linkage is NOT represented — say so wherever
                         the ratio or flux is reported.
      isomer-surrogate   a DIFFERENT, related species stands in (N1-acetylspermine ->
                         N1-acetylspermidine). NOT an identity: it changes what the flux means,
                         so it needs an explicit reason and it is reported as a substitution.
      model-scope-absent the compound is genuinely not in the model — a real answer, not a
                         failure. Pass no mam.
      id-gap             should be in the model but nothing has resolved it yet. Pass no mam.

    Anti-fabrication: every species id must exist in the model (checked against the model
    itself), mirroring `record_decision`'s rule that no id may be invented. Results are stored
    per model label and are protected from being overwritten by a later `gem_crosswalk`.
    """
    s = Session(workdir)
    e = s.get(feature_id)
    if e is None:
        return {"error": f"unknown feature_id {feature_id}"}
    if relation not in GEM_RELATIONS:
        return {"error": f"relation must be one of {sorted(GEM_RELATIONS)}"}
    mams = [mam] if isinstance(mam, str) else list(mam or [])
    mams = [m for m in (str(x).strip() for x in mams) if m]
    if relation in GEM_MAPPED_RELATIONS and not mams:
        return {"error": f"relation '{relation}' needs mam=[...]; use 'model-scope-absent' or "
                         "'id-gap' to record that no species applies"}
    if mams and relation not in GEM_MAPPED_RELATIONS:
        return {"error": f"relation '{relation}' must not carry a mam"}
    if len((rationale or "").strip()) < 8:
        return {"error": "rationale must state WHY this species (>=8 chars) — it is the "
                         "provenance for every flux number derived from it"}

    path = model_path or (s.data.get("models", {}) or {}).get(model_label or "") or GEM_MODEL
    label = model_label or _gem.model_label(path)
    warnings: list[str] = []
    names, resolved, pooled = [], [], []
    for m in mams:
        det = _gem.species_detail(m, model_path=path)
        if not det.get("found"):
            return {"error": f"species '{m}' is not in model {label} — run gem_search and use "
                             "a mam_base it returned (no invented accessions)"}
        resolved.append(det["mam_base"])
        names.append(det.get("gem_name") or "")
        if det.get("pool_species"):
            pooled.append(det["mam_base"])
    known = {c.get("mam_base") for c in e.get("gem_candidates", [])}
    unbacked = [m for m in resolved if known and m not in known]
    if unbacked:
        warnings.append(f"{unbacked} was not among this entry's gem_search candidates — "
                        "confirm you meant this species.")
    if pooled and relation == "exact":
        # Not automatically wrong: a glycosphingolipid species is defined by its HEADGROUP
        # while the model still writes the ceramide tail as an R group. It IS wrong when the
        # axis the marker depends on is the part the R group swallows.
        warnings.append(f"{pooled} carries an R-group (pool) formula, so chain length and "
                        "linkage of that part are not represented. 'exact' is only right if "
                        "the axis this marker depends on (e.g. the glycan headgroup) IS "
                        "specified; if the chain/tail is the point, use 'class-proxy'.")
    qname = _entry_query_name(e)
    for m, nm in zip(resolved, names):
        if not nm:
            continue
        rep = _isomer.compare(qname, nm)
        if rep["verdict"] == "conflict" and relation != "isomer-surrogate":
            warnings.append(f"{m} ('{nm}') CONFLICTS with the name on "
                            f"{[c['axis'] for c in rep['conflicts']]} — either it is the wrong "
                            "species or the relation is 'isomer-surrogate'.")
        if rep["verdict"] == "class-level" and relation == "exact":
            warnings.append(f"{m} ('{nm}') is a class-level entry for a species-level name — "
                            "relation should be 'class-proxy'.")
    if relation in ("isomer-surrogate", "class-proxy") and len(rationale.strip()) < 20:
        warnings.append("a surrogate/proxy needs a fuller rationale: state what differs and "
                        "how it changes the interpretation of the flux.")
    for other, fids in s.gem_index(label).items():
        if other in resolved and relation == "exact":
            clash = [f for f in fids if f != feature_id
                     and s.entries[f].get("gem_relation") == "exact"]
            if clash:
                warnings.append(f"species {other} is already assigned 'exact' to {clash} — "
                                "two compounds cannot BE the same species (run collision_check).")

    s.data.setdefault("models", {})[label] = str(path)
    _gem_write(e, label, str(path), resolved, relation, rationale, names, curated=True)
    s.add_decision(feature_id, "gem_assign", rationale, model=label, mam=resolved,
                   relation=relation, gem_names=names)
    s.save()
    out = {"feature_id": feature_id, "model_label": label, "mam": resolved,
           "gem_names": names, "relation": relation,
           "caveat": {"isomer-surrogate": "SUBSTITUTION — report it wherever this flux is used",
                      "class-proxy": "CLASS-level input — chain length/linkage not represented",
                      "model-scope-absent": "not in the model — a real answer, not a failure",
                      "id-gap": "unresolved — keep searching or mark model-scope-absent",
                      "exact": None}[relation]}
    if warnings:
        out["warnings"] = warnings
    return out


def backfill_hmdb(workdir: str, db: str | None = None) -> dict:
    """Bridge missing HMDB for every non-excluded entry lacking an accepted HMDB, using its
    available ids (InChIKey/KEGG/ChEBI/PubChem) via BridgeDb in ONE batched call, then accept.

    Reflects the run-2 finding that HMDB rides along with KEGG matching and is under-counted;
    this maximizes HMDB coverage. Does not touch KEGG assignments or the final_class/origin of
    KEGG-mapped or exogenous entries (HMDB is added as an extra xref); a structure-only entry
    that gains HMDB → HMDB-mapped. xenobiotic-excluded entries are skipped.
    """
    s = Session(workdir)
    targets = [(fid, e) for fid, e in s.entries.items()
               if e.get("final_class") in ("KEGG-mapped", "structure-only", "exogenous")
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
    # HMDB accessions already accepted by ANOTHER entry. Bridging is per-entry and blind to
    # this: in a verified run it re-created the Lc3Cer/nLc4Cer shared-accession collision that
    # had just been hand-corrected. Never hand one accession to a second compound here.
    taken = {str(v): f for (k, v), fids in s.accepted_index().items() if k == "hmdb"
             for f in fids}
    gained, skipped_conflicts = [], []
    for fid, _ in targets:
        e = s.get(fid)
        hmdbs = [c["hmdb"] for c in e["candidates"] if c.get("hmdb")]
        free = [h for h in hmdbs if taken.get(str(h), fid) == fid]
        if hmdbs and not free:
            other = taken[str(hmdbs[0])]
            skipped_conflicts.append({
                "feature_id": fid, "name": e["original_name"], "hmdb": hmdbs[0],
                "already_accepted_for": other,
                "other_name": s.entries[other]["original_name"],
                "isomer_verdict": _isomer.compare(e["original_name"],
                                                  s.entries[other]["original_name"])["verdict"]})
            if "hmdb_backfill_conflict" not in e.setdefault("_flags", []):
                e["_flags"].append("hmdb_backfill_conflict")
            s.save()
            continue
        hmdbs = free
        if not hmdbs:
            continue
        acc = dict(e.get("accepted", {})); acc["hmdb"] = hmdbs[0]
        # preserve class + origin: an exogenous compound stays exogenous when it gains HMDB;
        # only a structure-only entry is allowed to be promoted (final_class=None → derive).
        keep_fc = e.get("final_class") if e.get("final_class") == "exogenous" else None
        record_decision(workdir, fid, rationale="HMDB backfill (BridgeDb)",
                        accepted=acc, final_class=keep_fc, origin=e.get("origin"),
                        confidence=e.get("confidence") or "M2")
        gained.append({"feature_id": fid, "name": e["original_name"], "hmdb": hmdbs[0]})
        taken[str(hmdbs[0])] = fid
        s = Session(workdir)
    return {"targets_missing_hmdb": len(targets), "hmdb_gained": len(gained),
            "gained": gained[:40], "skipped_conflicts": skipped_conflicts,
            "counts": Session(workdir).class_counts(),
            "note": "run before coverage_summary to normalize HMDB coverage; structure-only "
                    "with no id-bridgeable HMDB remain (need name→HMDB, a separate step). "
                    "skipped_conflicts = the only HMDB the bridge offered is already accepted "
                    "for another compound — it was NOT assigned; find the species-specific "
                    "accession (or leave it without HMDB) rather than sharing one."}


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
    reason) and the two origin buckets. Writes kegg_recovered.tsv, hmdb_recovered.tsv,
    unmapped_harmonization.tsv, exogenous_kept.tsv (biological outside-host, KEPT + origin),
    and xenobiotic_excluded.tsv (non-biological contaminant, EXCLUDED + origin)."""
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

    def category_of(e, reason):
        # prefer the recorded origin tag; else fall back to the contaminant lexicon / reason text
        origin = e.get("origin")
        if origin and origin != "endogenous":
            return origin
        cat = _lex.xenobiotic_category(e["original_name"])
        if cat:
            return cat
        r = (reason or "").lower()
        if any(w in r for w in ("antibiotic", "drug", "pharmac", "sulfonamide", "appetite",
                                "antihistamine", "alzheimer")):
            return "drug"
        if any(w in r for w in ("diet", "food", "plant", "dietary")):
            return "diet"
        if any(w in r for w in ("microb", "gut", "bacter")):
            return "microbial"
        if any(w in r for w in ("surfactant", "detergent", "diethanolamide", "sulfonic",
                                "amidopropyl", "alkyl sulfate")):
            return "surfactant"
        if any(w in r for w in ("plasticizer", "phthalate", "phosphate ester", "flame")):
            return "plasticizer"
        if any(w in r for w in ("additive", "trifluoroacetic", "mobile phase", "ion-pair", "ion pairing")):
            return "additive"
        if any(w in r for w in ("antioxidant", "nitrophenol", "tert-butyl", "reagent",
                                "building block", "industrial", "thiocyanate", "synthetic")):
            return "industrial"
        return "unspecified"

    def reason_of(e):
        excl = [d.get("rationale", "") for d in e.get("decisions", [])
                if d.get("action") in ("exclude", "accept", "auto_accept")]
        return excl[-1] if excl else (e.get("decisions", [{}])[-1].get("rationale", ""))

    kegg_rows, hmdb_rows, unmapped_rows, exo_kept_rows, xeno_rows = [], [], [], [], []
    for fid, e in s.entries.items():
        a = e.get("accepted", {})
        m = ma(e)
        ids = ";".join(f"{k}:{a[k]}" for k in ("kegg", "hmdb", "chebi", "pubchem") if a.get(k))
        if e.get("final_class") == "xenobiotic-excluded":
            reason = reason_of(e)
            xeno_rows.append([fid, e["original_name"], category_of(e, reason), ids, reason])
        elif e.get("final_class") == "exogenous":
            reason = reason_of(e)
            xref = "in-model" if e.get("gem_mam") else (e.get("gem_cause") or "")
            exo_kept_rows.append([fid, e["original_name"], category_of(e, reason), ids, xref, reason])
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

    # GEM curation table: the resolution axis a flux design needs — which species each marker
    # maps to, whether that is the species itself, a class/pool proxy or a surrogate, and which
    # ones the model genuinely does not contain (absence is a result, not a gap in the table).
    gem_rows = []
    for fid, e in s.entries.items():
        rel = e.get("gem_relation")
        if not rel:
            if e.get("gem_mam"):
                rel = "unspecified"      # mapped by an older run, before relations existed
            elif e.get("final_class") in (PRIMARY_CLASSES | {"exogenous"}):
                rel = "not-evaluated"    # model-relevant but never put to the model
            else:
                continue
        rec = (e.get("gem") or {}).get(e.get("gem_model") or "", {})
        det = next((d for d in e.get("decisions", []) if d.get("action") == "gem_assign"), {})
        gem_rows.append([fid, e["original_name"], e.get("gem_model") or "",
                         ";".join(e.get("gem_mam", [])), rel or "",
                         ";".join(rec.get("gem_names", []) or []),
                         "yes" if e.get("gem_curated") else "no",
                         rec.get("rationale") or det.get("rationale") or ""])

    tables = {
        "gem_curation.tsv": (["feature_id", "original_name", "model", "mam", "relation",
                              "gem_name", "curated", "rationale"], gem_rows),
        "kegg_recovered.tsv": (["feature_id", "original_name", "harmonized_name", "kegg",
                                "kegg_name", "route", "logic"], kegg_rows),
        "hmdb_recovered.tsv": (["feature_id", "original_name", "harmonized_name", "hmdb",
                                "route", "logic"], hmdb_rows),
        "unmapped_harmonization.tsv": (["feature_id", "original_name", "harmonized_names_tried",
                                        "final_class", "structure_ids", "reason"], unmapped_rows),
        "exogenous_kept.tsv": (["feature_id", "original_name", "origin", "ids", "gem_xref",
                                "reason"], exo_kept_rows),
        "xenobiotic_excluded.tsv": (["feature_id", "original_name", "origin", "ids", "reason"],
                                    xeno_rows),
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
    rel_counts: dict[str, int] = {}
    for r in gem_rows:
        rel_counts[r[4] or "?"] = rel_counts.get(r[4] or "?", 0) + 1
    return {"kegg_recovered": len(kegg_rows), "hmdb_recovered": len(hmdb_rows),
            "unmapped": len(unmapped_rows), "exogenous_kept": len(exo_kept_rows),
            "xenobiotic_excluded": len(xeno_rows), "gem_curation": len(gem_rows),
            "gem_by_relation": rel_counts, "exports": paths,
            "note": "recovered = mapped BEYOND MetaboAnalyst 1st pass; route/logic per row. "
                    "exogenous_kept = biological outside-host signal (kept); xenobiotic_excluded "
                    "= non-biological contaminant (dropped). gem_curation = per-entry model "
                    "resolution (exact / class-proxy / isomer-surrogate / model-scope-absent)."}


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
    add = ["final_class", "origin", "confidence", "kegg", "hmdb", "chebi", "pubchem",
           "inchikey", "gem_mam", "gem_relation"]
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
        df.at[i, "ID_origin"] = e.get("origin") or ""
        df.at[i, "ID_confidence"] = e.get("confidence", "")
        for k in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
            df.at[i, "ID_" + k] = acc.get(k, "")
        df.at[i, "ID_gem_mam"] = ";".join(e.get("gem_mam", []))
        df.at[i, "ID_gem_relation"] = e.get("gem_relation") or ""
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
    exogenous = counts.get("exogenous", 0)
    xenobiotic_excluded = counts.get("xenobiotic-excluded", 0)
    gem_by_relation: dict[str, int] = {}
    for e in s.entries.values():
        rel = e.get("gem_relation")
        if rel:
            gem_by_relation[rel] = gem_by_relation.get(rel, 0) + 1
    # species-level model input = the model carries THIS compound (not a pool/surrogate)
    gem_species_level = gem_by_relation.get("exact", 0)
    summary = {"total": n, "counts": counts,
               "primary_usable": primary, "has_kegg": has_kegg, "has_hmdb": has_hmdb,
               "gem_mapped": gem_mapped, "gem_by_relation": gem_by_relation,
               "gem_species_level": gem_species_level, "exogenous": exogenous,
               "xenobiotic_excluded": xenobiotic_excluded}
    if export:
        led = s.workdir / "master_ledger.tsv"
        cols = ["feature_id", "original_name", "normalized", "final_class", "origin",
                "confidence", "kegg", "hmdb", "chebi", "pubchem", "inchikey",
                "gem_model", "gem_mam", "gem_relation", "gem_cause", "flags",
                "n_candidates", "n_decisions"]
        with open(led, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(cols)
            for fid, e in s.entries.items():
                acc = e.get("accepted", {})
                w.writerow([fid, e["original_name"],
                            e["normalized"].get("normalized", ""),
                            e.get("final_class", ""), e.get("origin") or "",
                            e.get("confidence", ""),
                            acc.get("kegg", ""), acc.get("hmdb", ""), acc.get("chebi", ""),
                            acc.get("pubchem", ""), acc.get("inchikey", ""),
                            e.get("gem_model") or "",
                            ";".join(e.get("gem_mam", [])), e.get("gem_relation") or "",
                            e.get("gem_cause") or "", ";".join(e.get("_flags", []) or []),
                            len(e.get("candidates", [])), len(e.get("decisions", []))])
        cov = s.workdir / "coverage_summary.tsv"
        with open(cov, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["metric", "value", "denom", "pct"])
            for c in sorted(ALL_CLASSES):
                w.writerow([f"class:{c}", counts[c], n,
                            round(100 * counts[c] / n, 1) if n else 0])
            for name, val in (("primary_usable", primary), ("has_kegg", has_kegg),
                              ("has_hmdb", has_hmdb), ("gem_mapped", gem_mapped),
                              ("gem_species_level", gem_species_level),
                              ("exogenous", exogenous),
                              ("xenobiotic_excluded", xenobiotic_excluded)):
                w.writerow([name, val, n, round(100 * val / n, 1) if n else 0])
            for rel, val in sorted(gem_by_relation.items()):
                w.writerow([f"gem_relation:{rel}", val, n,
                            round(100 * val / n, 1) if n else 0])
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
    id_name_check, isomer_guard, collision_check,
    screen_exogenous, record_decision, backfill_hmdb, gem_crosswalk, gem_search, gem_assign,
    mapping_provenance, annotate_source, plot_coverage, export_code,
    export_report_ppt, coverage_summary, harness_audit,
]
