"""Governance / process-completeness harness for a metabo-idmapper run (read-only).

This is NOT a mapping tool: it makes no identity call and touches nothing. It reads the
finished session ledger + workdir artifacts and checks that the reasoning layer (the driver)
actually HONORED metabo-idmapper's OWN operating contract — the rules written in
`guidance.CANONICAL` and `roles.DRIVER`/`roles.REVIEWER`. Rules that are "defined but not
followed" (an ID no tool produced, a mass-only weak hit used as primary, a locant-sensitive
name accepted via a non-exact route and never re-checked, a contaminant class excluded
inconsistently, a run that never reached the always-emit stage-7 artifacts) are exactly what
this catches. Emits a per-check pass/warn/fail scorecard; run it LAST.

Each check maps to a specific metabo-idmapper invariant — see the docstring on each `_check_*`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from .engine import lexicon as _lex
from .state import (CLASS_FOR_ORIGIN, EXCLUDED_CLASSES, PRIMARY_CLASSES,
                    VALID_ORIGINS, Session)

# classes that KEEP the compound in the analysable set (vs xenobiotic-excluded / unmapped)
KEPT_CLASSES = PRIMARY_CLASSES | {"structure-only", "exogenous"}

# confidence-tier vocabulary from guidance.CANONICAL
VALID_CONF = {"M1", "M2", "M3", "M4", "W", "X", "U"}
FUZZY_CONF = {"M2", "M3"}          # verified synonym / typo — evidence is REQUIRED first
_LEVELS = {"pass": 0, "warn": 1, "fail": 2}
_CAP = 25                           # cap offender lists so the summary stays small


def _decided(e: dict) -> bool:
    return bool(e.get("final_class"))


def _backed_ids(e: dict) -> set[tuple[str, str]]:
    """(db, id) pairs any tool actually produced for this entry (hmdb padded like _commit)."""
    seen: set[tuple[str, str]] = set()
    for c in e.get("candidates", []):
        for k in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
            if c.get(k):
                v = _lex.pad_hmdb(str(c[k])) if k == "hmdb" else str(c[k])
                seen.add((k, v))
    return seen


def _accept_decisions(e: dict) -> list[dict]:
    return [d for d in e.get("decisions", []) if d.get("action") in ("accept", "auto_accept")]


def _flags(e: dict) -> list[str]:
    return e.get("_flags", []) or []


# --------------------------------------------------------------------------- checks
# Each returns (level, summary, offenders). Level is worst-case over the entries it inspects.

def _check_anti_fabrication(s: Session) -> tuple[str, str, list]:
    """Contract: `record_decision` refuses any accepted ID no tool produced. Re-verify across
    the whole ledger (catches auto-accept bypasses / hand-edited ledgers)."""
    off = []
    for fid, e in s.entries.items():
        backed = _backed_ids(e)
        for k, v in (e.get("accepted") or {}).items():
            key = (k, _lex.pad_hmdb(str(v)) if k == "hmdb" else str(v))
            if key not in backed:
                off.append(f"{fid}:{k}={v}")
    return ("fail" if off else "pass",
            f"{len(off)} accepted id(s) not backed by any candidate" if off
            else "every accepted id is backed by a tool-produced candidate", off)


def _check_confidence_present(s: Session) -> tuple[str, str, list]:
    """Contract: every CALL carries a confidence tier (M1–U)."""
    off = [fid for fid, e in s.entries.items()
           if _decided(e) and e.get("confidence") not in VALID_CONF]
    return ("warn" if off else "pass",
            f"{len(off)} decided entr(y/ies) with missing/invalid confidence tier" if off
            else "all decided entries carry a valid confidence tier", off)


def _check_class_confidence_consistency(s: Session) -> tuple[str, str, list]:
    """Contract: final_class and confidence must agree. A primary-usable (KEGG/HMDB) mapping
    tagged W (mass-only), X (exogenous) or U (unresolved) is incoherent -> fail."""
    off = []
    for fid, e in s.entries.items():
        fc, conf = e.get("final_class"), e.get("confidence")
        if fc in PRIMARY_CLASSES and conf in {"W", "X", "U"}:
            off.append(f"{fid}:{fc}/{conf}")
    return ("fail" if off else "pass",
            f"{len(off)} primary-usable mapping(s) carry a weak/incoherent tier" if off
            else "final_class and confidence tiers are coherent", off)


def _check_w_tier_not_primary(s: Session) -> tuple[str, str, list]:
    """Contract (guidance): 'W mass/adduct-only weak candidate (do NOT use as primary).'"""
    off = [fid for fid, e in s.entries.items()
           if e.get("confidence") == "W" and e.get("final_class") in PRIMARY_CLASSES]
    return ("fail" if off else "pass",
            f"{len(off)} mass-only (W) candidate(s) used as primary-usable" if off
            else "no mass-only (W) candidate is used as primary", off)


def _check_verify_before_fuzzy(s: Session) -> tuple[str, str, list]:
    """Contract (driver): 'keggFind is noisy — NEVER accept a searched KEGG without
    verify_candidate passing.' A fuzzy-tier (M2/M3) accept needs recorded structure/formula
    evidence beyond the MetaboAnalyst pass."""
    off = []
    for fid, e in s.entries.items():
        if not _decided(e) or e.get("confidence") not in FUZZY_CONF:
            continue
        cands = e.get("candidates", [])
        has_evidence = any(c.get("source") != "metaboanalyst" and (c.get("formula") or c.get("mass"))
                           for c in cands)
        if not has_evidence:
            off.append(fid)
    return ("warn" if off else "pass",
            f"{len(off)} synonym/typo (M2/M3) accept(s) with no recorded formula/mass verification" if off
            else "fuzzy-tier accepts carry structure/formula verification", off)


def _check_isomer_locant_safety(s: Session) -> tuple[str, str, list]:
    """Contract (driver): a digit-locant/anomer-sensitive name accepted via a NON-exact route
    must be re-checked for the same locant. Surface any that slipped through unreviewed."""
    off = []
    for fid, e in s.entries.items():
        flags = e.get("normalized", {}).get("flags", [])
        if ("has_digit_locant" in flags
                and e.get("final_class") in KEPT_CLASSES
                and e.get("confidence") not in {"M1", None}):
            off.append(fid)
    return ("warn" if off else "pass",
            f"{len(off)} locant/anomer-sensitive name(s) accepted via a non-exact route — confirm same isomer" if off
            else "no unreviewed locant/anomer-sensitive acceptance", off)


def _check_contaminant_consistency(s: Session) -> tuple[str, str, list]:
    """Contract (driver): a NON-biological (xenobiotic-lexicon) name must be
    xenobiotic-excluded, applied CONSISTENTLY per class. Flags (a) xenobiotic-class names still
    kept in the analysable set (primary/structure-only/exogenous) and (b) a class where some
    members are excluded but others kept."""
    kept, by_class_excluded, by_class_kept = [], {}, {}
    for fid, e in s.entries.items():
        tags = _lex.xenobiotic_hits(e["original_name"])
        if not tags:
            continue
        is_kept = e.get("final_class") in KEPT_CLASSES
        if is_kept:
            kept.append(fid)
        for t in tags:
            (by_class_kept if is_kept else by_class_excluded).setdefault(t, []).append(fid)
    inconsistent = sorted(set(by_class_excluded) & set(by_class_kept))
    off = [f"{fid}:kept-analysable" for fid in kept]
    off += [f"class[{t}]:excluded+kept" for t in inconsistent]
    return ("warn" if off else "pass",
            f"{len(kept)} xenobiotic-class name(s) kept in the analysable set; "
            f"{len(inconsistent)} class(es) excluded inconsistently" if off
            else "xenobiotic-exclusion applied consistently", off)


def _check_origin_coherence(s: Session) -> tuple[str, str, list]:
    """Contract (new taxonomy): the origin tag must agree with final_class — biological
    exogenous origins (diet/drug/microbial/plant) go with 'exogenous' (KEPT); non-biological
    origins (contaminant/industrial/additive/...) go with 'xenobiotic-excluded' (EXCLUDED). An
    'exogenous' entry must carry an origin. Invalid or mismatched tags are governance failures."""
    off = []
    for fid, e in s.entries.items():
        fc, origin = e.get("final_class"), e.get("origin")
        if fc == "exogenous" and not origin:
            off.append(f"{fid}:exogenous-untagged")
            continue
        if origin is None:
            continue
        if origin not in VALID_ORIGINS:
            off.append(f"{fid}:invalid-origin[{origin}]")
            continue
        implied = CLASS_FOR_ORIGIN.get(origin)
        if implied and fc and implied != fc:
            off.append(f"{fid}:{origin}->{implied}!={fc}")
    return ("fail" if off else "pass",
            f"{len(off)} entr(y/ies) with a missing/invalid/mismatched origin↔class" if off
            else "origin tags cohere with final_class (exogenous kept, xenobiotic excluded)", off)


def _check_flagged_auto_accepts_reviewed(s: Session) -> tuple[str, str, list]:
    """Contract (driver+reviewer): exact_match 'auto_accept_review' flags (messy/trade-name
    M1 auto-accepts) are SUSPECT — route to the reviewer and re-verify. Treat an entry that
    still has only its single auto_accept decision as unreviewed."""
    off = []
    for fid, e in s.entries.items():
        if "auto_accept_review" not in _flags(e):
            continue
        decs = _accept_decisions(e)
        only_auto = decs and all(d.get("action") == "auto_accept" for d in decs) and len(decs) == 1
        if only_auto:
            off.append(fid)
    return ("warn" if off else "pass",
            f"{len(off)} messy/trade-name M1 auto-accept(s) never re-verified after flagging" if off
            else "flagged auto-accepts were reviewed/re-decided", off)


def _check_id_gap_generic_kegg(s: Session) -> tuple[str, str, list]:
    """Contract (driver): id_gap_try_generic_kegg = a KEGG that failed GEM crosswalk; search a
    generic KEGG and re-crosswalk. Flag entries still carrying a KEGG with no GEM MAM."""
    off = [fid for fid, e in s.entries.items()
           if "id_gap_try_generic_kegg" in _flags(e)
           and (e.get("accepted") or {}).get("kegg") and not e.get("gem_mam")]
    return ("warn" if off else "pass",
            f"{len(off)} anomer/stereo KEGG still failing GEM crosswalk — try a generic KEGG" if off
            else "no unresolved id-gap-with-KEGG entries", off)


def _check_rationale_quality(s: Session) -> tuple[str, str, list]:
    """Contract: the ledger is the provenance trail — every CALL needs a real one-line why."""
    trivial = {"", ".", "guess", "test", "n/a", "na", "-", "?"}
    off = []
    for fid, e in s.entries.items():
        for d in e.get("decisions", []):
            if d.get("action") not in ("accept", "auto_accept", "exclude"):
                continue
            r = (d.get("rationale") or "").strip().lower()
            if r in trivial or len(r) < 4:
                off.append(fid)
                break
    return ("warn" if off else "pass",
            f"{len(off)} decision(s) with an empty/trivial rationale" if off
            else "every decision carries a substantive rationale", off)


def _check_no_pending(s: Session) -> tuple[str, str, list]:
    """Contract: a finalized run has no entry left without a CALL."""
    off = s.pending_ids()
    return ("warn" if off else "pass",
            f"{len(off)} entr(y/ies) still pending — run not finalized" if off
            else "no pending entries", off[:_CAP])


def _check_gem_crosswalk_ran(s: Session) -> tuple[str, str, list]:
    """Contract: stage 6 maps accepted KEGG/HMDB/ChEBI -> Mouse-GEM MAM. If primary-usable
    entries exist but none carry gem_mam or a gem_cause, the crosswalk never ran."""
    primary = [fid for fid, e in s.entries.items() if e.get("final_class") in PRIMARY_CLASSES]
    if not primary:
        return ("pass", "no primary-usable entries yet — gem_crosswalk n/a", [])
    touched = any(e.get("gem_mam") or e.get("gem_cause") for e in s.entries.values())
    return ("pass" if touched else "warn",
            "gem_crosswalk has been run" if touched
            else f"{len(primary)} primary-usable entr(y/ies) but gem_crosswalk never ran", []
            if touched else primary[:_CAP])


def _check_finalize_artifacts(s: Session) -> tuple[str, str, list]:
    """Contract (stage 7 'ALWAYS emit'): master ledger + coverage + provenance tables + the two
    DB-matching figures + enriched_xref + raw-API reproduction code. Report which are missing."""
    wd = s.workdir
    required = [
        "master_ledger.tsv", "coverage_summary.tsv",
        "kegg_recovered.tsv", "hmdb_recovered.tsv",
        "unmapped_harmonization.tsv", "exogenous_kept.tsv", "xenobiotic_excluded.tsv",
        "enriched_xref.tsv",
        "figures/db_matching_upset.png", "figures/db_matching_improvement.png",
        "code/reproduce_mapping.py", "code/reproduce_mapping.ipynb",
    ]
    missing = [f for f in required if not (wd / f).exists()]
    return ("warn" if missing else "pass",
            f"{len(missing)} stage-7 always-emit artifact(s) missing (run coverage_summary)" if missing
            else "all stage-7 always-emit artifacts are present", missing)


_CHECKS: list[Callable[[Session], tuple[str, str, list]]] = [
    _check_anti_fabrication,
    _check_confidence_present,
    _check_class_confidence_consistency,
    _check_w_tier_not_primary,
    _check_verify_before_fuzzy,
    _check_isomer_locant_safety,
    _check_contaminant_consistency,
    _check_origin_coherence,
    _check_flagged_auto_accepts_reviewed,
    _check_id_gap_generic_kegg,
    _check_rationale_quality,
    _check_no_pending,
    _check_gem_crosswalk_ran,
    _check_finalize_artifacts,
]


def audit(workdir: str) -> dict[str, Any]:
    """Run the full governance scorecard over a session. Read-only."""
    s = Session(workdir)
    n = len(s.entries)
    if n == 0:
        return {"workdir": str(s.workdir), "n_entries": 0,
                "verdict": "warn", "score": {"pass": 0, "warn": 1, "fail": 0},
                "checks": [{"check": "ingest", "level": "warn",
                            "summary": "empty ledger — nothing to audit; run ingest_names first",
                            "offenders": []}],
                "scorecard": ["WARN  ingest — empty ledger, nothing to audit"],
                "note": "read-only governance audit; no ledger state was changed."}

    checks: list[dict] = []
    score = {"pass": 0, "warn": 0, "fail": 0}
    for fn in _CHECKS:
        level, summary, offenders = fn(s)
        score[level] += 1
        checks.append({"check": fn.__name__.removeprefix("_check_"), "level": level,
                       "summary": summary, "offenders": offenders[:_CAP]})

    verdict = max(checks, key=lambda c: _LEVELS[c["level"]])["level"]
    scorecard = [f"{c['level'].upper():5} {c['check']} — {c['summary']}" for c in checks]
    return {
        "workdir": str(s.workdir), "n_entries": n,
        "verdict": verdict, "score": score,
        "checks": checks, "scorecard": scorecard,
        "note": "read-only governance audit; no ledger state was changed. verdict = worst "
                "check level. fail = a hard contract violation (fix before finalizing); "
                "warn = review-and-confirm; pass = invariant honored.",
    }
