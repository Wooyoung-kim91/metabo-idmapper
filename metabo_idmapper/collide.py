"""Shared-identifier classification, used by both `collision_check` and the harness.

A single identifier standing for two entries has two very different meanings, and conflating
them would make the audit either useless or unusable:

  collision        two DIFFERENT compounds hold one id (Lc3Cer and nLc4Cer both landing on
                   PubChem CID 131770449 / HMDB0062485, because PubChem lists the two names as
                   synonyms of one entry). If those two are a ratio's numerator and denominator
                   the ratio collapses to exactly 1. This is a hard failure.
  class-id-shared  the id is a CLASS-level entry and the DB has no species entry, so several
                   species legitimately share it (LysoPC 16:0 / 18:2 / 20:4 all resolving to
                   KEGG C04230, "1-Acyl-sn-glycero-3-phosphocholine"). Not an error — but it
                   must be reported AS class-level, because a KEGG-keyed downstream analysis
                   still cannot tell the three apart.
  duplicate-name   the same compound listed twice under the same name. Usually harmless.

Evidence for class-level comes from the run itself: an `id_name_check` back-check that recorded
`db_resolution == "class"`, a `class_level_id` flag, or a candidate whose DB name resolves as a
generic class name. So back-checking an id is what downgrades a suspected collision to a
documented class-level share — the audit follows the evidence rather than a hardcoded list.
"""

from __future__ import annotations

from .engine import isomer as _isomer


def _class_level_evidence(entry: dict, db: str, val: str) -> str | None:
    """Why this entry's `db` id is believed to be a class-level DB entry (or None)."""
    for d in entry.get("decisions", []):
        if (d.get("action") == "backcheck" and d.get("db") == db
                and str(d.get("id")) == str(val)):
            if d.get("db_resolution") == "class" or d.get("verdict") == "class-level":
                return f"back-check: {db}:{val} is a class-level DB entry"
            for nm in d.get("db_names") or []:
                if _isomer.features(nm)["resolution"] == "class":
                    return f"back-check name '{nm}' is a class-level entry"
    for c in entry.get("candidates", []):
        if str(c.get(db) or "") == str(val) and c.get("note"):
            if _isomer.features(c["note"])["resolution"] == "class":
                return f"candidate name '{c['note']}' is a class-level entry"
    if "class_level_id" in (entry.get("_flags") or []):
        return "entry flagged class_level_id"
    return None


def classify_shared(entries: dict[str, dict], db: str, val: str,
                    fids: list[str]) -> dict:
    """Classify one shared (db, id) across the entries holding it."""
    names = [entries[f]["original_name"] for f in fids]
    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            rep = _isomer.compare(names[i], names[j])
            pairs.append({"a": names[i], "b": names[j], "verdict": rep["verdict"],
                          "conflicts": rep["conflicts"]})
    different = (any(p["verdict"] == "conflict" for p in pairs)
                 or len({n.strip().lower() for n in names}) > 1)
    evidence = [e for e in (_class_level_evidence(entries[f], db, val) for f in fids) if e]
    if not different:
        severity = "duplicate-name"
    elif evidence:
        severity = "class-id-shared"
    else:
        severity = "collision"
    return {"db": db, "id": val, "feature_ids": fids, "names": names, "pairwise": pairs,
            "severity": severity, "class_level_evidence": evidence,
            "advice": {
                "collision": "one id on two different compounds — find the species-specific id "
                             "for each before finalizing (a ratio between them is destroyed)",
                "class-id-shared": "the DB has only this class entry: keep it, but report these "
                                   "as class-level (a KEGG-keyed analysis cannot separate them)",
                "duplicate-name": "same compound listed twice — usually fine",
            }[severity]}


def scan(entries: dict[str, dict], accepted_index: dict) -> list[dict]:
    """Classify every accepted id held by more than one entry."""
    return [classify_shared(entries, db, val, fids)
            for (db, val), fids in sorted(accepted_index.items()) if len(fids) > 1]
