"""Review flags — raised with evidence, DERIVED open/closed, acknowledged with a reason.

A flag marks an entry the reasoning layer still has to answer for: an auto-accept nobody
re-checked, an id whose own DB record contradicts the name, an identifier shared with another
compound. The earlier design stored them as a bare list of strings that was appended to and
never cleared, so "is this still a problem?" had to be re-derived from decision history — the
audit ended up asking questions like "was there an accept AFTER the flag was raised?", which
is a proxy for the real condition and gets the answer wrong as soon as the flow changes.

Here a flag has two parts:

  raised      an append-only record (when, which tool, what evidence) — the provenance.
  still_open  a PREDICATE over the current entry + session — the truth.

`open_flags()` returns the flags whose predicate still holds, so resolving the underlying
problem closes the flag by itself and nothing goes stale. What cannot be derived is a human
judgement that a flag is acceptable as-is (a class-level id shared by three lyso-PCs is the
DB's limit, not an error) — that is `acknowledge()`, which records WHO decided and WHY and is
reported as an acknowledged flag, never as an absent one.

Severity separates work from context:
  action  something must change, or be acknowledged with a reason.
  info    a property of the mapping to carry into the report (e.g. this id is class-level).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable

from . import collide as _collide
from .engine import lexicon as _lex
from .state import PRIMARY_CLASSES

KEPT_CLASSES = PRIMARY_CLASSES | {"structure-only", "exogenous"}
_MAX_EVENTS = 5


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


# --------------------------------------------------------------------- predicates
def _re_decided(entry: dict) -> bool:
    """A deliberate call was made after the automatic one (auto_accept is its own action)."""
    return any(d.get("action") in ("accept", "exclude") for d in entry.get("decisions", []))


def _conflicting_id_still_accepted(entry: dict, _s: Any) -> bool:
    accepted = entry.get("accepted") or {}
    for d in entry.get("decisions", []):
        if d.get("action") == "backcheck" and d.get("verdict") == "conflict":
            if str(accepted.get(d.get("db")) or "") == str(d.get("id")):
                return True
    return False


def _class_level_id_accepted(entry: dict, _s: Any) -> bool:
    accepted = entry.get("accepted") or {}
    for d in entry.get("decisions", []):
        if d.get("action") == "backcheck" and (d.get("verdict") == "class-level"
                                               or d.get("db_resolution") == "class"):
            if str(accepted.get(d.get("db")) or "") == str(d.get("id")):
                return True
    return False


def _still_colliding(entry: dict, session: Any) -> bool:
    """Some accepted id of this entry is STILL held by a different compound."""
    if session is None:
        return True
    fid = entry.get("feature_id")
    index = session.accepted_index()
    for db, val in (entry.get("accepted") or {}).items():
        fids = index.get((db, str(val)), [])
        if len(fids) > 1:
            row = _collide.classify_shared(session.entries, db, str(val), fids)
            if row["severity"] == "collision" and fid in fids:
                return True
    return False


def _still_xenobiotic_cue(entry: dict, _s: Any) -> bool:
    return bool(_lex.xenobiotic_hits(entry.get("original_name", ""))
                and entry.get("final_class") in KEPT_CLASSES)


def _gem_gap_unresolved(entry: dict, _s: Any) -> bool:
    return bool((entry.get("accepted") or {}).get("kegg")
                and not entry.get("gem_mam")
                and entry.get("gem_relation") in (None, "id-gap"))


@dataclass(frozen=True)
class FlagSpec:
    name: str
    severity: str                       # "action" | "info"
    description: str
    resolution: str                     # what closes it
    still_open: Callable[[dict, Any], bool]


SPECS: dict[str, FlagSpec] = {f.name: f for f in (
    FlagSpec(
        "auto_accept_review", "action",
        "a messy / trade-name-like name was auto-accepted at M1 without review",
        "re-decide the entry with record_decision (or exclude it)",
        lambda e, s: not _re_decided(e)),
    FlagSpec(
        "isomer_token_conflict", "action",
        "the DB's own Match name disagrees with the query on a discriminant axis",
        "re-decide the entry after checking which isomer the id actually is",
        lambda e, s: not _re_decided(e)),
    FlagSpec(
        "id_name_conflict", "action",
        "a back-check found the accepted id resolves to a different compound/isomer",
        "accept a different id, or drop that id",
        _conflicting_id_still_accepted),
    FlagSpec(
        "class_level_id", "info",
        "an accepted id is a CLASS-level DB entry answering a species-level name",
        "keep it, but report the entry as class-level (or find a species-level id)",
        _class_level_id_accepted),
    FlagSpec(
        "shared_id_collision", "action",
        "an accepted id is also accepted for a different compound",
        "find the species-specific id for each entry",
        _still_colliding),
    FlagSpec(
        "hmdb_backfill_conflict", "action",
        "the only bridged HMDB belongs to another compound, so none was assigned",
        "find this compound's own accession, or accept that it has no HMDB",
        lambda e, s: not (e.get("accepted") or {}).get("hmdb")),
    FlagSpec(
        "id_gap_try_generic_kegg", "action",
        "a KEGG id failed the model crosswalk — it may be an anomer/stereo-specific id",
        "search the generic KEGG, or record the model relation (incl. model-scope-absent)",
        _gem_gap_unresolved),
    FlagSpec(
        "possible_xenobiotic", "info",
        "the name matches a NON-biological (contaminant) class but is kept as analysable",
        "reclassify as xenobiotic-excluded, or confirm it is a real metabolite",
        _still_xenobiotic_cue),
)}

ACTION_FLAGS = {n for n, f in SPECS.items() if f.severity == "action"}


# --------------------------------------------------------------------- api
def raise_flag(entry: dict, name: str, evidence: str, tool: str) -> None:
    """Record that `tool` observed the condition behind `name`, with its evidence."""
    rec = entry.setdefault("flags", {}).setdefault(name, {"raised": [], "acknowledged": None})
    rec["raised"] = (rec.get("raised") or [])[-(_MAX_EVENTS - 1):] + [
        {"ts": _now(), "tool": tool, "evidence": evidence}]


def acknowledge(entry: dict, name: str, why: str, who: str = "reasoning-layer") -> dict:
    """Mute a flag deliberately. The reason is kept and reported — this is a decision."""
    if name not in (entry.get("flags") or {}):
        return {"error": f"flag '{name}' was never raised for this entry"}
    if len((why or "").strip()) < 8:
        return {"error": "acknowledging a flag needs a reason (>=8 chars): it is a judgement "
                         "that the flagged state is acceptable, and it is reported as one"}
    entry["flags"][name]["acknowledged"] = {"ts": _now(), "why": why.strip(), "who": who}
    return {"flag": name, "acknowledged": True, "why": why.strip()}


def open_flags(entry: dict, session: Any = None, include_info: bool = True) -> list[str]:
    """Flags that were raised, are not acknowledged, and whose condition STILL holds."""
    out = []
    for name, rec in (entry.get("flags") or {}).items():
        if not (rec.get("raised") or []) or rec.get("acknowledged"):
            continue
        spec = SPECS.get(name)
        if spec is None:
            out.append(name)            # unknown flag: cannot prove it is resolved
            continue
        if not include_info and spec.severity == "info":
            continue
        if spec.still_open(entry, session):
            out.append(name)
    return sorted(out)


def acknowledged_flags(entry: dict) -> dict[str, dict]:
    return {n: r["acknowledged"] for n, r in (entry.get("flags") or {}).items()
            if r.get("acknowledged")}


def entries_with(session: Any, name: str) -> list[str]:
    """feature_ids whose `name` flag is currently open."""
    return [fid for fid, e in session.entries.items() if name in open_flags(e, session)]


def summary(session: Any) -> dict:
    """Counts of open / acknowledged / self-resolved flags across the session."""
    open_counts: dict[str, int] = {}
    ack_counts: dict[str, int] = {}
    resolved: dict[str, int] = {}
    for e in session.entries.values():
        opened = set(open_flags(e, session))
        for name, rec in (e.get("flags") or {}).items():
            if rec.get("acknowledged"):
                ack_counts[name] = ack_counts.get(name, 0) + 1
            elif name in opened:
                open_counts[name] = open_counts.get(name, 0) + 1
            elif rec.get("raised"):
                resolved[name] = resolved.get(name, 0) + 1
    return {"open": dict(sorted(open_counts.items())),
            "acknowledged": dict(sorted(ack_counts.items())),
            "self_resolved": dict(sorted(resolved.items())),
            "definitions": {n: {"severity": f.severity, "means": f.description,
                                "resolved_by": f.resolution}
                            for n, f in SPECS.items()
                            if n in open_counts or n in ack_counts}}
