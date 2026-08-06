"""Session ledger — the authoritative, on-disk state for one mapping run.

The ledger is a single JSON file in the session workdir. Every tool reads it, appends
evidence or a decision, and writes it back. This is the provenance trail: original name
-> normalized -> candidates (with source) -> accepted IDs -> final_class/confidence ->
GEM MAM, plus a timestamped decisions log.
"""

from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

LEDGER_NAME = "midmap_ledger.json"

# Ledger format version. Bump when the on-disk shape changes and add a step to migrate_ledger.
#   1  original: entries with accepted/final_class/origin/gem_mam, review cues in `_flags`
#   2  per-model `gem` results + gem_relation, flag records (flags.py), repaired model base ids
LEDGER_SCHEMA = 2

# identifier namespaces an entry can carry (order = report/column order)
ID_KEYS = ("kegg", "hmdb", "chebi", "pubchem", "inchikey")

# --- GEM relation axis: HOW an entry relates to the species recorded for a model ---
# exact             the model species IS this compound.
# class-proxy       the model carries only the generic class/pool (R-group) species: usable,
#                   but chain length / linkage is not represented — report it as class-level.
# isomer-surrogate  a DIFFERENT but closely related species stands in. Never an identity:
#                   the substitution must be stated wherever the flux result is used.
# model-scope-absent the compound is genuinely not in the model (not a mapping failure).
# id-gap            the compound should be in the model but no route resolved it yet.
GEM_RELATIONS = {"exact", "class-proxy", "isomer-surrogate", "model-scope-absent", "id-gap"}
GEM_MAPPED_RELATIONS = {"exact", "class-proxy", "isomer-surrogate"}
GEM_UNMAPPED_RELATIONS = {"model-scope-absent", "id-gap"}

# --- final_class vocabulary (ID-coverage / inclusion axis) ---
PRIMARY_CLASSES = {"KEGG-mapped", "HMDB-mapped"}      # endogenous, primary-usable (pathway/flux)
EXCLUDED_CLASSES = {"xenobiotic-excluded"}            # non-biological contaminant -> dropped
# 'exogenous' is KEPT (biologically real, from outside the host) and carries an origin tag.
ALL_CLASSES = PRIMARY_CLASSES | {"structure-only", "exogenous",
                                 "xenobiotic-excluded", "unmapped"}

# --- origin axis (fine provenance tag; orthogonal to ID coverage) ---
# Biological but from outside the host -> the compound is REAL signal; final_class 'exogenous'.
EXOGENOUS_ORIGINS = {"diet", "drug", "microbial", "plant"}
# Non-biological / technical -> not a metabolite; final_class 'xenobiotic-excluded'.
XENOBIOTIC_ORIGINS = {"contaminant", "industrial", "additive",
                      "surfactant", "plasticizer", "reagent"}
VALID_ORIGINS = {"endogenous"} | EXOGENOUS_ORIGINS | XENOBIOTIC_ORIGINS

# origin -> the final_class it implies (for validation / auto-routing)
CLASS_FOR_ORIGIN = ({o: "exogenous" for o in EXOGENOUS_ORIGINS}
                    | {o: "xenobiotic-excluded" for o in XENOBIOTIC_ORIGINS})

# legacy class name -> canonical (older ledgers used a single excluded bucket)
LEGACY_CLASSES = {"exogenous-excluded": "xenobiotic-excluded"}


@dataclass
class Entry:
    feature_id: str
    original_name: str
    normalized: dict[str, Any] = field(default_factory=dict)
    candidates: list[dict[str, Any]] = field(default_factory=list)
    accepted: dict[str, str] = field(default_factory=dict)  # kegg/hmdb/chebi/pubchem/inchikey
    final_class: str | None = None
    confidence: str | None = None
    origin: str | None = None  # endogenous | diet/drug/microbial/plant | contaminant/industrial/additive/...
    gem_mam: list[str] = field(default_factory=list)      # active model's species (base ids)
    gem_cause: str | None = None
    gem_relation: str | None = None   # exact | class-proxy | isomer-surrogate | absent | id-gap
    gem_model: str | None = None      # label of the model gem_mam/gem_relation refer to
    gem_curated: bool = False         # set by gem_assign; gem_crosswalk must not overwrite it
    gem: dict[str, Any] = field(default_factory=dict)     # per-model-label results
    gem_candidates: list[dict[str, Any]] = field(default_factory=list)  # gem_search evidence
    flags: dict[str, Any] = field(default_factory=dict)   # review flags (see flags.py)
    decisions: list[dict[str, Any]] = field(default_factory=list)


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


# --------------------------------------------------------------------------- migration
_TRAILING_SEP = re.compile(r"[_\-]+$")


def _repair_base_id(mam: str) -> str | None:
    """Repair a model accession left broken by the old compartment-stripping rule.

    That rule dropped ONE trailing lowercase letter, which turns the Recon3D-style
    'tdchola_c' into 'tdchola_' — a string no model resolves. Only a dangling separator is
    removed; anything else is left alone rather than guessed at.
    """
    fixed = _TRAILING_SEP.sub("", str(mam))
    return fixed if fixed and fixed != str(mam) else None


def migrate_ledger(data: dict[str, Any]) -> list[dict]:
    """Bring a loaded ledger up to LEDGER_SCHEMA **in memory**; the next save persists it.

    Every repair is appended to data['migrations'] so a changed value is visible in the
    provenance trail rather than appearing as if it had always been that way. Returns the
    steps applied in this load (empty when the ledger is already current).
    """
    version = int(data.get("schema_version") or 1)
    steps: list[dict] = []
    if version >= LEDGER_SCHEMA:
        data.setdefault("schema_version", LEDGER_SCHEMA)
        return steps

    if version < 2:
        repaired, moved = [], 0
        for fid, e in (data.get("entries") or {}).items():
            # v1 stored review cues as a bare list of names; v2 keeps a record per flag
            legacy = e.pop("_flags", None)
            if legacy:
                recs = e.setdefault("flags", {})
                for name in legacy:
                    rec = recs.setdefault(name, {"raised": [], "acknowledged": None})
                    rec["raised"].append({"ts": _now(), "tool": "migration",
                                          "evidence": "carried over from schema v1"})
                moved += 1
            e.setdefault("gem", {})
            e.setdefault("gem_candidates", [])
            e.setdefault("flags", {})
            e.setdefault("gem_relation", None)
            e.setdefault("gem_model", None)
            e.setdefault("gem_curated", False)
            # broken accessions written by the old base-id rule
            fixed = [(m, _repair_base_id(m)) for m in e.get("gem_mam") or []]
            if any(new for _old, new in fixed):
                e["gem_mam"] = [new or old for old, new in fixed]
                repaired += [{"feature_id": fid, "from": old, "to": new}
                             for old, new in fixed if new]
        steps.append({"schema": "1->2", "ts": _now(), "flag_records_created": moved,
                      "model_ids_repaired": repaired,
                      "note": "flag cues became records; model accessions with a dangling "
                              "separator (tdchola_ -> tdchola) were repaired — re-run "
                              "gem_crosswalk to confirm them against the model"})

    data["schema_version"] = LEDGER_SCHEMA
    if steps:
        data.setdefault("migrations", []).extend(steps)
    return steps


class LedgerConflict(RuntimeError):
    """The ledger changed on disk between load and save — refusing to overwrite.

    The ledger is the run's single source of truth, and a whole-file rewrite would silently
    discard whatever the other writer recorded. Reload and redo the operation.
    """


class Session:
    """Load/save the ledger and expose small helpers used by the tools.

    Write contract: a Session holds a whole-file snapshot, so `save()` REPLACES the file. It
    therefore refuses to write when the file changed underneath it (`LedgerConflict`, detected
    by the `revision` counter stored in the ledger) — two tool calls interleaving on one
    workdir would otherwise lose the earlier one's decisions entirely. Load once per tool call,
    mutate, save; do not re-open a Session mid-operation.
    """

    def __init__(self, workdir: str, *, migrate: bool = True):
        self.workdir = Path(workdir).expanduser().resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.path = self.workdir / LEDGER_NAME
        self.data: dict[str, Any] = {"created": _now(), "schema_version": LEDGER_SCHEMA,
                                     "revision": 0, "entries": {}}
        if self.path.exists():
            self.data = json.loads(self.path.read_text())
        self._rev = int(self.data.get("revision") or 0)
        self.migrations: list[dict] = migrate_ledger(self.data) if migrate else []

    # --- persistence ---
    def _disk_revision(self) -> int | None:
        """The revision counter currently on disk, or None when there is no ledger yet.

        A write counter INSIDE the file is what detects a concurrent write reliably: a
        timestamp+size stamp misses the case that matters most — two writes of equally sized
        content within one filesystem timestamp tick, which is exactly what two tool calls on
        the same ledger produce.
        """
        try:
            raw = self.path.read_text()
        except FileNotFoundError:
            return None
        try:
            return int(json.loads(raw).get("revision") or 0)
        except (ValueError, AttributeError):
            return 0

    def save(self, force: bool = False) -> None:
        """Write the ledger back. Raises LedgerConflict if another writer got there first."""
        on_disk = self._disk_revision()
        if not force and on_disk is not None and on_disk != self._rev:
            raise LedgerConflict(
                f"{self.path} changed since it was read (revision {on_disk}, this session "
                f"holds {self._rev}) — another tool call wrote to this workdir. Nothing was "
                "written; re-run this operation on the current ledger.")
        self._rev += 1
        self.data["schema_version"] = LEDGER_SCHEMA
        self.data["revision"] = self._rev
        tmp = self.path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.data, indent=2, ensure_ascii=False))
        tmp.replace(self.path)

    # --- entries ---
    @property
    def entries(self) -> dict[str, dict]:
        return self.data.setdefault("entries", {})

    def get(self, feature_id: str) -> dict | None:
        return self.entries.get(feature_id)

    def upsert(self, entry: Entry) -> None:
        self.entries[entry.feature_id] = asdict(entry)

    def add_candidate(self, feature_id: str, cand: dict) -> None:
        e = self.entries[feature_id]
        # de-dup on (source, kegg, hmdb, chebi, pubchem)
        key = (cand.get("source"), cand.get("kegg"), cand.get("hmdb"),
               cand.get("chebi"), cand.get("pubchem"))
        for c in e["candidates"]:
            if (c.get("source"), c.get("kegg"), c.get("hmdb"),
                    c.get("chebi"), c.get("pubchem")) == key:
                c.update({k: v for k, v in cand.items() if v})
                return
        e["candidates"].append(cand)

    def add_decision(self, feature_id: str, action: str, rationale: str, **payload) -> None:
        self.entries[feature_id]["decisions"].append(
            {"action": action, "rationale": rationale, "ts": _now(), **payload}
        )

    # --- summaries ---
    def class_counts(self) -> dict[str, int]:
        counts = {c: 0 for c in ALL_CLASSES}
        counts["pending"] = 0
        for e in self.entries.values():
            fc = e.get("final_class")
            if fc in counts:
                counts[fc] += 1
            else:
                counts["pending"] += 1
        return counts

    def pending_ids(self) -> list[str]:
        return [fid for fid, e in self.entries.items() if not e.get("final_class")]

    def accepted_index(self) -> dict[tuple[str, str], list[str]]:
        """{(db, id): [feature_id, ...]} over ACCEPTED ids — the basis for collision checks.

        Two different compounds holding one identifier is not a cosmetic problem: when the
        two are the numerator and denominator of a ratio, the ratio collapses to 1.
        """
        idx: dict[tuple[str, str], list[str]] = {}
        for fid, e in self.entries.items():
            for db, val in (e.get("accepted") or {}).items():
                if val:
                    idx.setdefault((db, str(val)), []).append(fid)
        return idx

    def gem_index(self, label: str | None = None) -> dict[str, list[str]]:
        """{model species base id: [feature_id, ...]} for the given model label (or active)."""
        idx: dict[str, list[str]] = {}
        for fid, e in self.entries.items():
            rec = (e.get("gem") or {}).get(label) if label else None
            mams = (rec or {}).get("mam") if rec else e.get("gem_mam") or []
            for m in mams or []:
                idx.setdefault(str(m), []).append(fid)
        return idx
