"""Session ledger — the authoritative, on-disk state for one mapping run.

The ledger is a single JSON file in the session workdir. Every tool reads it, appends
evidence or a decision, and writes it back. This is the provenance trail: original name
-> normalized -> candidates (with source) -> accepted IDs -> final_class/confidence ->
GEM MAM, plus a timestamped decisions log.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

LEDGER_NAME = "midmap_ledger.json"

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
    decisions: list[dict[str, Any]] = field(default_factory=list)


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


class Session:
    """Load/save the ledger and expose small helpers used by the tools."""

    def __init__(self, workdir: str):
        self.workdir = Path(workdir).expanduser().resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.path = self.workdir / LEDGER_NAME
        self.data: dict[str, Any] = {"created": _now(), "entries": {}}
        if self.path.exists():
            self.data = json.loads(self.path.read_text())

    # --- persistence ---
    def save(self) -> None:
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
