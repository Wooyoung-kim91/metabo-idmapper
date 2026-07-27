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

PRIMARY_CLASSES = {"KEGG-mapped", "HMDB-mapped"}
ALL_CLASSES = PRIMARY_CLASSES | {"structure-only", "exogenous-excluded", "unmapped"}


@dataclass
class Entry:
    feature_id: str
    original_name: str
    normalized: dict[str, Any] = field(default_factory=dict)
    candidates: list[dict[str, Any]] = field(default_factory=list)
    accepted: dict[str, str] = field(default_factory=dict)  # kegg/hmdb/chebi/pubchem/inchikey
    final_class: str | None = None
    confidence: str | None = None
    gem_mam: list[str] = field(default_factory=list)
    gem_cause: str | None = None
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
