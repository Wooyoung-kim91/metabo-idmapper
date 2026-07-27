"""Deterministic plumbing + evidence helpers surfaced by the reviewer's hardening notes.

- pad_hmdb: canonicalize HMDB accessions to the zero-padded 7-digit form (HMDB00349 ->
  HMDB0000349) so downstream joins are stable.
- xenobiotic_hits: flag names matching known NON-biological classes (LC-MS additives,
  surfactants/detergents, plasticizers). This EMITS evidence only — the endogenous-vs-
  xenobiotic CALL stays with the reasoning layer (drugs/dietary compounds are a judgement).
"""

from __future__ import annotations

import re

_HMDB_RE = re.compile(r"^HMDB0*([0-9]+)$", re.I)


def pad_hmdb(acc: str | None) -> str | None:
    """HMDB00349 / HMDB0000349 -> HMDB0000349 (7-digit). Non-HMDB strings pass through."""
    if not acc:
        return acc
    m = _HMDB_RE.match(str(acc).strip())
    if not m:
        return acc
    return "HMDB" + m.group(1).zfill(7)


# Patterns for clearly non-biological classes the reviewer flagged as inconsistently kept.
# Deliberately conservative — surfactants/plasticizers/reagent additives, not ambiguous
# dietary/plant compounds (those stay a reasoning-layer call).
_CONTAMINANT_PATTERNS = [
    (r"trifluoroacetic", "lcms-additive"),
    (r"\bphthalate\b", "plasticizer"),
    (r"benzenesulfonic", "surfactant"),
    (r"\blauryl sulfate\b|dodecyl sulfate|tetradecylsulfate|\balkyl sulfate", "surfactant"),
    (r"diethanolamide|amidopropyl betaine|lauramidopropyl", "surfactant"),
    (r"butoxyethyl\).*phosphate|tributyl phosphate|tricresyl", "plasticizer/flame-retardant"),
    (r"nitrophenol|di-?tert-butyl", "industrial-antioxidant"),
    (r"thiocyanate|erucamide|pentadecylbenzoic", "industrial"),
]


def xenobiotic_hits(name: str) -> list[str]:
    """Return the contaminant-class tags a name matches (empty = no lexicon hit)."""
    low = (name or "").lower()
    tags = []
    for pat, tag in _CONTAMINANT_PATTERNS:
        if re.search(pat, low):
            tags.append(tag)
    return sorted(set(tags))
