"""Deterministic plumbing + evidence helpers surfaced by the reviewer's hardening notes.

- pad_hmdb: canonicalize HMDB accessions to the zero-padded 7-digit form (HMDB00349 ->
  HMDB0000349) so downstream joins are stable.
- xenobiotic_hits / xenobiotic_category: flag names matching known NON-biological classes
  (LC-MS additives, surfactants/detergents, plasticizers, industrial reagents) and map them to
  an origin category (additive/surfactant/plasticizer/industrial). This EMITS evidence only —
  it detects the *non-biological* (xenobiotic → excluded) classes deterministically. The
  BIOLOGICAL-exogenous call (diet/drug/microbial/plant → kept as final_class 'exogenous') stays
  with the reasoning layer, because those are not lexicon-matchable without judgement.
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
# dietary/plant compounds (those stay a reasoning-layer call). Each tag carries an `origin`
# category (state.XENOBIOTIC_ORIGINS) so a hit routes straight to final_class
# 'xenobiotic-excluded'.
_CONTAMINANT_PATTERNS = [
    (r"trifluoroacetic", "lcms-additive", "additive"),
    (r"\bphthalate\b", "plasticizer", "plasticizer"),
    (r"benzenesulfonic", "surfactant", "surfactant"),
    (r"\blauryl sulfate\b|dodecyl sulfate|tetradecylsulfate|\balkyl sulfate", "surfactant", "surfactant"),
    (r"diethanolamide|amidopropyl betaine|lauramidopropyl", "surfactant", "surfactant"),
    (r"butoxyethyl\).*phosphate|tributyl phosphate|tricresyl", "plasticizer/flame-retardant", "plasticizer"),
    (r"nitrophenol|di-?tert-butyl", "industrial-antioxidant", "industrial"),
    (r"thiocyanate|erucamide|pentadecylbenzoic", "industrial", "industrial"),
]


def xenobiotic_hits(name: str) -> list[str]:
    """Return the contaminant-class tags a name matches (empty = no lexicon hit)."""
    low = (name or "").lower()
    tags = [tag for pat, tag, _origin in _CONTAMINANT_PATTERNS if re.search(pat, low)]
    return sorted(set(tags))


def xenobiotic_category(name: str) -> str | None:
    """Map a name to its dominant NON-biological origin category — one of
    state.XENOBIOTIC_ORIGINS (additive/surfactant/plasticizer/industrial) — or None if no
    lexicon class matches. Used to auto-suggest origin for a 'xenobiotic-excluded' call.
    A None result does NOT mean endogenous: drugs/diet/microbial/plant are exogenous origins
    that this deterministic lexicon deliberately leaves to the reasoning layer."""
    low = (name or "").lower()
    for pat, _tag, origin in _CONTAMINANT_PATTERNS:
        if re.search(pat, low):
            return origin
    return None
