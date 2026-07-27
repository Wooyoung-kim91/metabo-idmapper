"""Deterministic metabolite-name normalization.

Pure string transforms only — NO identity guessing. Emits a normalized form plus
FLAGS the reasoning LLM should look at (parenthetical abbreviation, combined name,
possible typo cue). It never decides what the compound is.
"""

from __future__ import annotations

import re

_GREEK = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon",
    "ω": "omega", "μ": "u", "λ": "lambda",
}

# Common salt / hydrate suffixes that are not part of the compound identity.
_SALTS = [
    "hydrochloride", "dihydrochloride", "hydrobromide", "sodium salt", "potassium salt",
    "hemisulfate", "sulfate salt", "acetate salt", "monohydrate", "dihydrate",
    "trihydrate", "hydrate",
]


def _strip_greek(s: str) -> str:
    for g, r in _GREEK.items():
        s = s.replace(g, r)
    return s


def normalize(name: str) -> dict:
    """Return {normalized, casefold, flags[], parenthetical, alternatives[]}.

    - normalized:  whitespace-collapsed, greek-expanded, salt-stripped display form.
    - casefold:    lowercase key for exact matching.
    - parenthetical: text inside the first (...) — often an abbreviation expansion
                     or the spelled-out name; surfaced for the LLM, never auto-picked.
    - alternatives: split forms of a combined "A/B" or "A=B" name.
    - flags:       cues (abbrev_parenthetical, combined_name, has_digit_locant, ...).
    """
    raw = name.strip()
    s = _strip_greek(raw)
    s = re.sub(r"\s+", " ", s).strip()

    flags: list[str] = []
    parenthetical = None
    alternatives: list[str] = []

    # Parentheticals are ONLY an abbreviation expansion when the OUTSIDE token is a short
    # all-caps/alnum acronym and the inside is a longer name — e.g. "TMAO(trimethylamine
    # N-oxide)", "SAM(S-adenosyl-L-methionine)". A chemical SUBSTITUENT in parens
    # ("4-(Dimethylamino)...", "Tri(butoxyethyl) phosphate", "9-(2,3-dihydroxypropoxy)...")
    # must NOT be stripped — doing so mangles the identity. Default: keep the full string.
    base = s
    m_abbr = re.fullmatch(r"([A-Za-z0-9]{2,7})\s*\((.{3,})\)", s)
    if m_abbr and m_abbr.group(1).isupper() and not m_abbr.group(1).isdigit():
        acronym, expansion = m_abbr.group(1), m_abbr.group(2).strip()
        base = expansion                       # the real compound name
        parenthetical = expansion
        alternatives.append(acronym)
        flags.append("abbrev_parenthetical")
    else:
        # surface a trailing/standalone abbrev-like parenthetical for the LLM, but do NOT strip
        m = re.search(r"\(([^)]+)\)", s)
        if m:
            inner = m.group(1).strip()
            if not re.fullmatch(r"[RSrs+\-,/ ]+", inner) and len(inner) > 2:
                parenthetical = inner

    # Always keep the un-stripped, cleaned original as an alternative query, so a decision
    # is never forced on an empty/mangled candidate set.
    if s != base:
        alternatives.append(s)

    # combined names: "A/B", "A=B" (e.g. "...=PEP", "2-PG/3-PG")
    for sep in ("=", "/"):
        if sep in base:
            parts = [p.strip() for p in base.split(sep) if p.strip()]
            if len(parts) > 1:
                alternatives.extend(parts)
                flags.append("combined_name")
    # salt stripping (display form only)
    display = base
    low = display.lower()
    for salt in _SALTS:
        if low.endswith(" " + salt):
            display = display[: -(len(salt) + 1)].strip()
            flags.append("salt_stripped")
            break

    if re.search(r"\b\d+[- ]", base):
        flags.append("has_digit_locant")  # isomer-sensitive; do not auto-substitute

    # Conjugate spacing variants — DBs list acyl-carnitines/CoA closed-up. Verified miss:
    # "Butyryl carnitine" resolves only as "Butyrylcarnitine". Collapse the separator
    # before the conjugate token (NOT all whitespace, which would garble multi-word names).
    for suf in ("carnitine", "coa"):
        if re.search(rf"[ -]{suf}$", display, re.I):
            alt = re.sub(rf"[ -]{suf}$", suf, display, flags=re.I)
            if alt != display:
                alternatives.append(alt)
                flags.append("conjugate_spacing")
    if re.search(r"\bcarnitine\b", display, re.I):
        flags.append("acylcarnitine")  # cue: also try the HMDB/KEGG synonym form

    return {
        "normalized": display,
        "casefold": display.casefold(),
        "parenthetical": parenthetical,
        "alternatives": sorted(set(alternatives)),
        "flags": sorted(set(flags)),
    }
