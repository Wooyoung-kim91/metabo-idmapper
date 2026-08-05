"""Lipidomics/glycolipid SHORTHAND -> the systematic and trivial names databases use.

Why: the names a metabolomics vendor reports and the names a DB or a GEM stores are almost
never the same string, and they share no whole tokens, so no token-overlap search finds them.
Verified pairs from a real run:

    LysoPC(16:0)                          -> "1-Palmitoylglycerophosphocholine"   (pcholpalm_hs)
    LysoPC(20:4)                          -> "1-Arachidonoyl-Glycero-3-Phospho..."(pcholar_hs)
    nLc4Cer                               -> "Lactoneotetraosylceramide"          (gal14acglc...)
    Lc3Cer                                -> "Acetylglucosaminylgalactosyl..."    (acglcgal...)
    3'-sialyl-neolactotetraosylceramide   -> "Sialyl-3-paragloboside"             (acngal14...)
    Glycochenodeoxycholic acid            -> "Chenodeoxyglycocholate"             (dgchol)
    15-HETE                               -> "15-hydroxyeicosatetraenoic acid"

This module turns one reported name into the DETERMINISTIC set of alternative query strings
that reach those records: chain shorthand -> trivial/systematic acyl names, head-group
abbreviations -> spelled-out class names, glyco shorthand -> its systematic alias, conjugated
bile-acid word order swaps, eicosanoid abbreviation expansion.

It proposes QUERIES ONLY — it never asserts an identity. Every returned string still has to
match something and survive `isomer.compare`.
"""

from __future__ import annotations

import re

# chain shorthand -> (trivial acyl, systematic acyl, systematic alkyl, systematic alkenyl)
_ACYL = {
    "14:0": ("myristoyl", "tetradecanoyl", "tetradecyl", "tetradecenyl"),
    "16:0": ("palmitoyl", "hexadecanoyl", "hexadecyl", "1Z-hexadecenyl"),
    "16:1": ("palmitoleoyl", "hexadecenoyl", "hexadecenyl", "hexadecadienyl"),
    "18:0": ("stearoyl", "octadecanoyl", "octadecyl", "1Z-octadecenyl"),
    "18:1": ("oleoyl", "octadecenoyl", "octadecenyl", "octadecadienyl"),
    "18:2": ("linoleoyl", "octadecadienoyl", "octadecadienyl", ""),
    "18:3": ("linolenoyl", "octadecatrienoyl", "", ""),
    "20:3": ("dihomo-gamma-linolenoyl", "eicosatrienoyl", "", ""),
    "20:4": ("arachidonoyl", "eicosatetraenoyl", "", ""),
    "20:5": ("eicosapentaenoyl", "eicosapentaenoyl", "", ""),
    "22:0": ("behenoyl", "docosanoyl", "", ""),
    "22:6": ("docosahexaenoyl", "docosahexaenoyl", "", ""),
    "24:0": ("lignoceroyl", "tetracosanoyl", "", ""),
    "24:1": ("nervonoyl", "tetracosenoyl", "", ""),
}

# head-group abbreviation -> (spelled-out class, systematic backbone)
_HEAD = {
    "lysopc": ("lysophosphatidylcholine", "sn-glycero-3-phosphocholine"),
    "lpc": ("lysophosphatidylcholine", "sn-glycero-3-phosphocholine"),
    "lysope": ("lysophosphatidylethanolamine", "sn-glycero-3-phosphoethanolamine"),
    "lpe": ("lysophosphatidylethanolamine", "sn-glycero-3-phosphoethanolamine"),
    "lysops": ("lysophosphatidylserine", "sn-glycero-3-phosphoserine"),
    "lysopi": ("lysophosphatidylinositol", "sn-glycero-3-phosphoinositol"),
    "pc": ("phosphatidylcholine", "sn-glycero-3-phosphocholine"),
    "pe": ("phosphatidylethanolamine", "sn-glycero-3-phosphoethanolamine"),
    "sm": ("sphingomyelin", "sphingomyelin"),
    "cer": ("ceramide", "ceramide"),
    "dhcer": ("dihydroceramide", "dihydroceramide"),
}

# Glyco-sphingolipid shorthand -> the systematic / trivial / cerebroside aliases KEGG and
# GEMs actually store, plus the SUGAR-COMPOSITION notation Recon3D uses for glycolipids
# ("(Gal)1 (Glc)1 (GlcNAc)1 (Cer)1"), which no name-similarity route can reach otherwise.
# NOTE the composition string is series-BLIND: Lc4Cer and nLc4Cer have identical compositions,
# so a composition hit must still be disambiguated (Recon3D encodes the linkage in the id:
# gal14acglc… = beta1-4 = neolacto vs galacglc… = beta1-3 = lacto).
_GLYCO = {
    "lc3cer": ["lactotriaosylceramide", "acetylglucosaminylgalactosylglucosylceramide",
               "GlcNAc-beta1-3-Gal-beta1-4-Glc-ceramide", "(Gal)1 (Glc)1 (GlcNAc)1 (Cer)1"],
    "lc4cer": ["lactotetraosylceramide",
               "galactosylacetylglucosaminylgalactosylglucosylceramide",
               "(Gal)2 (Glc)1 (GlcNAc)1 (Cer)1"],
    "nlc4cer": ["lactoneotetraosylceramide", "neolactotetraosylceramide", "paragloboside",
                "Gal-beta1-4-GlcNAc-beta1-3-Gal-beta1-4-Glc-ceramide",
                "(Gal)2 (Glc)1 (GlcNAc)1 (Cer)1"],
    "nlc5cer": ["lactoneopentaosylceramide"],
    "laccer": ["lactosylceramide", "galactosylglucosylceramide", "galactosyl glucosyl ceramide",
               "(Gal)1 (Glc)1 (Cer)1"],
    "lactosylceramide": ["galactosyl glucosyl ceramide", "galactosylglucosylceramide",
                         "(Gal)1 (Glc)1 (Cer)1"],
    "glccer": ["glucosylceramide", "glucocerebroside", "(Glc)1 (Cer)1"],
    "glucosylceramide": ["glucocerebroside", "glucosyl ceramide", "(Glc)1 (Cer)1"],
    "galcer": ["galactosylceramide", "galactocerebroside", "(Gal)1 (Cer)1"],
    "galactosylceramide": ["galactocerebroside", "galactosyl ceramide", "(Gal)1 (Cer)1"],
    "gb3": ["globotriaosylceramide", "(Gal)2 (Glc)1 (Cer)1"],
    "gb4": ["globotetraosylceramide"],
}

_CHAIN_RE = re.compile(r"(?<![\d.:])([dtmoOpP]?-?)(\d{1,2}:\d)(?![\d:])")
_HEAD_RE = re.compile(r"^\s*(lyso ?p[cesi]|lpc|lpe|lps|lpi|pc|pe|sm|dhcer|cer)\b", re.I)
_HETE_RE = re.compile(r"\b(\d{1,2})\s*[-(]?\s*(?:\(?[sr]\)?)?\s*-?\s*(h[ep]?ete|hpete|hode)\b",
                      re.I)
_BILE_RE = re.compile(r"^(glyco|tauro)(.*?)(cholic acid|cholate|cholic)$", re.I)


def _acyl_names(chain: str, ether: str | None) -> list[str]:
    """Acyl/alkyl/alkenyl name forms for one chain shorthand, honouring the ether prefix."""
    trivial, systematic, alkyl, alkenyl = _ACYL.get(chain, ("", "", "", ""))
    if ether == "plasmanyl":
        return [n for n in (alkyl,) if n]
    if ether == "plasmenyl":
        return [n for n in (alkenyl,) if n]
    return [n for n in (trivial, systematic) if n]


def variants(name: str, max_out: int = 14) -> list[str]:
    """Alternative query strings for one reported name (deterministic, de-duplicated).

    The input name itself is NOT included — these are the strings to try IN ADDITION to it.
    """
    raw = (name or "").strip()
    low = raw.lower()
    out: list[str] = []

    def add(s: str) -> None:
        s = re.sub(r"\s+", " ", (s or "").strip())
        if s and s.lower() != low and s not in out:
            out.append(s)

    # --- glyco shorthand (nLc4Cer, Lc3Cer, LacCer, GlcCer, Gb3) ---
    for key, aliases in _GLYCO.items():
        if re.search(rf"\b{key}\b", low) or low.replace(" ", "") == key:
            for a in aliases:
                add(a)

    # --- sialylated glyco series: linkage-preserving trivial names ---
    m = re.search(r"([36])'?\s*-?\s*sialyl", low)
    if m and re.search(r"neolacto|nlc|paraglobo", low):
        n = m.group(1)
        add(f"Sialyl-{n}-paragloboside")
        add(f"{n}'-LM1")
        add(f"alpha2,{'3' if n == '3' else '6'}-sialylneolactotetraosylceramide")
        add(f"{n}'-sialylparagloboside")

    # --- eicosanoid abbreviations (15-HETE -> 15-hydroxyeicosatetraenoic acid) ---
    for m in _HETE_RE.finditer(raw):
        loc, kind = m.group(1), m.group(2).lower()
        stem = {"hete": "hydroxyeicosatetraenoic acid",
                "hpete": "hydroperoxyeicosatetraenoic acid",
                "hode": "hydroxyoctadecadienoic acid"}.get(kind, "")
        if stem:
            add(f"{loc}-{stem}")
            add(f"{loc}-{stem.replace('eicosa', 'icosa')}")
            add(f"{loc}(S)-{kind.upper()}")

    # --- conjugated bile acids: DBs store the conjugate token next to 'chol' ---
    m = _BILE_RE.match(raw)
    if m:
        conj, stem, tail = m.group(1).lower(), m.group(2).strip("- "), m.group(3).lower()
        if stem:
            add(f"{stem}{conj}cholate")
            add(f"{stem}{conj}cholic acid")
            add(f"{stem}cholyl{'glycine' if conj == 'glyco' else 'taurine'}")
        add(f"{conj}{stem}{tail}")

    # --- glycerophospholipid / sphingolipid shorthand -> species names ---
    mh = _HEAD_RE.match(raw)
    chains = _CHAIN_RE.findall(raw)
    if mh:
        key = mh.group(1).lower().replace(" ", "")
        spelled, backbone = _HEAD.get(key, ("", ""))
        lyso = key.startswith(("lyso", "lp"))
        if key in ("cer", "dhcer", "sm"):
            # In "Cer(d18:1/24:0)" / "SM(d18:1/24:1)" the d18:1 is the SPHINGOID BASE, not an
            # acyl chain — naming it "1-oleoyl-…" would be wrong. Keep only the N-acyl chain.
            chains = [(p, c) for p, c in chains if p.strip("-").lower() not in ("d", "t")]
        glycero = "glycero" in backbone
        for pref, chain in chains:
            p = pref.strip("-").lower()
            ether = ("plasmenyl" if p == "p" else "plasmanyl" if p == "o" else None)
            for acyl in _acyl_names(chain, ether):
                if glycero and lyso:
                    add(f"1-{acyl}-{backbone}")
                    add(f"1-{acyl}glycero-3-phospho{backbone.split('phospho')[-1]}")
                elif glycero:
                    add(f"1-{acyl}-2-acyl-{backbone}")
                elif spelled:                       # sphingolipids: N-acyl naming
                    add(f"N-{acyl} {spelled}")
                if spelled:
                    add(f"{acyl} {spelled}")
        if spelled:
            for _pref, chain in chains:
                add(f"{spelled} {chain}")
            add(spelled)
        if any(p.strip('-').lower() == "p" for p, _c in chains):
            add("plasmenylcholine" if "c" in key else "plasmenylethanolamine")
            add("1-alkenyl-2-acyl-sn-glycero-3-phosphocholine")
            add("1-alkenyl-sn-glycero-3-phosphocholine")
        if any(p.strip('-').lower() == "o" for p, _c in chains):
            add("plasmanylcholine")
            add("1-alkyl-sn-glycero-3-phosphocholine")

    # --- oxidized lipids ---
    if re.search(r"\box|oxidi[sz]ed", low) and re.search(r"lpc|lysopc|pc\b", low):
        add("oxidized lysophosphatidylcholine")
        add("oxidized phosphatidylcholine")

    # --- polyamine acetylation ---
    if "spermine" in low or "spermidine" in low:
        if re.search(r"n1,\s*n12-?diacetyl|diacetyl", low):
            add("diacetylspermine")
            add("N1,N12-diacetylspermine")
        elif re.search(r"n1-?acetyl|monoacetyl", low):
            add("monoacetylspermine" if "spermine" in low else "monoacetylspermidine")
            add("N1-acetylspermine" if "spermine" in low else "N1-acetylspermidine")

    # --- generic string forms: drop the parenthetical, swap Lyso<->L prefix ---
    stripped = re.sub(r"\s*\([^)]*\)", "", raw).strip()
    add(stripped)
    if low.startswith("lyso"):
        add("L" + raw[4:].lstrip())
    elif re.match(r"^l[pP][ceEsSiI]", raw):
        add("Lyso" + raw[1:])

    return out[:max_out]
