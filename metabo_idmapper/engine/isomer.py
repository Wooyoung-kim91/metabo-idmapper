"""Isomer / class-vs-species discriminant screen (deterministic, offline, no network).

WHY THIS EXISTS. In verified runs the dangerous failure was never a MISSING id — it was a
SILENTLY WRONG one that looked right: a bile acid carrying a cortisol-metabolite KEGG, a
plasmenyl (vinyl-ether) lyso-PC carrying the 1-ACYL class id, a neolacto (beta1-4) glycan
carrying the lacto (beta1-3) isomer, two different glycolipids sharing one PubChem/HMDB
entry. None of those are caught by comparing IDs, and none are caught by formula/mass —
the isomers are formula-identical. They ARE caught by comparing the DISCRIMINANT TOKENS of
the two names.

This module extracts those tokens (`features`) and compares a query name against a
candidate DB/model entry name (`compare`). It emits evidence only: a verdict plus the
axes that disagree. The identity CALL stays with the reasoning layer.

Axes, each of which was a real mis-assignment:
  skeleton        bile-acid vs corticosteroid vs eicosanoid vs polyamine ...  (TCDCA/C05472)
  ether_linkage   plasmenyl (P-, 1-alkenyl) vs plasmanyl (O-, 1-alkyl) vs acyl  (LPC P-16:0)
  glycan_series   neolacto (nLc, paragloboside) vs lacto (Lc)                (nLc4Cer/C04910)
  glycan_length   Lc3 (tri) vs nLc4 (tetra)                                  (Lc3/nLc4 share)
  sialyl_linkage  alpha2,3 vs alpha2,6                                       (a2,6 arm)
  chains          16:0 / 20:4 / d18:1 shorthand
  db_positions    Delta5,8,11,14 (arachidonoyl) vs Delta8,11,14,17 (n-3)     (LPC 20:4)
  omega           n-3 vs n-6 vs n-9
  oxidation       oxidized/hydroxy/hydroperoxy/epoxy present vs absent       (oxLPC, HETE)
  hydroxy_count   mono-HETE vs di-HETE                                       (15-HETE vs DiHETE)
  acetyl_count    N1-acetyl (1) vs N1,N12-diacetyl (2)
  polyamine       spermine vs spermidine vs putrescine backbone              (N1-acetylspermine)

`resolution` additionally reports whether a name is SPECIES-level (chain- or headgroup-
defined) or CLASS-level (generic acyl/alkyl/R-group pool) — the axis a flux model needs,
because a class-level id is usable but must be reported as such, not as the species.
"""

from __future__ import annotations

import re

# --------------------------------------------------------------------- token vocabularies

# Core skeletons. A candidate whose skeleton set is DISJOINT from the query's is a
# different compound class — the hardest possible conflict.
_SKELETONS = [
    ("bile-acid", r"cholic|cholate|cholanoic|cholan-|deoxychol|lithochol|ursodeoxychol|"
                  r"chenodeoxy|taurochol|glycochol|muricholic|\bbile acid"),
    ("corticosteroid", r"cortisol|cortisone|corticosterone|cortolone|\bcortol\b|pregnan|"
                       r"pregnenolone|androst|\bestra|\bestrone|estradiol|aldosterone"),
    # [a-z0-9]cer(?![a-z]) catches the glycolipid shorthand suffix (nLc4Cer, LacCer, dhCer)
    ("sphingolipid", r"ceramide|sphingo|sphingomyelin|globoside|globotri|paraglobo|"
                     r"ganglioside|\bgb[3-5]\b|cerebroside|[a-z0-9]cer(?![a-z])|\bcer\b"),
    ("glycerophospholipid", r"phosphocholine|phosphoethanolamine|phosphoserine|phosphoinositol|"
                            r"phosphatidyl|glycero-3-phospho|glycerophospho|lyso ?p[cesi]\b|"
                            r"lyso ?pc|lysope|\bgp[ce]\b|plasmenylcholine|plasmanylcholine|"
                            r"\bp[ce]\(|\blp[ce]\b"),
    ("eicosanoid", r"\bh[ep]?ete\b|hpete|\bhode\b|dihete|prostagland|thromboxane|leukotriene|"
                   r"epoxy(?:e)?icosa|(?:e)?icosa(?:tri|tetra|penta)eno(?:ic|ate)|"
                   r"octadecadieno(?:ic|ate)|lipoxin|resolvin"),
    ("polyamine", r"spermidine|spermine|putrescine|cadaverine"),
    ("acylcarnitine", r"carnitine"),
    ("nucleoside", r"adenosine|guanosine|cytidine|uridine|thymidine|inosine|\badenine\b|"
                   r"\bguanine\b|\bcytosine\b|\buracil\b|\bthymine\b"),
    ("amino-acid", r"\b(?:l-|d-|dl-)?(?:alanine|arginine|asparagine|aspartate|cysteine|"
                   r"glutamine|glutamate|glycine|histidine|isoleucine|leucine|lysine|"
                   r"methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|"
                   r"valine)\b"),
    ("acyl-CoA", r"\-coa\b|coenzyme a"),
    ("sugar", r"glucose|fructose|galactose|mannose|ribose|sucrose|maltose|trehalose"),
]

_ETHER = [
    ("plasmenyl", r"\bp-\d{1,2}:\d|plasmenyl|plasmalogen|1z-\w*enyl|alk-1-enyl|1-alkenyl|"
                  r"alkenyl-glycero|alkenylglycero"),
    ("plasmanyl", r"\bo-\d{1,2}:\d|plasmanyl|1-alkyl|alkyl-glycero|alkylglycero|"
                  r"alkyl-sn-glycero|o-alkyl"),
    ("acyl", r"\d-acyl|\bacyl-sn-glycero|acylglycero|acyl-glycero|diacyl|monoacyl|\bacyl\b"),
]

# glyco-sphingolipid series. neolacto patterns are tested FIRST because "lactoneo..." and
# "neolacto..." both contain the lacto stem.
_NEOLACTO = r"neolacto|lactoneo|\bn-?lc\d?|paraglobo"
_LACTO = r"\blc\d|lactotetraosyl|lactotriaosyl|\blacto(?!neo|se|syl|bio)"
_GLOBO = r"globo(?:tri|tetra|penta)|\bgb[3-5]\b|globoside"

# Glycosidic linkages as systematic names write them: "1,3-beta-D-Galactosyl…",
# "beta-D-Galactosyl-(1->4)-…", "beta1-4". KEGG states the SERIES only as a linkage — C04910's
# own name is "1,3-beta-D-Galactosyl-N-acetyl-D-glucosaminyl-…" with no "Lc4Cer" anywhere — so
# the linkage is the only place the lacto/neolacto distinction is written down.
_SUGAR = r"(?:galactosyl|glucosyl|glucosaminyl|galactosaminyl|mannosyl|fucosyl|n-acetyl|sialyl)"
_LINKAGE_RE = re.compile(
    rf"\((\d)\s*->\s*(\d)\)|\bbeta\s*(\d)\s*-\s*(\d)\b"
    rf"|(\d)\s*,\s*(\d)\s*-\s*(?:alpha|beta)?-?\s*[dl]?-?\s*{_SUGAR}")
# the Gal->GlcNAc linkage that DEFINES the series: beta1-4 = neolacto, beta1-3 = lacto
_SERIES_LINK_RE = re.compile(
    r"galactosyl\s*-?\s*\(1\s*->\s*(\d)\)\s*-?\s*n-acetyl"
    r"|1\s*,\s*(\d)\s*-\s*beta-d-galactosyl-n-acetyl"
    r"|galactosyl\s*-?\s*1\s*,\s*(\d)\s*-\s*n-acetyl")

_OSYL_LEN = {"mono": 1, "di": 2, "tri": 3, "tetra": 4, "penta": 5, "hexa": 6}
_OSYL_RE = re.compile(r"(mono|di|tri|tetra|penta|hexa)a?osyl")
_SERIES_LEN_RE = re.compile(r"\bn?-?lc(\d)|\bgb(\d)")

_SIALYL = [
    ("2,3", r"alpha\s*2\s*[,\-–]\s*3|\ba2\s*[,\-]\s*3|\b3'?\s*-?\s*sialyl|sialyl\s*-?\s*3|"
            r"neu5?ac?a?2\s*-?\s*3|\(2->3\)"),
    ("2,6", r"alpha\s*2\s*[,\-–]\s*6|\ba2\s*[,\-]\s*6|\b6'?\s*-?\s*sialyl|sialyl\s*-?\s*6|"
            r"neu5?ac?a?2\s*-?\s*6|\(2->6\)"),
]

_OXIDATION = (r"oxidi[sz]ed|\boxlpc|\box-?(?:lpc|pc|pe|pl|ldl)\b|hydroxy|hydroperoxy|epoxy|"
              r"\d-oxo|\boxo-|\bh[ep]?ete\b|hpete|\bhode\b|dihete|peroxide")

# Trivial acyl names -> (chain shorthand, omega, canonical double-bond positions).
_ACYL_NAMES = {
    "arachidonoyl": ("20:4", "n-6", {5, 8, 11, 14}),
    "arachidoyl": ("20:0", None, set()),
    "linoleoyl": ("18:2", "n-6", {9, 12}),
    "linolenoyl": ("18:3", "n-3", {9, 12, 15}),
    "oleoyl": ("18:1", "n-9", {9}),
    "elaidoyl": ("18:1", "n-9", {9}),
    "palmitoleoyl": ("16:1", "n-7", {9}),
    "palmitoyl": ("16:0", None, set()),
    "stearoyl": ("18:0", None, set()),
    "myristoyl": ("14:0", None, set()),
    "lauroyl": ("12:0", None, set()),
    "behenoyl": ("22:0", None, set()),
    "lignoceroyl": ("24:0", None, set()),
    "nervonoyl": ("24:1", "n-9", {15}),
    "docosahexaenoyl": ("22:6", "n-3", {4, 7, 10, 13, 16, 19}),
    "eicosapentaenoyl": ("20:5", "n-3", {5, 8, 11, 14, 17}),
}

_CHAIN_RE = re.compile(r"(?<![\d.:])(\d{1,2}):(\d)(?![\d:])")
_OMEGA_RE = re.compile(r"\b(?:n|omega)\s*-?\s*([369])\b")
_HETE_LOC_RE = re.compile(r"(\d{1,2})\s*\(?[sr]?\)?\s*-?\s*h[ep]?ete")
_DIHETE_RE = re.compile(r"\bdi-?h[ep]?ete|\bdihydroxy")
_DBPOS_RE = re.compile(r"(\d{1,2}(?:,\d{1,2}){1,5})\s*-?\s*(?:all-)?(?:cis-|trans-)?[a-z]*"
                       r"(?:enoyl|enoic|enyl)")
_GENERIC_RE = re.compile(r"\b\d?-?acyl\b|\balkyl\b|\bany\b|\br\d?-group|\bunspecified\b")
_HEADGROUP_RE = re.compile(r"gluco?syl|galactosyl|lactosyl|neolacto|lactoneo|\bn?-?lc\d|"
                           r"globo|paraglobo|sialyl|fucosyl|gangli|\bgb[3-5]\b|\bgm[1-3]\b")


def _norm(name: str) -> str:
    """Lowercase + strip stereo letters inside double-bond position lists (5Z,8Z -> 5,8)."""
    low = (name or "").lower()
    low = low.replace("α", "alpha").replace("β", "beta").replace("→", "->")
    return re.sub(r"(\d+)[ze](?=[,)\s\-])", r"\1", low)


# Systematic acyl stems -> carbon count / unsaturation, so a model or KEGG entry written
# systematically ("1-(8,11,14,17-eicosatetraenoyl)-glycero-3-phosphocholine") is comparable
# with lipidomics shorthand ("LPC 20:4").
_CARBONS = {"octan": 8, "decan": 10, "dodec": 12, "tetradec": 14, "pentadec": 15,
            "hexadec": 16, "heptadec": 17, "octadec": 18, "nonadec": 19, "eicos": 20,
            "icos": 20, "heneicos": 21, "docos": 22, "tricos": 23, "tetracos": 24,
            "hexacos": 26}
_UNSAT = {"anoyl": 0, "anoic": 0, "enoyl": 1, "enoic": 1, "dienoyl": 2, "dienoic": 2,
          "trienoyl": 3, "trienoic": 3, "tetraenoyl": 4, "tetraenoic": 4,
          "pentaenoyl": 5, "pentaenoic": 5, "hexaenoyl": 6, "hexaenoic": 6}
_SYSTEMATIC_RE = re.compile(
    "(" + "|".join(sorted(_CARBONS, key=len, reverse=True)) + ")"
    "(" + "|".join(sorted(_UNSAT, key=len, reverse=True)) + ")")


def _systematic_chains(low: str) -> set[str]:
    """Chain shorthand implied by systematic acyl names (eicosatetraenoyl -> 20:4)."""
    return {f"{_CARBONS[stem]}:{_UNSAT[suf]}" for stem, suf in _SYSTEMATIC_RE.findall(low)}


def _matched(low: str, table: list[tuple[str, str]]) -> set[str]:
    return {tag for tag, pat in table if re.search(pat, low)}


# --------------------------------------------------------------------- feature extraction
def features(name: str) -> dict:
    """Extract the discriminant features of ONE name. Pure string logic, no lookups."""
    low = _norm(name)
    skeleton = _matched(low, _SKELETONS)
    ether = _matched(low, _ETHER)
    chains = {f"{a}:{b}" for a, b in _CHAIN_RE.findall(low)}
    omega = None
    dbpos: set[int] = set()

    for trivial, (chain, om, pos) in _ACYL_NAMES.items():
        if trivial not in low:
            continue
        chains.add(chain)
        omega = omega or om
        dbpos |= pos
        if re.search(r"glycero|phosphocholine|phosphoethanolamine", low):
            ether.add("acyl")
    chains |= _systematic_chains(low)
    m = _OMEGA_RE.search(low)
    if m:
        omega = f"n-{m.group(1)}"
    for grp in _DBPOS_RE.findall(low):
        dbpos |= {int(x) for x in grp.split(",")}

    # a chain-bearing glycerophospholipid with no ether marker is an ESTER (1-acyl) species
    if chains and "glycerophospholipid" in skeleton and not ether:
        ether.add("acyl")

    linkages = ["-".join(g for g in m.groups() if g) for m in _LINKAGE_RE.finditer(low)]
    linkages = ["1-" + lk if len(lk) == 1 else lk for lk in linkages]

    series = None
    if re.search(_NEOLACTO, low):
        series = "neolacto"
    elif re.search(_LACTO, low):
        series = "lacto"
    elif re.search(_GLOBO, low):
        series = "globo"
    else:
        m = _SERIES_LINK_RE.search(low)
        if m:
            link = next(g for g in m.groups() if g)
            series = {"4": "neolacto", "3": "lacto"}.get(link)

    glycan_len = None
    m = _SERIES_LEN_RE.search(low)
    if m:
        glycan_len = int(m.group(1) or m.group(2))
    else:
        m = _OSYL_RE.search(low)
        if m:
            glycan_len = _OSYL_LEN[m.group(1)]
        elif "paraglobo" in low:
            glycan_len = 4  # paragloboside IS neolactotetraosylceramide (4 sugars)

    sialyl = None
    for tag, pat in _SIALYL:
        if re.search(pat, low):
            sialyl = tag
            break

    oxidation = bool(re.search(_OXIDATION, low))
    # hydroPEROXY vs hydroXY is a different compound, not a spelling variant: Recon3D carries
    # 15HPET (the hydroperoxide precursor) but no 15-HETE, so the two must never merge.
    peroxidation = (bool(re.search(r"hydroperoxy|hpete|\bhpet\b|\bhpode\b|peroxide", low))
                    if oxidation else None)
    hydroxy_count = None
    if re.search(r"h[ep]?ete|hode", low):
        hydroxy_count = 2 if _DIHETE_RE.search(low) else 1
    elif "eicosanoid" in skeleton and "hydroxy" in low:
        # systematic eicosanoid spelling ("(15S)-15-Hydroxy-…-eicosatetraenoate")
        hydroxy_count = 2 if "dihydroxy" in low else 1
    hete_locants = {int(x) for x in _HETE_LOC_RE.findall(low)}

    polyamine = None
    for tag in ("spermidine", "spermine", "putrescine", "cadaverine"):
        if tag in low:
            polyamine = tag
            break

    # Acetylation count is a discriminant only where it distinguishes the PARENT compound
    # (mono- vs di-acetylspermine). N-acetylhexosamines inside a glycan name are structural
    # parts of the headgroup, not an acetylation of the parent — counting them would fire a
    # false conflict on every glycosphingolipid.
    acetyl = None
    if polyamine or ("amino-acid" in skeleton):
        if re.search(r"di-?acetyl|\bbis-?acetyl", low):
            acetyl = 2
        elif re.search(r"monoacetyl|n\d*-acetyl|\bacetyl", low):
            acetyl = 1

    headgroup = bool(_HEADGROUP_RE.search(low))
    generic = bool(_GENERIC_RE.search(low))
    lipid = bool(skeleton & {"sphingolipid", "glycerophospholipid"})
    if chains or headgroup:
        resolution = "species"
    elif lipid or generic:
        resolution = "class"
    else:
        resolution = "species"

    return {"name": name, "skeleton": sorted(skeleton), "ether_linkage": sorted(ether),
            "chains": sorted(chains), "omega": omega, "db_positions": sorted(dbpos),
            "glycan_series": series, "glycan_length": glycan_len, "sialyl_linkage": sialyl,
            "linkages": linkages,
            "oxidation": oxidation, "peroxidation": peroxidation,
            "hydroxy_count": hydroxy_count,
            "hete_locants": sorted(hete_locants), "acetyl_count": acetyl,
            "polyamine_backbone": polyamine, "resolution": resolution,
            "generic_token": generic}


# --------------------------------------------------------------------- comparison
_SET_AXES = ("skeleton", "ether_linkage")
_SCALAR_AXES = ("glycan_series", "glycan_length", "sialyl_linkage", "omega",
                "peroxidation", "hydroxy_count", "acetyl_count", "polyamine_backbone")


def compare(query: str, candidate: str) -> dict:
    """Compare a query metabolite name against a candidate DB/model entry name.

    Returns {verdict, conflicts[], ambiguities[], shared[], query/candidate features}.

    verdict:
      conflict     — an axis disagrees: the candidate is a DIFFERENT compound/isomer.
                     Do NOT accept it as the same identity (it may still be a documented
                     surrogate, which must be recorded as such).
      class-level  — the query is species-level (chains/headgroup) but the candidate is a
                     generic class/pool entry. Usable, but report it as class-level.
      ambiguous    — the query specifies an axis the candidate leaves open (or vice versa),
                     so the two cannot be distinguished from names alone: resolve it with
                     another evidence source before accepting.
      ok           — no axis disagrees.
    """
    q, c = features(query), features(candidate)
    conflicts, ambiguities, shared = [], [], []

    for ax in _SET_AXES:
        qa, ca = set(q[ax]), set(c[ax])
        if qa and ca:
            if qa & ca:
                shared.append(ax)
            else:
                conflicts.append({"axis": ax, "query": sorted(qa), "candidate": sorted(ca)})
        elif qa or ca:
            ambiguities.append({"axis": ax, "query": sorted(qa), "candidate": sorted(ca)})

    for ax in _SCALAR_AXES:
        qv, cv = q[ax], c[ax]
        if qv is None and cv is None:
            continue
        if qv is None or cv is None:
            ambiguities.append({"axis": ax, "query": qv, "candidate": cv})
        elif qv == cv:
            shared.append(ax)
        else:
            conflicts.append({"axis": ax, "query": qv, "candidate": cv})

    qch, cch = set(q["chains"]), set(c["chains"])
    if qch and cch:
        if qch == cch:
            shared.append("chains")
        elif qch & cch:
            ambiguities.append({"axis": "chains", "query": sorted(qch), "candidate": sorted(cch),
                                "note": "overlapping but not identical chain sets"})
        else:
            conflicts.append({"axis": "chains", "query": sorted(qch), "candidate": sorted(cch)})

    qp, cp = set(q["db_positions"]), set(c["db_positions"])
    if qp and cp:
        if qp == cp:
            shared.append("db_positions")
        else:
            conflicts.append({"axis": "db_positions", "query": sorted(qp),
                              "candidate": sorted(cp)})
    elif (qp or cp) and (qch & cch):
        # same chain length/unsaturation, but only ONE side pins the double-bond positions:
        # a positional isomer (e.g. 20:4 n-6 arachidonoyl vs 20:4 n-3) cannot be told apart
        # from the names alone — this is exactly the LPC 20:4 / pcholar_hs vs pcholn204_hs case.
        ambiguities.append({"axis": "db_positions", "query": sorted(qp), "candidate": sorted(cp),
                            "note": "same chain shorthand, positions specified on one side "
                                    "only — positional isomers are indistinguishable by name"})

    if q["linkages"] and c["linkages"]:
        if q["linkages"] == c["linkages"]:
            shared.append("linkages")
        else:
            # Reported as an AMBIGUITY, never a conflict: two spellings of the same compound
            # state their linkages with different completeness ("beta-D-Galactosyl-(1->4)-…"
            # vs "…-1,4-D-glucosylceramide"), so a sequence mismatch is a prompt to read the
            # linkages, not proof of a different isomer. The series axis carries the verdict.
            ambiguities.append({"axis": "linkages", "query": q["linkages"],
                                "candidate": c["linkages"],
                                "note": "glycosidic linkages as written differ — check them "
                                        "(beta1-3 = lacto vs beta1-4 = neolacto)"})

    ql, cl = set(q["hete_locants"]), set(c["hete_locants"])
    if ql and cl and ql != cl:
        conflicts.append({"axis": "hete_locants", "query": sorted(ql), "candidate": sorted(cl)})

    if q["oxidation"] != c["oxidation"] and (q["skeleton"] and c["skeleton"]):
        conflicts.append({"axis": "oxidation", "query": q["oxidation"],
                          "candidate": c["oxidation"]})

    if conflicts:
        verdict = "conflict"
    elif q["resolution"] == "species" and c["resolution"] == "class":
        verdict = "class-level"
    elif ambiguities:
        verdict = "ambiguous"
    else:
        verdict = "ok"

    return {"verdict": verdict, "query": query, "candidate": candidate,
            "conflicts": conflicts, "ambiguities": ambiguities,
            "shared_axes": [s for s in shared if isinstance(s, str)],
            "candidate_resolution": c["resolution"], "query_resolution": q["resolution"],
            "query_features": q, "candidate_features": c}


def screen(query: str, candidates: list[dict], name_key: str = "name") -> list[dict]:
    """Annotate a candidate list with the isomer verdict against `query`.

    Each candidate dict gains `isomer_verdict` + `isomer_conflicts` (+ `resolution`).
    Candidates with no name are passed through untouched (nothing to compare).
    Sorted safest-first: ok < class-level < ambiguous < conflict.
    """
    rank = {"ok": 0, "class-level": 1, "ambiguous": 2, "conflict": 3}
    out = []
    for cand in candidates:
        nm = cand.get(name_key) or ""
        if not nm:
            out.append({**cand, "isomer_verdict": "unknown",
                        "isomer_note": "candidate has no name to compare"})
            continue
        rep = compare(query, nm)
        out.append({**cand, "isomer_verdict": rep["verdict"],
                    "isomer_conflicts": rep["conflicts"],
                    "resolution": rep["candidate_resolution"]})
    out.sort(key=lambda c: rank.get(c.get("isomer_verdict"), 4))
    return out


CONFLICT_VERDICTS = {"conflict"}
REVIEW_VERDICTS = {"conflict", "class-level", "ambiguous"}
