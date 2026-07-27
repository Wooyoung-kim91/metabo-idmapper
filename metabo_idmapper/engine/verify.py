"""Formula / monoisotopic-mass verification via molmass, plus xref cross-checks.

Given a proposed (name, formula, mass) and any observed formula/mass, report whether
they are consistent. This is the deterministic gate the LLM must pass before accepting
a fuzzy / synonym / mass-only candidate.
"""

from __future__ import annotations

from molmass import Formula

# adduct -> (charge, delta mass in Da) for [M+adduct] singly-charged ions
ADDUCTS = {
    "[M+H]+": 1.007276,
    "[M-H]-": -1.007276,
    "[M+Na]+": 22.989218,
    "[M+K]+": 38.963158,
    "[M+NH4]+": 18.033823,
    "[M+H-H2O]+": 1.007276 - 18.010565,
    "[M+Cl]-": 34.969402,
    "[M+FA-H]-": 44.998201,  # formic acid adduct
}


def formula_mass(formula: str) -> dict:
    """Monoisotopic + average mass for a chemical formula. {ok, mono, average}."""
    try:
        f = Formula(formula)
        return {"ok": True, "mono": float(f.monoisotopic_mass), "average": float(f.mass),
                "formula": str(f.formula)}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}"}


def verify(proposed_formula: str | None = None,
           proposed_mass: float | None = None,
           observed_formula: str | None = None,
           observed_mass: float | None = None,
           mass_tol_ppm: float = 15.0) -> dict:
    """Cross-check a proposed identity against observed data.

    Returns {formula_match, mass_match, ppm_error, notes[]}. Any field that cannot be
    evaluated (missing input) is returned as None rather than a fabricated pass.
    """
    notes: list[str] = []
    formula_match = None
    if proposed_formula and observed_formula:
        a = formula_mass(proposed_formula)
        b = formula_mass(observed_formula)
        if a["ok"] and b["ok"]:
            formula_match = a["formula"] == b["formula"]
            if not formula_match:
                notes.append(f"formula {a['formula']} != {b['formula']}")
        else:
            notes.append("formula parse failed")

    # resolve the proposed mass from formula if not given
    pmass = proposed_mass
    if pmass is None and proposed_formula:
        fm = formula_mass(proposed_formula)
        if fm["ok"]:
            pmass = fm["mono"]

    mass_match = None
    ppm_error = None
    if pmass is not None and observed_mass is not None:
        ppm_error = (observed_mass - pmass) / pmass * 1e6
        mass_match = abs(ppm_error) <= mass_tol_ppm
        notes.append(f"ppm={ppm_error:.1f} (tol {mass_tol_ppm})")

    return {
        "formula_match": formula_match,
        "mass_match": mass_match,
        "ppm_error": ppm_error,
        "notes": notes,
    }


def mass_candidates(mz: float, adducts: list[str] | None = None,
                    tol_ppm: float = 15.0) -> list[dict]:
    """Given an observed m/z, return neutral-mass hypotheses per adduct.

    This does NOT search a compound DB (that is the LLM's job via search_synonym);
    it deterministically converts m/z -> neutral monoisotopic mass for each adduct so
    the LLM can hand a mass window to a DB search. Weak evidence by construction.
    """
    adducts = adducts or ["[M+H]+", "[M-H]-", "[M+Na]+", "[M+NH4]+"]
    out = []
    for a in adducts:
        if a not in ADDUCTS:
            continue
        neutral = mz - ADDUCTS[a]
        window = neutral * tol_ppm / 1e6
        out.append({
            "adduct": a, "neutral_mono_mass": round(neutral, 5),
            "mass_low": round(neutral - window, 5), "mass_high": round(neutral + window, 5),
        })
    return out
