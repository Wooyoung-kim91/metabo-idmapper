"""Reporting — the artifacts a run is read through, separated from the tools that decide.

Everything here READS the finished ledger and writes files: the master table, coverage, the
recovery/provenance tables, the model-resolution table, the figures, the raw-API reproduction
code, and the original data file with the ID columns appended. Nothing here makes an identity
call or changes a mapping.

Sequencing lives in `finalize_run`, not inside a computation. `coverage_summary` used to emit
the figures, the provenance tables, the reproduction code and the annotated source as side
effects, which meant "compute the coverage numbers" could not be asked for on its own and the
stage-7 order was an implicit pipeline hidden inside one tool — in a registry whose whole
premise is that ordering is the reasoning layer's job, not the code's.
"""

from __future__ import annotations

import csv
from pathlib import Path

from . import errors as _err
from . import flags as _flags
from .engine import lexicon as _lex
from .state import ALL_CLASSES, PRIMARY_CLASSES, Session


_ROUTE = {
    "kegg-search": "KEGG name/synonym search (keggFind) + verify",
    "pubchem-search": "PubChem name search",
    "chebi-ols4": "ChEBI OLS4 name search",
    "bridgedb": "xref bridge (BridgeDb)",
    "pubchem": "PubChem structure lookup",
    "metaboanalyst": "MetaboAnalyst exact",
}


def mapping_provenance(workdir: str, export: bool = True) -> dict:
    """Recovery provenance for the report: what was mapped to KEGG / HMDB BEYOND the
    MetaboAnalyst 1st pass and by which logic, plus the unmapped (harmonized-name-tried,
    reason) and the two origin buckets. Writes kegg_recovered.tsv, hmdb_recovered.tsv,
    unmapped_harmonization.tsv, exogenous_kept.tsv (biological outside-host, KEPT + origin),
    and xenobiotic_excluded.tsv (non-biological contaminant, EXCLUDED + origin)."""
    s = Session(workdir)

    def ma(e):
        return next((c for c in e.get("candidates", []) if c.get("source") == "metaboanalyst"), {})

    def cand_for(e, db, val):
        return next((c for c in e.get("candidates", []) if str(c.get(db)) == str(val)
                     and c.get("source") != "metaboanalyst"), None)

    def rationales(e):
        return [d.get("rationale", "") for d in e.get("decisions", [])
                if d.get("action") in ("accept", "auto_accept")]

    def harmon_name(e, c):
        # the query that matched, unless it is a bridge-style "src:id" -> use the compound name
        q = (c or {}).get("query", "")
        if not q or (":" in q and q.split(":", 1)[0] in
                     ("kegg", "hmdb", "chebi", "pubchem", "inchikey")):
            return e["normalized"].get("normalized", "")
        return q

    def category_of(e, reason):
        # prefer the recorded origin tag; else fall back to the contaminant lexicon / reason text
        origin = e.get("origin")
        if origin and origin != "endogenous":
            return origin
        cat = _lex.xenobiotic_category(e["original_name"])
        if cat:
            return cat
        r = (reason or "").lower()
        if any(w in r for w in ("antibiotic", "drug", "pharmac", "sulfonamide", "appetite",
                                "antihistamine", "alzheimer")):
            return "drug"
        if any(w in r for w in ("diet", "food", "plant", "dietary")):
            return "diet"
        if any(w in r for w in ("microb", "gut", "bacter")):
            return "microbial"
        if any(w in r for w in ("surfactant", "detergent", "diethanolamide", "sulfonic",
                                "amidopropyl", "alkyl sulfate")):
            return "surfactant"
        if any(w in r for w in ("plasticizer", "phthalate", "phosphate ester", "flame")):
            return "plasticizer"
        if any(w in r for w in ("additive", "trifluoroacetic", "mobile phase", "ion-pair", "ion pairing")):
            return "additive"
        if any(w in r for w in ("antioxidant", "nitrophenol", "tert-butyl", "reagent",
                                "building block", "industrial", "thiocyanate", "synthetic")):
            return "industrial"
        return "unspecified"

    def reason_of(e):
        excl = [d.get("rationale", "") for d in e.get("decisions", [])
                if d.get("action") in ("exclude", "accept", "auto_accept")]
        return excl[-1] if excl else (e.get("decisions", [{}])[-1].get("rationale", ""))

    kegg_rows, hmdb_rows, unmapped_rows, exo_kept_rows, xeno_rows = [], [], [], [], []
    for fid, e in s.entries.items():
        a = e.get("accepted", {})
        m = ma(e)
        ids = ";".join(f"{k}:{a[k]}" for k in ("kegg", "hmdb", "chebi", "pubchem") if a.get(k))
        if e.get("final_class") == "xenobiotic-excluded":
            reason = reason_of(e)
            xeno_rows.append([fid, e["original_name"], category_of(e, reason), ids, reason])
        elif e.get("final_class") == "exogenous":
            reason = reason_of(e)
            xref = "in-model" if e.get("gem_mam") else (e.get("gem_cause") or "")
            exo_kept_rows.append([fid, e["original_name"], category_of(e, reason), ids, xref, reason])
        if a.get("kegg") and not m.get("kegg"):
            c = cand_for(e, "kegg", a["kegg"])
            route = _ROUTE.get(c["source"], c["source"]) if c else "?"
            harmon = harmon_name(e, c)
            logic = next((r for r in rationales(e) if "re-bridge" not in r and "backfill" not in r), "")
            kegg_rows.append([fid, e["original_name"], harmon, a["kegg"],
                              (c or {}).get("note", ""), route, logic])
        if a.get("hmdb") and not m.get("hmdb"):
            c = cand_for(e, "hmdb", a["hmdb"])
            route = _ROUTE.get(c["source"], c["source"]) if c else "?"
            harmon = harmon_name(e, c)
            logic = next((r for r in rationales(e) if "hmdb" in r.lower() or "re-bridge" in r
                          or "backfill" in r.lower()), rationales(e)[0] if rationales(e) else "")
            hmdb_rows.append([fid, e["original_name"], harmon, a["hmdb"], route, logic])
        if e.get("final_class") in ("structure-only", "unmapped"):
            tried = sorted({c.get("query") for c in e.get("candidates", []) if c.get("query")})
            sid = ";".join(f"{k}:{a[k]}" for k in ("chebi", "pubchem", "inchikey") if a.get(k))
            reason = next((r for r in rationales(e)), "") or "no KEGG/HMDB returned by any source"
            unmapped_rows.append([fid, e["original_name"], " | ".join(tried),
                                  e.get("final_class"), sid, reason])

    # GEM curation table: the resolution axis a flux design needs — which species each marker
    # maps to, whether that is the species itself, a class/pool proxy or a surrogate, and which
    # ones the model genuinely does not contain (absence is a result, not a gap in the table).
    gem_rows = []
    for fid, e in s.entries.items():
        rel = e.get("gem_relation")
        if not rel:
            if e.get("gem_mam"):
                rel = "unspecified"      # mapped by an older run, before relations existed
            elif e.get("final_class") in (PRIMARY_CLASSES | {"exogenous"}):
                rel = "not-evaluated"    # model-relevant but never put to the model
            else:
                continue
        rec = (e.get("gem") or {}).get(e.get("gem_model") or "", {})
        det = next((d for d in e.get("decisions", []) if d.get("action") == "gem_assign"), {})
        gem_rows.append([fid, e["original_name"], e.get("gem_model") or "",
                         ";".join(e.get("gem_mam", [])), rel or "",
                         ";".join(rec.get("gem_names", []) or []),
                         "yes" if e.get("gem_curated") else "no",
                         rec.get("rationale") or det.get("rationale") or ""])

    tables = {
        "gem_curation.tsv": (["feature_id", "original_name", "model", "mam", "relation",
                              "gem_name", "curated", "rationale"], gem_rows),
        "kegg_recovered.tsv": (["feature_id", "original_name", "harmonized_name", "kegg",
                                "kegg_name", "route", "logic"], kegg_rows),
        "hmdb_recovered.tsv": (["feature_id", "original_name", "harmonized_name", "hmdb",
                                "route", "logic"], hmdb_rows),
        "unmapped_harmonization.tsv": (["feature_id", "original_name", "harmonized_names_tried",
                                        "final_class", "structure_ids", "reason"], unmapped_rows),
        "exogenous_kept.tsv": (["feature_id", "original_name", "origin", "ids", "gem_xref",
                                "reason"], exo_kept_rows),
        "xenobiotic_excluded.tsv": (["feature_id", "original_name", "origin", "ids", "reason"],
                                    xeno_rows),
    }
    paths = []
    if export:
        for fname, (cols, rows) in tables.items():
            p = Session(workdir).workdir / fname
            with open(p, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(cols)
                w.writerows(rows)
            paths.append(str(p))
    rel_counts: dict[str, int] = {}
    for r in gem_rows:
        rel_counts[r[4] or "?"] = rel_counts.get(r[4] or "?", 0) + 1
    return {"kegg_recovered": len(kegg_rows), "hmdb_recovered": len(hmdb_rows),
            "unmapped": len(unmapped_rows), "exogenous_kept": len(exo_kept_rows),
            "xenobiotic_excluded": len(xeno_rows), "gem_curation": len(gem_rows),
            "gem_by_relation": rel_counts, "exports": paths,
            "note": "recovered = mapped BEYOND MetaboAnalyst 1st pass; route/logic per row. "
                    "exogenous_kept = biological outside-host signal (kept); xenobiotic_excluded "
                    "= non-biological contaminant (dropped). gem_curation = per-entry model "
                    "resolution (exact / class-proxy / isomer-surrogate / model-scope-absent)."}


def annotate_source(workdir: str, source: str | None = None, sheet: "str | int" = 0,
                    header: int = 0, name_column: str | None = None,
                    out: str | None = None) -> dict:
    """Append the final ID-result columns to the ORIGINAL data file — the intensity matrix
    keeps all its columns and gains ID_final_class / ID_confidence / ID_kegg / ID_hmdb /
    ID_chebi / ID_pubchem / ID_inchikey / ID_gem_mam per compound (matched by name).

    Uses the source recorded by ingest_names if `source` is omitted. `header` is the 0-based
    header row (some vendor sheets have a preamble row, e.g. header=1). Writes <source>_annotated
    .xlsx + .tsv into the workdir."""
    import pandas as pd

    s = Session(workdir)
    rec = s.data.get("source", {})
    source = source or rec.get("xlsx")
    if sheet == 0 and rec.get("sheet") is not None:
        sheet = rec["sheet"]
    name_column = name_column or rec.get("column")
    if not source:
        return _err.fail("invalid_argument",
                         "no source file to annotate — pass source=... "
                         "(ingest_names records it when called with xlsx=)")
    df = pd.read_excel(source, sheet_name=sheet, header=header)
    col = name_column
    if col is None or col not in df.columns:
        # pick the column whose values best overlap the ledger names
        names = {e["original_name"].strip().lower() for e in s.entries.values()}
        best, bestn = df.columns[0], -1
        for c in df.columns:
            hit = sum(str(v).strip().lower() in names for v in df[c].dropna())
            if hit > bestn:
                best, bestn = c, hit
        col = best
    idmap = {e["original_name"].strip().lower(): e for e in s.entries.values()}
    add = ["final_class", "origin", "confidence", "kegg", "hmdb", "chebi", "pubchem",
           "inchikey", "gem_mam", "gem_relation"]
    for a in add:
        df["ID_" + a] = ""
    matched = 0
    for i, v in df[col].items():
        e = idmap.get(str(v).strip().lower())
        if not e:
            continue
        matched += 1
        acc = e.get("accepted", {})
        df.at[i, "ID_final_class"] = e.get("final_class", "")
        df.at[i, "ID_origin"] = e.get("origin") or ""
        df.at[i, "ID_confidence"] = e.get("confidence", "")
        for k in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
            df.at[i, "ID_" + k] = acc.get(k, "")
        df.at[i, "ID_gem_mam"] = ";".join(e.get("gem_mam", []))
        df.at[i, "ID_gem_relation"] = e.get("gem_relation") or ""
    stem = str(Path(source).name).rsplit(".", 1)[0]
    out = out or str(Session(workdir).workdir / f"{stem}_annotated.xlsx")
    df.to_excel(out, index=False)
    tsv = out.rsplit(".", 1)[0] + ".tsv"
    df.to_csv(tsv, sep="\t", index=False)
    return {"source": source, "name_column": col, "rows": len(df), "matched": matched,
            "columns_added": ["ID_" + a for a in add], "annotated": [out, tsv],
            "note": "intensity matrix + final ID columns; unmatched rows have empty ID_* fields."}


def export_code(workdir: str) -> dict:
    """Write standalone reproduction code that reproduces the run using the ORIGINAL library
    APIs (MetaboAnalystR/BridgeDbR/KEGGREST via Rscript, PubChem PUG-REST, molmass, COBRApy,
    matplotlib) — NOT the tool wrappers. Emits both code/reproduce_mapping.py (flow-based,
    reads the saved ledger) and code/reproduce_mapping.ipynb (same flow unrolled into linear
    cells, no def, with detailed input/output/reuse comments per cell)."""
    from . import scriptgen
    try:
        script = scriptgen.generate(workdir)
        notebook = scriptgen.generate_notebook(workdir)
    except Exception as e:
        return _err.fail("export_failed", f"scriptgen failed: {type(e).__name__}: {e}")
    return {"script": script, "notebook": notebook,
            "note": "raw-API reproduction; no metabo_idmapper import."}


def export_report_ppt(workdir: str, out: str | None = None) -> dict:
    """Build a PPTX report from the run's outputs + figures (Title · Coverage KPI · Methods ·
    Pipeline · UpSet · Improvement · Recovery cause→fix · KEGG/HMDB recovered · Unmapped ·
    Exogenous · Outputs). Reads coverage_summary.tsv + provenance tsvs + figures; run
    coverage_summary first. Requires python-pptx + pillow."""
    from . import report_ppt
    try:
        path = report_ppt.generate(workdir, out=out)
    except Exception as e:
        return _err.fail("export_failed", f"ppt export failed: {type(e).__name__}: {e}")
    return {"pptx": path, "note": "values read from run artifacts; no fabrication."}


def plot_coverage(workdir: str) -> dict:
    """Always-on DB-matching figures: figures/db_matching_upset.png (5-DB UpSet + enriched_xref.tsv)
    and figures/db_matching_improvement.png (MetaboAnalyst baseline vs current logic)."""
    from .engine import plots as _plots
    try:
        up = _plots.upset(workdir)
        imp = _plots.improvement(workdir)
    except Exception as e:  # never let a plotting error break the run
        return _err.fail("export_failed", f"plot failed: {type(e).__name__}: {e}")
    return {"upset": up, "improvement": imp}


def coverage_summary(workdir: str, export: bool = True) -> dict:
    """Compute the run's coverage numbers and write the two tables that ARE the run's result:
    `master_ledger.tsv` (one row per feature: ids, class, origin, confidence, model species +
    relation, open flags) and `coverage_summary.tsv` (per-class and per-metric percentages,
    including the GEM relation split).

    Computation only. The figures, provenance tables, reproduction code and annotated source
    are separate tools; `finalize_run` is the one that orders them.
    """
    s = Session(workdir)
    counts = s.class_counts()
    n = len(s.entries)
    primary = sum(counts[c] for c in PRIMARY_CLASSES)
    has_kegg = sum(1 for e in s.entries.values() if e.get("accepted", {}).get("kegg"))
    has_hmdb = sum(1 for e in s.entries.values() if e.get("accepted", {}).get("hmdb"))
    gem_mapped = sum(1 for e in s.entries.values() if e.get("gem_mam"))
    exogenous = counts.get("exogenous", 0)
    xenobiotic_excluded = counts.get("xenobiotic-excluded", 0)
    gem_by_relation: dict[str, int] = {}
    for e in s.entries.values():
        rel = e.get("gem_relation")
        if rel:
            gem_by_relation[rel] = gem_by_relation.get(rel, 0) + 1
    # species-level model input = the model carries THIS compound (not a pool/surrogate)
    gem_species_level = gem_by_relation.get("exact", 0)
    summary = {"total": n, "counts": counts,
               "primary_usable": primary, "has_kegg": has_kegg, "has_hmdb": has_hmdb,
               "gem_mapped": gem_mapped, "gem_by_relation": gem_by_relation,
               "gem_species_level": gem_species_level, "exogenous": exogenous,
               "xenobiotic_excluded": xenobiotic_excluded}
    if export:
        led = s.workdir / "master_ledger.tsv"
        cols = ["feature_id", "original_name", "normalized", "final_class", "origin",
                "confidence", "kegg", "hmdb", "chebi", "pubchem", "inchikey",
                "gem_model", "gem_mam", "gem_relation", "gem_cause", "flags",
                "n_candidates", "n_decisions"]
        with open(led, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(cols)
            for fid, e in s.entries.items():
                acc = e.get("accepted", {})
                w.writerow([fid, e["original_name"],
                            e["normalized"].get("normalized", ""),
                            e.get("final_class", ""), e.get("origin") or "",
                            e.get("confidence", ""),
                            acc.get("kegg", ""), acc.get("hmdb", ""), acc.get("chebi", ""),
                            acc.get("pubchem", ""), acc.get("inchikey", ""),
                            e.get("gem_model") or "",
                            ";".join(e.get("gem_mam", [])), e.get("gem_relation") or "",
                            e.get("gem_cause") or "", ";".join(_flags.open_flags(e, s)),
                            len(e.get("candidates", [])), len(e.get("decisions", []))])
        cov = s.workdir / "coverage_summary.tsv"
        with open(cov, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["metric", "value", "denom", "pct"])
            for c in sorted(ALL_CLASSES):
                w.writerow([f"class:{c}", counts[c], n,
                            round(100 * counts[c] / n, 1) if n else 0])
            for name, val in (("primary_usable", primary), ("has_kegg", has_kegg),
                              ("has_hmdb", has_hmdb), ("gem_mapped", gem_mapped),
                              ("gem_species_level", gem_species_level),
                              ("exogenous", exogenous),
                              ("xenobiotic_excluded", xenobiotic_excluded)):
                w.writerow([name, val, n, round(100 * val / n, 1) if n else 0])
            for rel, val in sorted(gem_by_relation.items()):
                w.writerow([f"gem_relation:{rel}", val, n,
                            round(100 * val / n, 1) if n else 0])
        summary["exports"] = [str(led), str(cov)]
    summary["note"] = ("coverage numbers + master_ledger.tsv + coverage_summary.tsv. For the "
                       "full stage-7 output (figures, provenance tables, reproduction code, "
                       "annotated source) call finalize_run.")
    return summary


def finalize_run(workdir: str, what: list[str] | None = None) -> dict:
    """Stage 7: emit the artifacts a finished run is read through — ONE tool, explicit order.

    `what` selects the steps (default: everything except the slide deck):

      coverage          master_ledger.tsv + coverage_summary.tsv
      figures           figures/db_matching_upset.png + _improvement.png + enriched_xref.tsv
      provenance        kegg_recovered / hmdb_recovered / unmapped_harmonization /
                        exogenous_kept / xenobiotic_excluded / gem_curation .tsv
      code              code/reproduce_mapping.py + .ipynb + the raw R and Python engines
      annotated_source  the ORIGINAL data file with the ID columns appended (needs an xlsx
                        ingest); skipped silently when this run has no source file
      ppt               metabolite_id_mapping_report.pptx (opt in: `what=[..., "ppt"]`)

    These were five separate tools; as one they cost the reasoning layer one name instead of
    five, and the ORDER is written down rather than implied. A step that fails is reported in
    place and does not abort the rest — a broken figure must not cost you the ledger table.
    Run `harness_audit` after this.
    """
    s = Session(workdir)
    steps = {
        "coverage": lambda: coverage_summary(workdir),
        "figures": lambda: plot_coverage(workdir),
        "provenance": lambda: mapping_provenance(workdir),
        "code": lambda: export_code(workdir),
        "annotated_source": lambda: annotate_source(workdir),
        "ppt": lambda: export_report_ppt(workdir),
    }
    default = ["coverage", "figures", "provenance", "code", "annotated_source"]
    want = list(what) if what else default
    unknown = [w for w in want if w not in steps]
    if unknown:
        return _err.fail("invalid_argument",
                         f"unknown step(s) {unknown}; choose from {sorted(steps)}")
    if "annotated_source" in want and not s.data.get("source", {}).get("xlsx"):
        want = [w for w in want if w != "annotated_source"]   # nothing was ingested from a file

    out: dict = {"workdir": str(s.workdir)}
    failed = []
    for name in want:
        try:
            out[name] = steps[name]()
        except Exception as ex:                      # one broken step must not lose the rest
            out[name] = _err.fail("export_failed", f"{type(ex).__name__}: {ex}")
        # the emitters catch their own failures too, so read the result rather than the raise
        if isinstance(out[name], dict) and out[name].get("error_code"):
            failed.append(name)
    out["steps_run"] = want
    out["steps_failed"] = failed
    out["note"] = ("stage-7 artifacts emitted. Run harness_audit last: it checks that these "
                   "exist AND that the reasoning layer honored the contract."
                   + (f" FAILED: {out['steps_failed']} — the mapping is unaffected; re-run "
                      "those steps with what=[...]." if failed else ""))
    return out
