"""Raw-API engine for the reproduction artifacts — ONE source, two destinations.

This module implements the pipeline with the ORIGINAL libraries (MetaboAnalystR / KEGGREST /
BridgeDbR via Rscript, PubChem PUG-REST, molmass, COBRApy, matplotlib) and has NO dependency
on the package that ships it. `scriptgen` inlines this file verbatim into
`code/reproduce_mapping.py` and ships it beside `code/reproduce_mapping.ipynb` as a sidecar,
exactly like `ma.R` / `bridge.R` / `kegg.R` / `keggent.R`.

Why it exists as a real module: the same logic used to be written out three times — once in
the engine, once in the notebook's cell strings, once in the script template — and a single
bug (the model base-id rule that turned "tdchola_c" into "tdchola_") had to be fixed in all
three. Code that is stringly-duplicated cannot be tested and will drift. This can be imported
and unit-tested, and the emitters copy it rather than restate it.

The R engines are passed in as SOURCE TEXT so this file stays free of build-time
substitution: the script inlines them as constants, the notebook reads the sidecar .R files.
"""

import csv
import json
import os
import re
import subprocess
import tempfile
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

DBS = ["kegg", "hmdb", "chebi", "pubchem", "inchikey"]
GEM_XREF_KEYS = {"kegg": "kegg.compound", "hmdb": "hmdb", "chebi": "chebi"}


# --------------------------------------------------------------------- R bridge
def run_r(script_source, payload, rscript="Rscript", timeout=1200):
    """Run an R engine (given as source text) with a JSON payload; return its parsed JSON.

    The R side writes to the file named by `out` and never to stdout — MetaboAnalystR and
    curl both pollute stdout, so stdout parsing is unreliable. Returns None on any failure;
    the caller decides what a missing stage means.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as sf:
        sf.write(script_source)
        spath = sf.name
    ofile = tempfile.mktemp(suffix=".json")
    try:
        subprocess.run([rscript, "--vanilla", spath],
                       input=json.dumps({**payload, "out": ofile}),
                       capture_output=True, text=True, timeout=timeout)
        if not os.path.exists(ofile) or os.path.getsize(ofile) == 0:
            return None
        return json.loads(Path(ofile).read_text())
    finally:
        for f in (spath, ofile):
            try:
                os.unlink(f)
            except OSError:
                pass


def metaboanalyst_map(names, ma_r, workdir, rscript="Rscript"):
    """RAW MetaboAnalystR 1st pass: name -> KEGG/HMDB/PubChem/ChEBI. {query: row}."""
    res = run_r(ma_r, {"names": list(names), "workdir": str(workdir)}, rscript) or []
    return {r["query"]: r for r in res}


def kegg_find(terms, kegg_r, rscript="Rscript", max_per_term=8):
    """RAW KEGGREST keggFind over query terms -> [{term, kegg, kegg_name}]."""
    return run_r(kegg_r, {"terms": list(terms), "max_per_term": max_per_term}, rscript) or []


def kegg_entries(ids, keggent_r, rscript="Rscript"):
    """RAW KEGGREST keggGet: what each id ITSELF is. {kegg_id: {names, formula, ...}}.

    This is the back-check. Searching name->id can never reveal a wrong id, because the wrong
    id arrives attached to the right name; the id's own NAME field can (C05472 is
    "Urocortisol", not a bile acid).
    """
    res = run_r(keggent_r, {"ids": sorted({str(i) for i in ids if i})}, rscript) or []
    return {r["kegg"]: r for r in res}


def bridge_map(queries, bridge_r, bridge_db, rscript="Rscript"):
    """RAW BridgeDbR xref bridging, batched so the 2 GB database loads once per call."""
    return run_r(bridge_r, {"db": str(bridge_db), "queries": list(queries)}, rscript) or []


# --------------------------------------------------------------------- structure / mass
def pubchem_by_name(name, timeout=30):
    """RAW PubChem PUG-REST: name -> {cid, formula, inchikey, mono_mass} (None on a miss)."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(str(name))
           + "/property/MolecularFormula,InChIKey,MonoisotopicMass/JSON")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            props = json.loads(r.read().decode())["PropertyTable"]["Properties"][0]
    except Exception:
        return None
    return {"cid": str(props.get("CID")), "formula": props.get("MolecularFormula"),
            "inchikey": props.get("InChIKey"),
            "mono_mass": float(props["MonoisotopicMass"]) if props.get("MonoisotopicMass")
            else None}


def mono_mass(formula):
    """RAW molmass monoisotopic mass, or None when the formula cannot be parsed."""
    try:
        from molmass import Formula
        return float(Formula(str(formula)).monoisotopic_mass)
    except Exception:
        return None


# --------------------------------------------------------------------- model (COBRApy)
def base_id(met_id, compartment=None):
    """Strip the compartment suffix from a model metabolite id.

    Two conventions exist and both must work: Mouse-GEM 'MAM01234c' and Recon3D 'tdchola_c'
    / 'met__L_c'. Dropping one trailing lowercase letter — the obvious rule — leaves
    'tdchola_', which is not an accession any model resolves, so it silently breaks the model
    input built from it.
    """
    mid = str(met_id)
    comp = str(compartment or "")
    if comp and mid.endswith(comp) and len(mid) > len(comp):
        return mid[: -len(comp)].rstrip("_") or mid
    m = re.match(r"^(.+?)_([a-z]{1,2})$", mid)
    if m:
        return m.group(1)
    return re.sub(r"[a-z]$", "", mid)


def build_model_index(model_path):
    """Load a COBRA model (SBML or JSON) and index it by xref, name and formula."""
    import cobra
    model = (cobra.io.load_json_model(str(model_path)) if str(model_path).endswith(".json")
             else cobra.io.read_sbml_model(str(model_path)))
    xref = {db: defaultdict(set) for db in GEM_XREF_KEYS}
    meta = {}
    for m in model.metabolites:
        base = base_id(m.id, getattr(m, "compartment", None))
        rec = meta.setdefault(base, {"name": m.name or "", "formula": getattr(m, "formula", None),
                                     "species": set()})
        rec["species"].add(m.id)
        for db, key in GEM_XREF_KEYS.items():
            v = m.annotation.get(key)
            if v is None:
                continue
            for x in ([v] if isinstance(v, str) else v):
                x = str(x)
                if db == "chebi":
                    x2 = x if x.upper().startswith("CHEBI:") else "CHEBI:" + x
                    xref[db][x2].add(base)
                    xref[db][x2.replace("CHEBI:", "")].add(base)
                else:
                    xref[db][x].add(base)
    for rec in meta.values():
        rec["species"] = sorted(rec["species"])
        rec["pool_species"] = bool(rec["formula"] and re.search(r"[RX]", str(rec["formula"])))
    return {"xref": xref, "meta": meta, "n_met": len(model.metabolites),
            "n_rxn": len(model.reactions)}


def gem_crosswalk(rows, index):
    """Map accepted KEGG/HMDB/ChEBI to model species. Sets row['gem_mam'] + row['gem_route']."""
    for r in rows:
        r["gem_mam"] = ""
        r["gem_route"] = ""
        for db in ("kegg", "hmdb", "chebi"):
            hit = index["xref"][db].get(str(r.get(db) or ""))
            if r.get(db) and hit:
                r["gem_mam"] = ";".join(sorted(hit))
                r["gem_route"] = "xref:" + db
                break
    return rows


def trigrams(text):
    flat = re.sub(r"[^a-z0-9]", "", str(text or "").lower())
    return {flat[i:i + 3] for i in range(max(len(flat) - 2, 0))}


def gem_name_search(rows, index, min_dice=0.4, top=5):
    """Candidate species for xref-missed rows, by character-trigram Dice over model names.

    Chemical names are agglutinative: "Glycochenodeoxycholic acid" and the model's
    "Chenodeoxyglycocholate" are the same compound and share no whole token, so token overlap
    finds nothing and n-gram similarity is the only route that reaches them. Candidates only —
    a pool_species (R-group formula) is a CLASS-level input, not the species.
    """
    model_grams = {b: trigrams(v["name"]) for b, v in index["meta"].items() if v["name"]}
    hits = {}
    for r in rows:
        if r.get("gem_mam") or not r.get("name"):
            continue
        qg = trigrams(r.get("harmonized") or r["name"])
        scored = []
        for b, g in model_grams.items():
            if not qg or not g:
                continue
            dice = 2 * len(qg & g) / (len(qg) + len(g))
            if dice > min_dice:
                scored.append((round(dice, 3), b))
        scored.sort(reverse=True)
        hits[r["feature_id"]] = [
            {"mam_base": b, "gem_name": index["meta"][b]["name"],
             "formula": index["meta"][b]["formula"], "dice": d,
             "pool_species": index["meta"][b]["pool_species"]}
            for d, b in scored[:top]]
    return hits


# --------------------------------------------------------------------- checks
def apply_kegg_backcheck(rows, entries):
    """Annotate rows with what their KEGG id resolves to, and whether it is a CLASS entry.

    A generic acyl/alkyl name with no chain spec means the id answers a species-level name
    with a class-level entry ('1-Acyl-sn-glycero-3-phosphocholine' for LysoPC(16:0)).
    """
    out = []
    for r in rows:
        e = entries.get(str(r.get("kegg") or ""))
        if not e or not e.get("found"):
            continue
        names = e["names"] if isinstance(e.get("names"), list) else [e.get("names")]
        names = [n for n in names if n]
        r["kegg_name"] = names[0] if names else ""
        r["kegg_resolution"] = "class" if any(
            re.search(r"\b\d?-?acyl\b|\balkyl\b", n, re.I) and not re.search(r"\d{1,2}:\d", n)
            for n in names) else "species"
        out.append({"feature_id": r["feature_id"], "name": r["name"], "kegg": r["kegg"],
                    "kegg_name": r["kegg_name"], "resolution": r["kegg_resolution"],
                    "formula": e.get("formula")})
    return out


def shared_id_scan(rows):
    """Identifiers held by more than one compound: {(db, id): [names]}.

    A shared id is silent and destroys any ratio between the two features — the verified case
    being two glycolipids merged into one PubChem entry as synonyms.
    """
    seen = defaultdict(list)
    for r in rows:
        for db in DBS:
            if r.get(db):
                seen[(db, str(r[db]))].append(r.get("name"))
    return {k: v for k, v in seen.items() if len(v) > 1}


# --------------------------------------------------------------------- outputs
MASTER_COLS = ["feature_id", "name", "harmonized", "final_class", "confidence",
               "kegg", "kegg_name", "kegg_resolution", "hmdb", "chebi", "pubchem",
               "inchikey", "gem_mam", "gem_route"]


def write_master(rows, out_path):
    """Write the master mapping table (the run's result as a flat file)."""
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(MASTER_COLS)
        for r in rows:
            w.writerow([r.get(c, "") for c in MASTER_COLS])
    return str(out_path)


def membership(rows):
    """Per-row DB membership (analysable rows only) — the UpSet source."""
    out = []
    for r in rows:
        if r.get("final_class") == "xenobiotic-excluded":
            continue
        out.append({"feature_id": r["feature_id"], "name": r.get("name", ""),
                    **{db: bool(r.get(db)) for db in DBS}})
    return out


def plot_upset(rows, fig_path):
    """DB identifier-coverage UpSet over analysable features (matplotlib only)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mem = membership(rows)
    per_db = {db: sum(m[db] for m in mem) for db in DBS}
    patterns = defaultdict(int)
    for m in mem:
        patterns[tuple(db for db in DBS if m[db])] += 1
    items = sorted(patterns.items(), key=lambda kv: kv[1], reverse=True)
    combos = [k for k, _ in items] or [()]
    sizes = [v for _, v in items] or [0]

    fig = plt.figure(figsize=(max(9, 1.15 * len(combos) + 3), 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 4.6], height_ratios=[2.5, 1.6],
                          wspace=0.32, hspace=0.06)
    top = fig.add_subplot(gs[0, 1])
    mat = fig.add_subplot(gs[1, 1], sharex=top)
    left = fig.add_subplot(gs[1, 0], sharey=mat)
    top.bar(range(len(combos)), sizes, color="#2a7db5", width=0.6)
    for i, sz in enumerate(sizes):
        top.text(i, sz, str(sz), ha="center", va="bottom", fontsize=9, fontweight="bold")
    top.set_ylabel("features in\nintersection")
    top.spines[["top", "right"]].set_visible(False)
    top.tick_params(axis="x", bottom=False, labelbottom=False)
    order = DBS[::-1]
    y_of = {db: i for i, db in enumerate(order)}
    for i, combo in enumerate(combos):
        for db in DBS:
            mat.plot(i, y_of[db], "o", ms=11,
                     color="#2a7db5" if db in combo else "#dcdcdc", zorder=3)
        ys = [y_of[db] for db in combo]
        if len(ys) > 1:
            mat.plot([i, i], [min(ys), max(ys)], color="#2a7db5", lw=2.2, zorder=2)
    for db in DBS:
        mat.text(-0.9, y_of[db], db.upper(), ha="right", va="center", fontsize=11,
                 fontweight="bold", clip_on=False)
    mat.set_yticks(range(len(DBS)))
    mat.set_yticklabels([])
    mat.tick_params(axis="both", left=False, bottom=False, labelbottom=False)
    for sp in mat.spines.values():
        sp.set_visible(False)
    tot = [per_db[db] for db in order]
    left.barh(range(len(DBS)), tot, color="#4f7ea8", height=0.55)
    left.set_xlim(max(tot) * 1.14 if tot and max(tot) else 1, 0)
    left.set_xlabel("set size")
    left.set_yticks(range(len(DBS)))
    left.set_yticklabels([])
    left.spines[["top", "left", "right"]].set_visible(False)
    left.tick_params(axis="y", left=False)
    fig.suptitle(f"Metabolite identifier coverage across {len(DBS)} DBs "
                 f"({len(mem)} analysable features)", fontsize=12.5, fontweight="bold")
    Path(fig_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"figure": str(fig_path), "per_db": per_db, "n_analysable": len(mem)}


# --------------------------------------------------------------------- ledger reading
def load_ledger(workdir):
    """The saved reasoning output: name -> normalized -> candidates -> accepted -> class."""
    return json.loads((Path(workdir) / "midmap_ledger.json").read_text())["entries"]


def harmonized_queries(entry):
    """The rename-rule query strings the reasoning layer actually used, per source."""
    out = {"kegg": [], "pubchem": [], "chebi": []}
    for c in entry.get("candidates", []):
        src, q = c.get("source", ""), c.get("query")
        if not q or (":" in q and q.split(":", 1)[0] in DBS):
            continue                      # skip bridge-style "src:id" queries
        for db in out:
            if db in src:
                out[db].append(q)
    nm = entry.get("normalized", {}).get("normalized")
    if nm:
        for db in out:
            out[db].append(nm)
    for db in out:
        out[db] = list(dict.fromkeys(out[db]))
    return out


def assemble_rows(ledger, stage1, recovered=None, backfilled=None):
    """Final rows: 1st-pass ids + recovered ids + the ledger's recorded decisions."""
    recovered = recovered or {}
    backfilled = backfilled or {}
    rows = []
    for fid, e in ledger.items():
        ids = {**(stage1.get(fid) or {}), **(recovered.get(fid) or {})}
        acc = e.get("accepted", {})
        row = {"feature_id": fid, "name": e["original_name"],
               "harmonized": e.get("normalized", {}).get("normalized", ""),
               "final_class": e.get("final_class", "") or "",
               "confidence": e.get("confidence", "") or "",
               "kegg": ids.get("kegg") or acc.get("kegg", "") or "",
               "hmdb": ids.get("hmdb") or backfilled.get(fid) or acc.get("hmdb", "") or "",
               "chebi": acc.get("chebi", "") or "",
               "pubchem": ids.get("pubchem") or acc.get("pubchem", "") or "",
               "inchikey": ids.get("inchikey") or acc.get("inchikey", "") or ""}
        if row["final_class"] == "xenobiotic-excluded":
            for db in DBS:
                row[db] = ""
        rows.append(row)
    return rows
