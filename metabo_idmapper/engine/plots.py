"""DB-matching coverage figures, always emitted by coverage_summary.

Two publication figures (design carried over from the project's hand-authored scripts):
  - upset:       5-DB identifier coverage UpSet over analysable (non-exogenous) features,
                 membership = UNION of all candidate + accepted ids. Writes enriched_xref.tsv.
  - improvement: MetaboAnalyst 1st-pass baseline vs current logic (coverage bars + GEM,
                 plus resolution composition).
These are pure visualizations of the ledger — no values are invented.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DBS = ["kegg", "hmdb", "chebi", "pubchem", "inchikey"]
DB_LABEL = {"kegg": "KEGG", "hmdb": "HMDB", "chebi": "ChEBI",
            "pubchem": "PubChem", "inchikey": "InChIKey"}


def _entries(workdir: str) -> dict:
    return json.loads((Path(workdir) / "midmap_ledger.json").read_text())["entries"]


def _union(e: dict) -> dict:
    ids = {db: set() for db in DBS}
    for c in [e.get("accepted", {})] + e.get("candidates", []):
        for db in DBS:
            if c.get(db):
                ids[db].add(str(c[db]).strip())
    return ids


def upset(workdir: str) -> dict:
    """5-DB UpSet over non-exogenous features. Writes figures/db_matching_upset.png +
    enriched_xref.tsv. Returns per-DB totals and intersection combos."""
    E = _entries(workdir)
    figdir = Path(workdir) / "figures"; figdir.mkdir(exist_ok=True)
    rows = []
    for fid, e in E.items():
        if e.get("final_class") == "exogenous-excluded":
            continue
        ids = _union(e)
        rows.append({"feature_id": fid, "name": e.get("original_name", ""),
                     "final_class": e.get("final_class", ""),
                     **{db: (sorted(ids[db])[0] if ids[db] else "") for db in DBS},
                     **{f"has_{db}": bool(ids[db]) for db in DBS}})
    N = len(rows)
    # enriched_xref.tsv
    with open(Path(workdir) / "enriched_xref.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["feature_id"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    per_db = {db: sum(r[f"has_{db}"] for r in rows) for db in DBS}
    patterns: dict = {}
    for r in rows:
        key = tuple(db for db in DBS if r[f"has_{db}"])
        patterns[key] = patterns.get(key, 0) + 1
    items = sorted(patterns.items(), key=lambda kv: kv[1], reverse=True)
    combos = [k for k, _ in items]; sizes = [v for _, v in items]
    n_dbs, n_col = len(DBS), max(1, len(combos))

    fig = plt.figure(figsize=(max(9, 1.15 * n_col + 3), 6.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 4.6], height_ratios=[2.5, 1.6],
                          wspace=0.32, hspace=0.06)
    ax_top = fig.add_subplot(gs[0, 1]); ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_top)
    ax_left = fig.add_subplot(gs[1, 0], sharey=ax_mat)
    x = range(n_col)
    ax_top.bar(x, sizes, color="#2a7db5", width=0.6)
    for i, sz in enumerate(sizes):
        ax_top.text(i, sz + max(sizes) * 0.01, str(sz), ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
    ax_top.set_ylabel("features in\nintersection"); ax_top.set_ylim(0, max(sizes) * 1.15)
    ax_top.spines[["top", "right"]].set_visible(False)
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    row_order = DBS[::-1]; y_of = {db: i for i, db in enumerate(row_order)}
    for i, combo in enumerate(combos):
        present = set(combo)
        for db in DBS:
            ax_mat.plot(i, y_of[db], "o", ms=11,
                        color="#2a7db5" if db in present else "#dcdcdc", zorder=3)
        ys = [y_of[db] for db in present]
        if len(ys) > 1:
            ax_mat.plot([i, i], [min(ys), max(ys)], color="#2a7db5", lw=2.2, zorder=2)
    ax_mat.set_yticks(range(n_dbs)); ax_mat.set_yticklabels([])
    for db in DBS:
        ax_mat.text(-0.9, y_of[db], DB_LABEL[db], ha="right", va="center",
                    fontsize=11, fontweight="bold", clip_on=False)
    ax_mat.tick_params(axis="y", left=False); ax_mat.set_ylim(-0.6, n_dbs - 0.4)
    ax_mat.set_xlim(-0.6, n_col - 0.4)
    ax_mat.tick_params(axis="x", bottom=False, labelbottom=False)
    for s in ax_mat.spines.values():
        s.set_visible(False)
    tot = [per_db[db] for db in row_order]
    ax_left.barh(range(n_dbs), tot, color="#4f7ea8", height=0.55)
    for i, t in enumerate(tot):
        ax_left.text(t, i, f"{t} ", ha="right", va="center", color="#1b3a52",
                     fontsize=9, fontweight="bold", clip_on=False)
    ax_left.set_xlim(max(tot) * 1.14 if tot and max(tot) else 1, 0)
    ax_left.set_xlabel("set size"); ax_left.set_yticks(range(n_dbs)); ax_left.set_yticklabels([])
    ax_left.spines[["top", "left", "right"]].set_visible(False)
    ax_left.tick_params(axis="y", left=False)
    fig.suptitle(f"Metabolite identifier coverage across 5 DBs (union of candidates; "
                 f"{N} analysable features)", fontsize=12.5, fontweight="bold", y=0.99)
    out = figdir / "db_matching_upset.png"
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    return {"figure": str(out), "n_analysable": N,
            "per_db": {DB_LABEL[k]: v for k, v in per_db.items()},
            "all_five": patterns.get(tuple(DBS), 0)}


def improvement(workdir: str) -> dict:
    """MetaboAnalyst 1st-pass baseline vs current logic. Writes figures/db_matching_improvement.png."""
    E = _entries(workdir)
    figdir = Path(workdir) / "figures"; figdir.mkdir(exist_ok=True)

    def ma(e):
        for c in e.get("candidates", []):
            if c.get("source") == "metaboanalyst":
                return c
        return {}

    N = len(E)
    base = {db: 0 for db in DBS}; after = {db: 0 for db in DBS}
    base_any = after_prim = struct = exog = gem_base = gem_after = 0
    for e in E.values():
        c = ma(e)
        for db in DBS:
            base[db] += bool(c.get(db))
        en = _union(e) if e.get("final_class") != "exogenous-excluded" else {db: set() for db in DBS}
        for db in DBS:
            after[db] += bool(en[db])
        base_any += bool(c.get("kegg") or c.get("hmdb"))
        fc = e.get("final_class")
        if fc in ("KEGG-mapped", "HMDB-mapped"):
            after_prim += 1
        elif fc == "structure-only":
            struct += 1
        elif fc == "exogenous-excluded":
            exog += 1
        if e.get("gem_mam"):
            gem_after += 1
            if c.get("kegg") or c.get("hmdb"):
                gem_base += 1

    metrics = [("primary\n(KEGG|HMDB)", base_any, after_prim), ("KEGG", base["kegg"], after["kegg"]),
               ("HMDB", base["hmdb"], after["hmdb"]), ("ChEBI", base["chebi"], after["chebi"]),
               ("PubChem", base["pubchem"], after["pubchem"]),
               ("InChIKey", base["inchikey"], after["inchikey"]),
               ("Mouse-GEM\nMAM", gem_base, gem_after)]
    labels = [m[0] for m in metrics]; b = [m[1] for m in metrics]; a = [m[2] for m in metrics]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                   gridspec_kw={"width_ratios": [2.5, 1.15]})
    x = np.arange(len(labels)); w = 0.38
    r1 = ax1.bar(x - w/2, b, w, label="MetaboAnalyst only (baseline)", color="#b9c6d1")
    r2 = ax1.bar(x + w/2, a, w, label="+ current logic", color="#2a7db5")
    for rects in (r1, r2):
        for rc in rects:
            h = rc.get_height()
            ax1.text(rc.get_x() + rc.get_width()/2, h + 1.5, str(int(h)), ha="center",
                     va="bottom", fontsize=9, fontweight="bold")
    for i, (bb, aa) in enumerate(zip(b, a)):
        d = aa - bb
        if d != 0:
            ax1.text(i, max(bb, aa) + 9, f"+{d}" if d > 0 else str(d), ha="center",
                     va="bottom", fontsize=9, color="#1a7a3a", fontweight="bold")
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=9.5)
    ax1.set_ylabel(f"features (of {N})"); ax1.set_ylim(0, N * 1.12)
    ax1.set_title("A. Coverage: MetaboAnalyst baseline vs. + current logic",
                  fontweight="bold", loc="left")
    ax1.legend(loc="lower right", frameon=False)
    ax1.spines[["top", "right"]].set_visible(False)
    cats_base = [("mapped", base_any, "#2a7db5"), ("unresolved", N - base_any, "#d0d0d0")]
    cats_after = [("primary\n(KEGG/HMDB)", after_prim, "#2a7db5"),
                  ("structure-only", struct, "#e0a13a"), ("exogenous", exog, "#b0524a"),
                  ("unresolved", N - after_prim - struct - exog, "#d0d0d0")]
    for xi, cats in [(0, cats_base), (1, cats_after)]:
        bottom = 0
        for name, v, col in cats:
            if v == 0:
                continue
            ax2.bar(xi, v, bottom=bottom, width=0.62, color=col)
            ax2.text(xi, bottom + v/2, f"{name}\n{v}", ha="center", va="center",
                     color="white" if col != "#d0d0d0" else "#555", fontsize=8.5,
                     fontweight="bold")
            bottom += v
    ax2.set_xticks([0, 1]); ax2.set_xticklabels(["MetaboAnalyst\nonly", "+ current\nlogic"], fontsize=10)
    ax2.set_ylabel(f"names (of {N})"); ax2.set_ylim(0, N * 1.05)
    ax2.set_title("B. Resolution of all names", fontweight="bold", loc="left")
    ax2.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Metabolite ID mapping — improvement over MetaboAnalyst 1st pass",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = figdir / "db_matching_improvement.png"
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    return {"figure": str(out), "baseline_primary": base_any, "after_primary": after_prim,
            "baseline_hmdb": base["hmdb"], "after_hmdb": after["hmdb"]}
