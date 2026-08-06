"""Emit standalone reproduction artifacts for a run, using the ORIGINAL library APIs.

Two artifacts, ONE implementation:

  code/reproduce_mapping.py     `codegen/raw_engine.py` inlined verbatim + a narrated main()
  code/reproduce_mapping.ipynb  the same flow unrolled into linear cells (no def/lambda) that
                                call `code/raw_engine.py`, shipped beside them like the .R
                                engines

Neither imports metabo_idmapper or any tool wrapper: the flow runs on MetaboAnalystR /
KEGGREST / BridgeDbR through Rscript, PubChem PUG-REST, molmass, COBRApy and matplotlib. The
reasoning layer's harmonized names, disambiguations and exclusions are READ from the saved
`midmap_ledger.json` rather than re-derived, and the stages are executed, not described:

  1 MetaboAnalystR 1st pass             5 back-check every accepted KEGG against its own
  2 extract the unmapped                  keggGet record + scan for shared identifiers
  3 rename-rule re-search with the       6 model crosswalk: xref AND a search over the
    SAVED harmonized names                 model's own names
  4 BridgeDb cross-check + HMDB          7 master table + coverage figure
    backfill

Everything computational lives in `codegen/raw_engine.py` — a real, unit-testable module —
because logic written out once per emitter is logic that drifts: one model base-id bug
previously had to be fixed in three places (engine, notebook cell string, script template).
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import BRIDGE_DB, GEM_MODEL, RSCRIPT_DIR
from .state import Session

_CODEGEN = Path(__file__).parent / "codegen"
_R_SIDECARS = {"ma.R": "metaboanalyst_map.R", "bridge.R": "bridge_xref.R",
               "kegg.R": "kegg_search.R", "keggent.R": "kegg_entry.R"}

# the raw-API surface the narrative calls; bound as `raw.*` in both artifacts
_ENGINE_API = sorted({
    "load_ledger", "metaboanalyst_map", "harmonized_queries", "kegg_find", "pubchem_by_name",
    "mono_mass", "bridge_map", "assemble_rows", "kegg_entries", "apply_kegg_backcheck",
    "shared_id_scan", "build_model_index", "gem_crosswalk", "gem_name_search", "write_master",
    "plot_upset"})


def _md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text}


def _code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text}


def _raw_engine_source() -> str:
    return (_CODEGEN / "raw_engine.py").read_text()


def _write_sidecars(outdir: Path) -> None:
    """The raw engines the artifacts call: the four R scripts + the raw-API Python module."""
    outdir.mkdir(parents=True, exist_ok=True)
    for name, src in _R_SIDECARS.items():
        (outdir / name).write_text((RSCRIPT_DIR / src).read_text())
    (outdir / "raw_engine.py").write_text(_raw_engine_source())


# --------------------------------------------------------------------- the shared narrative
# Written ONCE and reused by both emitters: the script runs it inside main(), the notebook
# runs the same calls cell by cell.
_STAGES = [
    ("01_SETUP", """
# INPUT   : midmap_ledger.json <- this run's ledger (원본 이름 + normalize + candidates
#           [source별 harmonized query] + accepted + decisions). 작성자: driver 세션.
#           BRIDGE_DB(.bridge) / GEM_MODEL(SBML|JSON) / code/*.R / code/raw_engine.py.
# OUTPUT  : LEDGER, names — 이후 모든 STAGE의 입력.
# REUSED  : 전 STAGE.
LEDGER = raw.load_ledger(OUT)
names = [e['original_name'] for e in LEDGER.values()]
print('entries:', len(LEDGER))"""),

    ("02_STAGE 1 — MetaboAnalystR 1st pass (raw)", """
# INPUT   : names(01) -> code/ma.R (raw MetaboAnalystR).
# DOES    : 이름 배치 매핑 -> 화합물당 KEGG/HMDB/PubChem/ChEBI 동시 반환.
# OUTPUT  : stage1[fid] = {kegg,hmdb,pubchem,chebi}.
# REUSED  : STAGE 2(미매핑 판정), STAGE 4(bridge source), assemble.
ma = raw.metaboanalyst_map(names, MA_R, OUT / '_ma_work', RSCRIPT)
stage1 = {}
for fid, e in LEDGER.items():
    row = ma.get(e['normalized'].get('normalized') or e['original_name'], {})
    stage1[fid] = {db: (row.get(db) if row.get(db) not in (None, '', 'NA', 'null') else None)
                   for db in ('kegg', 'hmdb', 'pubchem', 'chebi')}
print('1st-pass primary (KEGG|HMDB):',
      sum(1 for v in stage1.values() if v['kegg'] or v['hmdb']), '/', len(LEDGER))"""),

    ("03_STAGE 2 — extract the unmapped", """
# INPUT   : stage1(02), LEDGER(01).
# DOES    : 1차에서 KEGG·HMDB 둘 다 없는 항목만 추림(xenobiotic 제외).
# OUTPUT  : unmapped (feature id 리스트).
# REUSED  : STAGE 3(rename 재검색 대상).
unmapped = [fid for fid, v in stage1.items()
            if not (v['kegg'] or v['hmdb'])
            and LEDGER[fid].get('final_class') != 'xenobiotic-excluded']
print('unmapped after 1st pass:', len(unmapped))"""),

    ("04_STAGE 3 — rename-rule ID extraction (핵심)", """
# INPUT   : unmapped(03), LEDGER의 candidate.query(= 저장된 harmonized 이름)
#           -> code/kegg.R(raw KEGGREST), PubChem PUG-REST, molmass.
# DOES    : 저장된 rename 이름으로 실제 재검색(최종 ID 붙여넣기가 아님). keggFind 다중 hit는
#           ledger가 고른 KEGG로 확정하고, 질량은 molmass로 대조.
# OUTPUT  : recovered[fid] = {kegg?, pubchem?, inchikey?, mono_mass?}.
# REUSED  : STAGE 4(bridge source), assemble.
recovered = {}
for fid in unmapped:
    e = LEDGER[fid]
    qs = raw.harmonized_queries(e)
    picked = e.get('accepted', {})
    found = {}
    hits = raw.kegg_find(qs['kegg'], KEGG_R, RSCRIPT) if qs['kegg'] else []
    keggs = [h['kegg'] for h in hits]
    if picked.get('kegg') and picked['kegg'] in keggs:
        found['kegg'] = picked['kegg']          # the ledger's disambiguation, re-found here
    elif len(set(keggs)) == 1:
        found['kegg'] = keggs[0]
    for q in qs['pubchem']:
        hit = raw.pubchem_by_name(q)
        if hit:
            found.update(pubchem=hit['cid'], inchikey=hit['inchikey'],
                         mono_mass=raw.mono_mass(hit['formula']))
            break
    if found:
        recovered[fid] = found
print('recovered by rename-rule re-search:', len(recovered))"""),

    ("05_STAGE 4 — cross-check across DBs (BridgeDb) + HMDB backfill", """
# INPUT   : stage1 + recovered 의 InChIKey/KEGG/ChEBI/PubChem -> code/bridge.R (raw BridgeDbR;
#           2 GB derby DB는 호출당 1회 로드되므로 전부 한 번에 배치).
# DOES    : 확보된 ID를 다른 시스템으로 승격하고, 없는 HMDB를 채운다.
# OUTPUT  : backfilled[fid] = HMDB.
# REUSED  : assemble.
SYS = {'inchikey': 'InChIKey', 'kegg': 'KEGG Compound', 'chebi': 'ChEBI',
       'pubchem': 'PubChem-compound'}
queries = []
for fid, e in LEDGER.items():
    ids = {**stage1[fid], **recovered.get(fid, {})}
    if ids.get('hmdb') or e.get('accepted', {}).get('hmdb'):
        continue
    for src in ('inchikey', 'kegg', 'chebi', 'pubchem'):
        val = ids.get(src) or e.get('accepted', {}).get(src)
        if val:
            queries.append({'feature_id': fid, 'id': str(val), 'source': SYS[src],
                            'targets': ['HMDB']})
            break
backfilled = {}
for q, res in zip(queries, raw.bridge_map(queries, BRIDGE_R, BRIDGE_DB, RSCRIPT) or []):
    hits = (res.get('mappings') or {}).get('HMDB') or []
    if hits:
        backfilled[q['feature_id']] = hits[0] if isinstance(hits, list) else hits
print('bridged ids:', len(queries), '| HMDB backfilled:', len(backfilled))"""),

    ("06_ASSEMBLE — 최종 행 구성", """
# INPUT   : stage1 + recovered + backfilled + LEDGER(accepted/final_class).
# DOES    : feature별 최종 ID + class 확정(xenobiotic-excluded는 ID를 비움).
# OUTPUT  : rows — master_ledger.tsv의 원본.
# REUSED  : 백체크, GEM crosswalk, 표/그림.
rows = raw.assemble_rows(LEDGER, stage1, recovered, backfilled)
print('rows:', len(rows))"""),

    ("07_BACK-CHECK — 확정 KEGG를 그 id의 KEGG 레코드와 대조", """
# INPUT   : rows(06) -> code/keggent.R (raw KEGGREST keggGet).
# DOES    : id -> 그 id 자신의 NAME/FORMULA를 받아 이름과 대조한다. name->id 경로로는 절대
#           안 잡히는 오류(맞는 이름에 붙어 온 틀린 id)가 여기서만 드러난다. 검증된 사례:
#           C05472의 KEGG 이름은 'Urocortisol'이며 담즙산이 아니다.
#           같은 id를 두 화합물이 나눠 갖는지도 함께 본다(둘 사이의 비를 1로 붕괴시킨다).
# OUTPUT  : rows[*]['kegg_name','kegg_resolution'], backcheck_rows, shared.
# REUSED  : 리포트의 판정 근거, 종/클래스 해상도 표기.
entries = raw.kegg_entries([r['kegg'] for r in rows if r.get('kegg')], KEGGENT_R, RSCRIPT)
backcheck_rows = raw.apply_kegg_backcheck(rows, entries)
shared = raw.shared_id_scan(rows)
print('back-checked:', len(backcheck_rows),
      '| class-level ids:', sum(1 for b in backcheck_rows if b['resolution'] == 'class'))
print('identifiers held by >1 compound:', shared or 'none')"""),

    ("08_GEM CROSSWALK — xref + 모델 이름 직접 조회", """
# INPUT   : rows(06), GEM_MODEL(SBML|JSON) via COBRApy.
# DOES    : (a) KEGG/HMDB/ChEBI xref 매핑, (b) xref 미스는 모델의 이름/분자식을 직접 조회.
#           GEM은 지질 species에 xref를 거의 달지 않으므로 (b)가 없으면 '모델에 있는데
#           없다'고 오판한다(검증 런: xref 8/27, 실제 보유 21/28).
#           주의: 여기의 (b)는 문자 유사도만 쓴다. MCP의 gem_search는 여기에 shorthand 확장
#           (LysoPC(16:0) -> 1-palmitoyl-...)과 isomer 판정을 더해 순위를 정하므로, 이 목록의
#           1위를 답으로 읽지 말고 후보로 읽을 것 — 채택은 사람/LLM의 판단이다.
# OUTPUT  : rows[*]['gem_mam','gem_route'], gem_name_candidates.json.
# REUSED  : 대사모델 flux(E-Flux/GIMME) 입력.
index = raw.build_model_index(GEM_MODEL)
raw.gem_crosswalk(rows, index)
gem_name_hits = raw.gem_name_search(rows, index)
print('model:', index['n_met'], 'metabolites | xref-mapped:',
      sum(1 for r in rows if r['gem_mam']),
      '| name-search candidates for:', sum(1 for v in gem_name_hits.values() if v))
(OUT / 'gem_name_candidates.json').write_text(json.dumps(gem_name_hits, indent=1))"""),

    ("09_WRITE — master table + coverage figure", """
# INPUT   : rows(06,07,08).
# OUTPUT  : master_ledger.tsv, figures/db_matching_upset.png.
# REUSED  : 대사체 통계/네트워크의 ID 조인, joint pathway(KEGG/HMDB), 리포트.
print('wrote', raw.write_master(rows, OUT / 'master_ledger.tsv'))
print(raw.plot_upset(rows, OUT / 'figures' / 'db_matching_upset.png'))"""),
]


def _setup_block(for_notebook: bool) -> str:
    """Paths, R sources and the engine binding — the only part that differs per artifact."""
    if for_notebook:
        return (
            "import json, os, sys\n"
            "from pathlib import Path\n\n"
            "OUT = Path.cwd()\n"
            "if OUT.name == 'code': OUT = OUT.parent      # run from the workdir or its code/\n"
            "CODE = OUT / 'code'\n"
            "sys.path.insert(0, str(CODE))                # raw_engine.py ships beside the .R\n"
            "import raw_engine as raw\n"
            f"BRIDGE_DB = {BRIDGE_DB!r}\n"
            f"GEM_MODEL = {GEM_MODEL!r}\n"
            "RSCRIPT = os.environ.get('METABO_IDMAP_RSCRIPT', 'Rscript')\n"
            "MA_R = (CODE / 'ma.R').read_text()\n"
            "BRIDGE_R = (CODE / 'bridge.R').read_text()\n"
            "KEGG_R = (CODE / 'kegg.R').read_text()\n"
            "KEGGENT_R = (CODE / 'keggent.R').read_text()\n"
            "print('engines:', sorted(p.name for p in CODE.glob('*.R')))")
    return (
        "import json, os, sys\n"
        "from pathlib import Path\n\n"
        "HERE = Path(__file__).resolve().parent\n"
        "OUT = HERE.parent if HERE.name == 'code' else HERE\n"
        f"BRIDGE_DB = {BRIDGE_DB!r}\n"
        f"GEM_MODEL = {GEM_MODEL!r}\n"
        "RSCRIPT = os.environ.get('METABO_IDMAP_RSCRIPT', 'Rscript')")


def generate(workdir: str) -> str:
    """Write <workdir>/code/reproduce_mapping.py (self-contained) and return its path."""
    s = Session(workdir)
    outdir = Path(s.workdir) / "code"
    _write_sidecars(outdir)

    r_consts = "\n".join(
        f'{var} = r"""{(RSCRIPT_DIR / src).read_text()}"""'
        for var, src in (("MA_R", "metaboanalyst_map.R"), ("BRIDGE_R", "bridge_xref.R"),
                         ("KEGG_R", "kegg_search.R"), ("KEGGENT_R", "kegg_entry.R")))

    body = "\n\n".join(
        "    # ══ {} ═══\n".format(title)
        + "\n".join(("    " + ln) if ln.strip() else "" for ln in code.strip("\n").splitlines())
        for title, code in _STAGES)

    # the inlined engine is at module level; bind it as `raw.*` so the narrative is identical
    # in both artifacts and neither drifts from the other
    binding = ("class raw:\n"
               + "\n".join(f"    {n} = staticmethod({n})" for n in _ENGINE_API))

    parts = [
        '#!/usr/bin/env python\n'
        '"""AUTO-GENERATED reproduction of a metabo-idmapper run — RAW LIBRARY APIs ONLY.\n\n'
        "Runs the real flow: 1) MetaboAnalystR 1st pass  2) extract the unmapped  3) re-search\n"
        "with the harmonized names SAVED in midmap_ledger.json  4) BridgeDb cross-check + HMDB\n"
        "backfill  5) back-check every accepted KEGG against its own keggGet record + scan for\n"
        "identifiers shared by two compounds  6) model crosswalk by xref AND by the model's own\n"
        "names  7) master table + coverage figure.\n\n"
        "No metabo_idmapper import and no tool-function names. The engine below is the same\n"
        "source the package ships as codegen/raw_engine.py, copied in verbatim.\n\n"
        'Run:  python reproduce_mapping.py\n"""',
        _setup_block(for_notebook=False),
        "# ═══ raw R engines (verbatim) ═════════════════════════════════════════════════════\n"
        + r_consts,
        "# ═══ raw-API engine — verbatim copy of metabo_idmapper/codegen/raw_engine.py ══════\n"
        "# Edit that file, not this copy: it is unit-tested and this is generated from it.\n"
        + _raw_engine_source(),
        binding,
        "# ═══ the flow ════════════════════════════════════════════════════════════════════\n"
        "def main():\n" + body + "\n\n\nif __name__ == '__main__':\n    main()",
    ]
    out = outdir / "reproduce_mapping.py"
    out.write_text("\n\n\n".join(parts) + "\n")
    return str(out)


def generate_notebook(workdir: str) -> str:
    """Write <workdir>/code/reproduce_mapping.ipynb — the same flow, one cell per stage.

    Cells stay linear (no `def`, no `lambda`): the computation lives in the `raw_engine.py`
    sidecar written next to the notebook, exactly like the .R engines it already calls.
    """
    s = Session(workdir)
    outdir = Path(s.workdir) / "code"
    _write_sidecars(outdir)

    cells = [_md(
        "# metabo-idmapper — reproduction (unrolled, one cell per stage)\n\n"
        "Raw library APIs only — **MetaboAnalystR / KEGGREST / BridgeDbR** via Rscript, "
        "**PubChem** PUG-REST, **molmass**, **COBRApy**, **matplotlib**. No tool wrappers and "
        "no `def` in the cells: the raw-API functions live in the `raw_engine.py` sidecar "
        "beside this notebook (written from the package's unit-tested `codegen/raw_engine.py`), "
        "the same way `ma.R` / `bridge.R` / `kegg.R` / `keggent.R` are.\n\n"
        "Flow: **1)** MetaboAnalyst 1st pass → **2)** extract the unmapped → **3)** re-search "
        "with the harmonized names READ from the saved ledger → **4)** BridgeDb cross-check + "
        "HMDB backfill → **5)** back-check every KEGG against its own keggGet record + "
        "shared-id scan → **6)** model crosswalk by xref AND by the model's own names → "
        "**7)** master table + coverage figure.\n\n"
        "Run it from the run workdir (or this `code/` directory); it reads "
        "`midmap_ledger.json`.\n\n"
        "## 파일 맵\n"
        "**입력** — `midmap_ledger.json`(driver 세션이 남긴 판정), `code/*.R` + "
        "`code/raw_engine.py`(원 API 엔진), `BRIDGE_DB`, `GEM_MODEL`.\n\n"
        "**산출물** — `master_ledger.tsv`(feature×ID×class×gem_mam), "
        "`gem_name_candidates.json`(xref 미스에 대한 모델 이름 후보), "
        "`figures/db_matching_upset.png`.\n\n"
        "**재사용처** — `gem_mam` → 대사모델 flux 입력; KEGG/HMDB → pathway 매핑; "
        "`master_ledger.tsv` → 통계·네트워크의 ID 조인.")]

    for i, (title, code) in enumerate(_STAGES):
        heading = title.split("_", 1)[1] if "_" in title else title
        cells.append(_md(f"## {i + 1:02d}_{heading}"))
        cell = f"# ══ {title} ═══" + code
        if title.startswith("01_"):
            cell = "# ══ 01_SETUP ═══\n" + _setup_block(for_notebook=True) + "\n" + code
        cells.append(_code(cell))

    nb = {"cells": cells,
          "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                      "name": "python3"},
                       "language_info": {"name": "python"}},
          "nbformat": 4, "nbformat_minor": 5}
    out = outdir / "reproduce_mapping.ipynb"
    out.write_text(json.dumps(nb, ensure_ascii=False, indent=1))
    return str(out)
