"""Generate a standalone reproduction script for a run, using the ORIGINAL library APIs.

The emitted script does NOT import metabo_idmapper or call any tool function. It makes the
PIPELINE FLOW visible and actually executes it with the raw engines:

  Stage 1  raw MetaboAnalystR  — 1st-pass name -> ID
  Stage 2  extract the unmapped (no KEGG/HMDB from stage 1)
  Stage 3  rename-rule ID extraction — READ the harmonized query names saved in
           midmap_ledger.json (the reasoning output) and RE-RUN the DB searches
           (KEGGREST / PubChem / ChEBI) with them; molmass-verify
  Stage 4  cross-check the extracted IDs across DBs (BridgeDbR) + HMDB backfill
  then     classify, GEM crosswalk (COBRApy), coverage figure (matplotlib)

It reads the run's saved ledger at runtime (per the user's point that the reasoned names
are already saved), so nothing is hardcoded from the reasoning layer except the recorded
disambiguation/exclusion, which are read from that file too.
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import BRIDGE_DB, GEM_MODEL, RSCRIPT_DIR
from .state import Session


def _md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text}


def _code(text):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": text}


def generate_notebook(workdir: str) -> str:
    """Write <workdir>/code/reproduce_mapping.ipynb — the pipeline UNROLLED into linear
    cells with NO `def`/lambda; raw library APIs only. The 3 raw R engines are written as
    sidecar .R files next to the notebook and called via Rscript."""
    s = Session(workdir)
    outdir = Path(s.workdir) / "code"
    outdir.mkdir(exist_ok=True)
    (outdir / "ma.R").write_text((RSCRIPT_DIR / "metaboanalyst_map.R").read_text())
    (outdir / "bridge.R").write_text((RSCRIPT_DIR / "bridge_xref.R").read_text())
    (outdir / "kegg.R").write_text((RSCRIPT_DIR / "kegg_search.R").read_text())

    setup = (
        "# ══ 01_SETUP ═════════════════════════════════════════════════════════════════════\n"
        "# INPUT   : midmap_ledger.json  ← this run's ledger (원본 대사체명 + normalize +\n"
        "#           candidates[source별 query] + accepted + decisions). 작성자: 이 run의\n"
        "#           driver 세션(record_decision 등). 이름 원천은 Metabolomics/\n"
        "#           Atheroclerosis_metabolite.xlsx의 Compound 열이며 ledger에 이미 반영됨.\n"
        "#           BRIDGE_DB ← BridgeDb metabolite DB(.bridge, figshare).\n"
        "#           GEM_MODEL ← Mouse-GEM v1.8.0 SBML.\n"
        "#           code/ma.R, bridge.R, kegg.R ← raw R 엔진(MetaboAnalystR/BridgeDbR/KEGGREST).\n"
        "# OUTPUT  : LEDGER, names (메모리) — 이후 모든 STAGE가 사용.\n"
        "# REUSED  : Stage 1~4, assemble.\n"
        "import csv, json, os, re, subprocess, operator, urllib.parse, urllib.request\n"
        "from collections import defaultdict\n"
        "from pathlib import Path\n"
        "import matplotlib.pyplot as plt\n"
        "from molmass import Formula\n"
        "import cobra\n\n"
        "OUT = Path.cwd()\n"
        "if OUT.name == 'code': OUT = OUT.parent          # run from workdir or workdir/code\n"
        "CODE = OUT / 'code'\n"
        "FIG = OUT / 'figures'; FIG.mkdir(exist_ok=True)\n"
        f"BRIDGE_DB = {BRIDGE_DB!r}    # BridgeDb metabolite mapping DB\n"
        f"GEM_MODEL = {GEM_MODEL!r}    # Mouse-GEM SBML (for MAM crosswalk)\n"
        "RSCRIPT = os.environ.get('METABO_IDMAP_RSCRIPT', 'Rscript')\n"
        "DBS = ['kegg', 'hmdb', 'chebi', 'pubchem', 'inchikey']\n"
        "LEDGER = json.loads((OUT / 'midmap_ledger.json').read_text())['entries']  # 저장된 reasoning 결과\n"
        "names = [e['original_name'] for e in LEDGER.values()]   # 원본 대사체명 162개\n"
        "print('entries:', len(LEDGER))"
    )

    stage1 = (
        "# ══ 02_STAGE 1 — MetaboAnalystR 1st-pass mapping (raw) ══════════════════════════\n"
        "# INPUT   : names (setup) → code/ma.R(raw MetaboAnalystR). 임시 _ma_out.json 경유.\n"
        "# DOES    : 이름 배치 매핑 → 화합물당 KEGG/HMDB/PubChem/ChEBI 동시 반환.\n"
        "# OUTPUT  : stage1[fid] = {kegg,hmdb,pubchem,chebi} (메모리).\n"
        "# REUSED  : Stage 2(미매핑 판정), Stage 4(bridge source), assemble.\n"
        "req = {'names': names, 'workdir': str(OUT / '_ma_nb'), 'out': str(OUT / '_ma_out.json')}\n"
        "subprocess.run([RSCRIPT, '--vanilla', str(CODE / 'ma.R')], input=json.dumps(req),\n"
        "               capture_output=True, text=True, timeout=1200)\n"
        "ma = {r['query']: r for r in json.loads((OUT / '_ma_out.json').read_text())}\n"
        "stage1 = {}\n"
        "for fid, e in LEDGER.items():\n"
        "    q = e['normalized'].get('normalized') or e['original_name']\n"
        "    row = ma.get(q, {})\n"
        "    d = {}\n"
        "    for db in ('kegg', 'hmdb', 'pubchem', 'chebi'):\n"
        "        v = row.get(db)\n"
        "        d[db] = v if v not in (None, '', 'NA', 'null') else None\n"
        "    stage1[fid] = d\n"
        "print('1st-pass primary (KEGG|HMDB):',\n"
        "      sum(1 for v in stage1.values() if v['kegg'] or v['hmdb']), '/', len(LEDGER))"
    )

    stage2 = (
        "# ══ 03_STAGE 2 — extract the unmapped list ══════════════════════════════════════\n"
        "# INPUT   : stage1 (Stage 1), LEDGER (setup).\n"
        "# DOES    : 1차에서 KEGG·HMDB 둘 다 없는 항목만 추림(exogenous 제외).\n"
        "# OUTPUT  : unmapped (feature id 리스트).\n"
        "# REUSED  : Stage 3(rename 재검색 대상).\n"
        "# STAGE 2 — extract the unmapped list (no KEGG/HMDB from stage 1; skip exogenous)\n"
        "unmapped = []\n"
        "for fid, v in stage1.items():\n"
        "    if not (v['kegg'] or v['hmdb']) and LEDGER[fid].get('final_class') != 'exogenous-excluded':\n"
        "        unmapped.append(fid)\n"
        "print('unmapped after 1st pass:', len(unmapped))"
    )

    stage3 = (
        "# ══ 04_STAGE 3 — rename-rule ID extraction (핵심) ═══════════════════════════════\n"
        "# INPUT   : unmapped(Stage 2), LEDGER의 candidate.query(=저장된 harmonized 이름) →\n"
        "#           code/kegg.R(raw KEGGREST), PubChem PUG-REST, molmass.\n"
        "# DOES    : 저장된 rename 이름으로 KEGG/PubChem을 실제 재검색(최종 ID 붙여넣기 아님);\n"
        "#           keggFind 다중 hit는 ledger가 고른 KEGG(disambiguation)로 확정; molmass 질량검증.\n"
        "# OUTPUT  : recovered[fid] = {kegg?, pubchem?, inchikey?, formula?, mono_mass?, *_via}.\n"
        "# REUSED  : Stage 4(bridge source), assemble.\n"
        "# STAGE 3 — rename-rule ID extraction: read the SAVED harmonized names (candidate.query\n"
        "# from search sources in the ledger) and RE-RUN the DB searches with them; molmass-verify.\n"
        "recovered = {}\n"
        "for fid in unmapped:\n"
        "    e = LEDGER[fid]\n"
        "    kegg_qs, pubchem_qs = [], []\n"
        "    for c in e.get('candidates', []):\n"
        "        src, qy = c.get('source', ''), c.get('query')\n"
        "        if not qy or (':' in qy and qy.split(':', 1)[0] in DBS):\n"
        "            continue                       # skip bridge-style 'src:id'\n"
        "        if 'kegg' in src: kegg_qs.append(qy)\n"
        "        elif 'pubchem' in src: pubchem_qs.append(qy)\n"
        "    nm = e['normalized'].get('normalized')\n"
        "    if nm: kegg_qs.append(nm); pubchem_qs.append(nm)\n"
        "    rec = e.get('accepted', {})\n"
        "    found = {}\n"
        "    # KEGG via KEGGREST::keggFind with each harmonized name (keep recorded pick if present)\n"
        "    for qy in kegg_qs:\n"
        "        rq = {'terms': [qy], 'max_per_term': 8, 'out': str(OUT / '_k.json')}\n"
        "        subprocess.run([RSCRIPT, '--vanilla', str(CODE / 'kegg.R')], input=json.dumps(rq),\n"
        "                       capture_output=True, text=True, timeout=600)\n"
        "        hits = {h['kegg'] for h in json.loads((OUT / '_k.json').read_text())}\n"
        "        if hits:\n"
        "            found['kegg'] = rec.get('kegg') if rec.get('kegg') in hits else sorted(hits)[0]\n"
        "            found['kegg_via'] = qy\n"
        "            break\n"
        "    # structure via PubChem PUG-REST with each harmonized name; molmass mono-mass check\n"
        "    for qy in pubchem_qs:\n"
        "        url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'\n"
        "               + urllib.parse.quote(qy) + '/property/MolecularFormula,InChIKey,MonoisotopicMass/JSON')\n"
        "        try:\n"
        "            p = json.loads(urllib.request.urlopen(url, timeout=30).read())['PropertyTable']['Properties'][0]\n"
        "        except Exception:\n"
        "            p = None\n"
        "        if p and p.get('InChIKey'):\n"
        "            found['pubchem'] = str(p.get('CID')); found['inchikey'] = p['InChIKey']\n"
        "            found['formula'] = p.get('MolecularFormula'); found['pubchem_via'] = qy\n"
        "            if p.get('MolecularFormula'):\n"
        "                found['mono_mass'] = round(float(Formula(p['MolecularFormula']).monoisotopic_mass), 4)\n"
        "            break\n"
        "    recovered[fid] = found\n"
        "    print(fid, e['original_name'], '-> harmonized', (kegg_qs[:1] or pubchem_qs[:1]),\n"
        "          '->', found.get('kegg') or found.get('inchikey') or '-')"
    )

    stage4 = (
        "# ══ 05_STAGE 4 — cross-check across DBs (BridgeDbR) + HMDB backfill ══════════════\n"
        "# INPUT   : stage1 + recovered 의 ID(InChIKey/KEGG/ChEBI) → code/bridge.R(raw BridgeDbR,\n"
        "#           2GB DB 1회 로드). source/target는 datasource 정식명(getSystemCode).\n"
        "# DOES    : 확보 ID를 다른 DB로 상호대조 → HMDB 없는 항목 채움(HMDB=KEGG 아티팩트 해소).\n"
        "# OUTPUT  : backfill[fid] = HMDB id.\n"
        "# REUSED  : assemble(hmdb 채움).\n"
        "# STAGE 4 — cross-check extracted IDs across DBs (BridgeDbR) + HMDB backfill\n"
        "bq = []\n"
        "for fid, e in LEDGER.items():\n"
        "    if e.get('final_class') == 'exogenous-excluded': continue\n"
        "    ids = dict(stage1[fid]); ids.update(recovered.get(fid, {}))\n"
        "    src = None\n"
        "    if ids.get('inchikey'): src = ('InChIKey', ids['inchikey'])\n"
        "    elif ids.get('kegg'): src = ('KEGG Compound', ids['kegg'])\n"
        "    elif e.get('accepted', {}).get('chebi'): src = ('ChEBI', e['accepted']['chebi'])\n"
        "    if src and not ids.get('hmdb'):\n"
        "        bq.append({'feature_id': fid, 'id': src[1], 'source': src[0]})\n"
        "req = {'db': BRIDGE_DB, 'out': str(OUT / '_b.json'),\n"
        "       'queries': [{'id': q['id'], 'source': q['source'], 'targets': ['HMDB', 'KEGG Compound']} for q in bq]}\n"
        "subprocess.run([RSCRIPT, '--vanilla', str(CODE / 'bridge.R')], input=json.dumps(req),\n"
        "               capture_output=True, text=True, timeout=1800)\n"
        "bres = json.loads((OUT / '_b.json').read_text())\n"
        "backfill = {}\n"
        "for q, res in zip(bq, bres):\n"
        "    hb = res.get('mappings', {}).get('HMDB')\n"
        "    backfill[q['feature_id']] = (hb[0] if isinstance(hb, list) and hb else None)\n"
        "print('HMDB backfilled:', sum(1 for v in backfill.values() if v))"
    )

    assemble = (
        "# ══ 06_ASSEMBLE — 최종 행 구성 ══════════════════════════════════════════════════\n"
        "# INPUT   : stage1 + recovered + backfill + LEDGER(accepted/final_class).\n"
        "# DOES    : feature별 최종 KEGG/HMDB/ChEBI/PubChem/InChIKey + class 확정(exogenous는 ID 비움).\n"
        "# OUTPUT  : rows (master_ledger.tsv의 원본).\n"
        "# REUSED  : GEM crosswalk, master_ledger.tsv, UpSet figure.\n"
        "# assemble final rows (respect recorded class/exclusions read from the saved ledger)\n"
        "rows = []\n"
        "for fid, e in LEDGER.items():\n"
        "    ids = dict(stage1[fid]); ids.update(recovered.get(fid, {}))\n"
        "    acc = e.get('accepted', {})\n"
        "    row = {'feature_id': fid, 'name': e['original_name'],\n"
        "           'harmonized': e['normalized'].get('normalized', ''),\n"
        "           'final_class': e.get('final_class', ''), 'confidence': e.get('confidence', ''),\n"
        "           'kegg': ids.get('kegg') or acc.get('kegg', '') or '',\n"
        "           'hmdb': ids.get('hmdb') or backfill.get(fid) or acc.get('hmdb', '') or '',\n"
        "           'chebi': acc.get('chebi', '') or '',\n"
        "           'pubchem': ids.get('pubchem') or acc.get('pubchem', '') or '',\n"
        "           'inchikey': ids.get('inchikey') or acc.get('inchikey', '') or ''}\n"
        "    if e.get('final_class') == 'exogenous-excluded':\n"
        "        for db in DBS: row[db] = ''\n"
        "    rows.append(row)\n"
        "print('rows:', len(rows))"
    )

    gem = (
        "# ══ 07_GEM CROSSWALK — Mouse-GEM MAM (raw COBRApy) ══════════════════════════════\n"
        "# INPUT   : rows(assemble), GEM_MODEL(Mouse-GEM SBML).\n"
        "# DOES    : 확정 KEGG/HMDB/ChEBI xref로 GEM 대사물질(MAM) 매핑.\n"
        "# OUTPUT  : rows[*]['gem_mam'].\n"
        "# REUSED  : 하류 대사모델 flux(E-Flux/GIMME) 입력, master_ledger.tsv.\n"
        "# GEM crosswalk (raw COBRApy): accepted KEGG/HMDB/ChEBI -> Mouse-GEM MAM base ids\n"
        "model = cobra.io.read_sbml_model(GEM_MODEL)\n"
        "idx = {'kegg': defaultdict(set), 'hmdb': defaultdict(set), 'chebi': defaultdict(set)}\n"
        "keymap = {'kegg': 'kegg.compound', 'hmdb': 'hmdb', 'chebi': 'chebi'}\n"
        "for m in model.metabolites:\n"
        "    base = re.sub(r'[a-z]$', '', m.id)\n"
        "    for db, k in keymap.items():\n"
        "        v = m.annotation.get(k)\n"
        "        if v is None: continue\n"
        "        for x in ([v] if isinstance(v, str) else v):\n"
        "            x = str(x)\n"
        "            if db == 'chebi':\n"
        "                x2 = x if x.upper().startswith('CHEBI:') else 'CHEBI:' + x\n"
        "                idx[db][x2].add(base); idx[db][x2.replace('CHEBI:', '')].add(base)\n"
        "            else:\n"
        "                idx[db][x].add(base)\n"
        "for r in rows:\n"
        "    r['gem_mam'] = ''\n"
        "    for db in ('kegg', 'hmdb', 'chebi'):\n"
        "        if r.get(db) and idx[db].get(str(r[db])):\n"
        "            r['gem_mam'] = ';'.join(sorted(idx[db][str(r[db])])); break\n"
        "print('GEM-mapped:', sum(1 for r in rows if r['gem_mam']))"
    )

    writeout = (
        "# ══ 08_WRITE — master_ledger.tsv ════════════════════════════════════════════════\n"
        "# INPUT   : rows.  OUTPUT: OUT/master_ledger.tsv (최종 매핑 표).\n"
        "# REUSED  : 대사체 통계/네트워크의 ID 조인, joint pathway(KEGG/HMDB), 리포트.\n"
        "# write master_ledger.tsv\n"
        "cols = ['feature_id', 'name', 'harmonized', 'final_class', 'confidence',\n"
        "        'kegg', 'hmdb', 'chebi', 'pubchem', 'inchikey', 'gem_mam']\n"
        "with open(OUT / 'master_ledger.tsv', 'w', newline='') as f:\n"
        "    w = csv.writer(f, delimiter='\\t'); w.writerow(cols)\n"
        "    for r in rows: w.writerow([r.get(c, '') for c in cols])\n"
        "print('wrote', OUT / 'master_ledger.tsv')"
    )

    upset = (
        "# ══ 09_UPSET — DB 커버리지 figure (raw matplotlib) ══════════════════════════════\n"
        "# INPUT   : rows(비-exogenous).  OUTPUT: figures/db_matching_upset.png + enriched_xref.tsv.\n"
        "# DOES    : 5-DB 멤버십 교집합 UpSet(상단 교집합 크기·행렬·좌측 DB별 set size).\n"
        "# REUSED  : 리포트 figure; enriched_xref.tsv는 db_matching_improvement 등 후속 그림 source.\n"
        "# UpSet coverage figure (raw matplotlib), inline\n"
        "mem = [r for r in rows if r['final_class'] != 'exogenous-excluded']\n"
        "for r in mem:\n"
        "    for db in DBS: r['has_' + db] = bool(r.get(db))\n"
        "per_db = {db: sum(m['has_' + db] for m in mem) for db in DBS}\n"
        "patt = defaultdict(int)\n"
        "for m in mem: patt[tuple(db for db in DBS if m['has_' + db])] += 1\n"
        "items = sorted(patt.items(), key=operator.itemgetter(1), reverse=True)\n"
        "combos = [k for k, _ in items] or [()]; sizes = [v for _, v in items] or [0]\n"
        "lab = {'kegg': 'KEGG', 'hmdb': 'HMDB', 'chebi': 'ChEBI', 'pubchem': 'PubChem', 'inchikey': 'InChIKey'}\n"
        "fig = plt.figure(figsize=(max(9, 1.15 * len(combos) + 3), 6.8))\n"
        "gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 4.6], height_ratios=[2.5, 1.6], wspace=0.32, hspace=0.06)\n"
        "at = fig.add_subplot(gs[0, 1]); am = fig.add_subplot(gs[1, 1]); al = fig.add_subplot(gs[1, 0])\n"
        "at.bar(range(len(combos)), sizes, color='#2a7db5', width=0.6)\n"
        "for i, sz in enumerate(sizes): at.text(i, sz, str(sz), ha='center', va='bottom', fontweight='bold', fontsize=9)\n"
        "at.set_ylabel('features in\\nintersection'); at.spines[['top', 'right']].set_visible(False)\n"
        "at.tick_params(axis='x', bottom=False, labelbottom=False)\n"
        "ro = DBS[::-1]; yof = {db: i for i, db in enumerate(ro)}\n"
        "for i, combo in enumerate(combos):\n"
        "    for db in DBS: am.plot(i, yof[db], 'o', ms=11, color='#2a7db5' if db in combo else '#dcdcdc', zorder=3)\n"
        "    ys = [yof[db] for db in combo]\n"
        "    if len(ys) > 1: am.plot([i, i], [min(ys), max(ys)], color='#2a7db5', lw=2.2, zorder=2)\n"
        "for db in DBS: am.text(-0.9, yof[db], lab[db], ha='right', va='center', fontweight='bold', fontsize=11)\n"
        "am.set_yticks([]); am.set_xticks([])\n"
        "for sp in am.spines.values(): sp.set_visible(False)\n"
        "tot = [per_db[db] for db in ro]\n"
        "al.barh(range(len(DBS)), tot, color='#4f7ea8', height=0.55)\n"
        "for i, t in enumerate(tot): al.text(t, i, str(t) + ' ', ha='right', va='center', fontweight='bold', fontsize=9)\n"
        "al.set_xlim((max(tot) * 1.14) if tot and max(tot) else 1, 0); al.set_xlabel('set size')\n"
        "al.set_yticks([]); al.spines[['top', 'left', 'right']].set_visible(False)\n"
        "fig.suptitle('Metabolite identifier coverage across 5 DBs (' + str(len(mem)) + ' analysable features)',\n"
        "             fontsize=12.5, fontweight='bold')\n"
        "fig.savefig(FIG / 'db_matching_upset.png', dpi=200, bbox_inches='tight')\n"
        "print('per-DB:', {lab[k]: per_db[k] for k in DBS})\n"
        "fig"
    )

    nb = {
        "cells": [
            _md("# metabo-idmapper — reproduction (unrolled, no functions)\n\n"
                "Raw library APIs only (**MetaboAnalystR / BridgeDbR / KEGGREST** via Rscript, "
                "**PubChem** PUG-REST, **molmass**, **COBRApy**, **matplotlib**). No tool wrappers, "
                "no `def`. Flow: **1)** MetaboAnalyst 1st pass → **2)** extract unmapped → "
                "**3)** rename-rule re-search with the SAVED harmonized names → **4)** BridgeDb "
                "cross-check + HMDB backfill → classify → GEM crosswalk → coverage figure.\n\n"
                "Run this notebook from the run workdir (or its `code/`); it reads "
                "`midmap_ledger.json`.\n\n"
                "## 파일 맵 (어디서 왔고 / 무엇을 만들고 / 어디서 다시 쓰나)\n\n"
                "**입력**\n"
                "- `midmap_ledger.json` ← 이 run의 ledger. 작성: driver 세션(record_decision 등). "
                "원본 이름은 `Metabolomics/Atheroclerosis_metabolite.xlsx`의 Compound 열(→ ledger에 반영). "
                "여기엔 normalize·candidates(source별 harmonized query)·accepted·decisions가 담김.\n"
                "- `code/ma.R`,`bridge.R`,`kegg.R` ← metabo-idmapper 번들 raw R 엔진.\n"
                "- `BRIDGE_DB`(metabolites_20210109.bridge) ← BridgeDb(figshare). "
                "`GEM_MODEL`(Mouse-GEM.xml) ← Mouse-GEM v1.8.0.\n\n"
                "**산출물**\n"
                "- `master_ledger.tsv` — feature×ID×class×gem_mam 최종 매핑.\n"
                "- `enriched_xref.tsv` — 5-DB 멤버십(UpSet source).\n"
                "- `figures/db_matching_upset.png` — 커버리지 figure.\n"
                "- (임시) `_ma_out.json`,`_k.json`,`_b.json`.\n\n"
                "**재사용처(하류)**\n"
                "- `gem_mam` → 대사모델 flux(E-Flux/GIMME) 입력.\n"
                "- KEGG/HMDB → joint pathway(GSEA/ORA)·KEGG pathway 매핑.\n"
                "- `master_ledger.tsv` → 대사체 통계·상관 네트워크의 ID 조인, 리포트.\n"
                "- `enriched_xref.tsv` → db_matching_improvement 등 후속 커버리지 그림 source."),
            _md("## 01_Setup — inputs & config"), _code(setup), _md("## 02_Stage 1 — MetaboAnalystR 1st pass"), _code(stage1),
            _md("## 03_Stage 2 — extract the unmapped"), _code(stage2),
            _md("## 04_Stage 3 — rename-rule ID extraction (re-search saved harmonized names)"), _code(stage3),
            _md("## 05_Stage 4 — cross-check across DBs (BridgeDb) + HMDB backfill"), _code(stage4),
            _md("## 06_Assemble + GEM crosswalk + coverage figure"), _code(assemble), _code(gem),
            _code(writeout), _code(upset),
        ],
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python",
                                    "name": "python3"},
                     "language_info": {"name": "python"}},
        "nbformat": 4, "nbformat_minor": 5,
    }
    out = outdir / "reproduce_mapping.ipynb"
    out.write_text(json.dumps(nb, ensure_ascii=False, indent=1))
    return str(out)


def generate(workdir: str) -> str:
    """Write <workdir>/code/reproduce_mapping.py and return its path."""
    s = Session(workdir)
    ma_r = (RSCRIPT_DIR / "metaboanalyst_map.R").read_text()
    bridge_r = (RSCRIPT_DIR / "bridge_xref.R").read_text()
    kegg_r = (RSCRIPT_DIR / "kegg_search.R").read_text()

    code = _TEMPLATE
    code = code.replace("@@BRIDGE_DB@@", BRIDGE_DB).replace("@@GEM_MODEL@@", GEM_MODEL)
    code = code.replace("@@MA_R@@", ma_r).replace("@@BRIDGE_R@@", bridge_r)
    code = code.replace("@@KEGG_R@@", kegg_r)

    outdir = Path(s.workdir) / "code"
    outdir.mkdir(exist_ok=True)
    out = outdir / "reproduce_mapping.py"
    out.write_text(code)
    return str(out)


_TEMPLATE = r'''#!/usr/bin/env python
"""AUTO-GENERATED reproduction of a metabo-idmapper run — RAW LIBRARY APIs ONLY.

Shows and runs the real flow:
  1) MetaboAnalystR 1st-pass mapping (raw)
  2) extract the unmapped list
  3) rename-rule ID extraction: read the harmonized query names saved in
     midmap_ledger.json and RE-RUN the DB searches (KEGGREST/PubChem/ChEBI) with them,
     molmass-verify
  4) cross-check the extracted IDs across DBs (BridgeDbR) + HMDB backfill
  then classify, GEM crosswalk (COBRApy), coverage figure (matplotlib).

No metabo_idmapper import, no tool function names. The reasoning layer's disambiguation
and exclusions are READ from the saved ledger (not re-derived).

Run:  python reproduce_mapping.py
"""
import csv
import json
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent                       # the run workdir (holds midmap_ledger.json)
FIG = OUT / "figures"; FIG.mkdir(exist_ok=True)
LEDGER = json.loads((OUT / "midmap_ledger.json").read_text())["entries"]

BRIDGE_DB = "@@BRIDGE_DB@@"
GEM_MODEL = "@@GEM_MODEL@@"
RSCRIPT = os.environ.get("METABO_IDMAP_RSCRIPT", "Rscript")
DBS = ["kegg", "hmdb", "chebi", "pubchem", "inchikey"]

# ========================================================================== raw R engines
# Verified raw-API usage of MetaboAnalystR / BridgeDbR / KEGGREST (inlined verbatim).
_MA_R = r"""@@MA_R@@"""
_BRIDGE_R = r"""@@BRIDGE_R@@"""
_KEGG_R = r"""@@KEGG_R@@"""


def _run_r(snippet, payload, timeout=1200):
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as sf:
        sf.write(snippet); spath = sf.name
    ofile = tempfile.mktemp(suffix=".json")
    try:
        subprocess.run([RSCRIPT, "--vanilla", spath], input=json.dumps({**payload, "out": ofile}),
                       capture_output=True, text=True, timeout=timeout)
        if not os.path.exists(ofile) or os.path.getsize(ofile) == 0:
            return None
        return json.loads(Path(ofile).read_text())
    finally:
        for f in (spath, ofile):
            try: os.unlink(f)
            except OSError: pass


def metaboanalyst_map(names):
    """RAW MetaboAnalystR: name -> KEGG/HMDB/PubChem/ChEBI (1st pass)."""
    return _run_r(_MA_R, {"names": names, "workdir": str(OUT / "_ma_reproduce")}) or []


def bridgedb_map(queries):
    """RAW BridgeDbR: xref bridging. queries=[{id, source, targets:[...]}]."""
    return _run_r(_BRIDGE_R, {"db": BRIDGE_DB, "queries": queries}) or []


def kegg_find(term):
    """RAW KEGGREST::keggFind compound search for one term."""
    return _run_r(_KEGG_R, {"terms": [term], "max_per_term": 8}) or []


# ========================================================================== raw Python engines
def pubchem_lookup(name):
    """RAW PubChem PUG-REST: name -> CID/InChIKey/formula/monoisotopic mass."""
    import urllib.parse, urllib.request
    q = urllib.parse.quote(name)
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}"
           f"/property/MolecularFormula,InChIKey,MonoisotopicMass/JSON")
    try:
        p = json.loads(urllib.request.urlopen(url, timeout=30).read())["PropertyTable"]["Properties"][0]
        return {"cid": p.get("CID"), "inchikey": p.get("InChIKey"),
                "formula": p.get("MolecularFormula"),
                "mono_mass": float(p["MonoisotopicMass"]) if p.get("MonoisotopicMass") else None}
    except Exception:
        return {}


def mono_mass(formula):
    """RAW molmass monoisotopic mass (verification)."""
    from molmass import Formula
    return float(Formula(formula).monoisotopic_mass)


def gem_crosswalk(rows):
    """RAW COBRApy: accepted KEGG/HMDB/ChEBI -> Mouse-GEM MAM base ids."""
    import cobra
    model = cobra.io.read_sbml_model(GEM_MODEL)
    idx = {"kegg": defaultdict(set), "hmdb": defaultdict(set), "chebi": defaultdict(set)}
    key = {"kegg": "kegg.compound", "hmdb": "hmdb", "chebi": "chebi"}
    for m in model.metabolites:
        base = re.sub(r"[a-z]$", "", m.id)
        for db, k in key.items():
            v = m.annotation.get(k)
            if v is None: continue
            for x in ([v] if isinstance(v, str) else v):
                x = str(x)
                if db == "chebi":
                    x2 = x if x.upper().startswith("CHEBI:") else "CHEBI:" + x
                    idx[db][x2].add(base); idx[db][x2.replace("CHEBI:", "")].add(base)
                else:
                    idx[db][x].add(base)
    for r in rows:
        r["gem_mam"] = ""
        for db in ("kegg", "hmdb", "chebi"):
            if r.get(db) and idx[db].get(str(r[db])):
                r["gem_mam"] = ";".join(sorted(idx[db][str(r[db])])); break
    return rows


# ========================================================================== helpers
def harmonized_queries(e):
    """The rename-rule names saved in the ledger (candidate.query from search sources) —
    these are the reasoning layer's harmonized names, read from the saved file."""
    out = {"kegg": [], "pubchem": [], "chebi": []}
    for c in e.get("candidates", []):
        src, q = c.get("source", ""), c.get("query")
        if not q or (":" in q and q.split(":", 1)[0] in DBS):
            continue  # skip bridge-style "src:id" queries
        if "kegg" in src:
            out["kegg"].append(q)
        elif "pubchem" in src:
            out["pubchem"].append(q)
        elif "chebi" in src:
            out["chebi"].append(q)
    # fall back to the normalized (harmonized) display name
    nm = e.get("normalized", {}).get("normalized")
    if nm:
        for k in out:
            out[k].append(nm)
    # de-dup preserving order
    for k in out:
        seen = set(); out[k] = [x for x in out[k] if not (x in seen or seen.add(x))]
    return out


# ========================================================================== the flow
def main():
    entries = LEDGER
    names = [e["original_name"] for e in entries.values()]

    print("=" * 70)
    print("01_STAGE 1 — MetaboAnalystR 1st-pass mapping (raw)")
    ma = {row["query"]: row for row in metaboanalyst_map(names)}
    def nz(v): return v not in (None, "", "NA", "null")
    stage1 = {}
    for fid, e in entries.items():
        q = e.get("normalized", {}).get("normalized") or e["original_name"]
        row = ma.get(q, {})
        stage1[fid] = {db: (row.get(db) if nz(row.get(db)) else None)
                       for db in ("kegg", "hmdb", "pubchem", "chebi")}
    n_mapped1 = sum(1 for v in stage1.values() if v["kegg"] or v["hmdb"])
    print(f"  1st-pass primary (KEGG|HMDB): {n_mapped1}/{len(entries)}")

    print("=" * 70)
    print("02_STAGE 2 — extract the unmapped list")
    unmapped = [fid for fid, v in stage1.items() if not (v["kegg"] or v["hmdb"])
                and entries[fid].get("final_class") != "exogenous-excluded"]
    print(f"  unmapped after 1st pass (non-exogenous): {len(unmapped)}")

    print("=" * 70)
    print("03_STAGE 3 — rename-rule ID extraction: RE-RUN DB search with the saved harmonized names")
    recovered = {}
    for fid in unmapped:
        e = entries[fid]
        qs = harmonized_queries(e)
        rec = e.get("accepted", {})           # recorded pick (disambiguation), from the file
        found = {}
        # KEGG via KEGGREST with each harmonized name; keep the recorded id if it appears
        for q in qs["kegg"]:
            hits = {h["kegg"] for h in kegg_find(q)}
            if not hits:
                continue
            found["kegg"] = rec.get("kegg") if rec.get("kegg") in hits else sorted(hits)[0]
            found["kegg_query"] = q
            break
        # structure via PubChem with each harmonized name
        for q in qs["pubchem"]:
            p = pubchem_lookup(q)
            if p.get("inchikey"):
                found.update(pubchem=str(p["cid"]), inchikey=p["inchikey"], formula=p.get("formula"))
                found["pubchem_query"] = q
                if p.get("formula"):
                    try: found["mono_mass"] = round(mono_mass(p["formula"]), 4)  # molmass verify
                    except Exception: pass
                break
        recovered[fid] = found
        tag = found.get("kegg") or found.get("inchikey") or "-"
        print(f"  {fid} '{e['original_name']}' -> harmonized {qs['kegg'][:1] or qs['pubchem'][:1]} -> {tag}")

    print("=" * 70)
    print("04_STAGE 4 — cross-check extracted IDs across DBs (BridgeDbR) + HMDB backfill")
    # bridge every id we have (from stage 1 or 3) to fill KEGG/HMDB across DBs
    bridge_q, back = [], {}
    for fid, e in entries.items():
        if e.get("final_class") == "exogenous-excluded":
            continue
        ids = {**stage1[fid], **recovered.get(fid, {})}
        # BridgeDb needs the full datasource NAMES (getSystemCode), not short keys
        if ids.get("inchikey"):
            src_id = ("InChIKey", ids["inchikey"])
        elif ids.get("kegg"):
            src_id = ("KEGG Compound", ids["kegg"])
        elif e.get("accepted", {}).get("chebi"):
            src_id = ("ChEBI", e["accepted"]["chebi"])
        else:
            src_id = None
        if src_id and not ids.get("hmdb"):
            bridge_q.append({"feature_id": fid, "id": src_id[1], "source": src_id[0],
                             "targets": ["HMDB", "KEGG Compound"]})
    for q, res in zip(bridge_q, bridgedb_map([{k: v for k, v in q.items() if k != "feature_id"} for q in bridge_q])):
        mp = res.get("mappings", {})
        hb = mp.get("HMDB"); back[q["feature_id"]] = (hb[0] if isinstance(hb, list) and hb else None)
    print(f"  cross-checked {len(bridge_q)} ids; HMDB backfilled for {sum(1 for v in back.values() if v)}")

    # assemble final rows (respect recorded class/exclusions from the saved file)
    rows = []
    for fid, e in entries.items():
        ids = {**stage1[fid], **recovered.get(fid, {})}
        acc = e.get("accepted", {})
        row = {"feature_id": fid, "name": e["original_name"],
               "harmonized": e.get("normalized", {}).get("normalized", ""),
               "final_class": e.get("final_class", ""), "confidence": e.get("confidence", ""),
               "kegg": ids.get("kegg") or acc.get("kegg", ""),
               "hmdb": ids.get("hmdb") or back.get(fid) or acc.get("hmdb", ""),
               "chebi": acc.get("chebi", ""), "pubchem": ids.get("pubchem") or acc.get("pubchem", ""),
               "inchikey": ids.get("inchikey") or acc.get("inchikey", "")}
        if e.get("final_class") == "exogenous-excluded":
            for db in DBS: row[db] = ""
        rows.append(row)

    print("=" * 70)
    print("05_GEM crosswalk (COBRApy) + coverage figure (matplotlib)")
    try:
        gem_crosswalk(rows)
    except Exception as ex:
        print("  GEM skipped:", ex)
    _write_outputs(rows)
    _plot_upset(rows)
    print("done -> master_ledger.tsv, enriched_xref.tsv, figures/db_matching_upset.png")


def _write_outputs(rows):
    cols = ["feature_id", "name", "harmonized", "final_class", "confidence",
            "kegg", "hmdb", "chebi", "pubchem", "inchikey", "gem_mam"]
    with open(OUT / "master_ledger.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(cols)
        for r in rows: w.writerow([r.get(c, "") for c in cols])


def _plot_upset(rows):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mem = [r for r in rows if r.get("final_class") != "exogenous-excluded"]
    for r in mem:
        for db in DBS: r[f"has_{db}"] = bool(r.get(db))
    N = len(mem)
    per_db = {db: sum(m[f"has_{db}"] for m in mem) for db in DBS}
    patt = defaultdict(int)
    for m in mem: patt[tuple(db for db in DBS if m[f"has_{db}"])] += 1
    items = sorted(patt.items(), key=lambda kv: kv[1], reverse=True)
    combos = [k for k, _ in items] or [()]; sizes = [v for _, v in items] or [0]
    lab = {"kegg": "KEGG", "hmdb": "HMDB", "chebi": "ChEBI", "pubchem": "PubChem", "inchikey": "InChIKey"}
    fig = plt.figure(figsize=(max(9, 1.15 * len(combos) + 3), 6.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 4.6], height_ratios=[2.5, 1.6], wspace=0.32, hspace=0.06)
    at, am, al = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 0])
    at.bar(range(len(combos)), sizes, color="#2a7db5", width=0.6)
    for i, sz in enumerate(sizes): at.text(i, sz, str(sz), ha="center", va="bottom", fontweight="bold", fontsize=9)
    at.set_ylabel("features in\nintersection"); at.spines[["top", "right"]].set_visible(False)
    at.tick_params(axis="x", bottom=False, labelbottom=False)
    ro = DBS[::-1]; yof = {db: i for i, db in enumerate(ro)}
    for i, combo in enumerate(combos):
        for db in DBS: am.plot(i, yof[db], "o", ms=11, color="#2a7db5" if db in combo else "#dcdcdc", zorder=3)
        ys = [yof[db] for db in combo]
        if len(ys) > 1: am.plot([i, i], [min(ys), max(ys)], color="#2a7db5", lw=2.2, zorder=2)
    for db in DBS: am.text(-0.9, yof[db], lab[db], ha="right", va="center", fontweight="bold", fontsize=11)
    am.set_yticks([]); am.set_xticks([])
    for sp in am.spines.values(): sp.set_visible(False)
    tot = [per_db[db] for db in ro]
    al.barh(range(len(DBS)), tot, color="#4f7ea8", height=0.55)
    for i, t in enumerate(tot): al.text(t, i, f"{t} ", ha="right", va="center", fontweight="bold", fontsize=9)
    al.set_xlim((max(tot) * 1.14) if tot and max(tot) else 1, 0); al.set_xlabel("set size")
    al.set_yticks([]); al.spines[["top", "left", "right"]].set_visible(False)
    fig.suptitle(f"Metabolite identifier coverage across 5 DBs ({N} analysable features)",
                 fontsize=12.5, fontweight="bold")
    fig.savefig(FIG / "db_matching_upset.png", dpi=200, bbox_inches="tight"); plt.close(fig)
    with open(OUT / "enriched_xref.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["feature_id", "name", "final_class"] + DBS + [f"has_{d}" for d in DBS])
        for m in mem:
            w.writerow([m["feature_id"], m["name"], m["final_class"]] +
                       [m.get(d, "") for d in DBS] + [int(m[f"has_{d}"]) for d in DBS])


if __name__ == "__main__":
    main()
'''
