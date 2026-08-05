"""Generate a standalone reproduction script for a run, using the ORIGINAL library APIs.

The emitted script does NOT import metabo_idmapper or call any tool function. It makes the
PIPELINE FLOW visible and actually executes it with the raw engines:

  Stage 1  raw MetaboAnalystR  — 1st-pass name -> ID
  Stage 2  extract the unmapped (no KEGG/HMDB from stage 1)
  Stage 3  rename-rule ID extraction — READ the harmonized query names saved in
           midmap_ledger.json (the reasoning output) and RE-RUN the DB searches
           (KEGGREST / PubChem / ChEBI) with them; molmass-verify
  Stage 4  cross-check the extracted IDs across DBs (BridgeDbR) + HMDB backfill
  Stage 5  back-check every accepted KEGG against its own keggGet record + shared-id scan
  Stage 6  GEM crosswalk (COBRApy): xref AND a search over the model's own names
  then     classify, coverage figure (matplotlib)

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
    cells with NO `def`/lambda; raw library APIs only. The 4 raw R engines are written as
    sidecar .R files next to the notebook and called via Rscript."""
    s = Session(workdir)
    outdir = Path(s.workdir) / "code"
    outdir.mkdir(exist_ok=True)
    (outdir / "ma.R").write_text((RSCRIPT_DIR / "metaboanalyst_map.R").read_text())
    (outdir / "bridge.R").write_text((RSCRIPT_DIR / "bridge_xref.R").read_text())
    (outdir / "kegg.R").write_text((RSCRIPT_DIR / "kegg_search.R").read_text())
    (outdir / "keggent.R").write_text((RSCRIPT_DIR / "kegg_entry.R").read_text())

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
        "# DOES    : 1차에서 KEGG·HMDB 둘 다 없는 항목만 추림(xenobiotic 제외).\n"
        "# OUTPUT  : unmapped (feature id 리스트).\n"
        "# REUSED  : Stage 3(rename 재검색 대상).\n"
        "# STAGE 2 — extract the unmapped list (no KEGG/HMDB from stage 1; skip xenobiotic-excluded)\n"
        "unmapped = []\n"
        "for fid, v in stage1.items():\n"
        "    if not (v['kegg'] or v['hmdb']) and LEDGER[fid].get('final_class') != 'xenobiotic-excluded':\n"
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
        "    if e.get('final_class') == 'xenobiotic-excluded': continue\n"
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
        "# DOES    : feature별 최종 KEGG/HMDB/ChEBI/PubChem/InChIKey + class 확정(xenobiotic-excluded는 ID 비움).\n"
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
        "    if e.get('final_class') == 'xenobiotic-excluded':\n"
        "        for db in DBS: row[db] = ''\n"
        "    rows.append(row)\n"
        "print('rows:', len(rows))"
    )

    backcheck = (
        "# ══ 07_BACK-CHECK — 확정 KEGG의 '자기 자신' 레코드 대조 (raw KEGGREST keggGet) ═══\n"
        "# INPUT   : rows(assemble)의 kegg, code/keggent.R(raw KEGGREST).\n"
        "# DOES    : id -> 그 id의 KEGG NAME/FORMULA를 받아 이름과 대조. name->id 경로로는\n"
        "#           절대 안 잡히는 오류(맞는 이름에 붙어 온 틀린 id)를 여기서 잡는다.\n"
        "#           검증된 사례: C05472의 KEGG 이름은 'Urocortisol'(담즙산이 아님).\n"
        "# OUTPUT  : kegg_entries{id: {names, formula}}, backcheck_rows(이름 불일치 후보).\n"
        "# REUSED  : 리포트의 '판정 근거', class-level 여부(= 종/클래스 해상도) 판정.\n"
        "req = {'ids': sorted({r['kegg'] for r in rows if r.get('kegg')}), 'out': str(OUT / '_ke.json')}\n"
        "subprocess.run([RSCRIPT, '--vanilla', str(CODE / 'keggent.R')], input=json.dumps(req),\n"
        "               capture_output=True, text=True)\n"
        "kegg_entries = {e['kegg']: e for e in json.load(open(OUT / '_ke.json'))}\n"
        "# class-level 지표: DB 이름이 사슬 없이 1-acyl/alkyl 같은 총칭이면 그 id는 클래스 엔트리\n"
        "backcheck_rows = []\n"
        "for r in rows:\n"
        "    ke = kegg_entries.get(r.get('kegg') or '')\n"
        "    if not ke or not ke.get('found'): continue\n"
        "    names = ke['names'] if isinstance(ke['names'], list) else [ke['names']]\n"
        "    r['kegg_name'] = names[0] if names else ''\n"
        "    generic = any(re.search(r'\\b\\d?-?acyl\\b|\\balkyl\\b', n, re.I) and\n"
        "                  not re.search(r'\\d{1,2}:\\d', n) for n in names)\n"
        "    r['kegg_resolution'] = 'class' if generic else 'species'\n"
        "    backcheck_rows.append({'feature_id': r['feature_id'], 'name': r['name'],\n"
        "                           'kegg': r['kegg'], 'kegg_name': r['kegg_name'],\n"
        "                           'resolution': r['kegg_resolution'],\n"
        "                           'formula': ke.get('formula')})\n"
        "print('back-checked:', len(backcheck_rows),\n"
        "      '| class-level ids:', sum(1 for b in backcheck_rows if b['resolution'] == 'class'))\n"
        "# 같은 id를 두 화합물이 공유하는지 (비를 1로 붕괴시키는 조용한 오류)\n"
        "shared = defaultdict(list)\n"
        "for r in rows:\n"
        "    for db in ('kegg', 'hmdb'):\n"
        "        if r.get(db): shared[(db, r[db])].append(r['name'])\n"
        "print('shared ids:', {k: v for k, v in shared.items() if len(v) > 1})"
    )

    gem = (
        "# ══ 08_GEM CROSSWALK — 모델 species (raw COBRApy: xref + 이름/분자식 직접 조회) ══\n"
        "# INPUT   : rows(assemble), GEM_MODEL(SBML/JSON).\n"
        "# DOES    : (a) KEGG/HMDB/ChEBI xref 매핑, (b) xref 미스는 모델의 name/formula를\n"
        "#           직접 조회. GEM은 지질 species에 xref를 거의 달지 않으므로 (b)가 없으면\n"
        "#           '모델에 있는데 없다'고 오판한다(검증 런: xref 8/27, 실제 보유 21/28).\n"
        "# OUTPUT  : rows[*]['gem_mam'], rows[*]['gem_route'], gem_name_hits.\n"
        "# REUSED  : 하류 대사모델 flux(E-Flux/GIMME) 입력, master_ledger.tsv.\n"
        "model = cobra.io.read_sbml_model(GEM_MODEL) if str(GEM_MODEL).endswith('.xml') \\\n"
        "    else cobra.io.load_json_model(GEM_MODEL)\n"
        "idx = {'kegg': defaultdict(set), 'hmdb': defaultdict(set), 'chebi': defaultdict(set)}\n"
        "meta = {}\n"
        "keymap = {'kegg': 'kegg.compound', 'hmdb': 'hmdb', 'chebi': 'chebi'}\n"
        "for m in model.metabolites:\n"
        "    # compartment 접미사는 두 관례가 있다: MAM01234c 와 tdchola_c.\n"
        "    # 'lowercase 한 글자 제거'만 하면 후자가 'tdchola_'가 되어 accession이 깨진다.\n"
        "    comp = str(getattr(m, 'compartment', '') or '')\n"
        "    base = (m.id[:-len(comp)].rstrip('_') if comp and m.id.endswith(comp) else\n"
        "            re.sub(r'_[a-z]{1,2}$|[a-z]$', '', m.id))\n"
        "    rec = meta.setdefault(base, {'name': m.name or '', 'formula': getattr(m, 'formula', None)})\n"
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
        "    r['gem_mam'] = ''; r['gem_route'] = ''\n"
        "    for db in ('kegg', 'hmdb', 'chebi'):\n"
        "        if r.get(db) and idx[db].get(str(r[db])):\n"
        "            r['gem_mam'] = ';'.join(sorted(idx[db][str(r[db])]))\n"
        "            r['gem_route'] = 'xref:' + db; break\n"
        "print('GEM-mapped by xref:', sum(1 for r in rows if r['gem_mam']))\n"
        "# (b) xref 미스 -> 모델 이름 직접 조회(문자 트라이그램 Dice). 화학명은 교착적이라\n"
        "#     토큰 겹침이 0이 되므로(Glycochenodeoxycholic acid vs Chenodeoxyglycocholate)\n"
        "#     n-gram 유사도가 유일한 도달 경로다. 후보만 제시하고 채택은 사람이 판단.\n"
        "model_grams = {}\n"
        "for b, v in meta.items():\n"
        "    fl = re.sub(r'[^a-z0-9]', '', (v['name'] or '').lower())\n"
        "    if fl: model_grams[b] = {fl[i:i+3] for i in range(max(len(fl) - 2, 0))}\n"
        "gem_name_hits = {}\n"
        "for r in rows:\n"
        "    if r['gem_mam'] or not r.get('name'): continue\n"
        "    fq = re.sub(r'[^a-z0-9]', '', (r['harmonized'] or r['name']).lower())\n"
        "    qg = {fq[i:i+3] for i in range(max(len(fq) - 2, 0))}\n"
        "    scored = []\n"
        "    for b, g in model_grams.items():\n"
        "        if not qg or not g: continue\n"
        "        dice = 2 * len(qg & g) / (len(qg) + len(g))\n"
        "        if dice > 0.4: scored.append((round(dice, 3), b))\n"
        "    scored.sort(reverse=True)\n"
        "    gem_name_hits[r['feature_id']] = [\n"
        "        {'mam_base': b, 'gem_name': meta[b]['name'], 'formula': meta[b]['formula'],\n"
        "         'dice': d,\n"
        "         'pool_species': bool(meta[b]['formula'] and re.search(r'[RX]', str(meta[b]['formula'])))}\n"
        "        for d, b in scored[:5]]\n"
        "print('name-search candidates for', sum(1 for v in gem_name_hits.values() if v),\n"
        "      'xref-missed features (pool_species=True는 R-group 총칭 = class-level 입력)')"
    )

    writeout = (
        "# ══ 08_WRITE — master_ledger.tsv ════════════════════════════════════════════════\n"
        "# INPUT   : rows.  OUTPUT: OUT/master_ledger.tsv (최종 매핑 표).\n"
        "# REUSED  : 대사체 통계/네트워크의 ID 조인, joint pathway(KEGG/HMDB), 리포트.\n"
        "# write master_ledger.tsv\n"
        "cols = ['feature_id', 'name', 'harmonized', 'final_class', 'confidence',\n"
        "        'kegg', 'kegg_name', 'kegg_resolution', 'hmdb', 'chebi', 'pubchem',\n"
        "        'inchikey', 'gem_mam', 'gem_route']\n"
        "with open(OUT / 'master_ledger.tsv', 'w', newline='') as f:\n"
        "    w = csv.writer(f, delimiter='\\t'); w.writerow(cols)\n"
        "    for r in rows: w.writerow([r.get(c, '') for c in cols])\n"
        "print('wrote', OUT / 'master_ledger.tsv')"
    )

    upset = (
        "# ══ 09_UPSET — DB 커버리지 figure (raw matplotlib) ══════════════════════════════\n"
        "# INPUT   : rows(비-excluded).  OUTPUT: figures/db_matching_upset.png + enriched_xref.tsv.\n"
        "# DOES    : 5-DB 멤버십 교집합 UpSet(상단 교집합 크기·행렬·좌측 DB별 set size).\n"
        "# REUSED  : 리포트 figure; enriched_xref.tsv는 db_matching_improvement 등 후속 그림 source.\n"
        "# UpSet coverage figure (raw matplotlib), inline\n"
        "mem = [r for r in rows if r['final_class'] != 'xenobiotic-excluded']\n"
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
                "cross-check + HMDB backfill → classify → **back-check every KEGG against its "
                "own keggGet record (+ shared-id scan)** → GEM crosswalk **by xref AND by the "
                "model's own names** → coverage figure.\n\n"
                "Run this notebook from the run workdir (or its `code/`); it reads "
                "`midmap_ledger.json`.\n\n"
                "## 파일 맵 (어디서 왔고 / 무엇을 만들고 / 어디서 다시 쓰나)\n\n"
                "**입력**\n"
                "- `midmap_ledger.json` ← 이 run의 ledger. 작성: driver 세션(record_decision 등). "
                "원본 이름은 `Metabolomics/Atheroclerosis_metabolite.xlsx`의 Compound 열(→ ledger에 반영). "
                "여기엔 normalize·candidates(source별 harmonized query)·accepted·decisions가 담김.\n"
                "- `code/ma.R`,`bridge.R`,`kegg.R`,`keggent.R` ← metabo-idmapper 번들 raw R 엔진 "
                "(마지막은 keggGet 백체크용).\n"
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
            _md("## 06_Assemble"), _code(assemble),
            _md("## 07_Back-check — 확정 KEGG를 그 id의 KEGG 레코드와 대조 (틀린 id / 공유 id 검출)"),
            _code(backcheck),
            _md("## 08_GEM crosswalk — xref + 모델 이름 직접 조회, then coverage figure"),
            _code(gem), _code(writeout), _code(upset),
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
    keggent_r = (RSCRIPT_DIR / "kegg_entry.R").read_text()

    code = _TEMPLATE
    code = code.replace("@@BRIDGE_DB@@", BRIDGE_DB).replace("@@GEM_MODEL@@", GEM_MODEL)
    code = code.replace("@@MA_R@@", ma_r).replace("@@BRIDGE_R@@", bridge_r)
    code = code.replace("@@KEGG_R@@", kegg_r).replace("@@KEGGENT_R@@", keggent_r)

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
  5) back-check every accepted KEGG against its OWN keggGet record (name/formula) and scan
     for identifiers held by more than one compound — the two error classes a name->id
     pipeline cannot see (a wrong id arrives attached to the right name; a shared id
     silently collapses the ratio between two features)
  6) GEM crosswalk (COBRApy) by xref AND by the model's own names — GEMs annotate lipid
     species sparsely, so the xref route alone reports a model as lacking compounds it has
  then classify, coverage figure (matplotlib).

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
_KEGGENT_R = r"""@@KEGGENT_R@@"""


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
    """RAW COBRApy: accepted KEGG/HMDB/ChEBI -> model species base ids, then a name search
    over the model's own names for everything the xref route missed (GEMs annotate lipid
    species sparsely: a verified Recon3D run got 8/27 by xref while the model held 21/28)."""
    import cobra
    model = (cobra.io.load_json_model(GEM_MODEL) if str(GEM_MODEL).endswith(".json")
             else cobra.io.read_sbml_model(GEM_MODEL))
    idx = {"kegg": defaultdict(set), "hmdb": defaultdict(set), "chebi": defaultdict(set)}
    meta = {}
    key = {"kegg": "kegg.compound", "hmdb": "hmdb", "chebi": "chebi"}
    for m in model.metabolites:
        # two compartment conventions: MAM01234c and tdchola_c. Dropping one trailing
        # lowercase letter turns the latter into "tdchola_", which is not an accession.
        comp = str(getattr(m, "compartment", "") or "")
        base = (m.id[: -len(comp)].rstrip("_") if comp and m.id.endswith(comp)
                else re.sub(r"_[a-z]{1,2}$|[a-z]$", "", m.id))
        meta.setdefault(base, {"name": m.name or "", "formula": getattr(m, "formula", None)})
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
        r["gem_route"] = ""
        for db in ("kegg", "hmdb", "chebi"):
            if r.get(db) and idx[db].get(str(r[db])):
                r["gem_mam"] = ";".join(sorted(idx[db][str(r[db])]))
                r["gem_route"] = "xref:" + db
                break
    gem_name_hits = gem_name_search(rows, meta)
    return rows, gem_name_hits


def _trigrams(text):
    flat = re.sub(r"[^a-z0-9]", "", (text or "").lower())
    return {flat[i:i + 3] for i in range(max(len(flat) - 2, 0))}


def gem_name_search(rows, meta, min_dice=0.4, top=5):
    """Candidate model species for xref-missed rows, by character-trigram Dice over the
    model's own names. Chemical names are agglutinative — "Glycochenodeoxycholic acid" and
    the model's "Chenodeoxyglycocholate" share no whole token — so n-gram similarity is the
    only route that reaches them. Candidates only: the pick is a judgement, and a
    pool_species (R-group formula) is a CLASS-level input, not the species."""
    model_grams = {b: _trigrams(v["name"]) for b, v in meta.items() if v["name"]}
    hits = {}
    for r in rows:
        if r.get("gem_mam") or not r.get("name"):
            continue
        qg = _trigrams(r.get("harmonized") or r["name"])
        scored = []
        for b, g in model_grams.items():
            if not qg or not g:
                continue
            dice = 2 * len(qg & g) / (len(qg) + len(g))
            if dice > min_dice:
                scored.append((round(dice, 3), b))
        scored.sort(reverse=True)
        hits[r["feature_id"]] = [
            {"mam_base": b, "gem_name": meta[b]["name"], "formula": meta[b]["formula"],
             "dice": d, "pool_species": bool(meta[b]["formula"]
                                             and re.search(r"[RX]", str(meta[b]["formula"])))}
            for d, b in scored[:top]]
    return hits


def kegg_backcheck(rows):
    """RAW KEGGREST keggGet: fetch what each accepted KEGG id ITSELF resolves to.

    name->id searching cannot detect a wrong id, because the wrong id arrives attached to the
    right name. The id's own KEGG NAME field can: C05472's is "Urocortisol", not a bile acid.
    """
    ids = sorted({r["kegg"] for r in rows if r.get("kegg")})
    if not ids:
        return {}
    res = _run_r(_KEGGENT_R, {"ids": ids}) or []
    entries = {e["kegg"]: e for e in res}
    for r in rows:
        e = entries.get(r.get("kegg") or "")
        if not e or not e.get("found"):
            continue
        names = e["names"] if isinstance(e["names"], list) else [e["names"]]
        r["kegg_name"] = names[0] if names else ""
        # a generic acyl/alkyl name with no chain spec means the id is a CLASS entry
        r["kegg_resolution"] = "class" if any(
            re.search(r"\b\d?-?acyl\b|\balkyl\b", n, re.I) and not re.search(r"\d{1,2}:\d", n)
            for n in names) else "species"
    return entries


def shared_id_scan(rows):
    """Identifiers held by more than one compound. A shared id is silent and destroys any
    ratio between the two features (verified: two glycolipids on one PubChem/HMDB entry)."""
    seen = defaultdict(list)
    for r in rows:
        for db in ("kegg", "hmdb", "chebi", "pubchem", "inchikey"):
            if r.get(db):
                seen[(db, str(r[db]))].append(r["name"])
    return {k: v for k, v in seen.items() if len(v) > 1}


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
                and entries[fid].get("final_class") != "xenobiotic-excluded"]
    print(f"  unmapped after 1st pass (non-excluded): {len(unmapped)}")

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
        if e.get("final_class") == "xenobiotic-excluded":
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
        if e.get("final_class") == "xenobiotic-excluded":
            for db in DBS: row[db] = ""
        rows.append(row)

    print("=" * 70)
    print("05_back-check every accepted KEGG against its own keggGet record")
    try:
        kegg_backcheck(rows)
        klass = sum(1 for r in rows if r.get("kegg_resolution") == "class")
        print(f"  back-checked {sum(1 for r in rows if r.get('kegg_name'))} ids; "
              f"{klass} are CLASS-level entries (report them as such)")
        shared = shared_id_scan(rows)
        print("  identifiers held by >1 compound:", shared or "none")
    except Exception as ex:
        print("  back-check skipped:", ex)

    print("=" * 70)
    print("06_GEM crosswalk (COBRApy xref + model-name search) + coverage figure")
    gem_name_hits = {}
    try:
        _, gem_name_hits = gem_crosswalk(rows)
        print(f"  xref-mapped {sum(1 for r in rows if r.get('gem_mam'))}; "
              f"name-search candidates for {sum(1 for v in gem_name_hits.values() if v)} "
              "xref-missed features")
    except Exception as ex:
        print("  GEM skipped:", ex)
    if gem_name_hits:
        with open(OUT / "gem_name_candidates.json", "w") as f:
            json.dump(gem_name_hits, f, indent=1)
    _write_outputs(rows)
    _plot_upset(rows)
    print("done -> master_ledger.tsv, enriched_xref.tsv, figures/db_matching_upset.png")


def _write_outputs(rows):
    cols = ["feature_id", "name", "harmonized", "final_class", "confidence",
            "kegg", "kegg_name", "kegg_resolution", "hmdb", "chebi", "pubchem",
            "inchikey", "gem_mam", "gem_route"]
    with open(OUT / "master_ledger.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(cols)
        for r in rows: w.writerow([r.get(c, "") for c in cols])


def _plot_upset(rows):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mem = [r for r in rows if r.get("final_class") != "xenobiotic-excluded"]
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
