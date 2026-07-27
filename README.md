# metabo-idmapper

Turn messy **metabolite names** into **input-usable IDs** — KEGG / HMDB / ChEBI /
PubChem / InChIKey, plus a **Mouse-GEM MAM crosswalk** for metabolic-model (flux) input —
exposed as a single **MCP** tool registry.

Same philosophy as scpilot: a **deterministic registry** where every tool *emits evidence*
(candidate IDs, xref bridges, formula/mass verification, coverage). The connecting **LLM is
the reasoning layer** that makes the identity *call* (which candidate is correct, endogenous
vs xenobiotic, confidence tier, final inclusion). Tools never fabricate an ID; `record_decision`
refuses any accepted ID that no tool produced.

## Code vs LLM-judgment split

| Deterministic **code** (MCP tools) | **LLM** judgment (the call) |
|---|---|
| normalize names; exact DB match; PubChem/KEGG/ChEBI search; BridgeDb xref; molmass formula/mass verify; m/z→mass windows; Mouse-GEM crosswalk; coverage | is a fuzzy/typo/synonym candidate correct? abbreviation expansion? endogenous vs **xenobiotic**? id-gap vs model-scope-absent? confidence tier; isomer disambiguation; final inclusion |

## Tools (20)

`midmap_guidance` · `detect_state` · `ingest_names` · `exact_match` · `structure_lookup` ·
`search_synonym` · `bridge_xref` · `verify_candidate` · `mass_match_candidates` ·
`screen_exogenous` · `record_decision` · `backfill_hmdb` · `gem_crosswalk` ·
`mapping_provenance` · `annotate_source` · `plot_coverage` · `export_code` ·
`export_report_ppt` · `coverage_summary` · `harness_audit`

Call `midmap_guidance` first for the canonical workflow, confidence tiers (M1–U), and gotchas.

### Governance harness (`harness_audit`)
A **read-only** self-audit — the metabo-idmapper counterpart to a run harness — that makes no
identity call and changes nothing. Run it **last** (after `coverage_summary`): it reads the
ledger + emitted artifacts and checks that the reasoning layer actually *honored this project's
own contract*, emitting a per-check **pass / warn / fail** scorecard. It catches rules that were
"defined but not followed": a fabricated ID (accepted but produced by no tool), a mass-only
**W**-tier candidate used as primary, an incoherent `final_class`↔confidence pair, a fuzzy
(M2/M3) accept with no recorded formula/mass verification, a locant/anomer-sensitive name
accepted via a non-exact route and never re-checked, a contaminant class excluded
inconsistently, a flagged trade-name auto-accept never reviewed, an id-gap KEGG never re-tried,
a decision with no rationale, entries left pending, a skipped `gem_crosswalk`, or missing
stage-7 always-emit artifacts. Fix every `fail`; review each `warn`.

### Report outputs (always emitted by `coverage_summary`)
- `master_ledger.tsv`, `coverage_summary.tsv`, `mapping_provenance.tsv`
- **Provenance tables** (`mapping_provenance`): `kegg_recovered.tsv` / `hmdb_recovered.tsv` —
  what was mapped BEYOND the MetaboAnalyst 1st pass, with the harmonized name, id, and the
  logic/route (typo fix, synonym search, xref bridge); `unmapped_harmonization.tsv` —
  structure-only entries with the names tried and why mapping failed; `exogenous_excluded.tsv`
  — non-metabolite contaminants/drugs excluded, with category + full reason.
- **PPTX report** (`export_report_ppt`): a slide deck built from the run artifacts —
  Title · Coverage KPIs · Methods · Pipeline · UpSet · Improvement · Recovery cause→fix ·
  KEGG/HMDB recovered · Unmapped · Exogenous · Outputs (full text, paginated).
- **Annotated source** (`annotate_source`): the ORIGINAL data file with the final ID columns
  appended (`<source>_annotated.xlsx/.tsv`: intensity matrix + ID_kegg / ID_hmdb / ID_chebi /
  ID_pubchem / ID_inchikey / ID_final_class / ID_gem_mam per compound). Auto-detects the name
  column; pass `header=` for vendor sheets with a preamble.
- **Figures** (`plot_coverage`): `figures/db_matching_upset.png` (5-DB coverage UpSet +
  `enriched_xref.tsv`) and `figures/db_matching_improvement.png` (MetaboAnalyst baseline
  vs current logic).
- **Reproduction code** (`export_code`): `code/reproduce_mapping.py` (+ `.ipynb`) — standalone
  reproductions using the ORIGINAL library APIs (MetaboAnalystR / BridgeDbR / KEGGREST via
  Rscript, PubChem PUG-REST, molmass, COBRApy, matplotlib), NOT the tool wrappers. They make
  the flow explicit and RUN it: Stage 1 MetaboAnalyst 1st pass → Stage 2 extract unmapped →
  Stage 3 re-run KEGG/PubChem searches with the harmonized names READ from the saved ledger →
  Stage 4 BridgeDb cross-check + HMDB backfill → GEM crosswalk → figure. The `.ipynb` is
  unrolled into linear cells (no `def`) with per-cell input-source / output / reuse comments;
  the 3 raw R engines ship as `code/{ma,bridge,kegg}.R` sidecars.

## Reasoning layer (embedded in the MCP)

The LLM reasoning layer is not an external subagent — it ships **inside** the MCP as two
prompts and two resources, placed at the stages where judgment actually lives:

| role | MCP prompt | resource | placed at |
|---|---|---|---|
| **driver** — drive the tools end-to-end, own the identity CALLs | `map_metabolites` | `metabo-idmapper://driver` | all stages (judgment in 3 & 5) |
| **reviewer** — adversarially verify accepted identities + exclusions | `review_mappings` | `metabo-idmapper://reviewer` | after `record_decision`, before `coverage_summary` |

Connect the server and invoke the `map_metabolites` prompt to become the driver; it invokes
`review_mappings` before finalizing. No config symlinks — the whole reasoning layer travels
with the server.

## Engine

Python (PubChem REST, ChEBI OLS4, molmass, COBRApy) + out-of-process **system R**
(`BridgeDbR`, `KEGGREST`, `MetaboAnalystR`) via bundled `rscripts/`. Reuses the verified
`metabolite-id-harmonization` skill logic.

## Install

```bash
# scientific stack lives in the conda env `cobragem`; don't let pip re-resolve it
/home/wykim/miniforge3/envs/cobragem/bin/pip install -e . --no-deps
```

Required external assets (overridable by env var):

| env var | default |
|---|---|
| `METABO_IDMAP_BRIDGE_DB` | `.../Omics/models/bridgedb/metabolites_20210109.bridge` |
| `METABO_IDMAP_GEM` | `.../Omics/models/Mouse-GEM/Mouse-GEM.xml` |
| `METABO_IDMAP_RSCRIPT` | `Rscript` |

## Run

```bash
python -m metabo_idmapper          # stdio MCP server
```

Register with Claude Code:

```bash
claude mcp add metabo-idmapper -- /home/wykim/miniforge3/envs/cobragem/bin/python -m metabo_idmapper
```

## Test

```bash
/home/wykim/miniforge3/envs/cobragem/bin/python -m pytest -q   # 7 offline smoke tests
```

Live end-to-end control (Taurine): `structure_lookup` → PubChem CID 1123 / `C2H7NO3S`;
`bridge_xref` InChIKey → KEGG `C00245`, HMDB, ChEBI; `gem_crosswalk` → Mouse-GEM MAM.
