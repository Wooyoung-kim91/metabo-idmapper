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
| normalize names; exact DB match; PubChem/KEGG/ChEBI search; BridgeDb xref; molmass formula/mass verify; m/z→mass windows; **ID→own-DB-record back-check**; **isomer/class token comparison**; **shared-ID detection**; GEM xref crosswalk + **model name/formula/mass search**; coverage | is a fuzzy/typo/synonym candidate correct? abbreviation expansion? endogenous vs **xenobiotic**? which isomer, when the names disagree? is a class-level ID acceptable here? is a model species the compound, a class proxy, a surrogate, or genuinely absent? confidence tier; final inclusion |

## Tools (27)

`midmap_guidance` · `detect_state` · `ingest_names` · `exact_match` · `structure_lookup` ·
`search_synonym` · `bridge_xref` · `verify_candidate` · `mass_match_candidates` ·
**`id_name_check`** · **`isomer_guard`** · **`collision_check`** · **`acknowledge_flag`** ·
`screen_exogenous` · `record_decision` · `backfill_hmdb` · `gem_crosswalk` ·
**`gem_search`** · **`gem_assign`** ·
`mapping_provenance` · `annotate_source` · `plot_coverage` · `export_code` ·
`export_report_ppt` · `coverage_summary` · **`finalize_run`** · `harness_audit`

Call `midmap_guidance` first for the canonical workflow, confidence tiers (M1–U), and gotchas.

### The wrong-ID problem (why the three verification tools exist)

A missing ID is visible. A **wrong** ID is not — it arrives attached to the right name. Every
serious error in verified runs was of that kind, and none of them is catchable by comparing
IDs, or by formula/mass, because the alternatives are formula-identical isomers:

| what happened | what catches it |
|---|---|
| a taurochenodeoxycholate entry carrying KEGG `C05472`, whose own KEGG name is **"Urocortisol"** (a cortisol metabolite) — a wrong BridgeDb link | `id_name_check`: ID → its own DB record → compare names |
| two different glycolipids (Lc3Cer, nLc4Cer) sharing one PubChem/HMDB entry, which **collapses the ratio between them to exactly 1** | `collision_check` |
| a plasmenyl (vinyl-ether) `LysoPC(P-16:0)` carrying the **1-acyl class** ID | `isomer_guard` / `id_name_check` (`ether_linkage`) |
| a **neolacto** (β1-4) glycan carrying the **lacto** (β1-3) isomer | `isomer_guard` (`glycan_series`, from KEGG's own linkage wording) |

`isomer_guard` compares two names on the axes that actually decide lipid and glycolipid
identity — skeleton · ether linkage (P- plasmenyl / O- plasmanyl / 1-acyl) · glycan series and
length · glycosidic linkages · sialyl linkage (α2,3 vs α2,6) · chains · double-bond positions ·
omega · oxidation vs **per**oxidation (HETE vs HPETE) · hydroxy count · acetyl count · polyamine
backbone — and reports whether a name is **species**- or **class**-level. It is pure string
logic: no network, no lookups. `id_name_check` runs it against the record an ID resolves to;
`harness_audit` **fails** a run where a searched or bridged KEGG was never back-checked, or
where an identifier stands for two different compounds.

### Genome-scale model side: search the model, then record HOW it maps

`gem_crosswalk` alone under-reports badly. GEMs annotate lipid species sparsely: a verified
Recon3D run crosswalked **8 of 27** compounds by xref while the model actually contained **21 of
28** — under names no identifier reaches (`dgchol` "Chenodeoxyglycocholate" for GCDCA,
"Sialyl-3-paragloboside", "Lactoneotetraosylceramide", and glycolipids written as sugar
compositions `(Gal)1 (Glc)1 (GlcNAc)1 (Cer)1`).

So every xref miss now returns `name_suggestions` from the model's **own** names, searched with
the entry name plus its deterministic lipid-shorthand variants (`LysoPC(16:0)` →
`1-palmitoyl-sn-glycero-3-phosphocholine`, `nLc4Cer` → `lactoneotetraosylceramide` /
`paragloboside`, `Glycochenodeoxycholic acid` → `chenodeoxyglycocholate`), each with an isomer
verdict. `gem_search` searches name / id text / formula / mass directly and returns
`same_formula_groups` — formula-identical isomer candidates that mass **cannot** separate
(`pcholar_hs` arachidonoyl Δ5,8,11,14 vs `pcholn204_hs` Δ8,11,14,17 for LPC 20:4).

`gem_assign` then commits the call with an explicit **relation**, per model label, protected
from being overwritten by a later crosswalk, and refusing any species not in the model:

| relation | meaning for a flux result |
|---|---|
| `exact` | the species IS this compound — species-level input |
| `class-proxy` | only the generic R-group pool species exists (`crm_hs`, `sphmyln_hs`): usable, but chain length/linkage is not represented — report as class-level |
| `isomer-surrogate` | a **different** species stands in (N1-acetylsperm**ine** → N1-acetylsperm**idine**): not an identity, the substitution travels with every number |
| `model-scope-absent` | genuinely not in the model — a **result**, so the untestable hypotheses are known |
| `id-gap` | should be there, nothing has resolved it yet |

Without this axis a curated mapping cannot be stored at all: in the verified run 13 hand-found
species lived only in a report while the ledger still said `id-gap`. `harness_audit` now warns
on any model-relevant entry with no recorded relation, and the run emits `gem_curation.tsv`.

### Origin taxonomy — exogenous vs xenobiotic
`record_decision` resolves two axes: **ID coverage** (`final_class`) and **provenance**
(`origin`). Non-endogenous compounds are split into two distinct classes instead of one
"excluded" bucket:

| final_class | origin | disposition |
|---|---|---|
| `KEGG-mapped` / `HMDB-mapped` / `structure-only` | `endogenous` (default) | host-produced, analysed |
| **`exogenous`** | `diet` · `drug` · `microbial` · `plant` | **KEPT + tagged** — a real outside-host signal; keeps its KEGG/HMDB IDs and can enter the GEM crosswalk |
| **`xenobiotic-excluded`** | `contaminant` · `industrial` · `additive` · `surfactant` · `plasticizer` · `reagent` | **EXCLUDED** — non-biological (LC-MS additive / surfactant / plasticizer / industrial) |

`screen_exogenous` detects only the **non-biological** classes deterministically (and
auto-suggests the origin); the biological-exogenous call (drug / diet / gut-microbial / plant)
is a reasoning-layer judgement. A drug or dietary metabolite is tagged `exogenous` and **kept** —
never dumped into the excluded bucket. The `harness_audit` `origin_coherence` check fails any
entry whose origin and class disagree (e.g. `origin='drug'` with `xenobiotic-excluded`).

### Ledger contract and review flags

The ledger (`midmap_ledger.json`) is the run's single source of truth and every tool call
**replaces the whole file**. So a call that finds the file changed since it read it refuses to
write and returns `error_code: "ledger_conflict"` with `wrote_nothing: true` — two interleaved
calls on one workdir would otherwise silently drop whichever decisions were recorded first.
The file carries a `schema_version`; an older ledger is migrated when loaded (persisted by the
next write) and the migration is reported, including **model accessions repaired** from an
earlier base-id bug (`tdchola_` → `tdchola`).

Review flags are **derived, not stored as stale strings**. Each flag has an append-only record
of when and why it was raised, plus a predicate that decides whether it is still true — so
fixing the underlying problem closes the flag by itself and the audit reads *state* rather than
replaying decision history. `detect_state.flags` splits them three ways:

| | meaning |
|---|---|
| `open` | still true right now — each carries what resolves it |
| `self_resolved` | raised earlier, no longer the case (nothing to do) |
| `acknowledged` | real but not fixable, deliberately accepted via `acknowledge_flag` **with a recorded reason** — reported, never silently dropped |

`action` flags (`id_name_conflict`, `shared_id_collision`, `auto_accept_review`,
`isomer_token_conflict`, `hmdb_backfill_conflict`, `id_gap_try_generic_kegg`) must all be
resolved or acknowledged before a run is finished. `info` flags (`class_level_id`,
`possible_xenobiotic`) are properties to carry into the report — e.g. the nine class-level KEGG
ids in the verified run.

### Governance harness (`harness_audit`)
A **read-only** self-audit — the metabo-idmapper counterpart to a run harness — that makes no
identity call and changes nothing. Run it **last** (after `finalize_run`): it reads the
ledger + emitted artifacts and checks that the reasoning layer actually *honored this project's
own contract*, emitting a per-check **pass / warn / fail** scorecard. It catches rules that were
"defined but not followed": a fabricated ID (accepted but produced by no tool), a mass-only
**W**-tier candidate used as primary, an incoherent `final_class`↔confidence pair, a fuzzy
(M2/M3) accept with no recorded formula/mass verification, a locant/anomer-sensitive name
accepted via a non-exact route and never re-checked, **a searched/bridged KEGG never
back-checked against its own DB record** (`fail`), **an ID that contradicts its own name on a
discriminant axis and is still accepted** (`fail`), **an identifier shared by two different
compounds** (`fail`; a shared *class-level* ID is a `warn` instead, since the DB has nothing
finer), an unresolved HMDB-backfill collision, a xenobiotic class excluded inconsistently, an
origin↔class mismatch, a flagged trade-name auto-accept never reviewed, an id-gap KEGG never
re-tried, **a model-relevant entry with no GEM relation recorded**, **a surrogate/proxy species
used without documenting what it changes**, a decision with no rationale, entries left pending,
a skipped `gem_crosswalk`, or missing stage-7 always-emit artifacts. Fix every `fail`; review
each `warn`.

### Report outputs (`finalize_run` emits all of these in order)
- `master_ledger.tsv`, `coverage_summary.tsv`, `mapping_provenance.tsv`
- **Provenance tables** (`mapping_provenance`): `kegg_recovered.tsv` / `hmdb_recovered.tsv` —
  what was mapped BEYOND the MetaboAnalyst 1st pass, with the harmonized name, id, and the
  logic/route (typo fix, synonym search, xref bridge); `unmapped_harmonization.tsv` —
  structure-only entries with the names tried and why mapping failed; `exogenous_kept.tsv` —
  biological outside-host metabolites kept + tagged; `xenobiotic_excluded.tsv` — non-biological
  contaminants excluded, each with its `origin` + full reason; **`gem_curation.tsv`** — per
  entry the model, species, **relation** (exact / class-proxy / isomer-surrogate /
  model-scope-absent), whether it was curated, and the rationale.
- **PPTX report** (`export_report_ppt`): a slide deck built from the run artifacts —
  Title · Coverage KPIs · Methods · Pipeline · UpSet · Improvement · Recovery cause→fix ·
  KEGG/HMDB recovered · Unmapped · Exogenous(kept) · Xenobiotic(excluded) · Outputs.
- **Annotated source** (`annotate_source`): the ORIGINAL data file with the final ID columns
  appended (`<source>_annotated.xlsx/.tsv`: intensity matrix + ID_kegg / ID_hmdb / ID_chebi /
  ID_pubchem / ID_inchikey / ID_final_class / ID_origin / ID_gem_mam / ID_gem_relation per
  compound). Auto-detects the name column; pass `header=` for vendor sheets with a preamble.
- **Figures** (`plot_coverage`): `figures/db_matching_upset.png` (5-DB coverage UpSet +
  `enriched_xref.tsv`) and `figures/db_matching_improvement.png` (MetaboAnalyst baseline
  vs current logic).
- **Reproduction code** (`export_code`): `code/reproduce_mapping.py` (+ `.ipynb`) — standalone
  reproductions using the ORIGINAL library APIs (MetaboAnalystR / BridgeDbR / KEGGREST via
  Rscript, PubChem PUG-REST, molmass, COBRApy, matplotlib), NOT the tool wrappers. They make
  the flow explicit and RUN it: 1 MetaboAnalyst 1st pass → 2 extract unmapped → 3 re-run
  KEGG/PubChem searches with the harmonized names READ from the saved ledger → 4 BridgeDb
  cross-check + HMDB backfill → **5 back-check every accepted KEGG against its own keggGet
  record + shared-id scan** → **6 model crosswalk by xref AND by the model's own names** →
  7 master table + coverage figure.

  Both artifacts are generated from **one** implementation, `codegen/raw_engine.py` — a real,
  unit-tested module with no dependency on this package. The script inlines it verbatim; the
  notebook ships it as a sidecar beside the four `code/*.R` engines and keeps its cells linear
  (no `def`) with per-cell input / output / reuse comments. That single source exists because
  the logic used to be written out once per emitter, and one model base-id bug then had to be
  fixed in three places; a test now asserts the emitted artifacts copy the engine rather than
  restate it.

## Reasoning layer (embedded in the MCP)

The LLM reasoning layer is not an external subagent — it ships **inside** the MCP as two
prompts and two resources, placed at the stages where judgment actually lives:

| role | MCP prompt | resource | placed at |
|---|---|---|---|
| **driver** — drive the tools end-to-end, own the identity CALLs | `map_metabolites` | `metabo-idmapper://driver` | all stages (judgment in 3 & 5) |
| **reviewer** — adversarially verify accepted identities + exclusions | `review_mappings` | `metabo-idmapper://reviewer` | after `record_decision`, before `finalize_run` |

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
