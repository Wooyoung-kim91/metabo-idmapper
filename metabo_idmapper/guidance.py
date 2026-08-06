"""Canonical workflow + operating rules, returned by the `midmap_guidance` tool.

This is the single source of truth the reasoning LLM reads at the start of a run
(mirrors scpilot_guidance): what the tools do, what the LLM must decide, and the
confidence-tier vocabulary.
"""

CANONICAL = r"""
# metabo-idmapper — canonical workflow

GOAL: turn a list of messy metabolite NAMES into INPUT-USABLE IDs
(KEGG / HMDB / ChEBI / PubChem / InChIKey) plus a genome-scale-model species crosswalk
(Mouse-GEM MAM, Recon3D, or any COBRA model you point it at) for flux input.

## The failure mode this engine is built against
A missing id is visible; a WRONG id is not. In verified runs every serious error was an id
that looked right: a taurochenodeoxycholate entry carrying KEGG C05472, whose own KEGG name
is "Urocortisol" (a cortisol metabolite); two different glycolipids sharing one PubChem/HMDB
entry, which collapses the ratio between them to exactly 1; a plasmenyl (vinyl-ether) lyso-PC
carrying the 1-ACYL class id; a neolacto (beta1-4) glycan carrying the lacto (beta1-3) isomer.
None of those are catchable by comparing IDs, and none by formula/mass — the isomers are
formula-identical. Three tools exist specifically for them, and they are NOT optional:
  id_name_check   id -> its OWN DB record (name/formula), compared on the isomer axes
  isomer_guard    name vs candidate name on the discriminant axes (offline, no lookups)
  collision_check one identifier standing for two different compounds

REASONING LAYER: the driver + reviewer roles ship INSIDE this MCP as prompts
(`map_metabolites`, `review_mappings`) and resources (`metabo-idmapper://driver`,
`metabo-idmapper://reviewer`). Invoke the `map_metabolites` prompt to become the driver;
before finalizing, invoke `review_mappings` for independent adversarial verification.

## Operating contract (summary-in -> decision-out)
- Tools EMIT evidence: candidate IDs, xref bridges, formula/mass checks, coverage.
- YOU (the LLM) make every identity CALL: which candidate is correct, whether a
  compound is endogenous vs xenobiotic, the confidence tier, and final inclusion.
- NEVER fabricate an ID a tool did not return. NEVER auto-substitute an isomer
  (e.g. Azetidine-1- vs 2-carboxylic acid): if the name says "1", do not accept "2".
- One judgement at a time; record it with `record_decision` so the ledger stays the
  authoritative provenance trail.

## Ledger + review flags (how state behaves)
- One ledger per workdir, and each tool call REPLACES it. A call that finds the file changed
  underneath it returns `error_code: "ledger_conflict"` and writes NOTHING — re-read the state
  and repeat that call. Do not run two writing tools on one workdir at the same time.
- The ledger carries a `schema_version`. An older ledger is migrated when it is loaded (the
  change lands in the file on the next write) and the migration is reported by `detect_state`,
  including model accessions repaired from an older base-id bug — re-run `gem_crosswalk` to
  confirm those against the model.
- Review flags are DERIVED, not stored as stale strings: `detect_state.flags.open` lists the
  ones whose condition is true RIGHT NOW, each with what resolves it. Fix the underlying
  problem and the flag disappears by itself; nothing needs clearing. `self_resolved` shows the
  ones that already closed. When a flag is real but cannot be fixed (KEGG only has a
  class-level entry for three lyso-PC species), `acknowledge_flag` records WHY it is
  acceptable, and it is then reported as acknowledged — never silently dropped.
- Finish a run with no open `action` flag: every one is either resolved or acknowledged.

## Stages (run only what the request needs; `detect_state` tells you what is pending)
1. ingest_names   — load + normalize raw names into the session ledger. Flags
                    parenthetical abbreviations (e.g. "TMAO(trimethylamine N-oxide)")
                    and combined "A/B" names for your attention.
2. exact_match    — MetaboAnalyst batch name->ID. Direct hits are high-confidence (M1).
2b. exact_match also returns `isomer_or_class_mismatch`: auto-accepted entries whose DB Match
                    name disagrees with the query on a discriminant axis, or is a class entry
                    for a species-level name. Re-decide those — an "exact" match is not an
                    identity when the matched name is a different isomer.
3. For each still-unmatched entry, gather EVIDENCE:
     - structure_lookup  (name -> PubChem CID / InChIKey / formula / monoisotopic mass)
     - search_synonym    (KEGG keggFind / PubChem / ChEBI OLS4 for YOUR proposed query
                          strings — abbreviation expansions, typo fixes, synonyms; every hit
                          comes back with an `isomer_verdict`, safest first. With
                          expand_shorthand=True the lipid/glycolipid shorthand variants of the
                          name are searched too: LysoPC(16:0) ->
                          1-palmitoyl-sn-glycero-3-phosphocholine, nLc4Cer ->
                          lactoneotetraosylceramide / paragloboside, GCDCA ->
                          chenodeoxyglycocholate)
     - mass_match_candidates (m/z + adduct -> formula candidates; weak evidence only)
   Then verify_candidate (molmass formula/mass consistency + xref cross-check) before
   accepting anything from a fuzzy/synonym/mass route. For LIPIDS and GLYCOLIPIDS formula/mass
   cannot decide anything (isomers are formula-identical) — use isomer_guard / id_name_check,
   and say in the rationale which axis settled it (linkage, ether bond, chain, locant).
4. bridge_xref    — promote a single accepted ID (e.g. PubChem CID) to the others
                    (KEGG/HMDB) via BridgeDb. A structure-only hit with no KEGG/HMDB
                    stays structure-only (not primary-usable for pathway/flux).
                    A bridged id is UNVERIFIED: BridgeDb links are wrong often enough that
                    C05472 was linked to taurochenodeoxycholate's HMDB. Back-check it.
4b. id_name_check — MANDATORY for any KEGG accepted through a search or a bridge: it fetches
                    the id's own DB record and compares the DB's name to yours. `conflict` =
                    wrong compound/isomer (find the right id, or drop the id); `class-level` =
                    a generic class entry answering a species name (keep it only if you report
                    it as class-level); `ok` = the DB itself states this compound.
5. record_decision — commit your CALL: accepted IDs, final_class, confidence, `origin`, or a
                     xenobiotic exclusion, with a one-line rationale. `origin` decides the class
                     for non-endogenous compounds (see origin taxonomy below).
5b. collision_check — no identifier may stand for two different compounds. `collision` = split
                    them (find the species-specific id for each); `class-id-shared` = the DB
                    only has a class entry (legitimate, but report those as class-level).
6. gem_crosswalk  — map accepted KEGG/HMDB/ChEBI IDs to model species. Expect the xref route to
                    MISS most lipids (a verified Recon3D run: 8 of 27 by xref, 21 of 28 in the
                    model). Every miss comes back as `name_suggestions` searched from the
                    model's own names, each with an isomer verdict.
6a. gem_search    — search the model's own name / id text / formula / mass. This is how the
                    species an xref cannot reach are found (GCDCA = `dgchol`
                    "Chenodeoxyglycocholate"; nLc4Cer = "Lactoneotetraosylceramide"; Recon3D
                    glycolipids as sugar compositions "(Gal)1 (Glc)1 (GlcNAc)1 (Cer)1").
                    `same_formula_groups` = formula-identical isomer candidates: decide those
                    on name/biology, never on mass.
6b. gem_assign    — commit the model CALL with an explicit relation: exact | class-proxy (the
                    model has only the generic R-group pool species — chain length/linkage not
                    represented) | isomer-surrogate (a DIFFERENT species stands in; state what
                    changes) | model-scope-absent (a real answer, not a failure) | id-gap.
                    Every model-relevant entry must end with one of these recorded.
6c. backfill_hmdb — bridge missing HMDB (KEGG/InChIKey/ChEBI/PubChem → HMDB) so HMDB
                    coverage is not under-counted relative to KEGG. Run before coverage. It
                    SKIPS an accession already accepted for another compound and reports it in
                    `skipped_conflicts` — resolve those, do not share an accession.
7. finalize_run   — stage 7 in ONE call, in order: coverage_summary (master ledger +
                    coverage tables) then ALWAYS emit:
                    DB-matching figures (db_matching_upset.png + enriched_xref.tsv;
                    db_matching_improvement.png); recovery-provenance tables
                    (mapping_provenance: kegg_recovered / hmdb_recovered / unmapped); and a
                    raw-API reproduction script (export_code → code/reproduce_mapping.py,
                    written with MetaboAnalystR/BridgeDbR/KEGGREST/PubChem/molmass/COBRApy/
                    matplotlib — NOT the tool wrappers). Each also runs standalone.
                    Each step is also its own tool and runs standalone; finalize_run is the
                    ORDER written down, not a computation with side effects — `coverage_summary`
                    alone computes the numbers and writes only its own two tables.
8. harness_audit  — READ-ONLY governance scorecard. Run LAST (after finalize_run):
                    verifies the reasoning layer honored THIS contract — no fabricated ID, no
                    mass-only (W) candidate used as primary, final_class↔confidence coherent,
                    fuzzy (M2/M3) accepts verified, locant-sensitive names re-checked, every
                    searched/bridged KEGG back-checked, no id contradicting its own name, no
                    identifier shared by two compounds, HMDB backfill collisions resolved,
                    xenobiotic-exclusion consistent, origin↔class coherent, flagged auto-accepts
                    reviewed, id-gap KEGG re-tried, a GEM relation recorded for every
                    model-relevant entry, surrogates documented, rationales present, nothing
                    pending, gem_crosswalk ran, all always-emit artifacts present. Fix every
                    `fail`; review each `warn`.

## final_class vocabulary
- KEGG-mapped / HMDB-mapped : endogenous, primary-usable (pathway + flux input).
- structure-only            : has PubChem/ChEBI/InChIKey but no KEGG/HMDB -> NOT primary.
- exogenous                 : biologically real but from OUTSIDE the host (diet/drug/microbial/
                              plant). KEPT + tagged; may carry KEGG/HMDB IDs and enter the GEM
                              crosswalk as boundary/exchange species. NOT the same as a
                              contaminant — this is a genuine signal, analysed separately.
- xenobiotic-excluded       : NON-biological (LC-MS additive / surfactant / plasticizer /
                              industrial / reagent contaminant) -> dropped from analysis.
- unmapped                  : no trustworthy identity found.

## origin taxonomy (the `origin` field on record_decision — decides exogenous vs xenobiotic)
Biological, outside-host  -> final_class 'exogenous' (KEPT):  diet · drug · microbial · plant
Non-biological / technical -> final_class 'xenobiotic-excluded' (EXCLUDED):
                             contaminant · industrial · additive · surfactant · plasticizer · reagent
Host-produced             -> origin 'endogenous' (default for KEGG/HMDB/structure-only).
screen_exogenous detects ONLY the non-biological classes (deterministic lexicon); the
diet/drug/microbial/plant call is YOUR judgement. Do not exclude a drug or dietary metabolite as
a contaminant — tag it 'exogenous' so it is kept and analysed as an outside-host signal.

## confidence tiers
- M1 exact/direct DB name match (MetaboAnalyst/KEGG/HMDB exact).
- M2 verified synonym or abbreviation expansion (structure + formula agree).
- M3 verified typo fix (edit-distance candidate, structure/formula agree).
- M4 structure-only (PubChem/ChEBI) with no biological DB xref.
- W  mass/adduct-only weak candidate (do NOT use as primary).
- X  reclassified non-endogenous (exogenous or xenobiotic-excluded).
- U  unresolved.

## GEM relation vocabulary (set by gem_assign; this is what a flux number is interpreted against)
- exact              : the model species IS this compound. Species-level input.
- class-proxy        : the model carries only the generic class / R-group pool species
                       (crm_hs, sphmyln_hs, ak2lgchol_hs). Usable, but chain length,
                       unsaturation and linkage are NOT represented — report it as class-level.
                       Beware mixed resolution: a ratio whose numerator is class-level and
                       denominator species-level is no longer the ratio the hypothesis stated.
- isomer-surrogate   : a DIFFERENT, related species stands in (N1-acetylspermine ->
                       N1-acetylspermidine, a different backbone). NOT an identity — state the
                       substitution wherever the flux is reported.
- model-scope-absent : the compound is genuinely not in the model. A RESULT, not a gap: record
                       it so the affected hypotheses are known to be unrepresentable.
- id-gap             : should be in the model but nothing has resolved it yet — keep working
                       (try the generic KEGG, then gem_search by name/formula/mass).

## isomer/class axes that decide lipid + glycolipid identity (isomer_guard / id_name_check)
skeleton (bile-acid vs corticosteroid …) · ether linkage (P- plasmenyl vs O- plasmanyl vs
1-acyl) · glycan series (neolacto/beta1-4 vs lacto/beta1-3) and length (Lc3 vs nLc4) ·
glycosidic linkages as written · sialyl linkage (a2,3 vs a2,6) · chains (16:0 / 20:4 / d18:1) ·
double-bond positions (arachidonoyl Delta5,8,11,14 vs n-3 Delta8,11,14,17) · omega · oxidation
vs PEROXidation (HETE vs HPETE) · hydroxy count (mono vs di) · acetyl count (mono vs di) ·
polyamine backbone (spermine vs spermidine). Formula and mass are BLIND to all of these.

## known gotchas (from verified runs)
- MetaboAnalystR needs InitDataObjects(default.dpi=72) or it errors.
- BridgeDb getDatabase() writes nothing; the DB must be downloaded already. loadDatabase
  on the 2 GB derby file is slow (~tens of seconds) — bridge_xref batches queries so the
  DB loads once per call.
- keggFind substring hits are noisy: always verify_candidate before accepting.
- PubChem REST rate-limits: calls are throttled; large batches take time.
- keggGet accepts at most 10 ids per request; id_name_check batches for you.
- A BridgeDb release can carry a wrong link (2021-01 links KEGG C05472 to
  taurochenodeoxycholate's HMDB). Bridged ids are leads, not identities.
- PubChem sometimes merges two compounds into one entry as synonyms (Lc3Cer and nLc4Cer on
  CID 131770449) — which is how two entries end up sharing an id. collision_check finds it.
- A DB entry's own synonym list can be internally inconsistent: KEGG C04910 lists BOTH
  "Lc4Cer" (beta1-3 lacto) and "Paragloboside" (which is neolacto). id_name_check therefore
  lets a conflicting synonym decide the verdict and returns both sides — when they disagree,
  read the systematic name/linkage, never the synonym that happens to match.
- GEMs annotate lipid species sparsely: expect the xref crosswalk to miss most lipids and
  plan on gem_search by name/formula. Recon3D stores most lipids as R-group POOL species and
  its glycolipids as sugar compositions.
- Model species ids carry a compartment suffix in two conventions (MAM01234c, tdchola_c);
  the accession to feed a model is the base id without it (`tdchola`, then `tdchola[c]`).
"""
