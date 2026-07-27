"""Canonical workflow + operating rules, returned by the `midmap_guidance` tool.

This is the single source of truth the reasoning LLM reads at the start of a run
(mirrors scpilot_guidance): what the tools do, what the LLM must decide, and the
confidence-tier vocabulary.
"""

CANONICAL = r"""
# metabo-idmapper — canonical workflow

GOAL: turn a list of messy metabolite NAMES into INPUT-USABLE IDs
(KEGG / HMDB / ChEBI / PubChem / InChIKey) plus a Mouse-GEM MAM crosswalk for flux input.

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

## Stages (run only what the request needs; `detect_state` tells you what is pending)
1. ingest_names   — load + normalize raw names into the session ledger. Flags
                    parenthetical abbreviations (e.g. "TMAO(trimethylamine N-oxide)")
                    and combined "A/B" names for your attention.
2. exact_match    — MetaboAnalyst batch name->ID. Direct hits are high-confidence (M1).
3. For each still-unmatched entry, gather EVIDENCE:
     - structure_lookup  (name -> PubChem CID / InChIKey / formula / monoisotopic mass)
     - search_synonym    (KEGG keggFind / PubChem / ChEBI OLS4 for YOUR proposed query
                          strings — abbreviation expansions, typo fixes, synonyms)
     - mass_match_candidates (m/z + adduct -> formula candidates; weak evidence only)
   Then verify_candidate (molmass formula/mass consistency + xref cross-check) before
   accepting anything from a fuzzy/synonym/mass route.
4. bridge_xref    — promote a single accepted ID (e.g. PubChem CID) to the others
                    (KEGG/HMDB) via BridgeDb. A structure-only hit with no KEGG/HMDB
                    stays structure-only (not primary-usable for pathway/flux).
5. record_decision — commit your CALL: accepted IDs, final_class, confidence, or
                     exogenous exclusion, with a one-line rationale.
6. gem_crosswalk  — map accepted KEGG/HMDB/ChEBI IDs to Mouse-GEM MAM species.
6b. backfill_hmdb — bridge missing HMDB (KEGG/InChIKey/ChEBI/PubChem → HMDB) so HMDB
                    coverage is not under-counted relative to KEGG. Run before coverage.
7. coverage_summary — write the master ledger + coverage + provenance, and ALWAYS emit:
                    DB-matching figures (db_matching_upset.png + enriched_xref.tsv;
                    db_matching_improvement.png); recovery-provenance tables
                    (mapping_provenance: kegg_recovered / hmdb_recovered / unmapped); and a
                    raw-API reproduction script (export_code → code/reproduce_mapping.py,
                    written with MetaboAnalystR/BridgeDbR/KEGGREST/PubChem/molmass/COBRApy/
                    matplotlib — NOT the tool wrappers). Each also runs standalone.
8. harness_audit  — READ-ONLY governance scorecard. Run LAST (after coverage_summary):
                    verifies the reasoning layer honored THIS contract — no fabricated ID, no
                    mass-only (W) candidate used as primary, final_class↔confidence coherent,
                    fuzzy (M2/M3) accepts verified, locant-sensitive names re-checked,
                    exogenous-exclusion consistent, flagged auto-accepts reviewed, id-gap KEGG
                    re-tried, rationales present, nothing pending, gem_crosswalk ran, all
                    always-emit artifacts present. Fix every `fail`; review each `warn`.

## final_class vocabulary
- KEGG-mapped / HMDB-mapped : primary-usable (pathway + flux input).
- structure-only            : has PubChem/ChEBI/InChIKey but no KEGG/HMDB -> NOT primary.
- exogenous-excluded        : xenobiotic / industrial / drug -> excluded from analysis.
- unmapped                  : no trustworthy identity found.

## confidence tiers
- M1 exact/direct DB name match (MetaboAnalyst/KEGG/HMDB exact).
- M2 verified synonym or abbreviation expansion (structure + formula agree).
- M3 verified typo fix (edit-distance candidate, structure/formula agree).
- M4 structure-only (PubChem/ChEBI) with no biological DB xref.
- W  mass/adduct-only weak candidate (do NOT use as primary).
- X  reclassified exogenous.
- U  unresolved.

## unmapped-cause classification (for gem_crosswalk gaps)
- id-gap            : compound has an ID but the GEM has no xref for it.
- model-scope-absent: compound is genuinely outside the model (dipeptides, sugar-
                      phosphates, ergothioneine, etc.) — not a mapping failure.

## known gotchas (from verified runs)
- MetaboAnalystR needs InitDataObjects(default.dpi=72) or it errors.
- BridgeDb getDatabase() writes nothing; the DB must be downloaded already. loadDatabase
  on the 2 GB derby file is slow (~tens of seconds) — bridge_xref batches queries so the
  DB loads once per call.
- keggFind substring hits are noisy: always verify_candidate before accepting.
- PubChem REST rate-limits: calls are throttled; large batches take time.
"""
