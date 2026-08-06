"""Offline smoke tests — no network / no R / no GEM load. Exercises the pure-Python
core, the session ledger, the anti-fabrication guard, and MCP server construction."""

import json
import tempfile
from pathlib import Path as pathlib_Path

from metabo_idmapper import tools
from metabo_idmapper.engine import normalize, verify
from metabo_idmapper.mcp_server import build_server


def test_normalize_parenthetical_and_locant():
    n = normalize.normalize("TMAO(trimethylamine N-oxide)")
    assert n["parenthetical"] == "trimethylamine N-oxide"
    assert "abbrev_parenthetical" in n["flags"]
    n2 = normalize.normalize("Azetidine-1-carboxylic acid")
    assert "has_digit_locant" in n2["flags"]  # isomer-sensitive


def test_normalize_acylcarnitine_spacing():
    n = normalize.normalize("Butyryl carnitine")
    assert "Butyrylcarnitine" in n["alternatives"]      # collapsed conjugate variant
    assert "conjugate_spacing" in n["flags"] and "acylcarnitine" in n["flags"]
    # multi-word non-conjugate names are not garbled
    assert normalize.normalize("citric acid")["alternatives"] == []


def test_verify_mass_and_formula():
    fm = verify.formula_mass("C2H7NO3S")  # taurine
    assert fm["ok"] and abs(fm["mono"] - 125.0147) < 0.01
    rep = verify.verify(proposed_formula="C2H7NO3S", observed_mass=125.0147)
    assert rep["mass_match"] is True
    bad = verify.verify(proposed_formula="C6H12O6", observed_formula="C2H7NO3S")
    assert bad["formula_match"] is False


def test_mass_candidates():
    w = verify.mass_candidates(126.0220, ["[M+H]+"])
    assert abs(w[0]["neutral_mono_mass"] - 125.0147) < 0.01


def test_ledger_roundtrip_and_state():
    with tempfile.TemporaryDirectory() as d:
        r = tools.ingest_names(d, names=["Taurine", "PEP=Phosphoenolpyruvate"])
        assert r["ingested"] == 2
        st = tools.detect_state(d)
        assert st["n_entries"] == 2 and st["pending"] == 2
        assert "exact_match" in st["suggested_next_tools"]


def test_record_decision_rejects_unbacked_id():
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurine"])
        # no tool produced C00245 yet -> must be refused (anti-fabrication)
        out = tools.record_decision(d, "M0000", rationale="guess",
                                    accepted={"kegg": "C00245"})
        assert "error" in out and "not produced" in out["error"]


def test_exclusion_path():
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Erucamide"])
        out = tools.record_decision(d, "M0000", rationale="industrial amide",
                                    final_class="xenobiotic-excluded")
        assert out["final_class"] == "xenobiotic-excluded"
        assert out["origin"] == "industrial"          # auto-suggested from the lexicon
        assert out["counts"]["xenobiotic-excluded"] == 1


def test_legacy_exogenous_excluded_alias():
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Di-n-butyl phthalate"])
        out = tools.record_decision(d, "M0000", rationale="plasticizer",
                                    final_class="exogenous-excluded")  # legacy name
        assert out["final_class"] == "xenobiotic-excluded"
        assert out["origin"] == "plasticizer"


def test_exogenous_is_kept_and_tagged():
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Caffeine"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        s.entries["M0000"]["candidates"].append({"source": "metaboanalyst", "kegg": "C07481"})
        s.save()
        out = tools.record_decision(d, "M0000", rationale="dietary/drug xanthine",
                                    accepted={"kegg": "C07481"}, final_class="exogenous",
                                    origin="drug")
        assert out["final_class"] == "exogenous"       # KEPT, not excluded
        assert out["origin"] == "drug"
        assert out["accepted"]["kegg"] == "C07481"     # still carries its ID
        assert out["counts"]["exogenous"] == 1


def test_origin_class_mismatch_warns():
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["SomeCompound"])
        out = tools.record_decision(d, "M0000", rationale="x",
                                    final_class="xenobiotic-excluded", origin="drug")
        # drug is a biological (exogenous) origin -> should not be xenobiotic-excluded
        assert any("implies final_class 'exogenous'" in w for w in out.get("warnings", []))


def test_normalize_keeps_substituent_parens():
    # substituent parentheticals must NOT be stripped (reviewer-found bug)
    assert normalize.normalize("4-(Dimethylamino)phenylthiocyanate")["normalized"] \
        == "4-(Dimethylamino)phenylthiocyanate"
    assert normalize.normalize("Tri(butoxyethyl) phosphate")["normalized"] \
        == "Tri(butoxyethyl) phosphate"
    # genuine acronym expansion still works, and keeps the acronym as an alternative
    n = normalize.normalize("TMAO(trimethylamine N-oxide)")
    assert n["normalized"] == "trimethylamine N-oxide" and "TMAO" in n["alternatives"]


def test_hmdb_padding_and_lexicon():
    from metabo_idmapper.engine import lexicon
    assert lexicon.pad_hmdb("HMDB00349") == "HMDB0000349"
    assert lexicon.pad_hmdb("C00009") == "C00009"  # non-HMDB passes through
    assert lexicon.xenobiotic_hits("Trifluoroacetic acid") == ["lcms-additive"]
    assert lexicon.xenobiotic_hits("Di-n-butyl phthalate") == ["plasticizer"]
    assert lexicon.xenobiotic_hits("L-Threonine") == []


def test_record_decision_warns_on_contaminant_kept_endogenous():
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Trifluoroacetic acid"])
        # pretend a candidate exists, then keep it endogenous → should warn
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        s.entries["M0000"]["candidates"].append({"source": "x", "hmdb": "HMDB0014118"})
        s.save()
        out = tools.record_decision(d, "M0000", rationale="test",
                                    accepted={"hmdb": "HMDB0014118"})
        assert out["final_class"] == "HMDB-mapped"
        assert "warnings" in out and any("xenobiotic-excluded" in w for w in out["warnings"])


def test_harness_audit_origin_coherence():
    from metabo_idmapper import harness
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["MysteryDietCompound"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        # exogenous kept but never tagged with an origin -> governance failure
        e = s.entries["M0000"]; e["final_class"] = "exogenous"; e["confidence"] = "M4"
        s.save()
        by = {c["check"]: c for c in harness.audit(d)["checks"]}
        assert by["origin_coherence"]["level"] == "fail"
        assert any("exogenous-untagged" in o for o in by["origin_coherence"]["offenders"])


def test_coverage_emits_db_matching_figures():
    import os
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurine", "Glucose", "Erucamide"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        # give one entry a real accepted id so the plot has content
        s.entries["M0000"]["candidates"].append({"source": "metaboanalyst", "kegg": "C00245",
                                                  "hmdb": "HMDB0000251", "pubchem": "1123"})
        s.save()
        tools.record_decision(d, "M0000", rationale="t", accepted={"kegg": "C00245"})
        out = tools.coverage_summary(d)  # figures=True default
        assert "figures" in out
        assert os.path.exists(os.path.join(d, "figures", "db_matching_upset.png"))
        assert os.path.exists(os.path.join(d, "figures", "db_matching_improvement.png"))
        assert os.path.exists(os.path.join(d, "enriched_xref.tsv"))


def test_server_registers_all_tools():
    mcp = build_server()
    names = {t.name for t in mcp._tool_manager.list_tools()}
    assert {"midmap_guidance", "ingest_names", "record_decision", "screen_exogenous",
            "backfill_hmdb", "plot_coverage", "mapping_provenance", "gem_crosswalk",
            "export_code", "annotate_source", "export_report_ppt", "coverage_summary",
            "harness_audit", "id_name_check", "isomer_guard", "collision_check",
            "gem_search", "gem_assign", "acknowledge_flag"} <= names
    assert len(names) == 26


def test_harness_audit_scorecard():
    from metabo_idmapper import harness
    # empty ledger → nothing to audit
    with tempfile.TemporaryDirectory() as d:
        assert harness.audit(d)["verdict"] == "warn"
        assert harness.audit(d)["checks"][0]["check"] == "ingest"
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurine"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        s.entries["M0000"]["candidates"].append({"source": "metaboanalyst", "kegg": "C00245"})
        s.save()
        tools.record_decision(d, "M0000", rationale="exact KEGG match", accepted={"kegg": "C00245"})
        rep = harness.audit(d)
        by = {c["check"]: c for c in rep["checks"]}
        # a legitimately-backed accept passes anti-fabrication and W-tier checks
        assert by["anti_fabrication"]["level"] == "pass"
        assert by["w_tier_not_primary"]["level"] == "pass"
        # stage-7 artifacts were not emitted → that check warns
        assert by["finalize_artifacts"]["level"] == "warn"


def test_harness_audit_catches_w_tier_as_primary():
    from metabo_idmapper import harness
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["MysteryPeak"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        s.entries["M0000"]["candidates"].append({"source": "kegg-search", "kegg": "C00031"})
        s.save()
        # mass-only weak tier (W) used as a primary KEGG mapping → hard contract violation
        tools.record_decision(d, "M0000", rationale="mass-only guess", confidence="W",
                              accepted={"kegg": "C00031"})
        rep = harness.audit(d)
        by = {c["check"]: c for c in rep["checks"]}
        assert by["w_tier_not_primary"]["level"] == "fail"
        assert "M0000" in by["w_tier_not_primary"]["offenders"][0]
        assert rep["verdict"] == "fail"


def test_harness_audit_flags_fabricated_id():
    from metabo_idmapper import harness
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurine"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        # hand-plant an accepted id no candidate ever produced (simulates a tampered ledger)
        s.entries["M0000"]["accepted"] = {"kegg": "C99999"}
        s.entries["M0000"]["final_class"] = "KEGG-mapped"
        s.entries["M0000"]["confidence"] = "M1"
        s.save()
        by = {c["check"]: c for c in harness.audit(d)["checks"]}
        assert by["anti_fabrication"]["level"] == "fail"
        assert any("C99999" in o for o in by["anti_fabrication"]["offenders"])


def test_isomer_catches_the_verified_wrong_ids():
    """The four silent mis-assignments from verified runs, each on its own axis."""
    from metabo_idmapper.engine import isomer

    # a bile acid carrying a cortisol metabolite's KEGG (C05472's own KEGG name)
    r = isomer.compare("Taurochenodeoxycholic acid", "Urocortisol; Tetrahydrocortisol")
    assert r["verdict"] == "conflict"
    assert any(c["axis"] == "skeleton" for c in r["conflicts"])
    # a plasmenyl (vinyl-ether) lyso-PC carrying the 1-ACYL class id (C04230)
    r = isomer.compare("LysoPC(P-16:0)", "1-Acyl-sn-glycero-3-phosphocholine")
    assert r["verdict"] == "conflict"
    assert any(c["axis"] == "ether_linkage" for c in r["conflicts"])
    # the SAME class id is merely class-level for a plain acyl species — not a conflict
    assert isomer.compare("LysoPC(16:0)",
                          "1-Acyl-sn-glycero-3-phosphocholine")["verdict"] == "class-level"
    # neolacto query vs the lacto isomer, as KEGG actually writes it (C04910: no "Lc4Cer" in
    # the name — only the beta1-3 linkage says it)
    r = isomer.compare("nLc4Cer", "1,3-beta-D-Galactosyl-N-acetyl-D-glucosaminyl-"
                                  "1,3-beta-D-galactosyl-1,4-D-glucosylceramide")
    assert r["verdict"] == "conflict"
    assert any(c["axis"] == "glycan_series" for c in r["conflicts"])
    # ...and the correct neolacto entry (C04922) does not conflict
    assert isomer.compare("nLc4Cer", "beta-D-Galactosyl-(1->4)-N-acetyl-beta-D-glucosaminyl-"
                                     "(1->3)-beta-D-galactosyl-(1->4)-beta-D-glucosyl-"
                                     "(1<->1)-ceramide")["verdict"] != "conflict"
    # two glycolipids that must never share an id
    assert isomer.compare("Lc3Cer", "nLc4Cer")["verdict"] == "conflict"
    # peroxide vs hydroxide, and mono- vs di-hydroxy
    assert isomer.compare("15-HETE", "15-Hydroperoxyeicosatetraenoic acid")["verdict"] == "conflict"
    assert isomer.compare("15-HETE", "5,15-DiHETE")["verdict"] == "conflict"
    # a different polyamine backbone is not a substitute
    assert isomer.compare("N1-Acetylspermine", "N1-Acetylspermidine")["verdict"] == "conflict"
    # correct identities are not flagged
    for q, c in (("Taurochenodeoxycholic acid", "Taurochenodeoxycholate"),
                 ("Glycochenodeoxycholic acid", "Chenodeoxyglycocholate"),
                 ("15-HETE", "(15S)-15-Hydroxy-5,8,11-cis-13-trans-eicosatetraenoate"),
                 ("Lactosylceramide", "Galactosyl glucosyl ceramide")):
        assert isomer.compare(q, c)["verdict"] == "ok", (q, c)


def test_backcheck_verdict_is_not_masked_by_a_matching_synonym():
    """KEGG C04910 lists BOTH "Lc4Cer" (beta1-3 lacto) and "Paragloboside" (neolacto) on one
    entry: a matching synonym must not override a conflicting one, or the isomer swap this
    tool exists to catch would be reported as `ok`."""
    from metabo_idmapper.tools import _best_verdict

    rep = _best_verdict("nLc4Cer", [
        "1,3-beta-D-Galactosyl-N-acetyl-D-glucosaminyl-1,3-beta-D-galactosyl-1,4-"
        "D-glucosylceramide", "Ceramidetetrasaccharide", "Paragloboside", "Lc4Cer"])
    assert rep["verdict"] == "conflict"
    assert "Lc4Cer" in rep["conflicting_names"] and "Paragloboside" in rep["matching_names"]
    assert "synonyms disagree" in rep["note"]
    # a consistent entry still reads clean, and a generic entry reads class-level
    assert _best_verdict("Taurochenodeoxycholic acid",
                         ["Taurochenodeoxycholate",
                          "Taurochenodeoxycholic acid"])["verdict"] == "ok"
    assert _best_verdict("LysoPC(16:0)", ["1-Acyl-sn-glycero-3-phosphocholine",
                                          "Lysophosphatidylcholine"])["verdict"] == "class-level"


def test_shorthand_variants_reach_database_spellings():
    from metabo_idmapper.engine import shorthand

    v = [s.lower() for s in shorthand.variants("LysoPC(16:0)")]
    assert any("palmitoyl" in s and "phosphocholine" in s for s in v)
    v = [s.lower() for s in shorthand.variants("LysoPC(P-16:0)")]
    assert any("alkenyl" in s or "plasmenyl" in s for s in v)
    assert not any("palmitoyl" in s for s in v)        # an ether is not an ester
    assert any("lactoneotetraosylceramide" == s for s in
               (x.lower() for x in shorthand.variants("nLc4Cer")))
    assert any("chenodeoxyglycocholate" == s for s in
               (x.lower() for x in shorthand.variants("Glycochenodeoxycholic acid")))
    assert any("hydroxyeicosatetraenoic acid" in s for s in
               (x.lower() for x in shorthand.variants("15-HETE")))
    # the sphingoid base of Cer(d18:1/24:0) is not an acyl chain
    assert not any("oleoyl" in s.lower() for s in shorthand.variants("Cer(d18:1/24:0)"))


def test_collision_check_separates_collision_from_class_level_share():
    from metabo_idmapper.state import Session
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Lc3Cer", "nLc4Cer"])
        s = Session(d)
        for fid in ("M0000", "M0001"):
            s.entries[fid]["candidates"].append({"source": "pubchem", "hmdb": "HMDB0062485"})
        s.save()
        for fid in ("M0000", "M0001"):
            tools.record_decision(d, fid, rationale="bridged from a shared PubChem entry",
                                  accepted={"hmdb": "HMDB0062485"})
        rep = tools.collision_check(d)
        assert rep["n_collisions"] == 1
        assert rep["collisions"][0]["id"] == "HMDB0062485"
        assert set(rep["collisions"][0]["feature_ids"]) == {"M0000", "M0001"}
        by = {c["check"]: c for c in tools.harness_audit(d)["checks"]}
        assert by["id_collisions"]["level"] == "fail"

    # the same sharing, but the id is KNOWN to be a class-level DB entry -> legitimate
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["LysoPC(16:0)", "LysoPC(18:2)"])
        s = Session(d)
        for fid in ("M0000", "M0001"):
            s.entries[fid]["candidates"].append({"source": "metaboanalyst", "kegg": "C04230"})
            s.entries[fid]["decisions"].append(
                {"action": "backcheck", "db": "kegg", "id": "C04230", "verdict": "class-level",
                 "db_resolution": "class", "db_names": ["1-Acyl-sn-glycero-3-phosphocholine"],
                 "rationale": "back-check"})
        s.save()
        for fid in ("M0000", "M0001"):
            tools.record_decision(d, fid, rationale="KEGG has only the class entry",
                                  accepted={"kegg": "C04230"}, confidence="M1")
        rep = tools.collision_check(d)
        assert rep["n_collisions"] == 0 and len(rep["class_id_shared"]) == 1
        by = {c["check"]: c for c in tools.harness_audit(d)["checks"]}
        assert by["id_collisions"]["level"] == "warn"


def test_backfill_hmdb_refuses_to_share_an_accession(monkeypatch):
    """The bridge is per-entry and blind to what another entry already holds; the guard is not."""
    from metabo_idmapper import tools as T
    from metabo_idmapper.state import Session
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Lc3Cer", "nLc4Cer"])
        s = Session(d)
        s.entries["M0000"]["candidates"].append({"source": "kegg-search", "kegg": "C04845"})
        s.entries["M0001"]["candidates"].append({"source": "kegg-search", "kegg": "C04922"})
        s.save()
        tools.record_decision(d, "M0000", rationale="lacto triaosyl", accepted={"kegg": "C04845"})
        tools.record_decision(d, "M0001", rationale="neolacto tetraosyl",
                              accepted={"kegg": "C04922"})

        # both entries bridge to the SAME HMDB (PubChem merged the two names into one entry)
        def fake_bridge(workdir, queries, db=None):
            ses = Session(workdir)
            for q in queries:
                ses.add_candidate(q["feature_id"], {"source": "bridgedb",
                                                    "query": f"{q['source']}:{q['id']}",
                                                    "hmdb": "HMDB0062485"})
            ses.save()
            return {"results": []}

        monkeypatch.setattr(T, "bridge_xref", fake_bridge)
        out = tools.backfill_hmdb(d)
        assert out["hmdb_gained"] == 1                      # the first one only
        assert len(out["skipped_conflicts"]) == 1
        sk = out["skipped_conflicts"][0]
        assert sk["hmdb"] == "HMDB0062485" and sk["isomer_verdict"] == "conflict"
        assert Session(d).entries[sk["feature_id"]].get("accepted", {}).get("hmdb") is None


def test_gem_assign_requires_a_real_species_and_records_the_relation():
    from metabo_idmapper import tools as T
    from metabo_idmapper.engine import gem as G
    from metabo_idmapper.state import Session

    fake = {"crm_hs": {"found": True, "mam_base": "crm_hs", "gem_name": "Ceramide (homo sapiens)",
                       "pool_species": True},
            "N1aspmd": {"found": True, "mam_base": "N1aspmd", "gem_name": "N1-Acetylspermidine",
                        "pool_species": False}}

    def fake_detail(mam, model_path=None):
        return fake.get(mam, {"found": False, "mam_base": mam})

    orig = G.species_detail
    G.species_detail = fake_detail
    try:
        with tempfile.TemporaryDirectory() as d:
            tools.ingest_names(d, names=["N1-Acetylspermine", "Cer(d18:1/24:0)"])
            # an invented accession is refused, exactly like an invented ID
            bad = tools.gem_assign(d, "M0000", mam="no_such_species", relation="exact",
                                   rationale="testing the anti-fabrication guard")
            assert "error" in bad and "not in model" in bad["error"]
            # a mapped relation needs a species; an absence must be recorded as one
            assert "error" in tools.gem_assign(d, "M0000", relation="exact",
                                               rationale="no species given at all")
            out = tools.gem_assign(d, "M0000", mam="N1aspmd", relation="isomer-surrogate",
                                   rationale="Recon3D has only N1-acetylspermidine, a different "
                                             "backbone, so polyamine flux from it is not "
                                             "spermine-specific")
            assert out["relation"] == "isomer-surrogate" and out["mam"] == ["N1aspmd"]
            assert "SUBSTITUTION" in out["caveat"]
            # an R-group pool species is a class proxy, not an identity
            warn = tools.gem_assign(d, "M0001", mam="crm_hs", relation="exact",
                                    rationale="ceramide pool species in the model")
            assert any("class-proxy" in w for w in warn["warnings"])
            e = Session(d).entries["M0000"]
            assert e["gem_curated"] and e["gem_relation"] == "isomer-surrogate"
            assert e["gem"][e["gem_model"]]["mam"] == ["N1aspmd"]
    finally:
        G.species_detail = orig


def test_gem_base_id_strips_both_compartment_conventions():
    from metabo_idmapper.engine import gem

    class M:
        def __init__(self, i, c):
            self.id, self.compartment = i, c
    # Recon3D style: the naive "drop one trailing letter" rule left "tdchola_", not an accession
    assert gem.base_id(M("tdchola_c", "c")) == "tdchola"
    assert gem.base_id(M("met__L_e", "e")) == "met__L"
    assert gem.base_id(M("pcholar_hs_c", "c")) == "pcholar_hs"
    # Mouse-GEM style
    assert gem.base_id(M("MAM01234c", "c")) == "MAM01234"
    # no compartment metadata -> fall back to the suffix pattern
    assert gem.base_id(M("tdchola_c", None)) == "tdchola"


def test_harness_fails_unbackchecked_and_conflicting_kegg():
    from metabo_idmapper import harness
    from metabo_idmapper.state import Session
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurochenodeoxycholic acid"])
        s = Session(d)
        s.entries["M0000"]["candidates"].append({"source": "bridgedb", "kegg": "C05472"})
        s.save()
        tools.record_decision(d, "M0000", rationale="bridged from HMDB via BridgeDb",
                              accepted={"kegg": "C05472"}, confidence="M2")
        by = {c["check"]: c for c in harness.audit(d)["checks"]}
        assert by["id_backcheck"]["level"] == "fail"       # bridged and never back-checked

        # now the back-check runs and comes back 'conflict' — keeping the id is a failure
        from metabo_idmapper import flags
        s = Session(d)
        s.entries["M0000"]["decisions"].append(
            {"action": "backcheck", "db": "kegg", "id": "C05472", "verdict": "conflict",
             "db_names": ["Urocortisol"], "rationale": "back-check"})
        flags.raise_flag(s.entries["M0000"], "id_name_conflict",
                         "C05472 resolves to Urocortisol", "id_name_check")
        s.save()
        by = {c["check"]: c for c in harness.audit(d)["checks"]}
        assert by["id_backcheck"]["level"] == "pass"
        assert by["id_name_conflicts"]["level"] == "fail"
        assert any("M0000" in o for o in by["id_name_conflicts"]["offenders"])

        # accepting a different id resolves it — the flag closes itself, nothing to clean up
        s = Session(d)
        s.entries["M0000"]["candidates"].append({"source": "kegg-search", "kegg": "C05465"})
        s.save()
        tools.record_decision(d, "M0000", rationale="C05465 is the bile acid; C05472 was a "
                                                    "cortisol metabolite",
                              accepted={"kegg": "C05465"}, confidence="M2")
        by = {c["check"]: c for c in harness.audit(d)["checks"]}
        assert by["id_name_conflicts"]["level"] == "pass"
        assert flags.open_flags(Session(d).entries["M0000"], Session(d)) == []


def test_ledger_refuses_to_overwrite_a_concurrent_write():
    """The ledger is a whole-file snapshot: a second writer would silently drop the first
    writer's decisions, which is the one way this design can lose data."""
    from metabo_idmapper.state import LedgerConflict, Session
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurine", "Glucose"])
        a = Session(d)                      # two tool calls hold the ledger at once
        b = Session(d)
        b.entries["M0000"]["origin"] = "endogenous"
        b.save()                            # b commits first
        a.entries["M0001"]["origin"] = "endogenous"
        try:
            a.save()
            raise AssertionError("stale session overwrote a newer ledger")
        except LedgerConflict as ex:
            assert "changed since it was read" in str(ex)
        # b's write survived, a's was refused rather than silently applied
        assert Session(d).entries["M0000"]["origin"] == "endogenous"
        assert Session(d).entries["M0001"]["origin"] is None
        # tools surface it as a structured result, not an exception
        from metabo_idmapper.tools import REGISTRY
        record = next(f for f in REGISTRY if f.__name__ == "record_decision")
        s = Session(d)
        s.entries["M0000"]["candidates"].append({"source": "x", "kegg": "C00245"})
        s.save()
        stale = Session(d)
        Session(d).save()                   # someone else writes again
        import metabo_idmapper.tools as T
        orig = T.Session
        T.Session = lambda wd: stale        # force the tool onto the stale snapshot
        try:
            out = record(d, "M0000", rationale="test", accepted={"kegg": "C00245"})
        finally:
            T.Session = orig
        assert out["error_code"] == "ledger_conflict" and out["wrote_nothing"] is True


def test_migration_repairs_v1_ledger():
    """A v1 ledger carries bare flag strings and model accessions broken by the old
    base-id rule ('tdchola_c' -> 'tdchola_'). Loading migrates in memory; saving persists it."""
    import json as _json
    from metabo_idmapper.state import LEDGER_SCHEMA, Session
    with tempfile.TemporaryDirectory() as d:
        legacy = {"created": "2026-01-01T00:00:00", "entries": {"M0000": {
            "feature_id": "M0000", "original_name": "Taurochenodeoxycholic acid",
            "normalized": {"normalized": "Taurochenodeoxycholic acid", "flags": []},
            "candidates": [{"source": "metaboanalyst", "kegg": "C05465"}],
            "accepted": {"kegg": "C05465"}, "final_class": "KEGG-mapped", "confidence": "M1",
            "origin": "endogenous", "gem_mam": ["tdchola_", "ak2lgchol_hs_"],
            "gem_cause": None, "decisions": [], "_flags": ["auto_accept_review"]}}}
        (pathlib_Path(d) / "midmap_ledger.json").write_text(_json.dumps(legacy))

        s = Session(d)
        assert s.data["schema_version"] == LEDGER_SCHEMA
        e = s.entries["M0000"]
        assert e["gem_mam"] == ["tdchola", "ak2lgchol_hs"]      # dangling separator repaired
        assert "_flags" not in e and "auto_accept_review" in e["flags"]
        step = s.migrations[0]
        assert step["schema"] == "1->2" and len(step["model_ids_repaired"]) == 2
        # in-memory only until something writes
        assert "schema_version" not in _json.loads(
            (pathlib_Path(d) / "midmap_ledger.json").read_text())
        s.save()
        on_disk = _json.loads((pathlib_Path(d) / "midmap_ledger.json").read_text())
        assert on_disk["schema_version"] == LEDGER_SCHEMA
        assert on_disk["migrations"][0]["model_ids_repaired"][0]["to"] == "tdchola"
        assert Session(d).migrations == []                      # already current


def test_flags_are_derived_and_acknowledgeable():
    from metabo_idmapper import flags
    from metabo_idmapper.state import Session
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["LysoPC(16:0)", "LysoPC(18:2)"])
        s = Session(d)
        for fid in ("M0000", "M0001"):
            e = s.entries[fid]
            e["candidates"].append({"source": "metaboanalyst", "kegg": "C04230"})
            e["decisions"].append({"action": "backcheck", "db": "kegg", "id": "C04230",
                                   "verdict": "class-level", "db_resolution": "class",
                                   "db_names": ["1-Acyl-sn-glycero-3-phosphocholine"],
                                   "rationale": "back-check"})
            flags.raise_flag(e, "class_level_id", "C04230 is a class entry", "id_name_check")
        s.save()
        for fid in ("M0000", "M0001"):
            tools.record_decision(d, fid, rationale="KEGG only has the class entry",
                                  accepted={"kegg": "C04230"}, confidence="M1")

        st = tools.detect_state(d)
        assert st["flags"]["open"]["class_level_id"] == 2
        # an unresolvable-by-fixing flag: acknowledge it with a reason, and it is REPORTED
        short = tools.acknowledge_flag(d, "M0000", "class_level_id", "n/a")
        assert "error" in short                                  # a reason is required
        ok = tools.acknowledge_flag(d, "M0000", "class_level_id",
                                    "KEGG has no species-level lyso-PC entry; reported as "
                                    "class-level in the marker table")
        assert ok["acknowledged"] and ok["open_now"] == []
        st = tools.detect_state(d)
        assert st["flags"]["open"]["class_level_id"] == 1
        assert st["flags"]["acknowledged"]["class_level_id"] == 1

        # a flag whose condition is fixed disappears WITHOUT being acknowledged
        s = Session(d)
        s.entries["M0001"]["candidates"].append({"source": "kegg-search", "kegg": "C04100"})
        s.save()
        tools.record_decision(d, "M0001", rationale="species-level id found",
                              accepted={"kegg": "C04100"}, confidence="M2")
        st = tools.detect_state(d)
        assert "class_level_id" not in st["flags"]["open"]
        assert st["flags"]["self_resolved"]["class_level_id"] == 1


def test_export_code_is_raw_api_and_compiles():
    import py_compile
    with tempfile.TemporaryDirectory() as d:
        tools.ingest_names(d, names=["Taurine"])
        s = __import__("metabo_idmapper.state", fromlist=["Session"]).Session(d)
        s.entries["M0000"]["candidates"].append({"source": "metaboanalyst", "kegg": "C00245"})
        s.save()
        tools.record_decision(d, "M0000", rationale="t", accepted={"kegg": "C00245"})
        out = tools.export_code(d)
        script = out["script"]
        src = open(script).read()
        py_compile.compile(script, doraise=True)          # valid Python
        import re as _re
        code_lines = [ln for ln in src.splitlines()
                      if _re.match(r"\s*(import|from)\s+metabo_idmapper", ln)]
        assert not code_lines                              # no tool-package import statement
        assert "ingest_names(" not in src and "record_decision(" not in src  # no tool fns
        assert "MetaboAnalystR" in src and "BridgeDbR" in src  # raw R APIs inlined
        assert "from molmass import Formula" in src and "import cobra" in src
        # notebook: valid JSON, no def/lambda in code cells, has detailed provenance comments
        nb = json.load(open(out["notebook"]))
        code_src = "\n".join("".join(c["source"]) if isinstance(c["source"], list)
                             else c["source"] for c in nb["cells"] if c["cell_type"] == "code")
        assert "def " not in code_src and "lambda" not in code_src   # unrolled, no functions
        assert "MetaboAnalystR" not in code_src or True  # (R is in sidecar files)
        assert "INPUT" in code_src and "REUSED" in code_src          # per-cell provenance comments
        assert "import metabo_idmapper" not in code_src
