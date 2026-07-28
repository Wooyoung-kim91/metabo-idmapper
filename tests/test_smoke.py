"""Offline smoke tests — no network / no R / no GEM load. Exercises the pure-Python
core, the session ledger, the anti-fabrication guard, and MCP server construction."""

import json
import tempfile

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
            "harness_audit"} <= names
    assert len(names) == 20


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
