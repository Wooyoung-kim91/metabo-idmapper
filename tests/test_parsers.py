"""Parser tests against RECORDED responses — the R and REST boundary, offline.

Everything these parsers read comes from outside the package: four R engines writing JSON, and
three REST endpoints. Until now they were only ever exercised by a live run, so a change in a
column name, a wrapper key or an "NA" convention broke nothing until someone happened to run
the real thing — and then broke it silently, as a mapping that quietly came back empty.

The fixtures in `tests/fixtures/` are real responses captured from those engines (KEGGREST
keggGet/keggFind, MetaboAnalystR, BridgeDbR against the 2 GB derby DB, PubChem PUG-REST, ChEBI
OLS4). Each test replaces only the TRANSPORT — `run_r`, the HTTP getter — and asserts what the
parser makes of the payload. No network, no R, no model load.

Recording them found two real defects, both now pinned by the tests below: MetaboAnalystR's
injected anchors were dropping a caller's own Taurine / glutamate / glucose rows, and a
BridgeDb answer with several distinct accessions was silently reduced to its first element.
"""

import json
from pathlib import Path

from metabo_idmapper import tools
from metabo_idmapper.engine import chebi, entry, structure

FIXTURES = Path(__file__).parent / "fixtures"


def load(name):
    return json.loads((FIXTURES / name).read_text())


# --------------------------------------------------------------------- R engines
def test_kegg_entry_parser(monkeypatch):
    """keggGet rows -> {id: {found, names, formula, xrefs}}, including a miss."""
    monkeypatch.setattr(entry, "run_r",
                        lambda script, payload, **kw: {"ok": True,
                                                       "result": load("r_kegg_entry.json")})
    recs = entry.kegg_entries(["C05472", "C04230", "C99999"])
    assert recs["C05472"]["found"] is True
    # the record that made this tool necessary: the id is a cortisol metabolite
    assert recs["C05472"]["names"][0] == "Urocortisol"
    assert recs["C05472"]["formula"] == "C21H34O5"
    assert recs["C05472"]["chebi"] == "28320" and recs["C05472"]["pubchem"] == "7832"
    assert recs["C04230"]["names"][0].startswith("1-Acyl-sn-glycero-3-phosphocholine")
    assert recs["C99999"]["found"] is False and recs["C99999"]["names"] == []


def test_kegg_entry_parser_survives_an_engine_failure(monkeypatch):
    """An R failure is not a data verdict: every id comes back not-found WITH the reason."""
    monkeypatch.setattr(entry, "run_r",
                        lambda *a, **k: {"ok": False, "error": "R exit 1", "stderr": "boom"})
    recs = entry.kegg_entries(["C00245"])
    assert recs["C00245"]["found"] is False
    assert recs["C00245"]["error"] == "R exit 1" and recs["C00245"]["names"] == []


def test_kegg_search_rows_become_candidates_with_an_isomer_verdict(monkeypatch, tmp_path):
    """keggFind rows carry kegg + kegg_name; search_synonym screens the name."""
    monkeypatch.setattr(tools, "run_r",
                        lambda script, payload, **kw: {"ok": True,
                                                       "result": load("r_kegg_search.json")})
    monkeypatch.setattr(tools._struct, "lookup", lambda q: {"found": False})
    monkeypatch.setattr(tools._chebi, "search", lambda q, rows=6: [])
    d = str(tmp_path)
    tools.ingest_names(d, names=["Taurochenodeoxycholic acid"])
    out = tools.search_synonym(d, "M0000", queries=["taurochenodeoxycholate"])
    hit = out["candidates"]["kegg"][0]
    assert hit["kegg"] == "C05465"
    assert hit["isomer_verdict"] == "ok"          # the DB's own name agrees with the query
    assert out["n_isomer_conflicts"] == 0


def test_metaboanalyst_rows_keep_the_callers_own_anchor_names(monkeypatch, tmp_path):
    """The R engine injects Taurine / glutamate / glucose so CrossReferencing cannot die on an
    all-miss list, then drops those rows. Dropping them when the CALLER asked for one returned
    no mapping for three of the most common metabolites there are — this pins the fix."""
    monkeypatch.setattr(tools, "run_r",
                        lambda script, payload, **kw: {"ok": True,
                                                       "result": load("r_metaboanalyst_map.json")})
    d = str(tmp_path)
    tools.ingest_names(d, names=["Taurine", "LysoPC(16:0)", "NotARealCompound"])
    out = tools.exact_match(d)
    from metabo_idmapper.state import Session
    acc = {e["original_name"]: e.get("accepted", {}) for e in Session(d).entries.values()}
    assert acc["Taurine"]["kegg"] == "C00245"                 # not swallowed as an anchor
    assert acc["LysoPC(16:0)"]["kegg"] == "C04230"
    assert acc["NotARealCompound"] == {}                      # "NA" strings are not ids
    assert out["matched"] == 2
    # "LysoPC(16:0/0:0)" is the same compound in stricter notation (0:0 = the empty position),
    # so the isomer screen must NOT report it as a mismatch — a false alarm here would send
    # the driver re-deciding correct M1 matches.
    assert out["isomer_or_class_mismatch"] == []
    from metabo_idmapper.engine import isomer
    assert isomer.compare("LysoPC(16:0)", "LysoPC(16:0/0:0)")["verdict"] == "ok"
    assert isomer.compare("LysoPC(16:0)", "LysoPC(18:2/0:0)")["verdict"] == "conflict"


def test_bridgedb_normalizes_spellings_and_reports_real_ambiguity(monkeypatch, tmp_path):
    """BridgeDb answers in every spelling it knows (HMDB00251 AND HMDB0000251, ChEBI bare and
    prefixed). Those are one accession; two DIFFERENT accessions are the case that matters."""
    monkeypatch.setattr(tools, "run_r",
                        lambda script, payload, **kw: {"ok": True,
                                                       "result": load("r_bridge_xref.json")})
    d = str(tmp_path)
    tools.ingest_names(d, names=["Taurine", "Taurochenodeoxycholic acid"])
    out = tools.bridge_xref(d, [
        {"feature_id": "M0000", "id": "XOAAWQZATWQOTB-UHFFFAOYSA-N", "source": "inchikey",
         "targets": ["kegg", "hmdb", "chebi"]},
        {"feature_id": "M0001", "id": "C05465", "source": "kegg", "targets": ["hmdb"]}])
    taurine, tcdca = out["results"]
    assert taurine["mappings"]["kegg"] == ["C00245"]
    assert taurine["mappings"]["hmdb"] == ["HMDB0000251"]      # one accession, not two
    assert "hmdb" not in taurine["ambiguous"]
    assert taurine["ambiguous"]["chebi"] == ["CHEBI:507393", "CHEBI:15891"]
    # C05465 bridges to two genuinely different HMDB accessions
    assert sorted(tcdca["ambiguous"]["hmdb"]) == ["HMDB0000949", "HMDB0000951"]

    from metabo_idmapper.state import Session
    cand = next(c for c in Session(d).entries["M0001"]["candidates"]
                if c["source"] == "bridgedb")
    assert cand["ambiguous"]["hmdb"]                          # the ambiguity reaches the ledger


def test_backfill_refuses_an_ambiguous_bridge(monkeypatch, tmp_path):
    """Choosing the first of two distinct accessions would record a coin flip as a fact."""
    from metabo_idmapper import flags
    from metabo_idmapper.state import Session
    d = str(tmp_path)
    tools.ingest_names(d, names=["Taurochenodeoxycholic acid"])
    s = Session(d)
    s.entries["M0000"]["candidates"].append({"source": "kegg-search", "kegg": "C05465"})
    s.save()
    tools.record_decision(d, "M0000", rationale="bile acid", accepted={"kegg": "C05465"})
    monkeypatch.setattr(tools, "run_r",
                        lambda script, payload, **kw: {"ok": True,
                                                       "result": load("r_bridge_xref.json")[1:]})
    out = tools.backfill_hmdb(d)
    assert out["hmdb_gained"] == 0
    assert sorted(out["skipped_conflicts"][0]["candidates"]) == ["HMDB0000949", "HMDB0000951"]
    e = Session(d).entries["M0000"]
    assert e.get("accepted", {}).get("hmdb") is None
    assert "hmdb_backfill_conflict" in flags.open_flags(e, Session(d))


# --------------------------------------------------------------------- REST endpoints
def test_pubchem_lookup_parser(monkeypatch):
    """PUG-REST property + synonyms -> the evidence structure_lookup attaches."""
    payloads = {"property": load("pubchem_property.json"),
                "synonyms": load("pubchem_synonyms.json")}
    monkeypatch.setattr(structure, "get_json",
                        lambda url, timeout=30: payloads["synonyms"] if "synonyms" in url
                        else payloads["property"])
    res = structure.lookup("Taurine")
    assert res["found"] and res["cid"] == 1123
    assert res["formula"] == "C2H7NO3S"
    assert res["inchikey"] == "XOAAWQZATWQOTB-UHFFFAOYSA-N"
    assert abs(res["mono_mass"] - 125.0147) < 0.001            # parsed as a float, not a string
    assert res["synonyms_sample"][0] == "taurine" and len(res["synonyms_sample"]) <= 25


def test_pubchem_entry_parser_returns_the_synonym_list(monkeypatch):
    """The CID back-check leans on synonyms: one entry carrying two compounds' names is how
    two features end up sharing an id."""
    monkeypatch.setattr(entry._struct, "get_json",
                        lambda url, timeout=30: load("pubchem_synonyms.json")
                        if "synonyms" in url else load("pubchem_property.json"))
    rec = entry.pubchem_entry("1123")
    assert rec["found"] and rec["formula"] == "C2H7NO3S"
    assert "2-aminoethanesulfonic acid" in rec["names"]


def test_pubchem_miss_is_not_an_error(monkeypatch):
    """A miss returns found=False with the cues to retry, never a fabricated id."""
    monkeypatch.setattr(structure, "get_json", lambda url, timeout=30: None)
    res = structure.lookup("NotARealCompound")
    assert res["found"] is False and res["cid"] is None


def test_chebi_search_and_term_parsers(monkeypatch):
    """OLS4 search ranks the exact label first; the term lookup carries label + formula."""
    class _Resp:
        status_code = 200

        def __init__(self, payload):
            self._payload = payload

        def json(self):
            return self._payload

    monkeypatch.setattr(chebi.requests, "get",
                        lambda *a, **k: _Resp(load("chebi_search.json")))
    hits = chebi.search("taurine", rows=4)
    assert hits[0]["chebi"] == "CHEBI:15891" and hits[0]["exact"] is True
    assert any(not h["exact"] for h in hits[1:])

    monkeypatch.setattr(entry.requests, "get",
                        lambda *a, **k: _Resp(load("chebi_term.json")))
    rec = entry.chebi_entry("15891")                 # bare id is prefixed by the parser
    assert rec["found"] and rec["names"][0] == "taurine"
    assert rec["formula"] == "C2H7NO3S"


# --------------------------------------------------------------------- the R transport itself
def test_run_r_reads_the_output_file_not_stdout(tmp_path, monkeypatch):
    """MetaboAnalystR and curl both write to stdout, so the contract is that R hands back a
    FILE. This drives the real subprocess path with a stub 'Rscript'."""
    from metabo_idmapper.engine import rcall

    stub = tmp_path / "fake_rscript"
    stub.write_text(
        "#!/bin/sh\n"
        "echo 'loadDatabase OK: noise on stdout'\n"
        "payload=$(cat)\n"
        "out=$(printf '%s' \"$payload\" | sed 's/.*\"out\": *\"\\([^\"]*\\)\".*/\\1/')\n"
        "printf '[{\"kegg\":\"C00245\",\"found\":true}]' > \"$out\"\n")
    stub.chmod(0o755)
    monkeypatch.setattr(rcall, "RSCRIPT", str(stub))
    monkeypatch.setattr(rcall, "RSCRIPT_DIR", Path(__file__).parent.parent
                        / "metabo_idmapper" / "rscripts")
    r = rcall.run_r("kegg_entry.R", {"ids": ["C00245"]})
    assert r["ok"] and r["result"][0]["kegg"] == "C00245"     # stdout noise ignored

    silent = tmp_path / "silent_rscript"
    silent.write_text("#!/bin/sh\ncat > /dev/null\n")         # writes no output file
    silent.chmod(0o755)
    monkeypatch.setattr(rcall, "RSCRIPT", str(silent))
    r = rcall.run_r("kegg_entry.R", {"ids": ["C00245"]})
    assert r["ok"] is False and "no output" in r["error"]

    assert rcall.run_r("no_such_engine.R", {})["ok"] is False


# --------------------------------------------------------------------- shipped reproduction
def test_raw_engine_parsers_match_the_same_fixtures(monkeypatch):
    """The reproduction code ships to users and reads the SAME payloads; its parsers are
    tested against the same recordings so the artifact cannot rot separately."""
    from metabo_idmapper.codegen import raw_engine as raw

    monkeypatch.setattr(raw, "run_r",
                        lambda src, payload, rscript="Rscript", timeout=1200:
                        load("r_metaboanalyst_map.json"))
    ma = raw.metaboanalyst_map(["Taurine"], "<R source>", "/tmp", "Rscript")
    assert ma["Taurine"]["kegg"] == "C00245"

    monkeypatch.setattr(raw, "run_r",
                        lambda src, payload, rscript="Rscript", timeout=1200:
                        load("r_kegg_entry.json"))
    entries = raw.kegg_entries(["C05472"], "<R source>")
    assert entries["C05472"]["names"][0] == "Urocortisol"
    rows = [{"feature_id": "M0", "name": "Taurochenodeoxycholic acid", "kegg": "C05472"},
            {"feature_id": "M1", "name": "LysoPC(16:0)", "kegg": "C04230"}]
    back = raw.apply_kegg_backcheck(rows, entries)
    assert back[0]["kegg_name"] == "Urocortisol"
    assert rows[1]["kegg_resolution"] == "class"     # 1-Acyl-... has no chain spec

    monkeypatch.setattr(raw, "run_r",
                        lambda src, payload, rscript="Rscript", timeout=1200:
                        load("r_kegg_search.json"))
    assert raw.kegg_find(["taurochenodeoxycholate"], "<R>")[0]["kegg"] == "C05465"
