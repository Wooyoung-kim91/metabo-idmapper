"""Structured tool failures — a code the reasoning layer can act on, not just a sentence.

Every tool here returns data, so a failure is a return value too. Until now each one wrote its
own `{"error": "..."}` string, which reads fine for a human and tells a model nothing: it
cannot distinguish "you skipped a prerequisite" (run it and retry) from "that identifier does
not exist" (stop and pick another) from "the R engine died" (an environment problem, retrying
will not help). scpilot — the sibling registry this project mirrors — solves it with an error
code plus the tools to run next, and this file is that vocabulary for metabo-idmapper.

    {"error": "<what happened>", "error_code": "<one of CODES>",
     "suggested_next_tools": [...], ...context}

Codes, and what the reasoning layer should DO with each:

  invalid_state     a prerequisite has not run (empty ledger, no candidates yet). Run the
                    suggested tool, then repeat this call.
  unknown_entry     that feature_id is not in this ledger. Check `detect_state`; do not invent.
  invalid_argument  the arguments cannot be satisfied as given (missing/contradictory/not in
                    the vocabulary). Fix the call — retrying it unchanged will fail identically.
  not_backed        an anti-fabrication refusal: an id or a model species that no tool
                    produced. Gather the evidence first; never hand-write the value.
  engine_failed     an external engine failed (Rscript / MetaboAnalystR / BridgeDb / a REST
                    endpoint). The run can continue by another route; this is not a data
                    verdict, so do NOT record it as "no id exists".
  ledger_conflict   another call wrote this workdir first and nothing was written here.
                    Re-read the state and repeat the call.
  export_failed     an artifact could not be written (plot, code, deck). The mapping is
                    unaffected — the ledger and its tables stand.
"""

from __future__ import annotations

CODES = {
    "invalid_state": "a prerequisite tool has not run yet",
    "unknown_entry": "no such feature_id in this ledger",
    "invalid_argument": "arguments missing, contradictory, or outside the vocabulary",
    "not_backed": "refused: no tool produced this id / species (anti-fabrication)",
    "engine_failed": "an external engine (R, REST) failed — not a data verdict",
    "ledger_conflict": "another call wrote this workdir first; nothing was written",
    "export_failed": "an artifact could not be written; the mapping is unaffected",
}


def fail(code: str, message: str, next_tools: list[str] | None = None, **context) -> dict:
    """Build a tool-level failure result. `code` must be one of CODES."""
    assert code in CODES, f"unknown error code {code!r}"
    out = {"error": message, "error_code": code, "means": CODES[code]}
    if next_tools:
        out["suggested_next_tools"] = next_tools
    out.update(context)
    return out


def unknown_entry(feature_id: str) -> dict:
    return fail("unknown_entry", f"unknown feature_id {feature_id}",
                ["detect_state", "ingest_names"])
