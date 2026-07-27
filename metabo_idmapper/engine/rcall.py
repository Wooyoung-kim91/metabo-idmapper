"""Invoke bundled R engine snippets (BridgeDb / KEGGREST / MetaboAnalystR).

Each snippet reads a JSON request on stdin and writes a JSON response on stdout; all
diagnostics go to stderr. We keep the R side thin and stateless so failures are
isolated and logged rather than crashing the server.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

from ..config import RSCRIPT, RSCRIPT_DIR


def run_r(script: str, payload: dict, timeout: int = 900) -> dict:
    """Run rscripts/<script>, passing `payload` (with an added `out` file path) as JSON
    stdin. The R side writes its JSON result to `out`; we read the FILE, never stdout —
    MetaboAnalystR/curl pollute stdout, so stdout parsing is unreliable.

    Returns {"ok": True, "result": <parsed>} or {"ok": False, "error": ..., "stderr": ...}.
    Never raises for a tool-level failure — the caller surfaces it to the LLM.
    """
    path = RSCRIPT_DIR / script
    if not path.exists():
        return {"ok": False, "error": f"missing R script: {path}"}
    out_file = Path(tempfile.mkstemp(suffix=".json", prefix="midmap_r_")[1])
    payload = {**payload, "out": str(out_file)}
    try:
        proc = subprocess.run(
            [RSCRIPT, "--vanilla", str(path)],
            input=json.dumps(payload),
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        out_file.unlink(missing_ok=True)
        return {"ok": False, "error": f"R timeout after {timeout}s ({script})"}
    try:
        if proc.returncode != 0:
            return {"ok": False, "error": f"R exit {proc.returncode}",
                    "stderr": proc.stderr[-2000:]}
        if not out_file.exists() or out_file.stat().st_size == 0:
            return {"ok": False, "error": "R produced no output file",
                    "stderr": proc.stderr[-2000:]}
        return {"ok": True, "result": json.loads(out_file.read_text())}
    except (ValueError, json.JSONDecodeError) as e:
        return {"ok": False, "error": f"unparseable R output: {e}",
                "stderr": proc.stderr[-1000:]}
    finally:
        out_file.unlink(missing_ok=True)
