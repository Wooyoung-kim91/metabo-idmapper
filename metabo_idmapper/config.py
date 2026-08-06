"""Runtime configuration — paths to the R engine, bridge DB, and GEM model.

Everything is overridable by environment variable. The convenience defaults point at the
assets this engine was verified against, but they are only used IF THEY EXIST on this host: a
path baked into the package would otherwise travel to another machine as a valid-looking
string that resolves to nothing, and the failure would surface deep inside COBRApy or the R
subprocess instead of at the call that needed it. Absent an asset the value is None, and the
tools say what to pass (`model_path=` / `db=`) or which variable to set.
"""

from __future__ import annotations

import os
from pathlib import Path

# Assets used to verify the engine; kept as a default only while they are present.
_VERIFIED = "/home/wykim/data/Atherosclerosis/Omics"


def _asset(env_var: str, default: str) -> str | None:
    """Env var if set (trusted as given), else the verified default only if it is there."""
    val = os.environ.get(env_var)
    if val:
        return val
    return default if Path(default).exists() else None


# BridgeDb metabolite mapping database (.bridge). Derby-backed; ~2 GB.
BRIDGE_DB = _asset("METABO_IDMAP_BRIDGE_DB",
                   f"{_VERIFIED}/models/bridgedb/metabolites_20210109.bridge")

# Default genome-scale model for the species crosswalk (any COBRA SBML/JSON). Per-call
# `model_path=` wins, and a run records the models it used per label in its ledger.
GEM_MODEL = _asset("METABO_IDMAP_GEM", f"{_VERIFIED}/models/Mouse-GEM/Mouse-GEM.xml")

# GEM annotation keys carrying the cross-references we crosswalk against.
GEM_XREF_KEYS = {
    "kegg": "kegg.compound",
    "hmdb": "hmdb",
    "chebi": "chebi",
}

# Rscript binary (system R with BridgeDbR/KEGGREST/MetaboAnalystR installed).
RSCRIPT = os.environ.get("METABO_IDMAP_RSCRIPT", "Rscript")

# Directory holding the bundled .R engine snippets.
RSCRIPT_DIR = Path(__file__).parent / "rscripts"

# PubChem REST politeness delay (seconds) — the public endpoint rate-limits.
PUBCHEM_DELAY = float(os.environ.get("METABO_IDMAP_PUBCHEM_DELAY", "0.25"))
