"""Runtime configuration — paths to the R engine, bridge DB, and GEM model.

All are overridable by environment variable so the package is not hardwired to one
project. Defaults point at the Atherosclerosis assets that were used to verify the
engine, so the server is usable out of the box on this host.
"""

from __future__ import annotations

import os
from pathlib import Path

_ATHERO = "/home/wykim/data/Atherosclerosis/Omics"

# BridgeDb metabolite mapping database (.bridge). Derby-backed; ~2 GB.
BRIDGE_DB = os.environ.get(
    "METABO_IDMAP_BRIDGE_DB", f"{_ATHERO}/models/bridgedb/metabolites_20210109.bridge"
)

# Mouse-GEM (or any COBRA SBML) used for the MAM crosswalk. Overridable per model.
GEM_MODEL = os.environ.get("METABO_IDMAP_GEM", f"{_ATHERO}/models/Mouse-GEM/Mouse-GEM.xml")

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
