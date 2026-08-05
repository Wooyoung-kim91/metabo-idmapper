"""metabo-idmapper — metabolite name -> input-usable ID (MCP tool registry).

A DETERMINISTIC registry: every tool emits EVIDENCE (candidate IDs, xref bridges,
formula/mass verification). The connecting LLM is the reasoning layer that makes the
identity CALL (which candidate is correct, endogenous vs xenobiotic, confidence tier).
Never hardcode the biological call; never fabricate an ID that a tool did not return.

Engine reuses the verified `metabolite-id-harmonization` skill logic (BridgeDb / KEGGREST
/ MetaboAnalystR / PubChem / molmass / COBRApy). See guidance.CANONICAL for the workflow.
"""

__version__ = "0.2.0"
