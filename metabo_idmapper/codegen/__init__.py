"""Source that is SHIPPED, not imported: the raw-API implementation the reproduction code is
built from. `raw_engine` is a real, testable module with no metabo_idmapper import, so the
emitted script can inline it verbatim and the notebook can ship it as a sidecar — the same
way the .R engines are shipped. Nothing in the package imports it at runtime."""
