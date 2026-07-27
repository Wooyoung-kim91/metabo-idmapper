"""FastMCP (stdio) server exposing the metabo-idmapper tool registry.

stdout carries ONLY MCP protocol JSON; all logging goes to stderr. Tool functions come
straight from `tools.REGISTRY` — their signatures/docstrings become the MCP schema.
"""

from __future__ import annotations

import logging
import os
import sys

from mcp.server.fastmcp import FastMCP

from . import roles, tools


def _log() -> logging.Logger:
    logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                        format="%(asctime)s [metabo-idmapper] %(message)s")
    return logging.getLogger("metabo-idmapper")


def _register_reasoning_layer(mcp: FastMCP) -> None:
    """Embed the reasoning-layer roles IN the MCP as prompts + resources, so the whole
    driver+reviewer layer ships with the server (no external subagent files)."""

    def map_metabolites(workdir: str = "") -> str:
        """Driver role: drive the metabo-idmapper tools as the reasoning layer to map a
        metabolite-name list to input-usable IDs + Mouse-GEM MAM. Optional workdir."""
        hdr = f"Session workdir: {workdir}\n\n" if workdir else ""
        return hdr + roles.DRIVER

    def review_mappings(workdir: str = "") -> str:
        """Reviewer role: independently, adversarially verify accepted metabolite
        identities + exclusions before finalizing. Optional workdir to read the ledger."""
        hdr = f"Review the ledger in workdir: {workdir}\n\n" if workdir else ""
        return hdr + roles.REVIEWER

    mcp.prompt(name="map_metabolites",
               description="Reasoning-layer DRIVER: map metabolite names -> IDs + GEM MAM.")(map_metabolites)
    mcp.prompt(name="review_mappings",
               description="Reasoning-layer REVIEWER: adversarially verify accepted identities.")(review_mappings)

    @mcp.resource("metabo-idmapper://driver")
    def _driver_res() -> str:
        return roles.DRIVER

    @mcp.resource("metabo-idmapper://reviewer")
    def _reviewer_res() -> str:
        return roles.REVIEWER


def build_server() -> FastMCP:
    lg = _log()
    mcp = FastMCP("metabo-idmapper")

    enable = {n.strip() for n in os.environ.get("MIDMAP_ENABLE_TOOLS", "").split(",") if n.strip()}
    disable = {n.strip() for n in os.environ.get("MIDMAP_DISABLE_TOOLS", "").split(",") if n.strip()}

    for fn in tools.REGISTRY:
        name = fn.__name__
        if enable and name not in enable:
            continue
        if name in disable:
            continue
        mcp.add_tool(fn, name=name, description=(fn.__doc__ or "").strip())

    _register_reasoning_layer(mcp)
    lg.info("registered %d tools + 2 reasoning-layer prompts",
            len(mcp._tool_manager.list_tools()))
    return mcp


def main() -> None:
    build_server().run()  # stdio transport


if __name__ == "__main__":
    main()
