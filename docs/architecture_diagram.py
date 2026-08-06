#!/usr/bin/env python3
"""Render the metabo-idmapper architecture diagram with matplotlib only."""

from __future__ import annotations

import os
import re
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metabo-idmapper-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


DOCS_DIR = Path(__file__).resolve().parent
OUTPUT = DOCS_DIR / "metabo_idmapper_architecture.png"
MIN_FONT_SIZE = 6.5

NAVY = "#17324D"
INK = "#17212B"
MUTED = "#4E5A65"
WHITE = "#FFFFFF"
PAPER = "#F8FAFC"
BLUE = "#DCEEFF"
BLUE_EDGE = "#2B6FA6"
AMBER = "#FFF0CF"
AMBER_EDGE = "#D98200"
RED = "#FBE1E1"
RED_EDGE = "#B53B3B"
GRAY = "#EFF2F5"
GRAY_EDGE = "#687580"

STAGE_Y = 0.495
STAGE_HEIGHT = 0.306
SPINE_CENTER_Y = 0.668
WRONG_ID_STATEMENTS = (
    "CATCH THE INVISIBLE ERROR",
    "A missing ID is visible. A WRONG ID arrives wearing the right name.",
)

TOOL_GROUPS = (
    ("ORIENT", ("midmap_guidance", "detect_state")),
    ("INGEST", ("ingest_names",)),
    ("1st PASS", ("exact_match",)),
    (
        "EVIDENCE",
        (
            "structure_lookup",
            "search_synonym",
            "mass_match_candidates",
            "verify_candidate",
            "bridge_xref",
        ),
    ),
    ("VERIFY THE ID", ("id_name_check", "isomer_guard", "collision_check")),
    ("THE CALL", ("record_decision", "screen_exogenous", "acknowledge_flag")),
    ("MODEL", ("gem_crosswalk", "gem_search", "gem_assign", "backfill_hmdb")),
    (
        "REPORT + AUDIT",
        ("coverage_summary", "finalize_run", "acknowledge_check", "harness_audit"),
    ),
)
TOOLS = tuple(tool for _, tools in TOOL_GROUPS for tool in tools)


def box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    face: str = WHITE,
    edge: str = GRAY_EDGE,
    linewidth: float = 1.2,
    radius: float = 0.008,
    zorder: int = 1,
):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0.004,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def label(
    ax,
    x: float,
    y: float,
    text: str,
    *,
    size: float = 7.0,
    color: str = INK,
    weight: str = "normal",
    ha: str = "center",
    va: str = "center",
    family: str = "DejaVu Sans",
    linespacing: float = 1.12,
    zorder: int = 5,
):
    return ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=size,
        color=color,
        fontweight=weight,
        ha=ha,
        va=va,
        family=family,
        linespacing=linespacing,
        zorder=zorder,
    )


def arrow(
    ax,
    start,
    end,
    *,
    color=GRAY_EDGE,
    width=1.3,
    style="-|>",
    connectionstyle="arc3",
    zorder=3,
):
    patch = FancyArrowPatch(
        start,
        end,
        transform=ax.transAxes,
        arrowstyle=style,
        mutation_scale=10,
        linewidth=width,
        color=color,
        connectionstyle=connectionstyle,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def tool_chip(ax, x, y, width, text, *, face=WHITE, edge=BLUE_EDGE, color=NAVY):
    height = 0.027
    box(ax, x, y, width, height, face=face, edge=edge, linewidth=0.8, radius=0.004, zorder=3)
    label(
        ax,
        x + width / 2,
        y + height / 2,
        text,
        size=6.7,
        color=color,
        weight="bold",
        family="DejaVu Sans Mono",
        zorder=4,
    )


def draw_llm_band(ax):
    box(ax, 0.018, 0.842, 0.964, 0.082, face=NAVY, edge=NAVY, linewidth=0, radius=0.010)
    label(ax, 0.036, 0.895, "LLM", size=11.5, color=WHITE, weight="bold", ha="left")
    label(
        ax,
        0.036,
        0.865,
        "makes every identity call",
        size=7.2,
        color="#D8E6F2",
        ha="left",
    )

    prompt_specs = (
        (0.223, 0.122, "DRIVER PROMPT", "map_metabolites\nmetabo-idmapper://driver"),
        (0.517, 0.122, "DECISION", "evidence in  →  call out"),
        (0.761, 0.164, "REVIEWER PROMPT", "review_mappings\nmetabo-idmapper://reviewer"),
    )
    for x, width, heading, body in prompt_specs:
        box(ax, x, 0.855, width, 0.056, face=WHITE, edge="#B8CDDD", linewidth=0.9, radius=0.006, zorder=2)
        label(ax, x + width / 2, 0.894, heading, size=6.7, color=BLUE_EDGE, weight="bold")
        label(ax, x + width / 2, 0.872, body, size=6.7, color=INK, linespacing=1.02)


def draw_stage(ax, x, width, number, title, tools, *, face=BLUE, edge=BLUE_EDGE, dominant=False):
    y = STAGE_Y
    height = STAGE_HEIGHT
    box(ax, x, y, width, height, face=face, edge=edge, linewidth=2.4 if dominant else 1.3, radius=0.009)

    ax.add_patch(
        Rectangle(
            (x, y + height - 0.047),
            width,
            0.047,
            transform=ax.transAxes,
            facecolor=edge,
            edgecolor="none",
            zorder=2,
        )
    )
    label(ax, x + 0.012, y + height - 0.0235, str(number), size=7.0, color=WHITE, weight="bold")
    label(
        ax,
        x + width / 2 + 0.004,
        y + height - 0.0235,
        title,
        size=7.2 if len(title) < 12 else 6.7,
        color=WHITE,
        weight="bold",
    )

    if dominant:
        label(ax, x + width / 2, y + 0.214, WRONG_ID_STATEMENTS[0], size=8.2, color=RED_EDGE, weight="bold")
        chip_start = y + 0.148
        gap = 0.012
    elif len(tools) >= 4:
        chip_start = y + 0.194
        gap = 0.011
    elif len(tools) == 3:
        chip_start = y + 0.180
        gap = 0.012
    else:
        chip_start = y + 0.186
        gap = 0.010

    chip_height = 0.027
    for index, tool in enumerate(tools):
        tool_chip(
            ax,
            x + 0.007,
            chip_start - index * (chip_height + gap),
            width - 0.014,
            tool,
            face="#FFFDFD" if dominant else WHITE,
            edge=edge,
            color=RED_EDGE if dominant else NAVY,
        )

    footer = {
        "ORIENT": "state + guidance",
        "INGEST": "names → entries",
        "1st PASS": "M1 auto-accept",
        "EVIDENCE": "tools emit evidence",
        "THE CALL": "LLM records decision",
        "MODEL": "flux relationship",
        "REPORT + AUDIT": "export + govern",
    }.get(title)
    if footer:
        label(ax, x + width / 2, y + 0.018, footer, size=6.7, color=MUTED, weight="bold")


def draw_spine(ax):
    starts = (0.018, 0.101, 0.174, 0.255, 0.399, 0.569, 0.686, 0.837)
    widths = (0.077, 0.067, 0.075, 0.138, 0.164, 0.111, 0.145, 0.145)
    colors = (
        (GRAY, GRAY_EDGE),
        (BLUE, BLUE_EDGE),
        (BLUE, BLUE_EDGE),
        (AMBER, AMBER_EDGE),
        (RED, RED_EDGE),
        (AMBER, AMBER_EDGE),
        (BLUE, BLUE_EDGE),
        (GRAY, GRAY_EDGE),
    )

    for index, ((title, tools), x, width, (face, edge)) in enumerate(
        zip(TOOL_GROUPS, starts, widths, colors), start=1
    ):
        draw_stage(ax, x, width, index, title, tools, face=face, edge=edge, dominant=(title == "VERIFY THE ID"))
        if index < len(TOOL_GROUPS):
            arrow(
                ax,
                (x + width + 0.001, SPINE_CENTER_Y),
                (starts[index] - 0.001, SPINE_CENTER_Y),
                color=NAVY,
                width=1.6,
                zorder=6,
            )

    label(ax, 0.018, 0.818, "DETERMINISTIC MCP TOOL REGISTRY", size=7.6, color=NAVY, weight="bold", ha="left")

    # LLM delegation and review gate.
    arrow(ax, (0.284, 0.842), (0.324, 0.803), color=AMBER_EDGE, width=1.5)
    arrow(ax, (0.578, 0.842), (0.625, 0.803), color=AMBER_EDGE, width=1.5)
    label(ax, 0.825, 0.816, "review gate", size=6.7, color=GRAY_EDGE, weight="bold", ha="right")
    arrow(
        ax,
        (0.842, 0.842),
        (0.834, 0.682),
        color=GRAY_EDGE,
        width=1.5,
        connectionstyle="angle,angleA=0,angleB=90,rad=0",
        zorder=7,
    )
    ax.add_patch(
        Rectangle(
            (0.8325, 0.656),
            0.003,
            0.024,
            transform=ax.transAxes,
            facecolor=NAVY,
            edgecolor=WHITE,
            linewidth=0.6,
            zorder=8,
        )
    )


def draw_wrong_id_panel(ax):
    box(ax, 0.216, 0.292, 0.568, 0.174, face="#FFF8F8", edge=RED_EDGE, linewidth=2.0, radius=0.010)
    label(
        ax,
        0.500,
        0.441,
        WRONG_ID_STATEMENTS[1],
        size=8.0,
        color=INK,
        weight="bold",
    )
    label(
        ax,
        0.236,
        0.411,
        "For lipids and glycans, formula + mass can be identical — identity lives in structure-aware name axes.",
        size=7.0,
        color=MUTED,
        ha="left",
    )

    cases = (
        ("EXACT CAN BE WRONG", "Taurochenodeoxycholate →\nKEGG C05472 ‘Urocortisol’"),
        ("LINKAGE CHANGES IDENTITY", "P-16:0 vinyl ether ≠\n1-acyl class ID"),
        ("ISOMERS HIDE IN PLAIN SIGHT", "neolacto β1-4 ≠ lacto β1-3\nformula and mass stay blind"),
    )
    case_x = (0.236, 0.419, 0.602)
    for x, (heading, body) in zip(case_x, cases):
        box(ax, x, 0.314, 0.162, 0.075, face=WHITE, edge="#E3B3B3", linewidth=0.9, radius=0.005)
        label(ax, x + 0.081, 0.369, heading, size=6.7, color=RED_EDGE, weight="bold")
        label(ax, x + 0.081, 0.339, body, size=6.7, color=INK, linespacing=1.08)

    arrow(ax, (0.481, STAGE_Y), (0.481, 0.466), color=RED_EDGE, width=2.0)


def fact_panel(ax, x, width, title, kicker, body, *, accent=BLUE_EDGE):
    box(ax, x, 0.100, width, 0.132, face=WHITE, edge="#C8D0D8", linewidth=1.0, radius=0.008)
    ax.add_patch(
        Rectangle((x, 0.205), width, 0.027, transform=ax.transAxes, facecolor=accent, edgecolor="none", zorder=2)
    )
    label(ax, x + 0.012, 0.2185, title, size=7.4, color=WHITE, weight="bold", ha="left")
    label(ax, x + 0.012, 0.181, kicker, size=7.1, color=accent, weight="bold", ha="left")
    label(ax, x + 0.012, 0.158, body, size=6.7, color=INK, ha="left", va="top", linespacing=1.12)


def draw_fact_panels(ax):
    fact_panel(
        ax,
        0.018,
        0.300,
        "LEDGER",
        "midmap_ledger.json = source of truth",
        "revision rejects stale writers\nledger_conflict → nothing written\nschema migration on load\nreview flags derived from state",
        accent=NAVY,
    )
    fact_panel(
        ax,
        0.350,
        0.300,
        "MODEL RELATION",
        "What does the flux number mean?",
        "exact  •  class-proxy\nisomer-surrogate  •  model-scope-absent\nid-gap\nxref found 8/27  •  model held 21/28",
        accent=BLUE_EDGE,
    )
    fact_panel(
        ax,
        0.682,
        0.300,
        "GOVERNANCE",
        "20 read-only audit checks",
        "pass / warn / fail\nwarn ACK needs a reason\nACK lapses when offenders change\nfail can never be acknowledged",
        accent=GRAY_EDGE,
    )

    label(
        ax,
        0.018,
        0.050,
        "CONFIDENCE",
        size=6.7,
        color=NAVY,
        weight="bold",
        ha="left",
    )
    label(
        ax,
        0.094,
        0.050,
        "M1 exact  ·  M2 verified synonym  ·  M3 verified typo  ·  M4 structure-only  ·  W mass-only  ·  X exogenous  ·  U unresolved",
        size=6.7,
        color=MUTED,
        ha="left",
    )

    legend_specs = (
        (0.638, 0.104, BLUE, BLUE_EDGE, "deterministic tool"),
        (0.748, 0.104, AMBER, AMBER_EDGE, "LLM judgment"),
        (0.858, 0.124, GRAY, GRAY_EDGE, "state + output"),
    )
    for x, width, face, edge, text in legend_specs:
        box(ax, x, 0.036, width, 0.028, face=face, edge=edge, linewidth=1.0, radius=0.004, zorder=2)
        label(ax, x + width / 2, 0.050, text, size=6.7, color=NAVY, weight="bold", zorder=3)


def build_figure():
    fig = plt.figure(figsize=(16, 9), facecolor=PAPER)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    label(ax, 0.018, 0.965, "metabo-idmapper", size=18.0, color=NAVY, weight="bold", ha="left")
    label(
        ax,
        0.982,
        0.965,
        "Evidence is deterministic. Identity is a judgment.",
        size=9.0,
        color=MUTED,
        ha="right",
    )
    label(
        ax,
        0.982,
        0.938,
        "tools emit evidence  •  code enforces state",
        size=6.8,
        color=MUTED,
        ha="right",
    )

    draw_llm_band(ax)
    draw_spine(ax)
    draw_wrong_id_panel(ax)
    draw_fact_panels(ax)
    return fig, ax


def _intersection_area(first, second) -> float:
    left = max(first.x0, second.x0)
    right = min(first.x1, second.x1)
    bottom = max(first.y0, second.y0)
    top = min(first.y1, second.y1)
    return max(0.0, right - left) * max(0.0, top - bottom)


def validate_layout(fig, ax) -> dict[str, object]:
    """Validate typography, tool inventory, and renderer-level text collisions."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [artist for artist in ax.texts if artist.get_visible() and artist.get_text().strip()]

    small_text = [(artist.get_text(), artist.get_fontsize()) for artist in texts if artist.get_fontsize() < MIN_FONT_SIZE]

    corpus = "\n".join(artist.get_text() for artist in texts)
    tool_counts = {
        tool: len(re.findall(rf"(?<![A-Za-z0-9_]){re.escape(tool)}(?![A-Za-z0-9_])", corpus))
        for tool in TOOLS
    }
    bad_tools = {tool: count for tool, count in tool_counts.items() if count != 1}
    wrong_id_counts = {statement: corpus.count(statement) for statement in WRONG_ID_STATEMENTS}

    bboxes = [(artist, artist.get_window_extent(renderer=renderer)) for artist in texts]
    overlaps = []
    for (first_artist, first_box), (second_artist, second_box) in combinations(bboxes, 2):
        # A sub-pixel intersection is normally antialiasing at a shared boundary, not
        # a readable glyph collision. Anything larger is a renderer-level overlap.
        area = _intersection_area(first_box, second_box)
        if area > 1.0:
            overlaps.append((first_artist.get_text(), second_artist.get_text(), area))

    errors = []
    if small_text:
        errors.append(f"{len(small_text)} text element(s) below {MIN_FONT_SIZE:.1f} pt")
    if bad_tools:
        errors.append(f"tool occurrence errors: {bad_tools}")
    if any(count != 1 for count in wrong_id_counts.values()):
        errors.append(f"wrong-ID thesis occurrence errors: {wrong_id_counts}")
    if overlaps:
        errors.append(f"{len(overlaps)} text overlap(s)")

    diagnostics = {
        "text_count": len(texts),
        "min_font": min(artist.get_fontsize() for artist in texts),
        "tool_count": sum(count == 1 for count in tool_counts.values()),
        "bad_tools": bad_tools,
        "wrong_id_count": sum(wrong_id_counts.values()),
        "overlaps": overlaps,
        "errors": errors,
    }
    if errors:
        for first, second, area in overlaps[:10]:
            print(f"  overlap ({area:.1f} px²): {first!r} <> {second!r}")
        raise RuntimeError("Layout validation failed: " + "; ".join(errors))
    return diagnostics


def main() -> None:
    fig, ax = build_figure()
    diagnostics = validate_layout(fig, ax)
    fig.savefig(OUTPUT, dpi=200, facecolor=fig.get_facecolor(), bbox_inches=None)
    plt.close(fig)

    print(f"Output: {OUTPUT}")
    print("Canvas: 3200 × 1800 px (16:9, 200 dpi)")
    print(f"Rendered text elements: {diagnostics['text_count']}")
    print(f"Tool-name check: PASS — {diagnostics['tool_count']}/{len(TOOLS)} appear exactly once")
    print(f"Wrong-ID thesis check: PASS — {diagnostics['wrong_id_count']} statements")
    print(
        f"Minimum font check: PASS — {diagnostics['min_font']:.1f} pt "
        f"(required ≥ {MIN_FONT_SIZE:.1f} pt)"
    )
    print(f"Text-overlap check: PASS — {len(diagnostics['overlaps'])} overlaps")
    print("Layout validator: PASS")


if __name__ == "__main__":
    main()
