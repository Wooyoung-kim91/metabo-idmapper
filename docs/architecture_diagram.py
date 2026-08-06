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
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


DOCS_DIR = Path(__file__).resolve().parent
OUTPUT = DOCS_DIR / "metabo_idmapper_architecture.png"
FIGSIZE = (12.0, 6.75)
DPI = 200
MIN_FONT_WIDTH_RATIO = 0.0125

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

WRONG_ID_STATEMENT = "A missing ID is visible. A WRONG ID arrives wearing the right name."


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
    size: float = 11.0,
    color: str = INK,
    weight: str = "normal",
    ha: str = "center",
    va: str = "center",
    family: str = "DejaVu Sans",
    linespacing: float = 1.08,
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


def label_in_patch(
    ax,
    patch,
    x: float,
    y: float,
    text: str,
    *,
    available_width: float,
    padding: float = 0.0,
    **kwargs,
):
    """Draw a label and register its intentional container and usable width."""
    artist = label(ax, x, y, text, **kwargs)
    if not hasattr(ax, "_intentional_text_patch_pairs"):
        ax._intentional_text_patch_pairs = set()
        ax._text_fit_constraints = []
    ax._intentional_text_patch_pairs.add((id(artist), id(patch)))
    ax._text_fit_constraints.append((artist, patch, available_width, padding))
    return artist


def arrow(ax, start, end, *, color=GRAY_EDGE, width=1.5, style="-|>", zorder=6):
    patch = FancyArrowPatch(
        start,
        end,
        transform=ax.transAxes,
        arrowstyle=style,
        mutation_scale=12,
        linewidth=width,
        color=color,
        fill=False,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def path_arrow(ax, vertices, *, color=NAVY, width=1.7, zorder=7):
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(vertices) - 1)
    patch = FancyArrowPatch(
        path=MplPath(vertices, codes),
        transform=ax.transAxes,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=width,
        color=color,
        fill=False,
        shrinkA=0,
        shrinkB=0,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def tool_chip(ax, x, y, width, text, *, face=WHITE, edge=BLUE_EDGE, color=NAVY):
    height = 0.032
    patch = box(ax, x, y, width, height, face=face, edge=edge, linewidth=0.9, radius=0.004, zorder=3)
    label_in_patch(
        ax,
        patch,
        x + width / 2,
        y + height / 2,
        text,
        available_width=width,
        padding=0.008,
        size=11.0,
        color=color,
        weight="bold",
        family="DejaVu Sans Mono",
        zorder=4,
    )


def draw_stage(ax, x, y, width, number, title, tools, *, face, edge, dominant=False):
    height = 0.240 if y > 0.6 else 0.225
    header_height = 0.045
    box(
        ax,
        x,
        y,
        width,
        height,
        face=face,
        edge=edge,
        linewidth=2.4 if dominant else 1.4,
        radius=0.008,
    )
    header_patch = Rectangle(
        (x, y + height - header_height),
        width,
        header_height,
        transform=ax.transAxes,
        facecolor=edge,
        edgecolor="none",
        zorder=2,
    )
    ax.add_patch(header_patch)
    label_in_patch(
        ax,
        header_patch,
        x + 0.018,
        y + height - header_height / 2,
        str(number),
        available_width=0.036,
        size=11.0,
        color=WHITE,
        weight="bold",
    )
    label_in_patch(
        ax,
        header_patch,
        x + width / 2 + 0.006,
        y + height - header_height / 2,
        title,
        available_width=width - 0.025,
        padding=0.004,
        size=14.0,
        color=WHITE,
        weight="bold",
    )

    if y > 0.6:  # phases 1–4
        if len(tools) == 5:
            chip_positions = [y + 0.157 - index * 0.036 for index in range(5)]
        elif len(tools) == 2:
            chip_positions = [y + 0.112, y + 0.060]
        else:
            chip_positions = [y + 0.086]
    elif dominant:
        chip_positions = [y + 0.103, y + 0.061, y + 0.019]
    elif len(tools) == 4:
        chip_positions = [y + 0.140 - index * 0.041 for index in range(4)]
    else:
        chip_positions = [y + 0.118 - index * 0.047 for index in range(len(tools))]

    for chip_y, tool in zip(chip_positions, tools):
        tool_chip(
            ax,
            x + 0.009,
            chip_y,
            width - 0.018,
            tool,
            face="#FFFDFD" if dominant else WHITE,
            edge=edge,
            color=RED_EDGE if dominant else NAVY,
        )


def draw_legend(ax):
    label(ax, 0.018, 0.915, "DETERMINISTIC MCP TOOL REGISTRY", size=11.0, color=NAVY, weight="bold", ha="left")
    legend_specs = (
        (0.520, 0.145, BLUE, BLUE_EDGE, "deterministic tool"),
        (0.675, 0.140, AMBER, AMBER_EDGE, "LLM judgment"),
        (0.825, 0.157, GRAY, GRAY_EDGE, "state + output"),
    )
    for x, width, face, edge, text in legend_specs:
        patch = box(ax, x, 0.900, width, 0.030, face=face, edge=edge, linewidth=1.0, radius=0.004, zorder=2)
        label_in_patch(
            ax,
            patch,
            x + width / 2,
            0.915,
            text,
            available_width=width,
            padding=0.008,
            size=11.0,
            color=NAVY,
            weight="bold",
            zorder=3,
        )


def draw_spine(ax):
    xs = (0.018, 0.262, 0.506, 0.750)
    top_y = 0.650
    bottom_y = 0.295
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

    # Phases 1–4 read left to right.
    for index in range(4):
        title, tools = TOOL_GROUPS[index]
        face, edge = colors[index]
        draw_stage(ax, xs[index], top_y, 0.232, index + 1, title, tools, face=face, edge=edge)
        if index < 3:
            arrow(ax, (xs[index] + 0.234, 0.755), (xs[index + 1] - 0.002, 0.755), color=NAVY)

    # Phase 5 starts the wrapped second row; phases 5–8 read left to right.
    bottom_specs = ((0.018, 0.300), (0.330, 0.216), (0.558, 0.206), (0.776, 0.206))
    for offset, (x, width) in enumerate(bottom_specs):
        group_index = 4 + offset
        title, tools = TOOL_GROUPS[group_index]
        face, edge = colors[group_index]
        draw_stage(
            ax,
            x,
            bottom_y,
            width,
            group_index + 1,
            title,
            tools,
            face=face,
            edge=edge,
            dominant=(group_index == 4),
        )
        if offset < 3:
            next_x, next_width = bottom_specs[offset + 1]
            arrow(ax, (x + width + 0.002, 0.407), (next_x - 0.002, 0.407), color=NAVY)

    path_arrow(
        ax,
        (
            (0.984, 0.755),
            (0.992, 0.755),
            (0.992, 0.637),
            (0.008, 0.637),
            (0.008, 0.407),
            (0.016, 0.407),
        ),
        color=NAVY,
        width=1.8,
    )


def draw_llm_band(ax):
    band_patch = box(ax, 0.018, 0.535, 0.964, 0.090, face=NAVY, edge=NAVY, linewidth=0, radius=0.009)
    label_in_patch(
        ax,
        band_patch,
        0.035,
        0.589,
        "LLM",
        available_width=0.170,
        size=15.0,
        color=WHITE,
        weight="bold",
        ha="left",
    )
    label_in_patch(
        ax,
        band_patch,
        0.035,
        0.553,
        "makes every identity call",
        available_width=0.170,
        size=11.0,
        color="#D8E6F2",
        ha="left",
    )

    prompt_specs = (
        (0.225, 0.185, "DRIVER", "map_metabolites · driver"),
        (0.430, 0.210, "DECISION", "evidence in → call out"),
        (0.660, 0.305, "REVIEWER", "review_mappings · gate"),
    )
    for x, width, heading, body in prompt_specs:
        prompt_patch = box(
            ax,
            x,
            0.548,
            width,
            0.064,
            face=WHITE,
            edge="#B8CDDD",
            linewidth=0.9,
            radius=0.005,
            zorder=2,
        )
        label_in_patch(
            ax,
            prompt_patch,
            x + width / 2,
            0.593,
            heading,
            available_width=width,
            padding=0.010,
            size=12.0,
            color=BLUE_EDGE,
            weight="bold",
        )
        label_in_patch(
            ax,
            prompt_patch,
            x + width / 2,
            0.562,
            body,
            available_width=width,
            padding=0.010,
            size=11.0,
            color=INK,
        )

    # Short, unobstructed links: driver → evidence; decision → call; reviewer → pre-report gate.
    path_arrow(ax, ((0.318, 0.625), (0.318, 0.632), (0.866, 0.632), (0.866, 0.650)), color=AMBER_EDGE, width=1.7)
    arrow(ax, (0.535, 0.535), (0.438, 0.520), color=AMBER_EDGE, width=1.7)
    path_arrow(
        ax,
        ((0.812, 0.535), (0.812, 0.528), (0.770, 0.528), (0.770, 0.421)),
        color=GRAY_EDGE,
        width=1.7,
    )
    ax.add_patch(
        Rectangle(
            (0.7685, 0.395),
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
    panel_patch = box(
        ax,
        0.018,
        0.145,
        0.964,
        0.133,
        face="#FFF8F8",
        edge=RED_EDGE,
        linewidth=2.0,
        radius=0.009,
    )
    label_in_patch(
        ax,
        panel_patch,
        0.500,
        0.255,
        WRONG_ID_STATEMENT,
        available_width=0.964,
        padding=0.030,
        size=15.0,
        color=INK,
        weight="bold",
    )
    label_in_patch(
        ax,
        panel_patch,
        0.500,
        0.224,
        "Formula + mass can be identical; identity lives in structure-aware names.",
        available_width=0.964,
        padding=0.030,
        size=11.0,
        color=MUTED,
    )

    examples = (
        ("EXACT CAN BE WRONG", "C05472 = Urocortisol"),
        ("LINKAGE CHANGES IDENTITY", "P-16:0 ≠ 1-acyl"),
        ("ISOMERS HIDE IN PLAIN SIGHT", "neolacto β1-4 ≠ lacto β1-3"),
    )
    for x, (heading, body) in zip((0.035, 0.351, 0.667), examples):
        example_patch = box(
            ax,
            x,
            0.153,
            0.298,
            0.060,
            face=WHITE,
            edge="#E3B3B3",
            linewidth=0.9,
            radius=0.004,
        )
        label_in_patch(
            ax,
            example_patch,
            x + 0.149,
            0.198,
            heading,
            available_width=0.298,
            padding=0.012,
            size=11.0,
            color=RED_EDGE,
            weight="bold",
        )
        label_in_patch(
            ax,
            example_patch,
            x + 0.149,
            0.166,
            body,
            available_width=0.298,
            padding=0.012,
            size=11.0,
            color=INK,
        )

    arrow(ax, (0.168, 0.295), (0.168, 0.278), color=RED_EDGE, width=2.0)


def fact_panel(ax, x, title, body, *, accent):
    width = 0.300
    panel_patch = box(ax, x, 0.018, width, 0.108, face=WHITE, edge="#C8D0D8", linewidth=1.0, radius=0.007)
    header_patch = Rectangle(
        (x, 0.091),
        width,
        0.035,
        transform=ax.transAxes,
        facecolor=accent,
        edgecolor="none",
        zorder=2,
    )
    ax.add_patch(header_patch)
    label_in_patch(
        ax,
        header_patch,
        x + 0.014,
        0.1085,
        title,
        available_width=width,
        padding=0.014,
        size=13.0,
        color=WHITE,
        weight="bold",
        ha="left",
    )
    label_in_patch(
        ax,
        panel_patch,
        x + 0.014,
        0.057,
        body,
        available_width=width,
        padding=0.014,
        size=11.0,
        color=INK,
        ha="left",
        linespacing=1.12,
    )


def draw_fact_panels(ax):
    fact_panel(
        ax,
        0.018,
        "LEDGER",
        "revision conflict → no write\nmigration on load · flags derived",
        accent=NAVY,
    )
    fact_panel(
        ax,
        0.350,
        "MODEL RELATION",
        "exact · class-proxy · isomer-surrogate\nmodel-scope-absent · id-gap",
        accent=BLUE_EDGE,
    )
    fact_panel(
        ax,
        0.682,
        "GOVERNANCE",
        "20 checks: pass / warn / fail\nwarn ACK expires · fail never ACK",
        accent=GRAY_EDGE,
    )


def build_figure():
    fig = plt.figure(figsize=FIGSIZE, facecolor=PAPER)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    label(ax, 0.018, 0.968, "metabo-idmapper", size=24.0, color=NAVY, weight="bold", ha="left")
    label(
        ax,
        0.982,
        0.968,
        "Evidence is deterministic. Identity is a judgment.",
        size=13.0,
        color=MUTED,
        ha="right",
    )

    draw_legend(ax)
    draw_spine(ax)
    draw_llm_band(ax)
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
    """Validate relative typography, tool inventory, and rendered text collisions."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [artist for artist in ax.texts if artist.get_visible() and artist.get_text().strip()]

    figure_width_points = fig.get_size_inches()[0] * 72.0
    font_ratios = {artist: artist.get_fontsize() / figure_width_points for artist in texts}
    smallest_artist = min(texts, key=font_ratios.get)
    small_text = [
        (artist.get_text(), artist.get_fontsize(), ratio)
        for artist, ratio in font_ratios.items()
        if ratio < MIN_FONT_WIDTH_RATIO
    ]

    corpus = "\n".join(artist.get_text() for artist in texts)
    tool_counts = {
        tool: len(re.findall(rf"(?<![A-Za-z0-9_]){re.escape(tool)}(?![A-Za-z0-9_])", corpus))
        for tool in TOOLS
    }
    bad_tools = {tool: count for tool, count in tool_counts.items() if count != 1}
    wrong_id_count = corpus.count(WRONG_ID_STATEMENT)

    bboxes = [(artist, artist.get_window_extent(renderer=renderer)) for artist in texts]
    overlaps = []
    for (first_artist, first_box), (second_artist, second_box) in combinations(bboxes, 2):
        area = _intersection_area(first_box, second_box)
        if area > 1.0:
            overlaps.append((first_artist.get_text(), second_artist.get_text(), area))

    intentional_pairs = getattr(ax, "_intentional_text_patch_pairs", set())
    patch_clips = []
    for artist, text_box in bboxes:
        for patch in ax.patches:
            if not patch.get_visible() or patch.get_zorder() < artist.get_zorder():
                continue
            if (id(artist), id(patch)) in intentional_pairs:
                continue
            facecolor = patch.get_facecolor()
            if not patch.get_fill() or len(facecolor) < 4 or facecolor[3] == 0:
                continue
            rendered_path = patch.get_path().transformed(patch.get_transform())
            if rendered_path.intersects_bbox(text_box, filled=True):
                patch_clips.append((artist.get_text(), type(patch).__name__, patch.get_zorder()))

    fit_violations = []
    for artist, _patch, available_width, padding in getattr(ax, "_text_fit_constraints", []):
        text_width_px = artist.get_window_extent(renderer=renderer).width
        available_width_px = (available_width - 2 * padding) * ax.bbox.width
        if text_width_px > available_width_px + 1.0:
            fit_violations.append((artist.get_text(), text_width_px, available_width_px))

    errors = []
    if small_text:
        errors.append(f"{len(small_text)} text element(s) below width ratio {MIN_FONT_WIDTH_RATIO:.4f}")
    if bad_tools:
        errors.append(f"tool occurrence errors: {bad_tools}")
    if wrong_id_count != 1:
        errors.append(f"wrong-ID thesis occurrence error: {wrong_id_count}")
    if overlaps:
        errors.append(f"{len(overlaps)} text overlap(s)")
    if patch_clips:
        errors.append(f"{len(patch_clips)} text/patch overlap(s)")
    if fit_violations:
        errors.append(f"{len(fit_violations)} text-fit violation(s)")

    diagnostics = {
        "text_count": len(texts),
        "min_font": smallest_artist.get_fontsize(),
        "min_ratio": font_ratios[smallest_artist],
        "tool_count": sum(count == 1 for count in tool_counts.values()),
        "wrong_id_count": wrong_id_count,
        "overlaps": overlaps,
        "patch_clips": patch_clips,
        "fit_violations": fit_violations,
    }
    if errors:
        for text, size, ratio in small_text[:10]:
            print(f"  subminimum ({ratio:.5f}, {size:.1f} pt): {text!r}")
        for first, second, area in overlaps[:10]:
            print(f"  overlap ({area:.1f} px²): {first!r} <> {second!r}")
        for text, patch_name, patch_zorder in patch_clips[:10]:
            print(f"  patch clip (patch={patch_name}, zorder={patch_zorder}): {text!r}")
        for text, text_width, available_width in fit_violations[:10]:
            print(
                f"  text does not fit: {text!r} "
                f"({text_width:.1f}px text vs {available_width:.1f}px available)"
            )
        raise RuntimeError("Layout validation failed: " + "; ".join(errors))
    return diagnostics


def main() -> None:
    fig, ax = build_figure()
    diagnostics = validate_layout(fig, ax)
    fig.savefig(OUTPUT, dpi=DPI, facecolor=fig.get_facecolor(), bbox_inches=None)
    plt.close(fig)

    pixel_width = round(FIGSIZE[0] * DPI)
    pixel_height = round(FIGSIZE[1] * DPI)
    print(f"Output: {OUTPUT}")
    print(f"Canvas: {pixel_width} × {pixel_height} px (16:9, {DPI} dpi)")
    print(f"Rendered text elements: {diagnostics['text_count']}")
    print(f"Tool-name check: PASS — {diagnostics['tool_count']}/{len(TOOLS)} appear exactly once")
    print(f"Wrong-ID thesis check: PASS — {diagnostics['wrong_id_count']} statement")
    print(
        "Relative font check: PASS — "
        f"minimum ratio {diagnostics['min_ratio']:.5f} ({diagnostics['min_ratio'] * 100:.3f}%) "
        f"at {diagnostics['min_font']:.1f} pt; required ≥ {MIN_FONT_WIDTH_RATIO:.4f} (1.250%)"
    )
    print(f"Text-overlap check: PASS — {len(diagnostics['overlaps'])} overlaps")
    print(f"Patch-clip check: PASS — {len(diagnostics['patch_clips'])} text/patch overlaps found")
    print(f"Text-fits-in-box check: PASS — {len(diagnostics['fit_violations'])} violations found")
    print("Layout validator: PASS")


if __name__ == "__main__":
    main()
