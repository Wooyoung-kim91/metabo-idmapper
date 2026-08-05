#!/usr/bin/env python3
"""Render the metabo-idmapper MCP architecture using matplotlib only."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/metabo-idmapper-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "metabo_idmapper_architecture.png"

NAVY = "#17324D"
BLUE = "#DCEEFF"
BLUE_EDGE = "#2B6FA6"
BLUE_DARK = "#0F5485"
AMBER = "#FFF0CF"
AMBER_EDGE = "#D98200"
AMBER_DARK = "#9B5400"
GRAY = "#EFF2F5"
GRAY_EDGE = "#687580"
GREEN = "#DDF4E4"
GREEN_EDGE = "#26834A"
RED = "#FBE1E1"
RED_EDGE = "#B53B3B"
INK = "#17212B"
MUTED = "#4E5A65"
WHITE = "#FFFFFF"


def rounded_box(
    ax,
    x,
    y,
    w,
    h,
    *,
    facecolor,
    edgecolor,
    linewidth=1.3,
    radius=0.012,
    zorder=2,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        transform=ax.transAxes,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def add_text(
    ax,
    x,
    y,
    text,
    *,
    size=7,
    color=INK,
    weight="normal",
    ha="center",
    va="center",
    zorder=4,
    linespacing=1.12,
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
        zorder=zorder,
        linespacing=linespacing,
        family="DejaVu Sans",
    )


def stage_box(
    ax,
    x,
    y,
    w,
    h,
    number,
    title,
    body,
    *,
    facecolor=BLUE,
    edgecolor=BLUE_EDGE,
    title_color=BLUE_DARK,
    body_size=6.35,
    title_size=8.2,
    footer=None,
    footer_face=WHITE,
    footer_edge=GRAY_EDGE,
    footer_color=INK,
):
    rounded_box(
        ax,
        x,
        y,
        w,
        h,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=2.2 if edgecolor == AMBER_EDGE else 1.4,
    )
    badge_r = min(0.014, w * 0.13)
    badge = plt.Circle(
        (x + 0.014, y + h - 0.022),
        badge_r,
        transform=ax.transAxes,
        facecolor=edgecolor,
        edgecolor=WHITE,
        linewidth=0.8,
        zorder=5,
    )
    ax.add_patch(badge)
    add_text(
        ax,
        x + 0.014,
        y + h - 0.022,
        str(number),
        size=6.6,
        color=WHITE,
        weight="bold",
        zorder=6,
    )
    add_text(
        ax,
        x + w / 2,
        y + h - 0.047,
        title,
        size=title_size,
        color=title_color,
        weight="bold",
        va="top",
    )
    # With a footer the body must clear the footer box AND the rounding pad the patch adds
    # around it (0.010 offset + 0.040 height + 0.006 pad), or the last body line is drawn
    # underneath it.
    body_bottom = y + (0.062 if footer else 0.018)
    body_top = y + h - 0.082
    add_text(
        ax,
        x + w / 2,
        (body_bottom + body_top) / 2,
        body,
        size=body_size,
        color=INK,
        va="center",
        linespacing=1.15,
    )
    if footer:
        fh = 0.040
        rounded_box(
            ax,
            x + 0.008,
            y + 0.010,
            w - 0.016,
            fh,
            facecolor=footer_face,
            edgecolor=footer_edge,
            linewidth=1.0,
            radius=0.005,
            zorder=4,
        )
        add_text(
            ax,
            x + w / 2,
            y + 0.010 + fh / 2,
            footer,
            size=5.65,
            color=footer_color,
            weight="bold",
            zorder=5,
        )


def arrow(
    ax,
    start,
    end,
    *,
    color=GRAY_EDGE,
    linewidth=1.4,
    style="-|>",
    linestyle="-",
    connectionstyle="arc3",
    mutation_scale=10,
    zorder=3,
):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=mutation_scale,
        linewidth=linewidth,
        linestyle=linestyle,
        color=color,
        connectionstyle=connectionstyle,
        transform=ax.transAxes,
        zorder=zorder,
        shrinkA=0,
        shrinkB=0,
    )
    ax.add_patch(patch)
    return patch


def lane(ax, y, h, label, description, color):
    rounded_box(
        ax,
        0.018,
        y,
        0.964,
        h,
        facecolor=WHITE,
        edgecolor="#CAD1D8",
        linewidth=1.0,
        radius=0.010,
        zorder=0,
    )
    ax.add_patch(
        Rectangle(
            (0.018, y + h - 0.036),
            0.964,
            0.036,
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
            zorder=1,
        )
    )
    add_text(
        ax,
        0.029,
        y + h - 0.018,
        label,
        size=7.6,
        color=NAVY,
        weight="bold",
        ha="left",
        zorder=2,
    )
    add_text(
        ax,
        0.096,
        y + h - 0.018,
        description,
        size=7.15,
        color=NAVY,
        ha="left",
        zorder=2,
    )


def evidence_group(ax, x, y, w, h):
    rounded_box(
        ax,
        x,
        y,
        w,
        h,
        facecolor=AMBER,
        edgecolor=AMBER_EDGE,
        linewidth=2.2,
        radius=0.012,
    )
    badge = plt.Circle(
        (x + 0.014, y + h - 0.022),
        0.014,
        transform=ax.transAxes,
        facecolor=AMBER_EDGE,
        edgecolor=WHITE,
        linewidth=0.8,
        zorder=5,
    )
    ax.add_patch(badge)
    add_text(ax, x + 0.014, y + h - 0.022, "3", size=6.6, color=WHITE, weight="bold", zorder=6)
    add_text(
        ax,
        x + w / 2,
        y + h - 0.043,
        "EVIDENCE & RECOVERY",
        size=8.2,
        color=AMBER_DARK,
        weight="bold",
        va="top",
    )

    inner_x = x + 0.009
    inner_w = w - 0.018
    rounded_box(
        ax,
        inner_x,
        y + 0.073,
        inner_w,
        h - 0.143,
        facecolor=BLUE,
        edgecolor=BLUE_EDGE,
        linewidth=1.0,
        radius=0.006,
        zorder=3,
    )
    evidence = (
        "structure_lookup\n"
        "PubChem -> CID / InChIKey / formula / mass\n\n"
        "search_synonym\n"
        "KEGG keggFind / PubChem / ChEBI OLS4\n"
        "each hit carries an isomer verdict\n\n"
        "mass_match_candidates\n"
        "m/z + adduct -> formula candidates (weak)"
    )
    add_text(
        ax,
        x + w / 2,
        y + 0.073 + (h - 0.143) / 2,
        evidence,
        size=5.15,
        color=INK,
        weight="normal",
        linespacing=1.08,
    )
    rounded_box(
        ax,
        inner_x,
        y + 0.012,
        inner_w,
        0.050,
        facecolor=WHITE,
        edgecolor=BLUE_EDGE,
        linewidth=1.2,
        radius=0.005,
        zorder=3,
    )
    add_text(
        ax,
        x + w / 2,
        y + 0.037,
        "verify_candidate: formula / mass REQUIRED\n"
        "(blind to isomers — see 4b)",
        size=5.3,
        color=BLUE_DARK,
        weight="bold",
        zorder=4,
    )
    add_text(
        ax,
        x + w / 2,
        y + h - 0.068,
        "JUDGMENT: LLM proposes recovery queries",
        size=5.8,
        color=AMBER_DARK,
        weight="bold",
    )


def build_diagram():
    fig, ax = plt.subplots(figsize=(19, 11), dpi=200, facecolor=WHITE)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_text(
        ax,
        0.5,
        0.972,
        "metabo-idmapper MCP — ID mapping pipeline & embedded reasoning-layer agents",
        size=18,
        color=NAVY,
        weight="bold",
    )
    add_text(
        ax,
        0.5,
        0.944,
        "Evidence is deterministic; identity, confidence, inclusion, origin, and model resolution are LLM calls. "
        "25 tools.",
        size=9.4,
        color=MUTED,
    )

    lane(
        ax,
        0.700,
        0.220,
        "LANE A",
        "LLM Reasoning Layer — agents embedded IN the MCP as prompts + resources; connect the server and they ship with it",
        "#FFF5DF",
    )
    lane(
        ax,
        0.335,
        0.340,
        "LANE B",
        "Deterministic tool pipeline — each tool EMITS evidence and returns a small JSON summary; the LLM makes the identity CALL",
        "#EAF4FC",
    )
    lane(
        ax,
        0.105,
        0.205,
        "LANE C",
        "State / provenance & outputs",
        "#F0F2F4",
    )

    # Embedded reasoning-layer agents.
    rounded_box(
        ax,
        0.190,
        0.735,
        0.330,
        0.135,
        facecolor=AMBER,
        edgecolor=AMBER_EDGE,
        linewidth=1.8,
        radius=0.014,
    )
    add_text(ax, 0.208, 0.842, "DRIVER agent", size=11, color=AMBER_DARK, weight="bold", ha="left")
    add_text(
        ax,
        0.208,
        0.812,
        "MCP prompt  map_metabolites\nresource  metabo-idmapper://driver",
        size=7.3,
        color=INK,
        ha="left",
        linespacing=1.25,
    )
    add_text(
        ax,
        0.208,
        0.762,
        "Drives stages 1-7 • judgment concentrated at stages 3 & 5",
        size=6.8,
        color=AMBER_DARK,
        weight="bold",
        ha="left",
    )

    rounded_box(
        ax,
        0.595,
        0.735,
        0.360,
        0.135,
        facecolor=AMBER,
        edgecolor=AMBER_EDGE,
        linewidth=1.8,
        radius=0.014,
    )
    add_text(ax, 0.613, 0.842, "REVIEWER agent", size=11, color=AMBER_DARK, weight="bold", ha="left")
    add_text(
        ax,
        0.613,
        0.812,
        "MCP prompt  review_mappings\nresource  metabo-idmapper://reviewer",
        size=7.3,
        color=INK,
        ha="left",
        linespacing=1.25,
    )
    add_text(
        ax,
        0.613,
        0.762,
        "Independent adversarial gate: confirmed / suspect / refuted",
        size=6.8,
        color=AMBER_DARK,
        weight="bold",
        ha="left",
    )

    # Pipeline geometry. Stage numbering follows guidance.CANONICAL, so the verification gate
    # is 4b (a sub-stage of bridging) rather than renumbering everything downstream.
    stage_y = 0.385
    stage_h = 0.238
    gap = 0.006
    widths = [0.076, 0.078, 0.150, 0.064, 0.124, 0.133, 0.122, 0.098, 0.056]
    xs = [0.026]
    for w in widths[:-1]:
        xs.append(xs[-1] + w + gap)

    stage_box(
        ax,
        xs[0],
        stage_y,
        widths[0],
        stage_h,
        1,
        "ingest_names",
        "Normalize names\n\nFlag:\nparenthetical-\nabbrev\ncombined-name\ndigit-locant\n(isomer-\nsensitive)\n\nXenobiotic\nlexicon cue",
        body_size=5.4,
        title_size=7.6,
    )
    stage_box(
        ax,
        xs[1],
        stage_y,
        widths[1],
        stage_h,
        2,
        "exact_match",
        "MetaboAnalyst\nbatch name -> ID\n\nExact hit:\nauto-accept M1\n\nFlagged: messy /\ntrade-name, and\nMatch name that\ndisagrees on a\ndiscriminant axis",
        body_size=5.35,
        title_size=7.6,
    )
    evidence_group(ax, xs[2], stage_y, widths[2], stage_h)
    stage_box(
        ax,
        xs[3],
        stage_y,
        widths[3],
        stage_h,
        4,
        "bridge_xref",
        "BridgeDb\n\nPromote one\naccepted ID\nto: KEGG /\nHMDB / ChEBI /\nPubChem /\nInChIKey\n\nA bridged ID\nis a LEAD,\nnot an\nidentity",
        body_size=5.2,
        title_size=7.2,
    )
    stage_box(
        ax,
        xs[4],
        stage_y,
        widths[4],
        stage_h,
        "4b",
        "VERIFY the ID",
        "id_name_check\nID -> its OWN DB record\n(keggGet / OLS4 / CID)\n\n"
        "isomer_guard   name-axis compare\nether P- / O- / acyl\nneolacto vs lacto • a2,3 vs a2,6\n"
        "chains • double-bond positions\nHETE vs HPETE • backbone\nspecies- vs CLASS-level\n\n"
        "collision_check\none ID on two compounds",
        facecolor=BLUE,
        edgecolor=RED_EDGE,
        title_color=RED_EDGE,
        body_size=4.45,
        title_size=7.8,
        footer="a wrong ID comes with\nthe RIGHT name — mass\nand formula are blind",
        footer_face=RED,
        footer_edge=RED_EDGE,
        footer_color=RED_EDGE,
    )
    stage_box(
        ax,
        xs[5],
        stage_y,
        widths[5],
        stage_h,
        5,
        "record_decision",
        "THE CALL — only tool that sets final_class\n\n"
        "ID coverage:\nKEGG-mapped • HMDB-mapped\nstructure-only • unmapped\n\n"
        "Origin:\nendogenous\n"
        "exogenous = KEPT + tagged\n"
        "diet / drug / microbial / plant\n"
        "xenobiotic = EXCLUDED\n"
        "contaminant / industrial / additive /\n"
        "surfactant / plasticizer / reagent\n\n"
        "screen_exogenous -> non-biological cues",
        facecolor=AMBER,
        edgecolor=AMBER_EDGE,
        title_color=AMBER_DARK,
        body_size=4.95,
    )
    stage_box(
        ax,
        xs[6],
        stage_y,
        widths[6],
        stage_h,
        6,
        "gem_crosswalk",
        "+ gem_search + gem_assign\n+ backfill_hmdb\n\n"
        "xref route MISSES most lipids:\nGEMs annotate them sparsely\n"
        "(verified: 8/27 by xref,\n21/28 present in the model)\n\n"
        "so every miss is searched by the\nmodel's OWN name / formula / mass\n"
        "+ lipid-shorthand expansion\n\n"
        "relation is the CALL:\nexact • class-proxy (R-group pool)\nisomer-surrogate • model-scope-absent",
        body_size=4.7,
        title_size=7.6,
    )
    stage_box(
        ax,
        xs[7],
        stage_y,
        widths[7],
        stage_h,
        7,
        "coverage_summary",
        "Master ledger +\ncoverage + provenance\n+ gem_curation\n\nFigures: UpSet +\nimprovement\n\nexport_code:\nraw-API reproduction\n\nannotate_source",
        body_size=5.35,
        title_size=7.6,
    )
    stage_box(
        ax,
        xs[8],
        stage_y,
        widths[8],
        stage_h,
        8,
        "harness_audit",
        "READ-ONLY\n\nRun LAST\n\n20 checks\npass / warn /\nfail\n\nGovernance\nscorecard\n\nChanges\nnothing",
        body_size=5.2,
        title_size=6.9,
    )

    # Straight pipeline arrows live only in the inter-stage gaps.
    for i in range(8):
        arrow(
            ax,
            (xs[i] + widths[i] + 0.001, stage_y + stage_h / 2),
            (xs[i + 1] - 0.001, stage_y + stage_h / 2),
            color=BLUE_EDGE,
            linewidth=1.5,
            mutation_scale=9,
            zorder=6,
        )

    # Driver judgment insertion points and review gate.
    arrow(
        ax,
        (0.310, 0.735),
        (xs[2] + widths[2] * 0.38, stage_y + stage_h + 0.002),
        color=AMBER_EDGE,
        linewidth=1.7,
        linestyle="--",
        connectionstyle="arc3,rad=0.08",
        mutation_scale=11,
        zorder=5,
    )
    arrow(
        ax,
        (0.430, 0.735),
        (xs[5] + widths[5] * 0.45, stage_y + stage_h + 0.002),
        color=AMBER_EDGE,
        linewidth=1.7,
        linestyle="--",
        connectionstyle="arc3,rad=-0.08",
        mutation_scale=11,
        zorder=5,
    )
    gate_x = xs[7] - gap / 2
    arrow(
        ax,
        (0.775, 0.735),
        (gate_x, stage_y + stage_h + 0.004),
        color=AMBER_EDGE,
        linewidth=1.8,
        linestyle="--",
        connectionstyle="arc3,rad=0.11",
        mutation_scale=11,
        zorder=5,
    )
    arrow(
        ax,
        (gate_x, stage_y + stage_h + 0.002),
        (gate_x, stage_y + stage_h / 2),
        color=AMBER_EDGE,
        linewidth=1.5,
        linestyle="--",
        connectionstyle="arc3",
        mutation_scale=9,
        zorder=7,
    )
    gate_marker = plt.Circle(
        (gate_x, stage_y + stage_h / 2),
        0.0045,
        transform=ax.transAxes,
        facecolor=WHITE,
        edgecolor=AMBER_EDGE,
        linewidth=1.5,
        zorder=8,
    )
    ax.add_patch(gate_marker)
    add_text(
        ax,
        gate_x,
        0.646,
        "BEFORE FINALIZE",
        size=5.8,
        color=AMBER_DARK,
        weight="bold",
        zorder=6,
    )
    add_text(
        ax,
        0.775,
        0.713,
        "suspect / refuted -> driver re-decides",
        size=6.2,
        color=AMBER_DARK,
        weight="bold",
    )

    # Lane C: authoritative ledger and generated artifacts.
    rounded_box(
        ax,
        0.045,
        0.135,
        0.475,
        0.120,
        facecolor=GRAY,
        edgecolor=GRAY_EDGE,
        linewidth=1.5,
        radius=0.012,
    )
    add_text(
        ax,
        0.065,
        0.229,
        "AUTHORITATIVE STORE  •  midmap_ledger.json",
        size=8.6,
        color=NAVY,
        weight="bold",
        ha="left",
    )
    add_text(
        ax,
        0.065,
        0.190,
        "name -> normalized -> candidates [with source] -> accepted -> final_class / confidence /\n"
        "origin -> back-checks -> gem_mam / gem_relation per model -> timestamped decisions log",
        size=7.15,
        color=INK,
        ha="left",
        linespacing=1.25,
    )
    add_text(
        ax,
        0.065,
        0.153,
        "All evidence and decisions remain traceable; candidate provenance gates record_decision.",
        size=6.3,
        color=MUTED,
        ha="left",
    )

    rounded_box(
        ax,
        0.545,
        0.135,
        0.410,
        0.120,
        facecolor=GRAY,
        edgecolor=GRAY_EDGE,
        linewidth=1.5,
        radius=0.012,
    )
    add_text(
        ax,
        0.565,
        0.229,
        "OUTPUT ARTIFACTS",
        size=8.6,
        color=NAVY,
        weight="bold",
        ha="left",
    )
    add_text(
        ax,
        0.565,
        0.193,
        "master_ledger.tsv  •  coverage_summary.tsv  •  gem_curation.tsv  •  exogenous_kept.tsv  •  xenobiotic_excluded.tsv",
        size=6.05,
        color=INK,
        ha="left",
    )
    add_text(
        ax,
        0.565,
        0.166,
        "figures/db_matching_upset.png  •  figures/db_matching_improvement.png",
        size=6.35,
        color=INK,
        ha="left",
    )
    add_text(
        ax,
        0.565,
        0.142,
        "code/reproduce_mapping.py (+ .ipynb)  •  annotated source",
        size=6.35,
        color=INK,
        ha="left",
    )

    # Anti-fabrication provenance gate: tool-produced ledger candidates feed the CALL.
    arrow(
        ax,
        (0.435, 0.257),
        (xs[5] + 0.020, 0.348),
        color=RED_EDGE,
        linewidth=1.35,
        connectionstyle="arc3,rad=-0.12",
        mutation_scale=10,
        zorder=6,
    )
    add_text(
        ax,
        0.372,
        0.300,
        "ledger candidates\n-> anti-fabrication gate",
        size=5.3,
        color=RED_EDGE,
        weight="bold",
        zorder=6,
        linespacing=1.15,
    )
    arrow(
        ax,
        (0.520, 0.195),
        (0.543, 0.195),
        color=GRAY_EDGE,
        linewidth=1.5,
        mutation_scale=10,
        zorder=5,
    )
    arrow(
        ax,
        (xs[7] + widths[7] / 2, stage_y - 0.004),
        (0.780, 0.257),
        color=GRAY_EDGE,
        linewidth=1.5,
        connectionstyle="arc3,rad=-0.06",
        mutation_scale=10,
        zorder=2,
    )

    # Anti-fabrication pill sits BELOW the CALL box (not as an in-box footer) so the two-axis
    # body has the full box height and never collides with it.
    pill_x = xs[5] + 0.004
    pill_w = widths[5] - 0.008
    pill_y = 0.349
    pill_h = 0.026
    rounded_box(
        ax,
        pill_x,
        pill_y,
        pill_w,
        pill_h,
        facecolor=RED,
        edgecolor=RED_EDGE,
        linewidth=1.1,
        radius=0.005,
        zorder=4,
    )
    add_text(
        ax,
        pill_x + pill_w / 2,
        pill_y + pill_h / 2,
        "anti-fabrication gate:\nno accepted ID without a backing tool candidate",
        size=4.7,
        color=RED_EDGE,
        weight="bold",
        zorder=5,
        linespacing=1.15,
    )

    # Bottom legends.
    rounded_box(
        ax,
        0.020,
        0.018,
        0.355,
        0.064,
        facecolor=WHITE,
        edgecolor="#CBD2D8",
        linewidth=1.0,
        radius=0.008,
    )
    add_text(ax, 0.033, 0.067, "ENCODING", size=7.3, color=NAVY, weight="bold", ha="left")
    legend_items = [
        (0.035, BLUE, BLUE_EDGE, "tools"),
        (0.086, AMBER, AMBER_EDGE, "LLM call"),
        (0.151, GRAY, GRAY_EDGE, "state / output"),
        (0.242, GREEN, GREEN_EDGE, "kept"),
        (0.292, RED, RED_EDGE, "excluded"),
    ]
    for x, fc, ec, label in legend_items:
        rounded_box(
            ax,
            x,
            0.031,
            0.013,
            0.018,
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.0,
            radius=0.002,
            zorder=3,
        )
        # +0.023, not +0.018: the swatch patch is 0.013 wide plus the 0.006 rounding pad on
        # each side, so a label at +0.018 lands inside the swatch it is labelling.
        add_text(ax, x + 0.023, 0.040, label, size=5.65, color=INK, ha="left")

    rounded_box(
        ax,
        0.385,
        0.018,
        0.595,
        0.064,
        facecolor=WHITE,
        edgecolor="#CBD2D8",
        linewidth=1.0,
        radius=0.008,
    )
    add_text(
        ax,
        0.398,
        0.067,
        "CONFIDENCE TIERS",
        size=7.3,
        color=NAVY,
        weight="bold",
        ha="left",
    )
    add_text(
        ax,
        0.398,
        0.040,
        "M1 exact  •  M2 verified synonym / abbrev  •  M3 verified typo  •  "
        "M4 structure-only  •  W mass-only (weak)  •  X non-endogenous  •  U unresolved",
        size=6.25,
        color=INK,
        ha="left",
    )
    add_text(
        ax,
        0.500,
        0.094,
        "Code (deterministic tools)  vs  LLM judgment (the CALL)",
        size=7.0,
        color=NAVY,
        weight="bold",
    )

    return fig, ax


def validate_layout(fig, ax):
    """Run basic renderer-level checks and return a concise diagnostics list."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    diagnostics = []

    width_px, height_px = fig.canvas.get_width_height()
    diagnostics.append(f"canvas={width_px}x{height_px}px")

    texts = [t for t in ax.texts if t.get_visible() and t.get_text().strip()]
    overlaps = []
    for i, left in enumerate(texts):
        lb = left.get_window_extent(renderer=renderer)
        for right in texts[i + 1 :]:
            rb = right.get_window_extent(renderer=renderer)
            x_overlap = min(lb.x1, rb.x1) - max(lb.x0, rb.x0)
            y_overlap = min(lb.y1, rb.y1) - max(lb.y0, rb.y0)
            if x_overlap > 2 and y_overlap > 2:
                # Deliberately stacked labels within a shared box never occupy the same y band;
                # any renderer-detected overlap is therefore actionable.
                overlaps.append((left.get_text().splitlines()[0], right.get_text().splitlines()[0]))
    diagnostics.append(f"text_overlap_pairs={len(overlaps)}")
    if overlaps:
        diagnostics.append(f"overlap_sample={overlaps[:8]}")

    return diagnostics, overlaps


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    fig, ax = build_diagram()
    diagnostics, overlaps = validate_layout(fig, ax)
    fig.savefig(
        OUTPUT,
        dpi=200,
        facecolor=WHITE,
        bbox_inches=None,
        pad_inches=0,
        metadata={"Software": "matplotlib", "Title": "metabo-idmapper MCP architecture"},
    )
    plt.close(fig)
    if overlaps:
        raise RuntimeError("layout validation found text overlaps: " + repr(overlaps[:8]))
    print("\n".join(diagnostics))
    print(f"{OUTPUT.resolve()} {OUTPUT.stat().st_size} bytes")


if __name__ == "__main__":
    main()
