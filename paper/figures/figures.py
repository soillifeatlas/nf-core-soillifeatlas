"""The figure set: which panels make up each main and supplementary figure.

Main text tells the resource story: the Atlas (Fig 1, unchanged here), what it adds to a soil
dataset (Fig 2), the workflow (Fig 3), interpreting single features (Fig 4) and the ClimGrass
demonstration (Fig 5). Analytical detail goes to the Supplement (Figs S1-S5).

A figure that needs an unconfirmed input is not drawn at all (no grey placeholder boxes): its
STATUS row says NOT BUILT and which switch in gates.yaml it waits for.
"""
from __future__ import annotations

import matplotlib.pyplot as plt

import panels_annotation as ann
import panels_application as app
import panels_atlas as atlas
import panels_interp as interp
import panels_qc as qc
import panels_quant as quant
import panels_workflow as wf
from common import DOUBLE, MM, Panels, fig_letter, gates, input_path, save, set_panel, show_image


def _not_built(figure, name, role, why, panels):
    P = Panels(figure, name, role)
    for p in panels:
        P.add(p, "NOT BUILT", why)
    return [], P.rows


# ------------------------------------------------------------------------------------ main
def fig2_atlas():
    P = Panels(2, "atlas_matching")
    fig = plt.figure(figsize=(DOUBLE, 150 * MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.75, 1.3], width_ratios=[1, 1.05, 0.95],
                          hspace=0.45, wspace=0.75, left=0.17, right=0.98, top=0.93, bottom=0.09)
    title = atlas.coverage(fig.add_subplot(gs[0, :]))
    fig_letter(fig, "A", 0.005, 0.985, title)
    P.add("A", "REAL", "POS; " + atlas.match_settings(), ["atlas_coverage"])
    atlas.breadth(fig.add_subplot(gs[1, 0]), "B")
    P.add("B", "REAL", "occurrence = phyla where the matched reference feature is detected",
          ["atlas_breadth"])
    atlas.per_phylum(fig.add_subplot(gs[1, 1]), "C")
    P.add("C", "REAL", "raw counts; phyla differ in reference sample numbers (in brackets)",
          ["atlas_per_phylum"])
    atlas.class_breadth(fig.add_subplot(gs[1, 2]), "D")
    P.add("D", "REAL", "10 most frequent classes among annotated matches",
          ["atlas_class_breadth"])
    return save(fig, 2, "atlas_matching"), P.rows


def fig3_workflow():
    P = Panels(3, "workflow")
    fig = plt.figure(figsize=(DOUBLE, 170 * MM))
    axA = fig.add_axes([0.0, 0.80, 1.0, 0.16])
    wf.stages(axA)
    fig_letter(fig, "A", 0.005, 0.995, "From raw files to Atlas-linked lipid profiles "
               "(dashed: run as scripts, not yet Nextflow processes)")
    P.add("A", "SCHEMATIC", "crop of the 2026-11-24 draft; stages 5–7 dashed",
          ["schematic_stages"])
    axB = fig.add_axes([0.02, 0.43, 0.44, 0.31])
    wf.tools_matrix(axB)
    fig_letter(fig, "B", 0.005, 0.775, "Where tools and Atlas layers enter")
    P.add("B", "SCHEMATIC", "crop of the 2026-11-24 draft, unchanged", ["schematic_tools"])
    gC = fig.add_gridspec(1, 2, left=0.58, right=0.98, top=0.735, bottom=0.47, wspace=0.12)
    qc.panel_B(fig, gC, tag="C", compact=True)
    P.add("C", "REAL", "pooled-QC RSD: peak-picked only vs including gap-filled", ["3B"])
    axD = fig.add_axes([0.08, 0.0, 0.84, 0.355])
    wf.report(axD)
    fig_letter(fig, "D", 0.005, 0.385, "One-page report (ClimGrass, positive mode)")
    P.add("D", "REAL", "top of report_pos.html; the 'lipids with an amount' tile is stage-5 "
          "output", ["report"])
    return save(fig, 3, "workflow"), P.rows


def fig4_interpretation():
    P = Panels(4, "interpretation")
    fig = plt.figure(figsize=(DOUBLE, 120 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], hspace=0.55, wspace=0.42, left=0.08,
                          right=0.98, top=0.93, bottom=0.1)
    interp.mirror(fig.add_subplot(gs[0, 0]), "broad", "A")
    P.add("A", "REAL", "example proposed, to confirm", ["mirror_broad"])
    interp.mirror(fig.add_subplot(gs[1, 0]), "restricted", "B")
    P.add("B", "REAL", "example proposed, to confirm", ["mirror_restricted"])
    interp.phylum_detection(fig.add_subplot(gs[:, 1]), "C")
    P.add("C", "REAL", "only phyla where either lipid is detected", ["phylum_detection"])
    return save(fig, 4, "interpretation"), P.rows


def fig5_application():
    g = gates()
    panels = ["A", "B", "C"] + (["D"] if g.get("fig7E_lipid") else [])
    if not g["stage5_fix_confirmed"]:
        return _not_built(5, "application", "main", "waits for stage5_fix_confirmed", panels)
    P = Panels(5, "application")
    fig = plt.figure(figsize=(DOUBLE, 150 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1.05, 1.05], height_ratios=[1, 1],
                          hspace=0.45, wspace=0.75, left=0.06, right=0.98, top=0.93,
                          bottom=0.1)
    top = gs[0, 0].subgridspec(2, 1, height_ratios=[1, 0.42], hspace=0.4)
    app.panel_A(fig.add_subplot(top[0]), tag="A")
    app.panel_A2(fig.add_subplot(top[1]))
    P.add("A", "REAL", "ClimGrass design; drying–rewetting only as a text box (#11)", ["7A"])
    app.panel_C(fig.add_subplot(gs[:, 1]), tag="B")
    P.add("B", "REAL", "drought-responsive lipids, q < 0.05", ["7C"])
    app.panel_D(fig.add_subplot(gs[:, 2]), tag="C")
    P.add("C", "REAL", "Atlas occurrence of the POS responsive lipids", ["7D"])
    if len(panels) == 4:
        app.panel_E(fig.add_subplot(gs[1, 0]), tag="D")
        P.add("D", "REAL", f"traced lipid: {g['fig7E_lipid']}", ["7E"])
    return save(fig, 5, "application"), P.rows


# ------------------------------------------------------------------------------ supplement
def figS1_quality():
    P = Panels("S1", "quality", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 165 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[0.55, 1.0, 1.35], hspace=0.62, left=0.09,
                             right=0.98, top=0.95, bottom=0.06)
    qc.panel_A(fig, outer[0].subgridspec(2, 1, hspace=1.5), tag="A")
    P.add("A", "REAL", "POS and NEG sequences, both blocks", ["3A"])
    qc.panel_B(fig, outer[1].subgridspec(1, 2, wspace=0.12), tag="B")
    P.add("B", "REAL", "ECDF, x cut at 150 %; n and % ≤ 30 % from rsd_summary.tsv", ["3B"])
    qc.panel_C(fig, outer[2].subgridspec(2, 2, hspace=0.12, wspace=0.15), tag="C")
    P.add("C", "REAL", "RT relative to recorded block median; tolerances from provenance.yaml",
          ["3C"])
    return save(fig, "S1", "quality"), P.rows


def figS2_annotation():
    P = Panels("S2", "annotation", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 185 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 0.95], hspace=0.5, left=0.08,
                             right=0.98, top=0.93, bottom=0.07)
    set_panel("4A")
    top = outer[0].subgridspec(1, 2, wspace=0.12)
    ann.upset(fig, top[0], "pos", True, tag="A")
    ann.upset(fig, top[1], "neg", False, tag="A")
    P.add("A", "REAL", "exclusive intersections, POS and NEG", ["4A"])
    mid = outer[1].subgridspec(1, 2, width_ratios=[1.1, 1], wspace=0.3)
    ann.panel_B(fig.add_subplot(mid[0]), tag="B")
    P.add("B", "REAL", "ClimGrass only", ["4B"])
    ann.panel_D(fig, mid[1], tag="D")
    P.add("D", "REAL", "PC; untested lipids counted, not drawn", ["4D"])
    ann.panel_C(fig, outer[2].subgridspec(1, 2, wspace=0.15), tag="C")
    P.add("C", "REAL", "12 largest classes per polarity", ["4C"])
    return save(fig, "S2", "annotation"), P.rows


def figS3_response():
    P = Panels("S3", "response_routes", "supplement")
    include_q3 = bool(gates()["stage5_fix_confirmed"])
    fig = plt.figure(figsize=(DOUBLE, 190 * MM))
    outer = fig.add_gridspec(3, 1, height_ratios=[1.05, 0.62, 1.0], hspace=0.42, left=0.08,
                             right=0.98, top=0.95, bottom=0.07)
    r1 = outer[0].subgridspec(1, 2, width_ratios=[1.55, 1], wspace=0.28)
    quant.panel_A(fig.add_subplot(r1[0]), tag="A")
    P.add("A", "REAL", "RIE relative to LPE", ["5A"])
    quant.panel_B(fig.add_subplot(r1[1]), tag="B")
    P.add("B", "REAL", "leave-one-class-out; statistics from ie_model_stats.tsv", ["5B"])
    set_panel("5D")
    axC = fig.add_subplot(outer[1])
    show_image(axC, input_path(quant.SCHEM), quant.CROP_5D)
    fig_letter(fig, "C", 0.01, axC.get_position().y1 + 0.035, "Quantification route per feature")
    P.add("C", "SCHEMATIC", "crop of the 2026-11-24 draft; Q2 not used in this run", ["5D"])
    r3 = outer[2].subgridspec(1, 2, width_ratios=[1, 1.2], wspace=0.3)
    quant.panel_E(fig.add_subplot(r3[0]), tag="D")
    P.add("D", "REAL", "Q1 and Q3" if include_q3 else "Q1 only; Q3 waits for "
          "stage5_fix_confirmed", ["5E"])
    return save(fig, "S3", "response_routes"), P.rows


def figS4_recovery():
    if not gates()["reconstitution_volume_confirmed"]:
        return _not_built("S4", "recovery", "supplement",
                          "waits for reconstitution_volume_confirmed", ["A", "B"])
    P = Panels("S4", "recovery", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 85 * MM))
    gs = fig.add_gridspec(1, 2, wspace=0.5, left=0.15, right=0.98, top=0.88, bottom=0.25)
    qc.panel_D(fig, gs[0, 0], tag="A")
    P.add("A", "REAL", "IS in soil vs calibration range", ["3D"])
    quant.panel_C(fig.add_subplot(gs[0, 1]), tag="B")
    P.add("B", "REAL", "apparent recovery", ["5C"])
    return save(fig, "S4", "recovery"), P.rows


def figS5_amounts():
    if not gates()["stage5_fix_confirmed"]:
        return _not_built("S5", "amounts", "supplement", "waits for stage5_fix_confirmed",
                          ["A", "B"])
    P = Panels("S5", "amounts", "supplement")
    fig = plt.figure(figsize=(DOUBLE, 85 * MM))
    gs = fig.add_gridspec(1, 2, wspace=0.3, left=0.08, right=0.98, top=0.88, bottom=0.2)
    quant.panel_F(fig.add_subplot(gs[0, 0]), tag="A")
    P.add("A", "REAL", "class amounts, common lipid set", ["5F"])
    app.panel_B(fig.add_subplot(gs[0, 1]), tag="B")
    P.add("B", "REAL", "class amounts by treatment (POS)", ["7B"])
    return save(fig, "S5", "amounts"), P.rows


FIGURES = [fig2_atlas, fig3_workflow, fig4_interpretation, fig5_application,
           figS1_quality, figS2_annotation, figS3_response, figS4_recovery, figS5_amounts]
