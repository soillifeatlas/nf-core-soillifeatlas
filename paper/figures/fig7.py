"""Figure 7 - ClimGrass application. A design (real); B-E gated on the stage-5 fix."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (DOUBLE, MM, OI, POL, POL_LABEL, Panels, gates, letter, pending, read_table,
                    save, set_panel)

POLS = ("pos", "neg")
TREAT = ["ambient", "future", "ambient_drought", "future_drought"]
TREAT_LABEL = {"ambient": "ambient", "future": "future", "ambient_drought": "ambient + drought",
               "future_drought": "future + drought"}
TREAT_COLOR = {"ambient": "#BBBBBB", "future": OI["orange"], "ambient_drought": OI["green"],
               "future_drought": OI["purple"]}
STAGE5 = "stage-5 quantification\nfix not yet confirmed"


def panel_A(ax):
    set_panel("7A")
    s = read_table("00_inputs/climgrass_samples.csv")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.set_axis_off()
    for i, wc in enumerate(["no", "yes"]):
        for j, dr in enumerate(["no", "yes"]):
            cell = s[(s.warming_co2 == wc) & (s.drought == dr)]
            t = cell.treatment.iloc[0]
            ax.add_patch(plt.Rectangle((i + 0.04, 1 - j + 0.04), 0.92, 0.92,
                                       facecolor=TREAT_COLOR[t], alpha=0.35, edgecolor="none"))
            lines = []
            for _, r in cell.iterrows():
                pols = [p for p, c in (("POS", r.used_pos), ("NEG", r.used_neg)) if c == "yes"]
                lines.append(f"plot {r['plot']}" + ("" if len(pols) == 2 else f" ({pols[0]} only)"))
            ax.text(i + 0.5, 1 - j + 0.72, TREAT_LABEL[t], ha="center", fontweight="bold",
                    fontsize=6.5)
            ax.text(i + 0.5, 1 - j + 0.36, "\n".join(lines), ha="center", va="center",
                    fontsize=6)
    ax.text(0.5, 2.06, "ambient climate", ha="center", fontsize=6)
    ax.text(1.5, 2.06, "+3 °C, +300 ppm CO$_2$", ha="center", fontsize=6)
    ax.text(-0.04, 1.5, "no drought", ha="right", va="center", fontsize=6, rotation=90)
    ax.text(-0.04, 0.5, "drought", ha="right", va="center", fontsize=6, rotation=90)
    n_pos, n_neg = (s.used_pos == "yes").sum(), (s.used_neg == "yes").sum()
    ax.text(1.0, -0.1, f"{n_pos} soils POS, {n_neg} NEG", ha="center", va="top", fontsize=6)
    letter(ax, "A", x=-0.1, y=1.08, title="ClimGrass 2 × 2 design")


def panel_A2(ax):
    ax.set_axis_off()
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="white",
                               edgecolor="#777777", lw=0.6))
    ax.text(0.5, 0.5, "Drying–rewetting validation\n\npre-registered and blinded;\n"
            "design shown only after unblinding\n(issue #11)", transform=ax.transAxes,
            ha="center", va="center", fontsize=6, color="#333333")


def panel_B(ax, top=8):
    set_panel("7B")
    d = {p: read_table(f"08_application/{p}/class_profiles.tsv") for p in POLS}
    if not gates()["stage5_fix_confirmed"]:
        pending(ax, STAGE5)
        letter(ax, "B", x=0.0, y=1.03, title="Class amounts by treatment")
        return "PENDING", "gated on stage5_fix_confirmed"
    p = d["pos"]
    order = p.groupby("lipid_class")["mean"].max().sort_values(ascending=False).head(top).index
    x = np.arange(len(order))
    w = 0.2
    for k, t in enumerate(TREAT):
        s = p[p.treatment == t].set_index("lipid_class").reindex(order)
        ax.bar(x + (k - 1.5) * w, s["mean"], w, yerr=s["std"], color=TREAT_COLOR[t],
               error_kw=dict(lw=0.4, capsize=1), label=TREAT_LABEL[t])
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=90)
    ax.set_ylabel(r"µg g$^{-1}$ (mean ± SD of plots)")
    ax.legend(fontsize=5.5, ncol=2)
    letter(ax, "B", x=-0.15, y=1.03, title="Class amounts by treatment (POS, common set)")
    return "REAL", "POS shown; NEG read"


def panel_C(ax):
    set_panel("7C")
    d = {p: read_table(f"08_application/{p}/responsive.tsv") for p in POLS}
    for p in POLS:
        read_table(f"08_application/{p}/test_summary.tsv")
    if not gates()["stage5_fix_confirmed"]:
        pending(ax, STAGE5)
        letter(ax, "C", x=0.0, y=1.03, title="Drought-responsive lipids")
        return "PENDING", "gated on stage5_fix_confirmed"
    rows = pd.concat([d[p].assign(pol=p) for p in POLS]).sort_values("log2fc_drought_vs_ambient")
    y = np.arange(len(rows))
    ax.barh(y, rows.log2fc_drought_vs_ambient, color=[POL[p] for p in rows.pol], height=0.7)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(rows.lipid_name, fontsize=5)
    ax.set_xlabel(r"log$_2$ fold change, drought vs ambient")
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=POL[p]) for p in POLS],
              [POL_LABEL[p] for p in POLS], loc="lower right", fontsize=5.5)
    letter(ax, "C", x=-0.45, y=1.03, title="Drought-responsive lipids (q < 0.05)")
    return "REAL", ""


def panel_D(ax):
    set_panel("7D")
    d = read_table("06_atlas_match/examples/responsive_atlas_occurrence.tsv")
    if not gates()["stage5_fix_confirmed"]:
        pending(ax, STAGE5)
        letter(ax, "D", x=0.0, y=1.03, title="Atlas occurrence of responsive lipids")
        return "PENDING", "gated on stage5_fix_confirmed"
    d = d.sort_values("n_ref_phyla", na_position="first")
    y = np.arange(len(d))
    ax.barh(y, d.n_ref_phyla.fillna(0), color=OI["green"], height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(d.lipid_name, fontsize=5)
    ax.set_xlabel("phyla with the matched reference feature")
    letter(ax, "D", x=-0.3, y=1.03, title="Atlas occurrence (POS responsive lipids)")
    return "REAL", "no bar = no Atlas match"


def panel_E(ax):
    set_panel("7E")
    resp = read_table("08_application/pos/responsive.tsv")
    occ = read_table("06_atlas_match/examples/responsive_atlas_occurrence.tsv")
    g = gates()
    if not g["stage5_fix_confirmed"] or not g.get("fig7E_lipid"):
        why = STAGE5 if not g["stage5_fix_confirmed"] else "lipid for the trace\nnot chosen"
        pending(ax, why)
        letter(ax, "E", x=0.0, y=1.03, title="One lipid traced end to end")
        return "PENDING", "gated on stage5_fix_confirmed and fig7E_lipid"
    r = resp[resp.lipid_name == g["fig7E_lipid"]].iloc[0]
    o = occ[occ.feature_id == r.feature_id].iloc[0]
    vals = [r[f"mean_ng_g_{t}"] for t in TREAT]
    ax.bar(range(4), vals, color=[TREAT_COLOR[t] for t in TREAT])
    ax.set_xticks(range(4))
    ax.set_xticklabels([TREAT_LABEL[t] for t in TREAT], rotation=30, ha="right")
    ax.set_ylabel(r"ng g$^{-1}$ (mean of plots)")
    ax.set_ylim(0, max(vals) * 1.6)
    ax.text(0.02, 0.98, f"{r.lipid_name} ({r.confidence}, route {r.quant_route})\n"
            f"q(drought) = {r.q_drought:.3g}\nAtlas: cosine {o.cosine:.2f}, "
            f"{int(o.n_ref_phyla)} phyla", transform=ax.transAxes, va="top", fontsize=5.8)
    letter(ax, "E", x=-0.2, y=1.03, title="One lipid traced end to end")
    return "REAL", f"{r.lipid_name}"


def build():
    P = Panels(7)
    fig = plt.figure(figsize=(DOUBLE, 150 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], hspace=0.5, wspace=0.45,
                          left=0.1, right=0.98, top=0.92, bottom=0.1)
    top = gs[0, 0].subgridspec(2, 1, height_ratios=[1, 0.42], hspace=0.35)
    panel_A(fig.add_subplot(top[0]))
    panel_A2(fig.add_subplot(top[1]))
    P.add("A", "REAL", "ClimGrass design from climgrass_samples.csv; drying–rewetting as a "
          "generic text box (blinded, #11)")
    for tag, fn, spot in (("B", panel_B, gs[0, 1:]), ("C", panel_C, gs[1, 0]),
                          ("D", panel_D, gs[1, 1]), ("E", panel_E, gs[1, 2])):
        st, note = fn(fig.add_subplot(spot))
        P.add(tag, st, note)
    return save(fig, 7, "application"), P.rows
