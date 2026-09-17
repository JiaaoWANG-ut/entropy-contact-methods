#!/usr/bin/env python3
"""Nature-format figures and source-data tables for metal–Ni wetting."""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, NullFormatter
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

from run_contact import (
    AREA_PER_ATOM,
    DATA_DIR,
    ELEMENTS,
    MET_CUT,
    OUT_DIR,
    crossing_time,
    load_results,
    series,
    value_at,
)
from run_entropy import SYSTEMS, TEMPERATURES

NATURE_DIR = os.path.join(OUT_DIR, "nature")

# Wong colour-blind palette
ELEM_COLOR = {
    "Fe": "#D55E00",
    "Zn": "#0072B2",
    "Cr": "#009E73",
    "Ru": "#CC79A7",
}
MM = 1.0 / 25.4
COL_DOUBLE = 183 * MM  # Nature double-column width
T_MIN_NS = 3.5e-4
T_MAX_NS = 6.0


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Nimbus Sans",
                                "Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.labelsize": 7,
            "axes.titlesize": 6.5,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.5,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.minor.width": 0.35,
            "ytick.minor.width": 0.35,
            "xtick.major.size": 2.6,
            "ytick.major.size": 2.6,
            "xtick.minor.size": 1.4,
            "ytick.minor.size": 1.4,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "lines.linewidth": 1.05,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(which="both", direction="in")


def panel_label(ax, letter, x=-0.22, y=1.12):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")


def shade_wetting(ax):
    ax.axvspan(T_MIN_NS, 1.0, color="#f0f0f0", lw=0, zorder=0)
    ax.axvline(1.0, color="0.62", lw=0.45, ls=":", zorder=1)


def t_ns(tab):
    t = tab["time_ps"] / 1000.0
    t = np.where(t <= 0, T_MIN_NS, t)
    return np.clip(t, T_MIN_NS, None)


def plot_element_curve(ax, tab, el, ykey):
    """Solid while A_c is a wetting patch; dashed after a late collapse."""
    t = t_ns(tab)
    y = tab[ykey]
    cr = tab[f"cr_{el}"]
    i_peak = int(np.argmax(cr))
    collapsed = (tab["time_ps"] > tab["time_ps"][i_peak]) & (
        cr < cr[i_peak] - 0.15)
    color = ELEM_COLOR[el]
    if collapsed.any():
        i0 = int(np.flatnonzero(collapsed)[0])
        ax.plot(t[:i0], y[:i0], color=color, lw=1.05, zorder=3, solid_capstyle="round")
        ax.plot(t[i0 - 1:], y[i0 - 1:], color=color, lw=1.05, ls=(0, (2.4, 1.15)),
                zorder=3)
    else:
        ax.plot(t, y, color=color, lw=1.05, zorder=3, solid_capstyle="round")


def save(fig, stem):
    os.makedirs(NATURE_DIR, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        path = os.path.join(NATURE_DIR, f"{stem}.{ext}")
        fig.savefig(path)
        print(f"Saved: {path}")
    plt.close(fig)


def fig_contact(results, ykey, ylabel, stem, ypad=None):
    apply_style()
    fig, axes = plt.subplots(
        2, 5, figsize=(COL_DOUBLE, 3.55),
        sharex=True, sharey=True,
    )
    letters = "abcdefghij"
    k = 0
    ymax = 0.0
    for row, system in enumerate(SYSTEMS):
        short = "NiO-6" if "NiO-6" in system else "NiO-8"
        for col, temp_k in enumerate(TEMPERATURES):
            ax = axes[row, col]
            tab = series(results, system, temp_k)
            shade_wetting(ax)
            for el in ELEMENTS:
                plot_element_curve(ax, tab, el, f"{ykey}_{el}")
                ymax = max(ymax, float(tab[f"{ykey}_{el}"].max()))
            ax.set_xscale("log")
            ax.set_xlim(T_MIN_NS, T_MAX_NS)
            ax.xaxis.set_major_locator(LogLocator(base=10, numticks=6))
            ax.xaxis.set_minor_formatter(NullFormatter())
            style_ax(ax)
            panel_label(ax, letters[k])
            k += 1
            if row == 0:
                ax.set_title(f"{temp_k} K", pad=3, fontsize=6.5)
            if col == 0:
                ax.set_ylabel(f"{short}\n{ylabel}", labelpad=2)
            if row == 1:
                ax.set_xlabel("Time (ns)", labelpad=2)
    if ykey == "cr":
        for ax in axes.ravel():
            ax.set_ylim(-0.02, 1.05)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    else:
        top = (ypad or ymax) * 1.08
        for ax in axes.ravel():
            ax.set_ylim(0, top)

    handles = [Line2D([0], [0], color=ELEM_COLOR[el], lw=1.35, label=el)
               for el in ELEMENTS]
    handles += [
        Line2D([0], [0], color="0.3", lw=1.1, label="wetting"),
        Line2D([0], [0], color="0.3", lw=1.1, ls=(0, (2.4, 1.15)),
               label="A$_c$ ill-defined"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=6, frameon=False,
               bbox_to_anchor=(0.5, 1.04), handlelength=1.6,
               columnspacing=1.15, borderpad=0)
    fig.tight_layout(h_pad=0.55, w_pad=0.28, rect=(0, 0, 1, 0.96))
    save(fig, stem)


def fig_kinetics(results):
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(COL_DOUBLE, 3.70))
    x = np.arange(len(ELEMENTS))
    width = 0.15
    temp_colors = {
        1500: "#5B8FA8",
        1800: "#4C78A8",
        2000: "#F2A93B",
        2500: "#E07A3D",
        3000: "#C44E52",
    }
    letters = [["a", "b"], ["c", "d"]]
    for row, system in enumerate(SYSTEMS):
        short = "NiO-6" if "NiO-6" in system else "NiO-8"
        t10 = np.zeros((len(TEMPERATURES), len(ELEMENTS)))
        r1 = np.zeros_like(t10)
        for i, temp_k in enumerate(TEMPERATURES):
            tab = series(results, system, temp_k)
            t = tab["time_ps"]
            for j, el in enumerate(ELEMENTS):
                t10[i, j] = crossing_time(t, tab[f"cr_{el}"], 0.10)
                r1[i, j] = value_at(t, tab[f"cr_{el}"], 1000.0)
        ax = axes[row, 0]
        for i, temp_k in enumerate(TEMPERATURES):
            ax.bar(x + (i - 2) * width, t10[i], width=width * 0.92,
                   color=temp_colors[temp_k], label=f"{temp_k} K", zorder=3)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(ELEMENTS)
        ax.set_ylabel(f"{short}\n$t(\\rho=0.1)$ (ps)", labelpad=2)
        ax.set_ylim(5, 2000)
        style_ax(ax)
        panel_label(ax, letters[row][0], x=-0.18, y=1.08)
        ax.set_xlabel("Element", labelpad=2)

        ax = axes[row, 1]
        for i, temp_k in enumerate(TEMPERATURES):
            ax.bar(x + (i - 2) * width, 100 * r1[i], width=width * 0.92,
                   color=temp_colors[temp_k], label=f"{temp_k} K", zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(ELEMENTS)
        ax.set_ylabel("Contact ratio at 1 ns (%)", labelpad=2)
        ax.set_ylim(0, 100)
        style_ax(ax)
        panel_label(ax, letters[row][1], x=-0.18, y=1.08)
        ax.set_xlabel("Element", labelpad=2)

    handles = [Rectangle((0, 0), 1, 1, color=temp_colors[T], label=f"{T} K")
               for T in TEMPERATURES]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, 1.03), handlelength=1.1,
               columnspacing=1.3, borderpad=0)
    fig.tight_layout(h_pad=0.85, w_pad=1.1, rect=(0, 0, 1, 0.95))
    save(fig, "Fig3_wetting_kinetics")


def write_source_xlsx(results):
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    header_font = Font(bold=True, name="Calibri", size=11)
    body = [
        ["Quantity", "Definition", "Units"],
        ["contact_ratio",
         "Fraction of X metal atoms (X = Fe, Zn, Cr, Ru) with ≥1 Ni neighbour "
         f"within {MET_CUT} Å; oxygen excluded. ρ = n_c / N_X.", "1"],
        ["n_contact", "Number of X atoms in contact with Ni", "atoms"],
        ["n_pairs", "Number of X–Ni neighbour pairs within the cutoff", "pairs"],
        ["area_A2",
         f"Contact-area proxy A_c = n_contact × {AREA_PER_ATOM:.4f} Å² "
         "(close-packed atomic area, d_nn = 2.50 Å).", "Å^2"],
        ["mean_cn_ni", "Mean Ni coordination of all X atoms", "1"],
        ["alloyed_frac", "Fraction of X atoms with CN_Ni ≥ 4 (embedded, not a "
         "surface contact)", "1"],
        ["attached_frac", "Fraction of X atoms in the Ni-containing metal "
         "connected component", "1"],
        ["com_dist_A", "Distance between X COM and Ni COM", "Å"],
        ["time_ps", "Dump Time attribute", "ps"],
        [],
        ["Figures", "Contents", ""],
        ["Fig1_contact_ratio",
         "2×5: rows NiO-6 / NiO-8, columns 1500–3000 K. y = contact ratio, "
         "x = log time (ns). Grey band = 0–1 ns wetting window. Dashed = "
         "A_c ill-defined after late collapse (3000 K evaporation).", ""],
        ["Fig2_contact_area", "Same layout, y = A_c (Å²).", ""],
        ["Fig3_wetting_kinetics",
         "a,c: t(ρ=0.1) on a log scale. b,d: contact ratio at 1 ns.", ""],
        [],
        ["Sampling",
         "t = 0 from model.xyz plus 122 logarithmically spaced dump frames "
         "in 0.4 ps–5 ns (denser below 1 ns).", ""],
        ["Protocol",
         "GPUMD NEP89, NVT-Langevin, dt = 1 fs, dump every 400 steps (0.4 ps).",
         ""],
    ]
    for row in body:
        ws.append(row)
    for cell in ws[1]:
        cell.font = header_font
    ws.column_dimensions["A"].width = 22
    ws.column_dimensions["B"].width = 92
    ws.column_dimensions["C"].width = 12

    ws = wb.create_sheet("summary")
    ws.append(["system", "T_K", "element", "N_X", "t10_ps", "t50_ps",
               "rho_1ns", "rho_5ns", "Ac_1ns_A2", "Ac_5ns_A2",
               "alloy_1ns", "alloy_5ns"])
    for system in SYSTEMS:
        for temp_k in TEMPERATURES:
            tab = series(results, system, temp_k)
            t = tab["time_ps"]
            for el in ELEMENTS:
                cr = tab[f"cr_{el}"]
                ws.append([
                    system, temp_k, el, int(tab[f"n_{el}"][0]),
                    crossing_time(t, cr, 0.10),
                    crossing_time(t, cr, 0.50),
                    value_at(t, cr, 1000.0),
                    value_at(t, cr, 5000.0),
                    value_at(t, tab[f"area_{el}"], 1000.0),
                    value_at(t, tab[f"area_{el}"], 5000.0),
                    value_at(t, tab[f"alloy_{el}"], 1000.0),
                    value_at(t, tab[f"alloy_{el}"], 5000.0),
                ])

    ws = wb.create_sheet("source_long")
    ws.append(["system", "temperature_K", "frame_index", "time_ps", "time_ns",
               "element", "N_X", "n_contact", "n_pairs", "contact_ratio",
               "area_A2", "mean_cn_ni", "alloyed_frac", "attached_frac",
               "com_dist_A"])
    for system in SYSTEMS:
        for temp_k in TEMPERATURES:
            tab = series(results, system, temp_k)
            for r in tab:
                for el in ELEMENTS:
                    ws.append([
                        system, temp_k, int(r["frame"]), float(r["time_ps"]),
                        float(r["time_ps"]) / 1000.0, el,
                        int(r[f"n_{el}"]), int(r[f"nc_{el}"]),
                        int(r[f"np_{el}"]), float(r[f"cr_{el}"]),
                        float(r[f"area_{el}"]), float(r[f"cn_{el}"]),
                        float(r[f"alloy_{el}"]), float(r[f"att_{el}"]),
                        float(r[f"dist_{el}"]),
                    ])

    for system in SYSTEMS:
        short = "NiO6" if "NiO-6" in system else "NiO8"
        for metric, tag in (("cr", "ratio"), ("area", "area")):
            ws = wb.create_sheet(f"{short}_{tag}"[:31])
            header = ["time_ns"]
            cols = []
            for temp_k in TEMPERATURES:
                tab = series(results, system, temp_k)
                t = t_ns(tab)
                cols.append(("time", t, temp_k, None))
                for el in ELEMENTS:
                    header.append(f"{temp_k}K_{el}")
            # Use NiO-6 1500 K times as the first time column; other T share
            # the same log grid to within dump spacing, so write per-T blocks.
            ws.append(["Note: each temperature has its own time column "
                       "because dump frames are log-sampled independently."])
            header = []
            arrays = []
            for temp_k in TEMPERATURES:
                tab = series(results, system, temp_k)
                header.append(f"{temp_k}K_time_ns")
                arrays.append(t_ns(tab))
                for el in ELEMENTS:
                    header.append(f"{temp_k}K_{el}")
                    arrays.append(tab[f"{metric}_{el}"].astype(float))
            ws.append(header)
            n = max(len(a) for a in arrays)
            for i in range(n):
                ws.append([float(a[i]) if i < len(a) else None for a in arrays])

    path = os.path.join(NATURE_DIR, "Fig1-3_source_data.xlsx")
    wb.save(path)
    print(f"Saved: {path}")
    return path


def write_captions():
    text = """Figure captions (Nature style)

Fig. 1 | Metal–Ni contact ratio versus logarithmic time.
Rows are the two NiO particles (NiO-6, 6,200 Ni; NiO-8, 14,672 Ni). Columns
are temperature. An X atom (X = Fe, Zn, Cr, Ru) is counted as contacting Ni
when it has at least one Ni neighbour within 3.0 Å; oxygen is excluded. The
contact ratio is ρ = n_c / N_X. The grey band is the 0–1 ns wetting window.
Solid lines: ρ traces a droplet/film contact patch. Dashed lines: A_c is no
longer a wetting area (3000 K, after the particle evaporates and ρ collapses).
Colour: Fe orange, Zn blue, Cr green, Ru purple (Wong palette).

Fig. 2 | Metal–Ni contact-area proxy versus logarithmic time.
Same layout as Fig. 1. A_c = n_c × 5.41 Å², the close-packed atomic area with
d_nn = 2.50 Å. At 2000–2500 K, A_c saturates near N_X × 5.41 Å² once the
metal has spread into a film, so the number ceases to be a 3D contact-patch
area. At 3000 K, A_c peaks inside the wetting window and then falls as the
condensed X/Ni interface disappears.

Fig. 3 | Wetting kinetics.
a, c, Time at which ρ first reaches 0.1 (log scale). b, d, Contact ratio at
1 ns, the end of the wetting window. Zn is the fastest-wetting metal at every
temperature; Ru forms the smallest contact patch. Raising T from 1500 K to
3000 K shortens t(ρ=0.1) by roughly two orders of magnitude.

Notes.
Bulk alloying is not reached in 5 ns: the fraction of X atoms with CN_Ni ≥ 4
remains below 0.22. The 3000 K collapse of ρ is evaporation of the particle,
consistent with the rise of S_config and the shrinkage of the largest metal
cluster in the companion entropy analysis.
"""
    path = os.path.join(NATURE_DIR, "figure_captions.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"Saved: {path}")
    return path


def write_email_body():
    text = """Dear all,

Metal-on-Ni wetting analysis for the FIG4 NEP89 runs (NiO-6 / NiO-8 + Fe, Zn,
Cr, Ru). Nature-format figures, captions, and the original plotting tables are
on GitHub Releases (link below). This note is the figure explanation.

Setup
  GPUMD NEP89, NVT-Langevin, dt = 1 fs, dump every 0.4 ps.
  Central NiO particle with four oxide particles (Fe, Zn, Cr, Ru) around it.
  10 trajectories: two sizes x 1500, 1800, 2000, 2500, 3000 K, first 5 ns.
  Sampling: t = 0 plus 122 logarithmically spaced frames, denser in 0-1 ns.

Definition (metal part only; oxygen excluded)
  An X atom contacts Ni if it has >=1 Ni neighbour within 3.0 A.
  Contact ratio  rho = n_c / N_X.
  Contact area   A_c = n_c x 5.41 A^2  (close-packed atomic area, d_nn = 2.50 A).
  When the droplet becomes a complete film (rho -> 1) or the particle evaporates,
  A_c is no longer a 3D contact-patch area; those segments are dashed in Fig. 1-2.

What the figures show

  Fig. 1  Contact ratio vs log time, 2 x 5 (rows = NiO-6/NiO-8, columns = T).
          Grey band = 0-1 ns wetting window.

  Fig. 2  Same layout for A_c.

  Fig. 3  a,c: time to rho = 0.1 (log scale). b,d: rho at 1 ns.

Wetting (0-1 ns)
  Ranking is Zn > Fe ~ Cr > Ru at every temperature.
  t(rho=0.1) is 140-900 ps at 1500 K and 7-12 ps at 3000 K.
  NiO-8 reaches a higher 1 ns rho than NiO-6 (more Ni surface).

After 1 ns, three regimes
  1500-1800 K: rho still rising at 5 ns. Distinct droplets, A_c is a valid
               contact patch.
  2000-2500 K: Zn/Fe/Cr saturate at rho = 0.73-0.92. Thin spread film, not a
               bulk alloy (fraction with CN_Ni >= 4 stays < 0.22). A_c has
               hit the N_X ceiling and no longer tracks a 3D contact line.
  3000 K:      rho peaks at 0.4-0.8 ns then collapses (NiO-6 Cr 67% -> 10%;
               Ru 54% -> 6%). Evaporation, not mixing: attached fraction -> 0,
               matching the S_config rise in the entropy analysis.

Source data
  Fig1-3_source_data.xlsx
    README, summary (t10, t50, rho/A_c at 1 and 5 ns), source_long (every
    analysed frame), and plot-ready wide tables per system.

Best regards,
Jiaao
"""
    path = os.path.join(NATURE_DIR, "email_body.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"Saved: {path}")
    return path


def write_release_notes(url_placeholder: str = ""):
    path = os.path.join(NATURE_DIR, "RELEASE_NOTES.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "Nature-format metal–Ni wetting figures for the FIG4 GPUMD runs "
            "(NEP89, NVT-Langevin, dt = 1 fs, dump every 0.4 ps).\n\n"
            "Metric\n"
            f"- Contact: X–Ni neighbours within {MET_CUT} Å, oxygen excluded\n"
            "- Contact ratio ρ = n_c / N_X for X = Fe, Zn, Cr, Ru\n"
            f"- Contact area A_c = n_c × {AREA_PER_ATOM:.2f} Å²\n"
            "- Log time axis; grey band = 0–1 ns wetting window\n"
            "- Dashed curves: A_c ill-defined after the 3000 K collapse\n\n"
            "Assets\n"
            "- Fig1_contact_ratio.pdf/.png/.svg\n"
            "- Fig2_contact_area.pdf/.png/.svg\n"
            "- Fig3_wetting_kinetics.pdf/.png/.svg\n"
            "- Fig1-3_source_data.xlsx  (README + long + wide tables)\n"
            "- figure_captions.txt\n"
            "- FIG4_wetting_nature.zip  (all of the above)\n"
        )
    print(f"Saved: {path}")
    return path


def build():
    os.makedirs(NATURE_DIR, exist_ok=True)
    results = load_results()
    if not results:
        raise SystemExit("Run run_contact.py first.")
    fig_contact(results, "cr", "Contact ratio  $\\rho$", "Fig1_contact_ratio")
    fig_contact(results, "area", "Contact area  $A_c$  (Å$^2$)",
                "Fig2_contact_area")
    fig_kinetics(results)
    write_source_xlsx(results)
    write_captions()
    write_email_body()
    write_release_notes()
    # copy the long CSV next to the Nature files for the zip
    src = os.path.join(OUT_DIR, "contact_raw_long.csv")
    dst = os.path.join(NATURE_DIR, "contact_raw_long.csv")
    if os.path.isfile(src) and os.path.abspath(src) != os.path.abspath(dst):
        import shutil
        shutil.copy2(src, dst)
        print(f"Copied: {dst}")


if __name__ == "__main__":
    build()
