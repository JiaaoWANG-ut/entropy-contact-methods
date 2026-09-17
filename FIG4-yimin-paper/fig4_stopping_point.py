#!/usr/bin/env python3
"""Nature Fig. 4: selecting a morphological stopping point."""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import label as cc_label

from run_contact import ELEMENTS, load_results, series
from run_entropy import TEMPERATURES

ROOT = os.path.dirname(os.path.abspath(__file__))
PNG_ROOT = "/public/home/jwang/test/testgpumd/ovito-figures/png"
ENT_CSV = os.path.join(ROOT, "entropy_analysis_5ns", "entropy_raw_long.csv")
OUT = os.path.join(ROOT, "contact_analysis", "nature")
SYSTEM = "NiO-8-Fe-Zn-Cr-Ru"

STAGE_C = {1: "#F6E6C4", 2: "#C5DCCF", 3: "#EBC4C1"}
STAGE_EDGE = {1: "#C9A45A", 2: "#5F8F78", 3: "#B36A66"}
T_SHOW = (1500, 2000, 3000)
T_COLOR = {1500: "#3D6B93", 2000: "#C9922A", 3000: "#B24745"}
MM = 1.0 / 25.4
T_MIN, T_MAX = 0.025, 5.4


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Nimbus Sans",
                                "Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6.5,
            "axes.labelsize": 7,
            "axes.titlesize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.45,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 500,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


def style(ax, box=False):
    ax.spines["top"].set_visible(box)
    ax.spines["right"].set_visible(box)
    for sp in ax.spines.values():
        sp.set_color("0.25")


def panel(ax, letter, x=-0.10, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left", color="0.08")


def load_entropy():
    return np.genfromtxt(
        ENT_CSV, delimiter=",", names=True, dtype=None,
        encoding="utf-8", invalid_raise=False,
    )


def entropy_series(rows, temp_k):
    m = (rows["system"] == SYSTEM) & (rows["temperature_K"] == temp_k)
    t = rows["time_ps"][m]
    order = np.argsort(t)
    return t[order], rows["S_atom"][m][order], rows["S_config"][m][order]


def mean_rho(tab):
    return np.mean([tab[f"cr_{el}"] for el in ELEMENTS], axis=0)


def crossing(t, y, thresh):
    hit = np.flatnonzero(y >= thresh)
    if len(hit) == 0:
        return np.nan
    i = int(hit[0])
    if i == 0:
        return float(t[0])
    y0, y1, t0, t1 = y[i - 1], y[i], t[i - 1], t[i]
    if y1 == y0:
        return float(t1)
    return float(t0 + (thresh - y0) / (y1 - y0) * (t1 - t0))


def stage_times(contact, ent_rows):
    t_wet, t_sep = [], []
    for T in TEMPERATURES:
        tab = series(contact, SYSTEM, T)
        tw = crossing(tab["time_ps"], mean_rho(tab), 0.45)
        t_wet.append(7.5e3 if not np.isfinite(tw) else tw)
        te, _, sc = entropy_series(ent_rows, T)
        imin = int(np.argmin(sc))
        rise = np.flatnonzero((te >= te[imin]) & (sc >= sc[imin] + 0.35))
        t_sep.append(9e3 if len(rise) == 0 else float(te[rise[0]]))
    return np.asarray(TEMPERATURES, float), np.asarray(t_wet), np.asarray(t_sep)


def snap_path(temp_k, t_ns):
    if abs(t_ns - 0.25) < 1e-9:
        tag = "0.25"
    elif abs(t_ns - 1) < 1e-9:
        tag = "1.0"
    else:
        tag = "5.0"
    return os.path.join(PNG_ROOT, SYSTEM, f"{temp_k}K", f"{tag}ns.png")


def _core_mask(body):
    mask = body.mean(axis=2) < 240
    lab, nlab = cc_label(mask)
    if nlab == 0:
        return mask
    counts = np.bincount(lab.ravel())
    counts[0] = 0
    thresh = max(int(counts.max() * 0.06), 400)
    keep = np.isin(lab, np.flatnonzero(counts >= thresh))
    return keep


def content_bbox(path):
    im = np.asarray(Image.open(path).convert("RGB")).copy()
    h, w = im.shape[:2]
    im[:int(0.13 * h), :] = 255
    y0c, y1c = int(0.13 * h), int(0.97 * h)
    x0c, x1c = int(0.03 * w), int(0.97 * w)
    body = im[y0c:y1c, x0c:x1c]
    keep = _core_mask(body)
    ys, xs = np.where(keep)
    if len(xs) < 30:
        return im, (0, im.shape[0], 0, im.shape[1])
    pad = 36
    y0 = y0c + max(0, int(ys.min()) - pad)
    y1 = y0c + min(body.shape[0], int(ys.max()) + pad)
    x0 = x0c + max(0, int(xs.min()) - pad)
    x1 = x0c + min(body.shape[1], int(xs.max()) + pad)
    return im, (y0, y1, x0, x1)


def crop_to(im, box, size=440):
    y0, y1, x0, x1 = box
    cy, cx = (y0 + y1) / 2.0, (x0 + x1) / 2.0
    half = max(y1 - y0, x1 - x0) / 2.0
    y0 = int(max(0, cy - half))
    y1 = int(min(im.shape[0], cy + half))
    x0 = int(max(0, cx - half))
    x1 = int(min(im.shape[1], cx + half))
    crop = im[y0:y1, x0:x1]
    side = max(crop.shape[0], crop.shape[1], 1)
    canvas = np.full((side, side, 3), 255, np.uint8)
    oy = (side - crop.shape[0]) // 2
    ox = (side - crop.shape[1]) // 2
    canvas[oy:oy + crop.shape[0], ox:ox + crop.shape[1]] = crop
    return np.asarray(Image.fromarray(canvas).resize(
        (size, size), Image.Resampling.LANCZOS))


def snap_axes(fig, spec, items, title, letter, stage_colors=None):
    n = len(items)
    inner = GridSpecFromSubplotSpec(
        2, n, subplot_spec=spec, height_ratios=[0.12, 1.0],
        hspace=0.05, wspace=0.045,
    )
    tax = fig.add_subplot(inner[0, :])
    tax.set_axis_off()
    tax.text(0.0, 0.25, letter, fontsize=8, fontweight="bold",
             transform=tax.transAxes, va="center", ha="left")
    tax.text(0.05, 0.25, title, fontsize=7, color="0.25",
             transform=tax.transAxes, va="center", ha="left")
    axes = []
    for i, item in enumerate(items):
        axi = fig.add_subplot(inner[1, i])
        axi.imshow(item["img"], interpolation="bilinear")
        axi.set_xticks([])
        axi.set_yticks([])
        for sp in axi.spines.values():
            sp.set_linewidth(0.45)
            sp.set_color("0.55")
        axi.set_xlabel(item["cap"], fontsize=6.2, labelpad=4, color="0.15")
        if stage_colors is not None:
            axi.plot([0.22, 0.78], [-0.07, -0.07], transform=axi.transAxes,
                     color=stage_colors[i], lw=1.8, solid_capstyle="round",
                     clip_on=False)
        axes.append(axi)
    return axes


def build():
    apply_style()
    os.makedirs(OUT, exist_ok=True)
    contact = load_results()
    ent = load_entropy()
    T_pts, t_wet, t_sep = stage_times(contact, ent)

    T_fine = np.linspace(1475, 3025, 240)
    tw = np.clip(PchipInterpolator(T_pts, t_wet)(T_fine) / 1000.0, T_MIN, T_MAX)
    ts = np.clip(PchipInterpolator(T_pts, t_sep)(T_fine) / 1000.0, T_MIN, T_MAX)

    # shared camera for the 2000 K time series (widest = 0.25 ns)
    ims_c, boxes_c = zip(*[content_bbox(snap_path(2000, t))
                           for t in (0.25, 1.0, 5.0)])
    y0 = min(b[0] for b in boxes_c)
    y1 = max(b[1] for b in boxes_c)
    x0 = min(b[2] for b in boxes_c)
    x1 = max(b[3] for b in boxes_c)
    shared = (y0, y1, x0, x1)
    imgs_time = [crop_to(im, shared) for im in ims_c]

    imgs_temp = []
    for T in (1500, 2000, 2500, 3000):
        im, box = content_bbox(snap_path(T, 5.0))
        imgs_temp.append(crop_to(im, box))

    fig = plt.figure(figsize=(183 * MM, 158 * MM))
    gs = GridSpec(
        3, 2, figure=fig,
        height_ratios=[1.05, 1.22, 1.02],
        width_ratios=[1.12, 1.0],
        hspace=0.28, wspace=0.22,
        left=0.065, right=0.985, top=0.93, bottom=0.045,
    )

    # ── a ──────────────────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, :])
    ax.fill_betweenx(T_fine, T_MIN, tw, color=STAGE_C[1], lw=0, zorder=0)
    ax.fill_betweenx(T_fine, tw, ts, color=STAGE_C[2], lw=0, zorder=0)
    ax.fill_betweenx(T_fine, ts, T_MAX, color=STAGE_C[3], lw=0, zorder=0)
    ax.plot(tw, T_fine, color="0.28", lw=0.9, ls=(0, (3.2, 1.6)), zorder=3)
    ax.plot(ts, T_fine, color="0.18", lw=1.05, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(T_MIN, T_MAX)
    ax.set_ylim(1475, 3025)
    ax.set_yticks(list(TEMPERATURES))
    ax.set_xlabel("Time (ns)", labelpad=2)
    ax.set_ylabel("Temperature (K)", labelpad=4)
    style(ax, box=True)
    panel(ax, "a", x=-0.055, y=1.12)

    ax.plot([0.22, 5.0], [2000, 2000], color="0.12", lw=0.85, zorder=4)
    ax.plot([5.0, 5.0], [1500, 3000], color="0.12", lw=0.85, zorder=4)
    ax.annotate("", xy=(5.0, 2000), xytext=(4.15, 2000),
                arrowprops=dict(arrowstyle="-|>", color="0.12", lw=0.85,
                                mutation_scale=7), zorder=4)
    ax.annotate("", xy=(5.0, 3000), xytext=(5.0, 2780),
                arrowprops=dict(arrowstyle="-|>", color="0.12", lw=0.85,
                                mutation_scale=7), zorder=4)

    stops = [(0.25, 2000, "1"), (1.0, 2000, "2"), (5.0, 2000, "3"),
             (5.0, 1500, "4"), (5.0, 2500, "5"), (5.0, 3000, "6")]
    for t, T, n in stops:
        ax.scatter([t], [T], s=86, facecolor="white", edgecolor="0.12",
                   linewidths=0.7, zorder=6)
        ax.text(t, T, n, ha="center", va="center", fontsize=6.2,
                fontweight="bold", zorder=7, color="0.08")

    ax.legend(
        handles=[
            Patch(facecolor=STAGE_C[1], edgecolor=STAGE_EDGE[1], lw=0.4,
                  label="I  wetting"),
            Patch(facecolor=STAGE_C[2], edgecolor=STAGE_EDGE[2], lw=0.4,
                  label="II  diffusion"),
            Patch(facecolor=STAGE_C[3], edgecolor=STAGE_EDGE[3], lw=0.4,
                  label="III  phase separation"),
        ],
        loc="lower left", frameon=False, ncol=3, handlelength=1.15,
        handleheight=0.8, columnspacing=1.4, borderpad=0,
        bbox_to_anchor=(0.08, 1.02),
    )

    # ── b ──────────────────────────────────────────────────────────
    sub = GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1, 0], hspace=0.06)
    axes_b = [fig.add_subplot(sub[i]) for i in range(3)]
    ylabels = [
        "Contact ratio  $\\rho$",
        "$S_\\mathrm{atom}$  (nats)",
        "$S_\\mathrm{config}$  (nats)",
    ]
    for T in T_SHOW:
        tab = series(contact, SYSTEM, T)
        t = np.clip(tab["time_ps"] / 1000.0, 2e-4, None)
        te, sa, sc = entropy_series(ent, T)
        te = np.clip(te / 1000.0, 2e-4, None)
        col = T_COLOR[T]
        axes_b[0].plot(t, mean_rho(tab), color=col, lw=1.2, label=f"{T} K",
                       solid_capstyle="round")
        axes_b[1].plot(te, sa, color=col, lw=1.2, solid_capstyle="round")
        axes_b[2].plot(te, sc, color=col, lw=1.2, solid_capstyle="round")
    for axb, ylab in zip(axes_b, ylabels):
        axb.set_xscale("log")
        axb.set_xlim(T_MIN, T_MAX)
        axb.set_ylabel(ylab, labelpad=3)
        style(axb)
        axb.axvline(1.0, color="0.78", lw=0.45, ls=":", zorder=1)
    axes_b[0].set_ylim(-0.02, 1.02)
    axes_b[0].set_yticks([0, 0.5, 1.0])
    axes_b[0].legend(frameon=False, loc="lower right", handlelength=1.25,
                     labelspacing=0.12, borderpad=0.15, ncol=1)
    axes_b[0].tick_params(labelbottom=False)
    axes_b[1].tick_params(labelbottom=False)
    axes_b[2].set_xlabel("Time (ns)", labelpad=2)
    panel(axes_b[0], "b", x=-0.16, y=1.08)

    # ── c ──────────────────────────────────────────────────────────
    snap_axes(
        fig, gs[1, 1],
        [
            {"img": imgs_time[0], "n": "1", "cap": "1   ·   0.25 ns"},
            {"img": imgs_time[1], "n": "2", "cap": "2   ·   1 ns"},
            {"img": imgs_time[2], "n": "3", "cap": "3   ·   5 ns"},
        ],
        "Time control at 2000 K",
        "c",
        stage_colors=[STAGE_EDGE[1], STAGE_EDGE[2], STAGE_EDGE[2]],
    )

    # ── d ──────────────────────────────────────────────────────────
    snap_axes(
        fig, gs[2, :],
        [
            {"img": imgs_temp[0], "n": "4", "cap": "4   ·   1500 K"},
            {"img": imgs_temp[1], "n": "3", "cap": "3   ·   2000 K"},
            {"img": imgs_temp[2], "n": "5", "cap": "5   ·   2500 K"},
            {"img": imgs_temp[3], "n": "6", "cap": "6   ·   3000 K"},
        ],
        "Temperature control at 5 ns",
        "d",
        stage_colors=[STAGE_EDGE[1], STAGE_EDGE[2], STAGE_EDGE[3], STAGE_EDGE[3]],
    )

    stem = os.path.join(OUT, "Fig4_stopping_point")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}")
        print("Saved", f"{stem}.{ext}")
    plt.close(fig)


if __name__ == "__main__":
    build()
