#!/usr/bin/env python3
"""
Metal-on-Ni wetting: per-element contact ratio and contact-area proxy.

Systems are a central NiO particle with Fe / Zn / Cr / Ru oxide particles.
Only metal–metal contacts are counted (O excluded).  An X atom (X = Fe, Zn,
Cr, Ru) is in contact with Ni if it has at least one Ni neighbour within
MET_CUT.  Contact ratio is n_contact / N_X; contact area is n_contact times
a close-packed atomic area.

When X dissolves into the Ni particle (high Ni coordination) the geometric
interface disappears and the contact-area number is no longer a wetting
metric — those frames are flagged as fused.

Time sampling is logarithmic and denser in 0–1 ns.
"""

from __future__ import annotations

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor
from collections import Counter

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

import traj_entropy as TE
import run_entropy as RE

MET_CUT = 3.0  # A, same metal-metal cutoff as the entropy analysis
ELEMENTS = ("Fe", "Zn", "Cr", "Ru")
# Close-packed area per contacting atom, d_nn = 2.50 A
AREA_PER_ATOM = 0.5 * np.sqrt(3.0) * (2.50 ** 2)  # ~5.41 A^2
ALLOY_CN = 4  # CN_Ni >= 4: embedded / alloyed, not a surface contact
DUMP_INTERVAL_PS = 0.4
MAX_PS = 5000.0

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT, "contact_analysis")
DATA_DIR = os.path.join(OUT_DIR, "data")
MODEL_DIR = os.path.join(ROOT, "FIG4-run", "models")


def log_frame_picks(n_frames: int, max_ps: float = MAX_PS) -> np.ndarray:
    """Log-spaced dump indices, denser below 1 ns."""
    t_lo = DUMP_INTERVAL_PS
    t_hi = min(max_ps, n_frames * DUMP_INTERVAL_PS)
    t1 = np.logspace(np.log10(t_lo), np.log10(min(1000.0, t_hi)), 90)
    if t_hi > 1000.0:
        t2 = np.logspace(np.log10(1000.0), np.log10(t_hi), 50)
        ts = np.concatenate([t1, t2])
    else:
        ts = t1
    idx = np.clip(np.round(ts / DUMP_INTERVAL_PS).astype(int) - 1, 0, n_frames - 1)
    return np.unique(idx)


def trajectory_from_cache(path: str, system: str, temp_k: int) -> TE.Trajectory:
    cache = RE.offsets_cache(system, temp_k)
    payload = np.load(cache)
    size, natoms, offsets = int(payload[0]), int(payload[1]), payload[2:]
    if int(os.path.getsize(path)) != size:
        raise RuntimeError(f"stale offset cache for {system} {temp_k}K")
    traj = TE.Trajectory.__new__(TE.Trajectory)
    traj.path = path
    traj.natoms = natoms
    traj._offsets = offsets
    traj.n_frames = max(len(offsets) - 1, 0)
    traj._fh = open(path, "rb", buffering=0)
    return traj


def read_xyz(path: str):
    with open(path) as fh:
        natoms = int(fh.readline())
        comment = fh.readline()
        species, pos = [], []
        for _ in range(natoms):
            parts = fh.readline().split()
            species.append(parts[0])
            pos.append((float(parts[1]), float(parts[2]), float(parts[3])))
    lat = TE._LATTICE_RE.search(comment.encode())
    cell = np.asarray(lat.group(1).split(), dtype=np.float64).reshape(3, 3)
    box = np.diag(cell).copy()
    return np.asarray(species, dtype="U3"), np.asarray(pos, dtype=np.float64), box


def frame_contact(species: np.ndarray, positions: np.ndarray, box: np.ndarray):
    """Return per-element contact metrics for one configuration."""
    wrapped = np.mod(positions, box)
    is_o = species == "O"
    is_ni = species == "Ni"
    metal = np.flatnonzero(~is_o)
    n_metal = len(metal)
    mpos = wrapped[metal]
    mspec = species[metal]
    m_is_ni = is_ni[metal]

    tree = cKDTree(mpos, boxsize=box)
    pairs = tree.query_pairs(MET_CUT, output_type="ndarray")

    cn_ni = np.zeros(n_metal, dtype=np.int32)
    if len(pairs):
        a, b = pairs[:, 0], pairs[:, 1]
        cn_ni += np.bincount(a[m_is_ni[b]], minlength=n_metal).astype(np.int32)
        cn_ni += np.bincount(b[m_is_ni[a]], minlength=n_metal).astype(np.int32)

    if len(pairs):
        graph = coo_matrix(
            (np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])),
            shape=(n_metal, n_metal),
        ).tocsr()
    else:
        graph = coo_matrix((n_metal, n_metal))
    n_clusters, labels = connected_components(graph, directed=False)
    ni_local = np.flatnonzero(m_is_ni)
    if len(ni_local):
        ni_lab, ni_cnt = np.unique(labels[ni_local], return_counts=True)
        ni_core = int(ni_lab[int(np.argmax(ni_cnt))])
    else:
        ni_core = -1
    in_ni_core = labels == ni_core if ni_core >= 0 else np.zeros(n_metal, dtype=bool)

    ni_com = wrapped[is_ni].mean(axis=0) if is_ni.any() else np.zeros(3)

    out = {}
    for el in ELEMENTS:
        sel = mspec == el
        n = int(sel.sum())
        if n == 0:
            out[el] = dict(n=0, n_contact=0, n_pairs=0, contact_ratio=0.0,
                           area_A2=0.0, mean_cn_ni=0.0, alloyed_frac=0.0,
                           attached_frac=0.0, com_dist=np.nan)
            continue
        cn = cn_ni[sel]
        n_contact = int((cn >= 1).sum())
        n_pairs = int(cn.sum())
        out[el] = dict(
            n=n,
            n_contact=n_contact,
            n_pairs=n_pairs,
            contact_ratio=n_contact / n,
            area_A2=n_contact * AREA_PER_ATOM,
            mean_cn_ni=float(cn.mean()),
            alloyed_frac=float((cn >= ALLOY_CN).sum()) / n,
            attached_frac=float(in_ni_core[sel].sum()) / n,
            com_dist=float(np.linalg.norm(mpos[sel].mean(axis=0) - ni_com)),
        )
    return out


ROW_DTYPE = [
    ("frame", "i8"), ("time_ps", "f8"),
    ("n_Fe", "i8"), ("n_Zn", "i8"), ("n_Cr", "i8"), ("n_Ru", "i8"),
]
for _el in ELEMENTS:
    ROW_DTYPE += [
        (f"cr_{_el}", "f8"), (f"nc_{_el}", "i8"), (f"np_{_el}", "i8"),
        (f"area_{_el}", "f8"), (f"cn_{_el}", "f8"),
        (f"alloy_{_el}", "f8"), (f"att_{_el}", "f8"), (f"dist_{_el}", "f8"),
    ]


def pack_row(frame, time_ps, metrics):
    rec = [int(frame), float(time_ps)]
    rec += [metrics[el]["n"] for el in ELEMENTS]
    for el in ELEMENTS:
        m = metrics[el]
        rec += [m["contact_ratio"], m["n_contact"], m["n_pairs"],
                m["area_A2"], m["mean_cn_ni"], m["alloyed_frac"],
                m["attached_frac"], m["com_dist"]]
    return tuple(rec)


def analyze_model(system: str):
    path = os.path.join(MODEL_DIR, f"{system}.xyz")
    species, pos, box = read_xyz(path)
    return pack_row(-1, 0.0, frame_contact(species, pos, box))


def worker(task):
    system, temp_k, out_path = task
    path = RE.dump_path(system, temp_k)
    t0 = time.time()
    traj = trajectory_from_cache(path, system, temp_k)
    n_dump = traj.n_frames
    n_window = min(n_dump, int(round(MAX_PS / DUMP_INTERVAL_PS)))
    picks = log_frame_picks(n_window)
    rows = [analyze_model(system)]
    try:
        for fidx in picks:
            species, positions, box, time_fs = traj.frame(int(fidx))
            metrics = frame_contact(species, positions, box)
            rows.append(pack_row(int(fidx), time_fs / 1000.0, metrics))
    finally:
        traj.close()
    table = np.array(rows, dtype=ROW_DTYPE)
    np.savez(
        out_path,
        table=table,
        n_frames_total=n_dump,
        n_frames_window=n_window,
        temp_k=temp_k,
        system=system,
        met_cut=MET_CUT,
        area_per_atom=AREA_PER_ATOM,
        alloy_cn=ALLOY_CN,
    )
    return {
        "system": system,
        "temp_k": temp_k,
        "n_dump": n_dump,
        "n_analyzed": len(table),
        "t_max_ps": float(table["time_ps"].max()),
        "wall_s": time.time() - t0,
    }


def run_all(workers: int = 10, force: bool = False):
    os.makedirs(DATA_DIR, exist_ok=True)
    tasks = [(s, t, os.path.join(DATA_DIR, f"{s}_{t}K.npz"))
             for s in RE.SYSTEMS for t in RE.TEMPERATURES]
    if not force:
        tasks = [t for t in tasks if not os.path.isfile(t[2])]
    if not tasks:
        print("All trajectories already analysed (use --force to redo).")
        return
    print(f"Analysing {len(tasks)} trajectories with {workers} workers "
          f"(log-spaced frames, first {MAX_PS:g} ps)")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for info in pool.map(worker, tasks):
            print("  {system} {temp_k}K: {n_dump} dump frames, "
                  "{n_analyzed} analysed (t_max={t_max_ps:.1f} ps), "
                  "{wall_s:.0f}s".format(**info), flush=True)


def load_results():
    results = {}
    for system in RE.SYSTEMS:
        for temp_k in RE.TEMPERATURES:
            path = os.path.join(DATA_DIR, f"{system}_{temp_k}K.npz")
            if os.path.isfile(path):
                results[(system, temp_k)] = np.load(path, allow_pickle=True)
    return results


def series(results, system, temp_k):
    payload = results.get((system, temp_k))
    if payload is None:
        return None
    table = payload["table"]
    return table[np.argsort(table["time_ps"])]


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

ELEM_COLOR = {
    "Fe": "#c45c26",
    "Zn": "#5a7aa5",
    "Cr": "#2a8a4a",
    "Ru": "#8b3d88",
}
FIG_W, FIG_H = 7.0, 3.15


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "lines.linewidth": 1.15,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", direction="in", length=3, width=0.6)


def temp_color(temp_k):
    norm = Normalize(vmin=min(RE.TEMPERATURES), vmax=max(RE.TEMPERATURES))
    return cm.coolwarm(norm(temp_k))


def save(fig, stem):
    for ext in ("png", "svg"):
        path = os.path.join(OUT_DIR, f"{stem}.{ext}")
        fig.savefig(path)
        print(f"Saved: {path}")
    plt.close(fig)


def shade_wetting(ax):
    ax.axvspan(4e-4, 1.0, color="0.92", lw=0, zorder=0)
    ax.axvline(1.0, color="0.55", lw=0.6, ls=":", zorder=1)


def plot_ratio_by_temp(results, system):
    """One panel per T: 4 elements, contact ratio vs log time."""
    apply_style()
    fig, axes = plt.subplots(1, 5, figsize=(11.2, 2.55), sharey=True)
    for ax, temp_k in zip(axes, RE.TEMPERATURES):
        tab = series(results, system, temp_k)
        shade_wetting(ax)
        t_ns = tab["time_ps"] / 1000.0
        t_ns = np.clip(t_ns, 3e-4, None)
        for el in ELEMENTS:
            cr = tab[f"cr_{el}"]
            alloy = tab[f"alloy_{el}"]
            fused = alloy >= 0.50
            ax.plot(t_ns, cr, color=ELEM_COLOR[el], lw=1.2, label=el, zorder=3)
            if fused.any():
                ax.plot(t_ns[fused], cr[fused], color=ELEM_COLOR[el],
                        lw=1.2, ls="--", zorder=4)
        ax.set_xscale("log")
        ax.set_xlim(3e-4, 6)
        ax.set_ylim(-0.03, 1.05)
        ax.set_title(f"{temp_k} K", loc="left", pad=4)
        style_axes(ax)
        ax.set_xlabel("Time (ns)", labelpad=2)
    axes[0].set_ylabel("Contact ratio  $n_\\mathrm{c}/N_X$", labelpad=3)
    handles = [Line2D([0], [0], color=ELEM_COLOR[el], lw=1.4, label=el)
               for el in ELEMENTS]
    handles.append(Line2D([0], [0], color="0.3", lw=1.2, ls="--",
                          label="fused ($f_\\mathrm{alloy}\\geq0.5$)"))
    axes[-1].legend(handles=handles, frameon=False, loc="upper left",
                    handlelength=1.5, borderpad=0.2, labelspacing=0.2)
    fig.suptitle(f"{system}  |  metal–Ni contact ratio (log time)",
                 x=0.01, ha="left", fontsize=8, y=1.04)
    fig.tight_layout(w_pad=0.45)
    save(fig, f"{system}_contact_ratio_byT")


def plot_ratio_by_element(results, system):
    """One panel per element: 5 temperatures."""
    apply_style()
    fig, axes = plt.subplots(1, 4, figsize=(10.4, 2.55), sharey=True)
    for ax, el in zip(axes, ELEMENTS):
        shade_wetting(ax)
        for temp_k in RE.TEMPERATURES:
            tab = series(results, system, temp_k)
            t_ns = np.clip(tab["time_ps"] / 1000.0, 3e-4, None)
            ax.plot(t_ns, tab[f"cr_{el}"], color=temp_color(temp_k),
                    lw=1.15, label=f"{temp_k} K", zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(3e-4, 6)
        ax.set_ylim(-0.03, 1.05)
        ax.set_title(el, loc="left", pad=4)
        style_axes(ax)
        ax.set_xlabel("Time (ns)", labelpad=2)
    axes[0].set_ylabel("Contact ratio  $n_\\mathrm{c}/N_X$", labelpad=3)
    axes[-1].legend(frameon=False, loc="upper left", handlelength=1.3,
                    borderpad=0.2, labelspacing=0.18, fontsize=6)
    fig.suptitle(f"{system}  |  contact ratio by element",
                 x=0.01, ha="left", fontsize=8, y=1.04)
    fig.tight_layout(w_pad=0.45)
    save(fig, f"{system}_contact_ratio_byElem")


def plot_area(results, system):
    apply_style()
    fig, axes = plt.subplots(1, 5, figsize=(11.2, 2.55), sharey=True)
    ymax = 0.0
    for ax, temp_k in zip(axes, RE.TEMPERATURES):
        tab = series(results, system, temp_k)
        shade_wetting(ax)
        t_ns = np.clip(tab["time_ps"] / 1000.0, 3e-4, None)
        for el in ELEMENTS:
            area = tab[f"area_{el}"]
            alloy = tab[f"alloy_{el}"]
            solid = alloy < 0.50
            ax.plot(t_ns[solid], area[solid], color=ELEM_COLOR[el],
                    lw=1.2, label=el, zorder=3)
            if (~solid).any():
                ax.plot(t_ns[~solid], area[~solid], color=ELEM_COLOR[el],
                        lw=1.15, ls="--", alpha=0.7, zorder=3)
            ymax = max(ymax, float(area.max()))
        ax.set_xscale("log")
        ax.set_xlim(3e-4, 6)
        ax.set_title(f"{temp_k} K", loc="left", pad=4)
        style_axes(ax)
        ax.set_xlabel("Time (ns)", labelpad=2)
    for ax in axes:
        ax.set_ylim(0, ymax * 1.08 if ymax > 0 else 1)
    axes[0].set_ylabel("Contact area  $A_c$  (Å$^2$)", labelpad=3)
    handles = [Line2D([0], [0], color=ELEM_COLOR[el], lw=1.4, label=el)
               for el in ELEMENTS]
    axes[-1].legend(handles=handles, frameon=False, loc="upper left",
                    handlelength=1.4, borderpad=0.2, labelspacing=0.2)
    fig.suptitle(
        f"{system}  |  metal–Ni contact area  "
        r"($A_c = n_c \times 5.41\,\mathrm{\AA}^2$; dashed = fused)",
        x=0.01, ha="left", fontsize=8, y=1.04,
    )
    fig.tight_layout(w_pad=0.45)
    save(fig, f"{system}_contact_area_byT")


def plot_overview(results):
    apply_style()
    fig, axes = plt.subplots(2, 5, figsize=(11.2, 4.6), sharex=True, sharey=True)
    for row, system in enumerate(RE.SYSTEMS):
        short = "NiO-6" if "NiO-6" in system else "NiO-8"
        for col, temp_k in enumerate(RE.TEMPERATURES):
            ax = axes[row, col]
            tab = series(results, system, temp_k)
            shade_wetting(ax)
            t_ns = np.clip(tab["time_ps"] / 1000.0, 3e-4, None)
            for el in ELEMENTS:
                ax.plot(t_ns, tab[f"cr_{el}"], color=ELEM_COLOR[el],
                        lw=1.1, zorder=3)
            ax.set_xscale("log")
            ax.set_xlim(3e-4, 6)
            ax.set_ylim(-0.03, 1.05)
            style_axes(ax)
            if row == 0:
                ax.set_title(f"{temp_k} K", loc="left", pad=3)
            if col == 0:
                ax.set_ylabel(f"{short}\ncontact ratio", labelpad=3)
            if row == 1:
                ax.set_xlabel("Time (ns)", labelpad=2)
    handles = [Line2D([0], [0], color=ELEM_COLOR[el], lw=1.4, label=el)
               for el in ELEMENTS]
    axes[0, -1].legend(handles=handles, frameon=False, loc="upper left",
                       handlelength=1.3, borderpad=0.15, labelspacing=0.15)
    fig.suptitle("Metal–Ni contact ratio  |  grey band = wetting window 0–1 ns",
                 x=0.01, ha="left", fontsize=8, y=1.01)
    fig.tight_layout(h_pad=0.6, w_pad=0.4)
    save(fig, "contact_ratio_overview")


def crossing_time(time_ps, y, thresh):
    """First time y reaches thresh; nan if never."""
    hit = np.flatnonzero(y >= thresh)
    if len(hit) == 0:
        return np.nan
    i = int(hit[0])
    if i == 0:
        return float(time_ps[0])
    y0, y1 = y[i - 1], y[i]
    t0, t1 = time_ps[i - 1], time_ps[i]
    if y1 == y0:
        return float(t1)
    frac = (thresh - y0) / (y1 - y0)
    return float(t0 + frac * (t1 - t0))


def value_at(time_ps, y, t_target):
    if t_target <= time_ps[0]:
        return float(y[0])
    if t_target >= time_ps[-1]:
        return float(y[-1])
    i = int(np.searchsorted(time_ps, t_target))
    t0, t1 = time_ps[i - 1], time_ps[i]
    w = (t_target - t0) / (t1 - t0) if t1 != t0 else 1.0
    return float(y[i - 1] * (1 - w) + y[i] * w)


def summary_table(results):
    rows = []
    header = [
        "system", "T_K", "element", "N_X",
        "t10_ps", "t50_ps",
        "cr_1ns", "cr_5ns", "area_1ns", "area_5ns",
        "alloy_1ns", "alloy_5ns", "fused",
    ]
    rows.append(header)
    for system in RE.SYSTEMS:
        for temp_k in RE.TEMPERATURES:
            tab = series(results, system, temp_k)
            t = tab["time_ps"]
            for el in ELEMENTS:
                cr = tab[f"cr_{el}"]
                alloy = tab[f"alloy_{el}"]
                fused = bool(alloy[-1] >= 0.50)
                rows.append([
                    system, temp_k, el, int(tab[f"n_{el}"][0]),
                    round(crossing_time(t, cr, 0.10), 2),
                    round(crossing_time(t, cr, 0.50), 2),
                    round(value_at(t, cr, 1000.0), 4),
                    round(value_at(t, cr, 5000.0), 4),
                    round(value_at(t, tab[f"area_{el}"], 1000.0), 1),
                    round(value_at(t, tab[f"area_{el}"], 5000.0), 1),
                    round(value_at(t, alloy, 1000.0), 4),
                    round(value_at(t, alloy, 5000.0), 4),
                    int(fused),
                ])
    return rows


def write_long_csv(results):
    path = os.path.join(OUT_DIR, "contact_raw_long.csv")
    fields = ["system", "temperature_K", "frame_index", "time_ps", "element",
              "N_X", "n_contact", "n_pairs", "contact_ratio", "area_A2",
              "mean_cn_ni", "alloyed_frac", "attached_frac", "com_dist_A"]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join(fields) + "\n")
        for system in RE.SYSTEMS:
            for temp_k in RE.TEMPERATURES:
                tab = series(results, system, temp_k)
                for r in tab:
                    for el in ELEMENTS:
                        fh.write(
                            f"{system},{temp_k},{int(r['frame'])},"
                            f"{float(r['time_ps']):.6f},{el},"
                            f"{int(r[f'n_{el}'])},{int(r[f'nc_{el}'])},"
                            f"{int(r[f'np_{el}'])},{float(r[f'cr_{el}']):.6f},"
                            f"{float(r[f'area_{el}']):.4f},"
                            f"{float(r[f'cn_{el}']):.4f},"
                            f"{float(r[f'alloy_{el}']):.6f},"
                            f"{float(r[f'att_{el}']):.6f},"
                            f"{float(r[f'dist_{el}']):.4f}\n"
                        )
    print(f"Saved: {path}")


def write_summary_csv(results):
    rows = summary_table(results)
    path = os.path.join(OUT_DIR, "contact_summary.csv")
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(",".join(str(x) for x in row) + "\n")
    print(f"Saved: {path}")
    return rows


def write_notes(rows):
    path = os.path.join(OUT_DIR, "NOTES.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(
            "Metal–Ni wetting contact analysis (FIG4 GPUMD / NEP89)\n\n"
            f"Cutoff: metal–metal < {MET_CUT} A, oxygen excluded.\n"
            "Contact ratio of element X: fraction of X metal atoms with ≥1 Ni neighbour.\n"
            f"Contact area A_c = n_contact × {AREA_PER_ATOM:.2f} Å² "
            "(close-packed atomic area, d_nn = 2.50 Å).\n"
            f"Fused flag: alloyed_frac = N(CN_Ni ≥ {ALLOY_CN}) / N_X ≥ 0.5 "
            "at 5 ns — the X/Ni interface is no longer a well-defined "
            "wetting contact area.\n"
            "Grey band in the figures is the 0–1 ns wetting window; "
            "the time axis is logarithmic.\n\n"
        )
        fh.write(
            "| System | T (K) | El | t(10%) ps | t(50%) ps | "
            "ρ(1 ns) | ρ(5 ns) | A_c(1 ns) Å² | fused |\n"
            "|---|---|---|---|---|---|---|---|---|\n"
        )
        for row in rows[1:]:
            t10 = row[4] if row[4] == row[4] else "—"
            t50 = row[5] if row[5] == row[5] else "—"
            fused = "yes" if row[12] else "no"
            fh.write(
                f"| {row[0]} | {row[1]} | {row[2]} | {t10} | {t50} | "
                f"{row[6]:.3f} | {row[7]:.3f} | {row[8]:.0f} | {fused} |\n"
            )
    print(f"Saved: {path}")


def build(results):
    os.makedirs(OUT_DIR, exist_ok=True)
    if not results:
        raise SystemExit("No results found; run the analysis first.")
    write_long_csv(results)
    rows = write_summary_csv(results)
    write_notes(rows)
    for system in RE.SYSTEMS:
        plot_ratio_by_temp(results, system)
        plot_ratio_by_element(results, system)
        plot_area(results, system)
    plot_overview(results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    if not args.plot_only:
        run_all(workers=args.workers, force=args.force)
    build(load_results())


if __name__ == "__main__":
    main()
