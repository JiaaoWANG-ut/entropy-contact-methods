#!/usr/bin/env python3
"""
Driver: atomic and configurational entropy vs time for all FIG4 GPUMD runs.

Systems   : NiO-6-Fe-Zn-Cr-Ru, NiO-8-Fe-Zn-Cr-Ru
Temperatures: 1500, 1800, 2000, 2500, 3000 K
Per trajectory 400 uniformly spaced frames are analysed (see traj_entropy.py).

Outputs (in entropy_analysis/):
  data/<system>_<T>K.npz          per-trajectory results + frame offset cache
  entropy_raw_long.csv            tidy raw data (one row per frame)
  <system>_Satom_wide.csv         plot-ready wide tables
  <system>_Sconfig_wide.csv
  entropy_tables.xlsx             all tables in one workbook
  <system>_Satom_time.png/.svg    figures
  <system>_Sconfig_time.png/.svg
  <system>_entropy_dual.png/.svg
  entropy_overview.png/.svg
"""

from __future__ import annotations

import argparse
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import traj_entropy as TE

ROOT = os.path.dirname(os.path.abspath(__file__))
RUN_ROOT = os.path.join(ROOT, "FIG4-run")
SYSTEMS = ("NiO-6-Fe-Zn-Cr-Ru", "NiO-8-Fe-Zn-Cr-Ru")
TEMPERATURES = (1500, 1800, 2000, 2500, 3000)
N_SAMPLE = 400
DUMP_INTERVAL_PS = 0.4  # dump_xyz every 400 steps at dt = 1 fs

# Frame byte-offset caches are independent of the analysis window and are
# therefore shared by all runs.
INDEX_DIR = os.path.join(ROOT, "entropy_frame_index")

# Set by configure(): analysis window and output location.
OUT_DIR = os.path.join(ROOT, "entropy_analysis")
DATA_DIR = os.path.join(OUT_DIR, "data")
MAX_PS = None  # None = use each trajectory in full
WINDOW_LABEL = "full trajectory"


def configure(tag: str = "", max_ps: float | None = None):
    """Select the analysis window and the matching output directory."""
    global OUT_DIR, DATA_DIR, MAX_PS, WINDOW_LABEL
    name = "entropy_analysis" + (f"_{tag}" if tag else "")
    OUT_DIR = os.path.join(ROOT, name)
    DATA_DIR = os.path.join(OUT_DIR, "data")
    MAX_PS = max_ps
    WINDOW_LABEL = ("full trajectory" if max_ps is None
                    else f"first {max_ps / 1000:g} ns")
    return OUT_DIR


def dump_path(system: str, temp_k: int) -> str:
    path = os.path.join(RUN_ROOT, system, f"{temp_k}K", "dump.xyz")
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    return path


def cache_path(system: str, temp_k: int) -> str:
    return os.path.join(DATA_DIR, f"{system}_{temp_k}K.npz")


def offsets_cache(system: str, temp_k: int) -> str:
    return os.path.join(INDEX_DIR, f"{system}_{temp_k}K_offsets.npy")


def load_or_build_offsets(system: str, temp_k: int, path: str):
    cache = offsets_cache(system, temp_k)
    size = os.path.getsize(path)
    if os.path.isfile(cache):
        payload = np.load(cache)
        if int(payload[0]) == size:
            return int(payload[1]), payload[2:]
    natoms, offsets = TE.index_frames(path)
    os.makedirs(INDEX_DIR, exist_ok=True)
    np.save(cache, np.concatenate(([size, natoms], offsets)).astype(np.int64))
    return natoms, offsets


def worker(task):
    system, temp_k, max_ps, out_path = task
    path = dump_path(system, temp_k)
    t0 = time.time()
    natoms, offsets = load_or_build_offsets(system, temp_k, path)
    n_frames_dump = len(offsets) - 1
    t_index = time.time() - t0

    n_frames = n_frames_dump
    if max_ps is not None:
        # frame i carries Time = (i + 1) * DUMP_INTERVAL_PS
        n_frames = min(n_frames, int(round(max_ps / DUMP_INTERVAL_PS)))
        if n_frames < 1:
            raise ValueError(f"{system} {temp_k}K shorter than {max_ps} ps")

    picks = TE.sample_indices(n_frames, N_SAMPLE)
    t1 = time.time()
    table = TE.analyze(path, picks)
    np.savez(
        out_path,
        table=table,
        natoms=natoms,
        n_frames_total=n_frames_dump,
        n_frames_window=n_frames,
        max_ps=-1.0 if max_ps is None else max_ps,
        temp_k=temp_k,
        system=system,
        dump_path=path,
        file_size=os.path.getsize(path),
    )
    return {
        "system": system,
        "temp_k": temp_k,
        "natoms": natoms,
        "n_frames_total": n_frames_dump,
        "n_frames_window": n_frames,
        "n_analyzed": len(table),
        "t_max_ps": float(table["time_ps"].max()),
        "index_s": t_index,
        "analyze_s": time.time() - t1,
    }


def run_all(workers: int = 10, force: bool = False):
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(INDEX_DIR, exist_ok=True)
    tasks = [(s, t, MAX_PS, cache_path(s, t))
             for s in SYSTEMS for t in TEMPERATURES]
    if not force:
        tasks = [t for t in tasks if not os.path.isfile(t[3])]
    if not tasks:
        print("All trajectories already analysed (use --force to redo).")
        return
    print(f"Analysing {len(tasks)} trajectories ({WINDOW_LABEL}) "
          f"with {workers} workers")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for info in pool.map(worker, tasks):
            print(
                "  {system} {temp_k}K: {n_frames_total} frames in dump, "
                "window {n_frames_window} frames ({t_max_ps:.0f} ps), "
                "{n_analyzed} analysed, index {index_s:.0f}s, "
                "entropy {analyze_s:.0f}s".format(**info),
                flush=True,
            )


def load_results():
    results = {}
    for system in SYSTEMS:
        for temp_k in TEMPERATURES:
            path = cache_path(system, temp_k)
            if os.path.isfile(path):
                payload = np.load(path, allow_pickle=True)
                results[(system, temp_k)] = payload
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--max-ps", type=float, default=None,
                        help="Analyse only the first MAX_PS ps of each run.")
    parser.add_argument("--tag", default="",
                        help="Suffix for the output directory "
                             "(entropy_analysis_<tag>).")
    args = parser.parse_args()

    configure(tag=args.tag, max_ps=args.max_ps)
    print(f"Output directory: {OUT_DIR}")

    if not args.plot_only:
        run_all(workers=args.workers, force=args.force)

    import entropy_report

    entropy_report.build(load_results(), out_dir=OUT_DIR,
                         window_label=WINDOW_LABEL)


if __name__ == "__main__":
    main()
