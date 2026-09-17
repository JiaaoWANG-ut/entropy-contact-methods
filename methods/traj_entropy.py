#!/usr/bin/env python3
"""
Entropy analysis of GPUMD extended-XYZ trajectories (numpy/scipy backend).

Two entropy definitions per frame:

  (1) Atomic entropy  S_atom  (states = discrete local environments)
        Each atom is labelled by the motif (species, CN_O, CN_M), where CN_O
        counts oxygen neighbours within O_CUT and CN_M counts metal
        (non-oxygen) neighbours within MET_CUT.
        S_atom = - sum_k p_k ln p_k   (nats),  p_k = N_k / N_atoms

  (2) Configurational entropy  S_config  (states = metal clusters)
        Metal (non-oxygen) atoms are decomposed into connected clusters
        (contact distance < MET_CUT, periodic).  Oxygen is excluded.
        S_config = - sum_c p_c ln p_c  (nats),  p_c = n_c / N_metal

The trajectory reader indexes frame byte offsets in a single streaming pass
(newline counting) so that arbitrary frames can be read by seeking, without
holding the multi-GB file in memory.
"""

from __future__ import annotations

import os
import re

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

O_CUT = 2.5    # A, first oxygen shell
MET_CUT = 3.0  # A, metal-metal contact / first-shell metal cutoff
OXYGEN = "O"

_LATTICE_RE = re.compile(rb'Lattice="([^"]+)"')
_TIME_RE = re.compile(rb"Time=([-\d.eE+]+)")


def index_frames(path: str, chunk: int = 1 << 26):
    """Return (natoms, offsets) where offsets[i] starts complete frame i.

    A trailing partial frame (job killed mid-write) is discarded.
    """
    with open(path, "rb") as fh:
        natoms = int(fh.readline())
    lines_per_frame = natoms + 2

    boundaries = [0]
    seen_lines = 0
    pos = 0
    with open(path, "rb", buffering=0) as fh:
        while True:
            buf = fh.read(chunk)
            if not buf:
                break
            nl = np.flatnonzero(np.frombuffer(buf, dtype=np.uint8) == 10)
            if nl.size:
                line_ids = seen_lines + np.arange(1, nl.size + 1)
                hits = nl[(line_ids % lines_per_frame) == 0]
                if hits.size:
                    boundaries.extend((pos + hits + 1).tolist())
                seen_lines += nl.size
            pos += len(buf)

    offsets = np.asarray(boundaries, dtype=np.int64)
    return natoms, offsets  # offsets[:-1] are starts of complete frames


class Trajectory:
    """Random-access reader for a fixed-composition extended-XYZ trajectory."""

    def __init__(self, path: str):
        self.path = path
        self.natoms, self._offsets = index_frames(path)
        self.n_frames = max(len(self._offsets) - 1, 0)
        self._fh = open(path, "rb", buffering=0)

    def close(self):
        if not self._fh.closed:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def frame(self, idx: int):
        """Return (species, positions, box_lengths, time_fs) for frame idx."""
        if not 0 <= idx < self.n_frames:
            raise IndexError(idx)
        start = int(self._offsets[idx])
        stop = int(self._offsets[idx + 1])
        self._fh.seek(start)
        raw = self._fh.read(stop - start)

        nl1 = raw.index(b"\n")
        nl2 = raw.index(b"\n", nl1 + 1)
        comment = raw[nl1 + 1 : nl2]

        tokens = raw[nl2 + 1 :].split()
        ncol = len(tokens) // self.natoms
        table = np.asarray(tokens).reshape(self.natoms, ncol)
        species = table[:, 0].astype("U3")
        positions = table[:, 1:4].astype(np.float64)

        lat = _LATTICE_RE.search(comment)
        cell = np.asarray(lat.group(1).split(), dtype=np.float64).reshape(3, 3)
        box = np.diag(cell).copy()

        tmatch = _TIME_RE.search(comment)
        time_fs = float(tmatch.group(1)) if tmatch else np.nan
        return species, positions, box, time_fs


def shannon(counts) -> float:
    counts = np.asarray(counts, dtype=np.float64)
    counts = counts[counts > 0]
    total = counts.sum()
    if total <= 0:
        return 0.0
    p = counts / total
    return float(-(p * np.log(p)).sum())


def frame_entropies(species_id: np.ndarray, is_oxygen: np.ndarray,
                    positions: np.ndarray, box: np.ndarray):
    """Compute (S_atom, S_config, n_motifs, n_clusters, n_metal, max_cluster)."""
    wrapped = np.mod(positions, box)
    tree = cKDTree(wrapped, boxsize=box)
    pairs = tree.query_pairs(MET_CUT, output_type="ndarray")

    natoms = len(species_id)
    cn_o = np.zeros(natoms, dtype=np.int32)
    cn_m = np.zeros(natoms, dtype=np.int32)

    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        delta = wrapped[i] - wrapped[j]
        delta -= box * np.round(delta / box)
        dist = np.linalg.norm(delta, axis=1)

        oi, oj = is_oxygen[i], is_oxygen[j]
        close = dist < O_CUT

        # oxygen coordination (within O_CUT)
        cn_o += np.bincount(i[oj & close], minlength=natoms).astype(np.int32)
        cn_o += np.bincount(j[oi & close], minlength=natoms).astype(np.int32)
        # metal coordination (within MET_CUT)
        cn_m += np.bincount(i[~oj], minlength=natoms).astype(np.int32)
        cn_m += np.bincount(j[~oi], minlength=natoms).astype(np.int32)

    motif = (species_id.astype(np.int64) * 10000
             + np.minimum(cn_o, 99) * 100
             + np.minimum(cn_m, 99))
    _, motif_counts = np.unique(motif, return_counts=True)
    s_atom = shannon(motif_counts)

    metal_idx = np.flatnonzero(~is_oxygen)
    n_metal = len(metal_idx)
    remap = np.full(natoms, -1, dtype=np.int64)
    remap[metal_idx] = np.arange(n_metal)

    if len(pairs):
        mm = pairs[(~is_oxygen[pairs[:, 0]]) & (~is_oxygen[pairs[:, 1]])]
        a, b = remap[mm[:, 0]], remap[mm[:, 1]]
    else:
        a = b = np.empty(0, dtype=np.int64)

    graph = coo_matrix(
        (np.ones(len(a), dtype=np.int8), (a, b)), shape=(n_metal, n_metal)
    ).tocsr()
    n_clusters, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels, minlength=n_clusters)
    s_config = shannon(sizes)

    return (s_atom, s_config, int(len(motif_counts)), int(n_clusters),
            int(n_metal), int(sizes.max()) if n_clusters else 0)


def sample_indices(n_frames: int, n_sample: int) -> np.ndarray:
    if n_frames <= n_sample:
        return np.arange(n_frames, dtype=int)
    return np.unique(np.round(np.linspace(0, n_frames - 1, n_sample)).astype(int))


def analyze(path: str, frame_indices, progress=None):
    """Analyze the requested frames of one trajectory."""
    with Trajectory(path) as traj:
        species0, _, _, _ = traj.frame(0)
        names = sorted(set(species0.tolist()))
        sid_of = {name: k for k, name in enumerate(names)}

        rows = []
        for count, fidx in enumerate(frame_indices, start=1):
            species, positions, box, time_fs = traj.frame(int(fidx))
            species_id = np.asarray([sid_of.get(s, -1) for s in species])
            is_oxygen = species == OXYGEN
            vals = frame_entropies(species_id, is_oxygen, positions, box)
            rows.append((int(fidx), time_fs / 1000.0) + vals)
            if progress is not None:
                progress(count, len(frame_indices), rows[-1])

    dtype = [("frame", "i8"), ("time_ps", "f8"), ("s_atom", "f8"),
             ("s_config", "f8"), ("n_motifs", "i8"), ("n_clusters", "i8"),
             ("n_metal", "i8"), ("max_cluster", "i8")]
    return np.array(rows, dtype=dtype)
