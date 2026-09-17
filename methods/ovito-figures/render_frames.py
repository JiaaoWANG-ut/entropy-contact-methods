#!/usr/bin/env python3
"""High-res shaded-sphere frames (OVITO-like Phong balls, core-framed camera)."""
from __future__ import annotations

import glob
import math
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from PIL import Image, ImageDraw, ImageFont

OUT = "/public/home/jwang/test/testgpumd/ovito-figures"
FRAME_ROOT = os.path.join(OUT, "frames")
PNG_ROOT = os.path.join(OUT, "png")

# Jmol colors, slightly lifted for print
COLORS = {
    "Ni": np.array([0.560, 0.560, 0.860], dtype=np.float32),
    "O":  np.array([0.980, 0.180, 0.160], dtype=np.float32),
    "Fe": np.array([0.920, 0.480, 0.180], dtype=np.float32),
    "Zn": np.array([0.460, 0.540, 0.780], dtype=np.float32),
    "Cr": np.array([0.420, 0.760, 0.520], dtype=np.float32),
    "Ru": np.array([0.160, 0.620, 0.640], dtype=np.float32),
}
RADII = {
    "O": 0.70,
    "Ni": 1.02,
    "Fe": 1.06,
    "Zn": 1.04,
    "Cr": 1.02,
    "Ru": 1.10,
}
DEFAULT_RGB = np.array([0.70, 0.70, 0.70], dtype=np.float32)
LIGHT = np.array([0.32, 0.48, 0.82], dtype=np.float32)
LIGHT /= np.linalg.norm(LIGHT)
FRAME_PX = 1800


def _font(size, bold=True):
    names = [
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for p in names:
        if os.path.isfile(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def read_xyz(path):
    with open(path) as f:
        first = f.readline().strip()
        if not first:
            raise IOError("incomplete xyz (empty): %s" % path)
        n = int(first)
        comment = f.readline()
        el = []
        xyz = np.empty((n, 3), dtype=np.float64)
        for i in range(n):
            line = f.readline()
            if not line:
                raise IOError("incomplete xyz at atom %d/%d: %s" % (i, n, path))
            sp = line.split()
            el.append(sp[0])
            xyz[i] = (float(sp[1]), float(sp[2]), float(sp[3]))
    return el, xyz, comment


def isometric_project(xyz):
    com = xyz.mean(axis=0)
    p = xyz - com
    a = math.radians(32.0)
    b = math.radians(38.0)
    Rx = np.array(
        [[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]],
        dtype=np.float64,
    )
    Rz = np.array(
        [[math.cos(b), -math.sin(b), 0], [math.sin(b), math.cos(b), 0], [0, 0, 1]],
        dtype=np.float64,
    )
    return p @ Rz.T @ Rx.T


def core_mask(xyz, cutoff=3.15, min_nb=7):
    """Atoms in the condensed nanoparticle(s), excluding vapor."""
    n = len(xyz)
    if n == 0:
        return np.zeros(0, dtype=bool)
    inv = 1.0 / cutoff
    keys = np.floor(xyz * inv).astype(np.int32)
    cells = defaultdict(list)
    for i, k in enumerate(map(tuple, keys)):
        cells[k].append(i)
    counts = np.zeros(n, dtype=np.int16)
    c2 = cutoff * cutoff
    offs = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for i, (x, y, z) in enumerate(xyz):
        kx, ky, kz = keys[i]
        acc = 0
        for dx, dy, dz in offs:
            for j in cells.get((kx + dx, ky + dy, kz + dz), ()):
                if j <= i:
                    continue
                ddx = xyz[j, 0] - x
                ddy = xyz[j, 1] - y
                ddz = xyz[j, 2] - z
                if ddx * ddx + ddy * ddy + ddz * ddz <= c2:
                    acc += 1
                    counts[j] += 1
        counts[i] += acc
    return counts >= min_nb


def view_window(proj, xyz):
    """Frame on condensed cores; keep full box if the system is still multi-particle."""
    mask = core_mask(xyz)
    frac = float(mask.mean()) if len(mask) else 1.0
    use = proj[mask] if 0.25 < frac < 0.92 else proj
    if len(use) < 32:
        use = proj
    lo = use.min(axis=0)
    hi = use.max(axis=0)
    pad = 7.5
    xmin, xmax = float(lo[0] - pad), float(hi[0] + pad)
    ymin, ymax = float(lo[1] - pad), float(hi[1] + pad)
    span = max(xmax - xmin, ymax - ymin, 20.0)
    cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    return cx, cy, span


def sphere_lut(r_px, rgb):
    """Anti-aliased Phong sphere stamp and its front-z offset."""
    r = max(float(r_px), 1.25)
    rad = int(math.ceil(r + 1.4))
    yy, xx = np.mgrid[-rad : rad + 1, -rad : rad + 1]
    # pixel centers
    fx = xx.astype(np.float32) + 0.0
    fy = yy.astype(np.float32) + 0.0
    rr2 = fx * fx + fy * fy
    r2 = r * r
    cover = np.clip(r + 0.65 - np.sqrt(np.maximum(rr2, 1e-8)), 0.0, 1.0).astype(np.float32)
    inside = rr2 <= (r + 0.65) ** 2
    nz = np.zeros_like(fx)
    nz[inside] = np.sqrt(np.maximum(r2 - rr2[inside], 0.0)) / r
    nx = np.zeros_like(fx)
    ny = np.zeros_like(fx)
    nx[inside] = fx[inside] / r
    ny[inside] = -fy[inside] / r  # image y down
    ndotl = np.clip(nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2], -1.0, 1.0)
    wrap = 0.5 * ndotl + 0.5  # hemisphere wrap, OVITO-like fill
    hx, hy, hz = LIGHT[0], LIGHT[1], LIGHT[2] + 1.0
    hn = math.sqrt(hx * hx + hy * hy + hz * hz)
    hx, hy, hz = hx / hn, hy / hn, hz / hn
    ndoth = np.clip(nx * hx + ny * hy + nz * hz, 0.0, 1.0)
    spec = ndoth ** 48 * 0.20
    shade = 0.40 + 0.52 * wrap
    rgb = rgb.astype(np.float32)
    stamp = np.zeros(xx.shape + (3,), dtype=np.float32)
    for c in range(3):
        stamp[..., c] = np.clip(rgb[c] * shade + spec, 0.0, 1.0)
    rim = np.clip(1.0 - nz, 0.0, 1.0) ** 1.6
    stamp *= (1.0 - 0.12 * rim)[..., None]
    zoff = (nz * r).astype(np.float32)
    return stamp, cover, zoff, rad


def raster_spheres(xs, ys, zs, el, cx, cy, span, size=FRAME_PX):
    n = len(xs)
    scale = size / span
    px = (xs - (cx - 0.5 * span)) * scale
    py = ((cy + 0.5 * span) - ys) * scale  # flip y for image
    img = np.ones((size, size, 3), dtype=np.float32)
    zbuf = np.full((size, size), -1e9, dtype=np.float32)
    cache = {}
    order = np.argsort(zs)  # far to near, z-buffer still used
    for i in order:
        e = el[i]
        rgb = COLORS.get(e, DEFAULT_RGB)
        rA = RADII.get(e, 0.85)
        r_px = rA * scale
        key = (e, round(r_px, 2))
        if key not in cache:
            cache[key] = sphere_lut(r_px, rgb)
        stamp, cover, zoff, rad = cache[key]
        xi = int(round(px[i]))
        yi = int(round(py[i]))
        x0, x1 = xi - rad, xi + rad + 1
        y0, y1 = yi - rad, yi + rad + 1
        if x1 <= 0 or y1 <= 0 or x0 >= size or y0 >= size:
            continue
        sx0 = 0 if x0 >= 0 else -x0
        sy0 = 0 if y0 >= 0 else -y0
        sx1 = stamp.shape[1] if x1 <= size else stamp.shape[1] - (x1 - size)
        sy1 = stamp.shape[0] if y1 <= size else stamp.shape[0] - (y1 - size)
        dx0 = max(x0, 0)
        dy0 = max(y0, 0)
        dx1 = dx0 + (sx1 - sx0)
        dy1 = dy0 + (sy1 - sy0)
        sl = (slice(sy0, sy1), slice(sx0, sx1))
        dl = (slice(dy0, dy1), slice(dx0, dx1))
        znew = zs[i] + zoff[sl]
        vis = (cover[sl] > 0.04) & (znew >= zbuf[dl])
        if not np.any(vis):
            continue
        a = cover[sl][vis][:, None]
        dst = img[dl]
        dst[vis] = stamp[sl][vis] * a + dst[vis] * (1.0 - a)
        zbuf[dl][vis] = znew[vis]
        img[dl] = dst
    img = apply_soft_ao(img, zbuf)
    return img


def apply_soft_ao(img, zbuf):
    """Screen-space ambient occlusion: soft contact shadows like OVITO default AO."""
    step = 4
    z = zbuf[::step, ::step]
    valid = z > -1e8
    zf = np.where(valid, z, np.nan)
    occ = np.zeros(z.shape, dtype=np.float32)
    shifts = [
        (0, 1), (0, -1), (1, 0), (-1, 0),
        (1, 1), (1, -1), (-1, 1), (-1, -1),
        (0, 2), (2, 0), (0, -2), (-2, 0),
        (2, 2), (2, -2), (-2, 2), (-2, -2),
        (0, 3), (3, 0), (-3, 0), (0, -3),
    ]
    for dy, dx in shifts:
        zs = np.roll(np.roll(zf, dy, 0), dx, 1)
        occ += np.clip(np.nan_to_num(zs - zf, nan=0.0) / 3.8, 0.0, 1.0)
    occ /= float(len(shifts))
    ao = np.clip(1.0 - 0.55 * occ, 0.52, 1.0)
    ao = np.where(valid, ao, 1.0).astype(np.float32)
    ao_u = np.repeat(np.repeat(ao, step, 0), step, 1)[: img.shape[0], : img.shape[1]]
    pad = np.pad(ao_u, 1, mode="edge")
    blur = (
        pad[0:-2, 0:-2] + 2 * pad[0:-2, 1:-1] + pad[0:-2, 2:]
        + 2 * pad[1:-1, 0:-2] + 4 * pad[1:-1, 1:-1] + 2 * pad[1:-1, 2:]
        + pad[2:, 0:-2] + 2 * pad[2:, 1:-1] + pad[2:, 2:]
    ) / 16.0
    out = img * blur[:, :, None]
    out[zbuf <= -1e8] = 1.0
    return out


def overlay_labels(img, title, subtitle):
    im = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8), "RGB")
    d = ImageDraw.Draw(im, "RGBA")
    ft = _font(42, True)
    fs = _font(46, True)
    # translucent chips
    def chip(xy, text, font, anchor):
        bbox = d.textbbox(xy, text, font=font, anchor=anchor)
        pad = 8
        box = [bbox[0] - pad, bbox[1] - pad, bbox[2] + pad, bbox[3] + pad]
        d.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 210))
        d.text(xy, text, font=font, fill=(18, 18, 20, 255), anchor=anchor)

    chip((28, 28), title, ft, "lt")
    chip((im.width - 28, 28), subtitle, fs, "rt")
    return im


def render_frame(xyz_path, png_path, title, subtitle, size=FRAME_PX):
    el, xyz, _ = read_xyz(xyz_path)
    proj = isometric_project(xyz)
    cx, cy, span = view_window(proj, xyz)
    img = raster_spheres(proj[:, 0], proj[:, 1], proj[:, 2], el, cx, cy, span, size=size)
    im = overlay_labels(img, title, subtitle)
    os.makedirs(os.path.dirname(png_path), exist_ok=True)
    tmp = png_path + ".tmp.png"
    im.save(tmp, "PNG", compress_level=3)
    os.replace(tmp, png_path)
    print("rendered", png_path, flush=True)


def elem_title(comp):
    if comp.startswith("NiO-6"):
        return "NiO-6  Fe/Zn/Cr/Ru"
    if comp.startswith("NiO-8"):
        return "NiO-8  Fe/Zn/Cr/Ru"
    return comp.replace("-", " ")


def _jobs(force=False):
    times = ["0.0", "0.25", "0.5", "0.75", "1.0", "2.0", "3.0", "4.0", "5.0"]
    labels = {
        "0.0": "0 ns",
        "0.25": "0.25 ns",
        "0.5": "0.5 ns",
        "0.75": "0.75 ns",
        "1.0": "1 ns",
        "2.0": "2 ns",
        "3.0": "3 ns",
        "4.0": "4 ns",
        "5.0": "5 ns",
    }
    jobs = []
    missing = 0
    for comp_dir in sorted(glob.glob(os.path.join(FRAME_ROOT, "*"))):
        if not os.path.isdir(comp_dir):
            continue
        comp = os.path.basename(comp_dir)
        title = elem_title(comp)
        for tdir in sorted(glob.glob(os.path.join(comp_dir, "*K"))):
            temp = os.path.basename(tdir)
            for ts in times:
                xyz = os.path.join(tdir, ts + "ns.xyz")
                png = os.path.join(PNG_ROOT, comp, temp, ts + "ns.png")
                if not os.path.isfile(xyz):
                    missing += 1
                    print("skip missing", xyz, flush=True)
                    continue
                if (
                    not force
                    and os.path.isfile(png)
                    and os.path.getmtime(png) >= os.path.getmtime(xyz)
                    and os.path.getsize(png) > 80000
                ):
                    continue
                jobs.append((xyz, png, title, labels[ts]))
    return jobs, missing


def _worker(job):
    xyz, png, title, subtitle = job
    try:
        render_frame(xyz, png, title, subtitle)
        return png, None
    except Exception as exc:
        return png, str(exc)


def render_all(force=False, workers=10):
    jobs, missing = _jobs(force=force)
    print("jobs=%d missing=%d workers=%d" % (len(jobs), missing, workers), flush=True)
    if not jobs:
        return
    if workers <= 1:
        for job in jobs:
            _worker(job)
    else:
        with Pool(processes=workers) as pool:
            for png, err in pool.imap_unordered(_worker, jobs):
                if err:
                    print("FAIL", png, err, flush=True)
    print("done render", flush=True)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        render_frame(sys.argv[1], sys.argv[2], title="test", subtitle="0 ns")
    elif "--force" in sys.argv:
        render_all(force=True, workers=10)
    else:
        render_all(force=False, workers=10)
