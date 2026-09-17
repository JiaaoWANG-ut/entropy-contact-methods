#!/usr/bin/env python3
"""Stitch per-composition T x time grids. Rows = temperature, cols = time."""
from __future__ import annotations

import os
import glob

import numpy as np
from PIL import Image, ImageDraw, ImageFont

OUT = "/public/home/jwang/test/testgpumd/ovito-figures"
PNG_ROOT = os.path.join(OUT, "png")
SUM_ROOT = os.path.join(OUT, "summary")

COMPS = ["NiO-6-Fe-Zn-Cr-Ru", "NiO-8-Fe-Zn-Cr-Ru"]
TEMPS = ["1500K", "1800K", "2000K", "2500K", "3000K"]
TIMES = ["0.0", "0.25", "0.5", "0.75", "1.0", "2.0", "3.0", "4.0", "5.0"]
TIME_LABELS = ["0 ns", "0.25 ns", "0.5 ns", "0.75 ns", "1 ns", "2 ns", "3 ns", "4 ns", "5 ns"]

# OVITO default element colors (ParticleType.load_defaults)
COLORS = {
    "Ni": (80, 208, 80),
    "O": (255, 13, 13),
    "Fe": (224, 102, 51),
    "Zn": (125, 128, 176),
    "Cr": (138, 153, 199),
    "Ru": (36, 143, 143),
}


def font(size, bold=False):
    candidates = [
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/liberation/LiberationSans-Bold.ttf" if bold else "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",
    ]
    for p in candidates:
        if os.path.isfile(p):
            return ImageFont.truetype(p, size)
    return ImageFont.load_default()


def elem_banner(comp):
    if comp.startswith("NiO-6"):
        return "NiO-6 + Fe / Zn / Cr / Ru   (NEP89 NVT)"
    if comp.startswith("NiO-8"):
        return "NiO-8 + Fe / Zn / Cr / Ru   (NEP89 NVT)"
    return comp


def load_cell(path, box):
    if path and os.path.isfile(path):
        im = Image.open(path).convert("RGB")
        im.thumbnail((box, box), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (box, box), (255, 255, 255))
        x = (box - im.width) // 2
        y = (box - im.height) // 2
        canvas.paste(im, (x, y))
        return canvas
    canvas = Image.new("RGB", (box, box), (245, 245, 245))
    d = ImageDraw.Draw(canvas)
    d.text((box // 2, box // 2), "not yet\n(5 ns pending)" if "5.0" in (path or "") else "missing",
           fill=(120, 120, 120), anchor="mm", font=font(36), align="center")
    return canvas


def stitch_one(comp, cell=1600, pad=24, left=260, top=200, bottom=140):
    os.makedirs(SUM_ROOT, exist_ok=True)
    ncol, nrow = len(TIMES), len(TEMPS)
    W = left + ncol * cell + (ncol + 1) * pad
    H = top + nrow * cell + (nrow + 1) * pad + bottom
    canvas = Image.new("RGB", (W, H), (255, 255, 255))
    d = ImageDraw.Draw(canvas)
    title = elem_banner(comp)
    d.text((W // 2, 36), title, fill=(20, 20, 20), anchor="mt", font=font(72, bold=True))

    for j, lab in enumerate(TIME_LABELS):
        x = left + pad + j * (cell + pad) + cell // 2
        d.text((x, top - 16), lab, fill=(20, 20, 20), anchor="ms", font=font(40, bold=True))

    for i, temp in enumerate(TEMPS):
        y = top + pad + i * (cell + pad) + cell // 2
        d.text((left - 28, y), temp.replace("K", " K"), fill=(20, 20, 20),
               anchor="rm", font=font(52, bold=True))
        for j, ts in enumerate(TIMES):
            png = os.path.join(PNG_ROOT, comp, temp, ts + "ns.png")
            cell_im = load_cell(png if os.path.isfile(png) else png, cell)
            x0 = left + pad + j * (cell + pad)
            y0 = top + pad + i * (cell + pad)
            canvas.paste(cell_im, (x0, y0))
            d.rectangle([x0, y0, x0 + cell - 1, y0 + cell - 1], outline=(210, 210, 210), width=1)

    # color legend
    lx = pad + 48
    ly = H - 88
    d.text((lx, ly - 6), "elements:", fill=(40, 40, 40), font=font(40), anchor="ls")
    lx += 220
    for name, rgb in COLORS.items():
        d.ellipse([lx, ly - 36, lx + 42, ly + 6], fill=rgb, outline=(40, 40, 40), width=3)
        d.text((lx + 54, ly - 12), name, fill=(20, 20, 20), font=font(40), anchor="ls")
        lx += 180

    out = os.path.join(SUM_ROOT, "summary_%s.png" % comp)
    canvas.save(out, compress_level=1)
    print("wrote", out, canvas.size, flush=True)
    return out


def main():
    outs = []
    for comp in COMPS:
        # only stitch if at least one png exists
        hits = glob.glob(os.path.join(PNG_ROOT, comp, "*", "*.png"))
        if not hits:
            print("no pngs for", comp, flush=True)
            continue
        outs.append(stitch_one(comp))
    return outs


if __name__ == "__main__":
    main()
