#!/usr/bin/env python3
"""OVITO Tachyon + AO frames, same style as the user snippet."""
import os
import glob
import sys

os.environ["OVITO_MODE"] = "script"
os.environ["QT_QPA_PLATFORM"] = "offscreen"

import numpy as np
import ovito
from ovito.data import ParticleType
from ovito.io import import_file
from ovito.modifiers import PythonModifier
from ovito.vis import TachyonRenderer, Viewport

OUT = "/public/home/jwang/test/testgpumd/ovito-figures"
FRAME_ROOT = os.path.join(OUT, "frames")
PNG_ROOT = os.path.join(OUT, "png")

TIMES = ["0.0", "0.25", "0.5", "0.75", "1.0", "2.0", "3.0", "4.0", "5.0"]
LABELS = {
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
SIZE = (3600, 3600)


def setup_types(frame, data):
    p = data.particles_
    names = {t.id: t.name for t in p.particle_types.types}
    types = np.asarray(p.particle_types)
    for prop in ("Color", "Radius"):
        if prop not in p:
            p.create_property(prop)
    colors, radii = np.asarray(p["Color"]).copy(), np.asarray(p["Radius"]).copy()
    for tid, name in names.items():
        pt = ParticleType()
        pt.name = name
        pt.load_defaults()
        m = types == tid
        colors[m] = pt.color
        radii[m] = pt.radius
    p["Color"][:] = colors
    p["Radius"][:] = radii
    p.vis.scaling = 1.0

    pos = np.asarray(p.positions)
    # keep all 5 clusters (NiO + Fe/Zn/Cr/Ru); drop only extreme outliers
    lo = np.percentile(pos, 0.2, axis=0) - 6.0
    hi = np.percentile(pos, 99.8, axis=0) + 6.0
    if data.cell is not None:
        data.cell_[:, 0] = (hi[0] - lo[0], 0, 0)
        data.cell_[:, 1] = (0, hi[1] - lo[1], 0)
        data.cell_[:, 2] = (0, 0, hi[2] - lo[2])
        data.cell_[:, 3] = lo
        data.cell_.vis.enabled = True
        data.cell_.vis.line_width = 0.45
        data.cell_.vis.rendering_color = (0.22, 0.22, 0.22)


def elem_title(comp):
    if comp.startswith("NiO-6"):
        return "NiO-6  Fe/Zn/Cr/Ru"
    if comp.startswith("NiO-8"):
        return "NiO-8  Fe/Zn/Cr/Ru"
    return comp


def _core_look(data):
    """Particle-core center and square ortho height (Å), ignoring the 300 Å box."""
    cell = np.asarray(data.cell)
    lo = cell[:, 3]
    hi = lo + cell[:, 0] + cell[:, 1] + cell[:, 2]
    look_at = 0.5 * (lo + hi)
    span = hi - lo
    # top-down xy: fit the four metal satellites around NiO
    fov = 1.08 * float(max(span[0], span[1]))
    return look_at, fov


def render_side(xyz, out, size=SIZE):
    """正交俯视：沿 -z 看，+y 朝上，无透视，5 个球同时可见。"""
    for pl in list(ovito.scene.pipelines):
        pl.remove_from_scene()
    pipe = import_file(xyz)
    pipe.modifiers.append(PythonModifier(function=setup_types))
    pipe.add_to_scene()
    look_at, fov = _core_look(pipe.compute())
    vp = Viewport(type=Viewport.Type.Ortho)
    vp.camera_dir = (0.0, 0.0, -1.0)
    vp.camera_up = (0.0, 1.0, 0.0)
    vp.fov = 0.5 * fov
    vp.camera_pos = tuple(np.asarray(look_at) - 80.0 * np.array(vp.camera_dir))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + ".tmp.png"
    vp.render_image(
        filename=tmp,
        size=size,
        frame=0,
        background=(1, 1, 1),
        renderer=TachyonRenderer(
            ambient_occlusion=True,
            shadows=False,
            antialiasing=True,
            antialiasing_samples=12,
        ),
    )
    os.replace(tmp, out)
    pipe.remove_from_scene()
    print("rendered", out, flush=True)


def label_png(path, title, subtitle):
    from PIL import Image, ImageDraw, ImageFont
    im = Image.open(path).convert("RGBA")
    d = ImageDraw.Draw(im, "RGBA")
    font_paths = [
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ]
    font = ImageFont.load_default()
    for p in font_paths:
        if os.path.isfile(p):
            font = ImageFont.truetype(p, 84)
            break

    def chip(xy, text, anchor):
        bbox = d.textbbox(xy, text, font=font, anchor=anchor)
        pad = 16
        box = [bbox[0] - pad, bbox[1] - pad, bbox[2] + pad, bbox[3] + pad]
        d.rounded_rectangle(box, radius=20, fill=(255, 255, 255, 215))
        d.text(xy, text, font=font, fill=(18, 18, 20, 255), anchor=anchor)

    chip((56, 56), title, "lt")
    chip((im.width - 56, 56), subtitle, "rt")
    im.convert("RGB").save(path, compress_level=1)


def all_jobs():
    jobs = []
    for comp_dir in sorted(glob.glob(os.path.join(FRAME_ROOT, "*"))):
        if not os.path.isdir(comp_dir):
            continue
        comp = os.path.basename(comp_dir)
        title = elem_title(comp)
        for tdir in sorted(glob.glob(os.path.join(comp_dir, "*K"))):
            temp = os.path.basename(tdir)
            for ts in TIMES:
                xyz = os.path.join(tdir, ts + "ns.xyz")
                png = os.path.join(PNG_ROOT, comp, temp, ts + "ns.png")
                if os.path.isfile(xyz):
                    jobs.append((xyz, png, title, LABELS[ts]))
                else:
                    print("skip missing", xyz, flush=True)
    return jobs


if __name__ == "__main__":
    jobs = all_jobs()
    if len(sys.argv) == 3 and not sys.argv[1].startswith("--"):
        jobs = [(sys.argv[1], sys.argv[2], "test", "0 ns")]
    elif len(sys.argv) == 3 and sys.argv[1] == "--shard":
        i, n = map(int, sys.argv[2].split("/"))
        jobs = [j for k, j in enumerate(jobs) if k % n == i]
    print("ovito jobs", len(jobs), flush=True)
    for xyz, png, title, subtitle in jobs:
        render_side(xyz, png)
        label_png(png, title, subtitle)
    print("done ovito render", flush=True)
