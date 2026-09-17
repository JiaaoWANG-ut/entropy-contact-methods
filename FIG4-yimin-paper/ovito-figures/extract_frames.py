#!/usr/bin/env python3
"""Stream GPUMD dump.xyz and write 0/1/2/3/4/5 ns frames (xyz, species+pos)."""
from __future__ import print_function
import os
import sys

ROOT = "/public/home/jwang/test/testgpumd/production-run/FIG4-run"
OUT = "/public/home/jwang/test/testgpumd/ovito-figures/frames"
COMPS = ["NiO-6-Fe-Zn-Cr-Ru", "NiO-8-Fe-Zn-Cr-Ru"]
TEMPS = ["1500K", "1800K", "2000K", "2500K", "3000K"]
# dump Time= is in fs; dump every 400 fs
# 0-1 ns: 5 frames, then 2/3/4/5 ns
TIMES_NS = [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0]
NS_FS = [int(round(t * 1e6)) for t in TIMES_NS]
TOL = 1.0  # fs


def frame_tag(ns):
    if abs(ns - round(ns)) < 1e-9:
        return "%.1f" % ns
    return ("%g" % ns)


def parse_time_fs(comment):
    text = comment.decode("ascii", "replace") if isinstance(comment, bytes) else comment
    for part in text.split():
        if part.startswith("Time="):
            return float(part.split("=", 1)[1])
    return None


def skip_lines(f, n):
    """Skip n newline-terminated records from a binary file."""
    left = n
    while left > 0:
        buf = f.read(1 << 20)
        if not buf:
            return False
        c = buf.count(b"\n")
        if c < left:
            left -= c
            continue
        idx = 0
        for _ in range(left):
            idx = buf.find(b"\n", idx) + 1
        f.seek(f.tell() - (len(buf) - idx))
        return True
    return True


def rewrite_properties(comment):
    """Keep Lattice/Time/pbc, but match the 4 written columns (species + pos)."""
    parts = comment.split()
    out = []
    for part in parts:
        if part.startswith("Properties="):
            out.append("Properties=species:S:1:pos:R:3")
        else:
            out.append(part)
    return " ".join(out)


def write_xyz(path, n, comment, atom_lines):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as g:
        g.write("%d\n" % n)
        comment = rewrite_properties(comment.rstrip("\n"))
        g.write(comment + "\n")
        for line in atom_lines:
            sp = line.split()
            g.write("%s %s %s %s\n" % (sp[0], sp[1], sp[2], sp[3]))
    os.replace(tmp, path)


def copy_zero(comp, temp):
    src = os.path.join(ROOT, "models", comp + ".xyz")
    dst = os.path.join(OUT, comp, temp, "0.0ns.xyz")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(src) as f:
        data = f.read()
    with open(dst, "w") as g:
        g.write(data)
    print("  wrote", dst, "from model.xyz", flush=True)
    return dst


def extract_one(comp, temp):
    dump = os.path.join(ROOT, comp, temp, "dump.xyz")
    outdir = os.path.join(OUT, comp, temp)
    os.makedirs(outdir, exist_ok=True)
    copy_zero(comp, temp)
    needed = {}
    for tfs, ns in zip(NS_FS, TIMES_NS):
        if ns <= 0:
            continue
        needed[float(tfs)] = ns
    found = {}
    for key, ns in list(needed.items()):
        dst = os.path.join(outdir, frame_tag(ns) + "ns.xyz")
        if os.path.isfile(dst) and os.path.getsize(dst) > 1000:
            found[ns] = dst
            needed.pop(key)
            print("  already", dst, flush=True)
    if not needed:
        return found
    if not os.path.isfile(dump):
        print("  MISSING dump", dump, flush=True)
        return found
    print("  scanning", dump, flush=True)
    last_needed = max(needed)
    with open(dump, "rb", buffering=8 * 1024 * 1024) as f:
        while needed:
            nline = f.readline()
            if not nline:
                break
            nline = nline.strip()
            if not nline:
                continue
            n = int(nline)
            comment = f.readline()
            if not comment:
                break
            tfs = parse_time_fs(comment)
            hit = None
            if tfs is not None:
                for key in list(needed):
                    if abs(tfs - key) <= TOL:
                        hit = key
                        break
            if hit is not None:
                atoms = [f.readline() for _ in range(n)]
                ns = needed.pop(hit)
                dst = os.path.join(outdir, frame_tag(ns) + "ns.xyz")
                write_xyz(dst, n, comment.decode("ascii", "replace"),
                          [a.decode("ascii", "replace") for a in atoms])
                found[ns] = dst
                print("  wrote", dst, "Time=%.0f fs" % tfs, flush=True)
                if not needed:
                    break
            else:
                if not skip_lines(f, n):
                    break
                if tfs is not None and tfs > last_needed + TOL:
                    break
    missing = sorted(needed.values())
    if missing:
        print("  missing ns:", missing, flush=True)
    return found


def main():
    only = sys.argv[1:]  # optional: Comp T
    if len(only) == 2:
        extract_one(only[0], only[1])
        return
    for comp in COMPS:
        for temp in TEMPS:
            print("=== %s %s ===" % (comp, temp), flush=True)
            extract_one(comp, temp)


if __name__ == "__main__":
    main()
