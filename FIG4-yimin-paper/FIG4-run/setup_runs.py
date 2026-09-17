#!/usr/bin/env python3
"""Convert FIG4 CIFs to GPUMD model.xyz and write per-T run directories."""
from __future__ import print_function
import os
from collections import Counter

ROOT = os.path.dirname(os.path.abspath(__file__))
POTENTIAL = (
    "/public/home/jwang/test/testgpumd/GPUMD/potentials/"
    "nep/nep89_20250409/nep89_20250409.txt"
)
TEMPS = [1500, 1800, 2000, 2500, 3000]
TIME_STEP_FS = 1.0
RUN_PS = 10000.0
DUMP_EVERY = 400
RESTART_EVERY = 100000
CIFS = [
    "NiO-6-Fe-Zn-Cr-Ru.cif",
    "NiO-8-Fe-Zn-Cr-Ru.cif",
]


def parse_cif(path):
    a = b = c = None
    atoms = []
    in_loop = False
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            key = s[0]
            if key == "_cell_length_a":
                a = float(s[1])
            elif key == "_cell_length_b":
                b = float(s[1])
            elif key == "_cell_length_c":
                c = float(s[1])
            elif key.startswith("_atom_site_"):
                in_loop = True
                continue
            if in_loop and len(s) >= 5 and s[1][0].isalpha() and not s[1].startswith("_"):
                elem, fx, fy, fz = s[1], float(s[2]), float(s[3]), float(s[4])
                atoms.append((elem, fx * a, fy * b, fz * c))
    if a is None or not atoms:
        raise RuntimeError("failed to parse %s" % path)
    return a, b, c, atoms


def write_xyz(path, a, b, c, atoms):
    with open(path, "w") as f:
        f.write("%d\n" % len(atoms))
        f.write(
            'Lattice="%.8f 0 0 0 %.8f 0 0 0 %.8f" '
            'Properties=species:S:1:pos:R:3 pbc="T T T"\n' % (a, b, c)
        )
        for elem, x, y, z in atoms:
            f.write("%s %.8f %.8f %.8f\n" % (elem, x, y, z))


def write_run_in(path, tag, temp):
    nsteps = int(round(RUN_PS * 1000.0 / TIME_STEP_FS))
    with open(path, "w") as f:
        f.write("# %s  NEP89  NVT-Langevin  %d K  %g ps\n" % (tag, temp, RUN_PS))
        f.write("# dt = %g fs; dump.xyz (exyz) every %d steps (%.3f ps)\n" % (
            TIME_STEP_FS, DUMP_EVERY, DUMP_EVERY * TIME_STEP_FS / 1000.0))
        f.write("potential  %s\n" % POTENTIAL)
        f.write("velocity   %d\n" % temp)
        f.write("time_step  %.1f\n" % TIME_STEP_FS)
        f.write("ensemble   nvt_lan %d %d 100\n" % (temp, temp))
        f.write("dump_thermo %d\n" % DUMP_EVERY)
        f.write("dump_xyz    %d dump.xyz velocity force potential\n" % DUMP_EVERY)
        f.write("dump_restart %d\n" % RESTART_EVERY)
        f.write("run         %d\n" % nsteps)


def write_submit(path, jobname, rundir):
    with open(path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=%s\n" % jobname)
        f.write("#SBATCH --partition=bw1000\n")
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --ntasks=1\n")
        f.write("#SBATCH --cpus-per-task=8\n")
        f.write("#SBATCH --gres=dcu:1\n")
        f.write("#SBATCH --mem=48G\n")
        f.write("#SBATCH --time=7-00:00:00\n")
        f.write("#SBATCH --output=slurm_%j.out\n")
        f.write("#SBATCH --error=slurm_%j.err\n")
        f.write("\n")
        f.write("set -e\n")
        f.write("source /public/software/dtk-25.04.1/env.sh\n")
        f.write("export LD_LIBRARY_PATH=/public/software/compiler/gnu/7.2.0/lib64:${LD_LIBRARY_PATH}\n")
        f.write("GPUMD=/public/home/jwang/test/testgpumd/GPUMD/src/gpumd\n")
        f.write("cd %s\n" % rundir)
        f.write("\n")
        f.write('echo "===== JOB INFO ====="\n')
        f.write('echo "host=$(hostname) job=${SLURM_JOB_ID} date=$(date)"\n')
        f.write('echo "node=${SLURM_JOB_NODELIST} partition=${SLURM_JOB_PARTITION}"\n')
        f.write('echo "pwd=$(pwd)"\n')
        f.write("hy-smi 2>/dev/null | head -12 || true\n")
        f.write("echo \"===== RUN =====\"\n")
        f.write("START=$(date +%s)\n")
        f.write('"${GPUMD}"\n')
        f.write("END=$(date +%s)\n")
        f.write('echo "===== DONE wall_seconds=$((END-START)) $(date) ====="\n')
    os.chmod(path, 0o755)


def short_tag(cif_name):
    # NiO-6-Fe-Zn-Cr-Ru.cif -> NiO-6
    return cif_name.replace(".cif", "").split("-Fe")[0]


def main():
    models_dir = os.path.join(ROOT, "models")
    os.makedirs(models_dir, exist_ok=True)
    nsteps = int(round(RUN_PS * 1000.0 / TIME_STEP_FS))
    all_runs = []

    for cif_name in CIFS:
        cif_path = os.path.join(ROOT, cif_name)
        tag = cif_name.replace(".cif", "")
        a, b, c, atoms = parse_cif(cif_path)
        counts = Counter(e for e, _, _, _ in atoms)
        xyz_path = os.path.join(models_dir, tag + ".xyz")
        write_xyz(xyz_path, a, b, c, atoms)
        print("=== %s ===" % tag)
        print("  box = %.4f x %.4f x %.4f A" % (a, b, c))
        print("  natoms = %d  %s" % (len(atoms), dict(counts)))
        print("  model = %s" % xyz_path)

        for temp in TEMPS:
            rundir = os.path.join(ROOT, tag, "%dK" % temp)
            os.makedirs(rundir, exist_ok=True)
            model_link = os.path.join(rundir, "model.xyz")
            if os.path.islink(model_link) or os.path.exists(model_link):
                os.remove(model_link)
            os.symlink(os.path.relpath(xyz_path, rundir), model_link)
            write_run_in(os.path.join(rundir, "run.in"), tag, temp)
            jobname = "f4-%s-%d" % (short_tag(cif_name), temp)
            write_submit(os.path.join(rundir, "submit.sh"), jobname, rundir)
            all_runs.append(rundir)

    submit_all = os.path.join(ROOT, "submit_all.sh")
    with open(submit_all, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Submit all 10 FIG4 NEP89 jobs. Review first; do not run blindly.\n")
        f.write("set -e\n")
        f.write("ROOT=%s\n" % ROOT)
        for rundir in all_runs:
            rel = os.path.relpath(rundir, ROOT)
            f.write('echo "sbatch %s"\n' % rel)
            f.write("( cd \"$ROOT/%s\" && sbatch submit.sh )\n" % rel)
    os.chmod(submit_all, 0o755)

    print("=== protocol ===")
    print("  potential = NEP89 (%s)" % POTENTIAL)
    print("  ensemble  = nvt_lan  T  T  100")
    print("  time_step = %g fs" % TIME_STEP_FS)
    print("  run       = %d steps = %g ps" % (nsteps, RUN_PS))
    print("  dump_xyz  = every %d steps (%.3f ps), %d frames" % (
        DUMP_EVERY, DUMP_EVERY * TIME_STEP_FS / 1000.0, nsteps // DUMP_EVERY))
    print("  jobs      = %d (not submitted)" % len(all_runs))
    for rundir in all_runs:
        print("    %s" % os.path.relpath(rundir, ROOT))


if __name__ == "__main__":
    main()
