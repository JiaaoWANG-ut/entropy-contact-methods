#!/bin/bash
# Submit all example NEP89 jobs. Review first; do not run blindly.
set -e
ROOT=/public/home/jwang/test/testgpumd/production-run/runs
echo "sbatch NiO-6-Fe-Zn-Cr-Ru/1500K"
( cd "$ROOT/NiO-6-Fe-Zn-Cr-Ru/1500K" && sbatch submit.sh )
echo "sbatch NiO-6-Fe-Zn-Cr-Ru/1800K"
( cd "$ROOT/NiO-6-Fe-Zn-Cr-Ru/1800K" && sbatch submit.sh )
echo "sbatch NiO-6-Fe-Zn-Cr-Ru/2000K"
( cd "$ROOT/NiO-6-Fe-Zn-Cr-Ru/2000K" && sbatch submit.sh )
echo "sbatch NiO-6-Fe-Zn-Cr-Ru/2500K"
( cd "$ROOT/NiO-6-Fe-Zn-Cr-Ru/2500K" && sbatch submit.sh )
echo "sbatch NiO-6-Fe-Zn-Cr-Ru/3000K"
( cd "$ROOT/NiO-6-Fe-Zn-Cr-Ru/3000K" && sbatch submit.sh )
echo "sbatch NiO-8-Fe-Zn-Cr-Ru/1500K"
( cd "$ROOT/NiO-8-Fe-Zn-Cr-Ru/1500K" && sbatch submit.sh )
echo "sbatch NiO-8-Fe-Zn-Cr-Ru/1800K"
( cd "$ROOT/NiO-8-Fe-Zn-Cr-Ru/1800K" && sbatch submit.sh )
echo "sbatch NiO-8-Fe-Zn-Cr-Ru/2000K"
( cd "$ROOT/NiO-8-Fe-Zn-Cr-Ru/2000K" && sbatch submit.sh )
echo "sbatch NiO-8-Fe-Zn-Cr-Ru/2500K"
( cd "$ROOT/NiO-8-Fe-Zn-Cr-Ru/2500K" && sbatch submit.sh )
echo "sbatch NiO-8-Fe-Zn-Cr-Ru/3000K"
( cd "$ROOT/NiO-8-Fe-Zn-Cr-Ru/3000K" && sbatch submit.sh )
