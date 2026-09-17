#!/bin/bash
#SBATCH --job-name=ovito-fig4
#SBATCH --partition=debug
#SBATCH --nodelist=dcu4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/public/home/jwang/test/testgpumd/ovito-figures/slurm_ovito_%j.out
#SBATCH --error=/public/home/jwang/test/testgpumd/ovito-figures/slurm_ovito_%j.err

set -e
echo "host=$(hostname) os=$(grep PRETTY /etc/os-release | head -1)"
ldd --version | head -1
echo "libstdc++: $(strings /lib64/libstdc++.so.6 | grep GLIBCXX_3.4.2 | tail -3)"

export OVITO_MODE=script
export QT_QPA_PLATFORM=offscreen
export PYTHONPATH=/public/home/jwang/.local/opt/ovito-py:/public/home/jwang/.local/opt/ovito-venv/lib/python3.11/site-packages
export LD_LIBRARY_PATH=/public/home/jwang/.local/opt/ovito-py/PySide6/Qt/lib:/public/home/jwang/.local/opt/ovito-py/ovito/plugins:${LD_LIBRARY_PATH:-}
export XDG_RUNTIME_DIR=/tmp/runtime-jwang-ovito
mkdir -p "$XDG_RUNTIME_DIR"

PY=/public/home/jwang/.local/opt/ovito-venv/bin/python
"$PY" - << 'PY'
import ovito
from ovito.vis import TachyonRenderer
print("OVITO_OK", ovito.version, TachyonRenderer)
PY

cd /public/home/jwang/test/testgpumd/ovito-figures
"$PY" render_ovito_tachyon.py
"$PY" stitch.py
echo "OVITO_RENDER_DONE"
ls -lh summary/*.png
