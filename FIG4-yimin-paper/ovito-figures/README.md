# FIG4 NiO–Fe/Zn/Cr/Ru figures

NEP89 NVT-Langevin snapshots for FIG4 NiO-6 / NiO-8 + Fe/Zn/Cr/Ru
(GPUMD on Hygon DCU). **Download PNGs from [Releases](https://github.com/JiaaoWANG-ut/FIG4-NiO-FeZnCrRu-figures/releases)** — this repo is for distribution only.

## Layout

- One summary grid per composition: **rows = T** (1500–3000 K), **cols = time** (0–5 ns)
- Each small frame is labeled with composition/elements and time
- OVITO Tachyon + AO, 3600×3600 HD frames, OVITO default colors; ortho top view (all 5 clusters)

## Coverage

| System | 0–5 ns |
|---|---|
| NiO-6-Fe-Zn-Cr-Ru | all 5 T |
| NiO-8-Fe-Zn-Cr-Ru | all 5 T |

Protocol: `dt = 1 fs`, `nvt_lan`, dump every 400 steps (0.4 ps), 10 ns target.
