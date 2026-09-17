# FIG4 NiO–Fe/Zn/Cr/Ru: entropy & contact analysis

Python tools for post-processing **GPUMD** extended-XYZ trajectories of
high-entropy oxide particle systems used in Nature Fig. 4–style wetting /
morphology analysis.

**Companion figure pack (OVITO PNGs):**  
[FIG4-NiO-FeZnCrRu-figures](https://github.com/JiaaoWANG-ut/FIG4-NiO-FeZnCrRu-figures)

| Systems | Temperatures | Protocol |
|---------|--------------|----------|
| `NiO-6-Fe-Zn-Cr-Ru`, `NiO-8-Fe-Zn-Cr-Ru` | 1500, 1800, 2000, 2500, 3000 K | NEP89, NVT-Langevin, `dt = 1 fs`, dump every 400 steps (0.4 ps) |

---

## Repository layout

```
FIG4-paper/
  traj_entropy.py          # Core: S_atom / S_config from dump.xyz
  run_entropy.py           # Batch entropy driver + CLI
  entropy_report.py        # Entropy tables & plots
  run_contact.py           # Metal–Ni contact ratio / area
  nature_contact_figures.py
  fig4_stopping_point.py   # Morphological stage / stopping-point figure
  send_report.py           # Optional: email result zip (needs Gmail app password)
  check_sent.py            # Optional: IMAP check of Sent folder
  FIG4-run/
    setup_runs.py          # CIF → model.xyz + per-T run.in / submit scripts
    submit_all.sh
    models/*.xyz           # Initial structures
    *.cif
  ovito-figures/           # Frame extract / Tachyon render / stitch scripts
```

Trajectory dumps (`FIG4-run/<system>/<T>K/dump.xyz`) are **not** shipped
(they are hundreds of GB). Place them locally in the layout above before running
analysis.

---

## Requirements

```bash
pip install -r requirements.txt
```

Python 3.10+ recommended. No GPUMD needed for analysis if `dump.xyz` files
already exist.

---

## Methods (short)

### Atomic entropy \(S_\mathrm{atom}\)

Each atom is labelled by motif `(species, CN_O, CN_M)`:

- `CN_O`: O neighbours within **2.5 Å**
- `CN_M`: metal neighbours within **3.0 Å**

\[
S_\mathrm{atom} = -\sum_k p_k \ln p_k \quad (\mathrm{nats}),\quad
p_k = N_k / N_\mathrm{atoms}
\]

### Configurational entropy \(S_\mathrm{config}\)

Metal atoms (O excluded) form clusters with contacts \< **3.0 Å** (PBC).

\[
S_\mathrm{config} = -\sum_c p_c \ln p_c,\quad
p_c = n_c / N_\mathrm{metal}
\]

### Contact / wetting metrics

For each metal X ∈ {Fe, Zn, Cr, Ru}: an X atom contacts Ni if it has ≥1 Ni
neighbour within 3.0 Å. Contact ratio \(\rho = n_\mathrm{contact}/N_X\);
contact area uses a close-packed atomic area proxy. Frames with high Ni
coordination (alloyed / fused) are flagged so area is not misread as wetting.

Entropy uses **400** uniformly spaced frames per trajectory. Contact uses
**log-spaced** frames denser in 0–1 ns (window 0–5 ns by default).

---

## Usage

All analysis commands are run from `FIG4-paper/`:

```bash
cd FIG4-paper
```

### 0. (Optional) Prepare GPUMD run directories

Edit the `POTENTIAL` path in `FIG4-run/setup_runs.py` for your machine, then:

```bash
python FIG4-run/setup_runs.py
# then submit jobs via FIG4-run/submit_all.sh (cluster-specific)
```

Expected dump path after MD:

```text
FIG4-run/<system>/<temp>K/dump.xyz
```

### 1. Entropy vs time

Full trajectories:

```bash
python run_entropy.py --workers 10
```

First 5 ns only (writes to `entropy_analysis_5ns/`):

```bash
python run_entropy.py --tag 5ns --max-ps 5000 --workers 10
```

Useful flags:

| Flag | Meaning |
|------|---------|
| `--workers N` | Parallel processes (default 10) |
| `--force` | Recompute even if `.npz` caches exist |
| `--plot-only` | Rebuild CSV/XLSX/figures from existing caches |
| `--max-ps X` | Analyse only the first `X` picoseconds |
| `--tag NAME` | Output dir `entropy_analysis_NAME/` |

**Outputs** (under `entropy_analysis/` or `entropy_analysis_<tag>/`):

- `data/<system>_<T>K.npz` — per-trajectory tables  
- `entropy_raw_long.csv`, `*_Satom_wide.csv`, `*_Sconfig_wide.csv`  
- `entropy_tables.xlsx`, PNG/SVG figures  

Frame byte-offset caches are shared in `entropy_frame_index/` (safe to delete;
they are rebuilt automatically).

### 2. Contact ratio & contact area

```bash
python run_contact.py --workers 10
# or
python run_contact.py --plot-only
```

**Outputs** under `contact_analysis/`:

- `data/*.npz`, `contact_raw_long.csv`, `contact_summary.csv`  
- PNG/SVG overview and per-system plots  

### 3. Nature-format contact figures

Requires step 2 results:

```bash
python nature_contact_figures.py
```

Writes to `contact_analysis/nature/` (Fig. 1–3 style plots, source Excel,
captions).

### 4. Morphological stopping-point figure (Fig. 4 panel)

Needs contact caches **and** entropy CSV (default:
`entropy_analysis_5ns/entropy_raw_long.csv`), plus OVITO snapshot PNGs.
Edit `PNG_ROOT` in `fig4_stopping_point.py` if your PNG path differs
(default points at a local `ovito-figures/png` tree).

```bash
python fig4_stopping_point.py
```

### 5. OVITO rendering (optional)

Scripts under `ovito-figures/` extract frames, render with Tachyon, and stitch
T–time grids. See `ovito-figures/README.md`. Pre-rendered packs are published
on the [figures release](https://github.com/JiaaoWANG-ut/FIG4-NiO-FeZnCrRu-figures/releases).

### 6. Email a result bundle (optional)

```bash
export GMAIL_APP_PASSWORD='xxxx'   # Google app password, not account password
python send_report.py --dir entropy_analysis_5ns
```

Uses IMAP to avoid duplicate sends on flaky egress. Inspect Sent with
`python check_sent.py`.

---

## Typical workflow

```text
GPUMD dumps ready
       │
       ▼
 run_entropy.py  ──►  entropy_analysis[_tag]/
       │
       ▼
 run_contact.py  ──►  contact_analysis/
       │
       ├──────────►  nature_contact_figures.py
       │
       └──────────►  fig4_stopping_point.py  (also needs entropy CSV + PNGs)
```

---

## Cutoffs & constants

| Symbol | Value | Role |
|--------|-------|------|
| `O_CUT` | 2.5 Å | First O shell (`traj_entropy`) |
| `MET_CUT` | 3.0 Å | Metal–metal contact / cluster link |
| `AREA_PER_ATOM` | ~5.41 Å² | Close-packed area for \(d_\mathrm{nn}=2.50\) Å |
| `ALLOY_CN` | 4 | Ni CN ≥ 4 → fused / alloyed flag |
| `N_SAMPLE` | 400 | Entropy frame subsample |

---

## License / citation

Analysis code for internal manuscript preparation. Cite the manuscript and the
NEP89 potential when publishing results derived from these trajectories.
