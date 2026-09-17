# General methods: entropy & metal–metal contact analysis

Reusable Python methods for post-processing **GPUMD** extended-XYZ
trajectories of multi-metal / oxide nanoparticle systems:

1. **Atomic entropy** $S_{\mathrm{atom}}$ and **configurational entropy**
   $S_{\mathrm{config}}$ vs time
2. **Contact ratio / contact-area** metrics for metal–metal wetting
3. Optional figure builders and OVITO helpers

The core algorithms live under `methods/`. Example structures and run-layout
helpers are included so you can reproduce an end-to-end workflow; swap in
your own systems by editing the `SYSTEMS` / temperature lists in the drivers.

---

## Repository layout

```
methods/
  traj_entropy.py          # Core: S_atom / S_config from dump.xyz
  run_entropy.py           # Batch entropy driver + CLI
  entropy_report.py        # Entropy tables & plots
  run_contact.py           # Metal–Ni contact ratio / area
  nature_contact_figures.py
  stopping_point.py        # Morphological stage / stopping-point figure
  send_report.py           # Optional: email result zip
  check_sent.py            # Optional: IMAP Sent check
  runs/
    setup_runs.py          # CIF → model.xyz + per-T run.in / submit scripts
    submit_all.sh
    models/*.xyz           # Example initial structures
    *.cif
  ovito-figures/           # Frame extract / Tachyon render / stitch scripts
```

Trajectory dumps (`runs/<system>/<T>K/dump.xyz`) are **not** shipped.
Place them locally in that layout (or point the drivers at your paths)
before running analysis.

---

## Requirements

```bash
pip install -r requirements.txt
```

Python 3.10+ recommended. GPUMD is not required once `dump.xyz` files exist.

---

## Methods

### Atomic entropy $S_{\mathrm{atom}}$

Each atom is labelled by motif `(species, CN_O, CN_M)`:

- `CN_O`: O neighbours within **2.5 Å**
- `CN_M`: metal neighbours within **3.0 Å**

$$
S_{\mathrm{atom}} = -\sum_k p_k \ln p_k \quad (\mathrm{nats}),\quad
p_k = N_k / N_{\mathrm{atoms}}
$$

### Configurational entropy $S_{\mathrm{config}}$

Metal atoms (O excluded) form clusters with contacts &lt; **3.0 Å** (PBC).

$$
S_{\mathrm{config}} = -\sum_c p_c \ln p_c,\quad
p_c = n_c / N_{\mathrm{metal}}
$$

### Contact / wetting metrics

For each secondary metal X (example set: Fe, Zn, Cr, Ru), an X atom contacts
Ni if it has ≥1 Ni neighbour within 3.0 Å. Contact ratio
$\rho = n_{\mathrm{contact}}/N_X$; contact area uses a close-packed atomic
area proxy. Frames with high Ni coordination (alloyed / fused) are flagged
so area is not misread as wetting.

Default sampling: **400** uniformly spaced frames for entropy; **log-spaced**
frames (denser in 0–1 ns) for contact, window 0–5 ns.

These cutoffs and element lists are constants at the top of the modules —
change them for other chemistries.

---

## Usage

```bash
cd methods
```

### 0. (Optional) Prepare GPUMD run directories

Edit the `POTENTIAL` path in `runs/setup_runs.py` for your machine, then:

```bash
python runs/setup_runs.py
# then submit via runs/submit_all.sh (cluster-specific)
```

Expected dump path after MD:

```text
runs/<system>/<temp>K/dump.xyz
```

### 1. Entropy vs time

```bash
python run_entropy.py --workers 10
```

First 5 ns only (writes to `entropy_analysis_5ns/`):

```bash
python run_entropy.py --tag 5ns --max-ps 5000 --workers 10
```

| Flag | Meaning |
|------|---------|
| `--workers N` | Parallel processes (default 10) |
| `--force` | Recompute even if `.npz` caches exist |
| `--plot-only` | Rebuild CSV/XLSX/figures from existing caches |
| `--max-ps X` | Analyse only the first `X` picoseconds |
| `--tag NAME` | Output dir `entropy_analysis_NAME/` |

**Outputs** under `entropy_analysis/` or `entropy_analysis_<tag>/`:

- `data/<system>_<T>K.npz`
- `entropy_raw_long.csv`, wide CSVs, `entropy_tables.xlsx`, PNG/SVG figures

Frame byte-offset caches live in `entropy_frame_index/` (safe to delete).

### 2. Contact ratio & contact area

```bash
python run_contact.py --workers 10
# or
python run_contact.py --plot-only
```

**Outputs** under `contact_analysis/`.

### 3. Publication-style contact figures

```bash
python nature_contact_figures.py
```

Writes plots, source Excel, and captions to `contact_analysis/nature/`.

### 4. Morphological stopping-point figure

Needs contact caches, an entropy CSV (default
`entropy_analysis_5ns/entropy_raw_long.csv`), and snapshot PNGs.
Edit `PNG_ROOT` in `stopping_point.py` if needed.

```bash
python stopping_point.py
```

### 5. OVITO rendering (optional)

See `ovito-figures/README.md` for extract / Tachyon / stitch scripts.

### 6. Email a result bundle (optional)

```bash
export GMAIL_APP_PASSWORD='xxxx'
python send_report.py --dir entropy_analysis_5ns
```

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
       └──────────►  stopping_point.py
```

---

## Cutoffs & constants

| Symbol | Value | Role |
|--------|-------|------|
| `O_CUT` | 2.5 Å | First O shell (`traj_entropy`) |
| `MET_CUT` | 3.0 Å | Metal–metal contact / cluster link |
| `AREA_PER_ATOM` | ~5.41 Å² | Close-packed area for $d_{\mathrm{nn}}=2.50$ Å |
| `ALLOY_CN` | 4 | Ni CN ≥ 4 → fused / alloyed flag |
| `N_SAMPLE` | 400 | Entropy frame subsample |

---

## Example systems shipped with the helpers

| Systems | Temperatures | Protocol |
|---------|--------------|----------|
| `NiO-6-Fe-Zn-Cr-Ru`, `NiO-8-Fe-Zn-Cr-Ru` | 1500–3000 K | NEP89, NVT-Langevin, `dt = 1 fs`, dump every 0.4 ps |

Adapt `SYSTEMS` and `TEMPERATURES` in `run_entropy.py` / `run_contact.py` for
your own runs.
