# HVSR — Horizontal-to-Vertical Spectral Ratio Processing Toolkit

A collection of Python (and MATLAB) scripts for computing, filtering, visualising, and inverting **HVSR (Horizontal-to-Vertical Spectral Ratio)** curves from seismic ambient-noise recordings. The workflow is built on top of [**hvsrpy**](https://github.com/jpvantassel/hvsrpy) and integrates with [**OpenHVSR-Inversion**](https://github.com/s-gaffet/OpenHVSR) for subsurface modelling.

> ⚠️ **Note from the author:** *"It is a bit messy because it contains a mix of scripts to calculate HV using hvsrpy (also filtering, plotting, etc) and some scripts to run OpenHVSR-Inversion."*

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Usage](#usage)
  - [Single-Day Processing](#single-day-processing)
  - [Multi-Day Processing](#multi-day-processing)
  - [Window Rejection & QC](#window-rejection--qc)
  - [Combining Results](#combining-results)
  - [Summary & Statistics](#summary--statistics)
  - [Mapping Fundamental Frequencies](#mapping-fundamental-frequencies)
  - [Dispersion / DFA Analysis](#dispersion--dfa-analysis)
  - [Batch Processing (HPC)](#batch-processing-hpc)
- [Processing Parameters](#processing-parameters)
- [Dependencies](#dependencies)
- [Data Format](#data-format)
- [Output](#output)

---

## Overview

This toolkit processes three-component seismic data (Z, E, N channels in **MiniSEED** format) to derive HVSR curves. It is designed for large seismic deployments running on HPC clusters, but can be adapted to any local dataset.

The two main use-cases are:

1. **HVSR computation** — using `hvsrpy` for traditional spectral ratio analysis with Konno–Ohmachi smoothing, Tukey windowing, frequency-domain window rejection, and geometric-mean horizontal combination.
2. **Inversion** — preparing and running **OpenHVSR-Inversion** to infer 1-D shear-wave velocity (Vs) profiles from the HVSR curves.

---

## Repository Structure

```
HVSR/
├── hvsr_oneday.py           # Process one station, one day at a time (memory-efficient)
├── hvsr_multiday.py         # Process one station across all available days (batched)
├── hvsr_fn_qc.py            # Quality control on the HVSR fundamental frequency (fn)
├── hvsr_window_rejection.py # Standalone window rejection pass
├── hvsr_combine.py          # Merge per-day CSV results into a single file
├── hvsr_summary.py          # Statistical summary and visualisation of merged results
├── hvsr_plot_combined.py    # Multi-station combined plot
├── hvsr_map_fn.py           # Spatial map of the fundamental frequency fn
├── hvsr_DFA.py              # Dispersion / DFA (Diffuse Field Assumption) analysis
├── array.sh                 # SLURM array job launcher for HPC clusters
├── hvsr.yml                 # Conda environment specification
├── messy/                   # Scratch / experimental scripts
└── README.md
```

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/miichenry/HVSR.git
cd HVSR
```

### 2. Create the Conda environment

```bash
conda env create -f hvsr.yml
conda activate hvsr
```

The environment pins all dependencies, including:

| Package | Version |
|---------|---------|
| Python  | 3.9 |
| hvsrpy  | 2.0.0 |
| ObsPy   | 1.4.1 |
| NumPy   | 1.26.4 |
| SciPy   | 1.13.1 |
| Matplotlib | 3.9.2 |
| Pandas  | 2.2.3 |
| Numba   | 0.60.0 |

---

## Usage

### Single-Day Processing

Processes one station day-by-day, saving one CSV per day. Peak memory is bounded to a single day of data — ideal for large deployments.

```bash
python hvsr_oneday.py <STATION_NAME>
# or pass full CSV row:
python hvsr_oneday.py "A6575,lat,lon,elev"
```

**Reads from:** `<base_path>/<STATION>/DPZ/`, `DPE/`, `DPN/`  
**Writes to:** `<output_dir>/perday/<STATION>/<STATION>_<YEAR>_<JDAY>.csv`

Days that already have an output CSV are skipped automatically.

---

### Multi-Day Processing

Loads all available days for a station and processes them in configurable batches (default `BATCH_SIZE = 3`). Produces window-rejection diagnostics and a combined HVSR output.

```bash
python hvsr_multiday.py <STATION_NAME>
```

**Writes to:** `<output_dir>/output_multiday/<STATION>_WR.csv` and `<STATION>_WR.png`

To reduce memory pressure, lower `BATCH_SIZE` in the script (e.g. to `1` for very large daily files).

> [!WARNING]
> This scrip not work when batch loading 3 days at a time for a total of 30 days, OOM error when starting hvsrpy preprocessing, for multiple day hvsr processing please use aforementioned hvsr_oneday.py``` 

---

### Window Rejection & QC

```bash
python hvsr_window_rejection.py
python hvsr_fn_qc.py
```

Applies frequency-domain window rejection (`hvsrpy.frequency_domain_window_rejection`) and produces pre/post rejection comparison figures.

---

### Combining Results

Merges all per-day CSVs for a station (or across stations) into a single consolidated result file.

```bash
python hvsr_combine.py
```

---

### Summary & Statistics

Generates statistical summaries (mean, standard deviation, median HVSR curves) and plots.

```bash
python hvsr_summary.py
```

---

### Mapping Fundamental Frequencies

Produces a spatial map of the HVSR fundamental frequency *fn* across all stations.

```bash
python hvsr_map_fn.py
```

---

### Dispersion / DFA Analysis

Runs Diffuse Field Assumption (DFA) or dispersion-based analysis.

```bash
python hvsr_DFA.py
```

---

### Batch Processing (HPC)

For processing many stations in parallel on a SLURM cluster:

```bash
# Submit an array job (one job per station)
sbatch array.sh
```

Edit `array.sh` to set the station list, memory, and time limits appropriate for your cluster.

---

## Processing Parameters

Default settings applied in both `hvsr_oneday.py` and `hvsr_multiday.py`:

| Parameter | Value |
|-----------|-------|
| Detrend | `constant` (mean removal) |
| Window length | 60 s |
| Orientation | 0° from North |
| Bandpass filter | 0.2 – 30 Hz |
| Window taper | Tukey (width = 0.01) |
| Smoothing operator | Konno & Ohmachi (bandwidth = 40) |
| Frequency axis | 256 log-spaced points, 0.2 – 50 Hz |
| Horizontal combination | Geometric mean |
| Resampling (dissimilar dt) | Frequency-domain resampling |
| Window rejection (n) | 2 |

---

## Dependencies

### Core

- **[hvsrpy](https://github.com/jpvantassel/hvsrpy)** `v2.0.0` — HVSR computation engine
- **[ObsPy](https://github.com/obspy/obspy)** `v1.4.1` — MiniSEED reading and seismic data handling
- **[NumPy](https://numpy.org/)** / **[SciPy](https://scipy.org/)** — numerical processing
- **[Matplotlib](https://matplotlib.org/)** — plotting (non-interactive `Agg` backend for HPC)
- **[Pandas](https://pandas.pydata.org/)** — CSV I/O

### Optional / Inversion

- **OpenHVSR-Inversion** (MATLAB) — scripts in `messy/` prepare inputs for this tool

---

## Data Format

Input data is expected as **MiniSEED** files organised by station and component:

```
<base_path>/
└── <STATION>/
    ├── DPZ/   ← vertical component   (*.mseed)
    ├── DPE/   ← east component       (*.mseed)
    └── DPN/   ← north component      (*.mseed)
```

File naming convention (parsed by regex):

```
VN.<STATION>..<COMP>.D.<YEAR>.<JDAY>.00.00.00.000.mseed
```

The scripts discover days automatically by matching `(YEAR, JDAY)` keys across all three components.

---

## Output

Each processed station/day produces a **CSV** file written by `hvsrpy.object_io.write_hvsr_object_to_file`, containing:

- Frequency axis
- Per-window HVSR amplitudes
- Window validity mask (post-rejection)
- Statistical summary (mean, std, median curves)

Figures (`.png`) are saved alongside each CSV showing:

- Pre- and post-rejection window comparison
- Single-panel HVSR curve with confidence bounds

## Acknowledgements

- [hvsrpy](https://github.com/jpvantassel/hvsrpy) by J.P. Vantassel
- [OpenHVSR-Inversion](https://github.com/s-gaffet/OpenHVSR) by S. Gaffet et al.

---
