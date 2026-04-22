"""
Quality-controlled global Fn CSV from combined HVSR CSVs.

Spread metric: σ_ln = std( ln(Fn,i accepted) )  — lognormal std, consistent
               with hvsrpy's own statistics framework.

The ±1σ_ln bounds in frequency space are:
    lower = exp( μ_ln - σ_ln ) = Fn,mu_geo * exp(-σ_ln)
    upper = exp( μ_ln + σ_ln ) = Fn,mu_geo * exp(+σ_ln)

Equivalently, σ_ln / ln(2) is the half-spread in octaves.

Performance
  The combined CSVs are ~264 MB each (one column per window).
  np.loadtxt on such files is very slow; this script uses two strategies:
    1. pandas.read_csv  — 5–10× faster than np.loadtxt for wide files
    2. .npz cache       — saves peak_frequencies + mean/std curves after the
                          first load; subsequent runs skip the 264 MB read entirely.
  Delete <combined_dir>/cache/<sta>.npz to force a re-read.

Quality flags
  "good"        σ_ln ≤ WARN_SIGMA_LN
  "wide_spread" WARN_SIGMA_LN < σ_ln ≤ MAX_SIGMA_LN
  "rejected"    σ_ln > MAX_SIGMA_LN  →  fn_hz set to NaN

Outputs
  <output_dir>/fn_qc_summary.csv        – global table with QC columns
  <output_dir>/fn_qc_histograms/        – one histogram per station
  <output_dir>/fn_qc_hvsr/approved/     – HVSR curves for passing stations
  <output_dir>/fn_qc_hvsr/rejected/     – HVSR curves for rejected stations

Usage:
    python hvsr_fn_qc.py [stations_csv]
"""

import sys
import csv
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hvsrpy
from hvsrpy.hvsr_traditional import HvsrTraditional
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
output_dir        = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")
combined_dir      = output_dir / "combined"
cache_dir         = combined_dir / "cache"
hist_dir          = output_dir / "fn_qc_histograms"
hvsr_approved_dir = output_dir / "fn_qc_hvsr" / "approved"
hvsr_rejected_dir = output_dir / "fn_qc_hvsr" / "rejected"
for d in (cache_dir, hist_dir, hvsr_approved_dir, hvsr_rejected_dir):
    d.mkdir(parents=True, exist_ok=True)

stations_csv = Path(sys.argv[1]) if len(sys.argv) > 1 else \
               Path("vulcano_stations.csv")

# ── Quality thresholds (σ_ln units) ──────────────────────────────────────────
WARN_SIGMA_LN = 0.3
MAX_SIGMA_LN  = 0.5

# ── Helpers ───────────────────────────────────────────────────────────────────

def _read_header(fname):
    """Return the JSON meta dict and n_data_cols from a hvsrpy CSV."""
    header_lines = []
    n_data_cols  = None
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line[2:])
            elif n_data_cols is None:
                n_data_cols = len(line.split(","))
    meta = json.loads("\n".join(header_lines[:-1]))
    return meta, n_data_cols


def _fast_read_hvsr(fname):
    """Load a hvsrpy CSV using pandas (much faster than np.loadtxt)."""
    meta, _ = _read_header(fname)
    # pandas is 5-10x faster than np.loadtxt on wide files
    array = pd.read_csv(fname, comment="#", header=None).values
    hvsr = HvsrTraditional(array[:, 0], array[:, 1:-2].T)
    hvsr.meta = meta
    hvsr.update_peaks_bounded(
        search_range_in_hz=tuple(meta["search_range_in_hz"]),
        find_peaks_kwargs=meta["find_peaks_kwargs"],
    )
    hvsr.valid_window_boolean_mask = np.array(meta.pop("valid_window_boolean_mask"))
    hvsr.valid_peak_boolean_mask   = np.array(meta.pop("valid_peak_boolean_mask"))
    return hvsr


def _load_or_cache(csv_path, cache_path):
    """
    Return (frequency, mean_curve, std_curve, fn_i, n_total, search_range).
    Uses the .npz cache when available; otherwise parses the full CSV and
    writes the cache for next time.
    """
    if cache_path.exists():
        d = np.load(cache_path, allow_pickle=True)
        return (d["frequency"], d["mean_curve"], d["std_curve"],
                d["fn_i"], int(d["n_total"]), tuple(d["search_range"]))

    # ── First run: parse full 264 MB CSV with pandas ──────────────────────────
    meta, n_cols = _read_header(csv_path)
    #print(f"    (cache miss — reading {csv_path.stat().st_size/1e6:.0f} MB CSV with pandas)")

    # Read only freq (col 0) and last 2 cols (mean, std) — fast, low memory
    light = pd.read_csv(csv_path, comment="#", header=None,
                        usecols=[0, n_cols - 2, n_cols - 1]).values
    frequency  = light[:, 0]
    mean_curve = light[:, 1]
    std_curve  = light[:, 2]

    # Read full array to compute per-window peaks (unavoidable once)
    full = pd.read_csv(csv_path, comment="#", header=None).values
    hvsr = HvsrTraditional(full[:, 0], full[:, 1:-2].T)
    hvsr.meta = meta
    search_range = tuple(meta.get("search_range_in_hz", (None, None)))
    hvsr.update_peaks_bounded(search_range_in_hz=search_range,
                              find_peaks_kwargs=meta.get("find_peaks_kwargs"))
    hvsr.valid_peak_boolean_mask = np.array(meta.get("valid_peak_boolean_mask", [True]*hvsr.n_curves))
    fn_i    = hvsr.peak_frequencies
    n_total = hvsr.n_curves

    np.savez_compressed(cache_path,
                        frequency=frequency,
                        mean_curve=mean_curve,
                        std_curve=std_curve,
                        fn_i=fn_i,
                        n_total=np.array(n_total),
                        search_range=np.array(list(search_range), dtype=object))

    return frequency, mean_curve, std_curve, fn_i, n_total, search_range


def _plot_hvsr_lightweight(frequency, mean_curve, std_curve,
                            fn_mc, fn_lower, fn_upper, sigma_ln, flag, sta):
    """Plot mean ± std curves without loading all individual windows."""
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogx(frequency, mean_curve, color="black", linewidth=1.5,
                label="Mean curve")
    ax.fill_between(frequency, mean_curve - std_curve, mean_curve + std_curve,
                    alpha=0.25, color="steelblue", label="± std")
    if np.isfinite(fn_mc):
        ax.axvline(fn_mc, color="red", linestyle="-", linewidth=1.5,
                   label=f"Fn,mc = {fn_mc:.3f} Hz")
    if np.isfinite(fn_lower) and np.isfinite(fn_upper):
        ax.axvspan(fn_lower, fn_upper, alpha=0.12, color="orange",
                   label=f"μ_ln ± σ_ln")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("HVSR amplitude")
    ax.set_title(f"{sta}  —  σ_ln={sigma_ln:.3f}  [{flag}]")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=7)
    fig.tight_layout()
    return fig


# ── Load station coordinates ──────────────────────────────────────────────────
coords = {}
if stations_csv.exists():
    with open(stations_csv) as f:
        for row in csv.reader(f):
            if len(row) >= 4:
                coords[row[0].strip()] = (float(row[1]), float(row[2]), float(row[3]))

# ── Process each station ──────────────────────────────────────────────────────
combined_csvs = sorted(combined_dir.glob("*_WR.csv"))
if not combined_csvs:
    raise FileNotFoundError(f"No combined CSVs found in {combined_dir}")

print(f"Found {len(combined_csvs)} station(s)   "
      f"thresholds: warn σ_ln>{WARN_SIGMA_LN}  reject σ_ln>{MAX_SIGMA_LN}\n")

global_rows = []

for csv_path in combined_csvs:
    sta        = csv_path.stem.replace("_WR", "")
    cache_path = cache_dir / f"{sta}.npz"

    frequency, mean_curve, std_curve, fn_i, n_total, search_range = \
        _load_or_cache(csv_path, cache_path)

    n_accepted = len(fn_i)

    # ── Lognormal spread ──────────────────────────────────────────────────────
    if n_accepted >= 2:
        log_fn    = np.log(fn_i[fn_i > 0])
        mu_ln     = log_fn.mean()
        sigma_ln  = log_fn.std(ddof=1)
        fn_mu_geo = np.exp(mu_ln)
        fn_lower  = fn_mu_geo * np.exp(-sigma_ln)
        fn_upper  = fn_mu_geo * np.exp(+sigma_ln)
        half_spread_oct = sigma_ln / np.log(2)
    else:
        fn_mu_geo = fn_lower = fn_upper = np.nan
        sigma_ln  = np.inf
        half_spread_oct = np.inf

    # ── Fn of mean curve ─────────────────────────────────────────────────────
    try:
        from hvsrpy.hvsr_curve import HvsrCurve
        fn_mc, fn_mc_amp = HvsrCurve._find_peak_bounded(
            frequency, mean_curve,
            search_range_in_hz=search_range,
            find_peaks_kwargs=None,
        )
        if fn_mc is None:
            fn_mc, fn_mc_amp = np.nan, np.nan
    except Exception:
        fn_mc, fn_mc_amp = np.nan, np.nan

    # ── Quality flag ─────────────────────────────────────────────────────────
    if sigma_ln <= WARN_SIGMA_LN:
        flag = "good"
    elif sigma_ln <= MAX_SIGMA_LN:
        flag = "wide_spread"
    else:
        flag = "rejected"

    fn_out = fn_mc if flag != "rejected" else np.nan

    Y, X, Z = coords.get(sta, (np.nan, np.nan, np.nan))
    global_rows.append(dict(
        station=sta, Y=Y, X=X, Z=Z,
        fn_hz=fn_out, fn_mc=fn_mc, fn_mc_amplitude=fn_mc_amp,
        fn_mu_geo=fn_mu_geo,
        fn_lower_1sigma=fn_lower, fn_upper_1sigma=fn_upper,
        sigma_ln=sigma_ln, half_spread_oct=half_spread_oct,
        n_accepted=n_accepted, n_total=n_total,
        quality=flag,
    ))

    sym = {"good": "✓", "wide_spread": "~", "rejected": "✗"}[flag]
    print(f"  {sym} {sta:8s}  fn_mc={fn_mc:.3f} Hz  "
          f"μ_geo={fn_mu_geo:.3f}  σ_ln={sigma_ln:.3f}  "
          f"(±{half_spread_oct:.2f} oct)  [{flag}]")

    # ── Histogram ─────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4))
    if n_accepted > 0:
        lo   = max(fn_i.min() * 0.7, 1e-3)
        hi   = fn_i.max() * 1.3
        bins = np.geomspace(lo, hi, 40)
        ax.hist(fn_i, bins=bins, color="steelblue", alpha=0.7,
                label=f"Fn,i accepted  (n={n_accepted}/{n_total})")
        ax.set_xscale("log")
    if np.isfinite(fn_lower) and np.isfinite(fn_upper):
        ax.axvspan(fn_lower, fn_upper, alpha=0.15, color="orange",
                   label=f"μ_ln ± σ_ln  (σ_ln={sigma_ln:.3f})")
    for val, lbl, color, ls in [
        (fn_mc,     "Fn,mc (mean curve peak)", "red",   "-"),
        (fn_mu_geo, "Fn,i geometric mean",     "green", "--"),
        (fn_lower,  "μ_ln − σ_ln",             "orange", ":"),
        (fn_upper,  "μ_ln + σ_ln",             "orange", ":"),
    ]:
        if np.isfinite(val):
            ax.axvline(val, color=color, linestyle=ls, linewidth=1.5, label=lbl)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Count")
    ax.set_title(f"{sta}  —  σ_ln={sigma_ln:.3f}  (±{half_spread_oct:.2f} oct)  [{flag}]")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(hist_dir / f"{sta}_fn_histogram.png", dpi=150)
    plt.close(fig)

    # ── HVSR plot (lightweight: mean ± std only, no individual curves) ────────
    hvsr_out_dir = hvsr_rejected_dir if flag == "rejected" else hvsr_approved_dir
    sfig = _plot_hvsr_lightweight(frequency, mean_curve, std_curve,
                                   fn_mc, fn_lower, fn_upper, sigma_ln, flag, sta)
    sfig.savefig(hvsr_out_dir / f"{sta}_hvsr.png", dpi=150, bbox_inches="tight")
    plt.close(sfig)

# ── Write global CSV ──────────────────────────────────────────────────────────
global_csv = output_dir / "fn_qc_summary.csv"
fieldnames = ["station", "Y", "X", "Z",
              "fn_hz", "fn_mc", "fn_mc_amplitude",
              "fn_mu_geo", "fn_lower_1sigma", "fn_upper_1sigma",
              "sigma_ln", "half_spread_oct",
              "n_accepted", "n_total", "quality"]

def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6f}" if np.isfinite(v) else "nan"
    return v

with open(global_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for row in global_rows:
        w.writerow({k: _fmt(v) for k, v in row.items()})

flags = [r["quality"] for r in global_rows]
print(f"\nSummary:  good={flags.count('good')}  "
      f"wide_spread={flags.count('wide_spread')}  "
      f"rejected={flags.count('rejected')}")
print(f"Global CSV  →  {global_csv}")
print(f"Histograms  →  {hist_dir}/")
print(f"HVSR plots  →  {output_dir}/fn_qc_hvsr/")
