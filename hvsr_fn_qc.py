"""
Quality-controlled global Fn CSV from combined HVSR CSVs.

Spread metric: σ_ln = std( ln(Fn,i accepted) )  — lognormal std, consistent
               with hvsrpy's own statistics framework.

The ±1σ_ln bounds in frequency space are:
    lower = exp( μ_ln - σ_ln ) = Fn,mu_geo * exp(-σ_ln)
    upper = exp( μ_ln + σ_ln ) = Fn,mu_geo * exp(+σ_ln)

Equivalently, σ_ln / ln(2) is the half-spread in octaves.

For S7232-style stations where Fn,i accepted spans 0.1–1 Hz (≈1 decade),
σ_ln ≈ 0.65, which triggers the rejection threshold below.

Quality flags
  "good"        σ_ln ≤ WARN_SIGMA_LN
  "wide_spread" WARN_SIGMA_LN < σ_ln ≤ MAX_SIGMA_LN
  "rejected"    σ_ln > MAX_SIGMA_LN  →  fn_hz set to NaN

Outputs
  <output_dir>/fn_qc_summary.csv    – global table with QC columns
  <output_dir>/fn_qc_histograms/    – one histogram per station

Usage:
    python hvsr_fn_qc.py [stations_csv]
"""

import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hvsrpy
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
output_dir = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")
combined_dir = output_dir / "combined"
hist_dir     = output_dir / "fn_qc_histograms"
hvsr_approved_dir = output_dir / "fn_qc_hvsr" / "approved"
hvsr_rejected_dir = output_dir / "fn_qc_hvsr" / "rejected"
for d in (hist_dir, hvsr_approved_dir, hvsr_rejected_dir):
    d.mkdir(parents=True, exist_ok=True)

stations_csv = Path(sys.argv[1]) if len(sys.argv) > 1 else \
               Path("vulcano_stations.csv")

# ── Quality thresholds (σ_ln units) ──────────────────────────────────────────
# σ_ln  →  half-spread in octaves (σ_ln / ln2)
#   0.3  →  0.43 oct   (tight)
#   0.5  →  0.72 oct   (moderate)
#   0.7  →  1.01 oct   (wide — S7232-like spanning ~1 decade triggers this)
WARN_SIGMA_LN = 0.3   # flag as "wide_spread" above this
MAX_SIGMA_LN  = 0.5   # flag as "rejected"   above this → fn_hz = NaN

# ── Load station coordinates (station, Y, X, Z — no header) ──────────────────
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

print(f"Found {len(combined_csvs)} station(s)\n"
      f"Thresholds: warn σ_ln>{WARN_SIGMA_LN:.2f}  "
      f"reject σ_ln>{MAX_SIGMA_LN:.2f}\n")

global_rows = []

for csv_path in combined_csvs:
    sta = csv_path.stem.replace("_WR", "")

    hvsr = hvsrpy.object_io.read_hvsr_object_from_file(str(csv_path))

    fn_i      = hvsr.peak_frequencies   # accepted individual peaks
    n_accepted = len(fn_i)
    n_total    = hvsr.n_curves

    # ── Lognormal statistics of Fn,i ─────────────────────────────────────────
    if n_accepted >= 2:
        # Use hvsrpy's own lognormal statistics for consistency:
        #   mean_fn_frequency → geometric mean = exp(μ_ln)
        #   std_fn_frequency  → σ_ln = std(ln(fn_i))  [no exp applied]
        fn_mu_geo = hvsr.mean_fn_frequency(distribution="lognormal")
        sigma_ln  = hvsr.std_fn_frequency(distribution="lognormal")
        fn_lower  = fn_mu_geo * np.exp(-sigma_ln)   # μ_ln − σ_ln
        fn_upper  = fn_mu_geo * np.exp(+sigma_ln)   # μ_ln + σ_ln
        half_spread_oct = sigma_ln / np.log(2)
    else:
        fn_mu_geo = fn_lower = fn_upper = np.nan
        sigma_ln  = np.inf
        half_spread_oct = np.inf

    # ── Fn of post-rejection mean curve ──────────────────────────────────────
    try:
        fn_mc, fn_mc_amp = hvsr.mean_curve_peak(distribution="lognormal")
    except ValueError:
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

    # ── Per-station histogram ─────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 4))

    if n_accepted > 0:
        lo   = max(fn_i.min() * 0.7, 1e-3)
        hi   = fn_i.max() * 1.3
        bins = np.geomspace(lo, hi, 40)
        ax.hist(fn_i, bins=bins, color="steelblue", alpha=0.7,
                label=f"Fn,i accepted  (n={n_accepted}/{n_total})")
        ax.set_xscale("log")

    # Shade the μ_ln ± σ_ln band
    if np.isfinite(fn_lower) and np.isfinite(fn_upper):
        ax.axvspan(fn_lower, fn_upper, alpha=0.15, color="orange",
                   label=f"μ_ln ± σ_ln  (σ_ln={sigma_ln:.3f})")

    # Reference vertical lines
    for val, label, color, ls in [
        (fn_mc,     "Fn,mc (mean curve peak)", "red",    "-"),
        (fn_mu_geo, "Fn,i  geometric mean",    "green",  "--"),
        (fn_lower,  "μ_ln − σ_ln",             "orange", ":"),
        (fn_upper,  "μ_ln + σ_ln",             "orange", ":"),
    ]:
        if np.isfinite(val):
            ax.axvline(val, color=color, linestyle=ls, linewidth=1.5, label=label)

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Count")
    ax.set_title(f"{sta}  —  σ_ln = {sigma_ln:.3f}  "
                 f"(±{half_spread_oct:.2f} oct)  [{flag}]")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(hist_dir / f"{sta}_fn_histogram.png", dpi=150)
    plt.close(fig)

    # ── Single-panel HVSR curve plot, sorted by quality flag ─────────────────
    hvsr_out_dir = hvsr_rejected_dir if flag == "rejected" else hvsr_approved_dir
    try:
        sfig, sax = hvsrpy.plot_single_panel_hvsr_curves(hvsr)
        # Annotate with QC info
        sax.set_title(f"{sta}  —  σ_ln={sigma_ln:.3f}  [{flag}]")
        if np.isfinite(fn_mc):
            sax.axvline(fn_mc, color="red", linestyle="-", linewidth=1.5,
                        label=f"Fn,mc = {fn_mc:.3f} Hz")
        if np.isfinite(fn_lower) and np.isfinite(fn_upper):
            sax.axvspan(fn_lower, fn_upper, alpha=0.12, color="orange",
                        label=f"μ_ln ± σ_ln")
        sax.get_legend().remove()
        sax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=7)
        sfig.savefig(hvsr_out_dir / f"{sta}_hvsr.png", dpi=150,
                     bbox_inches="tight")
        plt.close(sfig)
    except Exception as e:
        print(f"    [warn] could not plot HVSR for {sta}: {e}")

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
