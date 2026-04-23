"""Produce summary outputs from combined HVSR CSVs (hvsr_combine.py output).

For each station found in <combined_dir>:
  - Writes a simplified CSV: frequency | mean_curve | mean_curve_plus_1std
    (std is skipped for DFA when only 1 curve exists)

Writes one global CSV: station | Y | X | Z | fn_hz | fn_amplitude

Usage:
    python hvsr_summary.py [method] [stations_csv]

    method       : "traditional" (default) or "dfa"
    stations_csv : path to stations CSV (station,Y,X,Z, no header)
"""

import sys
import csv
import numpy as np
import hvsrpy
from pathlib import Path

args         = sys.argv[1:]
method       = "traditional"
stations_csv = Path("vulcano_stations.csv")

for a in args:
    if a in ("traditional", "dfa"):
        method = a
    else:
        stations_csv = Path(a)

output_dir   = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")
combined_dir = output_dir / f"combined_{method}"
summary_dir  = output_dir / f"summary_{method}"
summary_dir.mkdir(exist_ok=True)

print(f"Method       : {method}")
print(f"Combined dir : {combined_dir}")

# ── Station coordinates ───────────────────────────────────────────────────────
coords = {}
if stations_csv.exists():
    with open(stations_csv) as f:
        for row in csv.reader(f):
            if len(row) >= 4:
                coords[row[0].strip()] = (float(row[1]), float(row[2]), float(row[3]))

combined_csvs = sorted(combined_dir.glob("*_WR.csv"))
if not combined_csvs:
    raise FileNotFoundError(f"No combined CSVs found in {combined_dir}")

print(f"Found {len(combined_csvs)} station(s)\n")

global_rows = []

for csv_path in combined_csvs:
    sta = csv_path.stem.replace("_WR", "")
    print(f"  {sta}", end=" ... ")

    hvsr = hvsrpy.object_io.read_hvsr_object_from_file(str(csv_path))

    freq = hvsr.frequency
    mean = hvsr.mean_curve(distribution="lognormal")

    # std_curve requires > 1 window; for DFA combined it may be just ~30 days
    try:
        plus_1std = hvsr.nth_std_curve(n=1, distribution="lognormal")
        has_std = True
    except (ValueError, AttributeError):
        plus_1std = np.full_like(mean, np.nan)
        has_std = False

    # Per-station simplified CSV
    out_csv = summary_dir / f"{sta}_mean_curve.csv"
    header  = "frequency_hz,mean_curve,mean_curve_plus_1std"
    data    = np.column_stack([freq, mean, plus_1std])
    np.savetxt(out_csv, data, delimiter=",",
               header=header, comments="", fmt="%.6f")
    print(f"saved {out_csv.name}", end=" | ")

    # Peak of mean curve
    try:
        fn_hz, fn_amp = hvsr.mean_curve_peak(distribution="lognormal")
    except (ValueError, AttributeError):
        fn_hz, fn_amp = np.nan, np.nan

    Y, X, Z = coords.get(sta, (np.nan, np.nan, np.nan))
    global_rows.append([sta, Y, X, Z, fn_hz, fn_amp])
    print(f"fn={fn_hz:.3f} Hz{'  (no std)' if not has_std else ''}")

# ── Global CSV ────────────────────────────────────────────────────────────────
global_csv = output_dir / f"stations_fn_summary_{method}.csv"
with open(global_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["station", "Y", "X", "Z", "fn_hz", "fn_amplitude"])
    for row in global_rows:
        w.writerow([row[0],
                    f"{row[1]:.8f}", f"{row[2]:.8f}", f"{row[3]:.4f}",
                    f"{row[4]:.4f}", f"{row[5]:.4f}"])

print(f"\nGlobal CSV saved to {global_csv}")
print(f"Per-station curves saved to {summary_dir}/")
