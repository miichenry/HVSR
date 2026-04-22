"""Produce summary outputs from combined HVSR CSVs (hvsr_combine.py output).

For each station found in <combined_dir>:
  - Writes a simplified CSV: frequency | mean_curve | mean_curve_plus_1std

Writes one global CSV: station | Y | X | Z | fn_hz | fn_amplitude

Usage:
    python hvsr_summary.py [stations_csv]

    stations_csv : optional path to stations CSV (station,Y,X,Z, no header)
                   default: /home/users/h/henrymi/jectpro/vulcano/vulcano_stations.csv
"""

import sys
import csv
import numpy as np
import hvsrpy
from pathlib import Path

output_dir = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")
combined_dir = output_dir / "combined"
summary_dir  = output_dir / "summary"
summary_dir.mkdir(exist_ok=True)

stations_csv = Path(sys.argv[1]) if len(sys.argv) > 1 else \
               Path("vulcano_stations.csv")

# Load station coordinates keyed by station name
coords = {}
if stations_csv.exists():
    with open(stations_csv) as f:
        for row in csv.reader(f):
            if len(row) >= 4:
                coords[row[0].strip()] = (float(row[1]), float(row[2]), float(row[3]))

combined_csvs = sorted(combined_dir.glob("*_WR.csv"))
if not combined_csvs:
    raise FileNotFoundError(f"No combined CSVs found in {combined_dir}")

print(f"Found {len(combined_csvs)} station(s)")

global_rows = []

for csv_path in combined_csvs:
    sta = csv_path.stem.replace("_WR", "")
    print(f"  {sta}", end=" ... ")

    hvsr = hvsrpy.object_io.read_hvsr_object_from_file(str(csv_path))

    freq      = hvsr.frequency
    mean      = hvsr.mean_curve(distribution="lognormal")
    plus_1std = hvsr.nth_std_curve(n=1, distribution="lognormal")

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
    except ValueError:
        fn_hz, fn_amp = np.nan, np.nan

    Y, X, Z = coords.get(sta, (np.nan, np.nan, np.nan))
    global_rows.append([sta, Y, X, Z, fn_hz, fn_amp])
    print(f"fn={fn_hz:.3f} Hz")

print(f"Per-station curves saved to {summary_dir}/")
