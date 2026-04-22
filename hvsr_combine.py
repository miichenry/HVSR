"""Combine per-day HVSR CSVs for one station, run window rejection, save results.

Expects per-day CSVs produced by hvsr_oneday.py at:
    <output_dir>/perday/<sta>/<sta>_<year>_<jday>.csv

Usage:
    python hvsr_combine.py <station_arg>
    station_arg: station name or CSV row like "A6575,lat,lon,elev"
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import hvsrpy
from pathlib import Path
import sys

plt.style.use(hvsrpy.HVSRPY_MPL_STYLE)

sta = sys.argv[1].split(',')[0]
output_dir = Path("/srv/beegfs/scratch/users/h/henrymi/vulcano/output_test")
perday_dir = output_dir / "perday" / sta
combined_dir = output_dir / "combined"
os.makedirs(combined_dir, exist_ok=True)

print(f"Station: {sta}")
print(f"Starting time: {datetime.now().strftime('%H:%M:%S')}")

csv_files = sorted(perday_dir.glob(f"{sta}_*.csv"))
if not csv_files:
    raise FileNotFoundError(f"No per-day CSVs found in {perday_dir}")

print(f"Found {len(csv_files)} per-day CSV(s)")

# Load each per-day HVSR object and collect the amplitude windows.
# read_hvsr_object_from_file reconstructs a full HvsrTraditional from disk,
# so no seismic records are needed — memory stays low.
frequency = None
amplitude_parts = []
for csv in csv_files:
    h = hvsrpy.object_io.read_hvsr_object_from_file(str(csv))
    if frequency is None:
        frequency = h.frequency
        meta = h.meta
    amplitude_parts.append(h.amplitude)
    print(f"  {csv.name}: {h.n_curves} windows")

combined_amplitude = np.concatenate(amplitude_parts, axis=0)
print(f"\nTotal windows: {combined_amplitude.shape[0]}")

hvsr = hvsrpy.HvsrTraditional(frequency, combined_amplitude, meta=meta)

n = 2
search_range_in_hz = (None, None)
_ = hvsrpy.frequency_domain_window_rejection(hvsr, n=n, search_range_in_hz=search_range_in_hz)

fname_prefix = str(combined_dir / f"{sta}_WR")

fname = f"{fname_prefix}.csv"
hvsrpy.object_io.write_hvsr_object_to_file(hvsr, fname)
print(f"Combined CSV saved to {fname}")

print("\nStatistical Summary:")
print("-"*20)
hvsrpy.summarize_hvsr_statistics(hvsr)

(sfig, ax) = hvsrpy.plot_single_panel_hvsr_curves(hvsr)
ax.get_legend().remove()
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
fname = f"{fname_prefix}_single_panel.png"
sfig.savefig(fname)
plt.close()
print(f"Figure saved to {fname}")

print(f"\nDone: {datetime.now().strftime('%H:%M:%S')}")
