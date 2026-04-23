"""Combine per-day HVSR CSVs for one station, run window rejection, save results.

Supports both processing methods produced by hvsr_oneday.py:

  traditional  — per-day CSVs are HvsrTraditional (many windows/day).
                 All window amplitudes are concatenated into one combined
                 HvsrTraditional and window rejection is applied.

  dfa          — per-day CSVs are HvsrDiffuseField (one single curve/day).
                 Each daily curve becomes one row in a combined HvsrTraditional
                 (n_curves = n_days), giving mean/std across days.
                 Window rejection is skipped (too few points to be meaningful).

The processing_method is auto-detected from the per-day CSV headers.

Input:  <output_dir>/perday_<method>/<sta>/<sta>_<year>_<jday>.csv
Output: <output_dir>/combined_<method>/<sta>_WR.csv  +  _single_panel.png

Usage:
    python hvsr_combine.py <station_arg> [method]
    method: "traditional" (default) or "dfa"
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime
import hvsrpy
from hvsrpy.hvsr_traditional import HvsrTraditional
from pathlib import Path
import sys

plt.style.use(hvsrpy.HVSRPY_MPL_STYLE)

sta    = sys.argv[1].split(',')[0]
method = sys.argv[2].lower() if len(sys.argv) > 2 else "traditional"

if method not in ("traditional", "dfa"):
    raise ValueError(f"method must be 'traditional' or 'dfa', got '{method}'")

output_dir   = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")
perday_dir   = output_dir / f"perday_{method}" / sta
combined_dir = output_dir / f"combined_{method}"
os.makedirs(combined_dir, exist_ok=True)

print(f"Station : {sta}")
print(f"Method  : {method}")
print(f"Starting: {datetime.now().strftime('%H:%M:%S')}")

csv_files = sorted(perday_dir.glob(f"{sta}_*.csv"))
if not csv_files:
    raise FileNotFoundError(f"No per-day CSVs found in {perday_dir}")

print(f"Found {len(csv_files)} per-day CSV(s)")


def _detect_processing_method(csv_path):
    """Read the JSON meta from the CSV header and return processing_method."""
    header_lines = []
    with open(csv_path) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line[2:])
            else:
                break
    meta = json.loads("\n".join(header_lines[:-1]))
    return meta.get("processing_method", "traditional")


# ── Load and stack per-day amplitudes ─────────────────────────────────────────
detected = _detect_processing_method(csv_files[0])
print(f"Detected per-day format: {detected}\n")

frequency       = None
amplitude_parts = []
meta_ref        = None

for csv in csv_files:
    h = hvsrpy.object_io.read_hvsr_object_from_file(str(csv))
    if frequency is None:
        frequency = h.frequency
        meta_ref  = h.meta

    if detected == "diffuse_field":
        # HvsrDiffuseField: amplitude is 1D (n_freqs,) — one curve per day
        amplitude_parts.append(h.amplitude[np.newaxis, :])   # → (1, n_freqs)
        print(f"  {csv.name}: 1 daily curve")
    else:
        # HvsrTraditional: amplitude is 2D (n_windows, n_freqs)
        amplitude_parts.append(h.amplitude)
        print(f"  {csv.name}: {h.n_curves} windows")

combined_amplitude = np.concatenate(amplitude_parts, axis=0)  # (N, n_freqs)
print(f"\nTotal {'days' if detected == 'diffuse_field' else 'windows'}: "
      f"{combined_amplitude.shape[0]}")

# Always produce an HvsrTraditional for the combined object so that the
# rest of the pipeline (summary, fn_qc, map) is uniform.
hvsr = HvsrTraditional(frequency, combined_amplitude, meta=meta_ref)

# ── Window rejection ──────────────────────────────────────────────────────────
# Skip for DFA: with only ~30 daily curves the 2-sigma rejection is unreliable.
if detected == "diffuse_field":
    print("Skipping window rejection for DFA (n_curves = n_days, too few for "
          "reliable 2-sigma rejection).")
else:
    n = 2
    search_range_in_hz = (None, None)
    _ = hvsrpy.frequency_domain_window_rejection(
        hvsr, n=n, search_range_in_hz=search_range_in_hz)

# ── Save ──────────────────────────────────────────────────────────────────────
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
sfig.savefig(fname, bbox_inches="tight")
plt.close()
print(f"Figure saved to {fname}")

print(f"\nDone: {datetime.now().strftime('%H:%M:%S')}")
