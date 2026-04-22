import pathlib
import os
import re
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
base_path = Path("/srv/beegfs/scratch/shares/cdff/VulcaNODES/2021/VN")
output_dir = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")

# Maximum number of days loaded into memory at once.
# Reduce to 1 if individual daily files are very large (> ~1 GB each).
BATCH_SIZE = 3

print(f"Station: {sta}")
print(f"Starting time: {datetime.now().strftime('%H:%M:%S')}")

sta_path = base_path / sta
if not sta_path.exists():
    raise FileNotFoundError(f"Station directory not found: {sta_path}")


def index_by_day(component_dir):
    day_map = {}
    for f in sorted(component_dir.glob("*.mseed")):
        # VN.A0000..DPE.D.2021.290.00.00.00.000.mseed
        m = re.search(r'\.(\d{4})\.(\d{3})\.', f.name)
        if m:
            key = (m.group(1), m.group(2))
            day_map[key] = str(f)
    return day_map


z_map = index_by_day(sta_path / "DPZ")
e_map = index_by_day(sta_path / "DPE")
n_map = index_by_day(sta_path / "DPN")

common_days = sorted(set(z_map) & set(e_map) & set(n_map))
if not common_days:
    raise RuntimeError(f"No days with complete Z/E/N data found for station {sta}")

print(f"Days with complete data: {len(common_days)}")
for year, jday in common_days:
    print(f"  {year}.{jday}")

all_fnames = [[z_map[d], e_map[d], n_map[d]] for d in common_days]

preprocessing_settings = hvsrpy.settings.HvsrPreProcessingSettings()
preprocessing_settings.detrend = "constant"
preprocessing_settings.window_length_in_seconds = 60
preprocessing_settings.orient_to_degrees_from_north = 0.0
preprocessing_settings.filter_corner_frequencies_in_hz = (0.2, 30)
preprocessing_settings.ignore_dissimilar_time_step_warning = False

processing_settings = hvsrpy.settings.HvsrTraditionalProcessingSettings()
processing_settings.window_type_and_width = ("tukey", 0.01)
processing_settings.smoothing = dict(operator="konno_and_ohmachi",
                                     bandwidth=40,
                                     center_frequencies_in_hz=np.geomspace(0.2, 50, 256))
processing_settings.method_to_combine_horizontals = "geometric_mean"
processing_settings.handle_dissimilar_time_steps_by = "frequency_domain_resampling"

print("Preprocessing Summary")
print("-"*60)
preprocessing_settings.psummary()
print("Processing Summary")
print("-"*60)
processing_settings.psummary()

# Process in batches to cap peak memory usage.
# After each batch, only the compact frequency-domain HVSR object is kept;
# the raw seismic records are released before the next batch is loaded.
hvsr_batches = []
srecords_preprocessed_all = []  # kept for the rejection plot

for batch_start in range(0, len(all_fnames), BATCH_SIZE):
    batch = all_fnames[batch_start:batch_start + BATCH_SIZE]
    days_in_batch = common_days[batch_start:batch_start + BATCH_SIZE]
    print(f"\nBatch {batch_start // BATCH_SIZE + 1}: "
          f"days {days_in_batch[0]} → {days_in_batch[-1]} "
          f"({len(batch)} day(s))")

    print(f"  Reading:    {datetime.now().strftime('%H:%M:%S')}")
    srecords = hvsrpy.read(batch)

    print(f"  Preprocess: {datetime.now().strftime('%H:%M:%S')}")
    srecords_pre = hvsrpy.preprocess(srecords, preprocessing_settings)
    del srecords  # release raw records immediately

    print(f"  Processing: {datetime.now().strftime('%H:%M:%S')}")
    hvsr_batch = hvsrpy.process(srecords_pre, processing_settings)

    hvsr_batches.append(hvsr_batch)
    srecords_preprocessed_all.append(srecords_pre)
    # srecords_pre is kept only for plot_pre_and_post_rejection;
    # if memory is still tight you can del it here and skip that figure.

print(f"\nAll batches done: {datetime.now().strftime('%H:%M:%S')}")

# Merge HVSR windows from all batches.
# HvsrTraditional stores windows as .amplitude (shape: [n_windows, n_freqs])
# and .frequency. Re-construct via the constructor so all internal state
# (valid_window_boolean_mask, peak arrays, n_curves) is initialised cleanly.
combined_amplitude = np.concatenate([h.amplitude for h in hvsr_batches], axis=0)
hvsr = hvsrpy.HvsrTraditional(hvsr_batches[0].frequency, combined_amplitude,
                               meta=hvsr_batches[0].meta)

# Merge preprocessed records for the rejection plot.
# preprocess() returns a list of SeismicRecording3C — plain list concatenation.
srecords_preprocessed = []
for sp in srecords_preprocessed_all:
    srecords_preprocessed.extend(sp)

n = 2
search_range_in_hz = (None, None)
_ = hvsrpy.frequency_domain_window_rejection(hvsr, n=n, search_range_in_hz=search_range_in_hz)

mfig, axs = hvsrpy.plot_pre_and_post_rejection(srecords_preprocessed, hvsr)

save_figure = True
save_results = True
out_subdir = output_dir / "output_multiday"
fname_prefix = str(out_subdir / f"{sta}_WR")
os.makedirs(out_subdir, exist_ok=True)

if save_figure:
    fname = f"{fname_prefix}.png"
    mfig.savefig(fname)
    plt.close()
    print(f"Figure saved to {fname}")

if save_results:
    fname = f"{fname_prefix}.csv"
    hvsrpy.object_io.write_hvsr_object_to_file(hvsr, fname)
    print(f"Results saved to {fname}")

print("\nStatistical Summary:")
print("-"*20)
hvsrpy.summarize_hvsr_statistics(hvsr)

(sfig, ax) = hvsrpy.plot_single_panel_hvsr_curves(hvsr)
ax.get_legend().remove()
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

if save_figure:
    fname = f"{fname_prefix}_single_panel.png"
    sfig.savefig(fname)
    plt.close()
    print(f"Figure saved to {fname}")
