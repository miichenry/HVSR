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
output_dir = Path("/srv/beegfs/scratch/users/h/henrymi/vulcano/hvsr_output")

print(f"Station: {sta}")
print(f"Starting time: {datetime.now().strftime('%H:%M:%S')}")

sta_path = base_path / sta
if not sta_path.exists():
    raise FileNotFoundError(f"Station directory not found: {sta_path}")

# Collect files per component, keyed by (year, julianday)
def index_by_day(component_dir):
    day_map = {}
    for f in sorted(component_dir.glob("*.mseed")):
        # Filename: VN.A0000..DPE.D.2021.290.00.00.00.000.mseed
        m = re.search(r'\.(\d{4})\.(\d{3})\.', f.name)
        if m:
            key = (m.group(1), m.group(2))  # (year, julian_day)
            day_map[key] = str(f)
    return day_map

z_map = index_by_day(sta_path / "DPZ")
e_map = index_by_day(sta_path / "DPE")
n_map = index_by_day(sta_path / "DPN")

# Keep only days present in all three components
common_days = sorted(set(z_map) & set(e_map) & set(n_map))

if not common_days:
    raise RuntimeError(f"No days with complete Z/E/N data found for station {sta}")

print(f"Days with complete data: {len(common_days)}")
for year, jday in common_days:
    print(f"  {year}.{jday}")

fnames = [[z_map[d], e_map[d], n_map[d]] for d in common_days]

# Verify all files exist
for fname_set in fnames:
    for f in fname_set:
        if not pathlib.Path(f).exists():
            raise FileNotFoundError(f"File not found: {f}")
print(f"All files verified. ({len(fnames)} days)")

preprocessing_settings = hvsrpy.settings.HvsrPreProcessingSettings()
preprocessing_settings.detrend = "constant"
preprocessing_settings.window_length_in_seconds = 60
preprocessing_settings.orient_to_degrees_from_north = 0.0
preprocessing_settings.filter_corner_frequencies_in_hz = (0.2, 30)
preprocessing_settings.ignore_dissimilar_time_step_warning = False

print("Preprocessing Summary")
print("-"*60)
preprocessing_settings.psummary()

processing_settings = hvsrpy.settings.HvsrTraditionalProcessingSettings()
processing_settings.window_type_and_width = ("tukey", 0.01)
processing_settings.smoothing = dict(operator="konno_and_ohmachi",
                                     bandwidth=40,
                                     center_frequencies_in_hz=np.geomspace(0.2, 50, 256))
processing_settings.method_to_combine_horizontals = "geometric_mean"
processing_settings.handle_dissimilar_time_steps_by = "frequency_domain_resampling"

print("Processing Summary")
print("-"*60)
processing_settings.psummary()

print(f"Start reading: {datetime.now().strftime('%H:%M:%S')}")
srecords = hvsrpy.read(fnames)
srecords_preprocessed = hvsrpy.preprocess(srecords, preprocessing_settings)
hvsr = hvsrpy.process(srecords_preprocessed, processing_settings)
print(f"End processing: {datetime.now().strftime('%H:%M:%S')}")

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
