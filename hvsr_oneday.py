"""Process HVSR for one station, one day at a time.

Each day is processed independently and saved as a per-day CSV.
No srecords are accumulated between days, so peak memory is bounded
to a single day's worth of data.

Usage:
    python hvsr_oneday.py <station_arg>
    station_arg: station name or CSV row like "A6575,lat,lon,elev"
"""

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
base_path = Path("/srv/beegfs/scratch/users/h/henrymi/vulcano/miniseed_all")
output_dir = Path("/srv/beegfs/scratch/users/h/henrymi/vulcano/output_test")
perday_dir = output_dir / "perday" / sta

print(f"Station: {sta}")
print(f"Starting time: {datetime.now().strftime('%H:%M:%S')}")

sta_path = base_path
if not sta_path.exists():
    raise FileNotFoundError(f"Station directory not found: {sta_path}")


def index_by_day(component_dir):
    day_map = {}
    for f in sorted(component_dir.glob("*.mseed")):
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
os.makedirs(perday_dir, exist_ok=True)

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

for year, jday in common_days:
    out_csv = perday_dir / f"{sta}_{year}_{jday}.csv"
    if out_csv.exists():
        print(f"  {year}.{jday} — already done, skipping")
        continue

    print(f"\n  {year}.{jday}  [{datetime.now().strftime('%H:%M:%S')}]")
    try:
        fnames = [[z_map[(year, jday)], e_map[(year, jday)], n_map[(year, jday)]]]

        srecords = hvsrpy.read(fnames)
        srecords_pre = hvsrpy.preprocess(srecords, preprocessing_settings)
        del srecords

        hvsr = hvsrpy.process(srecords_pre, processing_settings)
        del srecords_pre

        hvsrpy.object_io.write_hvsr_object_to_file(hvsr, str(out_csv))
        del hvsr
        print(f"    Saved {out_csv.name}")

    except Exception as e:
        print(f"    ERROR on {year}.{jday}: {e}")

print(f"\nAll days done: {datetime.now().strftime('%H:%M:%S')}")
