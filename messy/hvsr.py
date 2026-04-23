import pathlib
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import hvsrpy
from hvsrpy import sesame

plt.style.use(hvsrpy.HVSRPY_MPL_STYLE)

# Set the directory path
directory_path = '/home/users/h/henrymi2/scratch/mseed_02'

# Initialize the filenames list
fnames = [["/srv/beegfs/scratch/shares/cdff/hautesorne/DATA_MSEED/01032024/453005265.2024.03.01.00.00.00.000.N.miniseed","/srv/beegfs/scratch/shares/cdff/hautesorne/DATA_MSEED/01032024/453005265.2024.03.01.00.00.00.000.E.miniseed","/srv/beegfs/scratch/shares/cdff/hautesorne/DATA_MSEED/01032024/453005265.2024.03.01.00.00.00.000.Z.miniseed"]]

# Loop through the files in the directory and add them to fnames directly
#for f in os.listdir(directory_path):
 #   file_path = os.path.join(directory_path, f)
  #  if os.path.isfile(file_path):
   #     fnames.append(file_path)

# Verify the structure of fnames and ensure each entry is a valid path
print(f"Number of recordings: {len(fnames)}")
for i, fname_set in enumerate(fnames):
    if not isinstance(fname_set, list):
        print(f"Error: Entry at index {i} in fnames is not a list.")
        continue
    
    for file in fname_set:
        # Check if file is indeed a path string
        if not isinstance(file, str):
            print(f"Error: Non-string file path detected: {file}")
            continue

        # Check if the file exists
        if not pathlib.Path(file).exists():
            raise FileNotFoundError(f"File '{file}' not found; check spelling and path.")
print(f"All files exist. (n={len(fnames)})")

preprocessing_settings = hvsrpy.settings.HvsrPreProcessingSettings()
preprocessing_settings.detrend = "linear"
preprocessing_settings.window_length_in_seconds = 50
preprocessing_settings.orient_to_degrees_from_north = 0.0
preprocessing_settings.filter_corner_frequencies_in_hz = (0.1, 30)
preprocessing_settings.ignore_dissimilar_time_step_warning = False

print("Preprocessing Summary")
print("-"*60)
preprocessing_settings.psummary()


processing_settings = hvsrpy.settings.HvsrTraditionalProcessingSettings()
processing_settings.window_type_and_width = ("tukey", 0.2)
processing_settings.smoothing=dict(operator="konno_and_ohmachi",
                                   bandwidth=40,
                                   center_frequencies_in_hz=np.geomspace(0.1, 40, 100))
processing_settings.method_to_combine_horizontals = "geometric_mean"
processing_settings.handle_dissimilar_time_steps_by = "frequency_domain_resampling"

print("Processing Summary")
print("-"*60)
processing_settings.psummary()

print(f"Start reading: {datetime.now().strftime('%H:%M:%S')}")
srecords = hvsrpy.read(fnames)
srecords = hvsrpy.preprocess(srecords, preprocessing_settings)
hvsr = hvsrpy.process(srecords, processing_settings)
print(f"End reading: {datetime.now().strftime('%H:%M:%S')}")
search_range_in_hz = (None, None)
verbose = 1
hvsr.update_peaks_bounded(search_range_in_hz=search_range_in_hz)

print("\nSESAME (2004) Clarity and Reliability Criteria:")
print("-"*47)
hvsrpy.sesame.reliability(
    windowlength=preprocessing_settings.window_length_in_seconds,
    passing_window_count=np.sum(hvsr.valid_window_boolean_mask),
    frequency=hvsr.frequency,
    mean_curve=hvsr.mean_curve(distribution="lognormal"),
    std_curve=hvsr.std_curve(distribution="lognormal"),
    search_range_in_hz=search_range_in_hz,
    verbose=verbose,
)
hvsrpy.sesame.clarity(
    frequency=hvsr.frequency,
    mean_curve=hvsr.mean_curve(distribution="lognormal"),
    std_curve=hvsr.std_curve(distribution="lognormal"),
    fn_std=hvsr.std_fn_frequency(distribution="normal"),
    search_range_in_hz=search_range_in_hz,
    verbose=verbose,
)

print("\nStatistical Summary:")
print("-"*20)
hvsrpy.summarize_hvsr_statistics(hvsr)

fig, ax = hvsrpy.plot_single_panel_hvsr_curves(hvsr)
ax.get_legend().remove()
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.show()

save_figure = True
save_results = True
fname_prefix = "hvsr_sesame"

if save_figure:
    fname = f"{fname_prefix}.png"
    fig.savefig(fname)
    plt.close()
    print(f"Figure saved successfully to {fname}!")

if save_results:
    fname = f"{fname_prefix}.csv"
    hvsrpy.object_io.write_hvsr_object_to_file(hvsr, fname)
    print(f"Results saved successfully to {fname}!")
