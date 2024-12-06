import pathlib
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import hvsrpy
from hvsrpy import sesame

import sys
from glob import glob

plt.style.use(hvsrpy.HVSRPY_MPL_STYLE)

sta = sys.argv[1]
path = os.path.join('/srv', 'beegfs', 'scratch', 'shares', 'cdff', 'hautesorne')

print(f"station:{sta}")
print(f"length:{len(sta)}")

current_time = datetime.now().strftime('%H:%M:%S')
print(f"Starting time: {current_time}")

Z = glob(os.path.join(path, 'DATA_MSEED', '01032024', f'*{sta}*Z*'))
E = glob(os.path.join(path, 'DATA_MSEED', '01032024', f'*{sta}*E*'))
N = glob(os.path.join(path, 'DATA_MSEED', '01032024', f'*{sta}*N*'))
print(Z)
def check_is_string(variable, name):
    if not isinstance(variable, str):
        raise TypeError(f"{name} must be a string, but received type: {type(variable)}")

# These checks will stop the program if any fail
check_is_string(Z[0], "Z")
check_is_string(E[0], "E")
check_is_string(N[0], "N")

# Code below will not execute if there's an error above
print("All checks passed. Z, E, and N are valid lists of strings.")

fnames = [[Z[0],E[0],N[0]]]
print(fnames)

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
preprocessing_settings.window_length_in_seconds = 60 #60s for coherency next run 
preprocessing_settings.orient_to_degrees_from_north = 0.0
preprocessing_settings.filter_corner_frequencies_in_hz = (0.1, 30)
preprocessing_settings.ignore_dissimilar_time_step_warning = False

print("Preprocessing Summary")
print("-"*60)
preprocessing_settings.psummary()


processing_settings = hvsrpy.settings.HvsrDiffuseFieldProcessingSettings()
processing_settings.window_type_and_width = ("tukey", 0.1)
processing_settings.smoothing=dict(operator="log_rectangular",
                                   bandwidth=0.1,
                                   center_frequencies_in_hz=np.geomspace(0.2, 50, 256))

print("Processing Summary")
print("-"*60)
processing_settings.psummary()

print(f"Start reading: {datetime.now().strftime('%H:%M:%S')}")
srecords = hvsrpy.read(fnames)
srecords = hvsrpy.preprocess(srecords, preprocessing_settings)
hvsr = hvsrpy.process(srecords, processing_settings)
print(f"End reading: {datetime.now().strftime('%H:%M:%S')}")

hvsrpy.summarize_hvsr_statistics(hvsr)
(fig, ax) = hvsrpy.plot_single_panel_hvsr_curves(hvsr,)
plt.show()

save_figure = True
save_results = True
fname_prefix = f"output_DFA_2/{sta}_DFA"

if save_figure:
    fname = f"{fname_prefix}.png"
    fig.savefig(fname)
    plt.close()
    print(f"Figure saved successfully to {fname}!")

if save_results:
    fname = f"{fname_prefix}.csv"
    hvsrpy.object_io.write_hvsr_object_to_file(hvsr, fname)
    print(f"Results saved successfully to {fname}!")
