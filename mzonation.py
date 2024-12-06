import time
import random
from obspy import read

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
#matplotlib.use('PDF')

import hvsrpy
from hvsrpy import utils
from obspy import Stream, read, Trace
from glob import glob
import os
from tqdm import tqdm
import multiprocessing
from scipy.interpolate import griddata

#---------------------------------------------------#PARAMETERS#--------------------------------------------------------
path = '/srv/beegfs/scratch/shares/cdff/hautesorne'
stations_file = path + '/HS_ANT_stations_locations_from_log_median_wName.csv'
out_folder = '/home/users/h/henrymi2/HVSR/hvfiles/'

# Window length in seconds. In general low frequency peaks require longer window lengths.
# See the SESAME guidelines for specific window length recommendations.
windowlength = 50

# Boolean to control whether Butterworth filter is applied.
# Geopsy does not apply a bandpass filter.
filter_bool = False
# Low-cut frequency for bandpass filter.
filter_flow = 0.1
# High-cut frequency for bandpass filter.
filter_fhigh = 30
# Filter order.
filter_order = 5

# Width of cosine taper {0. - 1.}. Geopsy default of 0.05 is equal to 0.1 -> 0.1 is recommended
width = 0.1

# Konno and Ohmachi smoothing constant. 40 is recommended.
bandwidth = 40

# Minimum frequency after resampling
resample_fmin = 0.1
# Maximum frequency after resampling
resample_fmax = 40
# Number of frequencies after resampling
resample_fnum = 100
# Type of resampling {'log', 'linear'}
resample_type = 'log'

# Upper and lower frequency limits to restrict peak selection. To use the entire range use `None`.
peak_f_lower = 0.1
peak_f_upper = 10.0

# Method for combining horizontal components {"squared-average", "geometric-mean", "single-azimuth"}.
# Geopsy's default is "squared-average" -> "geometric-mean" is recommended.
method = "geometric-mean"
# If method="single-azimuth", set azimuth in degree clock-wise from north. If method!="single-azimuth", value is ignored.
azimuth = 0

# Boolean to control whether frequency domain rejection proposed by Cox et al. (2020) is applied.
# Geopsy does not offer this functionality.
rejection_bool = True
# Number of standard deviations to consider during rejection. Smaller values will reject more windows -> 2 is recommended.
n = 2
# Maximum number of iterations to perform for rejection -> 50 is recommended.
max_iterations = 50

# Distribution of f0 {"lognormal", "normal"}. Geopsy default "normal" -> "lognormal" is recommended.
distribution_f0 = "lognormal"
# Distribution of mean curve {"lognormal", "normal"}. Geopsy default "lognormal" -> "lognormal" is recommended.
distribution_mc = "lognormal"

# Manually set the ylimits of the HVSR figures. Default is None so limits will be set automatically.
ymin, ymax = 0, 10

calculate_HVSR = False
#-----------------------------------------------------------------------------------------------------------------------

stations = np.loadtxt(stations_file, usecols=2, skiprows=1, delimiter=',', dtype=str)
lats = np.loadtxt(stations_file, usecols=5, skiprows=1, delimiter=',')
lons = np.loadtxt(stations_file, usecols=4, skiprows=1, delimiter=',')

def process_station(sta):
    if os.path.exists(out_folder+'/{:s}_geopsy.hv'.format(sta)) == True:
        print('The output file for the {:s} station already exists, or there are no data files for this seismic station.'.format(sta))
    else:

        stZ = Stream()
        stN = Stream()
        stE = Stream()

        for i in sorted(glob(path + '/DATA_MSEED/' + '*/' + '*{:s}*Z*'.format(sta))):
            st = read(i)
            fr = st[0].stats.sampling_rate
            stZ += st
        for i in sorted(glob(path + '/DATA_MSEED/' + '*/' + '*{:s}*E*'.format(sta))):
            st = read(i)
            stE += st
        for i in sorted(glob(path + '/DATA_MSEED/' + '*/' + '*{:s}*N*'.format(sta))):
            st = read(i)
            stN += st

       # Check if the streams are not empty
        if len(stZ) == 0 or len(stE) == 0 or len(stN) == 0:
            print(f"No data found for station {sta}. Skipping.")
            return

        stZ.merge()
        stE.merge()
        stN.merge()

        stZ.trim(starttime=stZ[0].stats.starttime, endtime=stZ[0].stats.endtime)
        stE.trim(starttime=stZ[0].stats.starttime, endtime=stZ[0].stats.endtime)
        stN.trim(starttime=stZ[0].stats.starttime, endtime=stZ[0].stats.endtime)

        lst = [len(stZ[0].data), len(stE[0].data), len(stN[0].data)]
        chk = all(element == lst[0] for element in lst)

        if not chk:
            axmin = np.argmin([len(stZ[0].data), len(stE[0].data), len(stN[0].data)])
            
            # Check if the data is a masked array before calling the 'filled' method
            traceZ_data = stZ[0].data.filled(np.nan) if np.ma.isMaskedArray(stZ[0].data) else stZ[0].data[:lst[axmin]]
            traceE_data = stE[0].data.filled(np.nan) if np.ma.isMaskedArray(stE[0].data) else stE[0].data[:lst[axmin]]
            traceN_data = stN[0].data.filled(np.nan) if np.ma.isMaskedArray(stN[0].data) else stN[0].data[:lst[axmin]]

            traceZ = Trace(traceZ_data, header=stZ[0].stats)
            traceE = Trace(traceE_data, header=stE[0].stats)
            traceN = Trace(traceN_data, header=stN[0].stats)

            # Convert masked arrays to normal arrays before writing
            St = Stream(traces=[traceZ, traceN, traceE])
            St.write('stream_file.miniseed', format='MSEED')
        else:
            # Check if the data is a masked array before calling the 'filled' method
            traceZ_data = stZ[0].data.filled(np.nan) if np.ma.isMaskedArray(stZ[0].data) else stZ[0].data
            traceE_data = stE[0].data.filled(np.nan) if np.ma.isMaskedArray(stE[0].data) else stE[0].data
            traceN_data = stN[0].data.filled(np.nan) if np.ma.isMaskedArray(stN[0].data) else stN[0].data

            traceZ = Trace(traceZ_data, header=stZ[0].stats)
            traceE = Trace(traceE_data, header=stE[0].stats)
            traceN = Trace(traceN_data, header=stN[0].stats)

            St = Stream(traces=[traceZ, traceN, traceE])
            St.write('stream_file.miniseed', format='MSEED')

        # Check if the file exists and is not empty before trying to read it
        file_name = os.getcwd() + '/stream_file.miniseed'
        if os.path.exists(file_name) and os.path.getsize(file_name) > 128:
            sensor = hvsrpy.Sensor3c.from_mseed(file_name)
        else:
            print(f"File {file_name} is empty or too small to be valid. Skipping station {sta}.")
            return
        
        fig = plt.figure(figsize=(6,6), dpi=150)
        gs = fig.add_gridspec(nrows=6,ncols=6)

        ax0 = fig.add_subplot(gs[0:2, 0:3])
        ax1 = fig.add_subplot(gs[2:4, 0:3])
        ax2 = fig.add_subplot(gs[4:6, 0:3])

        if rejection_bool:
            ax3 = fig.add_subplot(gs[0:3, 3:6])
            ax4 = fig.add_subplot(gs[3:6, 3:6])
        else:
            ax3 = fig.add_subplot(gs[0:3, 3:6])
            ax4 = False

        start = time.time()
        #sensor = hvsrpy.Sensor3c.from_mseed(file_name)
        bp_filter = {"flag":filter_bool, "flow":filter_flow, "fhigh":filter_fhigh, "order":filter_order}
        resampling = {"minf":resample_fmin, "maxf":resample_fmax, "nf":resample_fnum, "res_type":resample_type}
        hv = sensor.hv(windowlength, bp_filter, width, bandwidth, resampling, method, f_low=peak_f_lower, f_high=peak_f_upper, azimuth=azimuth)
        end = time.time()
        print(f"Elapsed Time: {str(end-start)[0:4]} seconds")

        individual_width = 0.3
        median_width = 1.3
        for ax, title in zip([ax3, ax4], ["Before Rejection", "After Rejection"]):
            # Rejected Windows
            if title=="After Rejection":
                if len(hv.rejected_window_indices):
                    label = "Rejected"
                    for amp in hv.amp[hv.rejected_window_indices]:
                        ax.plot(hv.frq, amp, color='#00ffff', linewidth=individual_width, zorder=2, label=label)
                        label=None

            # Accepted Windows
            label="Accepted"
            for amp in hv.amp[hv.valid_window_indices]:
                ax.plot(hv.frq, amp, color='#888888', linewidth=individual_width,
                        label = label if title=="Before Rejection" else "")
                label=None

            # Window Peaks
            ax.plot(hv.peak_frq, hv.peak_amp, linestyle="", zorder=2,
                    marker='o', markersize=2.5, markerfacecolor="#ffffff", markeredgewidth=0.5, markeredgecolor='k',
                    label="" if title=="Before Rejection" and rejection_bool else r"$f_{0,i}$")

            # Peak Mean Curve
            ax.plot(hv.mc_peak_frq(distribution_mc), hv.mc_peak_amp(distribution_mc), linestyle="", zorder=4,
                    marker='D', markersize=4, markerfacecolor='#66ff33', markeredgewidth=1, markeredgecolor='k',
                    label = "" if title=="Before Rejection" and rejection_bool else r"$f_{0,mc}$")

            # Mean Curve
            label = r"$LM_{curve}$" if distribution_mc=="lognormal" else "Mean"
            ax.plot(hv.frq, hv.mean_curve(distribution_mc), color='k', linewidth=median_width,
                    label="" if title=="Before Rejection" and rejection_bool else label)

            # Mean +/- Curve
            label = r"$LM_{curve}$"+" ± 1 STD" if distribution_mc=="lognormal" else "Mean ± 1 STD"
            ax.plot(hv.frq, hv.nstd_curve(-1, distribution_mc),
                    color='k', linestyle='--', linewidth=median_width, zorder=3,
                    label = "" if title=="Before Rejection" and rejection_bool else label)
            ax.plot(hv.frq, hv.nstd_curve(+1, distribution_mc),
                    color='k', linestyle='--', linewidth=median_width, zorder=3)

            # f0 +/- STD
            if ymin is not None and ymax is not None:
                ax.set_ylim((ymin, ymax))
            label = r"$LM_{f0}$"+" ± 1 STD" if distribution_f0=="lognormal" else "Mean f0 ± 1 STD"
            _ymin, _ymax = ax.get_ylim()
            ax.plot([hv.mean_f0_frq(distribution_f0)]*2, [_ymin, _ymax], linestyle="-.", color="#000000")
            ax.fill([hv.nstd_f0_frq(-1, distribution_f0)]*2 + [hv.nstd_f0_frq(+1, distribution_f0)]*2, [_ymin, _ymax, _ymax, _ymin],
                    color = "#ff8080",
                    label="" if title=="Before Rejection" and rejection_bool else label)
            ax.set_ylim((_ymin, _ymax))

            ax.set_xscale('log')
            ax.set_xlabel("Frequency (Hz)")
            ax.set_ylabel("HVSR Amplitude")
            if rejection_bool:
                if title=="Before Rejection":
                    print("\nStatistics before rejection:")
                    #hv.print_stats(distribution_f0)
                    c_iter = hv.reject_windows(n, max_iterations=max_iterations,
                                               distribution_f0=distribution_f0, distribution_mc=distribution_mc)
                elif title=="After Rejection":
                    fig.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.51, 0), columnspacing=2)

                    print("\nAnalysis summary:")
                    #display(pd.DataFrame(columns=[""], index=["Window length", "No. of windows", "Number of iterations to convergence", "No. of rejected windows"],
                    #        data=[f"{windowlength}s", str(sensor.ns.nseries), f"{c_iter} of {max_iterations} allowed", str(sum(hv.rejected_window_indices))]))
                    #print("\nStatistics after rejection:")
                    #hv.print_stats(distribution_f0)
            else:
                #display(pd.DataFrame(columns=[""], index=["Window length", "No. of windows"],
                #                 data=[f"{windowlength}s", str(sensor.ns.nseries)]))
                #hv.print_stats(distribution_f0)
                fig.legend(loc="upper center", bbox_to_anchor=(0.77, 0.4))
                break
            ax.set_title(title)

        norm_factor = sensor.normalization_factor
        for ax, timerecord, name in zip([ax0,ax1,ax2], [sensor.ns, sensor.ew, sensor.vt], ["NS", "EW", "VT"]):
            ctime = timerecord.time
            amp = timerecord.amp/norm_factor
            ax.plot(ctime.T, amp.T, linewidth=0.2, color='#888888')
            ax.set_title(f"Time Records ({name})")
            ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            ax.set_xlim(0, windowlength*timerecord.nseries)
            ax.set_ylim(-1, 1)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Normalized Amplitude')
            ax.plot(ctime[hv.rejected_window_indices].T, amp[hv.rejected_window_indices].T, linewidth=0.2, color="cyan")

        if rejection_bool:
            axs = [ax0, ax3, ax1, ax4, ax2]
        else:
            axs = [ax0, ax3, ax1, ax2]

        for ax, letter in zip(axs, list("abcde")):
            ax.text(0.97, 0.97, f"({letter})", ha="right", va="top", transform=ax.transAxes, fontsize=12)
            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)


            fig.tight_layout(h_pad=1, w_pad=2, rect=(0,0.08,1,1))

            figure_name_out = out_folder + "{:s}_hvsr_figure.png".format(sta)

            fig.savefig(figure_name_out, dpi=300, bbox_inches='tight')
            plt.close()
            print("Figure saved successfully!")

            file_name_out = "{:s}_hvsrpy.hv".format(sta)

            hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="hvsrpy")
            print("Results saved successfully!")

            file_name_out = "{:s}_geopsy.hv".format(sta)

            hv.to_file(file_name_out, distribution_f0, distribution_mc, data_format="geopsy")
            print("Results saved successfully!")


if __name__ == '__main__':
    stations = np.loadtxt(stations_file, usecols=2, skiprows=1, delimiter=',', dtype=str)
    lats = np.loadtxt(stations_file, usecols=4, skiprows=1, delimiter=',')
    lons = np.loadtxt(stations_file, usecols=5, skiprows=1, delimiter=',')

    num_processes = multiprocessing.cpu_count()  # Number of CPU cores
    num_processes = 5  # Number of CPU cores
    pool = multiprocessing.Pool(processes=num_processes)

    for _ in tqdm(pool.imap_unordered(process_station, stations), total=len(stations)):
        pass

    pool.close()
    pool.join()



grid_x, grid_y = np.mgrid[np.min(lons):np.max(lons):100j, np.min(lats):np.max(lats):100j]
points = []
f0 = []
A0 = []
m = []

for ix, sta in enumerate(stations):
    print(sta)
    out = open(out_folder +'{:s}_hvsrpy.hv'.format(sta), 'r').readlines()
    points.append((lons[ix], lats[ix]))
    for r in out:
        if r.split(',')[0] == "# Median Curve Peak Frequency (Hz) [f0mc]":
            f0.append(r.split(',')[1])

        if r.split(',')[0] == "# Median Curve Peak Amplitude ()":
            A0.append(r.split(',')[1])

gridf0 = griddata(points, f0, (grid_x, grid_y), method='cubic')
gridA0 = griddata(points, A0, (grid_x, grid_y), method='cubic')

M = 96*pow(gridf0, -1.388)

fig = plt.figure(figsize=(12, 5))
gs = GridSpec(1, 4)

# Creating subplots within the GridSpec layout
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[3])

kg = (gridA0**2)/gridf0


np.savez(path + '/gridf0.npz', matrix=gridf0, x=grid_x, y=grid_y)
np.savez(path + '/gridA0.npz', matrix=gridA0, x=grid_x, y=grid_y)
np.savez(path + '/gridM.npz', matrix=gridM, x=grid_x, y=grid_y)
np.savez(path + '/gridkg.npz', matrix=kg, x=grid_x, y=grid_y)

im = ax1.pcolormesh(grid_x, grid_y, gridf0)
ax1.contourf(grid_x, grid_y, gridf0)
ax1.set_xlabel('Longitude (°)')
ax1.set_ylabel('Latitude (°)')
cbar = plt.colorbar(im, ax=ax1)
cbar.set_label('frequency (Hz)')

im = ax2.pcolormesh(grid_x, grid_y, gridA0)
ax2.contourf(grid_x, grid_y, gridA0)
ax2.set_xlabel('Longitude (°)')
ax2.set_yticklabels([])
cbar = plt.colorbar(im, ax=ax2)
cbar.set_label('Amplitude')

im = ax3.pcolormesh(grid_x, grid_y, kg)
ax3.contourf(grid_x, grid_y, kg)
ax3.set_xlabel('Longitude (°)')
ax3.set_yticklabels([])
cbar = plt.colorbar(im, ax=ax3)
cbar.set_label('Vulnerability index')

im = ax4.pcolormesh(grid_x, grid_y, M)
ax4.contourf(grid_x, grid_y, M)
ax4.set_xlabel('Longitude (°)')
ax4.set_yticklabels([])
cbar = plt.colorbar(im, ax=ax4)
cbar.set_label('Thickness (m)')

plt.show()

