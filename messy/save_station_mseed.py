import numpy as np
from obspy import Stream, read, Trace
from glob import glob
import os
from tqdm import tqdm
import multiprocessing
import pandas as pd
import sys
from datetime import datetime

#---------------------------------------------------# PARAMETERS #-----------------------------------------------------
num = sys.argv[1]
path = os.path.join('/srv', 'beegfs', 'scratch', 'shares', 'cdff', 'hautesorne')
stations_file = os.path.join('/home', 'users', 'h', 'henrymi2', 'HVSR', 'csv', f'fitlered_HS_ANT_stations_part_{num}.csv')
out_folder = os.path.join('/home', 'users', 'h', 'henrymi2', 'scratch', 'mseed_02')
print(num,stations_file)
#station = pd.read_csv(stations_file, header=0).iloc[0, 2]
#station = sys.argv[1]

stations = np.loadtxt(stations_file, usecols=2, skiprows=1, delimiter=',', dtype=str)
print(f"station:{stations}")
print(f"length:{len(stations)}")
current_time = datetime.now().strftime('%H:%M:%S')
print(f"Starting time: {current_time}")
def process_station(sta):
    stZ = Stream()
    stN = Stream()
    stE = Stream()

    for i in sorted(glob(path + '/DATA_MSEED/' + '*22024/' + f'*{sta}*Z*')):
        st = read(i)
        fr = st[0].stats.sampling_rate
        stZ += st
    for i in sorted(glob(path + '/DATA_MSEED/' + '*22024/' + f'*{sta}*E*')):
        st = read(i)
        stE += st
    for i in sorted(glob(path + '/DATA_MSEED/' + '*22024/' + f'*{sta}*N*')):
        st = read(i)
        stN += st

    # Check if the streams are not empty
    if len(stZ) == 0 or len(stE) == 0 or len(stN) == 0:
        print(f"No data found for station {sta}. Skipping.")
        return

    stZ.merge()
    stE.merge()
    stN.merge()

    # Display stats for Z, E, and N componentsi
    print(f"{datetime.now().strftime('%H:%M:%S')}")
    print(f"Station {sta} - Stream Z stats: {stZ}")
    print(f"Station {sta} - Stream E stats: {stE}")
    print(f"Station {sta} - Stream N stats: {stN}")

    stZ.trim(starttime=stZ[0].stats.starttime, endtime=stZ[0].stats.endtime)
    stE.trim(starttime=stZ[0].stats.starttime, endtime=stZ[0].stats.endtime)
    stN.trim(starttime=stZ[0].stats.starttime, endtime=stZ[0].stats.endtime)
    print(f"Trimming pass: {datetime.now().strftime('%H:%M:%S')}")

    # Check if the data is a masked array and fill with NaN if necessary
    stZ[0].data = stZ[0].data.filled(np.nan) if np.ma.isMaskedArray(stZ[0].data) else stZ[0].data
    stE[0].data = stE[0].data.filled(np.nan) if np.ma.isMaskedArray(stE[0].data) else stE[0].data
    stN[0].data = stN[0].data.filled(np.nan) if np.ma.isMaskedArray(stN[0].data) else stN[0].data
    print(f"Masked array pass: {datetime.now().strftime('%H:%M:%S')}")
    lst = [len(stZ[0].data), len(stE[0].data), len(stN[0].data)]
    chk = all(element == lst[0] for element in lst)
    print(f"lst+chk: {datetime.now().strftime('%H:%M:%S')}")
    
    try:
    # Attempt to write the file directly if memory allows
        if chk:
            St = Stream(traces=[stZ[0], stN[0], stE[0]])
            St.write(f'{out_folder}/{sta}.mseed', format='MSEED')
        else:
            axmin = np.argmin(lst)
            traceZ = Trace(stZ[0].data[:lst[axmin]])
            traceZ.stats.component = 'Z'
            traceZ.stats.starttime = stZ[0].stats.starttime
            traceZ.stats.sampling_rate = fr

            traceE = Trace(stE[0].data[:lst[axmin]])
            traceE.stats.component = 'E'
            traceE.stats.starttime = stZ[0].stats.starttime
            traceE.stats.sampling_rate = fr

            traceN = Trace(stN[0].data[:lst[axmin]])
            traceN.stats.component = 'N'
            traceN.stats.starttime = stZ[0].stats.starttime
            traceN.stats.sampling_rate = fr

            St = Stream(traces=[traceZ, traceN, traceE])
            St.write(f'{out_folder}/{sta}.mseed', format='MSEED')

        print(f'MSEED file saved successfully for station {sta}!')

    except Exception as e:
        print(f"Error writing MSEED file for station {sta}: {e}")
        sys.exit(1)  # Exit if a fatal error occurs during file writing

if __name__ == '__main__':
   #print(f"Processing station {station}")
   #process_station(station)
   num_processes = multiprocessing.cpu_count()
   pool = multiprocessing.Pool(processes=num_processes)

   for _ in tqdm(pool.imap_unordered(process_station, stations), total=len(stations)):
       pass

   pool.close()
   pool.join()
