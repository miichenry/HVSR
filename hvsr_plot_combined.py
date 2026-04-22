"""Generate single-panel HVSR plot from a combined CSV produced by hvsr_combine.py.

Usage:
    python hvsr_plot_combined.py <station_arg>
    station_arg: station name or CSV row like "A6575,lat,lon,elev"
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hvsrpy
from pathlib import Path
import sys

plt.style.use(hvsrpy.HVSRPY_MPL_STYLE)

#sta = sys.argv[1].split(',')[0]
#output_dir = Path("/srv/beegfs/scratch/users/h/henrymi/vulcano/output_test")
#combined_dir = output_dir / "combined"

#csv_path = combined_dir / f"L6767_WR.csv"
csv_path = "../L6767_WR.csv"
#if not csv_path.exists():
#    raise FileNotFoundError(f"Combined CSV not found: {csv_path}")

hvsr = hvsrpy.object_io.read_hvsr_object_from_file(str(csv_path))

(sfig, ax) = hvsrpy.plot_single_panel_hvsr_curves(hvsr)
ax.get_legend().remove()
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

#fname = str(combined_dir / f"{sta}_WR_single_panel.png")
fname = str(f"../L6767_WR_single_panel.png")
sfig.savefig(fname)
plt.close()
print(f"Figure saved to {fname}")
