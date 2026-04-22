"""
Plot interpolated Fn (fundamental frequency) from HVSR on a map.

Reads stations_fn_summary.csv (station, Y, X, Z, fn_hz, fn_amplitude)
and produces a single-panel map with interpolated Fn values.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import os
from pathlib import Path

# ── I/O ──────────────────────────────────────────────────────────────────────
output_dir = Path("/srv/beegfs/scratch/users/c/cabrerap/hvsr_output")
summary_csv  = output_dir / "stations_fn_summary.csv"
output_fig   = output_dir / "map_fn_hz.png"

# ── Interpolation ─────────────────────────────────────────────────────────────
grid_size             = 500
interpolation_method  = "linear"   # 'nearest', 'linear', or 'cubic'
apply_gaussian        = True
gaussian_sigma        = 0.5

# ── Map / plot ────────────────────────────────────────────────────────────────
tile_zoom_level  = 13
colormap         = "RdYlBu_r"     # low Fn = blue (soft sediment), high = red (rock)
map_padding      = 0.005          # degrees added around station extent
marker_shape     = "^"
marker_color     = "darkgrey"
marker_edge_color = "black"
marker_size      = 40
alpha_markers    = 0.7
fig_size         = (8, 7)
base_font_size   = 11

plt.rcParams.update({"font.size": base_font_size})

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(summary_csv)
df = df.dropna(subset=["fn_hz", "X", "Y"])

print(f"Stations loaded : {len(df)}")
print(f"fn_hz range     : {df['fn_hz'].min():.3f} – {df['fn_hz'].max():.3f} Hz")

lats = df["Y"].values
lons = df["X"].values
fn   = df["fn_hz"].values

# ── Map extent ────────────────────────────────────────────────────────────────
min_lat = lats.min() - map_padding
max_lat = lats.max() + map_padding
min_lon = lons.min() - map_padding
max_lon = lons.max() + map_padding

# ── Interpolation grid ────────────────────────────────────────────────────────
lon_grid, lat_grid = np.mgrid[
    min_lon:max_lon:complex(0, grid_size),
    min_lat:max_lat:complex(0, grid_size)
]
fn_grid = griddata((lons, lats), fn, (lon_grid, lat_grid),
                   method=interpolation_method)
if apply_gaussian:
    fn_grid = gaussian_filter(fn_grid, sigma=gaussian_sigma)

vmin, vmax = fn.min(), fn.max()

# ── Figure ────────────────────────────────────────────────────────────────────
url  = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}.jpg"
tile = cimgt.GoogleTiles(url=url)

fig, ax = plt.subplots(1, 1,
                        subplot_kw={"projection": tile.crs},
                        figsize=fig_size)

ax.add_image(tile, tile_zoom_level)
ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.Geodetic())

im = ax.pcolormesh(lon_grid, lat_grid, fn_grid,
                   cmap=colormap, shading="auto",
                   transform=ccrs.PlateCarree(),
                   vmin=vmin, vmax=vmax,
                   alpha=0.75)

ax.scatter(lons, lats,
           marker=marker_shape,
           color=marker_color,
           edgecolor=marker_edge_color,
           s=marker_size,
           transform=ccrs.PlateCarree(),
           label="Stations",
           zorder=5,
           alpha=alpha_markers)

# Station name labels
for _, row in df.iterrows():
    ax.text(row["X"], row["Y"], row["station"],
            transform=ccrs.PlateCarree(),
            fontsize=5, ha="center", va="bottom",
            color="black", zorder=6)

cbar = fig.colorbar(im, ax=ax, orientation="vertical",
                    pad=0.05, fraction=0.046)
cbar.set_label("Fn (Hz)")

gl = ax.gridlines(draw_labels=True, alpha=0.3)
gl.top_labels   = False
gl.right_labels = False

ax.set_title("Fundamental frequency Fn from HVSR", fontsize=base_font_size + 1)

plt.tight_layout()
plt.savefig(output_fig, dpi=300, bbox_inches="tight")
plt.close()
print(f"Figure saved to {output_fig}")
