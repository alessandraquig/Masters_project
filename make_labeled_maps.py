import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from itertools import product

import cartopy.crs as ccrs
from cartopy import feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes

import parameters as par

if __name__ == '__main__':
    plot_dir = "Maps_output/labelled_maps"

    this_rc_params = {
        "text.usetex": True,
        "font.family": "roman"
    }
    plt.rcParams.update(this_rc_params)
    plt.rcParams['figure.dpi'] = 400

    os.makedirs(plot_dir, exist_ok=True)

    # %%
    fig: plt.Figure = plt.figure(figsize=(5, 5))
    ax: GeoAxes = fig.add_subplot(1, 1, 1, projection=par.m)
    ax.set_extent(par.geo_bounds, ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)  # or LAND or OCEAN

    land_label_colour = "k"
    ocean_label_colour = "k"

    if par.HEMI == "north":
        import grid_set as gs

        GPathfinder = gs.grid_set(par.m)
        GPathfinder.load_grid(par.ID_grid)

        bg_xs = 125, 125, 180, 180
        bg_ys = 217, 257, 257, 217
        bg_points = np.array([[GPathfinder.xpts[x, y], GPathfinder.ypts[x, y]] for x, y in zip(bg_xs, bg_ys)])

        beaufort = mpatches.Polygon(bg_points, closed=True, edgecolor="tab:red", facecolor="none")
        ax.add_patch(beaufort)

        place_labels = [
            {
                "label": "USA\n(Alaska)",
                "x": 0.272,
                "y": 0.068,
                "c": land_label_colour
            },
            {
                "label": "Canada",
                "x": 0.075,
                "y": 0.226,
                "c": land_label_colour
            },
            {
                "label": "Greenland",
                "x": 0.323,
                "y": 0.720,
                "c": land_label_colour
            },
            {
                "label": "Russia",
                "x": 0.909,
                "y": 0.349,
                "c": land_label_colour
            },
            {
                "label": "Chukchi\nSea",
                "x": 0.45,
                "y": 0.125,
                "c": ocean_label_colour,
                # "fontsize": 10
            },
            {
                "label": "Beaufort\nSea",
                "x": 0.35,
                "y": 0.237,
                "c": ocean_label_colour
            },
            {
                "label": "Baffin\nBay",
                "x": 0.2,
                "y": 0.639,
                "c": ocean_label_colour
            },
            {
                "label": "Fram\nStrait",
                "x": 0.489,
                "y": 0.75,
                "c": ocean_label_colour
            },
            {
                "label": "Barents\nSea",
                "x": 0.7,
                "y": 0.75,
                "c": ocean_label_colour
            },
            {
                "label": "Laptev\nSea",
                "x": 0.721,
                "y": 0.35,
                "c": ocean_label_colour
            },
            {
                "label": "East\nSiberian\nSea",
                "x": 0.591,
                "y": 0.2,
                "c": ocean_label_colour
            },
            {
                "label": "Arctic\nOcean",
                "x": 0.5,
                "y": 0.5,
                "c": ocean_label_colour,
                "fontsize": 16
            }
        ]
    else:
        place_labels = [
            {
                "label": "Antarctica",
                "x": 0.5,
                "y": 0.5,
                "c": land_label_colour,
                "fontsize": 16
            },
            {
                "label": "Weddell Sea",
                "x": 0.3,
                "y": 0.317,
                "c": ocean_label_colour,
            },
            {
                "label": "Riiser-Larsen\nSea",
                "x": 0.642,
                "y": 0.22,
                "c": ocean_label_colour,
            },
            {
                "label": "Indian\nOcean",
                "x": 0.9,
                "y": 0.5,
                "c": ocean_label_colour,
            },
            {
                "label": "Ross Sea",
                "x": 0.441,
                "y": 0.708,
                "c": ocean_label_colour,
            },
            {
                "label": "Ad\\'{e}lie\nCoast",
                "x": 0.7,
                "y": 0.65,
                "c": land_label_colour,
            },
            {
                "label": "D'Urville Sea",
                "x": 0.85,
                "y": 0.76,
                "c": ocean_label_colour,
            },
            {
                "label": "Prydz Bay",
                "x": 0.87,
                "y": 0.413,
                "c": ocean_label_colour,
            },
            {
                "label": "Atlantic Ocean",
                "x": 0.244,
                "y": 0.036,
                "c": ocean_label_colour,
            },
            {
                "label": "Pacific Ocean",
                "x": 0.255,
                "y": 0.948,
                "c": ocean_label_colour,
            },
            {
                "label": "Southern Ocean",
                "x": 0.5,
                "y": 0.15,
                "c": ocean_label_colour,
                "fontsize": 16
            },
            {
                "label": "Southern Ocean",
                "x": 0.5,
                "y": 0.85,
                "c": ocean_label_colour,
                "fontsize": 16
            },
        ]

    for label in place_labels:
        ax.text(label["x"], 1 - label["y"], label["label"],
                fontsize=label.get("fontsize", 11), color=label["c"],
                ha="center", va="center", transform=ax.transAxes)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    fig.savefig(f"{plot_dir}/{par.HEMI}_map.png")
    fig.savefig(f"{plot_dir}/{par.HEMI}_map.pdf")
    fig.show()
    plt.close(fig)
