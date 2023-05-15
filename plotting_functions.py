import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter1d

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path
import parameters as par

path.insert(0, f"{par.path}")

import grid_set as gs
import data_classes as dc

### GPLOT SET UP
m = par.m
bounds = par.geo_bounds
f = plt.figure()
Gplot = gs.grid_set(m)
ax = f.add_subplot(1, 1, 1, projection=m)
ax.set_extent(bounds, ccrs.PlateCarree())
Gplot.set_grid_mn(30, 30, ax)
Gplot.get_grid_info(av_ang=False)
plt.close()

# ice drift
DRIFT = dc.Pathfinder(f'{par.path}Pathfinder/')
GPathfinder = gs.grid_set(m)
GPathfinder.load_grid(par.ID_grid)
GPathfinder2Gplot = gs.Gs2Gs(GPathfinder, Gplot, vectors=True)

# ice concentration
IConc = dc.NSIDC_nt(f'{par.path}NSIDC_nt/')
GIC = gs.grid_set(m)
GIC.load_grid(par.IC_grid)
GIC2GPathfinder = gs.Gs2Gs(GIC, GPathfinder, vectors=False)

# geo currents
GEO = dc.CPOM_geo(f'{par.path}CPOM_geo/')
GCPOM = gs.grid_set(m)
GCPOM.load_grid(par.GC_grid)
GCPOM2GPathfinder = gs.Gs2Gs(GCPOM, GPathfinder, vectors=True)
GCPOM2Gplot = gs.Gs2Gs(GCPOM, Gplot, vectors=True)

# winds
MWinds = dc.ERA5_months(f'{par.path}ERA5/')
GEmonth = gs.grid_set(m)
lonE = MWinds.f_nc.variables['longitude'][:].data
latE = MWinds.f_nc.variables['latitude'][:].data
lon, lat = np.meshgrid(lonE, latE)
if par.HEMI == "north":
    lon[lat < 60] = np.nan
    lat[lat < 60] = np.nan
if par.HEMI == "south":
    lon[lat > -55] = np.nan
    lat[lat > -55] = np.nan
GEmonth.set_grid_lon_lat(lon, lat)
GEmonth.blank_grid_info()
GEmonth.ang_c[:] = 1.0
GEmonth2GPathfinder = gs.Gs2Gs(GEmonth, GPathfinder, vectors=True, NaN_avoid=True)

geo_jan2018 = np.load(f"{par.path}Data_arrays/{par.HEMI}/geo_1979-2020.npy")[-2, 0]
land_mask = np.isfinite(geo_jan2018[..., 0]) & np.isfinite(geo_jan2018[..., 1])


# Here are the broadly generalizable functions
def mask_data(data):
    """
    Applies a land mask - sets all land locations to nan. Run this for everything except wind.
    :param data:
    :return:
    """
    data[~land_mask] = np.nan
    return data


def get_average(data, data_start, years=None, months=None):
    """
    Takes an array and a range of years and months
    Gives the average for each grid cell over that time period
    If you don't specify months and/or years, it will take the values from parameters.py

    :param data: Data array - your choices are drift, conc, geo, wind, ekman[..., 0], or pump[..., 3]
    :param data_start: The first year of the available (not necessarily used) dataset (1979 or 2011)
    :param years: Your range of years (does not include final year, e.g. put 2021 to include 2020)
    :param months: Your range of months (does include all months)
    :return: Data array with one average value for each grid cell
    """
    if months is None:
        months = par.MONTHS
    if years is None:
        years = par.YEARS
    months_idxs = [par.MONTHS.index(mon) for mon in months]
    print("Getting average")
    print(f"months_idxs = {months_idxs}")
    start_idx = years[0] - data_start
    end_idx = years[-1] - data_start
    print(f"start_idx, end_idx: {start_idx} = {years[0]}, {end_idx} = {years[-1]}")

    data_avg = np.nanmean(data[start_idx:end_idx + 1, months_idxs], axis=(0, 1))
    return data_avg


def colorplot(data, fig, title, rows=1, cols=1, pos=1, cmap="PiYG", norm=None):
    """
    This will plot a variable with no vectors, just a single value for each grid square
    Works for conc or pump

    :param data: This is an average calculated in get_average or a month and year
    :param fig: Define a figure before you call this function. This is so you can use this function to
    call subplots in case you're plotting several things
    :param rows: How many rows of subplots in your figure (default is 1)
    :param cols: How many columns of subplots in your figure (default is 1)
    :param pos: The position of this subplot within your figure (default is 1)
    :param title: Title of this subplot (not the whole figure - define that outside)
    :param cmap: Colormap - default is an ugly one so I remember to change it for each variable
    :return: Returns the axis for a subplot
    """
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=m)
    ax.set_extent(bounds, ccrs.PlateCarree())
    ax.set_title(title)
    img = ax.pcolormesh(GPathfinder.xpts, GPathfinder.ypts, data, cmap=cmap, norm=norm)
    ax.add_feature(cfeature.COASTLINE)
    return ax, img


def vectorplot(data, fig, title, rows=1, cols=1, pos=1, cmap="PiYG", arrow_every=10, norm=None):
    """
    This will plot a variable with vectors
    Works for drift, wind, geo, or ekman

    :param data: This is either an unprocessed data array (drift, conc, geo, wind, ekman[..., 0], or pump[..., 3])
    or an average calculated in get_average
    :param fig: Define a figure before you call this function. This is so you can use this function to
    call subplots in case you're plotting several things
    :param rows: How many rows of subplots in your figure (default is 1)
    :param cols: How many columns of subplots in your figure (default is 1)
    :param pos: The position of this subplot within your figure (default is 1)
    :param title: Title of this subplot (not the whole figure - define that outside)
    :param cmap: Colormap - default is an ugly one so I remember to change it for each variable
    :return: Returns the axis to be plotted
    """

    # calculate the magnitude of the averaged Ekman current field using np.hypot
    ur_avg, vr_avg = data[..., 0], data[..., 1]
    mag = np.hypot(ur_avg, vr_avg)

    # plot the magnitude of the averaged Ekman current field
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=m)
    ax.set_extent(bounds, ccrs.PlateCarree())
    img = ax.contourf(GPathfinder.xpts, GPathfinder.ypts, mag, cmap=cmap)

    if par.HEMI == "north":
        ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                  GPathfinder.ypts[::arrow_every, ::arrow_every],
                  ur_avg[::arrow_every, ::arrow_every],
                  vr_avg[::arrow_every, ::arrow_every], color="k")

    if par.HEMI == "south":
        ax.quiver(Gplot.xpts, Gplot.ypts,
                  *GPathfinder2Gplot.rg_vecs_to_plot(ur_avg, vr_avg),
                  color="k")
    ax.add_feature(cfeature.COASTLINE)
    ax.set_title(title)
    return ax, img


def vectorplot2(data1, data2, fig, title, rows=1, cols=1, pos=1, cmap="PiYG", arrow_every=10, norm=None):
    """
    This will plot a variable with vectors
    Works for drift, wind, geo, or ekman

    :param data: This is either an unprocessed data array (drift, conc, geo, wind, ekman[..., 0], or pump[..., 3])
    or an average calculated in get_average
    :param fig: Define a figure before you call this function. This is so you can use this function to
    call subplots in case you're plotting several things
    :param rows: How many rows of subplots in your figure (default is 1)
    :param cols: How many columns of subplots in your figure (default is 1)
    :param pos: The position of this subplot within your figure (default is 1)
    :param title: Title of this subplot (not the whole figure - define that outside)
    :param cmap: Colormap - default is an ugly one so I remember to change it for each variable
    :return: Returns the axis to be plotted
    """

    # calculate the magnitude of the averaged Ekman current field using np.hypot
    ur_avg = data1
    vr_avg = data2

    mag = np.hypot(ur_avg, vr_avg)

    # plot the magnitude of the averaged Ekman current field
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=m)
    ax.set_extent(bounds, ccrs.PlateCarree())
    img = ax.contourf(GCPOM.xpts[::arrow_every, ::arrow_every], GCPOM.ypts[::arrow_every, ::arrow_every],
                      mag[::arrow_every, ::arrow_every], cmap=cmap)
    ax.quiver(GCPOM.xpts[::arrow_every, ::arrow_every],
              GCPOM.ypts[::arrow_every, ::arrow_every],
              ur_avg[::arrow_every, ::arrow_every],
              vr_avg[::arrow_every, ::arrow_every], color="k")
    ax.add_feature(cfeature.COASTLINE)
    ax.set_title(title)
    return ax, img


def vectorplot_old(data, fig, title, rows=1, cols=1, pos=1, cmap="PiYG", arrow_every=10, norm=None):
    """
    This will plot a variable with vectors
    Works for drift, wind, geo, or ekman

    :param data: This is either an unprocessed data array (drift, conc, geo, wind, ekman[..., 0], or pump[..., 3])
    or an average calculated in get_average
    :param fig: Define a figure before you call this function. This is so you can use this function to
    call subplots in case you're plotting several things
    :param rows: How many rows of subplots in your figure (default is 1)
    :param cols: How many columns of subplots in your figure (default is 1)
    :param pos: The position of this subplot within your figure (default is 1)
    :param title: Title of this subplot (not the whole figure - define that outside)
    :param cmap: Colormap - default is an ugly one so I remember to change it for each variable
    :return: Returns the axis to be plotted
    """

    # calculate the magnitude of the averaged Ekman current field using np.hypot
    ur_avg = data[..., 0]
    vr_avg = data[..., 1]

    mag = np.hypot(ur_avg, vr_avg)

    # plot the magnitude of the averaged Ekman current field
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=m)
    ax.set_extent(bounds, ccrs.PlateCarree())
    img = ax.contourf(GPathfinder.xpts, GPathfinder.ypts, mag, cmap=cmap)
    ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
              GPathfinder.ypts[::arrow_every, ::arrow_every],
              ur_avg[::arrow_every, ::arrow_every],
              vr_avg[::arrow_every, ::arrow_every], color="k")
    ax.add_feature(cfeature.COASTLINE)
    ax.set_title(title)
    return ax, img


def smooth_timeseries(timeseries: np.ndarray, scale, window_type: str = "gaussian"):
    """
    Smooth input ``timeseries`` by a window with length ``scale``, using convolution.

    :param timeseries: 1D input to be smoothed
    :param scale: Length/stdev of the smoothing window.
    :param window_type: Type of window. Options: gaussian, box
    :return: Returns the smoothed input
    """
    if window_type == "box":
        window = np.ones(scale) / scale  # tophat function
        output = np.convolve(timeseries, window, mode="same")
    elif window_type == "gaussian":
        output = gaussian_filter1d(timeseries, sigma=scale, mode="nearest")
    else:
        raise ValueError(f"Unrecognised window {window_type}")
    return output
