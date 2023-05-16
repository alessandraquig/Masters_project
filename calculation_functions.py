import numpy as np
import matplotlib
from scipy.ndimage import gaussian_filter1d

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path
import parameters as par
import plotting_functions as pf

path.insert(0, f"{par.path}")

import grid_set as gs
import data_classes as dc

hemisphere = par.HEMI
model = par.MODEL

### GPLOT SET UP
# m = par.m
# bounds = par.geo_bounds
# f = plt.figure()
# Gplot = gs.grid_set(m)
# ax = f.add_subplot(1, 1, 1, projection=m)
# ax.set_extent(bounds, ccrs.PlateCarree())
# Gplot.set_grid_mn(30, 30, ax)
# Gplot.get_grid_info(av_ang=False)
# cell_area_tuple = Gplot.xdist * Gplot.ydist
# cell_area = np.array(cell_area_tuple).mean()
# # print(f"xdist = {Gplot.xdist}, ydist = {Gplot.ydist}")
# # print(f"cell area shape = {cell_area.shape()}")
# plt.close()
#
# # ice drift
# DRIFT = dc.Pathfinder(f'{par.path}Pathfinder/')
# GPathfinder = gs.grid_set(m)
# GPathfinder.load_grid(par.ID_grid)
# GPathfinder2Gplot = gs.Gs2Gs(GPathfinder, Gplot, vectors=True)
#
# # ice concentration
# IConc = dc.NSIDC_nt(f'{par.path}NSIDC_nt/')
# GIC = gs.grid_set(m)
# GIC.load_grid(par.IC_grid)
# GIC2GPathfinder = gs.Gs2Gs(GIC, GPathfinder, vectors=False)
#
# # geo currents
# GEO = dc.CPOM_geo(f'{par.path}CPOM_geo/')
# GCPOM = gs.grid_set(m)
# GCPOM.load_grid(par.GC_grid)
# GCPOM2GPathfinder = gs.Gs2Gs(GCPOM, GPathfinder, vectors=True)
#
# # winds
# MWinds = dc.ERA5_months(f'{par.path}ERA5/')
# GEmonth = gs.grid_set(m)
# lonE = MWinds.f_nc.variables['longitude'][:].data
# latE = MWinds.f_nc.variables['latitude'][:].data
# lon, lat = np.meshgrid(lonE, latE)
# if hemisphere == "north":
#     lon[lat < 60] = np.nan
#     lat[lat < 60] = np.nan
# if hemisphere == "south":
#     lon[lat > -55] = np.nan
#     lat[lat > -55] = np.nan
# GEmonth.set_grid_lon_lat(lon, lat)
# GEmonth.blank_grid_info()
# GEmonth.ang_c[:] = 1.0
# GEmonth2GPathfinder = gs.Gs2Gs(GEmonth, GPathfinder, vectors=True, NaN_avoid=True)

##### Loading arrays from saved files
drift = np.load(f"Data_arrays/{hemisphere}/drift_1979-2020.npy")
conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
geo = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")
wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")

if model == "model3" or model == "model7":
    ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_1979-2020.npy")
    tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_1979-2020.npy")
    pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_1979-2020.npy")
else:
    ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_2011-2020.npy")
    tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_2011-2020.npy")
    pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_2011-2020.npy")

# def get_ice_area(year, month, region=None):
#     if region is not None:
#         # TODO: support regions
#         raise ValueError("Regions not supported")
#
#     ice_area = 0
#     for i in range(par.hem_shape[0]):
#         for j in range(par.hem_shape[1]):
#             if np.isfinite(conc[year, month, i, j]) and conc[year, month, i, j] > 0:
#                 cell_ice_area = cell_area * conc[year, month, i, j]
#                 ice_area = ice_area + cell_ice_area
#                 ice_area = np.sum(ice_area)
#     print(f"ice area = {ice_area} km2")
#     return ice_area

def get_ice_area_old(year, month, conc=conc, region=None, cell_area=1.0):
    if region is not None:
        # TODO: support regions
        raise ValueError("Regions not supported")

    # if hemisphere == "south":
    #     cell_area = 254813.02741155517 * 260349.8460127463    # m, not km

    conc = conc[year, month, ...]  # 2D grid of ice concentrations
    cell_area = 25.067525*25.067525  # km2 value from Summary of NOAA/NASA Polar Pathfinder Grid Relationships

    ice_area = 0
    for i in range(conc.shape[0]):
        for j in range(conc.shape[1]):
            if np.isfinite(conc[i, j]) and conc[i, j] > 0:
                cell_ice_area = cell_area * conc[i, j]
                ice_area += cell_ice_area

    #ice_area_km2 = ice_area / 1e6  # convert square meters to square kilometers
    print(f"ice area = {ice_area} km2")
    return ice_area

def get_ice_area(conc):
    ice_area = 0
    cell_area = 25.067525*25.067525  # km2 value from Summary of NOAA/NASA Polar Pathfinder Grid Relationships

    for i in range(conc.shape[0]):
        for j in range(conc.shape[1]):
            if np.isfinite(conc[i, j]) and conc[i, j] > 0:
                cell_ice_area = cell_area * conc[i, j]
                ice_area += cell_ice_area

    #ice_area_km2 = ice_area / 1e6  # convert square meters to square kilometers
    print(f"ice area = {ice_area} km2")
    return ice_area

if __name__ == '__main__':
    print(f"Shape = {np.shape(pf.get_average(conc, 1979, years=np.arange(1980, 1989+1)))}")

    ice_1983_apr = get_ice_area_old(5, 3, conc=conc)
    ice_1980s = get_ice_area(pf.get_average(conc, 1979, years=np.arange(1980, 1989+1)))
    ice_2010s = get_ice_area(pf.get_average(conc, 1979, years=np.arange(2010, 2019+1)))

#%%
# Setting grid squares for Beaufort Gyre (taken from Olivia's code - talk to Michel)
xmax = 180
xmin = 125
ymax = 257
ymin = 217

pump_80s = pf.get_average(pump[..., xmin:xmax, ymin:ymax, 3], 1979, years=np.arange(1980, 1989+1))
pump_80s_avg = np.nanmean(pump_80s)
pump_10s = pf.get_average(pump[..., xmin:xmax, ymin:ymax, 3], 1979, years=np.arange(2010, 2019+1))
pump_10s_avg = np.nanmean(pump_10s)
pump_change = pump_10s - pump_80s
pump_change_avg = np.nanmean(pump_change)

pump1 = np.load(f"Data_arrays/north/model1/pump_2011-2020.npy")
avg_pump_mag = np.nanmean(pf.get_average(pump1[..., xmin:xmax, ymin:ymax, 3], 2011, years=np.arange(2011, 2020)))

pump4 = np.load(f"Data_arrays/north/model4/pump_2011-2020.npy")
avg_pump_no_ice = np.nanmean(pf.get_average(pump4[..., xmin:xmax, ymin:ymax, 3], 2011, years=np.arange(2011, 2020)))
print(f"No ice pumping = {avg_pump_no_ice}")

# ekman = np.load(f"Data_arrays/north/model7/ekman_1979-2020.npy")
# ekman_80s = pf.get_average(ekman[..., xmin:xmax, ymin:ymax, ..., 0], 1979, years=np.arange(1980, 1989+1))
# ekman_80s_avg = np.nanmean(ekman_80s)
# ekman_10s = pf.get_average(ekman[..., xmin:xmax, ymin:ymax, ..., 0], 1979, years=np.arange(2010, 2019+1))
# ekman_10s_avg = np.nanmean(ekman_10s)
# ekman_change = ekman_10s - ekman_80s
# ekman_change_avg = np.nanmean(ekman_change)

wind = np.abs(np.load(f"Data_arrays/north/wind_1979-2020.npy"))
winter_wind = np.nanmean(pf.get_average(wind[..., xmin:xmax, ymin:ymax, :], 1979, months=["sep", "oct", "nov", "dec", "jan", "feb"]))
summer_wind = np.nanmean(pf.get_average(wind[..., xmin:xmax, ymin:ymax, :], 1979, months=["mar", "apr", "may", "jun", "jul", "aug"]))
print(f"Winter wind = {winter_wind}, summer wind = {summer_wind}")