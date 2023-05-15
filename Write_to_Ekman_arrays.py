# Load packages
import numpy as np
import traceback
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import parameters as par
from sys import path
import time

# path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/')
path.insert(0, par.path)
import grid_set as gs
import Calculating_pumping as CP
import data_classes as dc

#
# Define variables
years = par.YEARS
months = par.MONTHS
# model = par.MODEL
models = ["model1", "model2", "model4", "model5", "model6"]

### GPLOT SET UP
f = plt.figure()
Gplot = gs.grid_set(par.m)
ax = f.add_subplot(1, 1, 1, projection=par.m)
ax.set_extent(par.geo_bounds, ccrs.PlateCarree())

### make a new grid
Gplot.set_grid_mn(30, 30, ax)
Gplot.get_grid_info(av_ang=False)
plt.close()

# opening files, getting native grids

# ice drift
GPathfinder = gs.grid_set(par.m)
GPathfinder.load_grid(par.ID_grid)
GPathfinder2Gplot = gs.Gs2Gs(GPathfinder, Gplot, vectors=True)

# ice concentration
GIC = gs.grid_set(par.m)
GIC.load_grid(par.IC_grid)
GIC2GPathfinder = gs.Gs2Gs(GIC, GPathfinder, vectors=False)

# geo currents
GCPOM = gs.grid_set(par.m)
GCPOM.load_grid(par.GC_grid)
GCPOM2GPathfinder = gs.Gs2Gs(GCPOM, GPathfinder, vectors=True)

# winds
MWinds = dc.ERA5_months(f'{par.path}ERA5/')
GEmonth = gs.grid_set(par.m)
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

# Initializing arrays of particular shapes
# ekman_values = np.ma.zeros((1, len(months), *par.hem_shape, 2, 2))  # 1 is for number of years, 2 is for ek and ekg0
# tau_values = np.ma.zeros(
#     (1, len(months), *par.hem_shape, 2, 5))  # 1 for num years, 5 is for tau_a, tau_i, tau_i0, tau_g, tau_all
# pump_values = np.ma.zeros((1, len(months), *par.hem_shape, 4))  # 1 for num years, 4 is for pump_a, pump_i, pump_g, pump

for i in range(len(models)):
    # model = models[i]
    par.MODEL = models[i]
    if par.MODEL == "model3" or par.MODEL == "model7":
        years = np.arange(1979, 2020 + 1)
    else:
        years = np.arange(2011, 2020 + 1)

    # Initializing arrays of particular shapes
    ekman_values = np.ma.zeros(
        (len(years), len(months), *par.hem_shape, 2, 2))  # 1 is for number of years, 2 is for ek and ekg0
    tau_values = np.ma.zeros((len(years), len(months), *par.hem_shape, 2,
                              5))  # 1 for num years, 5 is for tau_a, tau_i, tau_i0, tau_g, tau_all
    pump_values = np.ma.zeros(
        (len(years), len(months), *par.hem_shape, 4))  # 1 for num years, 4 is for pump_a, pump_i, pump_g, pump

    # Loads arrays: this is so I can collect the data and then calculate separately
    # Saves me from having to collect the data each time I run this
    drift = np.load(f"Data_arrays/{par.HEMI}/drift_1979-2020.npy")
    conc = np.load(f"Data_arrays/{par.HEMI}/conc_1979-2020.npy")
    geo = np.load(f"Data_arrays/{par.HEMI}/geo_1979-2020.npy")
    geo[..., 0] *= -1
    geo[..., 1] *= -1
    wind = np.load(f"Data_arrays/{par.HEMI}/wind_1979-2020.npy")

    land_mask = np.isfinite(geo[-1, -1, :, :, 0]) & np.isfinite(geo[-1, -1, :, :, 1])
    for year_idx, year in enumerate(years):
        for month_idx, month in enumerate(months):
            start_time = time.perf_counter()
            try:
                # TODO: replace with CP.current_fast and see if it works
                # TODO: mask the geostrophic array with land mask pre 2011
                ek, ekg0, tau_all, tau_a, tau_i, taui0, tau_g = CP.current(
                    wind[year_idx, month_idx, :, :, 0],  # 0 for u direction
                    wind[year_idx, month_idx, :, :, 1],  # 1 for v direction
                    drift[year_idx, month_idx, :, :, 0],  # 0 for u
                    drift[year_idx, month_idx, :, :, 1],  # 1 for v
                    geo[year_idx, month_idx, :, :, 0],  # 0 for u
                    geo[year_idx, month_idx, :, :, 1],  # 1 for v
                    conc[year_idx, month_idx, :, :]  # No u or v b/c not directional
                )
                pump_a, pump_i, pump_g, pump = CP.pump(tau_a, tau_i, tau_g, tau_all)
                # print(u.shape, v.shape, ur.shape, vr.shape)
                # Take u & v and stack them into a array of size (361, 361, 2) = (*SHAPE_ID_UV, 2)
                # Indexed as ekman_values[year, month, grid_u, grid_v, (u/v), (ek/ekg0)]
                ekman_values[year_idx, month_idx, :, :, :] = np.stack([ek.squeeze(), ekg0.squeeze()], axis=-1)
                tau_values[year_idx, month_idx, :, :, :] = np.stack(
                    [tau_all.squeeze(), tau_a.squeeze(), tau_i.squeeze(), taui0.squeeze(), tau_g.squeeze()], axis=-1)
                pump_values[year_idx, month_idx, :, :] = np.stack(
                    [pump_a.squeeze(), pump_i.squeeze(), pump_g.squeeze(), pump.squeeze()], axis=-1)
            except ValueError as e:
                print(f"Caught error '{e}' on: {year}, {month}")
                print(traceback.format_exc())
            end_time = time.perf_counter()
            print(f"{year} {month} complete for {par.MODEL}: {end_time - start_time:.2f}s")
    ekman_array = ekman_values.filled(np.nan)
    tau_array = tau_values.filled(np.nan)
    pump_array = pump_values.filled(np.nan)

    np.save(f"Data_arrays/{par.HEMI}/{par.MODEL}/ekman_{years[0]}-2020.npy", ekman_array)
    np.save(f"Data_arrays/{par.HEMI}/{par.MODEL}/tau_{years[0]}-2020.npy", tau_array)
    np.save(f"Data_arrays/{par.HEMI}/{par.MODEL}/pump_{years[0]}-2020.npy", pump_array)
    print(par.MODEL, "saved")
