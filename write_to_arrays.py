# Load packages
import datetime as dt
import traceback
from sys import path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

import parameters as par

path.insert(0, par.path)
# path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/')
import grid_set as gs
import data_classes as dc

# %%
# Define variables
hemisphere = par.HEMI
years = par.YEARS
months = par.MONTHS
model = par.MODEL
# %%
### GPLOT SET UP
f = plt.figure()
Gplot = gs.grid_set(par.m)
ax = f.add_subplot(1, 1, 1, projection=par.m)
ax.set_extent(par.geo_bounds, ccrs.PlateCarree())

### make a new grid
Gplot.set_grid_mn(30, 30, ax)
Gplot.get_grid_info(av_ang=False)
plt.close()
# %%
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
lonE = MWinds.f_nc.variables['longitude'][:].data
latE = MWinds.f_nc.variables['latitude'][:].data
lon, lat = np.meshgrid(lonE, latE)
if hemisphere == "north":
    lon[lat < 60] = np.nan
    lat[lat < 60] = np.nan
if hemisphere == "south":
    lon[lat > -55] = np.nan
    lat[lat > -55] = np.nan
GEmonth = gs.grid_set(par.m)
GEmonth.set_grid_lon_lat(lon, lat)
GEmonth.blank_grid_info()
GEmonth.ang_c[:] = 1.0
GEmonth2GPathfinder = gs.Gs2Gs(GEmonth, GPathfinder, vectors=True, NaN_avoid=True)
# %%
# Creating dictionaries to run Olivia's functions
# id_function_dict = {
#     "jan": ID.ID_av.jan,
#     "feb": ID.ID_av.feb,
#     "mar": ID.ID_av.mar,
#     "apr": ID.ID_av.apr,
#     "may": ID.ID_av.may,
#     "jun": ID.ID_av.jun,
#     "jul": ID.ID_av.jul,
#     "aug": ID.ID_av.aug,
#     "sep": ID.ID_av.sep,
#     "oct": ID.ID_av.oct1,
#     "nov": ID.ID_av.nov,
#     "dec": ID.ID_av.dec
# }
# ic_function_dict = {
#     "jan": IC.IC_av.jan,
#     "feb": IC.IC_av.feb,
#     "mar": IC.IC_av.mar,
#     "apr": IC.IC_av.apr,
#     "may": IC.IC_av.may,
#     "jun": IC.IC_av.jun,
#     "jul": IC.IC_av.jul,
#     "aug": IC.IC_av.aug,
#     "sep": IC.IC_av.sep,
#     "oct": IC.IC_av.oct1,
#     "nov": IC.IC_av.nov,
#     "dec": IC.IC_av.dec
# }
# %%
# Initializing arrays of particular shapes

# #ice drift
# SHAPE_ID_UVR = (30, 30) #bogus; residual from Olivia's code
# ice_drift_uv = np.ma.zeros((1, len(months), *par.hem_shape, 2)) #The 1 is for number of years and 2 is for u and v
# ice_drift_uvr = np.ma.zeros((1, len(months), *SHAPE_ID_UVR, 2))

# #ice concentration
# ice_conc_ic = np.ma.zeros((1, len(months), *par.ic_shape)) #1 is for number of years
# ice_conc_icr = np.ma.zeros((1, len(months), *par.hem_shape))

# geostrophic currents
GEO = dc.CPOM_geo(f'{par.path}CPOM_geo/')
gc_uvg = np.ma.zeros((1, len(months), *par.gc_shape, 2))  # 1 is for number of years, 2 is for u and v
gc_uvgr = np.ma.zeros((1, len(months), *par.hem_shape, 2))
print("Got to checkpoint 1")
#
# #winds
# SHAPE_UVW = (121, 1440)
# wind_uvw = np.ma.zeros((1, len(months), *SHAPE_UVW, 2)) #1 is for number of years, 2 is for u and v
# wind_uvwr = np.ma.zeros((1, len(months), *par.hem_shape, 2))
#
# #ekman
# ekman_values = np.ma.zeros((1, len(months), *par.hem_shape, 2, 2)) #1 is for number of years, 2 is for ek and ekg0
# tau_values = np.ma.zeros((1, len(months), *par.hem_shape, 2, 5)) #1 for num years, 5 is for tau_a, tau_i, tau_i0, tau_g, tau_all
# pump_values = np.ma.zeros((1, len(months), *par.hem_shape, 4)) #1 for num years, 4 is for pump_a, pump_i, pump_g, pump
# #%%
# #ice drift
# for year_idx, year in enumerate(years):
#     for month_idx, month in enumerate(months):
#         try:
#             print(f"{month} ice drift")
#             u, v, ur, vr = id_function_dict[month](year)
#             #print('try', u.shape, v.shape, ur.shape, vr.shape)
#             # Take u & v and stack them into a array of size (361, 361, 2) = (*SHAPE_ID_UV, 2)
#             #print("year_idx = ", year_idx, "month_idx = ", month_idx)
#             ice_drift_uv[0, month_idx, :, :] = np.stack([u.squeeze(), v.squeeze()], axis=-1)
#             ice_drift_uvr[0, month_idx, :, :] = np.stack([ur.squeeze(), vr.squeeze()], axis=-1)
#         except ValueError as e:
#             print(f"Caught error '{e}' on: {year}, {month}")
#             print(traceback.format_exc())
#             #print()
#             ice_drift_uv[0, month_idx, :, :] = np.ma.masked
#             ice_drift_uvr[0, month_idx, :, :] = np.ma.masked
#
#     iduv = ice_drift_uv.filled(np.nan)
#     iduvr = ice_drift_uvr.filled(np.nan)
#
#     np.save(f"Data_arrays/{hemisphere}/ice_drift_uv_{year}.npy", iduv)
#     np.save(f"Data_arrays/{hemisphere}/ice_drift_uvr_{year}.npy", iduvr)
#     #print(iduv.shape, iduvr.shape)
# %%
# #ice concentration
# for year_idx, year in enumerate(years):
#     print(f"ICE: year={year}")
#     for month_idx, month in enumerate(months):
#         print(f"\tmonth={month}")
#         try:
#             ic, icr = ic_function_dict[month](year)
#             ice_conc_ic[0, month_idx] = ic.squeeze()
#             ice_conc_icr[0, month_idx] = icr.squeeze()
#         except ValueError as e:
#             print(f"Caught error '{e}' on: {year}, {month}")
#             print(traceback.format_exc())
#             ice_conc_ic[0, month_idx] = np.ma.masked
#             ice_conc_icr[0, month_idx] = np.ma.masked
#
#     ic = ice_conc_ic.filled(np.nan)
#     icr = ice_conc_icr.filled(np.nan)
#
#     np.save(f"Data_arrays/{hemisphere}/ice_conc_ic_{year}.npy", ic)
#     np.save(f"Data_arrays/{hemisphere}/ice_conc_icr_{year}.npy", icr)
# #%%
# #This is a junk loop to create empty geostrophic arrays pre-2011
# #Turn to code when you want to run pre-2011
# #I was getting errors and this was the easiest way around them
#
# for year_idx, year in enumerate(years):
#     uvg = np.zeros((1, 12, *par.gc_shape, 2))
#     uvgr = np.zeros((1, 12, *par.hem_shape, 2))
#     np.save(f"Data_arrays/{hemisphere}/gc_uvg_{year}.npy", uvg)
#     np.save(f"Data_arrays/{hemisphere}/gc_uvgr_{year}.npy", uvgr)
# #%%
# geostrophic currents
print("Got to checkpoint Charlie")
print(years)
print(months)
for year_idx, year in enumerate(years):
    print(f"starting {year}")
    for month_idx, month in enumerate(months):
        print(year, month)
        try:
            ug, vg = GEO.get_vels([dt.datetime(year, month_idx + 1, 1)], verbos=True)
            ugr, vgr = GCPOM2GPathfinder.rg_vecs(ug, vg)
            gc_uvg[0, month_idx] = np.stack([ug.squeeze(), vg.squeeze()], axis=-1)  # 0 for year index
            gc_uvgr[0, month_idx] = np.stack([ugr.squeeze(), vgr.squeeze()], axis=-1)
            print(year, month)
        except ValueError as e:
            print(f"Caught error '{e}' on: {year}, {month}")
            print(traceback.format_exc())
            print()
            gc_uvg[0, month_idx] = np.ma.masked
            gc_uvgr[0, month_idx] = np.ma.masked

    uvg = gc_uvg.filled(np.nan)
    uvgr = gc_uvgr.filled(np.nan)

    np.save(f"Data_arrays/{hemisphere}/gc_uvg_{year}.npy", uvg)
    np.save(f"Data_arrays/{hemisphere}/gc_uvgr_{year}.npy", uvgr)
print(gc_uvg.shape, gc_uvgr.shape)
# #%%
# #wind
# for year_idx, year in enumerate(years):
#     for month_idx, month in enumerate(months):
#         print(f"{month} {year}")
#         try:
#             MWinds.get_dates(dt.datetime(year,month_idx+1,1))
#             uw,vw = MWinds.get_vels()
#             uwr,vwr = GEmonth2GPathfinder.rg_vecs(uw,vw)
#             # print(u.shape, v.shape, ur.shape, vr.shape)
#             # Take u & v and stack them into a array of size (361, 361, 2) = (*SHAPE_ID_UV, 2)
#             wind_uvw[0, month_idx, :, :] = np.stack([uw.squeeze(), vw.squeeze()], axis=-1) #0 for year index
#             wind_uvwr[0, month_idx, :, :] = np.stack([uwr.squeeze(), vwr.squeeze()], axis=-1)
#         except ValueError as e:
#             #print(f"Caught error '{e}' on: {year}, {month}")
#             #print(traceback.format_exc())
#             #print()
#             wind_uvw[0, month_idx, :, :] = np.ma.masked
#             wind_uvwr[0, month_idx, :, :] = np.ma.masked
#
#     uvw = wind_uvw.filled(np.nan)
#     uvwr = wind_uvwr.filled(np.nan)
#
#     np.save(f"Data_arrays/{hemisphere}/wind_uvw_{year}.npy", uvw)
#     np.save(f"Data_arrays/{hemisphere}/wind_uvwr_{year}.npy", uvwr)
#     #print(wind_uvw.shape, wind_uvwr.shape)
# #%%
