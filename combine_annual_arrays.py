import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
import pandas as pd
import traceback
import importlib
#import xarrary as xr
from sys import path
from os.path import exists
#from imp import reload
from importlib import reload
from netCDF4 import Dataset
# path.insert(0, '/Users/H/WAVES/geo_data_group/')
# path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project')
import grid_set as gs
import data_classes as dc
import parameters as par
#%%
hemisphere = "north"
years = np.arange(1979, 2020+1)
model = "model6"
months = ["jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]
path = par.path

### GPLOT SET UP

if hemisphere == "north":
    m = ccrs.NorthPolarStereo()

    #### setup plotting grid - Gplot
    f = plt.figure()
    Gplot= gs.grid_set(m)

    ax = f.add_subplot(1,1,1,projection=m)

    #for north
    ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())

    ### make a new grid
    Gplot.set_grid_mn(30,30,ax)
    Gplot.get_grid_info(av_ang=False)
    plt.close()

if hemisphere == "south":
    m = ccrs.SouthPolarStereo()

    #### setup plotting grid - Gplot
    f = plt.figure()
    Gplot= gs.grid_set(m)

    ax = f.add_subplot(1,1,1,projection=m)

    #for south
    ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())

    ### make a new grid
    Gplot.set_grid_mn(30,30,ax)
    Gplot.get_grid_info(av_ang=False)
    plt.close()
#%%
#opening files, getting native grids

#ice drift
DRIFT = dc.Pathfinder(f'{path}Pathfinder/')

GPathfinder = gs.grid_set(m)
GPathfinder.load_grid(par.ID_grid)
GPathfinder2Gplot = gs.Gs2Gs(GPathfinder,Gplot,vectors=True)


#ice concentration
IConc = dc.NSIDC_nt(f'{path}NSIDC_nt')

GIC = gs.grid_set(m)
if hemisphere == "north":
    GIC.load_grid(f'{path}NSIDC_gs.npz')
if hemisphere == "south":
    GIC.load_grid(f'{path}NSIDC_gs_SH.npz')
GIC2GPathfinder = gs.Gs2Gs(GIC,GPathfinder,vectors=False)


#geo currents
GEO = dc.CPOM_geo(f'{path}CPOM_geo/')

GCPOM = gs.grid_set(m)
if hemisphere == "north":
    GCPOM.load_grid(f'{path}PS_20km_gs2021.npz')
if hemisphere == "south":
    GCPOM.load_grid(f'{path}Polar_stereo_50km_SH.npz')
GCPOM2GPathfinder = gs.Gs2Gs(GCPOM,GPathfinder,vectors=True)


#winds
MWinds = dc.ERA5_months(f'{path}ERA5/')

if hemisphere == "north":
    lonE = MWinds.f_nc.variables['longitude'][:].data
    latE = MWinds.f_nc.variables['latitude'][:].data
    lon,lat = np.meshgrid(lonE,latE)
    lon[lat<60] = np.nan
    lat[lat<60] = np.nan
    GEmonth = gs.grid_set(m)
    GEmonth.set_grid_lon_lat(lon,lat)
    GEmonth.blank_grid_info()
    GEmonth.ang_c[:] = 1.0
if hemisphere == "south":
    lonE = MWinds.f_nc.variables['longitude'][:].data
    latE = MWinds.f_nc.variables['latitude'][:].data
    lon,lat = np.meshgrid(lonE,latE)
    lon[lat>-55] = np.nan
    lat[lat>-55] = np.nan
    GEmonth = gs.grid_set(m)
    GEmonth.set_grid_lon_lat(lon,lat)
    GEmonth.blank_grid_info()
    GEmonth.ang_c[:] = 1.0
GEmonth2GPathfinder = gs.Gs2Gs(GEmonth,GPathfinder,vectors=True,NaN_avoid=True)

#%% md
### MONTHLY AVERAGES
#%% md
##### Loading arrays from saved files
#%%
#Making whole drift array
# iduv_list = []
# for year_idx, year in enumerate(years):
#     #if year != 1993 and year != 1994 and year != 1995:
#     try:
#         iduv = np.load(f"Annual_data/{hemisphere.capitalize()}/ice_drift_uv_{year}.npy")
#         iduv_list.append(iduv)
#     except:
#         print("No data for ", year)
#         iduv_empty = np.zeros((1, 12, len(GPathfinder.xpts), len(GPathfinder.ypts), 2))
#         iduv_list.append(iduv_empty)
#
# drift = np.concatenate(iduv_list, axis=0)
# np.save(f"Data_arrays/{hemisphere}/drift_{years[0]}-{years[-1]}.npy", drift)
#%%
#Making whole concentration array
icr_list = []
for year_idx, year in enumerate(years):
    #if year != 1993 and year != 1994 and year != 1995:
    try:
        icr = np.load(f"Data_arrays/{hemisphere}/ice_conc_icr_{year}.npy")
        icr_list.append(icr)
    except:
        print("No data for ",year)
        icr_empty = np.zeros((1, 12, len(GPathfinder.xpts), len(GPathfinder.ypts)))
        icr_list.append(icr_empty)

concentration = np.concatenate(icr_list, axis=0)
np.save(f"Data_arrays/{hemisphere}/conc_{years[0]}-{years[-1]}.npy", concentration)
# #%%
#Making whole geostrophic array
#I've made arrays of zeros for 1979-2010. This allows me to run the code without it crashing, but keep in mind that only models not using geostrophic data work before then
# gcr_list = []
# for year_idx, year in enumerate(years):
#     #if year != 1993 and year != 1994 and year != 1995:
#     try:
#         gcr = np.load(f"Annual_data/{hemisphere.capitalize()}/gc_uvgr_{year}.npy")
#         gcr_list.append(gcr)
#     except:
#         print("No data for ",year)
#         gcr_empty = np.zeros((1, 12, len(GPathfinder.xpts), len(GPathfinder.ypts), 2))
#         gcr_list.append(gcr_empty)
#
# geostrophic = np.concatenate(gcr_list, axis=0)
# np.save(f"Data_arrays/{hemisphere}/geo_{years[0]}-{years[-1]}.npy", geostrophic)
# #%%
# #Making whole wind array
# windr_list = []
# for year_idx, year in enumerate(years):
#     #if year != 1993 and year != 1994 and year != 1995:
#     try:
#         windr = np.load(f"Annual_data/{hemisphere.capitalize()}/wind_uvwr_{year}.npy")
#         windr_list.append(windr)
#     except:
#         print("No data for ",year)
#         windr_empty = np.zeros((1, 12, len(GPathfinder.xpts), len(GPathfinder.ypts), 2))
#         windr_list.append(windr_empty)
#
# wind = np.concatenate(windr_list, axis=0)
# np.save(f"Data_arrays/{hemisphere}/wind_{years[0]}-{years[-1]}.npy", wind)
#%%
# #Making whole Ekman arrays
# ek_list = []
# tau_list = []
# pump_list = []
#
# for year_idx, year in enumerate(years):
#     #if year != 1993 and year != 1994 and year != 1995:
#     #try:
#     #These have a model number, because the model changes the calculation
#     ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_{year}.npy")
#     ek = ekman[:,:,:,:,0]
#     tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_{year}.npy")
#     t = tau[:,:,:,:,0]
#     pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_{year}.npy")
#     p = pump[:,:,:,:,0]
#     print(np.shape(ek))
#     ek_list.append(ek)
#     tau_list.append(t)
#     pump_list.append(p)
#     #except:
#     #    print("No data for ",year)
#     #    windr_empty = np.zeros((1, 12, len(GPathfinder.xpts), len(GPathfinder.ypts), 2))
#     #    windr_list.append(windr_empty)
#
# ekman = np.concatenate(ek_list, axis=0)
# tau = np.concatenate(tau_list, axis=0)
# pump = np.concatenate(pump_list, axis=0)
#
# np.save(f"Data_arrays/{hemisphere}/{model}/ekman_{years[0]}-{years[-1]}", ekman)
# np.save(f"Data_arrays/{hemisphere}/{model}/tau_{years[0]}-{years[-1]}", tau)
# np.save(f"Data_arrays/{hemisphere}/{model}/pump_{years[0]}-{years[-1]}", pump)
# #%%
# ek = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_{year}.npy")
# p = np.load(f"Data_arrays/{hemisphere}/{model}/pump_{year}.npy")
# t = np.load(f"Data_arrays/{hemisphere}/{model}/tau_{year}.npy")
# print(np.shape(ek), np.shape(p), np.shape(t))
#
# print(np.shape(ekman), np.shape(pump), np.shape(tau))