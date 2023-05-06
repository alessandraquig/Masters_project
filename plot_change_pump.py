import numpy as np
import matplotlib
from matplotlib.colors import CenteredNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

matplotlib.use('Agg')
matplotlib.rcParams['figure.dpi'] = 400
import matplotlib.pyplot as plt
from sys import path
import parameters as par

path.insert(0, f"{par.path}")
import grid_set as gs
import data_classes as dc
import plotting_functions as pf

hemisphere = par.HEMI
years = par.YEARS
model = "model7"
months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

# GPLOT SET UP
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

# winds
MWinds = dc.ERA5_months(f'{par.path}ERA5/')
GEmonth = gs.grid_set(m)
lonE = MWinds.f_nc.variables['longitude'][:].data
latE = MWinds.f_nc.variables['latitude'][:].data
lon, lat = np.meshgrid(lonE, latE)
if hemisphere == "north":
    lon[lat < 60] = np.nan
    lat[lat < 60] = np.nan
if hemisphere == "south":
    lon[lat > -55] = np.nan
    lat[lat > -55] = np.nan
GEmonth.set_grid_lon_lat(lon, lat)
GEmonth.blank_grid_info()
GEmonth.ang_c[:] = 1.0
GEmonth2GPathfinder = gs.Gs2Gs(GEmonth, GPathfinder, vectors=True, NaN_avoid=True)

# Loading arrays from saved files
pump = np.load(f"Data_arrays/{hemisphere}/model7/pump_1979-2020.npy")
# Assumes I'm using model 7, which I am bc it's the most complete model pre-2011

# Plot change in pumping (Average for 2010s - average for 1980s)
pump_80s = pf.get_average(pump[..., 3], 1979, np.arange(1980, 1989 + 1))
pump_10s = pf.get_average(pump[..., 3], 1979, np.arange(2010, 2019 + 1))
pump_change = pump_10s - pump_80s
pump_change = pf.mask_data(pump_change)

pump_80s_sep = pf.get_average(pump[..., 3], 1979, np.arange(1980, 1989 + 1), months=["sep"])
pump_10s_sep = pf.get_average(pump[..., 3], 1979, np.arange(2010, 2019 + 1), months=["sep"])
pump_change_sep = pump_10s_sep - pump_80s_sep
pump_change_sep = pf.mask_data(pump_change_sep)

pump_80s_mar = pf.get_average(pump[..., 3], 1979, np.arange(1980, 1989 + 1), months=["mar"])
pump_10s_mar = pf.get_average(pump[..., 3], 1979, np.arange(2010, 2019 + 1), months=["mar"])
pump_change_mar = pump_10s_mar - pump_80s_mar
pump_change_mar = pf.mask_data(pump_change_mar)

pump_norm = CenteredNorm(halfrange=25, vcenter=0)
fig_pump = plt.figure(figsize=(18, 9))
fig_pump.suptitle("Change in Average Pumping 1980s-2010s")
ax_pump, img_pump = pf.colorplot(data=pump_change, fig=fig_pump, rows=1, cols=3, pos=1, title="Whole Year", norm=pump_norm)
ax_pump_sep, _ = pf.colorplot(data=pump_change_sep, fig=fig_pump, rows=1, cols=3, pos=2, title="September", norm=pump_norm)
ax_pump_mar, _ = pf.colorplot(data=pump_change_mar, fig=fig_pump, rows=1, cols=3, pos=3, title="March", norm=pump_norm)

fig_pump.colorbar(img_pump, ax=[ax_pump, ax_pump_sep, ax_pump_mar], orientation="horizontal",
                  pad=0.05, aspect=25, norm=pump_norm)
fig_pump.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)

fig_pump.savefig(f"Maps_output/{hemisphere}/{model}/Pump_change_seasons.png")