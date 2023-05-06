import numpy as np
import matplotlib

matplotlib.use('Agg')
matplotlib.rcParams['figure.dpi'] = 400

import plotting_functions as pf
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import grid_set as gs
import data_classes as dc
from sys import path
import parameters as par

path.insert(0, f"{par.path}")

hemisphere = par.HEMI
years = np.arange(2011, 2020 + 1)
months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]

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

##### Loading arrays from saved files
conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")

# Get averages of wind and conc (over 2011 to 2020)
wind_avg = pf.mask_data(pf.get_average(wind, 1979, years))
conc_avg = pf.mask_data(pf.get_average(conc, 1979, years))
ur_avg = wind_avg[..., 0]
vr_avg = wind_avg[..., 1]

# Winter
wind_avg_jas = pf.mask_data(pf.get_average(wind, 1979, years, ["jul", "aug", "sep"]))
conc_avg_jas = pf.mask_data(pf.get_average(conc, 1979, years, ["jul", "aug", "sep"]))
ur_avg_jas = wind_avg[..., 0]
vr_avg_jas = wind_avg[..., 1]

# Summer
wind_avg_jfm = pf.mask_data(pf.get_average(wind, 1979, years, ["jan", "feb", "mar"]))
conc_avg_jfm = pf.mask_data(pf.get_average(conc, 1979, years, ["jan", "feb", "mar"]))
ur_avg_jfm = wind_avg[..., 0]
vr_avg_jfm = wind_avg[..., 1]

# plot vector field of wind over colorplot of ice concentration
fig = plt.figure(figsize=(18, 9))
fig.suptitle("Average Ice Concentration and Wind 2011-2020")
arrow_every = 10  # scales how close together arrows are
ice_cmap = "Blues"

# Whole year
ax: plt.Axes = fig.add_subplot(1, 3, 1, projection=m)
ax.set_extent(bounds, ccrs.PlateCarree())
s = ax.contourf(GPathfinder.xpts, GPathfinder.ypts, conc_avg, np.linspace(0, 1, 10 + 1),
                cmap=ice_cmap)
ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
          GPathfinder.ypts[::arrow_every, ::arrow_every],
          ur_avg[::arrow_every, ::arrow_every],
          vr_avg[::arrow_every, ::arrow_every], color="r")
ax.add_feature(cfeature.COASTLINE)
ax.set_title("Whole Year")

# Summer - jfm
ax_summer: plt.Axes = fig.add_subplot(1, 3, 2, projection=m)
ax_summer.set_extent(bounds, ccrs.PlateCarree())
s = ax_summer.contourf(GPathfinder.xpts, GPathfinder.ypts, conc_avg_jfm, np.linspace(0, 1, 10 + 1),
                       cmap=ice_cmap)
ax_summer.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                 GPathfinder.ypts[::arrow_every, ::arrow_every],
                 ur_avg_jfm[::arrow_every, ::arrow_every],
                 vr_avg_jfm[::arrow_every, ::arrow_every], color="r")
ax_summer.add_feature(cfeature.COASTLINE)
ax_summer.set_title("JFM")

# Winter - jas
ax_winter: plt.Axes = fig.add_subplot(1, 3, 3, projection=m)
ax_winter.set_extent(bounds, ccrs.PlateCarree())
s = ax_winter.contourf(GPathfinder.xpts, GPathfinder.ypts, conc_avg_jas, np.linspace(0, 1, 10 + 1),
                       cmap=ice_cmap)
ax_winter.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                 GPathfinder.ypts[::arrow_every, ::arrow_every],
                 ur_avg_jas[::arrow_every, ::arrow_every],
                 vr_avg_jas[::arrow_every, ::arrow_every], color="r")
ax_winter.add_feature(cfeature.COASTLINE)
ax_winter.set_title("JAS")

fig.colorbar(s, ax=[ax, ax_summer, ax_winter], orientation="horizontal",
             pad=0.05, aspect=25)

fig.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)

fig.savefig(f"Maps_output/{hemisphere}/Conc_and_wind_seasons.png")
