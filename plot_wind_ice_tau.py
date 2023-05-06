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
years = np.arange(1979, 2020+1)
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
tau = np.load(f"Data_arrays/{hemisphere}/model7/tau_1979-2020.npy")

# Get averages of wind and conc (over 2011 to 2020)
plot_years = np.arange(2011, 2020+1)
wind_avg = pf.get_average(wind, 1979, plot_years)
conc_avg = pf.get_average(conc, 1979, plot_years)
tau_avg = pf.get_average(tau[...,0], 1979, plot_years)
taui_avg = pf.get_average(tau[...,2], 1979, plot_years)
taua_avg = pf.get_average(tau[...,1], 1979, plot_years)

# Winter
wind_avg_jas = pf.get_average(wind, 1979, years, ["jul", "aug", "sep"])
conc_avg_jas = pf.get_average(conc, 1979, years, ["jul", "aug", "sep"])

# Summer
wind_avg_jfm = pf.get_average(wind, 1979, years, ["jan", "feb", "mar"])
conc_avg_jfm = pf.get_average(conc, 1979, years, ["jan", "feb", "mar"])


#Plotting
fig_stress: plt.Figure = plt.figure(figsize=(10, 5))
fig_stress.suptitle("Wind and Ice Stress over Time", fontsize=15)

# norm=CenteredNorm(halfrange=100, vcenter=0)
arrow_every = 10
ax_conc, img_conc = pf.colorplot(conc_avg, fig_stress, "Total", 1, 3, 1, cmap="Blues")
taua_u, taua_v = taua_avg[..., 0], taua_avg[..., 1]
ax_taua = ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
          GPathfinder.ypts[::arrow_every, ::arrow_every],
          taua_u[::arrow_every, ::arrow_every],
          taua_v[::arrow_every, ::arrow_every], color="r")

taui_u, taui_v = taui_avg[..., 0], taui_avg[..., 1]
ax_taui = ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
          GPathfinder.ypts[::arrow_every, ::arrow_every],
          taui_u[::arrow_every, ::arrow_every],
          taui_v[::arrow_every, ::arrow_every], color="k")


fig_stress.colorbar

fig_stress.savefig(f"Maps_output/{hemisphere}/Ice_and_wind_stress.png")
