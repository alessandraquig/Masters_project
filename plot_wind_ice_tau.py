import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')
matplotlib.rcParams['figure.dpi'] = 400
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)

import plotting_functions as pf
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
taua_avg = pf.mask_data(pf.get_average(tau[...,1], 1979, plot_years))

# March
conc_avg_mar = pf.get_average(conc, 1979, plot_years, months=["mar"])
tau_avg_mar = pf.get_average(tau[...,0], 1979, plot_years, months=["mar"])
taui_avg_mar = pf.get_average(tau[...,2], 1979, plot_years, months=["mar"])
taua_avg_mar = pf.mask_data(pf.get_average(tau[...,1], 1979, plot_years, months=["mar"]))

#September
conc_avg_sep = pf.get_average(conc, 1979, plot_years, months=["sep"])
tau_avg_sep = pf.get_average(tau[...,0], 1979, plot_years, months=["sep"])
taui_avg_sep = pf.get_average(tau[...,2], 1979, plot_years, months=["sep"])
taua_avg_sep = pf.mask_data(pf.get_average(tau[...,1], 1979, plot_years, months=["sep"]))

#Plotting
fig_stress: plt.Figure = plt.figure(figsize=(10, 5))
fig_stress.suptitle("Wind and Ice Stress over Time", fontsize=15)
arrow_every = 10

# Whole year
ax_conc, img_conc = pf.colorplot(conc_avg, fig_stress, "Total", 1, 3, 1, cmap="Blues")
taua_u, taua_v = taua_avg[..., 0], taua_avg[..., 1]
ax_taua = ax_conc.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
          GPathfinder.ypts[::arrow_every, ::arrow_every],
          taua_u[::arrow_every, ::arrow_every],
          taua_v[::arrow_every, ::arrow_every], color="r", label="Wind Stress (Pa)")

taui_u, taui_v = taui_avg[..., 0], taui_avg[..., 1]
ax_taui = ax_conc.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
          GPathfinder.ypts[::arrow_every, ::arrow_every],
          taui_u[::arrow_every, ::arrow_every],
          taui_v[::arrow_every, ::arrow_every], color="k", label="Ice Stress (Pa)")

# March
ax_conc_mar, img_conc_mar = pf.colorplot(conc_avg_mar, fig_stress, "March", 1, 3, 2, cmap="Blues")
taua_u_mar, taua_v_mar = taua_avg_mar[..., 0], taua_avg_mar[..., 1]
ax_taua_mar = ax_conc_mar.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                         GPathfinder.ypts[::arrow_every, ::arrow_every],
                         taua_u_mar[::arrow_every, ::arrow_every],
                         taua_v_mar[::arrow_every, ::arrow_every], color="r")

taui_u_mar, taui_v_mar = taui_avg_mar[..., 0], taui_avg_mar[..., 1]
ax_taui_mar = ax_conc_mar.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                         GPathfinder.ypts[::arrow_every, ::arrow_every],
                         taui_u_mar[::arrow_every, ::arrow_every],
                         taui_v_mar[::arrow_every, ::arrow_every], color="k")

# September
ax_conc_sep, img_conc_sep = pf.colorplot(conc_avg_sep, fig_stress, "September", 1, 3, 3, cmap="Blues")
taua_u_sep, taua_v_sep = taua_avg_sep[..., 0], taua_avg_sep[..., 1]
ax_taua_sep = ax_conc_sep.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                                 GPathfinder.ypts[::arrow_every, ::arrow_every],
                                 taua_u_sep[::arrow_every, ::arrow_every],
                                 taua_v_sep[::arrow_every, ::arrow_every], color="r")

taui_u_sep, taui_v_sep = taui_avg_sep[..., 0], taui_avg_sep[..., 1]
ax_taui_sep = ax_conc_sep.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                                 GPathfinder.ypts[::arrow_every, ::arrow_every],
                                 taui_u_sep[::arrow_every, ::arrow_every],
                                 taui_v_sep[::arrow_every, ::arrow_every], color="k")


cbar = fig_stress.colorbar(img_conc, ax=[ax_conc, ax_conc_mar, ax_conc_sep], orientation="horizontal",
                pad=0.05, aspect=25)
fig_stress.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)
ax_conc.legend(handles=[ax_taua, ax_taui], loc="upper left")
cbar.set_label("Ice Concentration")

fig_stress.savefig(f"Maps_output/{hemisphere}/Ice_and_wind_stress_seasons.png")
