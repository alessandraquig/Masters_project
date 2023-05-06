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

##### Loading arrays from saved files
conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")

# Find change in wind (Average for 2010s - average for 1980s)
wind_80s = pf.get_average(wind, 1979, np.arange(1980, 1989 + 1))
wind_10s = pf.get_average(wind, 1979, np.arange(2010, 2019 + 1))
mag_80s = np.hypot(wind_80s[..., 0], wind_80s[..., 1])
mag_10s = np.hypot(wind_10s[..., 0], wind_10s[..., 1])
mag_change = pf.mask_data(mag_10s - mag_80s)
wind_change = pf.mask_data(wind_10s - wind_80s)
ur_avg = wind_change[..., 0]
vr_avg = wind_change[..., 1]

wind_80s_sep = pf.get_average(wind, 1979, np.arange(1980, 1989 + 1), months=["sep"])
wind_10s_sep = pf.get_average(wind, 1979, np.arange(2010, 2019 + 1), months=["sep"])
mag_80s_sep = np.hypot(wind_80s_sep[..., 0], wind_80s_sep[..., 1])
mag_10s_sep = np.hypot(wind_10s_sep[..., 0], wind_10s_sep[..., 1])
mag_change_sep = pf.mask_data(mag_10s_sep - mag_80s_sep)
wind_change_sep = pf.mask_data(wind_10s_sep - wind_80s_sep)
ur_avg = wind_change_sep[..., 0]
vr_avg = wind_change_sep[..., 1]

wind_80s_mar = pf.get_average(wind, 1979, np.arange(1980, 1989 + 1), months=["mar"])
wind_10s_mar = pf.get_average(wind, 1979, np.arange(2010, 2019 + 1), months=["mar"])
mag_80s_mar = np.hypot(wind_80s_mar[..., 0], wind_80s_mar[..., 1])
mag_10s_mar = np.hypot(wind_10s_mar[..., 0], wind_10s_mar[..., 1])
mag_change_mar = pf.mask_data(mag_10s_mar - mag_80s_mar)
wind_change_mar = pf.mask_data(wind_10s_mar - wind_80s_mar)
ur_avg = wind_change_mar[..., 0]
vr_avg = wind_change_mar[..., 1]

#Find changes in conc
conc_80s = pf.get_average(conc, 1979, np.arange(1980, 1989 + 1))
conc_10s = pf.get_average(conc, 1979, np.arange(2010, 2019 + 1))
conc_change = conc_10s - conc_80s
conc_change = pf.mask_data(conc_change)

conc_80s_sep = pf.get_average(conc, 1979, np.arange(1980, 1989 + 1), months=["sep"])
conc_10s_sep = pf.get_average(conc, 1979, np.arange(2010, 2019 + 1), months=["sep"])
conc_change_sep = conc_10s_sep - conc_80s_sep
conc_change_sep = pf.mask_data(conc_change_sep)

conc_80s_mar = pf.get_average(conc, 1979, np.arange(1980, 1989 + 1), months=["mar"])
conc_10s_mar = pf.get_average(conc, 1979, np.arange(2010, 2019 + 1), months=["mar"])
conc_change_mar = conc_10s_mar - conc_80s_mar
conc_change_mar = pf.mask_data(conc_change_mar)

conc_norm = CenteredNorm(halfrange=5, vcenter=0)
wind_norm = CenteredNorm(halfrange=25, vcenter=0)
arrow_every = 10
fig = plt.figure(figsize=(18, 9))
fig.suptitle("Change in Average Wind and Ice Conc 1980s-2010s")

ax_year, img_year = pf.colorplot(data=conc_change, fig=fig, rows=1, cols=3, pos=1, title="Whole Year", norm=conc_norm)
ax_year.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
             GPathfinder.ypts[::arrow_every, ::arrow_every],
             ur_avg[::arrow_every, ::arrow_every],
             vr_avg[::arrow_every, ::arrow_every], color="k")

ax_sep, _ = pf.colorplot(data=conc_change_sep, fig=fig, rows=1, cols=3, pos=2, title="September", norm=conc_norm)
ax_sep.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
               GPathfinder.ypts[::arrow_every, ::arrow_every],
               ur_avg[::arrow_every, ::arrow_every],
               vr_avg[::arrow_every, ::arrow_every], color="k")

ax_mar, _ = pf.colorplot(data=conc_change_mar, fig=fig, rows=1, cols=3, pos=3, title="March", norm=conc_norm)
ax_mar.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
               GPathfinder.ypts[::arrow_every, ::arrow_every],
               ur_avg[::arrow_every, ::arrow_every],
               vr_avg[::arrow_every, ::arrow_every], color="k")

fig.colorbar(img_year, ax=[ax_year, ax_sep, ax_mar], orientation="horizontal",
                  pad=0.05, aspect=25, norm=conc_norm)
fig.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)

fig.savefig(f"Maps_output/{hemisphere}/{model}/Wind_conc_change_seasons.png")