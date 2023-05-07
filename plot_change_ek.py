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
ekman = np.load(f"Data_arrays/{hemisphere}/model7/ekman_1979-2020.npy")
# Assumes I'm using model 7, which I am bc it's the most complete model pre-2011

# Plot change in currents (Average for 2010s - average for 1980s)
ek_80s = pf.get_average(ekman[..., 0], 1979, np.arange(1980, 1989 + 1))
mag_80s = np.hypot(ek_80s[..., 0], ek_80s[..., 1])
ek_10s = pf.get_average(ekman[..., 0], 1979, np.arange(2010, 2019 + 1))
mag_10s = np.hypot(ek_10s[..., 0], ek_10s[..., 1])
mag_change = pf.mask_data(mag_10s - mag_80s)
ek_change = pf.mask_data(ek_10s - ek_80s)
ur_avg = ek_change[..., 0]
vr_avg = ek_change[..., 1]

# plot the magnitude of the averaged Ekman current field
ek_norm = CenteredNorm(halfrange=0.04, vcenter=0)
ek_levels = np.linspace(-0.03, 0.03, num=9)
fig_ek = plt.figure(figsize=(18, 9))
fig_ek.suptitle("Change in Average Currents 1980s-2010s", fontsize=15)

arrow_every = 10
ax_ek: plt.Axes = fig_ek.add_subplot(1, 3, 1, projection=par.m)
ax_ek.set_extent(bounds, ccrs.PlateCarree())
img_ek = ax_ek.contourf(GPathfinder.xpts, GPathfinder.ypts, mag_change,
                        cmap="PiYG", norm=ek_norm, levels=ek_levels)
ax_ek.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
          GPathfinder.ypts[::arrow_every, ::arrow_every],
          ur_avg[::arrow_every, ::arrow_every],
          vr_avg[::arrow_every, ::arrow_every], color="k")
ax_ek.add_feature(cfeature.COASTLINE)
ax_ek.set_title("Whole Year Average")

print(f"Max change = {np.nanmax(ek_change)}, min change = {np.nanmin(ek_change)}")

# September sea ice (minimum for Arctic)
# Plot change in currents (Average for 2010s - average for 1980s)
ek_80s_sep = pf.get_average(ekman[..., 0], 1979, np.arange(1980, 1989 + 1), months=["sep"])
mag_80s_sep = np.hypot(ek_80s_sep[..., 0], ek_80s_sep[..., 1])
ek_10s_sep = pf.get_average(ekman[..., 0], 1979, np.arange(2010, 2019 + 1), months=["sep"])
mag_10s_sep = np.hypot(ek_10s_sep[..., 0], ek_10s_sep[..., 1])
mag_change_sep = pf.mask_data(mag_10s_sep - mag_80s_sep)
ek_change_sep = pf.mask_data(ek_10s_sep - ek_80s_sep)
ur_avg_sep = ek_change_sep[..., 0]
vr_avg_sep = ek_change_sep[..., 1]

ax_ek_sep: plt.Axes = fig_ek.add_subplot(1, 3, 2, projection=par.m)
ax_ek_sep.set_extent(bounds, ccrs.PlateCarree())
img_ek_sep = ax_ek_sep.contourf(GPathfinder.xpts, GPathfinder.ypts, mag_change,
                                cmap="PiYG", norm=ek_norm, levels=ek_levels)
ax_ek_sep.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
             GPathfinder.ypts[::arrow_every, ::arrow_every],
             ur_avg_sep[::arrow_every, ::arrow_every],
             vr_avg_sep[::arrow_every, ::arrow_every], color="k")
ax_ek_sep.add_feature(cfeature.COASTLINE)
ax_ek_sep.set_title("September")

print(f"Max change = {np.nanmax(ek_change_sep)}, min change = {np.nanmin(ek_change_sep)}")

# March sea ice (maximum for Arctic)
# Plot change in currents (Average for 2010s - average for 1980s)
ek_80s_mar = pf.get_average(ekman[..., 0], 1979, np.arange(1980, 1989 + 1), months=["mar"])
mag_80s_mar = np.hypot(ek_80s_mar[..., 0], ek_80s_mar[..., 1])
ek_10s_mar = pf.get_average(ekman[..., 0], 1979, np.arange(2010, 2019 + 1), months=["mar"])
mag_10s_mar = np.hypot(ek_10s_mar[..., 0], ek_10s_mar[..., 1])
mag_change_mar = pf.mask_data(mag_10s_mar - mag_80s_mar)
ek_change_mar = pf.mask_data(ek_10s_mar - ek_80s_mar)
ur_avg_mar = ek_change_mar[..., 0]
vr_avg_mar = ek_change_mar[..., 1]

ax_ek_mar: plt.Axes = fig_ek.add_subplot(1, 3, 3, projection=par.m)
ax_ek_mar.set_extent(bounds, ccrs.PlateCarree())
img_ek_mar = ax_ek_mar.contourf(GPathfinder.xpts, GPathfinder.ypts, mag_change,
                                cmap="PiYG", norm=ek_norm, levels=ek_levels)
ax_ek_mar.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                 GPathfinder.ypts[::arrow_every, ::arrow_every],
                 ur_avg_mar[::arrow_every, ::arrow_every],
                 vr_avg_mar[::arrow_every, ::arrow_every], color="k")
ax_ek_mar.add_feature(cfeature.COASTLINE)
ax_ek_mar.set_title("March")

print(f"Max change = {np.nanmax(ek_change_mar)}, min change = {np.nanmin(ek_change_mar)}")

fig_ek.colorbar(img_ek, ax=[ax_ek, ax_ek_sep, ax_ek_mar], orientation="horizontal",
                pad=0.05, aspect=25, norm=ek_norm)
fig_ek.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)
fig_ek.savefig(f"Maps_output/{hemisphere}/{model}/Currents_change_seasons.png")

ek_comp = ek_change_mar - ek_change_sep
mag_comp = mag_change_mar - mag_change_sep
u_comp = ur_avg_mar - ur_avg_sep
v_comp = vr_avg_mar - vr_avg_sep

comp_norm = CenteredNorm(halfrange=0.01, vcenter=0)
fig_comp = plt.figure(figsize=(18, 9))
fig_comp.suptitle("Comparing March and September")
ax_comp: plt.Axes = fig_comp.add_subplot(1, 1, 1, projection=par.m)
ax_comp.set_extent(bounds, ccrs.PlateCarree())
img_comp = ax_comp.contourf(GPathfinder.xpts, GPathfinder.ypts, mag_comp, cmap="PiYG")
ax_comp.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
                 GPathfinder.ypts[::arrow_every, ::arrow_every],
                 u_comp[::arrow_every, ::arrow_every],
                 v_comp[::arrow_every, ::arrow_every], color="k")
ax_comp.add_feature(cfeature.COASTLINE)
ax_comp.set_title("March")

fig_comp.colorbar(img_comp, ax=[ax_comp], orientation="vertical",
                pad=0.05, aspect=25, norm=comp_norm)
fig_comp.savefig(f"Maps_output/{hemisphere}/{model}/Currents_change_between_seasons.png")

print(f"All close is {np.allclose(ek_change_mar, ek_change_sep, equal_nan=True)}")