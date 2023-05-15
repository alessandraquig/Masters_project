# %%
import matplotlib
from matplotlib.colors import CenteredNorm
from matplotlib.colorbar import Colorbar
import matplotlib.pyplot as plt
import numpy as np
this_rc_params = {
    "text.usetex": True,
    "font.family": "roman"
}
plt.rcParams.update(this_rc_params)
plt.rcParams['figure.dpi'] = 400
matplotlib.use('Agg')
import cartopy.crs as ccrs
from sys import path

import grid_set as gs
import data_classes as dc
import plotting_functions as pf
import parameters as par

path.insert(0, par.path)

hemisphere = par.HEMI
years = par.YEARS
# models = ["model1", "model2", "model3", "model4", "model5", "model6", "model7"]
models = ["model7"]
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

# geo currents
GEO = dc.CPOM_geo(f'{par.path}CPOM_geo/')
GCPOM = gs.grid_set(m)
GCPOM.load_grid(par.GC_grid)
GCPOM2GPathfinder = gs.Gs2Gs(GCPOM, GPathfinder, vectors=True)
GCPOM2Gplot = gs.Gs2Gs(GCPOM, GPathfinder, vectors=True)

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
drift = np.load(f"Data_arrays/{hemisphere}/drift_1979-2020.npy")
conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
geo = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")
wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")

# MASTER LOOPS: Plots wind, drift, concentration, and geo for each season over 10 years.
# Seasons are defined as winter (JFM), spring (AMJ), summer (JAS), and fall (OND)
# This is to be easily comparable with Ramadhan et al., 2022
# Second master loop plots Ekman currents with and without geostrophic included
# It's a separate loop in case it crashes at some point

year = 1981
jfm = ["jan", "feb", "mar"]
amj = ["apr", "may", "jun"]
jas = ["jul", "aug", "sep"]
ond = ["oct", "nov", "dec"]
seasons = [jfm, amj, jas, ond]

fig_wind = plt.figure(figsize=(15, 20))
fig_drift = plt.figure(figsize=(15, 20))
fig_conc = plt.figure(figsize=(15, 20))
axes_wind = []
axes_drift = []
axes_conc = []

# for pos in range(1, 16 + 1):  # 3 rows by 4 columns = 12 subfigures
plot_pos = 1

fig_wind: plt.Figure
fig_wind.suptitle("Average Wind By Season and Decade", fontsize=15)
fig_drift: plt.Figure
fig_drift.suptitle("Average Ice Drift By Season and Decade", fontsize=15)
fig_conc: plt.Figure
fig_conc.suptitle("Average Ice Concentration By Season and Decade", fontsize=15)

img_wind = img_drift = img_conc = None
while year < 2021:  # 80s, 90s, 00s, 10s (the collective era of Green Day)
    start_year = year
    end_year = year + 9  # This gives the average over 10 years
    year = year + 10

    for season in seasons:
        if season == jfm:
            if hemisphere == "north":
                season_name = "Winter"
            else:
                season_name = "Summer"
        elif season == amj:
            if hemisphere == "north":
                season_name = "Spring"
            else:
                season_name = "Autumn"
        elif season == jas:
            if hemisphere == "north":
                season_name = "Summer"
            else:
                season_name = "Winter"
        elif season == ond:
            if hemisphere == "north":
                season_name = "Autumn"
            else:
                season_name = "Spring"
        else:
            raise RuntimeError("lmao this should never happen")


        # Wind
        wind_avg = pf.mask_data(pf.get_average(wind, 1979, np.arange(start_year, end_year+1), months=season))
        ax_wind, img_wind = pf.vectorplot(wind_avg, fig_wind, f"{start_year} - {end_year} {season_name}", rows=4, cols=4, pos=plot_pos,
                      cmap="RdPu")
        axes_wind.append(ax_wind)

        # Drift
        drift_avg = pf.mask_data(pf.get_average(drift, 1979, np.arange(start_year, end_year+1), months=season))
        ax_drift, img_drift = pf.vectorplot(drift_avg, fig_drift, f"{start_year} - {end_year} {season_name}", rows=4, cols=4, pos=plot_pos,
                      cmap="viridis")
        axes_drift.append(ax_drift)

        # Concentration
        conc_avg = pf.mask_data(pf.get_average(conc, 1979, np.arange(start_year, end_year+1), months=season))
        ax_conc, img_conc = pf.colorplot(conc_avg, fig_conc, f"{start_year} - {end_year} {season_name}", rows=4, cols=4, pos=plot_pos,
                     cmap="Blues")
        axes_conc.append(ax_conc)

        plot_pos = plot_pos + 1

cbar_wind: Colorbar = fig_wind.colorbar(img_wind, ax=axes_wind,
                                   orientation="horizontal", pad=0.05, aspect=25)
cbar_wind.set_label("Wind Velocity (m s$^{-1}$)")

cbar_drift: Colorbar = fig_drift.colorbar(img_drift, ax=axes_drift,
                                   orientation="horizontal", pad=0.05, aspect=25)
cbar_drift.set_label("Drift Velocity (m s$^{-1}$)")

cbar_conc: Colorbar = fig_conc.colorbar(img_conc, ax=axes_conc,
                                   orientation="horizontal", pad=0.05, aspect=25)
cbar_conc.set_label(r"Ice Concentration (\%)")

fig_geo: plt.Figure = plt.figure(figsize=(17, 7))
fig_geo.suptitle("Average Geostrophic Currents By Season 2011-2020", fontsize=15)
axes_geo = []

fig_total: plt.Figure = plt.figure(figsize=(17, 7))
fig_total.suptitle("Average Surface Currents By Season 2011-2020", fontsize=15)
axes_total = []

ekman = np.load(f"Data_arrays/{hemisphere}/model1/ekman_2011-2020.npy")
if hemisphere == "south":
    geo[..., 0] *= -1
    geo[..., 1] *= -1
total = ekman[..., 0] + geo[-10:, ...]
print(f"Shape of total = {np.shape(total)} (should be (10, 12, 361, 361, 2) for NH)")
print(f"Valid values = {np.count_nonzero(np.isfinite(total))}")

img_geo = img_total = None
for s_idx, season in enumerate(seasons):
    if season == jfm:
        if hemisphere == "north":
            season_name = "Winter"
        else:
            season_name = "Summer"
    elif season == amj:
        if hemisphere == "north":
            season_name = "Spring"
        else:
            season_name = "Autumn"
    elif season == jas:
        if hemisphere == "north":
            season_name = "Summer"
        else:
            season_name = "Winter"
    elif season == ond:
        if hemisphere == "north":
            season_name = "Autumn"
        else:
            season_name = "Spring"
    else:
        raise RuntimeError("lmao this should never happen")
    
    geo_avg = pf.mask_data(pf.get_average(geo, 1979, np.arange(2011, 2020 + 1), months=season))
    ax_geo, img_geo = pf.vectorplot(geo_avg, fig_geo, f"2011 - 2020 {season_name}", rows=1,
                                    cols=4, pos=(s_idx + 1), cmap="Greens")
    axes_geo.append(ax_geo)

    total_avg = pf.mask_data(pf.get_average(total, 2011, np.arange(2011, 2020 + 1), months=season))
    print(np.shape(total_avg))
    ax_total, img_total = pf.vectorplot(total_avg, fig_total, f"2011 - 2020 {season_name}", rows=1,
                                    cols=4, pos=(s_idx + 1), cmap="viridis")
    axes_total.append(ax_total)

fig_geo.colorbar(img_geo, ax=axes_geo, orientation="horizontal",
                 pad=0.05, aspect=25)
fig_geo.set_label("Geostrophic Current Velocity (m s$^{-1}$)")
fig_geo.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.9, bottom=0.25)

fig_total.colorbar(img_total, ax=axes_total, orientation="horizontal",
                 pad=0.05, aspect=25)
fig_total.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.9, bottom=0.25)

fig_total.savefig(f"Maps_output/{hemisphere}/Currents_Seasonal_Avg_Decades.png")
fig_geo.savefig(f"Maps_output/{hemisphere}/Geo/Geo_Seasonal_Avg_Decades.png")
fig_drift.savefig(f"Maps_output/{hemisphere}/Drift/Drift_Seasonal_Avg_Decades.png")
fig_conc.savefig(f"Maps_output/{hemisphere}/Conc/Conc_Seasonal_Avg_Decades.png")
fig_wind.savefig(f"Maps_output/{hemisphere}/Wind/Wind_Seasonal_Avg_Decades.png")

print()
print("=" * 80)
print(f"{'EKMAN':.^80s}")
print("=" * 80)
print()

pump_cmap = plt.get_cmap("seismic").copy()
pump_cmap.set_bad(color='lightgrey')
for model in models:
    print(f"{model:.^80s}")
    if model == "model3" or model == "model7":
        data_start_year = 1979
        year = 1981
        rows = 4
        ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_1979-2020.npy")
        tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_1979-2020.npy")
        pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_1979-2020.npy")
    else:
        data_start_year = 2011
        year = 2011
        rows = 1
        ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_2011-2020.npy")
        tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_2011-2020.npy")
        pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_2011-2020.npy")

    fig_ek: plt.Figure = plt.figure(figsize=(12, 5*rows))
    fig_pump: plt.Figure = plt.figure(figsize=(12, 5*rows))
    fig_ek.suptitle(f"Average Ekman Currents By Season and Decade", fontsize=15)
    fig_pump.suptitle(f"Average Ekman Pumping By Season and Decade", fontsize=15)
    axes_ek = []
    axes_pump = []

    plot_pos = 1

    while year < 2021:  # 80s, 90s, 00s, 10s
        start_year = year
        end_year = year + 9  # This gives the average over 10 years
        year = year + 10

        for season in seasons:
            if season == jfm:
                if hemisphere == "north":
                    season_name = "Winter"
                else:
                    season_name = "Summer"
            elif season == amj:
                if hemisphere == "north":
                    season_name = "Spring"
                else:
                    season_name = "Autumn"
            elif season == jas:
                if hemisphere == "north":
                    season_name = "Summer"
                else:
                    season_name = "Winter"
            elif season == ond:
                if hemisphere == "north":
                    season_name = "Autumn"
                else:
                    season_name = "Spring"

            # Ekman Currents
            ek_avg = pf.get_average(ekman[..., 0], data_start_year, np.arange(start_year, end_year+1), months=season)
            ek_avg_masked = pf.mask_data(ek_avg)
            ax_ek, img_ek = pf.vectorplot(ek_avg, fig_ek, f"{start_year} - {end_year} {season_name}", rows=rows, cols=4, pos=plot_pos,
                          cmap="viridis")
            axes_ek.append(ax_ek)

            # Ekman Pumping
            pump_avg = pf.get_average(pump[..., 3], data_start_year, np.arange(start_year, end_year+1), months=season)
            pump_avg_masked = pf.mask_data(pump_avg)
            print(f"{np.count_nonzero(np.isfinite(pump_avg_masked))/pump_avg_masked.size:.2%} unmasked")
            ax_pump, img_pump = pf.colorplot(pump_avg, fig_pump,
                                             f"{start_year} - {end_year} {season_name}",
                                             rows=rows, cols=4, pos=plot_pos, cmap=pump_cmap,
                                             norm=CenteredNorm(halfrange=75, vcenter=0))
            axes_pump.append(ax_pump)

            plot_pos = plot_pos + 1

    cbar: Colorbar = fig_ek.colorbar(img_ek, ax=axes_ek,
                                     orientation="horizontal", pad=0.05, aspect=25)
    fig_ek.savefig(f"Maps_output/{hemisphere}/{model}/Ekman_Currents_Seasonal_Avg_Decades.png")

    cbar: Colorbar = fig_pump.colorbar(img_pump, ax=axes_pump,
                                       orientation="horizontal", pad=0.05, aspect=25)
    cbar.set_label("Ekman Pumping (m yr$^{-1}$)")
    fig_pump.savefig(f"Maps_output/{hemisphere}/{model}/Ekman_Pumping_Seasonal_Avg_Decades.png")


