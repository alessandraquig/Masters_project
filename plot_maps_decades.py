# %%
import numpy as np
import matplotlib
from matplotlib.colors import CenteredNorm

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path

import grid_set as gs
import data_classes as dc
import plotting_functions as pf
import parameters as par

path.insert(0, par.path)

hemisphere = par.HEMI
years = par.YEARS
models = ["model1", "model2", "model3", "model4", "model5", "model6", "model7"]
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

        print(start_year, end_year)

        fig_wind = plt.figure(figsize=(10, 10))
        wind_avg = pf.get_average(wind, 1979, np.arange(start_year, end_year+1), months=season)
        ax_wind, img_wind = pf.vectorplot(wind_avg, fig_wind, f"Average Wind {season_name} {start_year}-{end_year}", rows=1, cols=1, pos=1,
                      cmap="PiYG", arrow_every=10)
        fig_wind.colorbar(img_wind, ax=ax_wind)
        fig_wind.savefig(f"Maps_output/{hemisphere}/Wind/Wind_{start_year}-{end_year}_{season_name}.png")
        plt.close(fig_wind)

        # Drift
        fig_drift = plt.figure(figsize=(10, 10))
        drift_avg = pf.get_average(drift, 1979, np.arange(start_year, end_year+1), months=season)
        drift_avg_masked = pf.mask_data(drift_avg)
        ax_drift, img_drift = pf.vectorplot(drift_avg_masked, fig_drift, f"Average Ice Drift {season_name} {start_year}-{end_year}", rows=1, cols=1, pos=1,
                      cmap="viridis", arrow_every=10)
        fig_drift.colorbar(img_drift, ax=ax_drift)
        fig_drift.savefig(f"Maps_output/{hemisphere}/Drift/Drift_{start_year}-{end_year}_{season_name}.png")
        plt.close(fig_drift)

        # Concentration
        fig_conc = plt.figure(figsize=(10, 10))
        conc_avg = pf.get_average(conc, 1979, np.arange(start_year, end_year+1), months=season)
        conc_avg_masked = pf.mask_data(conc_avg)
        ax_conc, img_conc = pf.colorplot(conc_avg, fig_conc, f"Average Ice Concentration {season_name} {start_year}-{end_year}", rows=1, cols=1, pos=1,
                     cmap="Blues")
        fig_conc.colorbar(img_conc, ax=ax_conc)
        fig_conc.savefig(f"Maps_output/{hemisphere}/Conc/Conc_{start_year}-{end_year}_{season_name}.png")
        plt.close(fig_conc)

        # Geostrophic
        fig_geo = plt.figure(figsize=(10, 10))
        geo_avg = pf.get_average(geo, 1979, np.arange(start_year, end_year+1), months=season)
        geo_avg_masked = pf.mask_data(geo_avg)
        ax_geo, img_geo = pf.vectorplot(geo_avg, fig_geo, f"Average Geostrophic Currents {season_name} {start_year}-{end_year}", rows=1, cols=1, pos=1,
                      cmap="PiYG", arrow_every=10)
        fig_geo.colorbar(img_geo, ax=ax_geo)
        fig_geo.savefig(f"Maps_output/{hemisphere}/Geo/Geo_{start_year}-{end_year}_{season_name}.png")
        plt.close(fig_geo)

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
            fig_ek: plt.Figure = plt.figure(figsize=(10, 10))
            ek_avg = pf.get_average(ekman[..., 0], data_start_year, np.arange(start_year, end_year+1), months=season)
            ek_avg_masked = pf.mask_data(ek_avg)
            ax_ek, img_ek = pf.vectorplot(ek_avg, fig_ek, f"Average Ekman Currents {season_name} {start_year}-{end_year}", rows=1, cols=1, pos=1,
                          cmap="viridis", arrow_every=10)
            fig_ek.colorbar(img_ek, ax=ax_ek)
            fig_ek.savefig(f"Maps_output/{hemisphere}/{model}/Ekman_Currents_{start_year}-{end_year}_{season_name}.png")
            plt.close(fig_ek)

            # Ekman Pumping
            fig_pump: plt.Figure = plt.figure(figsize=(10, 10))
            pump_avg = pf.get_average(pump[..., 3], data_start_year, np.arange(start_year, end_year+1), months=season)
            pump_avg_masked = pf.mask_data(pump_avg)
            print(f"{np.count_nonzero(np.isfinite(pump_avg_masked))/pump_avg_masked.size:.2%} unmasked")
            ax_pump, img_pump = pf.colorplot(pump_avg, fig_pump, f"Average Ekman Pumping {season_name} {start_year}-{end_year}", rows=1, cols=1, pos=1,
                         cmap=pump_cmap, norm=CenteredNorm(halfrange=250, vcenter=0))
            fig_pump.colorbar(img_pump, ax=ax_pump)
            fig_pump.savefig(f"Maps_output/{hemisphere}/{model}/Ekman_Pumping_{start_year}-{end_year}_{season_name}.png")
            plt.close(fig_pump)



