import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 400
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path
import parameters as par

path.insert(0, f"{par.path}")
import grid_set as gs
import data_classes as dc
import plotting_functions as pf

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

##### Loading arrays from saved files
drift_data = np.load(f"Data_arrays/{hemisphere}/drift_1979-2020.npy")
conc_data = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
geo_data = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")
wind_data = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")

print(f"Drift shape = {drift_data.shape}")
print(f"Conc shape = {conc_data.shape}")
print(f"Geo shape = {geo_data.shape}")
print(f"Wind shape = {wind_data.shape}")

# MASTER LOOP: Plots wind, drift, concentration, and geo for each month in each year.
# Also plots Ekman surface currents and pumping for each month, year, and model

for year_idx, year in enumerate(years):
    # year_diff = years[0] - 1979  # since data arrays start in 1979 (geostrophic is empty before 2011)
    # data_year_idx = year_idx + year_diff
    # # year_idx will work for Ekman but won't always for for raw data
    print(f"Year index = {year_idx}")
    for month_idx, month in enumerate(months):
        print(f"Month index = {month_idx}")
        month_num = month_idx + 1
        if month_num < 10:
            month_num = f"0{month_num}"

        #Formatting data
        wind = pf.mask_data(pf.get_average(wind_data, 1979, years=[year], months=[month]))
        conc = pf.mask_data(pf.get_average(conc_data, 1979, years=[year], months=[month]))
        geo = pf.mask_data(pf.get_average(geo_data, 1979, years=[year], months=[month]))
        drift = pf.mask_data(pf.get_average(drift_data, 1979, years=[year], months=[month]))
        print(np.shape(wind))

        # Wind
        fig_wind = plt.figure(figsize=(10, 10))
        fig_wind.suptitle("Wind")
        ax_wind, img_wind = pf.vectorplot(wind, fig_wind,
                                          title=f"{month.capitalize()} {year}", rows=1,
                                          cols=1, pos=1, cmap="viridis", arrow_every=10, norm=None)
        fig_wind.colorbar(img_wind, ax=ax_wind)
        plt.savefig(f"Maps_output/{hemisphere}/Wind/Wind_{year}-{month_num}.png")
        plt.close()

        # Drift
        fig_drift = plt.figure(figsize=(10, 10))
        fig_drift.suptitle("Ice Drift")
        ax_drift, img_drift = pf.vectorplot(drift, fig_drift,
                                          title=f"{month.capitalize()} {year}", rows=1,
                                          cols=1, pos=1, cmap="viridis", arrow_every=10, norm=None)
        fig_drift.colorbar(img_drift, ax=ax_drift)
        plt.savefig(f"Maps_output/{hemisphere}/Drift/Drift_{year}-{month_num}.png")
        plt.close()

        # Concentration
        fig_conc = plt.figure(figsize=(10, 10))
        fig_conc.suptitle("Ice Concentration")
        ax_conc, img_conc = pf.colorplot(conc, fig_conc,
                                      title=f"{month.capitalize()} {year}", rows=1,
                                      cols=1, pos=1, cmap="seismic", norm=None)
        fig_conc.colorbar(img_conc, ax=ax_conc)
        fig_conc.savefig(f"Maps_output/{hemisphere}/Conc/Conc_{year}-{month_num}.png")
        plt.close()

        # Geostrophic
        fig_geo = plt.figure(figsize=(10, 10))
        fig_geo.suptitle("Geostrophic Currents")
        ax_geo, img_geo = pf.vectorplot(geo, fig_geo,
                                          title=f"{month.capitalize()} {year}", rows=1,
                                          cols=1, pos=1, cmap="viridis", arrow_every=10, norm=None)
        fig_geo.colorbar(img_geo, ax=ax_geo)
        plt.savefig(f"Maps_output/{hemisphere}/Geo/Geo_{year}-{month_num}.png")
        plt.close()

# Plotting Ekman currents and pumping
for year_idx, year in enumerate(years):
    print(f"Year: {year}")
    for month_idx, month in enumerate(months):
        print(f"\tmonth: {month}")
        month_num = month_idx + 1
        if month_num < 10:
            month_num = f"0{month_num}"
        for model in models:
            print(f"\t\tmodel: {model}")
            if model == "model3" or model == "model7":
                data_start = 1979
                ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_1979-2020.npy")
                tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_1979-2020.npy")
                pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_1979-2020.npy")
            else:
                data_start = 2011
                ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_2011-2020.npy")
                tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_2011-2020.npy")
                pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_2011-2020.npy")

            ekman = pf.mask_data(pf.get_average(ekman, data_start, years=[year], months=[month]))
            pump = pf.mask_data(pf.get_average(pump, data_start, years=[year], months=[month]))

            # Ekman Pumping
            fig_pump = plt.figure(figsize=(10, 10))
            fig_pump.suptitle("Ekman Pumping")
            ax_pump, img_pump = pf.colorplot(pump[..., 3], fig_pump,
                                      title=f"{month.capitalize()} {year}", rows=1,
                                      cols=1, pos=1, cmap="seismic", norm=None)
            fig_pump.colorbar(img_pump, ax=ax_pump)
            plt.savefig(f"Maps_output/{hemisphere}/{model}/Ekman_Pumping_{year}-{month_num}.png")
            plt.close()

            # Ekman Currents
            fig_ek = plt.figure(figsize=(10, 10))
            fig_ek.suptitle("ekstrophic Currents")
            ax_ek, img_ek = pf.vectorplot(ekman[..., 0], fig_ek,
                                        title=f"{month.capitalize()} {year}", rows=1,
                                        cols=1, pos=1, cmap="viridis", arrow_every=10, norm=None)
            fig_ek.colorbar(img_ek, ax=ax_ek)
            plt.savefig(f"Maps_output/{hemisphere}/{model}/Ekman_Currents_{year}-{month_num}.png")
            plt.close()
