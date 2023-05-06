# %%
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path

path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project')
import grid_set as gs
import data_classes as dc
import parameters as par
import plotting_functions as plot

hemisphere = par.HEMI
years = np.arange(2011, 2020+1)
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
#drift = np.load(f"Data_arrays/{hemisphere}/drift_1979-2020.npy")
#conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
geo = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")
#wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")

def plot_ek_annual(fig: plt.Figure, rows, cols, pos, title, year_index, months=None, arrow_every=10):
    """
   Takes a range of years and months (default months value is all) and gives the average over that grid cell for that time range

    :param rows: How many rows of subplots (for single plot, 1)
    :param cols: How many columns of subplots (for single plot, 1)
    :param pos: Position in the subplot (for single plot, 1)
    :param title: Title of subplot or plot
    :param start_year: The first year of the average
    :param end_year: Include up to and including this year
    :param arrow_every: Density of arrows (usually 10)
    """
    if months is None:
        months = par.MONTHS
    months_idxs = [par.MONTHS.index(mon) for mon in months]
    ur_avg = np.nanmean(ekman[year_index, ..., 0, 0], axis=0)
    vr_avg = np.nanmean(ekman[year_index, ..., 1, 0], axis=0)

    # calculate the magnitude of the averaged Ekman current field using np.hypot
    mag = np.hypot(ur_avg, vr_avg)

    # plot the magnitude of the averaged Ekman current field
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=m)
    ax.set_extent(bounds, ccrs.PlateCarree())
    s = ax.contourf(GPathfinder.xpts, GPathfinder.ypts, mag)
    ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
              GPathfinder.ypts[::arrow_every, ::arrow_every],
              ur_avg[::arrow_every, ::arrow_every],
              vr_avg[::arrow_every, ::arrow_every], color="k")
    ax.add_feature(cfeature.COASTLINE)
    fig.colorbar(s)
    ax.set_title(title)

def plot_geo_annual(fig: plt.Figure, rows, cols, pos, title, year_index, months=None, arrow_every=10):
    if months is None:
        months = par.MONTHS
    months_idxs = [par.MONTHS.index(mon) for mon in months]
    ur_avg = np.nanmean(geo[year_index, ..., 0], axis=0)
    vr_avg = np.nanmean(geo[year_index, ..., 1], axis=0)

    # calculate the magnitude of the averaged Ekman current field using np.hypot
    mag = np.hypot(ur_avg, vr_avg)

    # plot the magnitude of the averaged Ekman current field
    ax: plt.Axes = fig.add_subplot(rows, cols, pos, projection=m)
    ax.set_extent(bounds, ccrs.PlateCarree())
    s = ax.contourf(GPathfinder.xpts, GPathfinder.ypts, mag)
    ax.quiver(GPathfinder.xpts[::arrow_every, ::arrow_every],
              GPathfinder.ypts[::arrow_every, ::arrow_every],
              ur_avg[::arrow_every, ::arrow_every],
              vr_avg[::arrow_every, ::arrow_every], color="k")
    ax.add_feature(cfeature.COASTLINE)
    fig.colorbar(s)
    ax.set_title(title)

#I want to create a plot with a given season (or overall) for each of ten years (2011-2020)
#Column 1 is model 1, column 2 is model 7, column 3 is geostrophic
jfm = ["jan", "feb", "mar"]
amj = ["apr", "may", "jun"]
jas = ["jul", "aug", "sep"]
ond = ["oct", "nov", "dec"]
seasons = [jfm, amj, jas, ond]

fig_ek = plt.figure(figsize=(15, 20))
fig_pump = plt.figure(figsize=(15, 20))


# for decade in range(4):  # 80s, 90s, 00s, 10s (the collective era of Green Day)
#     start_year = year
#     end_year = year + 9  # This gives the average over 10 years
#     year = year + 10
#
#     for season in seasons:
#         if season == jfm:
#             if hemisphere == "north":
#                 season_name = "Winter"
#             else:
#                 season_name = "Summer"
#         elif season == amj:
#             if hemisphere == "north":
#                 season_name = "Spring"
#             else:
#                 season_name = "Autumn"
#         elif season == jas:
#             if hemisphere == "north":
#                 season_name = "Summer"
#             else:
#                 season_name = "Winter"
#         elif season == ond:
#             if hemisphere == "north":
#                 season_name = "Autumn"
#             else:
#                 season_name = "Spring"


fig_ek = plt.figure(figsize=(12, 24))
fig_ek.suptitle("Total, Ekman, and Geostrophic Currents")
fig_pump = plt.figure(figsize=(12, 24))

plot_pos = 1 #Goes up to 30 (inclusive)

# Plotting Ekman currents and pumping
for year_idx, year in enumerate(years):

    #fig_ek: plt.Figure
    year_index = year_idx
    ekman = np.load(f"Data_arrays/{hemisphere}/model1/ekman_2011-2020.npy")
    plot_ek_annual(fig_ek, 10, 3, plot_pos, f"Total {year}", year_index)
    plot_pos = plot_pos + 1

    year_diff = years[0] - 1979  # since data arrays start in 1979 (geostrophic is empty before 2011)
    year_index = year_idx + year_diff

    ekman = np.load(f"Data_arrays/{hemisphere}/model7/ekman_1979-2020.npy")
    plot_ek_annual(fig_ek, 10, 3, plot_pos, f"Ekman {year}", year_index)
    plot_pos = plot_pos + 1


    plot_geo_annual(fig_ek, 10, 3, plot_pos, f"Geostrophic {year}", year_index)
    plot_pos = plot_pos + 1


fig_ek.savefig(f"Maps_output/{hemisphere}/Surface_Ek_Geo_2011-2020.png")
#fig_pump.savefig(f"Maps_output/{hemisphere}/{model}/Pump_Geo_2011-2020.png")
