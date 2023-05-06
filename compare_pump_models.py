# %%
import numpy as np
import matplotlib
from matplotlib.colors import Normalize

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm

plt.rcParams['figure.dpi'] = 400
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path

path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project')
import grid_set as gs
import data_classes as dc
import parameters as par
import plotting_functions as pf

hemisphere = par.HEMI
years = np.arange(2011, 2020 + 1)
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
pump_all = np.load(f"Data_arrays/{hemisphere}/model1/pump_2011-2020.npy")
pump_wind_ice = np.load(f"Data_arrays/{hemisphere}/model7/pump_1979-2020.npy")
pump_wind = np.load(f"Data_arrays/{hemisphere}/model3/pump_1979-2020.npy")

# Getting averages over 2011-2020
pump_all_avg = pf.mask_data(pf.get_average(pump_all[..., 3], 2011, years=years))
pump_wind_avg = pf.mask_data(pf.get_average(pump_wind[..., 3], 1979, years=years))
pump_wind_ice_avg = pf.mask_data(pf.get_average(pump_wind_ice[..., 3], 1979, years=years))

print(pump_wind_ice_avg.shape, pump_wind_avg.shape, pump_all_avg.shape)

jfm = ["jan", "feb", "mar"]
amj = ["apr", "may", "jun"]
jas = ["jul", "aug", "sep"]
ond = ["oct", "nov", "dec"]
seasons = [jfm, amj, jas, ond]

fig_pump: plt.Figure = plt.figure(figsize=(10, 5))
fig_pump.suptitle("Ekman Pumping 2011-2020 Mean", fontsize=15)

#
# Plot
norm = CenteredNorm(halfrange=100, vcenter=0)
ax_total, img_total = pf.colorplot(pump_all_avg, fig_pump, "Total", 1, 3, 1, cmap="seismic", norm=norm)
ax_wind_ice, _ = pf.colorplot(pump_wind_ice_avg, fig_pump, "Wind and Ice", 1, 3, 2, cmap="seismic", norm=norm)
ax_wind, _ = pf.colorplot(pump_wind_avg, fig_pump, "Wind", 1, 3, 3, cmap="seismic", norm=norm)

fig_pump.colorbar(img_total, ax=[ax_total, ax_wind_ice, ax_wind], orientation="horizontal",
                  pad=0.05, aspect=25)

fig_pump.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)

# Are all of the subplots using the same colorbar, or have I only printed it for one?

fig_pump.savefig(f"Maps_output/{hemisphere}/Pump_model_comparison.png")
# fig_pump.savefig(f"Maps_output/{hemisphere}/{model}/Pump_Geo_2011-2020.png")
