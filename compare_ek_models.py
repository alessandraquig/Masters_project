# %%
import numpy as np
import matplotlib
from matplotlib.colors import Normalize

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from sys import path

import parameters as par
path.insert(0, par.path)
import grid_set as gs
import data_classes as dc
import plotting_functions as pf
plt.rcParams['figure.dpi'] = 400

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
ek_all = np.load(f"Data_arrays/{hemisphere}/model1/ekman_2011-2020.npy")
ek_geo = np.load(f"Data_arrays/{hemisphere}/model6/ekman_2011-2020.npy")
ek_wind_ice = np.load(f"Data_arrays/{hemisphere}/model7/ekman_1979-2020.npy")
geo = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")

# ek_all = pf.mask_data(ek_all)
# ek_geo = pf.mask_data(ek_geo)
# ek_wind_ice = pf.mask_data(ek_wind_ice)
# geo = pf.mask_data(geo)

# Getting averages over 2011-2020
ek_all_avg = pf.mask_data(pf.get_average(ek_all[..., 0], 2011, years=years))
ek_geo_avg = pf.mask_data(pf.get_average(ek_geo[..., 0], 2011, years=years))
ek_wind_ice_avg = pf.mask_data(pf.get_average(ek_wind_ice[..., 0], 1979, years=years))
geo_avg = pf.mask_data(pf.get_average(geo, 1979, years=years))


jfm = ["jan", "feb", "mar"]
amj = ["apr", "may", "jun"]
jas = ["jul", "aug", "sep"]
ond = ["oct", "nov", "dec"]
seasons = [jfm, amj, jas, ond]

fig_ek: plt.Figure = plt.figure(figsize=(10, 5))
fig_ek.suptitle("Surface Currents 2011-2020 Mean", fontsize=15)

#
# TODO: give them all the same norm
v_max = max(np.max(ek_all_avg), np.max(geo_avg), np.max(ek_wind_ice_avg))
norm = Normalize(vmin=0, vmax=v_max)
ax_total, img_total = pf.vectorplot(ek_all_avg, fig_ek, "Total", 1, 3, 1, cmap="viridis", norm=norm)
ax_wind_ice, _ = pf.vectorplot(ek_wind_ice_avg, fig_ek, "Wind and Ice", 1, 3, 2, cmap="viridis", norm=norm)
ax_geo, _ = pf.vectorplot(geo_avg, fig_ek, "Geostrophic", 1, 3, 3, cmap="viridis", norm=norm)

fig_ek.colorbar(img_total, ax=[ax_total, ax_wind_ice, ax_geo], orientation="horizontal",
                pad=0.05, aspect=25)

fig_ek.subplots_adjust(left=0.01, right=0.99, wspace=0.05, top=0.85, bottom=0.25)


#Are all of the subplots using the same colorbar, or have I only printed it for one?

fig_ek.savefig(f"Maps_output/{hemisphere}/Surface_Ek_Geo_2011-2020.png")
#fig_pump.savefig(f"Maps_output/{hemisphere}/{model}/Pump_Geo_2011-2020.png")