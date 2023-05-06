import cartopy.crs as ccrs
import numpy as np
from pathlib import Path

# path = '/home/aq/Masters_project/'
path = str(Path("../Masters_project/").resolve()) + "/"
print(path)

HEMI = "north"
YEARS = np.arange(1979, 2020 + 1)
MONTHS = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
MODEL = "model7"

# Explanation of models:
# model1: All variables included, turning angle = 45 degrees
# 2011 and later
# model2: All variables included, turning angle = 90 degrees
# 2011 and later
# model3:Wind only: no ice or geostrophic circulation
# 1979 and later
# model4: Ice-free: ice concentration = 0
# 2011 and later
# model5: Wind-free: ice and ocean
# 2011 and later
# model6: Ocean only, no ice or wind (just geostrophic balance)
# 2011 and later
# model7: Wind and ice (most complete model pre-2011 because no geostrophic data)
# 1979 and later

# We have no dynamic ocean topography data until 2011, so any models calling geostrophic data
# won't work before that year. We have other data dating back to 1979, so models that
# don't factor in geostrophy will work all the way back to then.

# Models get called in Calculating_pumping.py and Ekman_pumping.py to alter the calculations
# Models get called in Create_arrays_year_loop.ipynb and Load_arrays_plot.ipynb to organize files
# Each model has a folder (within each hemisphere folder) in Data_arrays and Maps_output

# Change these based on papers

D_E = 20.0  # Ekman depth in meters

m_AIDJEX = [1.25e-3,  # air-ocean drag coefficient
            5.5e-3,  # ice-ocean drag coefficient
            0.0,  # angle between wind velocity and surface stress
            np.deg2rad(24)]  # angle between ice velocity and surface stress, given in degrees and coverted to radians

m_AIDJEX_0 = [1.25e-3,
              5.5e-3,
              0.0,
              0.0]

# constants
rho_a = 1.25  # Air density
rho_i = 917.0  # Ice density
rho_o = 1026.0  # Ocean density

if HEMI == "north":
    geo_bounds = [-180, 180, 65, 90]
    m = ccrs.NorthPolarStereo()
    IC_grid = path + 'NSIDC_gs.npz'
    ID_grid = path + 'Pathfinder_gs.npz'
    GC_grid = path + 'PS_20km_gs2021.npz'
    hem_shape = (361, 361)
    ic_shape = (448, 304)
    gc_shape = (334, 334)
if HEMI == "south":
    geo_bounds = [-180, 180, -90, -55]
    m = ccrs.SouthPolarStereo()
    IC_grid = path + 'NSIDC_gs_SH.npz'
    ID_grid = path + 'Pathfinder_gs_SH.npz'
    GC_grid = path + 'Polar_stereo_50km_SH.npz'
    hem_shape = (321, 321)
    ic_shape = (332, 316)
    gc_shape = (212, 212)

# %%
