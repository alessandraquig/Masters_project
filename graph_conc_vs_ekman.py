#importing packages
import numpy as np
import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# import grid_set as gs
# import data_classes as dc
import parameters as par
from sys import path
path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project')

#Setting constants
hemisphere = par.HEMI
years = par.YEARS
model = par.MODEL
months = par.MONTHS

#Note: I haven't loaded any grids because I don't think I'll need them.
    #If this doesn't work, keep that in mind.

#Loading arrays
#drift = np.load(f"Data_arrays/{hemisphere}/drift_1979-2020.npy")
conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
#geo = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")
#wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")
ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_{years[0]}-{years[-1]}.npy")
#tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_{years[0]}-{years[-1]}.npy")
pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_{years[0]}-{years[-1]}.npy")

#Take average ice conc across the whole map
#Conc shape is (42, 12, 361, 361) for Northern hemisphere
#Since both conc arrays (north and south) go from 1979-2020, number of years = 42
conc_avg = np.zeros((42, 12)) #42 years (1979-2020 inclusive), 12 months
for year in range(42):
    for month in range(12):
        conc_avg[year,month] = np.nanmean(conc[year, month, :, :], axis=(0, 1))

#Take average magnitude of surface currents
#Ekman shape is (10, 12, 361, 361, 2, 2) years, months, hem shape, u/v, ek/ekg0
ur_avg = np.zeros((len(years), len(months)))
vr_avg = np.zeros((len(years), len(months)))
mag = np.zeros((len(years), len(months)))
for year in range(len(years)):
    for month in range(len(months)):
        ur_avg[year, month] = np.nanmean(ekman[year, month, ..., 0, 0], axis=(0, 1))
        vr_avg[year, month] = np.nanmean(ekman[year, month, ..., 1, 0], axis=(0, 1))
        # calculate the magnitude of the averaged Ekman current field using np.hypot
        mag[year, month] = np.hypot(ur_avg[year, month], vr_avg[year, month])

#Take average magnitude of Ekman pumping
#Pump shape is (10, 12, 361, 361, 4)
pump_avg = np.zeros((len(years), len(months)))
for year in range(len(years)):
    for month in range(len(months)):
        pump_avg[year,month] = np.nanmean(pump[year, month, :, :, 3], axis=(0, 1))

#Make a scatterplot with ice conc on x axis and magnitude of Ekman pumping on y axis

if model != "model3" and model != "model7":
    conc_avg = conc_avg[32:,:]

fig, ax = plt.subplots()
for m, conc, pump in zip(months, conc_avg, pump_avg):
    if m in ["sep", "oct", "nov"]:
        color = 'orange'
    elif m in ["dec", "jan", "feb"]:
        color = 'blue'
    elif m in ["mar", "apr", "may"]:
        color = 'green'
    else:
        color = 'red'
    ax.scatter(conc, pump, c=color)
    plt.xlabel("Ice Concentration Monthly Average over Arctic")
    plt.ylabel("Pumping Monthly Average over Arctic")
    plt.title("Pumping vs. Ice Concentration (Wind and Ice-Driven)")

    blue_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='blue', label='December-February')
    green_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='green', label='March-May')
    red_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='red', label='June-August')
    orange_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='orange', label='September-November')
    plt.legend(handles=[blue_patch[0], green_patch[0], red_patch[0], orange_patch[0]])

    plt.savefig(f"Maps_output/{hemisphere}/{model}/Pumping_vs_conc.png")

# #Make another scatterplot with ice conc on x and magnitude of surface currents on y
fig, ax = plt.subplots()
for m, conc, pump in zip(months, conc_avg, mag):
    if m in ["sep", "oct", "nov"]:
        color = 'orange'
    elif m in ["dec", "jan", "feb"]:
        color = 'blue'
    elif m in ["mar", "apr", "may"]:
        color = 'green'
    else:
        color = 'red'
    ax.scatter(conc, pump, c=color)
    plt.xlabel("Ice Concentration Monthly Average over Arctic")
    plt.ylabel("Surface Currents Monthly Average over Arctic")
    plt.title("Current Magnitude vs. Ice Concentration (Wind and Ice-Driven)")

    blue_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='blue', label='December-February')
    green_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='green', label='March-May')
    red_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='red', label='June-August')
    orange_patch = plt.plot([],[], marker="o", markersize=10, ls="", mec=None, color='orange', label='September-November')
    plt.legend(handles=[blue_patch[0], green_patch[0], red_patch[0], orange_patch[0]])

    plt.savefig(f"Maps_output/{hemisphere}/{model}/Surface_vs_conc.png")
