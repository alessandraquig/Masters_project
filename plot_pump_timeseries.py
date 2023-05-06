# importing packages
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
import parameters as par
from sys import path
import plotting_functions as pf

path.insert(0, par.path)

# Setting constants
hemisphere = par.HEMI
years = np.arange(2011, 2020+1)
months = par.MONTHS

# Loading arrays
pump1 = np.load(f"Data_arrays/{hemisphere}/model1/pump_2011-2020.npy")
pump7 = np.load(f"Data_arrays/{hemisphere}/model7/pump_1979-2020.npy")

print(pump1.shape, pump7.shape)

# Getting monthly averages across the whole map
# Drift, conc, geo, and wind are all loaded for 1979-2020 (geo is nan before 2011), so years=42

# Take average magnitude of Ekman pumping
# Pump shape is (10, 12, 361, 361, 4)
pump1_avg = np.zeros((len(years), len(months)))
pump7_avg = np.zeros((len(years), len(months)))
start_diff = 2011 - 1979
for year_idx, year in enumerate(years):
    for month_idx, month in enumerate(months):
        pump1_avg[year_idx, month_idx] = np.nanmean(pump1[year_idx, month_idx, ..., 3], axis=(0, 1))

        year_idx7 = year_idx + start_diff
        pump7_avg[year_idx, month_idx] = np.nanmean(pump7[year_idx7, month_idx, ..., 3], axis=(0, 1))

# Flatten arrays so they can be plotted linearly against time
pump1_flat = pump1_avg.flatten()
pump7_flat = pump7_avg.flatten()

print(np.nanmean(pump1_flat), np.nanmax(pump1_flat), np.nanmin(pump1_flat))
print(np.nanmean(pump7_flat), np.nanmax(pump7_flat), np.nanmin(pump7_flat))

window_size = 2 / 2  # standard deviation of window in months
#pump1_flat = pf.smooth_timeseries(pump1_flat, window_size)
#pump7_flat = pf.smooth_timeseries(pump7_flat, window_size)

print(np.nanmean(pump1_flat), np.nanmax(pump1_flat), np.nanmin(pump1_flat))
print(np.nanmean(pump7_flat), np.nanmax(pump7_flat), np.nanmin(pump7_flat))

# Plot timescale
total_months = np.arange(10 * 12)
years_ad = total_months / 12 + 2011

f, time_ax = plt.subplots(1, 1, figsize=(18, 4), dpi=200)
time_ax: plt.Axes


time_ax.plot(years_ad, pump1_flat, label="With Geostrophic", c='C0')
time_ax.plot(years_ad, pump7_flat, label="Without Geostrophic", c='C1')

time_ax.set_title("Ekman Pumping Over Time")

time_ax.xaxis.set_major_locator(MultipleLocator(5))
time_ax.xaxis.set_minor_locator(MultipleLocator(1))
# time_ax.set_xlim(2010, 2021)
# time_ax.set_ylim(-4,4)

time_ax.set_ylabel("Pumping ($m^{2}/yr$)")
# time_ax.set_ylim(bottom=0)

time_ax.legend()

f.subplots_adjust(left=0.05, right=0.95)

f.savefig(f"Maps_output/{hemisphere}/pump_timeseries.png")

# f.savefig('IceGovernor2020.pdf',dpu = 250)
