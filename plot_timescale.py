# importing packages
import numpy as np
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits import axisartist
import parameters as par
from sys import path

path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project')

# Setting constants
hemisphere = par.HEMI
years = par.YEARS
model = par.MODEL
months = par.MONTHS

# Loading arrays
drift = np.load(f"Data_arrays/{hemisphere}/drift_1979-2020.npy")
conc = np.load(f"Data_arrays/{hemisphere}/conc_1979-2020.npy")
geo = np.load(f"Data_arrays/{hemisphere}/geo_1979-2020.npy")
wind = np.load(f"Data_arrays/{hemisphere}/wind_1979-2020.npy")
ekman = np.load(f"Data_arrays/{hemisphere}/{model}/ekman_{years[0]}-{years[-1]}.npy")
# tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_{years[0]}-{years[-1]}.npy")
pump = np.load(f"Data_arrays/{hemisphere}/{model}/pump_{years[0]}-{years[-1]}.npy")

# Getting monthly averages across the whole map
# Drift, conc, geo, and wind are all loaded for 1979-2020 (geo is nan before 2011), so years=42

# Take average ice conc across the whole map
# Conc shape is (42, 12, 361, 361) for Northern hemisphere
# Since both conc arrays (north and south) go from 1979-2020, number of years = 42
conc_avg = np.zeros((42, 12))  # 42 years (1979-2020 inclusive), 12 months
for year in range(42):
    for month in range(12):
        conc_avg[year, month] = np.nanmean(conc[year, month, :, :], axis=(0, 1))

# Take average ice drift across the whole map
ur_avg = np.zeros((42, 12))
vr_avg = np.zeros((42, 12))
drift_avg = np.zeros((42, 12))  # 42 years (1979-2020 inclusive), 12 months
for year in range(42):
    for month in range(12):
        ur_avg[year, month] = np.nanmean(drift[year, month, ..., 0], axis=(0, 1))
        vr_avg[year, month] = np.nanmean(drift[year, month, ..., 1], axis=(0, 1))
        # calculate the magnitude of the averaged Ekman current field using np.hypot
        drift_avg[year, month] = np.hypot(ur_avg[year, month], vr_avg[year, month])

# Take average geostrophic across the whole map
ur_avg = np.zeros((42, 12))
vr_avg = np.zeros((42, 12))
geo_avg = np.zeros((42, 12))  # 42 years (1979-2020 inclusive), 12 months
for year in range(42):
    for month in range(12):
        ur_avg[year, month] = np.nanmean(geo[year, month, ..., 0], axis=(0, 1))
        vr_avg[year, month] = np.nanmean(geo[year, month, ..., 1], axis=(0, 1))
        # calculate the magnitude of the averaged Ekman current field using np.hypot
        geo_avg[year, month] = np.hypot(ur_avg[year, month], vr_avg[year, month])

# Take average wind across the whole map
ur_avg = np.zeros((42, 12))
vr_avg = np.zeros((42, 12))
wind_avg = np.zeros((42, 12))  # 42 years (1979-2020 inclusive), 12 months
for year in range(42):
    for month in range(12):
        ur_avg[year, month] = np.nanmean(wind[year, month, ..., 0], axis=(0, 1))
        vr_avg[year, month] = np.nanmean(wind[year, month, ..., 1], axis=(0, 1))
        # calculate the magnitude of the averaged Ekman current field using np.hypot
        wind_avg[year, month] = np.hypot(ur_avg[year, month], vr_avg[year, month])

# Take average magnitude of Ekman currents
# Ekman shape is (10, 12, 361, 361, 2, 2) years, months, hem shape, u/v, ek/ekg0
ur_avg = np.zeros((len(years), len(months)))
vr_avg = np.zeros((len(years), len(months)))
ek_avg = np.zeros((len(years), len(months)))
for year in range(len(years)):
    for month in range(len(months)):
        ur_avg[year, month] = np.nanmean(ekman[year, month, ..., 0, 0], axis=(0, 1))
        vr_avg[year, month] = np.nanmean(ekman[year, month, ..., 1, 0], axis=(0, 1))
        # calculate the magnitude of the averaged Ekman current field using np.hypot
        ek_avg[year, month] = np.hypot(ur_avg[year, month], vr_avg[year, month])

# Take average magnitude of Ekman pumping
# Pump shape is (10, 12, 361, 361, 4)
pump_avg = np.zeros((len(years), len(months)))
for year in range(len(years)):
    for month in range(len(months)):
        pump_avg[year, month] = np.nanmean(pump[year, month, :, :, 3], axis=(0, 1))

# Flatten arrays so they can be plotted linearly against time
drift_flat = drift_avg.flatten()
conc_flat = conc_avg.flatten()
geo_flat = geo_avg.flatten()
wind_flat = wind_avg.flatten()
ek_flat = ek_avg.flatten()
pump_flat = pump_avg.flatten()


def smooth_timeseries(timeseries: np.ndarray, scale, window_type: str = "gaussian"):
    """
    Smooth input ``timeseries`` by a window with length ``scale``, using convolution.

    :param timeseries: 1D input to be smoothed
    :param scale: Length/stdev of the smoothing window.
    :param window_type: Type of window. Options: gaussian, box
    :return:
    """
    if window_type == "box":
        window = np.ones(scale) / scale  # tophat function
        output = np.convolve(timeseries, window, mode="same")
    elif window_type == "gaussian":
        output = gaussian_filter1d(timeseries, sigma=scale, mode="nearest")
    else:
        raise ValueError(f"Unrecognised window {window_type}")
    return output


window_size = 3 / 2  # standard deviation of window in months
drift_flat = smooth_timeseries(drift_flat, window_size)
conc_flat = smooth_timeseries(conc_flat, window_size)
geo_flat = smooth_timeseries(geo_flat, window_size)
wind_flat = smooth_timeseries(wind_flat, window_size)
ek_flat = smooth_timeseries(ek_flat, window_size)
pump_flat = smooth_timeseries(pump_flat, window_size)

# Getting all flattened arrays to the same length
if len(years) == 10:  # All models go from either 1979 or 2011, so this is if it's one from 2011
    nan_years = np.zeros(32 * 12)
    ek_flat = np.concatenate([nan_years, ek_flat])
    pump_flat = np.concatenate([nan_years, pump_flat])
# Plot timescale

total_months = np.arange(42 * 12)
f, time_ax = plt.subplots(1, 1, figsize=(18, 4), dpi=200)
time_ax: plt.Axes

# wind_ax: plt.Axes = time_ax.twinx()
#
# years_ad = total_months / 12 + 1979
# p1, = time_ax.plot(years_ad[-10 * 12:], ek_flat[-10 * 12:], label="Ekman Currents", c='C0')
# p2, = time_ax.plot(years_ad, drift_flat, label="Drift Speed", c='C1')
# p3, = wind_ax.plot(years_ad, wind_flat, label="Wind Speed", c='C2')
# p4, = time_ax.plot(years_ad[-10 * 12:], geo_flat[-10 * 12:], label="Geostrophic Currents", c='C3')

conc_ax: plt.Axes = time_ax.twinx()

years_ad = total_months / 12 + 1979
wind_factor = 25
# p1, = time_ax.plot(years_ad[-10 * 12:], ek_flat[-10 * 12:], label="Ekman Currents", c='C0')
p2, = time_ax.plot(years_ad, drift_flat * wind_factor,
                   label=fr"Drift Speed $\times {wind_factor}$", c='C1')
p3, = time_ax.plot(years_ad, wind_flat, label=r"Wind Speed", c='C2')
# p4, = time_ax.plot(years_ad[-10 * 12:], geo_flat[-10 * 12:], label="Geostrophic Currents", c='C3')
pConc = conc_ax.fill_between(years_ad, 100 * conc_flat,
                             color='tab:cyan', alpha=0.3, label="Ice Concentration")

time_ax.set_title("Wind Speed, Ice Drift Speed, and Ice Concentration Over Time")

time_ax.xaxis.set_major_locator(MultipleLocator(5))
time_ax.xaxis.set_minor_locator(MultipleLocator(1))
time_ax.set_xlim(1979, 2021)

time_ax.set_ylabel("Wind / Ice Speed ($m/s$)")
time_ax.set_ylim(bottom=0)
# wind_ax.set_ylabel("Wind speed ($m/s$)")
# wind_ax.yaxis.label.set_color(p3.get_color())
# wind_ax.set_ylim(bottom=0)

conc_ax.set_ylabel("Ice concentration (%)")
# wind_ax.yaxis.label.set_color(p3.get_color())
conc_ax.set_ylim(bottom=0)

# time_ax.legend(handles=[p1, p2, p3, p4])
time_ax.legend(handles=[p2, p3, pConc])

plt.subplots_adjust(left=0.05, right=0.95)

plt.savefig(f"Maps_output/{hemisphere}/sample_conc_timeseries.png")



# f.savefig('IceGovernor2020.pdf',dpu = 250)


# %%
