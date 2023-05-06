# importing packages
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import parameters as par
from sys import path
import plotting_functions as pf

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
tau = np.load(f"Data_arrays/{hemisphere}/{model}/tau_{years[0]}-{years[-1]}.npy")
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

# Take average magnitude of tau (air and ice)
# Tau shape for model 1 North is (10, 12, 361, 361, 2, 5) years, months, hem shape, u/v, tau_a/tau_i/tau_i0/tau_g/tau_all
tau_a_u_avg = np.zeros((len(years), len(months)))
tau_a_v_avg = np.zeros((len(years), len(months)))
tau_i_u_avg = np.zeros((len(years), len(months)))
tau_i_v_avg = np.zeros((len(years), len(months)))
tau_a_avg = np.zeros((len(years), len(months)))
tau_i_avg = np.zeros((len(years), len(months)))
for year in range(len(years)):
    for month in range(len(months)):
        tau_a_u_avg[year, month] = np.nanmean(tau[year, month, ..., 0, 0], axis=(0, 1))
        tau_a_v_avg[year, month] = np.nanmean(tau[year, month, ..., 1, 0], axis=(0, 1))
        tau_i_u_avg[year, month] = np.nanmean(tau[year, month, ..., 0, 1], axis=(0, 1))
        tau_i_v_avg[year, month] = np.nanmean(tau[year, month, ..., 1, 1], axis=(0, 1))

        # calculate the magnitude of the averaged Ekman current field using np.hypot
        tau_a_avg[year, month] = np.hypot(tau_a_u_avg[year, month], tau_a_v_avg[year, month])
        tau_i_avg[year, month] = np.hypot(tau_i_u_avg[year, month], tau_i_v_avg[year, month])

# Flatten arrays so they can be plotted linearly against time
conc_flat = conc_avg.flatten()
tau_a_flat = tau_a_avg.flatten()
tau_i_flat = tau_i_avg.flatten()


window_size = 3 / 2  # standard deviation of window in months
conc_flat = pf.smooth_timeseries(conc_flat, window_size)
tau_a_flat = pf.smooth_timeseries(tau_a_flat, window_size)
tau_i_flat = pf.smooth_timeseries(tau_i_flat, window_size)

# Getting all flattened arrays to the same length
#if len(years) == 10:  # All models go from either 1979 or 2011, so this is if it's one from 2011
#nan_years = np.zeros(32 * 12)
#tau_a_flat = np.concatenate([nan_years, tau_a_flat])
#tau_i_flat = np.concatenate([nan_years, tau_i_flat])
print(tau_a_flat.shape, tau_i_flat.shape)
print(conc_flat.shape)

# Plot timescale
total_months = np.arange(42 * 12)
f, time_ax = plt.subplots(1, 1, figsize=(18, 4), dpi=200)
time_ax: plt.Axes
conc_ax: plt.Axes = time_ax.twinx()

years_ad = total_months / 12 + 1979

# p1, = time_ax.plot(years_ad[-10 * 12:], ek_flat[-10 * 12:], label="Ekman Currents", c='C0')
p1, = time_ax.plot(years_ad, tau_flat, label=fr"Total Stress", c='C1')
p2, = time_ax.plot(years_ad, tau_a_flat,
                   label=fr"Wind Stress", c='C2')
p3, = time_ax.plot(years_ad, tau_i_flat, label=r"Ice Stress", c='C3')
# p4, = time_ax.plot(years_ad[-10 * 12:], geo_flat[-10 * 12:], label="Geostrophic Currents", c='C3')
pConc = conc_ax.fill_between(years_ad, 100 * conc_flat,
                             color='tab:cyan', alpha=0.3, label="Ice Concentration")

time_ax.set_title("Wind and Ice Stress and Ice Concentration Over Time")

time_ax.xaxis.set_major_locator(MultipleLocator(5))
time_ax.xaxis.set_minor_locator(MultipleLocator(1))
#time_ax.set_xlim(2011, 2021)

time_ax.set_ylabel("Wind / Ice Stress (Pa)")
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

plt.savefig(f"Maps_output/{hemisphere}/{model}/stress_timeseries.png")
