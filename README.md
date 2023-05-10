# Masters_project
The aim is to create plots of surface currents, Ekman pumping, and stress.

Step 1 (assuming you have all the data - see the bottom of this document if that's not the case): Create raw data arrays.

Change the parameters you want (years, months, hemisphere) in parameters.py. Every other file gets values from this script.
Run write_to_arrays.py create arrays for each variable by latitude and longitude.
The values will be monthly averages, and each array represents one year of data.
Run combine_annual_arrays.py to combine these annual arrays into arrays spanning the whole desired time period.
You should now have one array each for geo, drift, conc, and wind.


Step 2: Calculate surface currents, pumping, and stress from your data arrays.

Run write_to_Ekman_arrays.py to get a single array for your time period for each model of stress (tau), Ekman currents (ekman), and pumping (pump).
These are the models (all models besides 2 use turning angle of 45 degrees):
Model 1 - all data (2011-2020) (based on Meneghello et al. 2018)
Model 2 - all data, turning angle of 90 degrees (2011-2020) (based on Yang 2006)
Model 3 - wind only: sets ice concentration and geostrophic = 0 (1979-2020)
Model 4 - ice-free: sets ice concentration to 0 (2011-2020)
Model 5 - wind-free: sets wind to 0 (2011-2020)
Model 6 - ocean only: sets wind and ice concentration to 0 (2011-2020)
Model 7 - no geostrophic: sets geo to 0 (1979-2011)


Step 3: Plot stuff.

Make sure plotting_functions.py is imported as pf. Almost all other plotting scripts call it.

Beaufort_stress_timeseries.py plots wind, ice, and total stress and ice concentration over time.
I have it set to use model 7 for 1979-2020.
plot_avg_maps.py plots figures where the subplots are seasonal averages for each decade.
I have seasons defined as JFM, AMJ, JAS, OND.
plot_change_ek.py takes the average Ekman currents in the 2010s and subtracts it from the 1980s average, then plots that difference.
It has subplots for March and September to see ice min and max.
plot_change_pump.py and plot_change_wind.py are the same as plot_change_ek.py but for pumping and wind.
plot_maps.py plots whatever variable you want for a given month and year.
It's currently set to loop through everything, so run it overnight to wake up to a bunch of maps.
plot_maps_decades.py does the same thing as plot_avg_maps.py, but the subplots in the former are separate figures in this one.
plot_pump_timeseries.py plots a timeseries of Ekman pumping with and without geostrophy.
plot_timescale.py plots a timeseries of ice concentration, drift speed, and wind speed.
plot_wind_ice_stress.py is the same as Beaufort_stress_timeseries.py but for a whole grid, not a specific area.
plot_wind_ice_tau.py plots seasonal subplots of wind and ice stress.




	I also inherited these files and changed them minimally:
grid_set.py
Ice_Conc.py
Ice_Drift.py
data_classes.py (I technically made this by cobbling together classes that were called repeatedly,
but I didn't really edit the classes.)

I modified these files only enough to have them intake parameters from parameters.py,
so you shouldn't have to deal with any of these to run anything.

Below this is the important bits from the README of the Github folder I inherited.
https://github.com/owylie/MSci-Project-Ekman-Dynamics

All the files in this repository should be downloaded in the same folder,
and as well as this the following data sets need to be downloaded (they're too large to upload here),
and put in separate folders.


1) CPOM Geostrophic current data. http://www.cpom.ucl.ac.uk/dynamic_topography/.
   Request access, and download the file 'Full_DOT_data_Arco.nc'.
   This file should go in a folder called CPOM_geo.

ncinfo on Full DOT:

    CPOM CryoSat2 Arctic Oceanography data 2020. Written by H Heorton
    
    Geoid: GOCO03s
    
    Surface_elevation: CPOM threshold
    
    Smoothing: 100 km guassian
    
    Time_written: 2021-10-06
    
    publisher_name: UCL_CPOM
    
    publisher_type: institution
    
    publisher_email: a.muir@ucl.ac.uk
    
    publisher_url: http://www.cpom.ucl.ac.uk/dynamic_topography/
    
    Time_dimension: Days since 2000-01-01
    
    dimensions(sizes): time(120), x(334), y(334)
    
    variables(dimensions): float32 Sea_level_anomaly(time, x, y), float32 time(time), float32 DOT_smoothed(time, x, y), float32 DOT_unsmoothed(time, x, y), float32 DOT_uncertainty_estimate(time, x, y), float32 Geo_surf_current_x(time, x, y), float32 Geo_surf_current_y(time, x, y), float32 lons(x, y), float32 lats(x, y), float32 xdist(x, y), float32 ydist(x, y), float32 ang_c(x, y), float32 ang_s(x, y)


2) ERA5 wind data. This can be downloaded as an hourly or monthly data set. To work in accordance with these notebooks, it will have to be dowloaded in the following way.

   Monthly: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form.
   Select 'monthly averaged reanalysis', then '10m u-component of wind' and '10m v-component of wind', then the years 2011 to 2020, and all months, at hour 00:00 (the only option)
   For effiency, I only downloaded from 60 degrees to 90 degrees North, but a larger spatial region shouldn't affect the code.
   This will download one file with all the data, which should be named 'monthly_2011-21.nc', and put in a folder called 'ERA5'.

   Hourly: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form.
   Select the same 'reanalysis' and 10m wind components. I used 2020, all months, all hours, from 60 to 90 degrees North, but other years can be looked at by changing the code.
   This will download one file with all the data, which should be named '2020_all.nc' (or whatever year you chose), and put in a folder called 'ERA5'.

ncinfo on the file:

    <class 'netCDF4._netCDF4.Dataset'>

    root group (NETCDF3_64BIT_OFFSET data model, file format NETCDF3):

    Conventions: CF-1.6
    
    history: 2021-12-16 15:38:32 GMT 
    
    by grib_to_netcdf-2.23.0: /opt/ecmwf/mars-client/bin/grib_to_netcdf -S param -o /cache/data7/adaptor.mars.internal-1639669111.9061267-31108-1-4e62d4da-6097-40db-8e36-d4c777aeb562.nc /cache/tmp/4e62d4da-6097-40db-8e36-d4c777aeb562-adaptor.mars.internal-1639669111.0741432-31108-1-tmp.grib
    
    dimensions(sizes): longitude(1440), latitude(721), time(2)
    
    variables(dimensions): float32 longitude(longitude), float32 latitude(latitude), int32 time(time), int16 u10(time, latitude, longitude), int16 v10(time, latitude, longitude)


3) Pathfinder Ice Drift. https://nsidc.org/data/NSIDC-0116/versions/4. You will need to make an account to access this data.
   This should be downloaded in yearly .nc files (dont change the names), and put in a folder called 'Pathfinder'.

ncinfo on one of the yearly ice drift files:

    <class 'netCDF4._netCDF4.Dataset'>
    
    root group (NETCDF4 data model, file format HDF5):

    version: 4.1
    
    release_date: Apr 2021
    
    Conventions: CF-1.4
    
    dataset_doi: 10.5067/INAWUWO7QH7B
    
    dimensions(sizes): x(361), y(361), time(366)
    
    variables(dimensions): float64 x(x), float64 y(y), float64 time(time), int32 crs(), float32 u(time, y, x), float32 v(time, y, x), int16 icemotion_error_estimate(time, y, x), float32 latitude(y, x), float32 longitude(y, x)



4) NSIDC Ice Concentration. https://nsidc.org/data/NSIDC-0051/versions/1. You will have to manually download the data for each year (it takes a loooong time)
   In the 'filter by date' box, put in the year you want (will need to do at least 2011 to 2020). In 'filter spatially', range from 60 to 90 degrees North.
   Save each year of data files in a folder labelled by the year (eg. '2015'), and put all these folders in a folder called 'NSIDC_nt'.
