import numpy as np
import datetime as dt
import struct
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from os.path import exists
import glob
from scipy.interpolate import griddata
from scipy.interpolate.interpnd import _ndim_coords_from_arrays
from scipy.spatial import cKDTree
import parameters as par


class Pathfinder():
    """
    forcing class for the budget
    lets the forcing load efficiently

    """

    def __init__(self, ppath, grid=False):
        self.name = 'Pathfinder'
        self.path = ppath
        # self.hemi = hemi
        self.vyear_load = 0
        self.erryear_load = 0
        self.vels_loaded = False
        self.err_loaded = False
        if type(grid) == bool:
            self.check_grid = False
        else:
            self.check_grid = True
            self.grid = grid

    def get_dates(self, time_start, time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        dates = []
        d0 = dt.datetime(1970, 1, 1)
        n_yrs = (time_end.year - time_start.year) + 1
        for y in range(n_yrs):
            yu = time_start.year + y
            if par.HEMI == 'north':
                f_name = 'icemotion_daily_nh_25km_' + str(yu) + '0101_' + str(yu) + '1231_v4.1.nc'
            elif par.HEMI == 'south':
                ### icemotion_daily_sh_25km_20100101_20101231_v4.1.nc
                f_name = 'icemotion_daily_sh_25km_' + str(yu) + '0101_' + str(yu) + '1231_v4.1.nc'
            if exists(self.path + f_name):
                f_nc = Dataset(self.path + f_name)
                [dates.append(d0 + relativedelta(days=d))
                 for d in f_nc['time'][:]]
                f_nc.close()
        self.dates = dates
        print(self.name + ' Found ' + str(np.shape(dates)[0]) + ' dates')

    # daily points in yearly files

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a
    def get_vels(self, dates_u, verbos=False):
        d0 = dt.datetime(1970, 1, 1)
        # does dates_u cover one year or more
        if (dates_u[-1].year - dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.vyear_load != yu) or (not self.vels_loaded)):
                print('loading new year of data: ' + str(yu))
                if par.HEMI == 'north':
                    f_name = 'icemotion_daily_nh_25km_' + str(yu) + '0101_' + str(yu) + '1231_v4.1.nc'
                elif par.HEMI == 'south':
                    ### icemotion_daily_sh_25km_20100101_20101231_v4.1.nc
                    f_name = 'icemotion_daily_sh_25km_' + str(yu) + '0101_' + str(yu) + '1231_v4.1.nc'
                else:
                    print("Define a hemisphere in parameters.py. If you've already done that, check your spelling.")
                #                 f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                f_nc = Dataset(self.path + f_name)
                #         print(p0,p1)
                self.u = f_nc['u'][:]
                self.v = f_nc['v'][:]
                self.u[self.u.mask] = np.nan
                self.v[self.v.mask] = np.nan
                f_nc.close()
                self.vyear_load = yu
                self.vels_loaded = True
            p0 = dates_u[0].timetuple().tm_yday - 1
            p1 = dates_u[-1].timetuple().tm_yday
            datau = self.u[p0:p1, :, :].transpose((0, 2, 1)) / 100
            datav = self.v[p0:p1, :, :].transpose((0, 2, 1)) / 100
            if self.check_grid:
                for n in range(np.shape(datau)[0]):
                    datau[n][self.grid.lats > 88] = np.nan
                    datav[n][self.grid.lats > 88] = np.nan
            return datau, datav

    def get_err(self, dates_u, verbos=False):
        ## errs need to be dimensional
        ### always make sure you call get_vels frist
        d0 = dt.datetime(1970, 1, 1)
        # does dates_u cover one year or more
        if (dates_u[-1].year - dates_u[0].year) == 0:
            # one year, one file
            yu = dates_u[0].year
            if ((self.erryear_load != yu) or (not self.err_loaded)):
                print('loading new year of data: ' + str(yu))
                #                 f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                if par.HEMI == 'north':
                    f_name = 'icemotion_daily_nh_25km_' + str(yu) + '0101_' + str(yu) + '1231_v4.1.nc'
                elif par.HEMI == 'south':
                    ### icemotion_daily_sh_25km_20100101_20101231_v4.1.nc
                    f_name = 'icemotion_daily_sh_25km_' + str(yu) + '0101_' + str(yu) + '1231_v4.1.nc'
                f_nc = Dataset(self.path + f_name)
                #         print(p0,p1)
                self.err = f_nc['icemotion_error_estimate'][:]
                self.err = self.err / 100.0  # percent to fraction
                self.err[self.err.mask] = np.nan
                self.err[self.err < 0] = np.nan
                #                 self.uerr = self.err*self.u
                #                 self.verr = self.err*self.v
                f_nc.close()
                self.erryear_load = yu
                self.err_loaded = True
            p0 = dates_u[0].timetuple().tm_yday - 1
            p1 = dates_u[-1].timetuple().tm_yday
            data_err = self.err[p0:p1, :, :].transpose((0, 2, 1))
            #             data_uerr = self.uerr[p0:p1,:,:].transpose((0,2,1))
            #             data_verr = self.verr[p0:p1,:,:].transpose((0,2,1))
            if self.check_grid:
                for n in range(np.shape(datau)[0]):
                    data_err[n][self.grid.lats > 88] = np.nanmean
            #                     data_uerr[n][self.grid.lats>88] = np.nanmean
            #                     data_verr[n][self.grid.lats>88] = np.nanmean
            return data_err


#             return data_uerr,data_verr

class Pathfinder_weekly():
    """
    forcing class for the budget
    lets the forcing load efficiently
    assumes a single NRT nc file for all the Pathfinder
    """

    def __init__(self, ppath, grid=False):
        self.name = 'Pathfinder_weekly'
        self.path = ppath
        self.files = glob.glob(ppath + '*.nc')
        self.files.sort()
        self.vels_loaded = False
        self.vyear_load = 0

        if type(grid) == bool:
            self.check_grid = False
        else:
            self.check_grid = True
            self.grid = grid

    def get_dates(self, time_start, time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        self.ds = []
        self.de = []
        self.all_dates = []
        self.dates = []
        d0 = dt.datetime(1970, 1, 1)
        for f in self.files:
            Pnc = Dataset(f)
            ds = d0 + relativedelta(days=int(Pnc.variables['time'][0]))
            de = d0 + relativedelta(days=int(Pnc.variables['time'][-1]))
            #             print(ds.year,time_start.year,de.year,time_end.year)
            if ds.year >= time_start.year and de.year <= time_end.year:
                self.ds.append(ds)
                self.de.append(de)
                self.all_dates.append([d0 + relativedelta(days=int(d)) for d in Pnc.variables['time']])
                self.dates.extend([t for t in self.all_dates[-1] if t >= time_start and t <= time_end])
            Pnc.close()
        print(self.name + ' Found ' + str(np.shape(self.dates)[0]) + ' dates')

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a
    def get_vels(self, dates_u, verbos=False):
        set_var_now = False
        join_data = False

        u_arrays = []
        v_arrays = []
        ### loop the dates
        for d in dates_u:
            iload = np.searchsorted(np.array(self.ds), d, side='right') - 1
            if (d.year != self.vyear_load) or (not self.vels_loaded):
                self.load_vels([d], verbos=verbos)
            ### fill the out arrays
            load_point = self.all_dates[iload].index(d)
            u_arrays.append(self.u[load_point])
            v_arrays.append(self.v[load_point])

        datau = np.ma.array(u_arrays)
        datav = np.ma.array(v_arrays)
        datau[datau.mask] = np.nan
        datav[datav.mask] = np.nan
        if self.check_grid:
            for n in range(np.shape(datau)[0]):
                datau[n][self.grid.lats > 88] = np.nan
                datav[n][self.grid.lats > 88] = np.nan
        return datau, datav

    def load_vels(self, dates_u, verbos=False):
        iload = np.searchsorted(np.array(self.ds), dates_u[0], side='right') - 1
        Pnc = Dataset(self.files[iload])
        if verbos: print(self.files[iload])
        self.u = Pnc.variables['u'][:].transpose((0, 2, 1)) / 100
        self.v = Pnc.variables['v'][:].transpose((0, 2, 1)) / 100
        self.vels_loaded = True
        self.vyear_load = self.ds[iload].year
        Pnc.close()
        print('Filling buffer from ' + self.ds[iload].strftime('%Y-%m-%d'))


# ERA5 10m MONTHLY winds
class ERA5_months():
    """
    forcing class for the budget
    lets the forcing load efficiently
   
    """

    def __init__(self, ppath):
        self.name = 'ERA5_Monthly_Winds'
        self.path = ppath
        if par.HEMI == "north":
            self.file = self.path + 'monthly_north.nc'
        if par.HEMI == "south":
            self.file = self.path + 'monthly_south.nc'
        self.f_nc = Dataset(self.file)

    def get_dates(self, month_wanted):
        """
        returns the all encompassing date list for use with the forcing object
        """
        d0 = dt.datetime(2011, 1, 1)  # first month of dataset
        if par.HEMI == "north":
            self.file = self.path + 'monthly_north.nc'
        if par.HEMI == "south":
            self.file = self.path + 'monthly_south.nc'
        self.f_nc = Dataset(self.file)
        self.months = self.f_nc.variables['time'][:]
        # this makes sure we have the month that is asked for:
        self.months_all = [d0 + relativedelta(months=m) for m in range(len(self.months))]
        self.times = [d for d in self.months_all if d >= month_wanted and d < month_wanted + relativedelta(months=1)]
        # print(self.name+' Found '+str(np.shape(self.times)[0])+' data points')

    def get_vels(self):
        """
        NOTE: gives the velocities of the month specified in get_dates above, NOT in get_vels as usual
        """
        # finding the index:
        n = 0
        for n in range(0, len(self.months)):
            if self.months_all[n] in self.times:
                idx = n
            else:
                # raise RuntimeError("oh fuck")
                continue

        u10 = self.f_nc.variables['u10'][idx]
        v10 = self.f_nc.variables['v10'][idx]

        return u10, v10


class NSIDC_nt():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """

    def __init__(self, ppath):
        self.name = 'NSIDC_n'
        self.path = ppath

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a
    def get_aice(self, dates_u, verbos=False):
        # does dates_u cover one year or more
        # daily files
        changed_names = ["f17", "f08", "f11", "f13", "n07"]  # options for filename
        if par.HEMI == "north":
            dimY = 304
            dimX = 448
        if par.HEMI == "south":
            dimY = 316
            dimX = 332

        d_no = np.shape(dates_u)[0]
        # noinspection PyUnboundLocalVariable
        data = np.full([d_no, dimX, dimY], np.nan)
        for n, d in enumerate(dates_u):
            for name in changed_names:
                if par.HEMI == "north":
                    infile = self.path + d.strftime('/%Y/') + d.strftime('/%Y.%m.%d/') + "nt_" + d.strftime(
                        '%Y%m%d_') + name + "_v1.1_n.bin"
                    # Example: Masters_project/NSIDC_nt/2011/2011.01.01/nt_20110101_f17_v1.1_n.bin
                if par.HEMI == "south":
                    infile = self.path + d.strftime('/%Y/') + d.strftime('/%Y.%m.%d/') + "nt_" + d.strftime(
                        '%Y%m%d_') + name + "_v1.1_s.bin"
                # noinspection PyUnboundLocalVariable
                if exists(infile):
                    with open(infile, 'rb') as fr:
                        hdr = fr.read(300)
                        ice = np.fromfile(fr, dtype=np.uint8)

                    ice = ice.reshape(dimX, dimY)
                    ice = np.flipud(ice)
                    data[n] = ice / 250.
                    break
            else:
                # TODO: HANDLE THIS
                print(f"No file found on date {d}")
        data[data > 1.0] = np.nan
        return data

    def get_dates(self, time_start, time_end):
        # does dates_u cover one year or more
        # daily files
        dates_u = []
        d_no = (time_end - time_start).days + 3
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start + relativedelta(days=dn - 1)
            # if d>=dt.datetime(2020,11,1):
            #    infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            # else:
            #             if d.year>2019:
            if par.HEMI == "north":
                infile = self.path + d.strftime('/%Y/') + "nt_" + d.strftime('%Y%m%d') + "_f17_v1.1_n.bin"
            if par.HEMI == "south":
                infile = self.path + d.strftime('/%Y/') + "nt_" + d.strftime('%Y%m%d') + "_f17_v1.1_s.bin"
            # check infile exists 
            if exists(infile):
                dates_u.append(d)
            # if it does append dates_u
        self.dates = dates_u
        print(self.name + ' Found ' + str(np.shape(dates_u)[0]) + ' dates')


#### for geostropic vels
### find it here:
### http://www.cpom.ucl.ac.uk/dynamic_topography/
### one file for all of it
class CPOM_geo():
    """
    forcing class for the budget
    lets the forcing load efficiently
    """

    def __init__(self, ppath, grid=False):
        self.name = 'CPOM_Geostrophic_Currents'
        d0 = dt.datetime(2000, 1, 1)
        self.path = ppath
        if par.HEMI == "north":
            self.file = self.path + 'Full_DOT_data_Arco.nc'
        if par.HEMI == "south":
            self.file = self.path + 'Full_DOT_data_Anto.nc'
        self.f_nc = Dataset(self.file)
        self.time_vec = self.f_nc.variables['time'][:]
        self.dates_all = [d0 + relativedelta(days=t) for t in self.time_vec]

    def get_dates(self, time_start, time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        d0 = dt.datetime(2000, 1, 1)
        if par.HEMI == "north":
            self.file = self.path + 'Full_DOT_data_Arco.nc'
        if par.HEMI == "south":
            self.file = self.path + 'Full_DOT_data_Anto.nc'
        self.f_nc = Dataset(self.file)
        self.time_vec = self.f_nc.variables['time'][:]
        self.dates_all = [d0 + relativedelta(days=t) for t in self.time_vec]
        self.dates = [d for d in self.dates_all if time_start <= d <= time_end]
        print(self.name + ' Found ' + str(np.shape(self.dates)[0]) + ' dates')

    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self, dates_u, verbos=False):
        d0 = dt.datetime(2000, 1, 1)
        ### find the indices
        idx = [np.argwhere(np.array([d == du for d in self.dates_all]))[0, 0]
               for du in dates_u]
        ### little bit of checking
        if verbos:
            for i, du in zip(idx, dates_u):
                dcheck = dt.datetime(1, 1, 1)
                dcheck = dcheck + relativedelta(days=self.time_vec[i])
                dcheck = dcheck + relativedelta(years=-1)
                print(du.strftime('%Y%m%d-') + dcheck.strftime('%Y%m%d'))
        datau = self.f_nc.variables['Geo_surf_current_x'][idx]
        datav = -self.f_nc.variables['Geo_surf_current_y'][idx]
        return datau, datav
