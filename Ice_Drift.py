import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
from sys import path
from os.path import exists
from imp import reload
from netCDF4 import Dataset
path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/')
import grid_set as gs
import parameters as par


#### setup plotting grid - Gplot
f = plt.figure()
Gplot= gs.grid_set(par.m)
ax = f.add_subplot(1,1,1,projection=par.m)
ax.set_extent(par.geo_bounds, ccrs.PlateCarree())
### make a new grid
Gplot.set_grid_mn(30,30,ax)
Gplot.get_grid_info(av_ang=False)
plt.close()

#### for ice drift
### find it here:
### https://nsidc.org/data/NSIDC-0116/versions/4
### one file per year, put them all in one folder
class Pathfinder():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,grid=False):
        #if north
        self.name = 'Pathfinder'
        self.path = ppath
        self.vyear_load = 0
        self.vels_loaded = False
        if type(grid) == bool:
            self.check_grid = False
        else:
            self.check_grid = True
            self.grid = grid
        
    def get_dates(self,time_start,time_end):
        """
        returns the all encompassing date list for use with the forcing object
        """
        dates =[]
        d0 = dt.datetime(1970,1,1)
        n_yrs = (time_end.year - time_start.year)+1
        for y in range(n_yrs):
            yu = time_start.year + y
            if par.HEMI == "north":
                f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
            if par.HEMI == "south":
                f_name = 'icemotion_daily_sh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
                ### converting days_since_1970-01-01 to datetimes
            if exists(self.path+f_name):
                f_nc = Dataset(self.path+f_name)
                [dates.append(d0 + relativedelta(days = d))
                for d in f_nc['time'][:]]
                f_nc.close()
        self.dates = dates
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')


    # next function will take a list of dates and return an appropriately orientated arrays
    # give a 
    def get_vels(self,dates_u,verbos=False):
        d0 = dt.datetime(1970,1,1)
        # does dates_u cover one year or more
        if (dates_u[-1].year -dates_u[0].year) != 0:
            raise ValueError("Can only handle one year at a time")
        # one year, one file
        yu = dates_u[0].year
        if (self.vyear_load != yu) or (not self.vels_loaded):
            print('loading new year of data: '+str(yu))
            if par.HEMI == "north":
                f_name = 'icemotion_daily_nh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
            if par.HEMI == "south":
                f_name = 'icemotion_daily_sh_25km_'+str(yu)+'0101_'+str(yu)+'1231_v4.1.nc'
            # if exists(self.path+f_name):
            f_nc = Dataset(self.path+f_name)
            self.u = f_nc['u'][:]
            self.v = f_nc['v'][:]
            self.u[self.u.mask] = np.nan
            self.v[self.v.mask] = np.nan
            f_nc.close()
            self.vyear_load = yu
            self.vels_loaded= True
        p0 = dates_u[ 0].timetuple().tm_yday -1
        p1 = dates_u[-1].timetuple().tm_yday 
        datau = self.u[p0:p1,:,:].transpose((0,2,1))/100
        datav = self.v[p0:p1,:,:].transpose((0,2,1))/100
        if self.check_grid:
            for n in range(np.shape(datau)[0]):
                datau[n][self.grid.lats>88] = np.nanmean
                datav[n][self.grid.lats>88] = np.nanmean
        return datau,datav

#Loads Pathfinder grid
GPathfinder = gs.grid_set(par.m)
GPathfinder.load_grid(par.ID_grid)
GPathfinder2Gplot = gs.Gs2Gs(GPathfinder,Gplot,vectors=True)
DRIFT = Pathfinder('/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/Pathfinder/')

class ID_av():
    """
    - For this to work, you need to have DRIFT = Pathfinder(....)
    - Also need to have defined GPathfinder2Gplot regridding
    - Month functions return the monthly averaged velocities, on GPathfinder and Gplot
    """
     
    def self_init(self):
        self.name = 'Ice_Drift_Averaging'


    def jan(year):
        
        u1,v1 = DRIFT.get_vels([dt.datetime(year,1,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,1,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,1,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,1,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,1,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,1,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,1,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,1,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,1,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,1,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,1,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,1,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,1,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,1,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,1,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,1,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,1,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,1,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,1,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,1,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,1,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,1,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,1,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,1,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,1,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,1,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,1,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,1,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,1,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,1,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,1,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        jan_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,
                                   u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        jan_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20,
                                   v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        jan_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        jan_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )
        
        return jan_u, jan_v, jan_ur, jan_vr
    
    
    
    def feb(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,2,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,2,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,2,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,2,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,2,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,2,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,2,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,2,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,2,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,2,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,2,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,2,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,2,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,2,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,2,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,2,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,2,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,2,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,2,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,2,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,2,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,2,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,2,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,2,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,2,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,2,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,2,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,2,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        feb_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                           u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28]), axis=0 )
        feb_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                           v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28]), axis=0 )
        feb_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,
                            ur16,ur17,ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28]), axis=0 )
        feb_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,
                            vr16,vr17,vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28]), axis=0 )
        
        return feb_u, feb_v, feb_ur, feb_vr
    
    
    def mar(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,3,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,3,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,3,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,3,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,3,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,3,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,3,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,3,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,3,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,3,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,3,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,3,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,3,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,3,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,3,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,3,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,3,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,3,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,3,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,3,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,3,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,3,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,3,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,3,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,3,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,3,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,3,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,3,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,3,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,3,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,3,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        mar_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        mar_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        mar_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        mar_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )
        
        return mar_u, mar_v, mar_ur, mar_vr
    
    
    def apr(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,4,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,4,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,4,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,4,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,4,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,4,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,4,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,4,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,4,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,4,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,4,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,4,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,4,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,4,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,4,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,4,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,4,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,4,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,4,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,4,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,4,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,4,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,4,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,4,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,4,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,4,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,4,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,4,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,4,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,4,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        apr_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30]), axis=0 )
        apr_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30]), axis=0 )
        apr_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30]), axis=0 )
        apr_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30]), axis=0 )
        
        return apr_u, apr_v, apr_ur, apr_vr
    
    
    def may(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,5,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,5,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,5,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,5,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,5,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,5,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,5,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,5,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,5,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,5,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,5,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,5,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,5,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,5,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,5,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,5,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,5,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,5,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,5,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,5,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,5,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,5,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,5,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,5,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,5,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,5,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,5,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,5,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,5,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,5,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,5,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        may_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        may_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        may_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        may_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )
        
        return may_u, may_v, may_ur, may_vr
    
    
    def jun(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,6,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,6,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,6,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,6,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,6,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,6,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,6,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,6,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,6,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,6,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,6,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,6,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,6,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,6,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,6,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,6,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,6,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,6,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,6,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,6,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,6,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,6,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,6,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,6,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,6,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,6,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,6,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,6,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,6,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,6,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        jun_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30]), axis=0 )
        jun_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30]), axis=0 )
        jun_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30]), axis=0 )
        jun_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30]), axis=0 )
        
        return jun_u, jun_v, jun_ur, jun_vr
    
    
    def jul(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,7,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,7,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,7,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,7,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,7,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,7,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,7,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,7,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,7,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,7,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,7,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,7,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,7,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,7,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,7,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,7,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,7,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,7,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,7,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,7,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,7,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,7,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,7,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,7,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,7,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,7,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,7,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,7,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,7,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,7,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,7,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        jul_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        jul_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        jul_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        jul_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )
        
        return jul_u, jul_v, jul_ur, jul_vr
    
    
    def aug(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,8,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,8,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,8,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,8,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,8,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,8,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,8,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,8,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,8,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,8,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,8,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,8,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,8,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,8,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,8,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,8,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,8,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,8,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,8,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,8,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,8,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,8,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,8,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,8,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,8,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,8,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,8,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,8,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,8,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,8,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,8,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        aug_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        aug_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        aug_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        aug_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )
        
        return aug_u, aug_v, aug_ur, aug_vr
    
    
    def sep(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,9,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,9,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,9,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,9,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,9,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,9,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,9,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,9,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,9,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,9,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,9,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,9,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,9,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,9,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,9,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,9,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,9,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,9,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,9,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,9,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,9,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,9,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,9,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,9,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,9,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,9,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,9,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,9,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,9,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,9,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        sep_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30]), axis=0 )
        sep_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30]), axis=0 )
        sep_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30]), axis=0 )
        sep_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30]), axis=0 )
        
        return sep_u, sep_v, sep_ur, sep_vr
    
    
    def oct1(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,10,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,10,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,10,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,10,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,10,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,10,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,10,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,10,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,10,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,10,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,10,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,10,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,10,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,10,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,10,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,10,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,10,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,10,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,10,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,10,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,10,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,10,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,10,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,10,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,10,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,10,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,10,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,10,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,10,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,10,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,10,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        oct1_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                    u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        oct1_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                    v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        oct1_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                     ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        oct1_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                     vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )
        
        return oct1_u, oct1_v, oct1_ur, oct1_vr
    
    
    def nov(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,11,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,11,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,11,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,11,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,11,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,11,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,11,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,11,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,11,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,11,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,11,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)
    
        u12,v12 = DRIFT.get_vels([dt.datetime(year,11,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,11,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,11,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,11,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,11,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,11,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,11,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,11,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,11,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,11,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,11,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,11,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,11,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,11,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,11,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,11,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)
    
        u28,v28 = DRIFT.get_vels([dt.datetime(year,11,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,11,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,11,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        nov_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30]), axis=0 )
        nov_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30]), axis=0 )
        nov_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30]), axis=0 )
        nov_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30]), axis=0 )
        
        return nov_u, nov_v, nov_ur, nov_vr
    
    
    def dec(year):
        u1,v1 = DRIFT.get_vels([dt.datetime(year,12,1)],verbos=True)
        ur1,vr1 = GPathfinder2Gplot.rg_vecs(u1,v1)

        u2,v2 = DRIFT.get_vels([dt.datetime(year,12,2)],verbos=True)
        ur2,vr2 = GPathfinder2Gplot.rg_vecs(u2,v2)

        u3,v3 = DRIFT.get_vels([dt.datetime(year,12,3)],verbos=True)
        ur3,vr3 = GPathfinder2Gplot.rg_vecs(u3,v3)

        u4,v4 = DRIFT.get_vels([dt.datetime(year,12,4)],verbos=True)
        ur4,vr4 = GPathfinder2Gplot.rg_vecs(u4,v4)

        u5,v5 = DRIFT.get_vels([dt.datetime(year,12,5)],verbos=True)
        ur5,vr5 = GPathfinder2Gplot.rg_vecs(u5,v5)

        u6,v6 = DRIFT.get_vels([dt.datetime(year,12,6)],verbos=True)
        ur6,vr6 = GPathfinder2Gplot.rg_vecs(u6,v6)

        u7,v7 = DRIFT.get_vels([dt.datetime(year,12,7)],verbos=True)
        ur7,vr7 = GPathfinder2Gplot.rg_vecs(u7,v7)

        u8,v8 = DRIFT.get_vels([dt.datetime(year,12,8)],verbos=True)
        ur8,vr8 = GPathfinder2Gplot.rg_vecs(u8,v8)

        u9,v9 = DRIFT.get_vels([dt.datetime(year,12,9)],verbos=True)
        ur9,vr9 = GPathfinder2Gplot.rg_vecs(u9,v9)

        u10,v10 = DRIFT.get_vels([dt.datetime(year,12,10)],verbos=True)
        ur10,vr10 = GPathfinder2Gplot.rg_vecs(u10,v10)

        u11,v11 = DRIFT.get_vels([dt.datetime(year,12,11)],verbos=True)
        ur11,vr11 = GPathfinder2Gplot.rg_vecs(u11,v11)

        u12,v12 = DRIFT.get_vels([dt.datetime(year,12,12)],verbos=True)
        ur12,vr12 = GPathfinder2Gplot.rg_vecs(u12,v12)

        u13,v13 = DRIFT.get_vels([dt.datetime(year,12,13)],verbos=True)
        ur13,vr13 = GPathfinder2Gplot.rg_vecs(u13,v13)

        u14,v14 = DRIFT.get_vels([dt.datetime(year,12,14)],verbos=True)
        ur14,vr14 = GPathfinder2Gplot.rg_vecs(u14,v14)

        u15,v15 = DRIFT.get_vels([dt.datetime(year,12,15)],verbos=True)
        ur15,vr15 = GPathfinder2Gplot.rg_vecs(u15,v15)

        u16,v16 = DRIFT.get_vels([dt.datetime(year,12,16)],verbos=True)
        ur16,vr16 = GPathfinder2Gplot.rg_vecs(u16,v16)

        u17,v17 = DRIFT.get_vels([dt.datetime(year,12,17)],verbos=True)
        ur17,vr17 = GPathfinder2Gplot.rg_vecs(u17,v17)

        u18,v18 = DRIFT.get_vels([dt.datetime(year,12,18)],verbos=True)
        ur18,vr18 = GPathfinder2Gplot.rg_vecs(u18,v18)

        u19,v19 = DRIFT.get_vels([dt.datetime(year,12,19)],verbos=True)
        ur19,vr19 = GPathfinder2Gplot.rg_vecs(u19,v19)

        u20,v20 = DRIFT.get_vels([dt.datetime(year,12,20)],verbos=True)
        ur20,vr20 = GPathfinder2Gplot.rg_vecs(u20,v20)

        u21,v21 = DRIFT.get_vels([dt.datetime(year,12,21)],verbos=True)
        ur21,vr21 = GPathfinder2Gplot.rg_vecs(u21,v21)

        u22,v22 = DRIFT.get_vels([dt.datetime(year,12,22)],verbos=True)
        ur22,vr22 = GPathfinder2Gplot.rg_vecs(u22,v22)

        u23,v23 = DRIFT.get_vels([dt.datetime(year,12,23)],verbos=True)
        ur23,vr23 = GPathfinder2Gplot.rg_vecs(u23,v23)

        u24,v24 = DRIFT.get_vels([dt.datetime(year,12,24)],verbos=True)
        ur24,vr24 = GPathfinder2Gplot.rg_vecs(u24,v24)

        u25,v25 = DRIFT.get_vels([dt.datetime(year,12,25)],verbos=True)
        ur25,vr25 = GPathfinder2Gplot.rg_vecs(u25,v25)

        u26,v26 = DRIFT.get_vels([dt.datetime(year,12,26)],verbos=True)
        ur26,vr26 = GPathfinder2Gplot.rg_vecs(u26,v26)

        u27,v27 = DRIFT.get_vels([dt.datetime(year,12,27)],verbos=True)
        ur27,vr27 = GPathfinder2Gplot.rg_vecs(u27,v27)

        u28,v28 = DRIFT.get_vels([dt.datetime(year,12,28)],verbos=True)
        ur28,vr28 = GPathfinder2Gplot.rg_vecs(u28,v28)

        u29,v29 = DRIFT.get_vels([dt.datetime(year,12,29)],verbos=True)
        ur29,vr29 = GPathfinder2Gplot.rg_vecs(u29,v29)

        u30,v30 = DRIFT.get_vels([dt.datetime(year,12,30)],verbos=True)
        ur30,vr30 = GPathfinder2Gplot.rg_vecs(u30,v30)

        u31,v31 = DRIFT.get_vels([dt.datetime(year,12,31)],verbos=True)
        ur31,vr31 = GPathfinder2Gplot.rg_vecs(u31,v31)

        dec_u = np.nanmean( np.array([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15,
                                   u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31]), axis=0 )
        dec_v = np.nanmean( np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,
                                   v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31]), axis=0 )
        dec_ur = np.nanmean( np.array([ur1,ur2,ur3,ur4,ur5,ur6,ur7,ur8,ur9,ur10,ur11,ur12,ur13,ur14,ur15,ur16,ur17,
                                    ur18,ur19,ur20,ur21,ur22,ur23,ur24,ur25,ur26,ur27,ur28,ur29,ur30,ur31]), axis=0 )
        dec_vr = np.nanmean( np.array([vr1,vr2,vr3,vr4,vr5,vr6,vr7,vr8,vr9,vr10,vr11,vr12,vr13,vr14,vr15,vr16,vr17,
                                    vr18,vr19,vr20,vr21,vr22,vr23,vr24,vr25,vr26,vr27,vr28,vr29,vr30,vr31]), axis=0 )

        return dec_u, dec_v, dec_ur, dec_vr