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
# path.insert(0, '/Users/H/WAVES/geo_data_group/')
path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/')
import grid_set as gs
import parameters as par

class NSIDC_nt():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    
    def __init__(self,ppath):
        self.name = 'NSIDC_n'
        self.path = ppath

        
# next function will take a list of dates and return an appropriately orientated arrays
# This version will change filenames as NASA changed them over time
# Testing with constant name "f17"
    def get_aice(self,dates_u,verbos=False):
        # does dates_u cover one year or more
        #daily files
        changed_names = ["f17","f08","f11","f13","n07"] #options for filename
        if par.HEMI == "north":
            dimY = 304
            dimX = 448
        if par.HEMI == "south":
            dimY = 316
            dimX = 332

        d_no = np.shape(dates_u)[0]
        data =  np.empty([d_no, dimX, dimY])
        for n,d in enumerate(dates_u):
            for name in changed_names:
                if par.HEMI == "north":
                    infile = self.path+d.strftime('/%Y/')+d.strftime('/%Y.%m.%d/')+"nt_"+d.strftime('%Y%m%d_')+name+"_v1.1_n.bin"
                    # Example: Masters_project/NSIDC_nt/2011/2011.01.01/nt_20110101_f17_v1.1_n.bin
                if par.HEMI == "south":
                    infile = self.path+d.strftime('/%Y/')+d.strftime('/%Y.%m.%d/')+"nt_"+d.strftime('%Y%m%d_')+name+"_v1.1_s.bin"
                if exists(infile):
                    with open(infile, 'rb') as fr:
                        hdr = fr.read(300)
                        ice = np.fromfile(fr, dtype=np.uint8)

                    ice = ice.reshape(dimX,dimY)
                    ice = np.flipud(ice)
                    data[n] = ice / 250.
                data[data>1.0] = np.nan
        return data
    
#     def get_aice(self,dates_u,verbos=False):
#         # does dates_u cover one year or more
#         #daily files
#         dimY = 304
#         dimX = 448
#         d_no = np.shape(dates_u)[0]
#         data =  np.empty([d_no, dimX, dimY])
#         for n,d in enumerate(dates_u):
#             infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
#             infile = f"/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/NSIDC_nt/{year}/nt_{year}{month_index+1}02_7_v1.1_n.bin"
#             print(infile)
#             with open(infile, 'rb') as fr:
#                 hdr = fr.read(300)
#                 ice = np.fromfile(fr, dtype=np.uint8)

#             ice = ice.reshape(dimX,dimY)
#             ice = np.flipud(ice)
#             data[n] = ice / 250.
#         data[data>1.0] = np.nan
#         return data

    def get_dates(self,time_start,time_end):
        # does dates_u cover one year or more
        #daily files
        dates_u = []
        d_no = (time_end-time_start).days +3 
        # make sure we get the bracket points
        for dn in range(d_no):
            d = time_start+ relativedelta(days = dn - 1)
            #if d>=dt.datetime(2020,11,1):
            #    infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f18_nrt_n.bin"
            #else:
#             if d.year>2019:
            infile = self.path+d.strftime('/%Y/')+"nt_"+d.strftime('%Y%m%d')+"_f17_v1.1_n.bin"
            # check infile exists 
            if exists(infile):
                dates_u.append(d)
            #if it does append dates_u
        self.dates= dates_u
        print(self.name+' Found '+str(np.shape(dates_u)[0])+' dates')



        
        

#grid
m = par.m
GIC = gs.grid_set(m)
GIC.load_grid(par.IC_grid)

IC = NSIDC_nt('/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/NSIDC_nt')

GPathfinder = gs.grid_set(par.m)
GPathfinder.load_grid(par.ID_grid)
GIC2GPathfinder = gs.Gs2Gs(GIC,GPathfinder,vectors=False)


class IC_av():
    """
    - For this to work, you need to have IC = NSIDC_nt(....)
    - Also need to have defined GIC2GPathfinder regridding
    - Month functions return the monthly concentrations, on GPathfinder and GIC
    """


    # @staticmethod
    def jan(year):
        a1 = IC.get_aice([dt.datetime(year,1,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,1,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,1,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,1,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,1,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,1,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,1,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,1,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,1,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,1,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,1,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,1,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,1,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,1,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,1,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,1,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,1,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,1,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,1,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,1,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,1,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,1,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,1,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,1,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,1,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,1,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,1,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,1,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])

        a29 = IC.get_aice([dt.datetime(year,1,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,1,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,1,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        jan = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        jan_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return jan, jan_r
    
    
    def feb(year):
        a1 = IC.get_aice([dt.datetime(year,12,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,2,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,2,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,2,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,2,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,2,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,2,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,2,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,2,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,2,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,2,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,2,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,2,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,2,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,2,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,2,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,2,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,2,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,2,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,2,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,2,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,2,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,2,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,2,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,2,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,2,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,2,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,2,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])

        feb = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28]), axis=0 )
        feb_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28]), axis=0 )
        
        return feb, feb_r
    
    
    def mar(year):
        a1 = IC.get_aice([dt.datetime(year,3,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,3,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,3,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,3,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,3,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,3,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,3,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,3,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,3,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,3,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,3,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,3,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,3,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,3,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,3,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,3,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,3,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,3,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,3,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,3,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,3,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,3,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,3,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,3,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,3,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,3,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,3,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,3,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])
        
        a29 = IC.get_aice([dt.datetime(year,3,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,3,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,3,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        mar = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        mar_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return mar, mar_r
    
    
    def apr(year):
        a1 = IC.get_aice([dt.datetime(year,4,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,4,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,4,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,4,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,4,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,4,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,4,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,4,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,4,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,4,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,4,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,4,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,4,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,4,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,4,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,4,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,4,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,4,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,4,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,4,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,4,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,4,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,4,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,4,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,4,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,4,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,4,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,4,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])

        a29 = IC.get_aice([dt.datetime(year,4,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,4,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        apr = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30]), axis=0 )
        apr_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30]), axis=0 )
        
        return apr, apr_r
    
    
    def may(year):
        a1 = IC.get_aice([dt.datetime(year,5,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,5,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,5,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,5,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,5,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,5,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,5,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,5,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,5,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,5,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,5,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,5,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,5,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,5,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,5,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,5,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,5,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,5,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,5,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,5,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,5,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,5,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,5,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,5,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,5,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,5,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,5,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,5,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])
        
        a29 = IC.get_aice([dt.datetime(year,5,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,5,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,5,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        may = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        may_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return may, may_r
    
    
    def jun(year):
        a1 = IC.get_aice([dt.datetime(year,6,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,6,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,6,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,6,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,6,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,6,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,6,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,6,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,6,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,6,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,6,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,6,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,6,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,6,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,6,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,6,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,6,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,6,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,6,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,6,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,6,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,6,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,6,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,6,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,6,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,6,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,6,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,6,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])
        
        a29 = IC.get_aice([dt.datetime(year,6,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,6,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        jun = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30]), axis=0 )
        jun_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30]), axis=0 )
        
        return jun, jun_r
    
    
    def jul(year):
        a1 = IC.get_aice([dt.datetime(year,7,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,7,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,7,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,7,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,7,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,7,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,7,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,7,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,7,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,7,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,7,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,7,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,7,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,7,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,7,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,7,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,7,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,7,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])
        
        a19 = IC.get_aice([dt.datetime(year,7,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,7,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,7,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,7,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,7,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,7,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,7,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,7,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,7,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,7,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])
        
        a29 = IC.get_aice([dt.datetime(year,7,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,7,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,7,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        jul = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        jul_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return jul, jul_r
    
    
    def aug(year):
        a1 = IC.get_aice([dt.datetime(year,8,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,8,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,8,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,8,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,8,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,8,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,8,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,8,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,8,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,8,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,8,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,8,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,8,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,8,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,8,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,8,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,8,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,8,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,8,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,8,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,8,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,8,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,8,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,8,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,8,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,8,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,8,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])
        
        a28 = IC.get_aice([dt.datetime(year,8,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])
        
        a29 = IC.get_aice([dt.datetime(year,8,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,8,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,8,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        aug = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        aug_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return aug, aug_r
    
    
    def sep(year):
        a1 = IC.get_aice([dt.datetime(year,9,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,9,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,9,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,9,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,9,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,9,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,9,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,9,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,9,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,9,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,9,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,9,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,9,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,9,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,9,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,9,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,9,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,9,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,9,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,9,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,9,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,9,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,9,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,9,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,9,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,9,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,9,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,9,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])
        
        a29 = IC.get_aice([dt.datetime(year,9,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,9,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        sep = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30]), axis=0 )
        sep_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30]), axis=0 )
        
        return sep, sep_r
    
    
    def oct1(year):
        a1 = IC.get_aice([dt.datetime(year,10,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,10,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,10,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,10,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,10,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,10,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,10,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,10,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,10,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,10,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,10,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,10,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,10,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,10,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,10,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,10,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,10,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,10,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,10,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,10,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,10,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,10,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,10,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,10,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,10,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,10,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,10,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,10,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])

        a29 = IC.get_aice([dt.datetime(year,10,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,10,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,10,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        oct1 = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        oct1_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return oct1, oct1_r
    
    
    def nov(year):
        a1 = IC.get_aice([dt.datetime(year,11,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,11,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,11,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,11,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,11,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,11,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,11,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,11,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,11,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,11,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,11,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,11,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,11,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,11,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,11,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,11,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,11,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,11,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,11,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,11,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,11,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,11,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,11,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,11,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,11,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,11,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,11,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,11,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])

        a29 = IC.get_aice([dt.datetime(year,11,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,11,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        nov = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30]), axis=0 )
        nov_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30]), axis=0 )
        
        return nov, nov_r
    
    
    def dec(year):
        a1 = IC.get_aice([dt.datetime(year,12,1)],verbos=True)
        ar1 = GIC2GPathfinder.rg_array(a1[0])

        a2 = IC.get_aice([dt.datetime(year,12,2)],verbos=True)
        ar2 = GIC2GPathfinder.rg_array(a2[0])

        a3 = IC.get_aice([dt.datetime(year,12,3)],verbos=True)
        ar3 = GIC2GPathfinder.rg_array(a3[0])

        a4 = IC.get_aice([dt.datetime(year,12,4)],verbos=True)
        ar4 = GIC2GPathfinder.rg_array(a4[0])

        a5 = IC.get_aice([dt.datetime(year,12,5)],verbos=True)
        ar5 = GIC2GPathfinder.rg_array(a5[0])

        a6 = IC.get_aice([dt.datetime(year,12,6)],verbos=True)
        ar6 = GIC2GPathfinder.rg_array(a6[0])

        a7 = IC.get_aice([dt.datetime(year,12,7)],verbos=True)
        ar7 = GIC2GPathfinder.rg_array(a7[0])

        a8 = IC.get_aice([dt.datetime(year,12,8)],verbos=True)
        ar8 = GIC2GPathfinder.rg_array(a8[0])

        a9 = IC.get_aice([dt.datetime(year,12,9)],verbos=True)
        ar9 = GIC2GPathfinder.rg_array(a9[0])

        a10 = IC.get_aice([dt.datetime(year,12,10)],verbos=True)
        ar10 = GIC2GPathfinder.rg_array(a10[0])

        a11 = IC.get_aice([dt.datetime(year,12,11)],verbos=True)
        ar11 = GIC2GPathfinder.rg_array(a11[0])

        a12 = IC.get_aice([dt.datetime(year,12,12)],verbos=True)
        ar12 = GIC2GPathfinder.rg_array(a12[0])

        a13 = IC.get_aice([dt.datetime(year,12,13)],verbos=True)
        ar13 = GIC2GPathfinder.rg_array(a13[0])

        a14 = IC.get_aice([dt.datetime(year,12,14)],verbos=True)
        ar14 = GIC2GPathfinder.rg_array(a14[0])

        a15 = IC.get_aice([dt.datetime(year,12,15)],verbos=True)
        ar15 = GIC2GPathfinder.rg_array(a15[0])

        a16 = IC.get_aice([dt.datetime(year,12,16)],verbos=True)
        ar16 = GIC2GPathfinder.rg_array(a16[0])

        a17 = IC.get_aice([dt.datetime(year,12,17)],verbos=True)
        ar17 = GIC2GPathfinder.rg_array(a17[0])

        a18 = IC.get_aice([dt.datetime(year,12,18)],verbos=True)
        ar18 = GIC2GPathfinder.rg_array(a18[0])

        a19 = IC.get_aice([dt.datetime(year,12,19)],verbos=True)
        ar19 = GIC2GPathfinder.rg_array(a19[0])

        a20 = IC.get_aice([dt.datetime(year,12,20)],verbos=True)
        ar20 = GIC2GPathfinder.rg_array(a20[0])

        a21 = IC.get_aice([dt.datetime(year,12,21)],verbos=True)
        ar21 = GIC2GPathfinder.rg_array(a21[0])

        a22 = IC.get_aice([dt.datetime(year,12,22)],verbos=True)
        ar22 = GIC2GPathfinder.rg_array(a22[0])

        a23 = IC.get_aice([dt.datetime(year,12,23)],verbos=True)
        ar23 = GIC2GPathfinder.rg_array(a23[0])

        a24 = IC.get_aice([dt.datetime(year,12,24)],verbos=True)
        ar24 = GIC2GPathfinder.rg_array(a24[0])

        a25 = IC.get_aice([dt.datetime(year,12,25)],verbos=True)
        ar25 = GIC2GPathfinder.rg_array(a25[0])

        a26 = IC.get_aice([dt.datetime(year,12,26)],verbos=True)
        ar26 = GIC2GPathfinder.rg_array(a26[0])

        a27 = IC.get_aice([dt.datetime(year,12,27)],verbos=True)
        ar27 = GIC2GPathfinder.rg_array(a27[0])

        a28 = IC.get_aice([dt.datetime(year,12,28)],verbos=True)
        ar28 = GIC2GPathfinder.rg_array(a28[0])

        a29 = IC.get_aice([dt.datetime(year,12,29)],verbos=True)
        ar29 = GIC2GPathfinder.rg_array(a29[0])

        a30 = IC.get_aice([dt.datetime(year,12,30)],verbos=True)
        ar30 = GIC2GPathfinder.rg_array(a30[0])

        a31 = IC.get_aice([dt.datetime(year,12,31)],verbos=True)
        ar31 = GIC2GPathfinder.rg_array(a31[0])

        dec = np.nanmean( np.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,
                                 a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31]), axis=0 )
        dec_r = np.nanmean( np.array([ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10,ar11,ar12,ar13,ar14,ar15,ar16,ar17,
                                   ar18,ar19,ar20,ar21,ar22,ar23,ar24,ar25,ar26,ar27,ar28,ar29,ar30,ar31]), axis=0 )
        
        return dec, dec_r
