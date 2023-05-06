# functions for playing Ekman pumping
# dependancies
import numpy as np
# from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
from scipy.interpolate import griddata
from glob import glob
# from numba import jit
# from numba import :
import datetime
from dateutil.relativedelta import relativedelta
import parameters as par

# constants
rhoa = par.rho_a  # Air density
rhoi = par.rho_i  # Ice density
rhoo = par.rho_o  # Ocean density

# going to make a class that holds all the data
# first it will be defined by start/finish datetime objects
# this will depend on year/month only for the limits but compatible with more times for getting the data.
# will have different def to create the ekman pumping
# will have a flag to show if it has been calculated.
# want to pump.forcingcomp(several names) = np.arrays
# want pump.results[iterable].components(several names) = np.arrays
# will now load an individual ppump_result into a list (maybe)
class ppump_forcing:

    def __init__(self, start, finish, period):
        self.start = start
        self.finish = finish
        self.period = period
        if period.months > 0:
            self.n_t = int(((finish.month - start.month)
                            + (finish.year - start.year) * 12) / period.months)
        elif period.days > 0:
            self.n_t = int((finish - start).days / period.days)
        else:
            print("Period of months or days ONLY, currently supported")
        self.n_yrs = finish.year - start.year + 1
        self.forcing = False
        self.grid = False

    # grid will be in a m/n lon/lat (raidians) np array
    def load_grid_from_np(self, file, lon_name='lonsG', lat_name='latsG'):
        # loads the grid first
        if self.grid:
            print('Grid exists, check it first or regrid')
            return
        else:
            npzfile = np.load(file)
            lonsG = npzfile[lon_name]
            latsG = npzfile[lat_name]
            # somehow save the projection as a dict then load it
            # self.proj =
            npzfile.close()
            self.lons = lonsG
            self.lats = latsG
            m, n = np.shape(latsG)
            self.m = m
            self.n = n
            self.grid = True

    def load_forcing_from_np(self, file):
        # loads the data in default format. For saved data only
        if self.forcing:
            print('Forcing exists, check it first')
            return
        elif self.grid:
            print('Grid exists, check it first or regrid')
            return
        else:
            npzfile = np.load(file)
            self.uax = npzfile['uax']
            self.uay = npzfile['uay']
            self.uix = npzfile['uix']
            self.uiy = npzfile['uiy']
            self.uox = npzfile['uox']
            self.uoy = npzfile['uoy']
            lons = npzfile["lons"]
            lats = npzfile["lats"]
            npzfile.close()

            self.lons = lons
            self.lats = lats
            m, n = np.shape(lats)
            self.m = m
            self.n = n
            self.grid = True


# This function inputs the results from the class ppump_forcing above, which includes the Pathfinder grid.
# It returns the grid and period of forcing. As far as I understand, this function turns the class into an instance.

# Olivia's comments:
# now need to regrid some forcing 
# takes the lon lat of
def start_pump(ppump_forcing, copy_grid=False):
    # inits a ppump_results from a forcing
    # need start,finish,period
    start = ppump_forcing.start
    finish = ppump_forcing.finish
    period = ppump_forcing.period
    pr = ppump_results(start, finish, period)
    if ppump_results.grid & copy_grid:
        pr.lats = ppump_forcing.lats
        pr.lons = ppump_forcing.lons
        pr.m = ppump_forcing.m
        pr.n = ppump_forcing.n
        pr.grid = True
        print('Copied grid size ', m, n, ' to results')
    return pr


class ppump_results:
    # saved pumps - want to easily save and load a pump first check it's consistent with lat lon
    # for plotting - want to easily loop over the pumps loaded, with their names
    # change now - single time point
    # set by time_start and period

    def __init__(self, start, finish, period):
        self.start = start
        self.finish = finish
        self.period = period
        if period.months > 0:
            self.n_t = int(((finish.month - start.month)
                            + (finish.year - start.year) * 12) / period.months)
        elif period.days > 0:
            self.n_t = int((finish - start).days / period.days)
        else:
            print("Period of months or days ONLY, currently supported")
        self.n_yrs = finish.year - start.year + 1
        self.grid = False
        self.fcor = False

    # grid will be in a m/n lon/lat (raidians) np array
    def load_grid_from_np(self, file, lon_name='lonsG', lat_name='latsG'):
        # loads the grid first
        if self.grid:
            print('Grid exists, check it first or regrid')
            return
        else:
            with np.load(file) as npzfile:
                lonsG = npzfile[lon_name]
                latsG = npzfile[lat_name]
            # somehow save the projection as a dict then load it
            # self.proj =
            self.lons = lonsG
            self.lats = latsG
            m, n = np.shape(latsG)  # I think m is latitude and n is longitude for regridding?
            self.m = m
            self.n = n
            self.grid = True

    def set_proj(self, mplot):
        # puts in a projection mplot too
        #        # make make sure you keep hold of regridding projs
        self.mplot = mplot
        self.proj = True

    def set_grid_dxy(self, dxRes, dyRes):
        # creates a grid depending on wanted resolution
        if self.proj:
            nx = int((self.mplot.xmax - self.mplot.xmin) / dxRes) + 1
            ny = int((self.mplot.ymax - self.mplot.ymin) / dyRes) + 1
            lons, lats, xpts, ypts = self.mplot.makegrid(nx, ny, returnxy=True)
            self.lons = lons
            self.lats = lats
            self.xpts = xpts
            self.ypts = ypts
            self.dxRes = dxRes
            self.dyRes = dyRes
            self.grid = True
            self.m = nx
            self.n = ny
            print("Got a grid res = ", nx, " x ", ny)
            print("Note that all grid info is in ny x nx grids, whilst data is in nx x ny")
        else:
            print("Projection not defined yet, do that first")

    def set_grid_mn(self, nx, ny):
        # creates a grid depending on wanted no. of points
        if self.proj:
            lons, lats, xpts, ypts = self.mplot.makegrid(nx, ny, returnxy=True)
            self.lons = lons
            self.lats = lats
            self.xpts = xpts
            self.ypts = ypts
            self.grid = True
            self.dxRes = (self.mplot.xmax - self.mplot.xmin) / (nx - 1)
            self.dyRes = (self.mplot.ymax - self.mplot.ymin) / (ny - 1)
            self.m = nx
            self.n = ny
            print("Got a grid res = ", nx, " x ", ny)
        else:
            print("Projection not defined yet, do that first")

    def get_grid_info(self):
        # creates a grid depending on wanted no. of points
        # print( self.grid and (not self.gridinfo))
        if self.grid and (not self.gridinfo):
            # iterate over the grid to get dimensions and angles
            # first iterate all x dimensions - m-1/n array
            # then  iterate all y dimensions - m/n-1 array
            xdims = np.empty([self.n, self.m - 1])
            ydims = np.empty([self.n - 1, self.m])
            self.xdist = np.empty([self.n, self.m])
            self.ydist = np.empty([self.n, self.m])
            self.ang_c = np.empty([self.n, self.m])
            self.ang_s = np.empty([self.n, self.m])
            for i in range(self.m - 1):
                for j in range(self.n):
                    xdims[j, i] = ellipsoidal_distance(
                        self.lons[j, i], self.lats[j, i],
                        self.lons[j, i + 1], self.lats[j, i + 1], deg=True)
            for i in range(self.m):
                for j in range(self.n - 1):
                    ydims[j, i] = ellipsoidal_distance(
                        self.lons[j, i], self.lats[j, i],
                        self.lons[j + 1, i], self.lats[j + 1, i], deg=True)

            # then average the available distances i-1,i j-1,j
            for i in range(self.m):
                for j in range(self.n):
                    self.xdist[j, i] = np.nanmean(xdims[j, :i + 1][-2:])
                    self.ydist[j, i] = np.nanmean(ydims[:j + 1, i][-2:])
            print("Grid distances calculated: ", np.nanmean(self.xdist), " x ", np.nanmean(self.ydist))

            # then  iterate all angles - this is all points plus the extra possible angles
            # pad the lon lat arrays for iteration
            lon_pad = np.pad(self.lons, (1, 1), 'linear_ramp', end_values=(np.nan))
            lat_pad = np.pad(self.lats, (1, 1), 'linear_ramp', end_values=(np.nan))
            for i in range(self.m):
                for j in range(self.n):
                    # i + angle
                    xPlus_c, xPlus_s = lon_lat_angle(lon_pad[j + 1, i + 1], lat_pad[j + 1, i + 1],
                                                     lon_pad[j + 1, i + 2], lat_pad[j + 1, i + 2], return_trig=True,
                                                     deg=True)
                    xMins_c, xMins_s = lon_lat_angle(lon_pad[j + 1, i + 1], lat_pad[j + 1, i + 1],
                                                     lon_pad[j + 1, i], lat_pad[j + 1, i], return_trig=True, deg=True)
                    yPlus_c, yPlus_s = lon_lat_angle(lon_pad[j + 1, i + 1], lat_pad[j + 1, i + 1],
                                                     lon_pad[j + 2, i + 1], lat_pad[j + 2, i + 1], return_trig=True,
                                                     deg=True)
                    yMins_c, yMins_s = lon_lat_angle(lon_pad[j + 1, i + 1], lat_pad[j + 1, i + 1],
                                                     lon_pad[j, i + 1], lat_pad[j, i + 1], return_trig=True, deg=True)
                    # average all the components first checking the orientation
                    # if j == 20 and i ==12:
                    # print([xPlus_c,xMins_c,yPlus_c,yMins_c])
                    # print([xPlus_s,xMins_s,yPlus_s,yMins_s])
                    self.ang_c[j, i] = np.nanmean([-xPlus_s, xMins_s, yPlus_c, -yMins_c])
                    self.ang_s[j, i] = np.nanmean([xPlus_c, xMins_c, yPlus_s, -yMins_s])
            self.gridinfo = True
        else:
            print("Grid not defined yet, do that first")

    def save_grid(self, file):
        if self.grid and self.gridinfo:
            # save lat/lon pts 
            np.savez(file,
                     lats=self.lats,
                     lons=self.lons,
                     xpts=self.xpts,
                     ypts=self.ypts,
                     dxRes=self.dxRes,
                     dyRes=self.dyRes,
                     m=self.m,
                     n=self.n,
                     ang_c=self.ang_c,
                     ang_s=self.ang_s,
                     xdist=self.xdist,
                     ydist=self.ydist)
            print("Grid saved in " + file)
        else:
            print("No grid to save - run get_grid_info")

    def load_grid(self, file):  # Check which files get loaded
        with np.load(file) as npzfile:  # This closes the file after it's done
            self.lats = npzfile["lats"]
            self.lons = npzfile["lons"]
            self.xpts = npzfile["xpts"]
            self.ypts = npzfile["ypts"]
            self.dxRes = npzfile["dxRes"]
            self.dyRes = npzfile["dyRes"]
            self.m = npzfile["m"]
            self.n = npzfile["n"]
            self.ang_c = npzfile["ang_c"]
            self.ang_s = npzfile["ang_s"]
            self.xdist = npzfile["xdist"]
            self.ydist = npzfile["ydist"]
            self.grid = True
            self.gridinfo = True
            print("Loaded a grid: " + file)

    def check_grid(self):
        # makes sure the projection and loaded grid are consistent
        if self.proj and self.grid and self.gridinfo:
            proj_dim = self.mplot.xmax - self.mplot.xmin
            proj_dim = proj_dim / self.m
            print("Projection av xdim = ", proj_dim)
            print("dxRes              = ", self.dxRes)
            print("xdist av           = ", np.mean(self.xdist))

    def get_grid_mask(self, inflate=0.0):
        # makes a land mask for each point then inflates by a distance m
        # makes a land mask for each point then inflates by a distance m
        if self.masked:
            print("Already masked, do it again? set mask = False first")
        else:
            self.mask = np.ones([self.m, self.n])
            for i in range(self.m):  # going through each latitude
                for j in range(self.n):  # going through each longitude
                    if self.mplot.is_land(self.xpts[j, i], self.ypts[j, i]):
                        self.mask[i, j] = np.nan  # sets land values to nan
            inf_mask = np.ones([self.m, self.n])
            if (inflate > 0.0) and self.gridinfo:  # As far as I can see, there's no way to make inflate >0?
                self.mask_inflate = inflate
                for i in range(self.m):
                    for j in range(self.n):
                        if np.isnan(self.mask[i, j]):
                            inf_p = int(inflate / np.hypot(self.xdist[j, i], self.ydist[j, i]))
                            inf_mask[i - inf_p:i + inf_p + 1, j - inf_p:j + inf_p + 1] = np.nan
                self.mask = inf_mask
        self.masked = True

    #    def load_pump_from_np(self,file,name):
    #        if self.grid:
    #            print('Grid exists, check it first or regrid')
    #        else:
    def generate_dates(self):
        # makes an indexable array of datetimes
        # makes index for slecting across years/months
        # check what form to make selecting these indicies easier
        if ~self.dates:
            n_mnth_s = []
            n_mnth_f = []

        if self.yr_mnth:
            n_mth_s.append(self.start.month - 1)
            n_mth_f.append(12)
            yr = 1
            while yr < self.n_yrs - 1:
                n_mth_s.append(0)
                n_mth_f.append(12)
                yr += 1
            n_mth_s.append(0)
            n_mth_f.append(self.finish.month)
        else:
            n_mth_s.append(0)
            n_mth_f.append((self.finish.month - self.start.month)
                           + (self.finish.year - self.start.year) * 12)
        print(self.n_yrs, pf.n_yrs, n_mth_s, n_mnth_f)


def check_bounds(pump_list):
    d = np.shape(pump_list)
    return d


# takes a pumpp_forcing and generates a ppump_results from the given eckman_pumping and surface_stress
# def calc_ekman(self,forcing,name,stress,ekman):
# 

# 

# this class defines a pumping, the equations will do everything needed for each pumping - 
# eckman current from stresses
# eckman stress from current
# test to see if they return correct values
# test to confirm each is consistent
# when calling the function do i want to check the m is correct shape?
class eckman_pumping:

    def __init__(self, current_func, stress_func):
        self.current_func = current_func
        self.stress_func = stress_func

    def eck_currents(self, taux, tauy, fcor, m):
        ue, ve = self.current_func(taux, tauy, fcor, m)
        return ue, ve

    def eck_stress(self, ue, ve, fcor, m):
        taux, tauy = self.stress_func(ue, ve, fcor, m)
        return taux, tauy

    # now check that the two functions are coherent and return result
    def check_funcs(self, ue, ve, fcor, m, eps):
        taux, tauy = self.stress_func(ue, ve, fcor, m)
        ue_test, ve_test = self.current_func(taux, tauy, fcor, m)
        if (np.abs(ue_test - ue) < eps) and (np.abs(ve_test - ve) < eps):
            return True
        else:
            return False


# this class defines a surface stress type, the defs will give everything we need from it.
# open ocean stress (wnd only)
# total stress for ice - from wnd ice eck 
# test to see if they return correct values
# class surface_stress:
# 
#     def __init__(self,current_func,stress_func):
#         self.current_func = current_func
#         self.stress_func = stress_func
# 
#     
#     def eck_currents(self,taux,tauy,fcor,m):
#         ue, ve = self.current_func(taux,tauy,fcor,m)
#         return ue, ve
# 
#     def eck_stress(self,ue, ve ,fcor,m):
#         taux,tauy = self.stress_func(ue, ve ,fcor,m)
# 	return taux,tauy
# 
# 	# now check that the two functions are coherent and return result
#     def check_funcs(self,ue, ve ,fcor,m,eps):
#         taux,tauy = self.stress_func(ue, ve ,fcor,m)
#         ue_test, ve_test = self.current_func(taux,tauy,fcor,m)
# 	if (np.abs(ue_test - ue) < eps) and (np.abs(ve_test - ve) < eps):
#             return True
# 	else:
#             return False


# now need something that takes an eckman pumping and 

# models of the form
#    Ca
#    Co
#    thA
#    thO
#    eck_rot

# takes the input vectors given, 
# returns the stress components,
# different function for each stress form considered

# the stress outputs to the ocean are:
# from the air tau_a
# from the ice tau_i
# from the ice with still ocean tau_i0
# something clever with eckman currents?

# this is where calc_pump went.

def taus_quad_open(m, uwnd, vwnd):
    wnd_mag = np.hypot(uwnd, vwnd)  # calculating wind vector magnitude

    tau_aox = m[0] * rhoa * wnd_mag * (uwnd * np.cos(m[2]) - vwnd * np.sin(m[2]))
    tau_aoy = m[0] * rhoa * wnd_mag * (vwnd * np.cos(m[2]) + uwnd * np.sin(m[2]))

    # Tau_ao is the surface stress between air and ocean. It takes the (m), air density (rhoa), wind magnitude (wnd_mag)

    return tau_aox, tau_aoy


# @jit
def taus_quad(m, uwnd, vwnd, uice, vice, ueck, veck, uocn, vocn, conc):
    ui_o = uice - uocn - ueck
    vi_o = vice - vocn - veck

    wnd_mag = np.hypot(uwnd, vwnd)
    i_o_mag = np.hypot(ui_o, vi_o)

    tau_aox = (1.0 - conc) * m[0] * rhoa * wnd_mag * (uwnd * np.cos(m[2]) - vwnd * np.sin(m[2]))
    tau_aoy = (1.0 - conc) * m[0] * rhoa * wnd_mag * (vwnd * np.cos(m[2]) + uwnd * np.sin(m[2]))

    tau_iox = conc * m[1] * rhoo * i_o_mag * (ui_o * np.cos(m[3]) - vi_o * np.sin(m[3]))
    tau_ioy = conc * m[1] * rhoo * i_o_mag * (vi_o * np.cos(m[3]) + ui_o * np.sin(m[3]))

    return tau_aox, tau_aoy, tau_iox, tau_ioy


# As far as I can tell, m[0] is the air-ocean drag coefficient, m[1] is the ic-ocean drag coefficient, m[2] is phi (), and m[3] is beta ()

# @jit
def ueck_quad(taux, tauy, fcor, ek_depth):
    #  note here that tau is 2 element vector
    mult = np.sqrt(2) / fcor / rhoo / ek_depth

    # note for here the default value of
    if par.HEMI == "north":
        eck_rot = -np.pi / 4
    elif par.HEMI == "south":
        eck_rot = np.pi / 4
    else:
        raise ValueError(f"Invalid hemisphere {par.HEMI}")

    ueck = mult * (taux * np.cos(eck_rot) - tauy * np.sin(eck_rot))
    veck = mult * (tauy * np.cos(eck_rot) + taux * np.sin(eck_rot))

    return ueck, veck


# this is the function to minimise in order to derive ueck
# @jit
def ueck_solve_quad(m, uwnd, vwnd, uice, vice, ueck, veck, uocn, vocn, conc, fcor, ek_depth):
    ui_o = uice - uocn - ueck  # ui_o and vi_o are relative veolcity
    vi_o = vice - vocn - veck  # between ocean and ice current speeds

    wnd_mag = np.hypot(uwnd, vwnd)  # wind speed
    i_o_mag = np.hypot(ui_o, vi_o)  # ice-ocean relative speed

    tau_aox = (1.0 - conc) * m[0] * rhoa * wnd_mag * (uwnd * np.cos(m[2]) - vwnd * np.sin(m[2]))  # air-ocean stress
    tau_aoy = (1.0 - conc) * m[0] * rhoa * wnd_mag * (vwnd * np.cos(m[2]) + uwnd * np.sin(m[2]))

    tau_iox = conc * m[1] * rhoo * i_o_mag * (ui_o * np.cos(m[3]) - vi_o * np.sin(m[3]))  # ice-ocean stress
    tau_ioy = conc * m[1] * rhoo * i_o_mag * (vi_o * np.cos(m[3]) + ui_o * np.sin(m[3]))

    mult = fcor * rhoo * ek_depth / np.sqrt(2)  # multiplier for ekman values

    if par.MODEL == "model2":  # Model 2 uses 90 degree turning angle (Yang et al., 2006)
        if par.HEMI == "north":
            turning_angle = np.pi/2
        if par.HEMI == "south":
            turning_angle = -np.pi/2
    else:  # Other models use 45 degree turning angle (Meneghello et al., 2018)
        if par.HEMI == "north":
            turning_angle = np.pi/4
        if par.HEMI == "south":
            turning_angle = -np.pi/4

    if par.HEMI == "north":
        tau_eckx = mult * (ueck * np.cos(turning_angle) - veck * np.sin(turning_angle))
        tau_ecky = mult * (veck * np.cos(turning_angle) + ueck * np.sin(turning_angle))
    elif par.HEMI == "south":
        tau_eckx = mult * (ueck * np.cos(-turning_angle) - veck * np.sin(-turning_angle))
        tau_ecky = mult * (veck * np.cos(-turning_angle) + ueck * np.sin(-turning_angle))
    else:
        raise ValueError(f"Invalid hemisphere {par.HEMI}")
        #This makes the turning angle correspond to the hemisphere

    tau_minx = tau_aox + tau_iox - tau_eckx
    tau_miny = tau_aoy + tau_ioy - tau_ecky

    return tau_minx, tau_miny


# @jit
def fb_residual_quad(m, uwnd, vwnd, uice, vice, ueck, veck, uocn, vocn, conc, hi, fcor):
    uo_i = uocn + ueck - uice
    vo_i = vocn + veck - vice

    wnd_mag = np.hypot(uwnd, vwnd)
    o_i_mag = np.hypot(uo_i, vo_i)

    tau_aix = conc * m[0] * rhoa * wnd_mag * (uwnd * np.cos(m[2]) - vwnd * np.sin(m[2]))
    tau_aiy = conc * m[0] * rhoa * wnd_mag * (vwnd * np.cos(m[2]) + uwnd * np.sin(m[2]))

    tau_oix = conc * m[1] * rhoo * o_i_mag * (uo_i * np.cos(m[3]) - vo_i * np.sin(m[3]))
    tau_oiy = conc * m[1] * rhoo * o_i_mag * (vo_i * np.cos(m[3]) + uo_i * np.sin(m[3]))

    cor_tltx = -rhoi * hi * fcor * vo_i
    cor_tlty = rhoi * hi * fcor * uo_i

    fb_r_x = tau_aix + tau_oix + cor_tltx
    fb_r_y = tau_aiy + tau_oiy + cor_tlty

    return fb_r_x, fb_r_y


# @jit
def calcCurlSq2dXYGradient(x_vel, y_vel, dx_res):
    # thanks Alek Petty!!!!
    # calculate the curl (squared) of a field given some input vector grid
    # CALCULATE THE CURL OF AN X/Y VEXTOR FIELD. DX_RES IF THE GRID RESOLUTION OF THIS REGULARLY SPACED GRID.
    # MULTIPLY BY MAGNITUDE TO GET THIS IN SQUARED UNITS ANALAGOUS TO THE WIND STRESS CURL.
    # USE GRADIENT FUNCTION WHICH USES CENTRAL DIFFERENCES IN MIDDLE CELLS AND FIRST DIFFERENCES AT THE BOUNDARIES (GIVES SAME SHAPE AS INPUT)

    mag = np.sqrt((x_vel ** 2) + (y_vel ** 2))

    x_vel_mag = x_vel * mag
    y_vel_mag = y_vel * mag
    x_vel_mag = x_vel
    y_vel_mag = y_vel

    # gradient [0] returns row divs (y direction) then [1] gives column divs (x direction)
    dvelydx = np.gradient(y_vel_mag, dx_res)[1]
    dvelxdy = np.gradient(x_vel_mag, dx_res)[0]

    zeta = dvelydx - dvelxdy
    # MASK ARRAY WHERE VALUES ARE NAN
    zeta = ma.array(zeta, mask=np.isnan(zeta))

    return zeta


# @jit
def gridvectors(m, xvel, yvel, lon, xptsG, yptsG, lonsG, latsG, xptsM, yptsM):
    # Comes in xy coordinates so need to rotate to UV
    alpha = lon * np.pi / 180.
    uvel = yvel * np.sin(alpha) + xvel * np.cos(alpha)
    vvel = yvel * np.cos(alpha) - xvel * np.sin(alpha)

    # Re-grid data
    uvelG = griddata((xptsM, yptsM), uvel, (xptsG, yptsG), method='linear')
    vvelG = griddata((xptsM, yptsM), vvel, (xptsG, yptsG), method='linear')

    # Rotate data onto new grid
    xvelG, yvelG = m.rotate_vector(uvelG, vvelG, lonsG, latsG)
    xvelG = ma.masked_invalid(xvelG)
    yvelG = ma.masked_invalid(yvelG)

    return xvelG, yvelG


# @jit(float64[:,:],float64[:,:],cache = True,nopython = True)
# @jit(cache = True,nopython = True)
# @jit
def get_vec_month_from_list(files, nm):
    file = files[0]
    f = np.genfromtxt(file)
    xvel = np.empty(np.shape(f)[0])
    yvel = np.empty(np.shape(f)[0])
    xvel[:] = f[:, 0]
    yvel[:] = f[:, 1]
    # xvel = f[:,0]
    # yvel = f[:,1]
    for i in range(1, nm):
        file = files[i]
        f = np.genfromtxt(file)
        f[np.isnan(f)] = 0
        xvel[:] += f[:, 0]
        yvel[:] += f[:, 1]
    xvel = xvel / nm
    yvel = yvel / nm
    return xvel, yvel


def ellipsoidal_distance(lat1, long1, lat2, long2):
    # Note from Alessandra: Sorry to whomever inherits this code.
    # I also inherited this and couldn't be bothered to change all the magic numbers.
    # This is completely illegible, but it seems to work.

    a = 6378137.0  # equatorial radius in meters
    f = 1 / 298.257223563  # ellipsoid flattening
    b = (1 - f) * a
    tolerance = 1e-11  # to stop iteration

    phi1, phi2 = lat1, lat2
    U1 = np.arctan((1 - f) * np.tan(phi1))
    U2 = np.arctan((1 - f) * np.tan(phi2))
    L1, L2 = long1, long2
    L = L2 - L1

    lambda_old = L + 0

    while True:

        t = (np.cos(U2) * np.sin(lambda_old)) ** 2
        t += (np.cos(U1) * np.sin(U2) - np.sin(U1) * np.cos(U2) * np.cos(lambda_old)) ** 2
        sin_sigma = t ** 0.5
        cos_sigma = np.sin(U1) * np.sin(U2) + np.cos(U1) * np.cos(U2) * np.cos(lambda_old)
        sigma = np.arctan2(sin_sigma, cos_sigma)

        sin_alpha = np.cos(U1) * np.cos(U2) * np.sin(lambda_old) / sin_sigma
        cos_sq_alpha = 1 - sin_alpha ** 2
        cos_2sigma_m = cos_sigma - 2 * np.sin(U1) * np.sin(U2) / cos_sq_alpha
        C = f * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha)) / 16

        t = sigma + C * sin_sigma * (cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m ** 2))
        lambda_new = L + (1 - C) * f * sin_alpha * t
        if np.abs(lambda_new - lambda_old) <= tolerance:
            break
        else:
            lambda_old = lambda_new

    u2 = cos_sq_alpha * ((a ** 2 - b ** 2) / b ** 2)
    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    t = cos_2sigma_m + 0.25 * B * (cos_sigma * (-1 + 2 * cos_2sigma_m ** 2))
    t -= (B / 6) * cos_2sigma_m * (-3 + 4 * sin_sigma ** 2) * (-3 + 4 * cos_2sigma_m ** 2)
    delta_sigma = B * sin_sigma * t
    s = b * A * (sigma - delta_sigma)

    return s

# %%
