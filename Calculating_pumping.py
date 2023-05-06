import numpy as np
from sys import path

# path.insert(0, '/Users/H/WAVES/geo_data_group/')
path.insert(0, '/Users/hp/OneDrive - University College London/Year 4 UCL/Master\'s project/Masters_project/')
import grid_set as gs
import parameters as par
import Ekman_pumping as ep


# need Pathfinder grid points
GPathfinder = gs.grid_set(par.m)
GPathfinder.load_grid(par.ID_grid)

# m_use = ep.m_AIDJEX_0 #(has angles = 0)
m_use = par.m_AIDJEX  # check parameters.py to see what constants this includes
[Ca, Co] = m_use[0:2]  # this just gets the values of Ca and Co, can change to regular values if you want

ek_depth = par.D_E  # Ekman depth
omega = 2. * np.pi / 24. / 60. / 60.  # Earths rotation
fcor = 2. * omega * np.sin(np.deg2rad(GPathfinder.lats.T))  # coriolis force


# # setting up empty arrays with the size of GPathfinder
# ue = np.zeros([len(GPathfinder.xpts), len(GPathfinder.xpts), 2])
# ueg0 = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
# tau_a = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
# tau_i = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
# tau_i0 = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
# tau_g = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
# tau_all = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
#
# np_ice = 0
# np_ocn = 0


### FUNCTION THAT CALCULATES THE TAUS AND THE EKMAN CURRENTS
def current_old(uwp, vwp, upp, vpp, ugp, vgp, alpha):
    """

    :param uwp: Wind u
    :param vwp: Wind v
    :param upp: Ice drift u
    :param vpp: Ice drift v
    :param ugp: Geostrophic u
    :param vgp: Geostrophic v
    :param alpha: Ice concentration
    :return: ?
    """
    # uwp and vwp are wind, upp and vpp are drift, ugp and vgp are geostrophic, and alpha is conc

    # setting up empty arrays with the size of GPathfinder
    ue = np.zeros([len(GPathfinder.xpts), len(GPathfinder.xpts), 2])
    ueg0 = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_a = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_i = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_i0 = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_g = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_all = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])

    # np_ice = 0
    # np_ocn = 0

    for i in range(len(GPathfinder.xpts)):
        for j in range(len(GPathfinder.ypts)):

            # if model == "model1":
            # All variables included, turning angle = 45 degrees
            # if model == "model2":
            # All variables included, turning angle = 90 degrees
            # See ep.ueck_solve_quad
            if par.MODEL == "model3":
                # Wind only: no ice or geostrophic circulation
                ugp[i, j] = 0.0
                vgp[i, j] = 0.0
                alpha[i, j] = 0.0
            if par.MODEL == "model4":
                # Ice-free
                alpha[i, j] = 0.0
            if par.MODEL == "model5":
                # Wind-free: ice and ocean
                uwp[i, j] = 0.0
                vwp[i, j] = 0.0
            if par.MODEL == "model6":
                # Ocean only, no ice or wind (just geostrophic balance)
                uwp[i, j] = 0.0
                vwp[i, j] = 0.0
                alpha[i, j] = 0.0
            if par.MODEL == "model7":
                # Wind and ice (most complete model pre-2011 because no geostrophic data)
                ugp[i, j] = 0.0
                vgp[i, j] = 0.0

            # ICE-COVERED
            if np.isfinite(ugp[i, j]) and np.isfinite(vgp[i, j]) and np.isfinite(alpha[i, j]) and alpha[i, j] > 0.0001:

                # lambda is basically creating a small function, solving x = Ek_u and y = Ek_v
                g = lambda x, y: ep.ueck_solve_quad(m_use, uwp[i, j], vwp[i, j], upp[i, j], vpp[i, j], x, y,
                                                    ugp[i, j], vgp[i, j], alpha[i, j], fcor[i, j], ek_depth)

                x0 = g(0, 0)  # 0,0 are Ekman current u&v vels

                # x and y components of minimized pumping force
                ue[i, j, 0] = x0[0]
                ue[i, j, 1] = x0[1]

                # x and y components of air and ice stress
                tau_a[i, j, 0], tau_a[i, j, 1], tau_i[i, j, 0], tau_i[i, j, 1] = ep.taus_quad(
                    m_use, uwp[i, j], vwp[i, j], upp[i, j], vpp[i, j], ue[i, j, 0], ue[i, j, 1], ugp[i, j], vgp[i, j],
                    alpha[i, j])

                # solve again without geostrophic currents - to see the difference
                g = lambda x, y: ep.ueck_solve_quad(m_use, uwp[i, j], vwp[i, j], upp[i, j], vpp[i, j], x, y,
                                                    0.0, 0.0, alpha[i, j], fcor[i, j], ek_depth)
                x0 = g(0, 0)
                ueg0[i, j, 0] = x0[0]
                ueg0[i, j, 1] = x0[1]

                _, _, tau_i0[i, j, 0], tau_i0[i, j, 1] = ep.taus_quad(m_use, uwp[i, j], vwp[i, j],
                                                                      upp[i, j], vpp[i, j], ueg0[i, j, 0],
                                                                      ueg0[i, j, 1],
                                                                      0.0, 0.0, alpha[i, j])

                # x and y components of geostrophic current forcing and total forcing
                tau_g[i, j, 0] = tau_i[i, j, 0] - tau_i0[i, j, 0]
                tau_g[i, j, 1] = tau_i[i, j, 1] - tau_i0[i, j, 1]
                tau_all[i, j, 0] = tau_a[i, j, 0] + tau_i[i, j, 0]
                tau_all[i, j, 1] = tau_a[i, j, 1] + tau_i[i, j, 1]

            # OPEN OCEAN
            elif np.isfinite(ugp[i, j]) and np.isfinite(
                    vgp[i, j]):  # and np.isfinite(uwp[i,j]) and np.isfinite(vwp[i,j]):
                # np_ocn += 1

                # x and y components of air stress (= total stress)
                tau_a[i, j, 0], tau_a[i, j, 1] = ep.taus_quad_open(m_use, uwp[i, j], vwp[i, j])
                tau_all[i, j, :] = tau_a[i, j, :]

                ue[i, j, 0], ue[i, j, 1] = ep.ueck_quad(tau_all[i, j, 0], tau_all[i, j, 1], fcor[i, j], ek_depth)

                ueg0[i, j, :] = np.nan
                tau_i[i, j, :] = np.nan
                tau_i0[i, j, :] = np.nan
                tau_g[i, j, :] = np.nan


            # OVER LAND
            else:
                ue[i, j, :] = np.nan
                ueg0[i, j, :] = np.nan
                tau_all[i, j, :] = np.nan
                tau_a[i, j, :] = np.nan
                tau_i[i, j, :] = np.nan
                tau_i0[i, j, :] = np.nan
                tau_g[i, j, :] = np.nan

            # print(tau_all[i,j])

    # print("Pump calc complete,",np_ice,"ice points,",np_ocn,"ocn points")

    return ue, ueg0, tau_all, tau_a, tau_i, tau_i0, tau_g


def current(uwp, vwp, upp, vpp, ugp, vgp, alpha):
    """

    :param uwp: Wind u
    :param vwp: Wind v
    :param upp: Ice drift u
    :param vpp: Ice drift v
    :param ugp: Geostrophic u
    :param vgp: Geostrophic v
    :param alpha: Ice concentration
    :return: ?
    """
    # uwp and vwp are wind, upp and vpp are drift, ugp and vgp are geostrophic, and alpha is conc

    # setting up empty arrays with the size of GPathfinder

    ue_ice = np.zeros([len(GPathfinder.xpts), len(GPathfinder.xpts), 2])
    ueg0_ice = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_a_ice = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_i_ice = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_i0_ice = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])

    ue_open = np.zeros([len(GPathfinder.xpts), len(GPathfinder.xpts), 2])
    ueg0_open = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_a_open = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_i_open = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])
    tau_i0_open = np.zeros([len(GPathfinder.xpts), len(GPathfinder.ypts), 2])

    # if model == "model1":
    # All variables included, turning angle = 45 degrees
    # if model == "model2":
    # All variables included, turning angle = 90 degrees
    # See ep.ueck_solve_quad
    if par.MODEL == "model3":
        # Wind only: no ice or geostrophic circulation
        ugp[:] = 0.0
        vgp[:] = 0.0
        alpha[:] = 0.0
    if par.MODEL == "model4":
        # Ice-free
        alpha[:] = 0.0
    if par.MODEL == "model5":
        # Wind-free: ice and ocean
        uwp[:] = 0.0
        vwp[:] = 0.0
    if par.MODEL == "model6":
        # Ocean only, no ice or wind (just geostrophic balance)
        uwp[:] = 0.0
        vwp[:] = 0.0
        alpha[:] = 0.0
    if par.MODEL == "model7":
        # Wind and ice (most complete model pre-2011 because no geostrophic data)
        ugp[:] = 0.0
        vgp[:] = 0.0

    # ICE-COVERED
    land_mask = ~(np.isfinite(ugp) & np.isfinite(vgp))  # True where land
    ice_mask = np.isfinite(alpha) & (alpha > 0.0001)  # True where ice.
    valid_ice_mask = ice_mask & ~land_mask
    valid_ocean_mask = ~ice_mask & ~land_mask

    # solving x = Ek_u and y = Ek_v
    # x and y components of minimized pumping force
    ue_ice[..., 0],  ue_ice[..., 1] = ep.ueck_solve_quad(m_use, uwp, vwp, upp, vpp,
                                0, 0, ugp, vgp, alpha, fcor, ek_depth)

    # x and y components of air and ice stress
    tau_a_ice[..., 0], tau_a_ice[..., 1], tau_i_ice[..., 0], tau_i_ice[..., 1] = ep.taus_quad(
        m_use, uwp, vwp, upp, vpp, ue_ice[..., 0], ue_ice[..., 1],
        ugp, vgp, alpha)

    # solve again without geostrophic currents - to see the difference
    ueg0_ice[..., 0], ueg0_ice[..., 1] = ep.ueck_solve_quad(m_use, uwp, vwp, upp, vpp,
                                0, 0, 0.0, 0.0, alpha, fcor, ek_depth)

    _, _, tau_i0_ice[..., 0], tau_i0_ice[..., 1] = ep.taus_quad(m_use, uwp, vwp,
                                                                upp, vpp, ueg0_ice[..., 0],
                                                                ueg0_ice[..., 1],
                                                                0.0, 0.0, alpha)

    # geostrophic current forcing and total forcing
    tau_g_ice = tau_i_ice - tau_i0_ice
    tau_all_ice = tau_a_ice + tau_i_ice

    # OPEN OCEAN

    # x and y components of air stress (= total stress)
    tau_a_open[..., 0], tau_a_open[..., 1] = ep.taus_quad_open(m_use, uwp, vwp)
    tau_all_open = tau_a_open.copy()

    ue_open[..., 0], ue_open[..., 1] = ep.ueck_quad(tau_all_open[..., 0], tau_all_open[..., 1],
                                                    fcor, ek_depth)

    # OUTPUT ARRAYS
    ue = np.where(ice_mask[..., np.newaxis], ue_ice, ue_open)  # like "ue_ice if ice_mask else ue_open"
    ueg0 = np.where(ice_mask[..., np.newaxis], ueg0_ice, np.nan)
    tau_a = np.where(ice_mask[..., np.newaxis], tau_a_ice, tau_a_open)
    tau_i = np.where(ice_mask[..., np.newaxis], tau_i_ice, np.nan)
    tau_i0 = np.where(ice_mask[..., np.newaxis], tau_i0_ice, np.nan)
    tau_g = np.where(ice_mask[..., np.newaxis], tau_g_ice, np.nan)
    tau_all = np.where(ice_mask[..., np.newaxis], tau_all_ice, tau_all_open)

    # LAND MASKS
    ueg0[land_mask] = np.nan
    tau_i[land_mask] = np.nan
    tau_i0[land_mask] = np.nan
    tau_g[land_mask] = np.nan
    ue[land_mask] = np.nan
    tau_all[land_mask] = np.nan
    tau_a[land_mask] = np.nan

    # print(tau_all[i,j])

    # print("Pump calc complete,",np_ice,"ice points,",np_ocn,"ocn points")

    return ue, ueg0, tau_all, tau_a, tau_i, tau_i0, tau_g


### FUNCTION THAT CALCULATES THE PUMPING
def pump(tau_a, tau_i, tau_g, tau_all):
    pump_scale = 60 * 60 * 24 * 365  # m/year
    pump_a = - gs.geo_curl(tau_a[:, :, 0], tau_a[:, :, 1], GPathfinder) * pump_scale / 1027.5 / fcor
    pump_i = - gs.geo_curl(tau_i[:, :, 0], tau_i[:, :, 1], GPathfinder) * pump_scale / ep.rhoo / fcor
    pump_g = - gs.geo_curl(tau_g[:, :, 0], tau_g[:, :, 1], GPathfinder) * pump_scale / ep.rhoo / fcor
    pump_all = - gs.geo_curl(tau_all[:, :, 0], tau_all[:, :, 1], GPathfinder) * pump_scale / ep.rhoo / fcor

    return pump_a, pump_i, pump_g, pump_all


def weighted_pump(tau_all):
    pump_scale = 60 * 60 * 24 * 365  # m/year
    # pump_a = gs.geo_curl(tau_a[:,:,0],tau_a[:,:,1],GPathfinder)*pump_scale/1027.5/fcor
    # pump_i = gs.geo_curl(tau_i[:,:,0],tau_i[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    # pump_g = gs.geo_curl(tau_g[:,:,0],tau_g[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    pump_all = - gs.geo_curl(tau_all[:, :, 0], tau_all[:, :, 1], GPathfinder) * pump_scale / ep.rhoo / fcor
    # pump_all2 = pump_a + pump_i

    return pump_all
