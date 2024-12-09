'''

_get_geq.py is a module managing quantities from the geqdsk

cjperks

06/01/23

'''

import numpy as np
from scipy.interpolate import interp1d
from omfit_classes import omfit_eqdsk
import matplotlib.pyplot as plt
from portals.gs_tools import GEQmodule
import scipy.constants as cnt
import os

__all__ = [
    'get_geq',
    'get_harm_res',
]

##########################################################
#
#                     Main
#
##########################################################

def get_geq(
    path_input = None,
    path_gfile = None,
    plt_all = None,
    ):

    dout = {}
    dout['paths'] = {}
    dout['paths']['input'] = path_input
    dout['paths']['gfile'] = path_gfile

    # sq. norm. tor. flux basis
    dout['rhot'] = np.linspace(0,1,101)

    # Output data from gfile
    gfile = omfit_eqdsk.OMFITgeqdsk(
        os.path.join(path_input, path_gfile)
        )

    # gfile scalar values
    torfluxa_g   = -1*gfile['fluxSurfaces']['geo']['phi'][-1]/(2*np.pi) # [Wb/rad], edge toroidal flux divided by 2pi (Note: direction dependent)
    rcentr_g     = gfile['RCENTR']          # [m], radius of the plasma center
    bcentr_g     = -1*gfile['BCENTR']      # [T], on-axis B-field (Note: direction dependent)
    current_g    = -1*gfile['CURRENT']/1e6  # [MA], plasma current (Note: direction dependent)

    # gfile geometry profiles
    psin_g = gfile['fluxSurfaces']['geo']['psin']    # radial flux coordinate mesh (norm. pol flux)
    sq_phin_g = np.sqrt(gfile['fluxSurfaces']['geo']['phi']/gfile['fluxSurfaces']['geo']['phi'][-1]) # sqrt of norm. toroidal flux
    rmin_g = gfile['fluxSurfaces']['geo']['a']                # [m], minor radius
    polflux_g = -1*(gfile['AuxQuantities']['PSI']-gfile['AuxQuantities']['PSI'][0])  # [Wb/rad] poloidal flux divided by 2pi
    q_g = gfile['QPSI']                                # [], q profile
    rmaj_g = gfile['fluxSurfaces']['geo']['R']         # [m], major radius
    zmag_g = gfile['fluxSurfaces']['geo']['Z']         # [m] height
    kappa_g = gfile['fluxSurfaces']['geo']['kap']      # [], elongation
    delta_g = gfile['fluxSurfaces']['geo']['delta']    # [], triangularity
    zeta_g = gfile['fluxSurfaces']['geo']['zeta']      # [], squareness
    johm_g = gfile['fluxSurfaces']['avg']['Jt']/1e6    # [MA/m^2], ohmic current density
    vol_g = gfile['fluxSurfaces']['geo']['vol']        # [m^3], volume of each flux surface
    Bt_g = gfile['fluxSurfaces']['avg']['Bt']          # [T], toroidal magnetic field

    shape_cos0_g = np.zeros(np.size(q_g)); # [], higher flux surface tilt moment
    shape_cos1_g = np.zeros(np.size(q_g)); # [], higher flux surface tilt moment
    shape_cos2_g = np.zeros(np.size(q_g)); # [], higher flux surface tilt moment
    shape_cos3_g = np.zeros(np.size(q_g)); # [], higher flux surface tilt moment
    shape_sin3_g = np.zeros(np.size(q_g)); # [], higher flux surface tilt moment
    jbs_g = np.zeros(np.size(q_g));        # [MA/m^2], parallel bootstrap current density
    jbstor_g = np.zeros(np.size(q_g));     # [MA/m^2], toroidal bootstrap current density
    jrf_g = np.zeros(np.size(q_g));        # [MA/m^2], RF-driven current density
    jnb_g = np.zeros(np.size(q_g));        # [MA/m^2], beam-driven current density

    # Interpolates gfile quantity into GACODE radial flux coordinate, \sqrt(\phi)
    dout['rmin_m'] = interp1d(sq_phin_g, rmin_g)(dout['rhot'])
    dout['polflux_Wb/rad'] = interp1d(sq_phin_g, polflux_g)(dout['rhot'])
    dout['q'] = interp1d(sq_phin_g, q_g)(dout['rhot'])
    dout['rmaj_m'] = interp1d(sq_phin_g, rmaj_g)(dout['rhot'])
    dout['zmag_m'] = interp1d(sq_phin_g, zmag_g)(dout['rhot'])
    dout['kappa'] = interp1d(sq_phin_g, kappa_g)(dout['rhot'])
    dout['delta'] = interp1d(sq_phin_g, delta_g)(dout['rhot'])
    dout['zeta'] = interp1d(sq_phin_g, zeta_g)(dout['rhot'])
    dout['shape_cos0'] = interp1d(sq_phin_g, shape_cos0_g)(dout['rhot'])
    dout['shape_cos1'] = interp1d(sq_phin_g, shape_cos1_g)(dout['rhot'])
    dout['shape_cos2'] = interp1d(sq_phin_g, shape_cos2_g)(dout['rhot'])
    dout['shape_cos3'] = interp1d(sq_phin_g, shape_cos3_g)(dout['rhot'])
    dout['shape_sin3'] = interp1d(sq_phin_g, shape_sin3_g)(dout['rhot'])
    dout['johm_MA/m2'] = interp1d(sq_phin_g, johm_g)(dout['rhot'])
    dout['jbs_MA/m2'] = interp1d(sq_phin_g, jbs_g)(dout['rhot'])
    dout['jbstor_MA/m2'] = interp1d(sq_phin_g, jbstor_g)(dout['rhot'])
    dout['jrf_MA/m2'] = interp1d(sq_phin_g, jrf_g)(dout['rhot'])
    dout['jnb_MA/m2'] = interp1d(sq_phin_g, jnb_g)(dout['rhot'])
    dout['vol_m3'] = interp1d(sq_phin_g, vol_g)(dout['rhot'])
    dout['Bt_T'] = interp1d(sq_phin_g, Bt_g)(dout['rhot'])

    dout['rhop'] = np.sqrt(dout['polflux_Wb/rad']/dout['polflux_Wb/rad'][-1])

    # Store scalar value
    dout['torfluxa_Wb/rad'] = "{:1.7E}".format(torfluxa_g)
    dout['bcentr_T'] = "{:1.7E}".format(bcentr_g)
    dout['rcentr_m'] = " {:1.7E}".format(rcentr_g)
    dout['current_MA'] = "{:1.7E}".format(current_g)

    # Plotting
    if plt_all:
        g = GEQmodule.PORTALSgeqdsk(
            os.path.join(path_input, path_gfile)
            )
        g.plot()

        plt.show()

    return dout

##########################################################
#
#                     Utilities
#
##########################################################

# Finds the (R,Z) contour of a cyclotron resonance layer
def get_harm_res(
    ddata = None,
    # Files
    path_input = None,
    path_gfile = None,
    # Wave
    freq = 80e6, # [Hz]
    order = 2,
    # Species
    cs = 16,
    amu = 40,
    ):

    if ddata is None:
        if path_gfile is not None:
            # Output data from gfile
            gfile = omfit_eqdsk.OMFITgeqdsk(
                os.path.join(path_input, path_gfile)
                )

            # Mesh
            R = gfile['AuxQuantities']['R'] # [m]
            Z = gfile['AuxQuantities']['Z'] # [m]

            # B-field
            Bt = gfile['AuxQuantities']['Bt'] # [T]

    else:
        R = ddata['R'] # [m], dim(nR,)
        Z = ddata['Z'] # [m], dim(nZ,)
        Bt = ddata['Bt'] # [T], dim(nZ, nR)

    # Get rid of helicity sign convenctions
    if np.max(Bt) <0:
        Bt *= -1

    # Initializes output
    Zout = Z.copy()
    Rout = np.zeros(len(Zout))

    # Fine Mesh
    Rfine = np.linspace(min(R), max(R), 100)

    # Loop over height
    for ii in np.arange(len(Z)):
        # Fine B-field
        Bfine = interp1d(R,Bt[ii,:])(Rfine)

        # Index of harmonic resonance
        try:
            ind = np.nanargmin(abs(
                Bfine - (freq*2*np.pi/order)*(cnt.m_p/cnt.e) *amu/cs
                ))
        except:
            ind = 0

        # Radial position
        Rout[ii] = Rfine[ind]

    # Output, [m]
    return Rout, Zout