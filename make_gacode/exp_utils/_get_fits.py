'''

_get_fits.py is a module managing fitted kinetic profiles

cjperks

06/01/23

'''

import numpy as np
from uncertainties import unumpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

__all__ = [
    'get_fits',
]

##########################################################
#
#                     Main
#
##########################################################

def get_fits(
    path_input = None,
    path_kin = None,
    t0 = None,
    dt = None,
    plt_all = None,
    ):
    dout = {}
    dout['paths'] = {}
    dout['paths']['input'] = path_input
    dout['paths']['kin'] = path_kin
    dout['t0_s'] = t0

    # Opens .npy file of saved profiles from QuickFit
    with open(path_input+path_kin, 'rb') as f: 
        rho_t = np.load(f) # norm. sq. tor. flux
        time = np.load(f)/1e3 # [s]
        ne_t = np.load(f, allow_pickle=True)/1e19 # dim(t, rho); [1e19 m^-3]
        Te_t = np.load(f, allow_pickle=True)/1e3 # dim(t, rho); [keV]
        Ti_t = np.load(f, allow_pickle=True)/1e3 # dim(t, rho); [keV]

        ne_t = unumpy.nominal_values(np.array(ne_t))
        Te_t = unumpy.nominal_values(np.array(Te_t))
        Ti_t = unumpy.nominal_values(np.array(Ti_t))

    # Finds the time slice of interest
    t_ind = np.where((time >= t0-dt) & (time <= t0+dt))

    # Takes kinetic profile data at time slice of interest
    ne = ne_t[t_ind].T # dim(rho, t); [1e19 m^-3]
    Te = Te_t[t_ind].T # dim(rho,t); [keV]
    Ti = Ti_t[t_ind].T # dim(rho,t); [keV]

    # Removes negative values from fits
    Te[Te<0.2] = 0.2
    Ti[Ti<0.2] = 0.2
    ne[ne<0.5] = 0.5

    # Interpolates onto rho = 0->1
    rho = np.linspace(0,1,101)
    ne = interp1d(rho_t, ne, axis = 0)(rho)
    Te = interp1d(rho_t, Te, axis = 0)(rho)
    Ti = interp1d(rho_t, Ti, axis = 0)(rho)

    # Takes the average in time
    dout['ne_19m3'] = np.mean(ne, axis=-1)
    dout['Te_keV'] = np.mean(Te, axis=-1)
    dout['Ti_keV'] = np.mean(Ti, axis=-1)
    dout['rhot'] = rho

    if plt_all:
        _plot(
            dout = dout,
            ne = ne,
            Te = Te,
            Ti = Ti,
            t0=t0,
            dt=dt,
            )

    return dout


##########################################################
#
#                     Extra
#
##########################################################

def _plot(
    dout=None,
    ne=None,
    Te=None,
    Ti=None,
    t0=None,
    dt=None,
    ):

    fig, ax = plt.subplots(1,3)

    # Plots ne
    if ne is not None:
        ax[0].plot(dout['rhot'], ne, label='fit')
    ax[0].plot(dout['rhot'], dout['ne_19m3'], 'k', label='avg')

    # Plots Te
    if Te is not None:
        ax[1].plot(dout['rhot'], Te)
    ax[1].plot(dout['rhot'], dout['Te_keV'], 'k')

    # Plots Ti
    if Ti is not None:
        ax[2].plot(dout['rhot'], Ti)
    ax[2].plot(dout['rhot'], dout['Ti_keV'], 'k')

    ax[0].set_xlabel(r'$\rho_p$')
    ax[0].set_ylabel(r'$n_e$ [1e19 $m^{-3}$]')
    ax[0].grid('on')
    leg = ax[0].legend()
    leg.set_draggable('on')

    ax[1].set_xlabel(r'$\rho_p$')
    ax[1].set_ylabel(r'$T_e$ [$keV$]')
    ax[1].grid('on')

    ax[2].set_xlabel(r'$\rho_p$')
    ax[2].set_ylabel(r'$T_i$ [$keV$]')
    ax[2].grid('on')

    fig.suptitle('Fitted Profiles t='+str(t0-dt)+'-'+str(t0+dt)+'s')


