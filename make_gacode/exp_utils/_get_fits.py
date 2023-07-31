'''

_get_fits.py is a module managing fitted kinetic profiles

cjperks

06/01/23

'''

import numpy as np
from uncertainties import unumpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from transport_world.make_gacode import read_tokamak as rTok

__all__ = [
    'get_fits',
]

##########################################################
#
#                     Main
#
##########################################################

def get_fits(
    dout = None,
    path_input = None,
    path_kin = None,
    shot = None,
    t0 = None,
    dt = None,
    plt_all = None,
    ):
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
    ne = interp1d(rho_t, ne, axis = 0)(dout['rhot'])
    Te = interp1d(rho_t, Te, axis = 0)(dout['rhot'])
    Ti = interp1d(rho_t, Ti, axis = 0)(dout['rhot'])

    # Takes the average in time
    dout['ne_19m3'] = np.mean(ne, axis=-1)
    dout['Te_keV'] = np.mean(Te, axis=-1)
    dout['Ti_keV'] = np.mean(Ti, axis=-1)

    if plt_all:
        _plot(
            dout = dout,
            shot = shot,
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

def _conv_rad(
    dout=None,
    diag=None,
    ):
    '''
    Converts various radial bases used for diags to sq. norm. tor. flux
    '''

    # rmin -> rhot
    if 'r_m' in diag.keys():
        rhot_all = interp1d(
            dout['rmin_m'],
            dout['rhot'],
            bounds_error=False,
            fill_value = (0, dout['rmin_m'][-1])
            )(diag['r_m']-float(dout['rcentr_m']))

    # psin -> rhot
    elif 'psin' in diag.keys():
        rhot_all = interp1d(
            dout['rhop']**2,
            dout['rhot'],
            bounds_error=False,
            fill_value = (0,1)
            )(diag['psin'])

    return rhot_all

def _plot(
    dout=None,
    shot=None,
    ne=None,
    Te=None,
    Ti=None,
    t0=None,
    dt=None,
    ):

    fig, ax = plt.subplots(1,3)

    # Plots ne
    if ne is not None:
        ax[0].plot(dout['rhot'], ne)
    ax[0].plot(dout['rhot'], dout['ne_19m3'], 'k', label='fit time-avg')

    # Plots Te
    if Te is not None:
        ax[1].plot(dout['rhot'], Te)
    ax[1].plot(dout['rhot'], dout['Te_keV'], 'k', label='fit time-avg')

    # Plots Ti
    if Ti is not None:
        ax[2].plot(dout['rhot'], Ti)
    ax[2].plot(dout['rhot'], dout['Ti_keV'], 'k', label='fit time-avg')
  
    # -----------------
    # Experimental data
    # -----------------

    if shot is not None:
        # Obtains experimental kinetic profiles
        dkin = rTok.profs_cmod(shot=shot)

        # Loop over kinetric profiles
        for prof in dkin.keys():
            if prof == 'ne':
                vals = 'val_1e20m3'
                num = 0
            else:
                vals = 'val_keV'
                if prof == 'Te':
                    num = 1
                elif prof == 'Ti':
                    num = 2

            # Plots electron density data
            for diag in dkin[prof].keys():
                # Skip unavailable diagnostics
                if len(dkin[prof][diag]) == 0:
                    continue

                # Finds time range
                t_ind = np.where(
                    (dkin[prof][diag]['t_s'] >=t0-dt) 
                    & (dkin[prof][diag]['t_s'] <=t0+dt) 
                    )[0]

                if diag == 'TCI':
                    nebar = np.mean(
                        dkin[prof][diag][vals][t_ind]
                        )/1e19
                    nebar_std = np.std(
                        dkin[prof][diag][vals][t_ind]
                        )/1e19

                    ax[num].plot([0,1], [nebar,nebar], 'g--', label=diag)

                else:
                    val = np.mean(
                        dkin[prof][diag][vals][:,t_ind],
                        axis = 1
                        )
                    if 'err_keV' in dkin[prof][diag].keys():
                        val_std = np.mean(
                            dkin[prof][diag]['err_keV'][:,t_ind],
                            axis = 1
                            )
                    else:
                        val_std = np.std(
                            dkin[prof][diag][vals][:,t_ind],
                            axis = 1
                            )

                    if prof == 'ne':
                        val /= 1e19
                        val_std /= 1e19

                    rhot_all = _conv_rad(
                        dout = dout,
                        diag = dkin[prof][diag]
                        )

                    rhot = np.mean(
                        rhot_all[:,t_ind],
                        axis = 1
                        )
                    rhot_std = np.std(
                        rhot_all[:,t_ind],
                        axis = 1
                        )

                    ax[num].errorbar(
                        rhot, 
                        val, 
                        yerr = val_std,
                        xerr = rhot_std,
                        fmt='*', label=diag
                        )


    ax[0].set_xlabel(r'$\rho_t$')
    ax[0].set_ylabel(r'$n_e$ [1e19 $m^{-3}$]')
    ax[0].grid('on')
    ax[0].set_ylim(5,15)
    leg0 = ax[0].legend()
    leg0.set_draggable('on')

    ax[1].set_xlabel(r'$\rho_t$')
    ax[1].set_ylabel(r'$T_e$ [$keV$]')
    ax[1].grid('on')
    ax[1].set_ylim(0,3)
    leg1 = ax[1].legend()
    leg1.set_draggable('on')

    ax[2].set_xlabel(r'$\rho_t$')
    ax[2].set_ylabel(r'$T_i$ [$keV$]')
    ax[2].grid('on')
    ax[2].set_ylim(0,2)
    leg2 = ax[2].legend()
    leg2.set_draggable('on')

    fig.suptitle('Fitted Profiles t='+str(t0-dt)+'-'+str(t0+dt)+'s')

    plt.show()