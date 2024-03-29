'''

_get_fits.py is a module managing fitted kinetic profiles

cjperks

06/01/23

'''

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from transport_world.make_gacode import read_tokamak as rTok

plt.rcParams.update({'font.size': 16})

__all__ = [
    'get_fits',
    'get_imprad',
    'conv_rad',
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
    plt_coor='rhot',
    source = 'quickfit',
    ):
    dout['paths']['kin'] = path_kin
    dout['t0_s'] = t0


    if source == 'quickfit':
        dout = _read_quickfit(
            dout=dout,
            path_input=path_input,
            path_kin=path_kin,
            t0=t0,
            dt=dt,
            )

    elif source == 'profiletools':
        dout = _read_profiletools(
            dout=dout,
            path_input=path_input,
            path_kin=path_kin,
            )

    if plt_all:
        _plot(
            dout = dout,
            shot = shot,
            ne = None,
            Te = None,
            Ti = None,
            t0=t0,
            dt=dt,
            plt_coor=plt_coor,
            )

    return dout


def get_imprad(
    fimprad = None,
    ind_key = 0,
    rhop = None,
    # Rescaling controls
    rescl = 1,
    con = None,
    ne_cm3 = None,
    vol = None,
    ):

    with open(fimprad, 'rb') as f: 
        ddata = np.load(f, allow_pickle=True)[()] # data dictionary

    if rhop is None:
        return ddata

    # Returns impurity density profile
    else:
        key = list(ddata.keys())[ind_key]

        nz_cm3 = interp1d(
            ddata[key]['result']['rhop'],
            ddata[key]['result']['n_imp'],
            axis = 0
            )(rhop) # dim(rhop, cs), [cm^-3]

        # If user wants to rescale to a concentration
        if con is not None:
            rescl = con/ (np.trapz(
                nz_cm3.sum(1)/ne_cm3,
                vol
                )/vol[-1])
            print(rescl)
            
            return ddata, nz_cm3*rescl, rescl

        # Else if user wants to rescale to known value
        else:
            return ddata, nz_cm3*rescl


##########################################################
#
#                     Utilities
#
##########################################################

# Reads fit data prepared by QuickFit within OMFIT
def _read_quickfit(
    dout=None,
    path_input=None,
    path_kin=None,
    t0=None,
    dt=None,
    ):

    # Modules
    from uncertainties import unumpy

    # Opens .npy file of saved profiles from QuickFit
    with open(path_input+path_kin, 'rb') as f: 
        rho_t = np.load(f) # norm. sq. tor. flux
        time = np.load(f)/1e3 # [s]
        ne_t = np.load(f, allow_pickle=True)/1e19 # dim(t, rho); [1e19 m^-3]
        Te_t = np.load(f, allow_pickle=True)/1e3 # dim(t, rho); [keV]
        try:
            Ti_t = np.load(f, allow_pickle=True)/1e3 # dim(t, rho); [keV]
            copied = False
        except:
            Ti_t = Te_t.copy()
            copied = True

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
    #Te[Te<0.2] = 0.2
    #Ti[Ti<0.2] = 0.2
    #ne[ne<0.5] = 0.5

    # Interpolates onto rho = 0->1
    ne = interp1d(rho_t, ne, axis = 0)(dout['rhot'])
    Te = interp1d(rho_t, Te, axis = 0)(dout['rhot'])

    if copied:
        Ti = interp1d(rho_t, Ti, axis = 0)(dout['rhot'])
    else:
        ### NOTE: I found an error under /home/sciortino/quickfit/CMOD/fetch_data.py
        # Assumes HIREX's rho = rho_t, but it's actually psin
        #Ti = interp1d(rho_t, Ti, axis = 0)(dout['rhot'])
        Ti = interp1d(rho_t, Ti, axis=0)(dout['rhop']**2)

    # Takes the average in time
    dout['ne_19m3'] = np.mean(ne, axis=-1)
    dout['Te_keV'] = np.mean(Te, axis=-1)
    dout['Ti_keV'] = np.mean(Ti, axis=-1)

    # Output
    return dout

# Reads fit data prepared by profiletools
def _read_profiletools(
    dout=None,
    path_input=None,
    path_kin=None,
    ):

    # Modules
    import netCDF4    

    # Assumes files are in the form
    #   path_input + '/ne' + path_kin
    #   path_input + '/Te' + path_kin
    ne = netCDF4.Dataset(path_input+'/ne'+path_kin, 'r')
    Te = netCDF4.Dataset(path_input+'/Te'+path_kin, 'r')

    # Radial coordinate
    if 'r/a' in Te.variables.keys():
        rad_Te = dout['rmin_m']/dout['rmin_m'][-1]
        key_Te = 'r/a'

    if 'r/a' in ne.variables.keys():
        rad_ne = dout['rmin_m']/dout['rmin_m'][-1]
        key_ne = 'r/a'

    # Interpolates onto grid
    dout['ne_19m3'] = interp1d(
        np.asarray(ne.variables[key_ne][:]),
        np.asarray(ne.variables['n_e'])*10,
        bounds_error=False,
        fill_value=(float(ne.variables['n_e'][0])*10, 0)
        )(rad_ne)
    dout['Te_keV'] = interp1d(
        np.asarray(Te.variables[key_Te][:]),
        np.asarray(Te.variables['T_e']),
        bounds_error=False,
        fill_value=(float(Te.variables['T_e'][0]), 0)
        )(rad_Te)
    
    # profiletools doesn't do Ti
    dout['Ti_keV'] = dout['Te_keV'].copy()

    # Output
    return dout

##########################################################
#
#                     Extra
#
##########################################################

def conv_rad(
    dout=None,
    diag=None,
    xout='rhot',
    ):
    '''
    Converts various radial bases used for diags to sq. norm. tor. flux
    '''

    if xout == 'rhot':
        xvar = dout['rhot']
    elif xout == 'rhop':
        xvar = dout['rhop']
    elif xout == 'rmin_m':
        xvar = dout['rmin_m']

    # rmin -> output x-variable
    if 'r_m' in diag.keys():
        xnew_all = interp1d(
            dout['rmin_m'],
            xvar,
            bounds_error=False,
            fill_value = (0, dout['rmin_m'][-1])
            )(diag['r_m']-float(dout['rcentr_m']))

    # psin -> output x-variable
    elif 'psin' in diag.keys():
        xnew_all = interp1d(
            dout['rhop']**2,
            xvar,
            bounds_error=False,
            fill_value = (0,1)
            )(diag['psin'])

    return xnew_all

def _plot(
    dout=None,
    shot=None,
    ne=None,
    Te=None,
    Ti=None,
    t0=None,
    dt=None,
    plt_coor = 'rhot',
    ):

    # x-axis variable
    if plt_coor == 'rhot':
        xvar = dout['rhot']
        xlab = r'$\rho_t$'
    elif plt_coor == 'rhop':
        xvar = dout['rhop']
        xlab = r'$\rho_p$'
    elif plt_coor == 'rmin_m':
        xvar = dout['rmin_m']
        xlab = r'$r$ [$cm$]'

    fig, ax = plt.subplots(1,3)

    # Plots ne
    if ne is not None:
        ax[0].plot(xvar, ne)
    ax[0].plot(xvar, dout['ne_19m3'], 'k', label='fit time-avg')

    # Plots Te
    if Te is not None:
        ax[1].plot(xvar, Te)
    ax[1].plot(xvar, dout['Te_keV'], 'k', label='fit time-avg')

    # Plots Ti
    if Ti is not None:
        ax[2].plot(xvar, Ti)
    ax[2].plot(xvar, dout['Ti_keV'], 'k', label='fit time-avg')
  
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
                err = 'err_1e20m3'
                num = 0
            else:
                vals = 'val_keV'
                err = 'err_keV'
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
                #t_ind = np.where(
                #    (dkin[prof][diag]['t_s'] >=t0-dt) 
                #    & (dkin[prof][diag]['t_s'] <=t0+dt) 
                #    )[0]
                t_ind = np.argmin(
                    abs(
                        dkin[prof][diag]['t_s'] - t0
                        )
                    )

                if diag == 'TCI':
                    nebar = dkin[prof][diag][vals][t_ind]/1e19

                    ax[num].plot([0,1], [nebar,nebar], 'g--', label=diag)

                else:
                    val = dkin[prof][diag][vals][:,t_ind]

                    if err in dkin[prof][diag].keys():
                        val_std = dkin[prof][diag][err][:,t_ind]
                    else:
                        val_std = dkin[prof][diag][vals][:,t_ind] * 0.1

                    if prof == 'ne':
                        val /= 1e19
                        val_std /= 1e19

                    xnew_all = conv_rad(
                        dout = dout,
                        diag = dkin[prof][diag],
                        xout = plt_coor,
                        )
                    
                    if xnew_all.ndim > 1:
                        xnew = xnew_all[:,t_ind]
                    else:
                        xnew = xnew_all
                    
                    ax[num].errorbar(
                        xnew, 
                        val, 
                        yerr = val_std,
                        fmt='*', label=diag
                        )


    ax[0].set_xlabel(xlab)
    ax[0].set_ylabel(r'$n_e$ [1e19 $m^{-3}$]')
    ax[0].grid('on')
    ax[0].set_ylim(0,20)
    leg0 = ax[0].legend()
    leg0.set_draggable('on')

    ax[1].set_xlabel(xlab)
    ax[1].set_ylabel(r'$T_e$ [$keV$]')
    ax[1].grid('on')
    ax[1].set_ylim(0,3)
    leg1 = ax[1].legend()
    leg1.set_draggable('on')

    ax[2].set_xlabel(xlab)
    ax[2].set_ylabel(r'$T_i$ [$keV$]')
    ax[2].grid('on')
    ax[2].set_ylim(0,3)
    leg2 = ax[2].legend()
    leg2.set_draggable('on')

    fig.suptitle('Fitted Profiles t='+str(t0-dt)+'-'+str(t0+dt)+'s')

    plt.show()