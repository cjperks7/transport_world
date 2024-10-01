'''

_get_fits.py is a module managing fitted kinetic profiles

cjperks

06/01/23

'''

import numpy as np
import sys, os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from transport_world.make_gacode import read_tokamak as rTok
sys.path.insert(0,'/home/cjperks/usr/python3modules/eqtools3')
import eqtools
sys.path.pop(0)

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

    # Init
    profs = ['ne_19m3', 'Te_keV', 'Ti_keV'] 
    ddata = {}

    for prof in profs:
        if prof == 'ne_19m3':
            prof_key = 'n_e'
            scale = 10
        elif prof == 'Te_keV':
            prof_key = 'T_e'
            scale = 1
        elif prof == 'Ti_keV':
            prof_key = 'T_i'
            scale = 1

        # Assumes files are in the form
        #   path_input + '/ne' + path_kin
        #   path_input + '/Te' + path_kin
        ff = netCDF4.Dataset(
            os.path.join(
                path_input,
                prof.split('_')[0].lower() + path_kin
                ),
            'r'
            )

        # Error checks
        if prof_key not in ff.variables.keys():
            others = [', CTS', ', HIREXSR, He-like Ar']
            for oth in others:
                if prof_key + oth in ff.variables.keys():
                    prof_key += oth

        # Radial coordinate
        if 'r/a' in ff.variables.keys():
            rad_val = dout['rmin_m']/dout['rmin_m'][-1]
            rad_key = 'r/a'
        elif 'sqrt{psi_n}' in ff.variables.keys():
            rad_val = dout['rhop']
            rad_key = 'sqrt{psi_n}'


        # Interpolates onto grid
        dout[prof] = interp1d(
            np.asarray(ff.variables[rad_key][:]),
            np.asarray(ff.variables[prof_key])*scale,
            bounds_error=False,
            fill_value=(float(ff.variables[prof_key][0])*scale, 0)
            )(rad_val)

        dout[prof+'_err'] = interp1d(
            np.asarray(ff.variables[rad_key][:]),
            np.asarray(ff.variables['err_'+prof_key])*scale,
            bounds_error=False,
            fill_value=(float(ff.variables['err_'+prof_key][0])*scale, 0)
            )(rad_val)


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
    edr = None,
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
        psin_tmp = np.zeros(diag['r_m'].shape)

        for ii in np.arange(diag['r_m'].shape[1]):
            psin_tmp[:,ii] = np.sqrt(edr.rmid2psinorm(
                list(diag['r_m'][:,ii]),
                diag['t_s'][ii],
                sqrt = False,
                ))**2

    # psin -> output x-variable
    elif 'psin' in diag.keys():
        psin_tmp = diag['psin']

    # Interpolates
    xnew_all = interp1d(
        dout['rhop']**2,
        xvar,
        bounds_error=False,
        fill_value = (0,1)
        )(psin_tmp)

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
    plt_q = True,
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

    # Plots profile fits
    if ne is not None:
        ax[0].plot(xvar, ne)
    if Te is not None:
        ax[1].plot(xvar, Te)
    if Ti is not None:
        ax[2].plot(xvar, Ti)

    keys = ['ne_19m3', 'Te_keV', 'Ti_keV']
    colors = ['r', 'b', 'g']
    for ii, key in enumerate(keys):
        ax[ii].plot(
            xvar,
            dout[key],
            label = 'fit tme-avg',
            color = colors[ii],
            )
        ax[ii].fill_between(
            xvar,
            dout[key]+dout[key+'_err'],
            dout[key]-dout[key+'_err'],
            alpha = 0.6,
            color = colors[ii],
            )
        ax[ii].fill_between(
            xvar,
            dout[key]+2*dout[key+'_err'],
            dout[key]-2*dout[key+'_err'],
            alpha = 0.3,
            color = colors[ii],
            )

    # -----------------
    # Experimental data
    # -----------------

    if shot is not None:
        # Obtains experimental kinetic profiles
        dkin = rTok.profs_cmod(shot=shot)

        # Gets eqdsk
        edr = eqtools.EqdskReader(
            gfile=os.path.join(
                dout['paths']['input'],
                dout['paths']['gfile']
                ),
            afile=None
            )

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
                elif np.isnan(dkin[prof][diag][vals]).all():
                    continue
                elif diag in ['GPC', 'GPC2', 'FRC']:
                    continue

                # Finds time range
                t_ind = np.where(
                    (dkin[prof][diag]['t_s'] >=dt[0]) 
                    & (dkin[prof][diag]['t_s'] <=dt[1]) 
                    )[0]

                if diag == 'TCI':
                    nebar = np.mean(dkin[prof][diag][vals][t_ind])/1e19

                    ax[num].plot([0,1], [nebar,nebar], 'g--', label=diag)

                else:
                    # Uncertainty weighted average
                    if err in dkin[prof][diag].keys():
                        weights = 1/dkin[prof][diag][err][:,t_ind]**2
                        weights[np.isinf(weights)] = 0
                        weights[np.isnan(weights)] = 0

                        val = (
                            (dkin[prof][diag][vals][:,t_ind] *weights).sum(axis=1)
                            /weights.sum(axis=1)
                            )
                        val_std = np.sqrt(1/weights.sum(axis=1))
                    else:
                        weights = None
                        val = np.mean(dkin[prof][diag][vals][:,t_ind], axis=1)
                        val_std = val * 0.1

                    if prof == 'ne':
                        val /= 1e19
                        val_std /= 1e19

                    xnew_all = conv_rad(
                        dout = dout,
                        diag = dkin[prof][diag],
                        xout = plt_coor,
                        edr = edr,
                        )
                    
                    if xnew_all.ndim > 1:
                        if weights is None:
                            xnew = np.mean(xnew_all[:,t_ind],axis=1)
                        else:
                            xnew = (
                                (xnew_all[:,t_ind]*weights).sum(axis=1)/
                                weights.sum(axis=1)
                                )
                    else:
                        xnew = xnew_all
                    
                    ax[num].errorbar(
                        xnew, 
                        val, 
                        yerr = val_std,
                        fmt='*', label=diag
                        )

    if plt_q:
        ind = np.argmin(abs(dout['q'] - 1))

        ax[0].plot(
            [xvar[ind], xvar[ind]],
            [0, 20],
            color = 'm',
            label = 'q =1'
            )
        ax[1].plot(
            [xvar[ind], xvar[ind]],
            [0, 3],
            color = 'm',
            label = 'q =1'
            )
        ax[2].plot(
            [xvar[ind], xvar[ind]],
            [0, 3],
            color = 'm',
            label = 'q =1'
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

    

    fig.suptitle('Fitted Profiles; %i; t=[%0.2f, %0.2f] s'%(shot, dt[0], dt[1]))

    plt.show()