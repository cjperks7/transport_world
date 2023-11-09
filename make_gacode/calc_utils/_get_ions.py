'''

_get_ions.py is a module for faciliating ion modeling

cjperks

06/01/23

'''

import numpy as np
import scipy.constants as cnt
import matplotlib.pyplot as plt
from transport_world.make_gacode import exp_utils as utz
import os

plt.rcParams.update({'font.size': 16})

__all__ = [
    'get_ions',
]

##########################################################
#
#                     Main
#
##########################################################

def get_ions(
    dout = None,
    dmodel = None,
    plt_all = None,
    fimprad = None,
    ):

    # Initializes entry
    dout['ions'] = {}
    dout['e'] = {}

    # Initializes scalars
    dout['ions']['name'] = ""
    dout['ions']['itype'] = ""
    dout['ions']['mass'] = ""
    dout['ions']['charge'] = ""

    # Store electron values
    dout['e']['mass'] = " {:1.7E}".format(5.4488741e-4)
    dout['e']['charge']  = "{:1.7E}".format(-1)

    # Determines modeling
    if dmodel['option'] == 'concentration':
        dout = _calc_con(
            dout=dout,
            dmodel=dmodel,
            fimprad=fimprad,
        )

    else:
        print('NOT IMPLEMENTED YET')

    # Calculates the total pressure and Zeff
    dout['ptot_Pa'] = dout['ne_19m3'] * dout['Te_keV'] # dim(rhot,)
    dout['Zeff'] = {}
    dout['Zeff']['prof'] = 0
    for ion in list(filter(None,dout['ions']['name'].split(' '))):
        dout['ptot_Pa'] += dout['ions'][ion]['ni_tot_19m3'] * dout['Ti_keV'] # dim(rhot,)
        dout['Zeff']['prof'] += (
            dout['ions'][ion]['ni_tot_19m3'] /dout['ne_19m3'] 
            * dout['ions'][ion]['Z']**2 
            ) # dim(rhot,)
    dout['Zeff']['avg'] = np.trapz(dout['Zeff']['prof'], dout['vol_m3'])/dout['vol_m3'][-1]
    dout['ptot_Pa'] *= 1000*cnt.e*1e19 # [Pa]

    # Rotation modeling
    dout = _get_rot(
        dout=dout,
        dmodel=dmodel,
        )

    # Plotting
    if plt_all:
        _plot(dout=dout)

    return dout


##########################################################
#
#                     Transport modeling
#
##########################################################

# Uniform concentration model -> ni = ci * ne
def _calc_con(
    dout=None,
    dmodel=None,
    fimprad=None,
    ):

    # lumped impurity per QN
    nLumped = dout['ne_19m3'].copy()

    # counts number of ions
    nion = 0

    # Loop over ions
    for ion in dmodel['ions'].keys():
        dout['ions'][ion] = {}

        if (fimprad is not None) & (ion == 'Ar'):
            _, nz_cm3, _ = utz.get_imprad(
                fimprad=os.path.join(dout['paths']['input'],fimprad),
                rhop=dout['rhop'],
                con = dmodel['ions'][ion]['con'],
                ne_cm3 = dout['ne_19m3']*1e13,
                vol = dout['vol_m3'],
                )

            dout['ions'][ion]['nz_19m3'] = nz_cm3/1e13 # [1e19 m3], dim(rhop, cs)
            dout['ions'][ion]['ni_tot_19m3'] = dout['ions'][ion]['nz_19m3'].sum(1) # dim(rhop,)

        else:
            dout['ions'][ion]['ni_tot_19m3'] = dmodel['ions'][ion]['con'] * dout['ne_19m3']

        # Obtains charge, mass, type
        dout = _ion_scalars(
            dout=dout,
            ion=ion,
            )

        nLumped -= dout['ions'][ion]['ni_tot_19m3'] * dout['ions'][ion]['Z']

        # increases counter
        nion +=1

    # Stores lumped impurity
    if any(nLumped <0):
        print('TOO MUCH IMPURITY FOR QN!!')
    else:
        dout['ions']['LUMPED'] = {}
        dout = _ion_scalars(
            dout=dout,
            ion='LUMPED',
            )
        dout['ions']['LUMPED']['ni_tot_19m3'] = nLumped / dout['ions']['LUMPED']['Z']

        # increases counter
        nion +=1

    # Stores number of ion
    dout['ions']['nion'] = str(nion)
    dout['nexp'] = str(len(dout['rhot']))

    return dout


##########################################################
#
#           Rotation modeling
#
##########################################################

def _get_rot(
    dout=None,
    dmodel=None,
    ):

    # If user wishes to use static plasma
    if dmodel['rotation']['option'] == 'zero':
        dout['omega0_rad/s'] = 0*dout['rhot']
        dout['vtor_m/s'] = 0*dout['rhot']
        dout['vpol_m/s'] = 0*dout['rhot']

    return dout

##########################################################
#
#           Extra
#
##########################################################

def _ion_scalars(
    dout =None,
    ion = None,
    ):
    # masses
    if ion == 'D':
        m = 2.0147294
        z = 1
    elif ion == 'H':
        m = 1.00784
        z = 1
    elif ion == 'Ar':
        m = 39.948
        z = 16
    elif ion == 'Mo':
        m = 95.95
        z = 30
    elif ion == 'LUMPED' or ion == 'B':
        m = 10.811
        z = 5
    elif ion == 'He4':
        m = 4.002603
        z = 2

    # Stores values
    dout['ions'][ion]['Z'] = z
    dout['ions'][ion]['M'] = m

    dout['ions']['name'] += (str(ion)+" ")
    dout['ions']['itype'] += ('[therm] ')
    dout['ions']['mass'] += (" {:1.7E}".format(m))
    dout['ions']['charge'] += (" {:1.7E}".format(z))

    return dout

##########################################################
#
#           Plotting
#
##########################################################

def _plot(
    dout=None,
    ):

    fig, ax = plt.subplots()

    # Loop over ions
    for ion in list(filter(None,dout['ions']['name'].split(' '))):
        Zf = (
            dout['ions'][ion]['Z']
            * np.trapz(
                dout['ions'][ion]['ni_tot_19m3']/dout['ne_19m3'],
                dout['vol_m3'])/dout['vol_m3'][-1]
        )

        ax.plot(
            dout['rhot'], 
            dout['ions'][ion]['ni_tot_19m3']/dout['ions'][ion]['ni_tot_19m3'][0],
            label = (
                r'$Z_{' +ion+r'}f_{'+ion+r'}$='+'{:.2f}'.format(Zf)
                )
            )

    ax.set_xlabel(r'$\rho_t$')
    ax.set_ylabel(r'$n_i$ [norm]')
    ax.grid('on')
    ax.legend()

    plt.show()