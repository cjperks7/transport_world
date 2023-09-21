'''

_get_cmod.py is a module to faciliate obtaining 
CMOD data useful for transport studies
    1) P_rad
    2) DD neutron rate
    3) Zeff
    4) Prf_in
    5) P_ohm
    6) Ar gas injection

TBD:
    1) include break-in-slope?
    2) spectroscopy?


'''

import numpy as np
import MDSplus
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})

__all__ = [
    'get_cmod',
    'plt_cmod',
    'rscl_cmod',
]


##########################################################
#
#                     Main
#
##########################################################

def get_cmod(
    dout = None,
    shot = None,
    quants = None,
    ):

    # If empty, load in everything
    if quants is None:
        quants = [
            'rad',
            'ohm',
            'rf',
            'fus',
            'Zeff',
            'Ar_gas',
            'Te/Ti',
            'H/D'
        ]

    # Initializes output
    dout['exp'] = {}
    dout['device'] = 'CMOD'
    dout['shot'] = int(shot)

    # Loads Bolometry data
    if 'rad' in quants:
        dout = _get_Prad(
            dout = dout,
            )

    # Loads DD neutron rate data
    if 'fus' in quants:
        dout = _get_DD(
            dout=dout,
            )

    # Loads injected RF power
    if 'rf' in quants:
        dout = _get_rf(
            dout=dout,
            )

    # Loads Z effective
    if 'Zeff' in quants:
        dout = _get_zeff(
            dout=dout,
            )

    # Loads ohmic power
    if 'ohm' in quants:
        dout = _get_ohmic(
            dout=dout,
            )

    # Loads Ar gas injection
    if 'Ar_gas' in quants:
        dout = _get_Ar_gas(
            dout=dout,
            )

    # Calculates Te/Ti
    if 'Te/Ti' in quants:
        dout = _get_TeoTi(
            dout=dout,
            )

    # Loads hydrogen-to-deuterium ratio
    if 'H/D' in quants:
        dout = _get_HoD(
            dout=dout,
            )

    # Loads stored energy
    if 'W_MHD' in quants:
        dout = _get_WMHD(
            dout=dout,
            )

    return dout


##########################################################
#
#                     Plotting
#
##########################################################

def plt_cmod(
    dout=None,
    ):

    # Number of subplots
    nsubs = len(dout['exp'].keys())

    # Initializes figure
    fig, ax  = plt.subplots(2, int(np.ceil(nsubs/2)))
    fig.tight_layout(pad=1.0)

    fig.suptitle(dout['shot'])

    # Loop over quanity
    for ii, quant in enumerate(dout['exp'].keys()):
        # Special case if Te/Ti plot
        if quant == 'Te/Ti':
            ax[int(ii%2), int(np.floor(ii/2))].plot(
                dout['rhot'],
                dout['Te_keV'],
                label = 'Te'
                )
            ax[int(ii%2), int(np.floor(ii/2))].plot(
                dout['rhot'],
                dout['Ti_keV'],
                label = 'Ti'
                )

            ax[int(ii%2), int(np.floor(ii/2))].set_xlim(
                [0.7,1]
                )
            ax[int(ii%2), int(np.floor(ii/2))].set_ylim(
                [0,1]
                )

            ax[int(ii%2), int(np.floor(ii/2))].set_xlabel(
                r'$\rho_t$'
                )
            ax[int(ii%2), int(np.floor(ii/2))].set_ylabel(
                r'$T_s$ [$keV$]'
                )
            leg = ax[int(ii%2), int(np.floor(ii/2))].legend()
            leg.set_draggable('on')
            ax[int(ii%2), int(np.floor(ii/2))].grid('on')

            continue

        # Loop over diags
        for jj in np.arange(len(dout['exp'][quant]['diags'])):
            ax[int(ii%2), int(np.floor(ii/2))].plot(
                dout['exp'][quant]['time'][jj],
                dout['exp'][quant]['val'][jj],
                label = dout['exp'][quant]['diags'][jj]
                )

        # Plots calculated value
        if quant in dout['powers'].keys():
            ax[int(ii%2), int(np.floor(ii/2))].plot(
                dout['t0_s'],
                dout['powers'][quant]['tot_MW']*1e3, # [MW] - > [kW]
                'r*',
                label = 'sim.',
                )
        elif quant == 'Zeff':
            ax[int(ii%2), int(np.floor(ii/2))].plot(
                dout['t0_s'],
                dout[quant]['avg'],
                'r*',
                label = 'sim.',
                )

        elif quant == 'H/(H+D)':
            ax[int(ii%2), int(np.floor(ii/2))].plot(
                dout['t0_s'],
                dout['exp'][quant]['sim'],
                'r*',
                label = 'sim.',
                )

        leg = ax[int(ii%2), int(np.floor(ii/2))].legend()
        leg.set_draggable('on')
        ax[int(ii%2), int(np.floor(ii/2))].grid('on')

        ax[int(ii%2), int(np.floor(ii/2))].set_xlim(0.5,1.5)
        ax[int(ii%2), int(np.floor(ii/2))].set_xlabel('time [s]')

        ax[int(ii%2), int(np.floor(ii/2))].set_ylabel(
            quant + ' [' + dout['exp'][quant]['units'] + ']'
            )

##########################################################
#
#              Rescale to experiment
#
##########################################################

def rscl_cmod(
    dout=None,
    dquants=None, # 'ohm', 'rad', 'rf', 'fus', 'Zeff'
    replot = None,
    ):

    # Loop over quantities to rescale
    for qnt in dquants.keys():
        if qnt == 'Zeff':
            # Rescales
            dout['Zeff']['prof'] *= dquants[qnt]
            dout['Zeff']['avg'] *= dquants[qnt]

            # Documents
            if 'rescale' not in dout['Zeff'].keys():
                dout['Zeff']['rescale'] = 1
            dout['Zeff']['rescale'] *= dquants[qnt]

        elif qnt == 'Te/Ti':
            # Rescales
            dout['Ti_keV'] += dquants[qnt]/1e3

            # Documents
            if 'Ti_rescale_eV' not in dout.keys():
                dout['Ti_rescale_eV'] = 0
            dout['Ti_rescale_eV'] += dquants[qnt]

        else:
            dout['powers'][qnt]['tot_MW'] *=dquants[qnt]
            dout['powers'][qnt]['prof_MW/m3'] *= dquants[qnt]

            # Documents
            if 'rescale' not in dout['powers'][qnt].keys():
                dout['powers'][qnt]['rescale'] = 1
            dout['powers'][qnt]['rescale'] *= dquants[qnt]

            # Adjusts other related components
            if qnt == 'fus':
                for qqq in ['alpi', 'alpe']:
                    dout['powers'][qqq]['prof_MW/m3'] *= dquants[qnt]

                    # Documents
                    if 'rescale' not in dout['powers'][qqq].keys():
                        dout['powers'][qqq]['rescale'] = 1
                    dout['powers'][qqq]['rescale'] *= dquants[qnt]

            elif qnt == 'rad':
                for qqq in ['sync', 'cont', 'line']:
                    dout['powers'][qqq]['tot_MW'] *=dquants[qnt]
                    dout['powers'][qqq]['prof_MW/m3'] *= dquants[qnt]

                    # Documents
                    if 'rescale' not in dout['powers'][qqq].keys():
                        dout['powers'][qqq]['rescale'] = 1
                    dout['powers'][qqq]['rescale'] *= dquants[qnt]

            elif qnt == 'rf':
                for qqq in ['rfe', 'rfi']:
                    dout['powers'][qqq]['tot_MW'] *=dquants[qnt]
                    dout['powers'][qqq]['prof_MW/m3'] *= dquants[qnt]

                    # Documents
                    if 'rescale' not in dout['powers'][qqq].keys():
                        dout['powers'][qqq]['rescale'] = 1
                    dout['powers'][qqq]['rescale'] *= dquants[qnt]

    # If one wishes to replot
    if replot:
        plt_cmod(dout=dout)

    return dout

##########################################################
#
#             Specialty functions
#
##########################################################

# Loads radiated power from Bolometry
def _get_Prad(
    dout = None,
    ):

    # Empirical rescale value
    fact = 3 # Jerry say 3, Bob says 4.5

    # MDSplus tree
    spec = MDSplus.Tree('spectroscopy', dout['shot'])

    # Radiated power from twopi diode
    bolo_nd_dio = spec.getNode(r'\spectroscopy::top.bolometer:twopi_diode') 
    bolo_dio = bolo_nd_dio.data()*fact # [kW]
    t_bolo_dio = bolo_nd_dio.dim_of(0).data() # [s]

    # Radiated power from twopi foil
    bolo_nd_foi = spec.getNode(r'\spectroscopy::top.bolometer:twopi_foil') 
    bolo_foil = bolo_nd_foi.data()/1e3 # [kW]
    t_bolo_foil = bolo_nd_foi.dim_of(0).data() # [s]

    # Stores output
    dout['exp']['rad'] = {}
    dout['exp']['rad']['val'] = [bolo_dio, bolo_foil]
    dout['exp']['rad']['time'] = [t_bolo_dio, t_bolo_foil]
    dout['exp']['rad']['units'] = 'kW'
    dout['exp']['rad']['diags'] = ['twopi_diode', 'twopi_foil']

    return dout


# Loads DD neutron rate
def _get_DD(
    dout = None,
    ):

    # MDSplus tree
    part = MDSplus.Tree('particles', dout['shot'])

    # Nodes for different diags
    neu_glo_nd = part.getNode(r'\particles::top.neutrons.global.results:neut_rate')
    neu_he3_nd = part.getNode(r'\particles::top.neutrons.he_3_bank.results:he3_nrate')

    # Loads in data
    neu_he3 = neu_he3_nd.data() # [1/s]
    t_neu_he3 = neu_he3_nd.dim_of(0).data() # [s]
    neu_glo = neu_glo_nd.data() # [1/s]
    t_neu_glo = neu_glo_nd.dim_of(0).data() # [s]

    neu_he3 = neu_he3[t_neu_he3 > 0.4]
    t_neu_he3 = t_neu_he3[t_neu_he3 > 0.4]
    neu_glo = neu_glo[t_neu_glo > 0.4]
    t_neu_glo = t_neu_glo[t_neu_glo > 0.4]

    # Conversion from rate to power, [kJ]
    E_DD = (0.5*(1.01+3.02) + 0.5*(0.82+2.45)) *1e3*1.6021e-19

    neu_he3 *= 2*E_DD
    neu_glo *= 2*E_DD

    # Store DD fusion rate data
    dout['exp']['fus'] = {}
    dout['exp']['fus']['val'] = [neu_he3, neu_glo]
    dout['exp']['fus']['time'] = [t_neu_he3, t_neu_glo]
    dout['exp']['fus']['units'] = 'kW'
    dout['exp']['fus']['diags'] = ['he3', 'global']

    return dout

# Loads Ohmic power
def _get_ohmic(
    dout = None,
    ):

    # MDSplus tree
    ana = MDSplus.Tree('analysis', dout['shot'])

    # Boundary flux
    ssibry_nd = ana.getNode(r'\analysis::efit_ssibry')
    ssibry= ssibry_nd.data()
    t_ssibry = ssibry_nd.dim_of(0).data()

    # Surface voltage
    vsurf = 2*np.pi * np.gradient(ssibry, t_ssibry) # [V]

    # Plasma current
    ip = abs(ana.getNode(r'\analysis::efit_aeqdsk:cpasma').data()) # [A]

    # Internal inductance
    li = ana.getNode(r'\analysis::efit_aeqdsk:ali').data() # internal inductance

    # Inductance
    L = li*6.28*67*1e-9 # [H]

    # Inductive Voltage
    vi = L * np.gradient(ip, t_ssibry) # [V], inductive voltage

    # Stores ohmic power
    dout['exp']['ohm'] = {}
    dout['exp']['ohm']['val'] = [ip*(vsurf-vi)/1e3] # [kW], dim(t,)
    dout['exp']['ohm']['time'] = [t_ssibry] # [s]
    dout['exp']['ohm']['units'] = 'kW'
    dout['exp']['ohm']['diags'] = ['ssibry']

    return dout

# Loads injected RF power
def _get_rf(
    dout = None,
    ):

    # MDSplus tree
    rf_nd = MDSplus.Tree('rf', dout['shot'])

    # Gets measured injected power
    rf = rf_nd.getNode(r'\rf::RF_POWER_NET').data()*1e3 # dim(t,); [kW]
    t_rf = np.array([ii/1e4 for ii in np.arange(rf.shape[0])]) # [s]

    # Loads in the LH tree for this shot
    try:
        lh_nd = MDSplus.Tree('lh', dout['shot'])
        lhnet = lh_nd.getNode(r'\LH::TOP.RESULTS:NETPOW') 
        lh = lhnet.data() # [kW]
        t_lh = lhnet.dim_of(0).data() # dim (t,); [s]
    except:
        lh = 0
        t_lh = 0

    # Stores injected RF power
    dout['exp']['rf'] = {}
    dout['exp']['rf']['val'] = [rf, lh]
    dout['exp']['rf']['time'] = [t_rf, t_lh]
    dout['exp']['rf']['units'] = 'kW'
    dout['exp']['rf']['diags'] = ['ICRH', 'LH']

    return dout

# Loads Z effective
def _get_zeff(
    dout = None,
    ):
    # Zeff from neoclassical resistivity calc
    from portals.experiment_tools.CMODtools import getZeff_neo

    # MDSplus tree
    spec = MDSplus.Tree('spectroscopy', dout['shot'])

    # Zeff from visible bremsstrahlung
    Zeff_nd = spec.getNode(r'\spectroscopy::z_ave')
    Zeff_VB = Zeff_nd.data()
    t_Zeff_VB = Zeff_nd.dim_of(0).data() # [s]

    # Calculates Zeff utilizing zeff_neo.pro
    Zeff_neo, t_Zeff_neo = getZeff_neo(dout['shot'])

    # Store Zeff
    dout['exp']['Zeff'] = {}
    dout['exp']['Zeff']['val'] = [Zeff_VB, Zeff_neo]
    dout['exp']['Zeff']['time'] = [t_Zeff_VB, t_Zeff_neo]
    dout['exp']['Zeff']['units'] = ''
    dout['exp']['Zeff']['diags'] = ['vis. brem.', 'neo. resis.']

    return dout


def _get_Ar_gas(
    dout = None,
    ):

    # MDSplus tree
    spec = MDSplus.Tree('spectroscopy', dout['shot'])

    # Loads Ar gas feed data
    gas_nd = spec.getNode(r'\spectroscopy::top.x_ray_pha:incaa16:gas_b_sd_lo')
    gas = gas_nd.data()
    t_gas = gas_nd.dim_of(0).data()

    # Stores data
    dout['exp']['Ar_gas'] = {}
    dout['exp']['Ar_gas']['val'] = [gas]
    dout['exp']['Ar_gas']['time'] = [t_gas]
    dout['exp']['Ar_gas']['units'] = ''
    dout['exp']['Ar_gas']['diags'] = ['incaa16']

    return dout

# Loads hydrogen-to-deuterium ratio
def _get_HoD(
    dout = None,
    ):

    # MDSplus tree
    spec = MDSplus.Tree('spectroscopy', dout['shot'])

    # Loads Ar gas feed data
    balmer_nd = spec.getNode(r'\spectroscopy::balmer_h_to_d')
    balmer = balmer_nd.data()
    t_balmer = balmer_nd.dim_of(0).data()

    # Stores data
    dout['exp']['H/(H+D)'] = {}
    dout['exp']['H/(H+D)']['val'] = [balmer/(1+balmer)]
    dout['exp']['H/(H+D)']['time'] = [t_balmer]
    dout['exp']['H/(H+D)']['units'] = ''
    dout['exp']['H/(H+D)']['diags'] = ['balmer']

    # Value used in simulation
    dout['exp']['H/(H+D)']['sim'] = np.trapz(
        (
            dout['ions']['H']['ni_tot_19m3']
            /(dout['ions']['H']['ni_tot_19m3']
            +dout['ions']['D']['ni_tot_19m3'])
            ),
        dout['vol_m3']
        )/dout['vol_m3'][-1]

    return dout

# Calculates Te/Ti
def _get_TeoTi(
    dout = None,
    ):

    # Calculates Te/Ti
    dout['exp']['Te/Ti'] = dout['Te_keV']/dout['Ti_keV'] # dim(rhop,)

    return dout

def _get_WMHD(
    dout = None,
    ):

    # MDSplus tree
    spec = MDSplus.Tree('analysis', dout['shot'])

    # Loads Ar gas feed data
    WMHD_nd = spec.getNode(r'\analysis::EFIT_AEQDSK:wplasm')/1e3
    WMHD = WMHD_nd.data()
    t_WMHD = WMHD_nd.dim_of(0).data()

    # Stores data
    dout['exp']['W_MHD'] = {}
    dout['exp']['W_MHD']['val'] = [WMHD]
    dout['exp']['W_MHD']['time'] = [t_WMHD]
    dout['exp']['W_MHD']['units'] = 'kJ'
    dout['exp']['W_MHD']['diags'] = ['EFIT']

    return dout