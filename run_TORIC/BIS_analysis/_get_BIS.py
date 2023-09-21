'''

_get_BIS.py is a script meant to estimate the coupled RF power 
given the break-in-slope of the time derivative of the total 
stored energy

cjperks
Sept. 20, 2023

'''

# Modules
import numpy as np
from transport_world.make_gacode import read_tokamak as rTok
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})

__all__ = [
    'get_BIS'
    ]

#####################################################
#
#                   Main
#
#####################################################

def get_BIS(
    shot=None,
    t1 = [0.75, 1.0],
    t2 = [1.0, 1.5],
    plt_all = True,
    ):

    # Loads total power time traces
    dout = rTok.get_cmod(
        dout = {},
        shot=shot,
        quants = ['rad', 'ohm', 'rf', 'W_MHD'],
        )

    # Finds average powers in time range of interest
    kpwr = ['rad', 'ohm']
    ddata = {}

    for pwr in kpwr:
        ddata[pwr] = _time_avg(
            ddata=dout['exp'][pwr],
            t1=t1,
            t2=t2,
            )

    # Performs linear fit of W_MHD in time
    ddata['W_MHD'] = _line_fit(
        ddata=dout['exp']['W_MHD'],
        t1=t1,
        t2=t2
        )

    # Performs break-in-slope analysis
    ddata['BIS'] = _calc_BIS(
        ddata=ddata
        )

    # Plotting
    if plt_all:
        _plot(
            dout=dout,
            ddata=ddata,
            t1 = t1,
            t2 = t2,
            )

    # Output
    return ddata


#####################################################
#
#                   Utilities
#
#####################################################

# Performs time average over window
def _time_avg(
    ddata=None,
    t1=None,
    t2=None,
    ):

    # Initialize data
    val1 = []
    val2 = []

    # Loop over diags
    for jj in np.arange(len(ddata['diags'])):

        # Finds time windows of interest
        ind1 = np.where((ddata['time'][jj] >= t1[0]) & (ddata['time'][jj] <= t1[-1]))[0]
        ind2 = np.where((ddata['time'][jj] >= t2[0]) & (ddata['time'][jj] <= t2[-1]))[0]

        val1.append(np.mean(ddata['val'][jj][ind1]))
        val2.append(np.mean(ddata['val'][jj][ind2]))

    # Outputs time average of window
    return [np.mean(val1), np.mean(val2)]

# Performs linear fit over time window
def _line_fit(
    ddata=None,
    t1=None,
    t2=None,
    ):

    # Finds time windows of interest
    ind1 = np.where((ddata['time'][0] >= t1[0]) & (ddata['time'][0] <= t1[-1]))[0]
    ind2 = np.where((ddata['time'][0] >= t2[0]) & (ddata['time'][0] <= t2[-1]))[0]

    # Performs linear fit
    m1, b1 = np.polyfit(ddata['time'][0][ind1], ddata['val'][0][ind1], 1)
    m2, b2 = np.polyfit(ddata['time'][0][ind2], ddata['val'][0][ind2], 1)

    # Outputs linear fit
    return [
        [m1, b1, ddata['time'][0][ind1]],
        [m2, b2, ddata['time'][0][ind2]]
        ]

# Performs simple break-in-slope analysis
def _calc_BIS(
    ddata=None,
    ):

    # Calculates transported power when RF is off (assumed t1)
    Ptrans = ddata['ohm'][0] -ddata['rad'][0] - ddata['W_MHD'][0][0]

    # Calculates the coupled RF power
    # NOTE: Assumes constant transported power
    Prf = ddata['W_MHD'][1][0] + ddata['rad'][1] + Ptrans - ddata['ohm'][1]

    return [Ptrans, Prf]

#####################################################
#
#                   Plotting
#
#####################################################

def _plot(
    dout = None,
    ddata = None,
    t1 = None,
    t2 = None,
    ):

    # Plot of time traces
    fig,ax=plt.subplots(1,2)

    # Power time traces
    kpwr = ['rad', 'ohm', 'rf']

    for pwr in kpwr:
        # Loop over diags
        for jj in np.arange(len(dout['exp'][pwr]['diags'])):
            ax[0].plot(
                dout['exp'][pwr]['time'][jj],
                dout['exp'][pwr]['val'][jj],
                label = pwr +'_'+dout['exp'][pwr]['diags'][jj]
                )
            if pwr == 'rf':
                if dout['exp'][pwr]['diags'][jj] == 'ICRH':
                    ax[0].plot(
                        [t1[0], t2[-1]],
                        [ddata['BIS'][0], ddata['BIS'][0]],
                        'b--',
                        label='transport'
                        )

                    ax[0].plot(
                        t2,
                        [ddata['BIS'][1], ddata['BIS'][1]],
                        'k--',
                        label='RF coupled'
                        )

            else:
                ax[0].plot(
                    [np.mean(t1), np.mean(t2)],
                    ddata[pwr],
                    '*',
                    label = pwr + '_time_avg'
                    )

    ax[0].set_xlim(t1[0], t2[-1])
    ax[0].set_ylim(0,2000)

    ax[0].set_xlabel('time [s]')
    ax[0].set_ylabel('power [kW]')

    leg = ax[0].legend()
    leg.set_draggable('on')

    ax[0].grid('on')

    ax[0].plot(
        [t1[0], t1[0]],
        [0,2000],
        'r--'
        )
    ax[0].plot(
        [t1[-1], t1[-1]],
        [0,2000],
        'r--'
        )

    ax[0].plot(
        [t2[0], t2[0]],
        [0,2000],
        'r--'
        )
    ax[0].plot(
        [t2[-1], t2[-1]],
        [0,2000],
        'r--'
        )

    # Stored energy time traces
    ax[1].plot(
        dout['exp']['W_MHD']['time'][0],
        dout['exp']['W_MHD']['val'][0],
        label = 'W_MHD'
        )

    ax[1].plot(
        ddata['W_MHD'][0][-1],
        ddata['W_MHD'][0][0] * ddata['W_MHD'][0][-1] + ddata['W_MHD'][0][1],
        'k--'
        )
    ax[1].plot(
        ddata['W_MHD'][1][-1],
        ddata['W_MHD'][1][0] * ddata['W_MHD'][1][-1] + ddata['W_MHD'][1][1],
        'k--'
        )
    
    ax[1].set_xlim(t1[0], t2[-1])
    ax[1].set_ylim(3e1,5e1)

    ax[1].set_xlabel('time [s]')
    ax[1].set_ylabel('W_MHD [kJ]')

    ax[1].grid('on')

    ax[1].plot(
        [t1[0], t1[0]],
        [0,5e4],
        'r--'
        )
    ax[1].plot(
        [t1[-1], t1[-1]],
        [0,5e4],
        'r--'
        )

    ax[1].plot(
        [t2[0], t2[0]],
        [0,5e4],
        'r--'
        )
    ax[1].plot(
        [t2[-1], t2[-1]],
        [0,5e4],
        'r--'
        )
    

    