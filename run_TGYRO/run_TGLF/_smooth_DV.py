'''

_smooth_DV.py is a function is facilitate
smoothing TGLF output

'''

# Modules
import numpy as np
from . import _smearr as _sm

__all__ = [
    'smooth_DV',
    ]


##########################################################
#
#                   Main
#
##########################################################

def smooth_DV(
    # Dictionary of D,V data
    ddata=None,
    # Where to chop data
    rho_accept = None,
    # Algorithm controls
    alfa = 8,
    alg = 'simple',
    # Plotting
    plt_all = False,
    ):

    # Smooths TGLF output
    ind_accept = np.argmin(abs(ddata['rhot']-rho_accept))

    # Runs smoothing function
    VZ_smooth = _sm.smearr(
        alfa = alfa,
        f_in = ddata['VZ'][ind_accept:],
        rho = ddata['rhot'][ind_accept:],
        alg = 'simple',
        )

    DZ_smooth = _sm.smearr(
        alfa = alfa,
        f_in = ddata['DZ'][ind_accept:],
        rho = ddata['rhot'][ind_accept:],
        alg = 'simple',
        )

    # Output dictionary
    dout = {
        'rhot': np.array(ddata['rhot'][ind_accept:]),
        'DZ': np.array(DZ_smooth),
        'VZ': np.array(VZ_smooth),
        }

    # Plotting
    if plt_all:
        _plot(
            ddata=ddata,
            dout=dout,
            )

    # Output
    return dout


##########################################################
#
#                   Main
#
##########################################################

def _plot(ddata=None,dout=None):
    R_CMOD = 0.68 # [m]

    fig, ax = plt.subplots(3,1)

    ax[0].plot(
        ddata['rhot'],
        ddata['DZ'],
        '*')
    ax[1].plot(
        ddata['rhot'],
        ddata['VZ'],
        '*'
        )

    ax[0].plot(
        dout['rhot'],
        dout['DZ'],
        '*'
        )

    ax[1].plot(
        dout['rhot'],
        dout['VZ'],
        '*'
        )
    ax[2].plot(
        dout['rhot'],
        R_CMOD*dout['VZ']/dout['DZ'],
        '*',
        )


    ax[0].set_xlim(0,1)
    ax[1].set_xlim(0,1)
    ax[2].set_xlim(0,1)

    ax[2].set_ylim(-10,5)

    ax[2].set_xlabel('rhot')

    ax[0].grid('on')
    ax[1].grid('on')
    ax[2].grid('on')

    ax[0].set_ylabel('D [m^2/s]')
    ax[1].set_ylabel('V [m/s]')
    ax[2].set_ylabel('RV/D')

    fig.show()