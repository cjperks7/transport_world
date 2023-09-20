'''

_read_idl_fits.py is a function to faciliate obtaining profile
raw and fitted data from CMOD using the IDL widgets w_dens 
and w_temp

NOTE: Assumes you have a folder  /IDL_fits/ within /path_input/
that then hold folders /ne/, /Te/, /Ti/ that stores the IDL
savesets v. time outputted by the widgets

cjperks
Sept. 20, 2023

'''

# Modules
import numpy as np
from scipy.io.idl import readsav

__all__ = [
    'read_idl',
    ]

####################################################################
#
#                       Main
#
####################################################################

def read_idl(
    path_input=None,
    shot = None,
    t0 = None,
    dt = None,
    profs = ['ne', 'Te'],
    quant = 'fit',
    ):

    # Reads IDL saveset (specific time) as python recarray
    data = readsav('dens_1140221012_0.819.dat')

    # Shot
    data['data']['shot']

    # Fit arrays
    data['data']['qfit_ne'] # [1e20 m^-3]
    data['data']['qfit_r'] # [cm]

