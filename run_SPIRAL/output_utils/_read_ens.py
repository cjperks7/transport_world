'''

Functions to read and process ensemble output data
from SPIRAL packaged by /u/cperks2/SPIRAL/python/read_ens.py

cjperks
Dec 4th, 2024

'''

# Modules
import numpy as np
import sys, os

__all__ = [
    'read_ens'
    ]

###############################################
#
#           Main
#
###############################################

# Reads process .txt file
def read_ens(
    nmark = None,
    ens = None
    ):

    # Init
    ddata = {}

    # Opens file
    with open(ens, 'r') as f:
        # Header
        line = f.readline()

        # Header data
        keys = []
        for kk in line.split():
            ddata[kk] = np.zeros(int(nmark))
            keys.append(kk)

        # Counter
        cnt = 0

        # Loop over markers
        for line in f.readlines():
            for ii, ll in enumerate(line.split()):
                ddata[keys[ii]][cnt] = float(ll)
            cnt += 1

    # Misc
    ddata['R_START'] = np.sqrt(
        ddata['X_START']**2
        + ddata['Y_START']**2
        )
    ddata['R_END'] = np.sqrt(
        ddata['X_END']**2
        + ddata['Y_END']**2
        )
        
    # Output
    return ddata