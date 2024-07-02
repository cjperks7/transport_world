'''

Script to read netCDF4 output files made by profiletools/gpfit

cjperks
June 28, 2024

'''

# Modules
import netCDF4
import numpy as np
import sys, os

__all__ = [
    'read',
    'write2ascii'
    ]

########################################
#
#           Main
#
########################################

# Reads netCDF output from GPfit
def read(
    file=None,
    ):
    # Init
    out = {}

    # Loads netCDF file
    nc = netCDF4.Dataset(file)

    # Loop over variables
    for var in nc.variables.keys():
        out[var] = {}

        # Get data
        out[var]['units'] = nc.variables[var].units
        out[var]['data'] = np.asarray(nc.variables[var][:])

    # Other data
    out['shot'] = nc.shot

    # Output
    return out

def write2ascii(
    fol=None,
    profs=None,
    ):

    # Opens file
    f.open(fol+'/fits_profiletools.txt', 'w')

    if (
        not np.array_equal(profs['ne']['sqrt{psi_n}']['data'], profs['Te']['sqrt{psi_n}']['data'])
        or not np.array_equal(profs['Ti']['sqrt{psi_n}']['data'], profs['Te']['sqrt{psi_n}']['data'])
        or not np.array_equal(profs['ne']['sqrt{psi_n}']['data'], profs['Ti']['sqrt{psi_n}']['data'])
        ):
        print('rho bases are not equal! Breaks logic')
        sys.exit(1)

    # Ignores values below cutoff
    cutoff = 1e-3
    profs['ne']['n_e, CTS']['data'][profs['ne']['n_e, CTS']['data']<cutoff] = cutoff
    profs['Te']['T_e, CTS']['data'][profs['Te']['T_e, CTS']['data']<cutoff] = cutoff
    profs['Ti']['T_i']['data'][profs['Ti']['T_i']['data']<cutoff] = cutoff

    # Writes lines
    for ii in range(len(profs['ne']['sqrt{psi_n}']['data'])):
        f.write((
            '%0.5f\t%0.5e\t%0.5e\t%0.5e\n'%(
                profs['ne']['sqrt{psi_n}']['data'][ii],
                profs['ne']['n_e, CTS']['data'][ii]*1e20,
                profs['Te']['T_e, CTS']['data'][ii]*1e3,
                profs['Ti']['T_i']['data'][ii]*1e3
                )
            ))

    # Closes file
    f.close()

