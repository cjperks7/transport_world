'''

Script to create the TXT files needed to use profile data
in SPIRAL

cjperks
Nov 27, 2024

'''

# Modules
import numpy as np
from omfit_classes import omfit_gapy

__all__ = [
    'p2s'
]

#############################################
#
#           Main
#
#############################################

# Manages density/temperature profiles to make
def p2s(
    outpath = None,
    bulk_ion = 'D_ga',
    fast_ion = 'He3_ga',
    fgacode = None,
    ext_bulk = None, # external bulk density array
    ext_fast = None, # external fast ion density array
    ):

    # Loads data from gacode
    ga = omfit_gapy.OMFITgacode(fgacode)

    # Writes the electron temperature data, [keV]
    with open(outpath+'/Te.txt', 'w') as f:
        for ii in np.arange(len(ga['Te'])):
            f.write(
                '%0.4f %0.4f\n'%(
                    abs(ga['polflux'][ii]/ga['polflux'][-1]),
                    ga['Te'][ii]
                    )
                )

    # Writes the ion temperature data, [keV]
    with open(outpath+'/Ti.txt', 'w') as f:
        for ii in np.arange(len(ga['Te'])):
            f.write(
                '%0.4f %0.4f\n'%(
                    abs(ga['polflux'][ii]/ga['polflux'][-1]),
                    ga['Ti_1'][ii]
                    )
                )

    # Writes the electron density data, [1e19 /m3]
    with open(outpath+'/ne.txt', 'w') as f:
        for ii in np.arange(len(ga['Te'])):
            f.write(
                '%0.4f %0.4f\n'%(
                    abs(ga['polflux'][ii]/ga['polflux'][-1]),
                    ga['ne'][ii]
                    )
                )

    # Writes the bulk ion density data, [1/m3]
    if bulk_ion.split('_')[-1] == 'ga':
        ga_ions = [ga['IONS'][kk][0] for kk in ga['IONS'].keys()]
        ind = ga_ions.index(bulk_ion.split('_')[0])

        dens = ga['ni_%i'%(ind+1)]*1e19 # [1/m3]

    elif bulk_ion.split('_')[-1] == 'ext':
        dens = ext_bulk # [1/m3]

    with open(outpath+'/n%s.txt'%(bulk_ion.split('_')[0]), 'w') as f:
        for ii in np.arange(len(ga['Te'])):
            f.write(
                '%0.4f %0.4e\n'%(
                    abs(ga['polflux'][ii]/ga['polflux'][-1]),
                    dens[ii]
                    )
                )

    # Writes the fast ion density data, [1/m3]
    if fast_ion.split('_')[-1] == 'ga':
        ga_ions = [ga['IONS'][kk][0] for kk in ga['IONS'].keys()]
        ind = ga_ions.index(fast_ion.split('_')[0])

        dens = ga['ni_%i'%(ind+1)]*1e19 # [1/m3]

    elif fast_ion.split('_')[-1] == 'ext':
        dens = ext_fast # [1/m3]

    with open(outpath+'/n%s.txt'%(fast_ion.split('_')[0]), 'w') as f:
        for ii in np.arange(len(ga['Te'])):
            f.write(
                '%0.4f %0.4e\n'%(
                    abs(ga['polflux'][ii]/ga['polflux'][-1]),
                    dens[ii]
                    )
                )
