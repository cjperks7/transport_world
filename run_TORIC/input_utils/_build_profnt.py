'''

_build_profnt.py is a module meant to take an input.gacode
file and reformat the density/temperature data as an 
ASCII file appropriate for TORIC5

cjperks
Sept. 19, 2023

'''

# Modules
import numpy as np
import os
from omfit_classes import omfit_gapy

__all__ = [
    'build_profnt',
    ]

###################################################
#
#                   Main
#
###################################################

def build_profnt(
    # path/to/input.gaocde
    in_path = None,
    fgacode = None,
    # Ions to include
    ions = ['D', 'H', 'Ar16', 'Ar17'],
    charges = [1, 1, 16, 17],
    masses = [2, 1, 40, 40],
    # Ion temperature controls
    diff_temp = 0,
    # Optional dictioanry to overwrite nz profile with external
    ext_nz = {
        'ion': ['xAr16', 'xAr17'],
        'data': [
            np.array([0]), # [cm^-3], dim(nrho,)
            np.array([0]) # [cm^-3], dim(nrho,)
            ]
        }
    ):

    # Reads input.gacode file
    inputga = omfit_gapy.OMFITgacode(
        os.path.join(
            in_path,
            fgacode,
            )
        )

    # Assumes main ion species is always deuterium
    mainsp = np.where('D' in ions)[0][0] +1

    # Number of grid points
    nrho = len(inputga['ne'])

    # Number of ions
    nions = len(ions)

    # Opens an ASCII file to write in
    f = open(in_path+'/profnt_TORIC.dat', 'w')

    # Header
    f.write(
        'Profiles'.rjust(10, ' ')                   # Variable name
        + str(nrho).rjust(4, ' ')     # Number of grid points
        + str(nions).rjust(4, ' ')              # Number of ions
        + str(mainsp).rjust(4, ' ')                 # Index of main ion
        + str(1).rjust(4, ' ')                      # Flag if different ion density shapes
        + str(diff_temp).rjust(4, ' ')              # Flag if different ion temperatures
        + "\n"
        )

    # Writes ion charge & mass
    for ii in np.arange(len(ions)):
        f.write(
            str(masses[ii]).rjust(4, ' ')
            + str(charges[ii]).rjust(4, ' ')
            + "\n"
            )
    
    # Writes mesh header
    f.write(
        'rho=sq.psi'.rjust(10, ' ')
        + "\n"
        )

    # Writes mesh data
    _write_block(
        f=f,
        nrho=nrho,
        data = np.sqrt(inputga['polflux']/inputga['polflux'][-1]),
        )

    # Writes header for electron density [cm^-3]
    f.write(
        'ne, cc'.rjust(10, ' ')
        + "\n"
        )

    # Writes electron density data
    _write_block(
        f=f,
        nrho=nrho,
        data = inputga['ne']*1e13,
        )

    # Writes header for electron temperature [keV]
    f.write(
        'Te, keV'.rjust(10, ' ')
        + "\n"
        )

    # Writes electron temperature data
    _write_block(
        f=f,
        nrho=nrho,
        data = inputga['Te'],
        )

    # Loop over ions
    Ti_row = False
    for ii,ion in enumerate(ions):
        for kion in inputga['IONS'].keys():
            if inputga['IONS'][kion][0] == ion:
                # Writes header for ion density [cm^-3]
                f.write(
                    ('ni_'+str(ii+1)+', cc').rjust(10, ' ')
                    + "\n"
                    )

                # Ion density data to use
                try ext_nz['ion'].index(ion):
                    data_nz = ext_nz['data'][ext_nz['ion'].index(ion)]
                except:
                    data_nz = inputga['ni_'+str(kion)]*1e13

                # Write ion temperature data
                _write_block(
                    f=f,
                    nrho = nrho,
                    data = data_nz,
                    )

                # If ion temperature not yet written
                if not Ti_row:
                    # Writes header for ion density [cm^-3]
                    f.write(
                        ('Ti_'+str(ii+1)+', keV').rjust(10, ' ')
                        + "\n"
                        )

                    # Write ion temperature data
                    _write_block(
                        f=f,
                        nrho = nrho,
                        data = inputga['Ti_'+str(kion)],
                        )

                    # Switchs flag to stop writing ion temperature block
                    if diff_temp == 0:
                        Ti_row = True

    # Finish file writing         
    f.close()


def _write_block(
    f=None,
    nrho=None,
    data = None,
    ):

    # Writes block data
    for rr in np.arange(int(np.ceil(nrho/5))):
        row = ''

        for ii in np.arange(5):
            try:
                row += (
                    "{:1.9E}".format(
                        data[5*rr + ii]
                        ).rjust(16, ' ')
                    )
            except:
                blah = 0

        f.write(row+"\n")

