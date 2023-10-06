'''

Writes the input.plasma_1d ASCOT4 input file

cjperks
Oct 2nd, 2023

'''

# Modules
from omfit_classes import omfit_gapy
from datetime import date
import numpy as np
import os

__all__ = [
    'write_plasma1d'
    ]

###########################################
#
#           Main
#
###########################################

def write_plasma1d(
    # path/to/input.gaocde
    in_path = None,
    fgacode = None,
    # Ions to include
    ions = ['D', 'H', 'Ar16', 'Ar17'],
    charges = [1, 1, 16, 17],
    masses = [2, 1, 40, 40],
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

    # Opens an ASCII file to write in
    f = open(in_path+'/input_ASCOT/input.plasma_1d', 'w')

    # Header
    f.write('# Input file for ASCOT containing radial 1D information of plasma temperature,density and toroidal rotation\n')
    f.write('# range must cover [0,1] of normalised poloidal rho. It can exceed 1.')
    f.write("# Created: "+date.today().strftime("%m/%d/%Y")+" for CMOD\n")
    f.write("# Author: Conor Perks\n")

    # Mesh data
    f.write(
        str(len(inputga['ne'])) 
        + "\t"
        + str(len(ions)+1).rjust(2, ' ')
        + "\t"
        + "# Nrad,Nion"
        + "\n"
        )

    # Loop over ions
    Zs = ''
    ms = ''
    cs = '1 '
    ns = ''
    for ii in np.arange(len(ions)):
        Zs += str(charges[ii]).ljust(3, ' ')
        ms += str(masses[ii]).ljust(3, ' ')
        cs += str(1).ljust(2, ' ')
        ns += ('Ni'+str(ii+1)+ "(1/m3)").rjust(16, ' ')

    # F-like Lumped impurity fill for QN
    Zs += str(9).ljust(3, ' ')
    ms += str(19).ljust(3, ' ')
    cs += str(1).ljust(2, ' ')
    ns += ('Ni'+str(len(ions)+1)+ "(1/m3)").rjust(16, ' ')

    # Ion charge
    f.write(
        Zs
        + "\t"
        + "# ion Znum"
        + "\n"
        )

    # Ion mass
    f.write(
        ms
        + "\t"
        + "# ion Amass"
        + "\n"
        )

    # Collision mode
    f.write(
        cs
        + "\t"
        + "# collision mode (0= no colls, 1=Maxw colls, 2=binary colls, 3=both colls) 1st number is for electrons"
        + "\n"
        )

    # Data headers
    f.write(
        'RHO (pol)'.rjust(16, ' ')
        + 'Te (eV)'.rjust(16, ' ')
        + 'Ne (1/m3)'.rjust(16, ' ')
        + 'Vtor_I (rad/s)'.rjust(16, ' ')
        + 'Ti1 (eV)'.rjust(16, ' ')
        + ns
        + "\n"
        )

    # Data

    # Loop over radius
    for rr in np.arange(len(inputga['ne'])):
        ni = ''
        # Loop over ions
        for ion in ions:
            # Loop over ions in input.gacode file
            ion_in_ga = False
            for kion in inputga['IONS'].keys():
                if inputga['IONS'][kion][0] == ion:
                    # Flag that ion is present in input.gacode file
                    ion_in_ga = True
                    if ion == 'D':
                        kionD = str(kion)
                    elif ion == 'H':
                        kionH = str(kion)

                    # Ion density data to use
                    data_nz = inputga['ni_'+str(kion)][rr]*1e19

            # If externally supplied profile
            if not ion_in_ga:
                # Tries to see if an external profile was supplied
                try:
                    data_nz = ext_nz['data'][ext_nz['ion'].index(ion)][rr]*1e6

                # Error message
                except:
                    print('Missing data for ion = '+ion)

            # Write data
            ni += "{:1.7E}".format(data_nz).rjust(16, ' ')

        # Lumped impurity fill
        data_lump = (
            (
                inputga['ne'][rr] - inputga['ni_'+kionD][rr] - inputga['ni_'+kionH][rr]
            )*1e19
            / 9
            )
        ni += "{:1.7E}".format(data_lump).rjust(16, ' ')

        # Writes data
        f.write(
            "{:1.7E}".format(
                np.sqrt(inputga['polflux'][rr]/inputga['polflux'][-1])
                ).rjust(16, ' ')
            + "{:1.7E}".format(
                inputga['Te'][rr]*1e3
                ).rjust(16, ' ')
            + "{:1.7E}".format(
                inputga['ne'][rr]*1e19
                ).rjust(16, ' ')
            + "{:1.7E}".format(
                0.0
                ).rjust(16, ' ')
            + "{:1.7E}".format(
                inputga['Ti_1'][rr]*1e3
                ).rjust(16, ' ')
            + ni
            + "\n"
            )

    # Closes file
    f.close()