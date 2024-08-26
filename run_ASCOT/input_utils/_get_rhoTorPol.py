'''

Script to write input.rhoTorPol file for
ASCOT4-RFOF

cjperks
July 29, 2025

'''

# Modules
import sys, os
import numpy as np
from omfit_classes import omfit_gapy

__all__ = [
    'write_rho'
    ]

#############################
#
#           Main
#
#############################

def write_rho(
    # path/to/input.gaocde
    in_path = None,
    fgacode = None,
    ):

    # Reads input.gacode file
    inputga = omfit_gapy.OMFITgacode(
        os.path.join(
            in_path,
            fgacode,
            )
        )

    # Opens an ASCII file to write in
    f = open(in_path+'/input_ASCOT/input.rhoTorPol', 'w')

    # Header
    f.write('#\n')
    f.write(
        '#'
        + '\trho_pol'
        + '\t\t\trho_tor'
        + '\n'
        )

    # Write data
    for rr in np.arange(len(inputga['polflux'])):
        line = (
            '\t'
            + "{:1.7E}".format(
                np.sqrt(inputga['polflux'][rr]/inputga['polflux'][-1])
                )
            + '\t'
            + "{:1.7E}".format(inputga['rho'][rr])
            )
        if rr < len(inputga['polflux'])-1:
            line += '\n'
        f.write(line)

    # Close file
    f.close()