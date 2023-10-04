'''

Writes 2d wall data from eqdsk

cjperks
Oct 2nd, 2023

'''

# Modules
from omfit_classes import omfit_eqdsk
import os
import numpy as np

__all__ = [
    'write_wall2d'
    ]

###########################################
#
#           Main
#
###########################################

def write_wall2d(
    # path/to/input.gfile
    in_path = None,
    fgfile = None,
    ):

    # Loads gfile
    gfile = omfit_eqdsk.OMFITgeqdsk(
        os.path.join(in_path,fgfile)
        )

    # Opens an ASCII file to write in
    f = open(in_path+'/input_ASCOT/input.wall_2d', 'w')

    # Header
    f.write(
        str(len(gfile['RLIM']))
        + ' (R,z) wall points & divertor flag (1 = divertor, 0 = wall)'
        + '\n'
        )
    
    # Loop over mesh point
    for ii in np.arange(len(gfile['RLIM'])):
        f.write(
            "{:1.7E}".format(gfile['RLIM'][ii]).rjust(16, ' ')
            + "{:1.7E}".format(gfile['ZLIM'][ii]).rjust(16, ' ')
            + "0".rjust(4, ' ')
            + "\n"
        )

    # Closes file
    f.close()

    