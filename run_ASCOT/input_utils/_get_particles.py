'''

Writes particles distribution file

cjperks
Oct 2nd, 2023

'''

# Modules
import numpy as np
from omfit_classes import omfit_gapy, omfit_eqdsk
import aurora
import os

__all__ = [
    'write_particles'
    ]

#################################################
#
#               Main
#
#################################################

def write_particles(
    # Files
    in_path = None,
    fgacode = None,
    fgfile = None,
    # Phase space mesh
    nE = 1,
    npitch = 1,
    nrho = 1,
    ntheta = 1,
    # Maxwellian to init RFOF
    Max = True,
    # Ion data
    ion = 'Ar',
    cs = 16,
    mass = 40,
    nz = None,
    ):

    # Number of markers
    nmark = int(nrho*ntheta*nE*npitch)

    # Reads input files
    inputga = omfit_gapy.OMFITgacode(
        os.path.join(in_path,fgacode)
        )
    gfile = omfit_eqdsk.OMFITgeqdsk(
        os.path.join(in_path,fgfile)
        )

    # Poloidal angle mesh
    theta = np.linspace(0, 2*np.pi, ntheta)

    # Radial mesh
    rhop = np.linspace(0,1, nrho)

    # Pitch mesh
    pitch = np.linspace(-1,1,npitch)

    # Energy mesh
    if Max:
        energy = _get_Max(
            Ti = inputga['Ti_1'],
            mass = mass,
            nE = nE,
            rhop = rhop,
            ) # dim(nrho, nE)
    else:
        print('NOT YET IMPLEMENTED')


    # R,Z flux-surface contours
    RV, ZV = aurora.rhoTheta2RZ(gfile, rhop, theta, coord_in="rhop", n_line=201)
    RV, ZV = RV.T, ZV.T # dim(fm_rhop, fm_theta); [m]

#################################################
#
#               Utilities
#
#################################################

def _get_Max(
    Ti=None,
    mass=None,
    nE=None,
    rhop=None,
    ):


    blah = 0
