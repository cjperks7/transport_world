'''

Functions to manage my eqtools scripts

cjperks
Nov 27, 204

'''

# Modules
import numpy as np
import sys, os
from scipy.spatial import ConvexHull, Delaunay

sys.path.insert(0,'/home/cjperks/usr/python3modules/eqtools3')
import eqtools
sys.path.pop(0)

__all__ = [
    '_get_eq',
    '_trim_hull',
    ]

########################################
#
#           Main
#
########################################

# Loads eq
def _get_eq(
    gfile = None,
    afile = None,
    machine = 'CMOD',
    ):

    # Reads gfile
    edr = eqtools.EqdskReader(
        gfile=gfile,
        afile=None
        )

    # Gathers data
    dedr = {}

    dedr['rGrid_1d'] = edr.getRGrid()
    dedr['zGrid_1d'] = edr.getZGrid()
    dedr['wall_R'] = edr.getMachineCrossSection()[0]
    dedr['wall_Z'] = edr.getMachineCrossSection()[1]


    dedr['rGrid_2d'], dedr['zGrid_2d'] = np.meshgrid(dedr['rGrid_1d'], dedr['zGrid_1d'])

    dedr['RLCFS'] = edr.getRLCFS()[0]
    dedr['ZLCFS'] = edr.getZLCFS()[0]
    dedr['Rmag'] = edr.getMagR()[0]
    dedr['Zmag'] = edr.getMagZ()[0]

    # Limits
    dedr['Rlim'] = [0.9*np.min(dedr['RLCFS']), 1.1*np.max(dedr['RLCFS'])]
    dedr['Zlim'] = [1.1*np.min(dedr['ZLCFS']), 1.1*np.max(dedr['ZLCFS'])]

    edr_psiRZ = edr.getFluxGrid()[0]
    edr_psiLCFS = edr.getFluxLCFS()[0]
    edr_psi0 = edr.getFluxAxis()[0]

    print(min(edr_psiRZ.flatten()))
    print(edr_psi0)
    print(edr_psiLCFS)
    if machine == 'CMOD':
        dedr['rhop2D'] = np.sqrt(-1*(edr_psiRZ+edr_psi0)/(edr_psiLCFS-edr_psi0))
    elif machine == 'SPARC':
        dedr['rhop2D'] = np.sqrt(-1*edr_psiRZ/edr_psiLCFS)
    elif machine == 'WEST':
        dedr['rhop2D'] = np.sqrt((edr_psiRZ-edr_psi0)/(edr_psiLCFS-edr_psi0))

    dedr['psiLCFS'] = edr_psiLCFS

    # Prepare interpolating 1D rho grid to 2D (R,Z) grid
    xx = np.linspace(0,1,edr_psiRZ.shape[1])
    yy = np.linspace(0,1,edr_psiRZ.shape[0])
    XX, YY = np.meshgrid(xx,yy)
    points = np.vstack((XX.ravel(), YY.ravel())).T
    dedr['rhop1D'] = dedr['rhop2D'].ravel()

    # Work around without afile
    if afile is None:
        edr._time = np.r_[1.]
        edr._RmidLCFS = np.r_[dedr['RLCFS'][np.argmin(abs(dedr['ZLCFS']))]]

    # Output
    return dedr, edr

########################################
#
#           Utilities
#
########################################

# Trims (R,Z) data outside rho>1
def _trim_hull(
    dedr = None,
    RR_2d = None,
    ZZ_2d = None,
    data_2d = None,
    ):

    # Exclude points outside the Convex hull of the LCFS
    LCFS = np.c_[dedr['RLCFS'], dedr['ZLCFS']]
    hull = ConvexHull(LCFS)
    hull_delaunay = Delaunay(LCFS[hull.vertices])

    for ii in np.arange(RR_2d.shape[0]):
        for jj in np.arange(RR_2d.shape[1]):
            if not hull_delaunay.find_simplex(np.c_[RR_2d[ii,jj], ZZ_2d[ii,jj]]) >= 0:

                data_2d[ii,jj] = np.nan

    # Output
    return data_2d

