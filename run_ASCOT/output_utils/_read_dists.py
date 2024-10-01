'''

Script to read distribution data from ASCOT

cjperks
Sep 26, 2024


'''

# Modules
import os, sys
import numpy as np

import h5py


__all__ = [
    '_read_dist'
    ]

##########################################
#
#               Main
#
##########################################

# Function to open .h5 data and format it
def _read_dist(
    h5_file = None,
    ):

    # Init
    dout = {}

    # Reads file
    ff = h5py.File(h5_file, 'r')

    # Loads distribution data
    try:
        dst = ff['distributions']
    except:
        try:
            dst = ff['dists']
        except:
            print('NO DISTRIBUTION DATA')
            return dout

    # (R,Z) x(mu, E) distribution
    if 'rzPitchEdist' in dst.keys():
        dout['rzDist'] = _read_rzDist(
            rzPE = dst['rzPitchEdist']
            )

    # rho distributions
    if 'rhoDist' in dst.keys():
        dout['rhoDist'] = _read_rhoDist(
            rhoD = dst['rhoDist']
            )

    # B-field data
    dout['bfield'] = _read_bfield(
        bfield = ff['bfield']
        )

    # Init state data
    if 'inistate' in ff.keys():
        dout['inistate'] = _read_inistate(
            inist = ff['inistate']
            )

    # End state data
    if 'endstate' in ff.keys():
        dout['endstate'] = _read_endstate(
            endst = ff['endstate']
            )

    # Output
    return dout

##########################################
#
#               Utils
#
##########################################

# Reads end states
def _read_endstate(
    endst = None,
    ):

    # Init
    dout = {}

    # Loop over keys
    for key in endst.keys():
        dout[key] = {}
        dout[key]['vals'] = endst[key][:]
        dout[key]['units'] = endst[key].attrs['unit']
        dout[key]['quant'] = endst[key].name

    # Ouput
    return dout

# Reads initial states
def _read_inistate(
    inist = None,
    ):

    # Init
    dout = {}

    # Loop over keys
    for key in inist.keys():
        dout[key] = {}
        dout[key]['vals'] = inist[key][:]
        dout[key]['units'] = inist[key].attrs['unit']
        dout[key]['quant'] = inist[key].name

    # Ouput
    return dout

# Reads B-field data
def _read_bfield(
    bfield = None,
    ):

    # Init
    dout = {}

    psiAxis = np.min(bfield['2d']['psi'][:].flatten())
    psiSep = 0
    dout['rhop'] = np.sqrt(
        (bfield['2d']['psi'][:] - psiAxis)
        /(psiSep - psiAxis)
        ) # dim(nR, nZ)
    dout['RR'] = bfield['r'][:] # dim(nR,)
    dout['ZZ'] = bfield['z'][:] # dim(nZ,)

    # Output
    return dout

# Reads rho distributions
def _read_rhoDist(
    rhoD = None,
    ):

    # Init
    dout = {}

    # Volume of each rho bin
    dout['shellVolume'] = {}
    dout['shellVolume']['vals'] = rhoD['shellVolume'][:]
    dout['shellVolume']['quant'] = 'shell Volume'
    dout['shellVolume']['units'] = 'm^3'
    dout['shellVolume']['size'] = len(dout['shellVolume']['vals'])

    # Cross-section area of each rho bin
    dout['shellArea'] = {}
    dout['shellArea']['vals'] = rhoD['shellArea'][:]
    dout['shellArea']['quant'] = 'shell Area'
    dout['shellArea']['units'] = 'm^2'
    dout['shellArea']['size'] = len(dout['shellArea']['vals'])

    # Attributes
    dout['n_1d_dists'] = rhoD.attrs['number_of_1d_dists']
    dout['n_bkg_species'] = rhoD.attrs['number_of_background_species']
    dout['n_deposition_dists'] = rhoD.attrs['number_of_deposition_dists']

    #### --- Edge dimensions --- ###
    dout['edges'] = {}

    # Loop over dimension
    for ii in np.arange(1,3+1):
        dimx = 'dim%i'%(ii)
        dout['edges'][dimx] = {}
        dout['edges'][dimx]['vals'] = rhoD['abscissae'][dimx][:] 
        dout['edges'][dimx]['quant'] = rhoD['abscissae'][dimx].attrs['name']
        dout['edges'][dimx]['units'] = rhoD['abscissae'][dimx].attrs['unit']
        dout['edges'][dimx]['size'] = len(dout['edges'][dimx]['vals'])
        dout['edges'][dimx]['dbin'] = abs(
            np.mean(dout['edges'][dimx]['vals'][1:] - dout['edges'][dimx]['vals'][:-1])
            )

    #### ---  Centers --- ####
    dout['cents'] = {}

    for dimx in dout['edges'].keys():
        dout['cents'][dimx] = {}

        dout['cents'][dimx]['vals'] = (
            dout['edges'][dimx]['vals'][1:] + dout['edges'][dimx]['vals'][:-1]
            )/2
        dout['cents'][dimx]['quant'] = dout['edges'][dimx]['quant']
        dout['cents'][dimx]['units'] = dout['edges'][dimx]['units']
        dout['cents'][dimx]['size'] = len(dout['cents'][dimx]['vals'])
        dout['cents'][dimx]['dbin'] = dout['edges'][dimx]['dbin']


    #### --- Distributions --- ####
    '''
    1 -- density [1/m3]
    2 -- energy density [J/m3]
    3 --- parallel current density [A/m2]
    4 -- toroidal current density [A/m2]
    5 -- jxB torque [N/m2]
    6 -- jxB torque (from ini/endstate) [N/m2]
    7 -- CX ion source [1/m3]
    8 -- CX ion energy source [J/m3]
    9 -- CX neutral source [1/m3]
    10 -- CX neutral energy source [J/m3]
    11 -- Pphi torque [N/m2]
    12 -- parallel energy density [J/m3]
    13 -- total toroidal current density [A/m2]
    14 -- total pressure [N/m2]
    15 -- parallel pressure [N/m2]
    16 -- perpendicular pressure [N/m2]
    17 -- finite Larmor radius torque [N/m2]
    18 -- Thermalized particle density [1/m3]
    19 -- Thermalized particle energy density [J/m3]
    20 -- Thermalized particle torque [N/m2]
    21 -- Absorbed ICRH power density [W/m3]
    22 -- J.B [A*T/m2]
    23 -- Power deposition to electrons [W/m3]
    24 -- Power deposition to background species 1 [W/m3]
    25 -- Power deposition to background species 2 [W/m3]
    26 -- Power deposition to background species 3 [W/m3]
    27 -- Power deposition to background species 3 [W/m3]
    28 -- Collisional torque deposition to electrons [N/m2]
    29 -- Collisional torque deposition to background species 1 [N/m2]
    30 -- Collisional torque deposition to background species 2 [N/m2]
    31 -- Collisional torque deposition to background species 3 [N/m2]
    32 -- Collisional torque deposition to background species 4 [N/m2]
    '''
    dout['dists'] = {}

    for ii in np.arange(int(len(rhoD['ordinates'].keys())/2)):
        dout['dists'][ii+1] = {}

        dout['dists'][ii+1]['vals'] = np.transpose(
            rhoD['ordinate'][0,0,0,:,:,:,ii],
            axes = (2,1,0)
            ) # dim(nrho, nt, nspecies)

        dout['dists'][ii+1]['quant'] = rhoD['ordinates']['name_'+('%i'%(ii+1)).rjust(6, '0')][()]
        dout['dists'][ii+1]['units'] = rhoD['ordinates']['unit_'+('%i'%(ii+1)).rjust(6, '0')][()]

    # Output
    return dout

# Reads (R,Z) x(mu, E) distribution
def _read_rzDist(
    rzPE = None,
    ):

    # Init
    dout = {}

    #### ---  Bin edges --- ####
    dout['edges'] = {}

    # Loop over dimensions
    for ii in np.arange(1,6+1):
        dimx = 'dim%i'%(ii)
        dout['edges'][dimx] = {}
        dout['edges'][dimx]['vals'] = rzPE['abscissae'][dimx][:] 
        dout['edges'][dimx]['quant'] = rzPE['abscissae'][dimx].attrs['name']
        dout['edges'][dimx]['units'] = rzPE['abscissae'][dimx].attrs['unit']
        dout['edges'][dimx]['size'] = len(dout['edges'][dimx]['vals'])
        dout['edges'][dimx]['dbin'] = abs(
            np.mean(dout['edges'][dimx]['vals'][1:] - dout['edges'][dimx]['vals'][:-1])
            )

    #### ---  Centers --- ####
    dout['cents'] = {}

    for dimx in dout['edges'].keys():
        dout['cents'][dimx] = {}

        dout['cents'][dimx]['vals'] = (
            dout['edges'][dimx]['vals'][1:] + dout['edges'][dimx]['vals'][:-1]
            )/2
        dout['cents'][dimx]['quant'] = dout['edges'][dimx]['quant']
        dout['cents'][dimx]['units'] = dout['edges'][dimx]['units']
        dout['cents'][dimx]['size'] = len(dout['cents'][dimx]['vals'])
        dout['cents'][dimx]['dbin'] = dout['edges'][dimx]['dbin']


    #### --- Distribution --- ####

    dout['dist'] = {}
    dout['dist']['vals'] = np.transpose( # Note python's h5py reads data row-major (C-style) whereas the data was stored as column-major (Fortran-style)
        rzPE['ordinate'][:],
        axes = np.arange(np.ndim(rzPE['ordinate'][:]), 0, -1)-1
        )[0][...,-1] # dim(nR, nZ, npitch, nE, nt)

    dout['dist']['quant'] = rzPE['ordinates']['name_000001'][()]
    dout['dist']['units'] = rzPE['ordinates']['unit_000001'][()]
    dout['dist']['ndim'] = rzPE['ordinate'].attrs['number_of_dims']
    dout['dist']['nslots'] = rzPE['ordinate'].attrs['number_of_slots']
    dout['dist']['vec_length'] = rzPE['ordinate'].attrs['vector_length']

    # Output
    return dout