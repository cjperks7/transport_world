'''

Script to read .cdf files prepared by ACTION_profiles in SPIRAL

cjperks
Nov 21, 2024

'''

# Modules
import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4
import scipy.constants as cnt

plt.rcParams.update({'font.size': 16})

sys.path.insert(0,'/home/cjperks/usr/python3modules/eqtools3')
import eqtools
sys.path.pop(0)

__all__ = [
    'read_prf',
    'cmp_prf_inputs'
    ]

###############################################
#
#           Main
#
###############################################

# Reads and organizes .cdf data
def read_prf(
    prf = None,
    ):

    # Init
    ddata = {}
    ff = netCDF4.Dataset(prf, 'r')

    # Gets mesh
    ddata['r_1d'] = {}
    ddata['r_1d']['data'] = ff.variables['profile_r_axis'][:] # dim(nr,)
    ddata['r_1d']['units'] = 'm'

    ddata['z_1d'] = {}
    ddata['z_1d']['data'] = ff.variables['profile_z_axis'][:] # dim(nz,)
    ddata['z_1d']['units'] = 'm'

    # Calculate the bin edges
    ddata['r_edges'] = {}
    bin_width = ddata['r_1d']['data'][1] - ddata['r_1d']['data'][0]
    ddata['r_edges']['data'] = np.concatenate((
        [ddata['r_1d']['data'][0] - bin_width / 2],
        ddata['r_1d']['data'] + bin_width / 2,
        ))
    ddata['z_edges'] = {}
    bin_width = ddata['z_1d']['data'][1] - ddata['z_1d']['data'][0]
    ddata['z_edges']['data'] = np.concatenate((
        [ddata['z_1d']['data'][0] - bin_width / 2],
        ddata['z_1d']['data'] + bin_width / 2,
        ))


    ddata['r_2d'] = {}
    ddata['z_2d'] = {}
    ddata['r_2d']['units'] = ddata['z_2d']['units'] = 'm'
    ddata['z_2d']['data'], ddata['r_2d']['data'] = np.meshgrid(
        ddata['z_1d']['data'], ddata['r_1d']['data']
        )

    ddata['phi_1d'] = {}
    ddata['phi_1d']['data'] = ff.variables['profile_phi_axis'][:] # dim(nphi,)
    ddata['phi_1d']['units'] = 'deg'

    # Gets fast ion density
    ddata['density'] = {}
    ddata['density']['data'] = ff.variables['profile_fast-ion_density'][:] # dim(nphi, nz,nr)
    ddata['density']['units'] = '1/m^3'

    # Gets fast ion effective temperature
    ddata['T_eff'] = {}
    ddata['T_eff']['data'] = ff.variables['profile_fast-ion_effective_temperature'][:]/cnt.e/1e3 # dim(nphi, nz,nr)
    ddata['T_eff']['units'] = 'keV'

    # Gets fast ion parallel temperature
    ddata['T_para'] = {}
    ddata['T_para']['data'] = ff.variables['profile_fast-ion_parallel_effective_temperature'][:]/cnt.e/1e3 # dim(nphi, nz,nr)
    ddata['T_para']['units'] = 'keV'

    # Gets fast ion perpendicular temperature
    ddata['T_perp'] = {}
    ddata['T_perp']['data'] = ff.variables['profile_fast-ion_perpendicular_effective_temperature'][:]/cnt.e/1e3 # dim(nphi, nz,nr)
    ddata['T_perp']['units'] = 'keV'

    # Gets voxel volume
    ddata['volume'] = {}
    ddata['volume']['data'] = ff.variables['profile_volume'][:] # dim(nphi, nz,nr)
    ddata['volume']['units'] = 'm^3'

    # Gets mode exchange energy
    ddata['mode_energy'] = {}
    ddata['mode_energy']['data'] = ff.variables['profile_mode_energy'][:] # dim(nphi, nz,nr)
    ddata['mode_energy']['units'] = 'W'

    # Gets magnetic axis data
    ddata['mag_axis'] = {}
    ddata['mag_axis']['R'] = float(ff.variables['r_axis'][:])
    ddata['mag_axis']['Z'] = float(ff.variables['z_axis'][:])

    # Gets total energy density
    ddata['energy_density'] = {}
    ddata['energy_density']['data'] = ff.variables['profile_fast-ion_pressure'][:] # dim(nphi, nz,nr)
    ddata['energy_density']['units'] = 'J/m^3'

    # Simulated time
    #ddata['dt'] = 3e-3 # [s]
    ddata['dt'] = 1e-4 # [s]

    # Gets heat flux
    ddata['heat_flux'] = {}
    ddata['heat_flux']['data'] = np.zeros_like(ddata['mode_energy']['data'])
    #ddata['heat_flux']['data'] = ddata['energy_density']['data']/ddata['dt']
    ddata['heat_flux']['units'] = 'W/m^3'

    
    # Output
    return ddata


###############################################
#
#           Plotting
#
###############################################

# Compares ACTION_profiles data against inputs
def cmp_prf_inputs(
    ddata=None,
    ):

    a=9