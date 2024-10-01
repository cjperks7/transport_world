'''

_write_gacode.py is a script that prepares a input file
in the GACODE format

'''

import numpy as np
from datetime import date
import os

__all__ = [
    'write_ga',
]

##########################################################
#
#                     Main
#
##########################################################

def write_ga(
    dout=None,
    name=None,
    ):

    # Opens an ASCII file to write in
    if name is None:
        f = open(
            os.path.join(
                dout['paths']['input'],
                'input_%s.gacode'%(dout['t0_s'])
                ),
            'w'
            )
    else:
        f=open(name, 'w')

    # Header
    f.write("# \n")
    f.write("# Created: "+date.today().strftime("%m/%d/%Y")+" for CMOD\n")
    f.write("# Author: Conor Perks\n")
    f.write("# \n")

    # Number of radial grid points
    f.write("# nexp\n")
    f.write(dout['nexp'] + "\n")

    # Number of ions
    f.write("# nion\n")
    f.write(dout['ions']['nion'] + "\n")

    # Shot number
    f.write("# shot\n")
    f.write(str(dout['shot']) + "\n")

    # Name of ions
    f.write("# name\n")
    f.write(dout['ions']['name'] + "\n")

    # Type of ions
    f.write("# type\n")
    f.write(dout['ions']['itype'] + "\n")

    # Mass of an electron in amu
    f.write("# masse\n")
    f.write(dout['e']['mass'] + "\n")

    # Mass of ions in amu
    f.write("# mass\n")
    f.write(dout['ions']['mass'] + "\n")

    # Charge of an electron
    f.write("# ze\n")
    f.write(dout['e']['charge'] + "\n")

    # Charge of ions
    f.write("# z\n")
    f.write(dout['ions']['charge'] + "\n")

    # Toroidal flux at the edge divided by 2pi
    f.write("# torfluxa | Wb/radian\n")
    f.write(dout['torfluxa_Wb/rad'] + "\n") 

    # Radial location of the plasma center
    f.write("# rcentr | m\n")
    f.write(dout['rcentr_m'] + "\n")

    # B-field on axis
    f.write("# bcentr | T\n")
    f.write(dout['bcentr_T'] + "\n")

    # Plasma current
    f.write("# current | MA\n")
    f.write(dout['current_MA'] + "\n")

    # Radial flux coordinate, sq. tor. flux
    f.write("# rho | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['rhot'][ii]).rjust(15, ' ') 
            + "\n"
            )
    
    # Minor radius
    f.write("# rmin | m\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['rmin_m'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Poloidal flux divided by 2pi
    f.write("# polflux | Wb/radian\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['polflux_Wb/rad'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # q profile
    f.write("# q | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['q'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Rotation frequency
    f.write("# omega0 | rad/s\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['omega0_rad/s'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Major radius
    f.write("# rmaj | m\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['rmaj_m'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Elevation
    f.write("# zmag | m\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['zmag_m'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Elongation
    f.write("# kappa | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['kappa'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Triangularity
    f.write("# delta | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['delta'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Squareness
    f.write("# zeta | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['zeta'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Higher flux-surface tilt moment
    f.write("# shape_cos0 | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['shape_cos0'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Higher flux-surface tilt moment
    f.write("# shape_cos1 | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['shape_cos1'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Higher flux-surface tilt moment
    f.write("# shape_cos2 | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['shape_cos2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Higher flux-surface tilt moment
    f.write("# shape_cos3 | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['shape_cos3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Higher flux-surface tilt moment
    f.write("# shape_sin3 | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['shape_sin3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Electron density
    f.write("# ne | 10^19/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['ne_19m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Ion Density
    f.write("# ni | 10^19/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        for jj,ion in enumerate(list(filter(None,dout['ions']['name'].split(' ')))):
            if jj == 0:
                f.write(
                    str(ii+1).rjust(3, ' ') 
                    + "{:1.7E}".format(dout['ions'][ion]['ni_tot_19m3'][ii]).rjust(15, ' ')
                    )
            elif jj == (int(dout['ions']['nion']) -1):
                f.write(
                    "{:1.7E}".format(dout['ions'][ion]['ni_tot_19m3'][ii]).rjust(15, ' ') 
                    + "\n"
                    )
            else:
                f.write(
                    "{:1.7E}".format(dout['ions'][ion]['ni_tot_19m3'][ii]).rjust(15, ' ')
                    )
    
    # Electron temperature
    f.write("# te | keV\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['Te_keV'][ii]).rjust(15, ' ') 
            + "\n"
            )
    
    # Ion temperature
    f.write("# ti | keV\n")
    for ii in np.arange(int(dout['nexp'])):
        for jj,ion in enumerate(list(filter(None,dout['ions']['name'].split(' ')))):
            if jj == 0:
                f.write(
                    str(ii+1).rjust(3, ' ') 
                    + "{:1.7E}".format(dout['Ti_keV'][ii]).rjust(15, ' ')
                    )
            elif jj == (int(dout['ions']['nion']) -1):
                f.write(
                    "{:1.7E}".format(dout['Ti_keV'][ii]).rjust(15, ' ') 
                    + "\n"
                    )
            else:
                f.write(
                    "{:1.7E}".format(dout['Ti_keV'][ii]).rjust(15, ' ')
                    )

    # Toroidal rotation
    f.write("# vtor | m/s\n")
    for ii in np.arange(int(dout['nexp'])):
        for jj,ion in enumerate(list(filter(None,dout['ions']['name'].split(' ')))):
            if jj == 0:
                f.write(
                    str(ii+1).rjust(3, ' ') 
                    + "{:1.7E}".format(dout['vtor_m/s'][ii]).rjust(15, ' ')
                    )
            elif jj == (int(dout['ions']['nion']) -1):
                f.write(
                    "{:1.7E}".format(dout['vtor_m/s'][ii]).rjust(15, ' ') 
                    + "\n"
                    )
            else:
                f.write(
                    "{:1.7E}".format(dout['vtor_m/s'][ii]).rjust(15, ' ')
                    )

    # Poloidal rotation
    f.write("# vpol | m/s\n")
    for ii in np.arange(int(dout['nexp'])):
        for jj,ion in enumerate(list(filter(None,dout['ions']['name'].split(' ')))):
            if jj == 0:
                f.write(
                    str(ii+1).rjust(3, ' ') 
                    + "{:1.7E}".format(dout['vpol_m/s'][ii]).rjust(15, ' ')
                    )
            elif jj == (int(dout['ions']['nion']) -1):
                f.write(
                    "{:1.7E}".format(dout['vpol_m/s'][ii]).rjust(15, ' ') 
                    + "\n"
                    )
            else:
                f.write(
                    "{:1.7E}".format(dout['vpol_m/s'][ii]).rjust(15, ' ')
                    )

    # Plasmas pressure
    f.write("# ptot | Pa\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['ptot_Pa'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Ohmic current density
    f.write("# johm | MA/m^2\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['johm_MA/m2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Parallel bootstrap current density
    f.write("# jbs | MA/m^2\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['jbs_MA/m2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # RF-driven current density
    f.write("# jrf | MA/m^2\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['jrf_MA/m2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Beam-driven current density
    f.write("# jnb | MA/m^2\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['jnb_MA/m2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Toroidal bootstrap current density
    f.write("# jbstor | MA/m^2\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['jbstor_MA/m2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Z_eff
    f.write("# z_eff | -\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['Zeff']['prof'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Ohmic power density
    f.write("# qohme | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['ohm']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Beam power density to electron
    f.write("# qbeame | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['beame']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Beam power density to ions
    f.write("# qbeami | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['beami']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # RF power density to electrons
    f.write("# qrfe | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['rfe']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # RF power density to ions
    f.write("# qrfi | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['rfi']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Fusion alpha power density to electrons
    f.write("# qfuse | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['alpe']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Fusion alpha power density to ions
    f.write("# qfusi | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['alpi']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Synchrotron radiation power density
    f.write("# qsync | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['sync']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Bremsstrahlung radiation power density
    f.write("# qbrem | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['cont']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Line radiation power density
    f.write("# qline | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['line']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Electron-ion energy exchange power density
    f.write("# qei | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['ei']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Recombination power density to electrons
    f.write("# qione | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['ione']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Recombination power density to ions
    f.write("# qioni | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['ioni']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Charge exchange power density to ions
    f.write("# qcxi | MW/m^3\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['cxi']['prof_MW/m3'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Beam particle source density
    f.write("# qpar_beam | 1/m^3/s\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['par_beam']['prof_1/m3/s'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Beam particle source density
    f.write("# qpar_wall | 1/m^3/s\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['par_wall']['prof_1/m3/s'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Momentum source density
    f.write("# qmom | N/m^2\n")
    for ii in np.arange(int(dout['nexp'])):
        f.write(
            str(ii+1).rjust(3, ' ') 
            + "{:1.7E}".format(dout['powers']['mom']['prof_N/m2'][ii]).rjust(15, ' ') 
            + "\n"
            )

    # Closes file
    f.close()

