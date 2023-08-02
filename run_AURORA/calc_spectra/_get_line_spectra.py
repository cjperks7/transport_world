'''

_get_spectra.py is a module that facilitates running AURORA to calculate
photon radiation emissiion spectra

NOTE: Outputs flux-surface-averaged spectra


'''

import numpy as np
import aurora
import os
from omfit_classes import omfit_gapy
from scipy.interpolate import interp1d

__all__ = [
    'get_line_spectra',
    ]

#########################################################
#
#               Main
#
#########################################################

def get_line_spectra(
    ddiag = None,   # Dictionary storing diagnostic information
    fgacode = None, # Background kinetic profiles
    dnz = None,     # Dictionary of trace impurity densities
    imp_res = False,    # Impurity-resovled spectra
    plt_all = None,
    ):
    '''
    NOTE: The philosophy used here is that the user would externally
    loop over impurity species when given density profiles

    '''

    # Loads background kinetic profiles
    inputga =  omfit_gapy.OMFITgacode(fgacode)

    # Master wavelength domain
    if ddiag['spec_rng']['units'] == 'AA':
        wave_master = np.linspace(
            ddiag['spec_rng']['rng'][0],
            ddiag['spec_rng']['rng'][1],
            int(ddiag['spec_rng']['nlamb'])
            )
    else:
        print('OTHER DIAG UNITS NOT IMPLEMENTED')

    # Initializes outs
    spec = np.zeros((len(inputga['rho']), len(wave_master))) # dim(rhop, lambda_fm)
    spec_imp = {}

    # Loop over ion species of interest
    for ion in ddiag['imps'].keys():
        spec_imp[ion] = {}

        # Loop over charge states of interest
        for cs in ddiag['imps'][ion]:
            # Obtains the synthetic spectra for this species/charge state
            spec_tmp = _get_line_emis(
                inputga=inputga,
                wave_master=wave_master,
                ion = ion,
                cs = cs,
                nz = dnz[ion],
                ) # [1/cm3/s/AA], dim(rhop,lambda)

            # Store output
            spec += spec_tmp
            spec_imp[ion][cs] = spec_tmp

    # Plotting
    if plt_all:
        from transport_world.plot_utils import plt_slider

        if imp_res:
            dzzz = spec_imp
        else:
            dzzz = spec

        plt_slider(
            xxx=wave_master,
            yyy=inputga['rho'],
            dzzz=dzzz,
            xxlabel=r"$\lambda$ [$\AA$]",
            yylabel=r"$\rho_t$",
            zzlabel=r"Emissivity [$1/cm^3/s/\AA$]",
            plt_sum=imp_res
            )

    # Output
    if imp_res:
        return {
            'lamb_AA': wave_master,
            'tot': spec,
            'imp_res': spec_imp,
            }
    else:
        return {
            'lamb_AA': wave_master,
            'tot': spec,
            }


#########################################################
#
#               Utilities
#
#########################################################

def _get_line_emis(
    inputga=None,
    wave_master=None,
    ion=None,
    cs=None,
    nz=None,
    ):

    # Obtains file location of PECs
    filepath = _get_PEC_loc(
        ion=ion,
        cs=cs,
        )

    # Intializes array to store spectrum data
    spec_fm = np.zeros((len(inputga['rho']), len(wave_master))) # dim(rhop, lambda_fm)

    print(ion)
    print(cs)
    # Loop over radial location
    for rho in np.arange(len(inputga['rho'])):
        # If possible, includes ionization, excitation, and recombination
        print(rho)
        if cs != 0:
            # Calculates the synthetic line radiation spectrum for the jth charge state
            out = aurora.radiation.get_local_spectrum(
                filepath, 
                inputga['ne'][rho]*1e13, 
                inputga['Te'][rho]*1e3, 
                ion_exc_rec_dens=[
                    nz[rho, cs-1],
                    nz[rho,cs],
                    nz[rho,cs+1]
                    ],
                Ti_eV = inputga['Ti_1'][rho]*1e3, 
                n0_cm3=0.0,
                dlam_A=0.0, 
                plot_spec_tot=False
                ) 

        else:
            # Calculates the synthetic line radiation spectrum for the jth charge state
            out = aurora.radiation.get_local_spectrum(
                filepath, 
                inputga['ne'][rho]*1e13, 
                inputga['Te'][rho]*1e3,
                ion_exc_rec_dens=[
                    0,
                    nz[rho,cs],
                    nz[rho,cs+1]
                    ],
                Ti_eV = inputga['Ti_1'][rho]*1e3,
                n0_cm3=0.0,
                dlam_A=0.0, 
                plot_spec_tot = False
                )
                
        # Sums over all interaction types
        spec_tmp = np.zeros(len(out[0]))
        for nn in range(1,6):
            if out[nn] is not None:
                spec_tmp += out[nn]

        # Forward Model total spectrum interpolated on master wavelength domain
        spec_fm[rho,:] = interp1d(
            out[0], 
            spec_tmp,
            bounds_error = False, fill_value= 0.0
            )(wave_master) # [1/cm^3/s/AA]; dim(,lambda)
        # Sums of ion, exc, rad rec, diel rec, and cx spectra

    # Outputs spectra for this species/charge state
    return spec_fm

def _get_PEC_loc(
    ion=None,
    cs=None
    ):

    # path to the common folder
    path_data = os.path.join(
        '/home/cjperks/tofu_sparc/'
        'atomic_data',
        'ADAS_PEC_files'
        )

    # File path to the PEC ADF15 data file for the iith  charge state
    if ion == 'W':
        #filepath = aurora.adas_files.get_adas_file_loc('pec40#w_ic#w' + str(cs) +'.dat', filetype='adf15')
        filepath = os.path.join(
            path_data,
            'Tungsten',
            'pec40#w_ic#w' + str(cs) +'.dat',
            )
    elif ion == 'Xe':
        filepath = os.path.join(
            path_data,
            'Xenon',
            'fs#0.50A_8.50A#xe' + str(cs) +'.dat',
            )
    elif ion == 'Kr':
        filepath = os.path.join(
            path_data,
            'Krypton',
            'fs#0.50A_8.50A#kr' + str(cs) +'.dat',
            )
    elif ion == 'H':
        filepath = aurora.adas_files.get_adas_file_loc('pec96#h_pju#h' + str(cs) +'.dat', filetype='adf15')
        #filepath = '/home/cjperks/XraySpec/PEC_files/Hydrogen/pec12#h_pju#h' + str(ii) +'.dat'
    elif ion == 'He':
        filepath = aurora.adas_files.get_adas_file_loc('pec96#he_pju#he' + str(cs) +'.dat', filetype='adf15')
        #filepath = '/home/cjperks/XraySpec/PEC_files/Helium/pec96#he_pju#he' + str(ii) +'.dat'
    elif ion == 'Ar':
        filepath = os.path.join(
            path_data,
            'fsciortino_share',
            'atomdb',
            'pec#ar' + str(cs) +'.dat',
            )
    elif ion == 'Mo':
        filepath = os.path.join(
            path_data,
            'Molybdenum',
            'transport_llu#mo' + str(cs) +'ic.dat',
            )
    #elif ion == 'F' and ii != 9:
    #    filepath = '/home/cjperks/XraySpec/PEC_files/Flourine/pec12#h_pju#h' + str(ii) +'.dat'
    elif ion == 'Ne':
        filepath = aurora.adas_files.get_adas_file_loc('pec96#ne_pju#ne' + str(cs) +'.dat', filetype='adf15')

    return filepath