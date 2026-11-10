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
    'get_PEC',
    'get_PEC_all',
    'get_sticks'
    ]

#########################################################
#
#               Main
#
#########################################################

'''

import os
from transport_world.run_AURORA import calc_spectra as csp

#Te_eV = 3.0e3
#Te_eV = 5e2
Te_eV = 0.4e3
ne_cm3 = 1e14

#Te_eV = 9e3
#ne_cm3 = 1e14

#sp = 'W'
#Znuc = 74
sp = 'Mo'
Znuc = 42
#sp = 'Ag'
#Znuc = 47
#sp = 'Se'
#Znuc = 34
#sp = 'Ge'
#Znuc = 32
#sp = 'Ni'
#Znuc = 28
#sp = 'Ar'
#Znuc = 18
#sp = 'Xe'
#Znuc = 54
#sp = 'Kr'
#Znuc = 36

#fol = '/home/cjperks/tofu_sparc/atomic_data/FAC_PEC_files/W/XRSLR_3D'
#fol = '/home/cjperks/tofu_sparc/atomic_data/FAC_PEC_files/W/XRSHRXe'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/XRS_LR'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/XRS_HR'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/LoWEUS/%s'%(sp)
fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/XEUS/%s'%(sp)
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/HIREXSR'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/EBIT'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/WEST/%s'%(sp)
#fol = '/home/cjperks/tofu_sparc/atomic_data/ADAS_PEC_files/Tungsten'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/EBIT_Lshell'
#fol = '/home/cjperks/work/papers/EBIT26/issues/neonlike/Se_nmax10'
#fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/all_Lshell'

files = [
    #os.path.join(fol, 'FAC_W%02i.dat'%(ii))
    #for ii in range(6,10)#,19)
    #os.path.join(fol, 'FAC_Mo%02i.dat'%(ii))
    #for ii in range(11,13)
    #os.path.join(fol, 'FAC_Ar%02i.dat'%(ii))
    #for ii in range(4,7)
    #os.path.join(fol, 'FAC_Ar%02i.dat'%(ii))
    #for ii in range(1,4)
    #os.path.join(fol, 'FAC_W%02i.dat'%(ii))
    #for ii in range(13,15)
    os.path.join(fol, 'FAC_%s%02i.dat'%(sp,ii))
    #os.path.join(fol, 'pec40#%s_ic#%s%02i.dat'%(sp.lower(),sp.lower(),Znuc-ii))
    #for ii in range(1,7)
    #for ii in range(11,13)
    #for ii in range(29,31)
    for ii in range(20,24)
    #for ii in range(13,14+1)
    #for ii in range(10,18+1)
    #for ii in range(10,10+1)
    #for ii in range(1,3+1)
    #for ii in range(1,3+1)
    #for ii in range(10,11)
    #for ii in range(2,3+1)
    #for ii in range(6,12+1)
    #for ii in [3,7]
    ]

csp.get_sticks(files=files, Te_eV=Te_eV, ne_cm3=ne_cm3)

'''

# Get PEC sticks (because I always forget)
def get_sticks(
    files=None, # string or list of strings
    Te_eV=None, # [scalar]
    ne_cm3=None, # [scalar]
    ):

    aurora.adf15_line_identification(files, Te_eV=Te_eV, ne_cm3=ne_cm3)

# Get PEC curve of all PECs in dataset
def get_PEC_all(
    # Ion data
    sp = None,
    cs = None,
    path_PEC = None,
    # Plasma data
    Te_eV = None, # [eV], dim(nrho,), electron temperature
    ne_cm3 = None, # [cm^-3], dim(nrho,), electron density
    # Controls
    wvl_lim = None,
    ):

    # Obtains file location of PECs
    if path_PEC is None:
        path_PEC = _get_PEC_loc(
            ion=sp,
            cs=cs,
            )

    # Reads ADF15
    trs = aurora.read_adf15(path_PEC)

    # Finds every unique wavelength
    wvls = np.unique(trs['lambda [A]'])
    if wvl_lim is not None:
        wvls = wvls[(wvls >= wvl_lim[0]) & (wvls <= wvl_lim[1])]

    # Init
    out = np.zeros((
        len(wvls), len(Te_eV), 3
        )) # dim(ntrans, nrho, ion/exc/rec)

    # Loop over PECs
    for ii, wvl in enumerate(wvls):
        # Finds all transitions at this wavelength
        sub_trs = trs.loc[trs['lambda [A]'] == wvl]

        for ll in sub_trs['isel']:
            if sub_trs.loc[ll - 1]['type'] in ['ion', 'ioniz', 'ionis']:
                xx = 0
            elif sub_trs.loc[ll - 1]['type'] in ['exc', 'excit']:
                xx = 1
            elif sub_trs.loc[ll - 1]['type'] in ['rec', 'recom', 'drsat']:
                xx = 2

            out[ii,:,xx] += 10**sub_trs.loc[ll - 1]['log10 PEC fun'].ev(np.log10(ne_cm3), np.log10(Te_eV))

    # Output, [cm^3/s], dim(ntrans, nrhop, ion/exc/rec)
    return out, wvls

'''

import sys, os
from transport_world.run_AURORA import calc_spectra as csp

# Plasma
Te_long = np.logspace(
    np.log10(1e-1),
    np.log10(4e0),
    51
    )*1e3 # [eV]
ne_long = 1e14*np.ones_like(Te_long) # [cm^-3]

# File management
fol = '/home/cjperks/work/2008_SPARCxray/Alt_atom_data/HIREXSR'

'''

# Get PEC curve
def get_PEC(
    # Ion data
    sp = None,
    cs = None,
    path_PEC = None,
    # Plasma data
    Te_eV = None, # [eV], dim(nrho,), electron temperature
    ne_cm3 = None, # [cm^-3], dim(nrho,), electron density
    # Line data
    lamb0 = None,   # [AA], central wavelength
    dlamb = None,   # [AA], interval to search over
    ):

    # Obtains file location of PECs
    if path_PEC is None:
        path_PEC = _get_PEC_loc(
            ion=sp,
            cs=cs,
            )

    # Reads ADF15
    trs = aurora.read_adf15(path_PEC)

    # Line of interest
    if dlamb is None:
        dlamb = 0.01*lamb0
    inds = np.where(
        (trs['lambda [A]'] >= lamb0-dlamb)
        & (trs['lambda [A]'] <= lamb0+dlamb)
        )[0]

    out = np.zeros((len(Te_eV), 3)) # dim(nrho, ion/exc/rec)

    # Loop over PECs
    for ind in inds:
        if trs['type'][ind] in ['ion', 'ioniz', 'ionis']:
            xx = 0
        elif trs['type'][ind] in ['exc', 'excit']:
            xx = 1
        elif trs['type'][ind] in ['rec', 'recom', 'drsat']:
            xx = 2

        out[:,xx] += 10**trs['log10 PEC fun'][ind].ev(np.log10(ne_cm3), np.log10(Te_eV))

    # Output, [cm^3/s], dim(nrhop, ion/exc/rec)
    return out

# Full spectrum forward modeling
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
#               Calculation
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
        # Sums of ion, exc, rad rec, diel rec, and cx spectra

        # Forward Model total spectrum interpolated on master wavelength domain
        spec_fm[rho,:] = _remesh(
            wave_fm = out[0],
            spec_fm = spec_tmp,
            wave_mstr = wave_master,
            ) # [1/cm^3/s/AA]; dim(,lambda)

    # Outputs spectra for this species/charge state
    return spec_fm


#########################################################
#
#               Utilities
#
#########################################################


# Remesh simulation energy grid onto global grid
# NOTE: benchmarked to conserve photons upon energy-integration
def _remesh(
    wave_fm=None,       # Model mesh
    spec_fm=None,       # Model output
    wave_mstr=None,     # Master mesh
    ):

    spec_mstr = np.zeros(len(wave_mstr)) # [1/cm^3/s/AA], dim(lambda,)

    # Loop over wavelength
    for yy in np.arange(len(wave_mstr)):
        # Calculates the edges of the wavelength mesh
        if yy == 0:
            y_max = wave_mstr[yy]
            y_min = (wave_mstr[yy] + wave_mstr[yy+1])/2
        elif yy == len(wave_mstr)-1:
            y_max = (wave_mstr[yy-1] + wave_mstr[yy])/2
            y_min = wave_mstr[yy]
        else:
            y_max = (wave_mstr[yy-1] + wave_mstr[yy])/2
            y_min = (wave_mstr[yy] + wave_mstr[yy+1])/2

        # Finds data within new bin
        ind = np.where(
            (wave_fm >= y_min)
            & (wave_fm <= y_max) 
            )[0] # dim(wave_fm,)

        # Interpolates data onto new bin
        if ind.size == 0:
            # If extrpolating, fill with zeros
            if wave_mstr[yy] > np.max(wave_fm) or wave_mstr[yy] < np.min(wave_fm):
                spec_mstr[yy] = 0.0

            # If interpolating between two points
            else:
                ind1 = np.argmin(abs(wave_fm-wave_mstr[yy]))

                if wave_mstr[yy] - wave_fm[ind1] > 0:
                    spec_mstr[yy] = (
                        spec_fm[ind1]
                        + (wave_mstr[yy] - wave_fm[ind1])
                        * (spec_fm[ind1+1] - spec_fm[ind1])
                        / (wave_fm[ind1+1] - wave_fm[ind1])
                        )
                else:
                    spec_mstr[yy] = (
                        spec_fm[ind1]
                        + (wave_mstr[yy] - wave_fm[ind1])
                        * (spec_fm[ind1-1] - spec_fm[ind1])
                        / (wave_fm[ind1-1] - wave_fm[ind1])
                        )

        
        else:
            # Handles edge points of bin
            if ind[0] == 0:
                eps_l = 0.0
            else:
                m_l = (
                    (spec_fm[ind[0]-1] - spec_fm[ind[0]])
                    /(wave_fm[ind[0]-1] - wave_fm[ind[0]])
                    )
                b_l = (
                    spec_fm[ind[0]] - m_l*wave_fm[ind[0]]
                    )
                eps_l = m_l*y_min + b_l

            if ind[-1] == len(wave_fm)-1:
                eps_r = 0.0
            else:
                m_r = (
                    (spec_fm[ind[-1]] - spec_fm[ind[-1]+1])
                    /(wave_fm[ind[-1]] - wave_fm[ind[-1]+1])
                    )
                b_r = (
                    spec_fm[ind[-1]] - m_r*wave_fm[ind[-1]]
                    )
                eps_r = m_r*y_max + b_r

            # Rebins the data while conserving photons
            eps = np.append(spec_fm[ind], eps_r)
            eps = np.append(eps_l, eps)

            lamb = np.append(wave_fm[ind], y_max)
            lamb = np.append(y_min, lamb)

            spec_mstr[yy] = np.trapz(eps,lamb)/(y_max-y_min)

    # Output, [1/cm^3/s/AA], dim(lambda,)
    return spec_mstr


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