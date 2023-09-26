'''

_calc_imp_DV is a function to calculate impurity D, V profiles
from TGYRO flux-matched profiles

'''

# Modules
from portals.gacode_tools 	import TGLFmodule
import numpy as np
from transport_world.plot_utils import plt_slider
from scipy.interpolate  import interp1d
import os

__all__ = [
    'calc_imp_turbDV',
    'cmp_tglf_spectra',
    ]

############################################
#
#               Main
#
############################################

def calc_imp_turbDV(
    # File management
    fgacode = None,
    folder = None,
    subfolder = None,
    # Impurity settings
    cs = 16,
    amu = 40,
    # Simulation settings
    rhos = list(np.linspace(0.25,0.85,31)), # sq. tor. flux
    restart = True,     # True = from scratch, False = already have output
    TGLFsettings = 4,   # 1 -> SAT1 old, 2 -> SAT0, 3 -> SAT1geo, 4 -> SAT2, 5 -> SAT2em
    # PLotting
    plt_all = False,
    plt_spec = False,
    plt_type = 'freq',
    plt_lmts = {
        0:{
            'xlim':None,
            'ylim':None,
            },
        1:{
            'xlim':None,
            'ylim':None,
            },
        }
    ):

    # Initialize TGLF class at the locations
    tglf = TGLFmodule.TGLF(rhos=rhos)

    # Prepare input files to TGLF (will run dummy iteration of TGYRO to populate files for TGLF)
    cdf = tglf.prep(
        folder,
        restart=restart,
        inputgacode=fgacode
        )

    # Runs PORTALS analysis function to determine impurity D,V profiles
    tglf.runAnalysis(
        subFolderTGLF = subfolder,
        analysisType='Z',
        TGLFsettings = TGLFsettings,
        restart = restart,
        trace = [cs,amu],
        )

    # Plots growth rate spectra
    if plt_spec:
        _plot(
            tglfs=[tglf],
            plt_type=plt_type,
            plt_lmts=plt_lmts
            )
        # Plots the TGLF analysis results
        if plt_all:
            tglf.plotAnalysis(
                labels=['analysis1'],
                analysisType='Z',
                )

    return {
        'rhot': tglf.rhos, # sq. norm. tor. flux, dim(tglf_rhot,)
        'DZ': tglf.scans['analysis1']['DZ'], # [m^2/s], dim(tglf_rhot,)
        'VZ': tglf.scans['analysis1']['VZ'], # [m/s], dim(tglf_rhot,)
        'tglf': tglf
        }

############################################
#
#               Compare TGLF spectra
#
############################################

def cmp_tglf_spectra(
    in_path=None,
    fgacodes=None,
    name_sims=None,
    folder=None,
    # Plotting
    plt_type = 'freq',
    plt_lmts = {
        0:{
            'xlim':None,
            'ylim':None,
            },
        1:{
            'xlim':None,
            'ylim':None,
            },
        }
    ):

    # Initializes
    tglfs = []

    # Loop over simulations
    for ii, name_sim in enumerate(name_sims):
        dout = calc_imp_turbDV(
            fgacode = os.path.join(
                in_path,
                fgacodes[ii]
                ),
            folder = os.path.join(
                in_path,
                'output_TGLF',
                folder,
                name_sim
                ),
            subfolder = name_sim,
            restart=False,
            )

        tglfs.append(dout['tglf'])

    # Plotting
    _plot(
        tglfs=tglfs,
        plt_type=plt_type,
        plt_lmts=plt_lmts
        )



############################################
#
#               Plotting
#
############################################

def _plot(
    tglfs=None,
    plt_type = 'freq',
    plt_lmts = None,
    ):

    # Initializes data to plot
    freq_mstr = np.zeros((len(tglfs), 31, 28)) # dim(ntglf, nrho, nky)
    gamma_mstr = np.zeros((len(tglfs), 31, 28)) # dim(ntglf, nrho, nky)

    # Unique labels
    dzlabels = []
    ky_mstr = []

    # Loop over TGLF results
    for tglf in tglfs:
        # Key to results
        keyr = list(tglf.results.keys())[1]

        dzlabels.append(
            keyr.split('_')[1]
            + '_'
            + keyr.split('_')[3].split('R')[0]
            )

        # Obtains wavenumber mesh
        ky = tglf.results[keyr]['ky'] # dim(nrho,nky)

        # Interpolates onto uniform ky mesh
        ky_mstr.append(np.mean(ky, axis=0)) # dim(nky,)

    ky_mstr = np.mean(ky_mstr, axis=0) # dim(nky,)

    # Loop over TGLF results
    for jj, tglf in enumerate(tglfs):

        # Key to reuslts
        keyr = list(tglf.results.keys())[1]

        # Obtains the radial mesh
        rhot = tglf.results[keyr]['x'] # dim(nrho,)

        # Obtains wavenumber mesh
        ky = tglf.results[keyr]['ky'] # dim(nrho,nky)

        # Obtains frequency and growth rate spectra
        freq = tglf.results[keyr]['freq'] # dim(nrho,nky)
        gamma = tglf.results[keyr]['gamma'] # dim(nrho, nky)

        # Loop over nrho
        for ii in np.arange(freq.shape[0]):

            freq_mstr[jj,ii,:] = interp1d(
                ky[ii,:],
                freq[ii,:],
                bounds_error = False,
                fill_value = (freq[ii,0], freq[ii,-1])
                )(ky_mstr) # dim(nrho, nky)

            gamma_mstr[jj,ii,:] = interp1d(
                ky[ii,:],
                gamma[ii,:],
                bounds_error = False,
                fill_value = (gamma[ii,0], gamma[ii,-1])
            )   (ky_mstr) # dim(nrho, nky)

    if plt_type == 'freq':
        data1 = freq_mstr
        data2 = gamma_mstr
        lab1 = r"$\omega$ [$C_s/a$]"
        lab2 = r"$\gamma$ [$C_s/a$]"
    elif plt_type == 'mixing':
        data1 = freq_mstr/ky_mstr[None,None,:]
        data2 = gamma_mstr/ky_mstr[None,None,:]
        lab1 = r"$\omega/(k_\theta\rho_s)$ [$C_s/a$]"
        lab2 = r"$\gamma/(k_\theta/\rho_s)$ [$C_s/a$]"

    # Plots slider
    plt_slider(
        xxx=ky_mstr,
        yyy=rhot,
        dzzz=data1,
        xxlabel=r"$k_\theta\rho_s$",
        yylabel=r"$\rho_t$",
        zzlabel= lab1,
        dzlabels=dzlabels,
        plt_sum=False,
        y_line = 0.0,
        xscale='log',
        yscale='symlog',
        xlim=plt_lmts['omega']['xlim'],
        ylim=plt_lmts['omega']['ylim'],
        )

    plt_slider(
        xxx=ky_mstr,
        yyy=rhot,
        dzzz=data2,
        xxlabel=r"$k_\theta\rho_s$",
        yylabel=r"$\rho_t$",
        zzlabel=lab2,
        dzlabels=dzlabels,
        plt_sum=False,
        xscale='log',
        yscale='log',
        xlim=plt_lmts['gamma']['xlim'],
        ylim=plt_lmts['gamma']['ylim'],
        )


