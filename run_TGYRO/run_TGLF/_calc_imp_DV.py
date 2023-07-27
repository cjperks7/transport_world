'''

_calc_imp_DV is a function to calculate impurity D, V profiles
from TGYRO flux-matched profiles

'''

# Modules
from portals.gacode_tools 	import TGLFmodule

__all__ = [
    'calc_imp_DV',
    ]

############################################
#
#               Main
#
############################################

def calc_imp_DV(
    # File management
    fgacode = None,
    folder = None,
    subfolder = None,
    # Impurity settings
    cs = 16,
    amu = 40,
    # Simulation settings
    rhos = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], # sq. tor. flux
    restart = True,     # True = from scratch, False = already have output
    TGLFsettings = 4,   # 1 -> SAT0, 2 -> SAT1, 3 -> SAT1geo, 4 -> SAT2, 5 -> SAT2em
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

    # Plots the TGLF analysis results
    tglf.plotAnalysis(
        labels=['analysis1'],
        analysisType='Z',
    )

    return {
        'rhot': tglf.rhos, # sq. norm. tor. flux, dim(tglf_rhot,)
        'DZ': tglf.scans['analysis1']['DZ'], # [m^2/s], dim(tglf_rhot,)
        'VZ': tglf.scans['analysis1']['VZ'] # [m/s], dim(tglf_rhot,)
        }
