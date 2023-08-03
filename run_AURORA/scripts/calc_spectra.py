'''

calc_spectra is an example script to use AURORA to
calculate photon radiation spectra

'''

from transport_world import run_AURORA as rA


# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2


#################################################
#
#           User-controls
#
#################################################

path = (
    '/home/cjperks/2201_Pumpout/CMOD/shots/'
    + '1140221013/'
    )
ga = 'input_t1.2s.gacode'
geq = 'g1140221013.01000'
fgacode = path+ga

# Diagnostic modeling controls
ddiag = {
    'spec_rng': {   # Spectral range of diag
        'units': 'AA',
        'rng': [3.94, 4.00],    # Wavelength mesh bounds
        'nlamb': 1e4,           # Wavelength mesh resolution   
    },
    'imps':{        # Impurities imaged by diag
        'Ar': [16], # Charge-states
    },
}

# Impurity density modeling controls
dmodel = {
    #'options': ['TGLF', 'FACIT'],
    'options': ['Flat', 'FACIT'],
    'imp':{
        'sym': 'Ar',
        'cs': 16,
        'amu': 40
        },
    'paths':{
        'fgacode': ga,
        'geqdsk': geq,
        'folder': path,
        'subfolder': 'xxx',
        },
    'Flat': {
        'D': 1e4, # [cm2/s]
        'V': -2e2 # [cm/s]
        },
    'TGLF': {
        'restart': True,
        'TGLFsettings': 4,
        },
    'FACIT': {
        'rotation_model': 2,
        'ICRH_option': None,
        },
    'AURORA':{
        'device': 'CMOD',
        'main_ion_A': 2,
        'times': np.array([0., 0.7]), # [s]
        'source': {
            'type': 'const',
            'rate': 1e20 # [particles/s]
            },
        'edge': {
            'source_cm_out_lcfs': 1.0,
            'recycling_flag': True,
            'wall_recycling': 1.0,
            },
        'target':{
            'c_imp': 1e-3,
            'P_rad': np.nan
            },
        },
    }


#################################################
#
#           Calculation
#
#################################################

# Obtains impurity density profile
dout = rA.calc_imp_dens(
    dmodel = dmodel,
    plt_all = False,
    )

# Obtains line radiation synthetic spectra
dspectra = rA.get_line_spectra(
    ddiag=ddiag,
    fgacode=fgacode,
    dnz={
        'Ar':dout['nz_ga'],
        },
    imp_res = True,
    plt_all = True,
    )