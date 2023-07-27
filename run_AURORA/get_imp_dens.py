'''

get_imp_dens is an example script to 
obtain impurity density profiles

'''

# Modules
from transport_world import run_AURORA as rA

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

#######################################################
#
#            User controls
#
#######################################################

dmodel = {
    'option': ['TGLF', 'FACIT'],
    #'option': ['Flat'],
    'imp':{
        'sym': 'Ar',
        'cs': 16,
        'amu': 40
        },
    'paths':{
        'fgacode': 'input.gacode.new',
        'geqdsk': 'input.geq',
        'folder': 'xxx/',
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
        },
    'AURORA':{
        'device': 'CMOD',
        'main_ion_A': 2,
        'times': np.array([0., 0.3]), # [s]
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


#######################################################
#
#            Calculation
#
#######################################################

nZ = rA.calc_imp_dens(
    dmodel = dmodel,
    plt_all = True,
    )