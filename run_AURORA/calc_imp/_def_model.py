'''
Stores default impurity modeling options

cjperks
Oct 2nd, 2023
'''

# Modules 
import numpy as np
import os


__all__ = [
    'def_dmodel',
    ]

################################################
#
# Main
#
#################################################

def def_dmodel(
    option=None,
    shot=None,
    fgacode=None,
    ):

    # Default models
    dmodel = {
        'Ar_TGLF': {
            'options': ['TGLF'],
            'imp':{
                'sym': 'Ar',
                'cs': 16,
                'amu': 40
                },
            'paths':{
                'in_path': os.path.join(
                    '/home/cjperks/'
                    '2201_Pumpout/CMOD/shots/',
                    shot
                    ),
                'fgacode': fgacode,
                'geqdsk': 'g'+shot+'.01000',
                'out_folder': os.path.join(
                    'output_TGLF',
                    'data_092623',
                    ),
                'name_sim': 'run_t0.8s_nrho31_SAT2',
                },
            'Flat': {
                'D': 1e4, # [cm2/s]
                'V': -5e2 # [cm/s]
                },
            'TGLF': {
                'restart': False,
                'TGLFsettings': 1,
                'rho_accept': 0.43,
                },
            'FACIT': {
                'rotation_model': 2,
                'ICRH_option': None,
                },
            'AURORA':{
                'device': 'CMOD',
                'main_ion_A': 2,
                'times': np.array([0., 1.0]), # [s]
                'source': {
                    'type': 'const',
                    'rate': 1e20 # [particles/s]
                    },
                'SOL': {
                    'decay': 1.0 # [cm]
                    },
                'edge': {
                    'source_cm_out_lcfs': 0.0,
                    'recycling_flag': True,
                    'wall_recycling': 1.0,
                    },
                'target':{
                    'c_imp': 0.5e-3,
                    'P_rad': np.nan,
                    #'BC': {'rhop': 0.98, 'val': 1.0},
                    },
                'acd': None,
                'scd': None,
                },
            },
        }

    # Output
    return dmodel[option]