'''

Script to take input.gacode and impurity density profile
estimations and create a TORIC5 input file

'''

# Modules
from transport_world.run_TORIC import input_utils as uts
import os
from transport_world import run_AURORA as rA

# Paths
shot = '1140221013'
in_path = os.path.join(
    '/home/cjperks',
    '2201_Pumpout',
    'CMOD/shots',
    shot
    )
fgacode = 'input_t0.8s.gacode'

#######################################################
#
#            Impurity Modeling
#
#######################################################

dmodel = {
    #'options': ['TGLF', 'FACIT'],
    #'options': ['Flat', 'FACIT'],
    'options': ['TGLF'],
    #'options': ['FACIT'],
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
        'geqdsk': 'g1140221013.01000',
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
            'c_imp': 1e-3,
            'P_rad': np.nan,
            #'BC': {'rhop': 0.98, 'val': 1.0},
            },
        },
    }

dimp = rA.calc_imp_dens(
    dmodel = dmodel,
    plt_all = False,
    )


#######################################################
#
#            File writing
#
#######################################################

# Writes input file
uts.build_profnt(
    # Files
    in_path=in_path,
    fgacode=fgacode,
    # Species
    ions = ['D', 'H', 'Ar16', 'Ar17'],
    charges = [1, 1, 16, 17],
    masses = [2, 1, 40, 40],
    # Ion temperature controls
    diff_temp = 0,
    # Optional dictioanry to overwrite nz profile with external
    ext_nz = {
        'ion': ['Ar16', 'Ar17'],
        'data': [
            dimp['nz_ga'][:,16], # [cm^-3], dim(nrho,)
            dimp['nz_ga'][:,17] # [cm^-3], dim(nrho,)
            ]
        }
    )
