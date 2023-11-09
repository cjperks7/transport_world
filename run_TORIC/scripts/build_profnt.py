'''

Script to take input.gacode and impurity density profile
estimations and create a TORIC5 input file

'''

# Modules
from transport_world.run_TORIC import input_utils as uts
import os


# Paths
shot = '1140408020'
in_path = os.path.join(
    '/home/cjperks',
    '2201_Pumpout',
    'CMOD/shots',
    shot
    )
#fgacode = 'input_t0.8s.gacode'
fgacode = 'input_t0.6s.gacode'

#######################################################
#
#            Impurity Modeling
#
#######################################################
'''
from transport_world import run_AURORA as rA

dmodel = rA.def_dmodel(
    option='Ar_TGLF',
    shot=shot,
    fgacode=fgacode,
    )

dimp = rA.calc_imp_dens(
    dmodel = dmodel,
    plt_all = True,
    )
'''

from omfit_classes import omfit_gapy
from transport_world.make_gacode import exp_utils as exp

inputga = omfit_gapy.OMFITgacode(
    os.path.join(
        in_path,
        fgacode,
        )
    )

_, dimp = exp.get_imprad(
    fimprad = os.path.join(in_path, "imprad_"+shot+".npy"),
    rhop = np.sqrt(inputga['polflux']/inputga['polflux'][-1]),
    rescl = 0.08657
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
            #dimp['nz_ga'][:,16], # [cm^-3], dim(nrho,)
            #dimp['nz_ga'][:,17] # [cm^-3], dim(nrho,)
            dimp[:,16], # [cm^-3], dim(nrho,)
            dimp[:,17] # [cm^-3], dim(nrho,)
            ]
        }
    )
