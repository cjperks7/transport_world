'''

build_inputs.py is an example script to write
input files for ASCOT4-RFOF

cjperks
Oct. 2nd 2023

'''

# Modules
from transport_world.run_ASCOT import input_utils as utz
from transport_world import run_AURORA as rA
import os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2


# Paths
shot = '1140221013'
in_path = os.path.join(
    '/home/cjperks',
    '2201_Pumpout',
    'CMOD/shots',
    shot
    )
fgacode = 'input_t0.8s.gacode'
fgfile = 'g'+shot+'.01000'

#######################################################
#
#            Impurity Modeling
#
#######################################################

dmodel = rA.def_dmodel(
    option='Ar_TGLF',
    shot=shot,
    fgacode=fgacode,
    )

dimp = rA.calc_imp_dens(
    dmodel = dmodel,
    plt_all = False,
    )

#######################################################
#
#            Writes ASCOT input files
#
#######################################################

# Writes standard ASCOT4 input files

# 1) Plasma n/T profiles
utz.write_plasma1d(
    # Files
    in_path=in_path,
    fgacode=fgacode,
    # Species
    ions = ['D', 'H', 'Ar16'],
    charges = [1, 1, 16],
    masses = [2, 1, 40],
    # Optional dictioanry to overwrite nz profile with external
    ext_nz = {
        'ion': ['Ar16', 'Ar17'],
        'data': [
            dimp['nz_ga'][:,16], # [cm^-3], dim(nrho,)
            dimp['nz_ga'][:,17] # [cm^-3], dim(nrho,)
            ]
        }
    )

# 2) Magnetic background
#### Use MATLAB scripts 

# 3) 2D wall geometry
utz.write_wall2d(
    # Files
    in_path=in_path,
    fgfile= fgfile,
    )


# Writes RFOF input files
utz.edit_waves(
    in_path=in_path,
    RF_abs = 4.0e3, # [W]
    )
