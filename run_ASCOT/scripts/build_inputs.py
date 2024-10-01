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
    'work/2201_Pumpout',
    'CMOD/shots',
    shot
    )
fgacode = 'input_t0.8s.gacode'
fgfile = 'g'+shot+'.01000'

# Species
ions = ['D', 'H', 'Ar16']
zn = [1, 1, 18]
zion = [1, 1, 16]
amn = [2, 1, 40]

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
    charges = zion,
    masses = amn,
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

# Locate X-point on plot
from transport_world.run_ASCOT import input_utils as utz
in_path = ''
fgfile = 'g1140221013.01000'
Xpt_RZ = [0.55334, -0.3939]
#utz.find_Xpt(
#    in_path=in_path,
#    fgfile=fgfile,
#    Xpt_RZ = Xpt_RZ,
#    )

# Writes input data
utz.write_magn(
    in_path=in_path,
    fgfile=fgfile,
    file_type = 'eqdsk',
    B_type = '2d',
    Xpt_RZ = Xpt_RZ,
    )

# 3) 2D wall geometry
utz.write_wall2d(
    # Files
    in_path=in_path,
    fgfile= fgfile,
    )


# 4) Writes RFOF input files
zion.append(9)
zn.append(9)
amn.append(19)
utz.edit_waves(
    in_path=in_path,
    RF_abs = 4.0e3, # [W]
    zion = zion,
    zn = zn, 
    amn = amn,
    nspec = 3, # matlab indexing
    )

# 5) Write rho_tor, rho_pol conversion
utz.write_rho(
    # Files
    in_path=in_path,
    fgacode=fgacode
    )