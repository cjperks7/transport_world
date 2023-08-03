'''

scope_impDV.py is an example script for running TGLF
to calculate trace impurity D,V profiles

'''

# Modules
from transport_world.run_TGYRO import run_TGLF as rT
import os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2


#################################################
#
#               User-controls
#
#################################################

# File management
in_path = os.path.join(
    '/home/cjperks/',
    'tofu_sparc/background_plasma',
    'PRD_plasma',
    'run1',
    )

fgacode = os.path.join(
    in_path,
    'input.gacode'
    )

folder = os.path.join(
    in_path,
    'output_TGYRO',
    'test_TGLF_080223'
    )

name_sim = 'run1'

# Trace impurity
#cs = 16
#amu = 40
cs = 34
amu = 83.8

# Simulation settings
rhos = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

restart = True

TGLFsettings = 4

# Plotting
plt_all = True



#################################################
#
#               Running
#
#################################################

# Output dictionary
dout = rT.calc_imp_turbDV(
    # File management
    fgacode=fgacode,
    folder=folder,
    subfolder=name_sim,
    # Trace impurity
    cs=cs,
    amu=amu,
    # Simulation settings
    rhos=rhos,
    restart=restart,
    TGLFsettings=TGLFsettings,
    # Plotting
    plt_all=plt_all,
    )

