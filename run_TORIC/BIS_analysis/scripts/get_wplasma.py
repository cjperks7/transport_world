'''

Script to obtain Wplasma from a local EFIT run

cjperks
Aug 7th, 2024

'''

# Modules
import os
import MDSplus

# Shot of interest
shot = 1140221013

# Sets environment to recognize local EFIT trees
path2efit = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD',
    'shots',
    str(shot),
    'BIS_analysis/trees'
    )
os.environ['efit_cjp_path'] = path2efit

# Opens local EFIT
efit = MDSplus.Tree('efit_cjp', shot)

# Gets stored energy
wplasma = efit.RESULTS.A_EQDSK.WPLASM.data() # [kJ]
t_wplasma = efit.RESULTS.A_EQDSK.WPLASM.dim_of().data() # [s]

