'''

Example script to turn radial profile output
into SPIRAL input files

cjperks
Nov 27, 2024

'''

# Module
import sys, os

from transport_world.run_SPIRAL import input_utils as iu

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# File management
#shot = 1140221012
outpath = os.path.join(
    '/home/cjperks/work',
    #'2201_Pumpout/CMOD',
    #'shots/%i/SPIRAL'%(shot),
    #'test_toric.cdf'
    '2201_Pumpout/SPARC',
    'runs/PRD',
    'T0_6keV_W47ss62/input_SPIRAL',
    )

fgacode = os.path.join(
    '/home/cjperks',
    'tofu_sparc/background_plasma',
    'PRD_plasma/run1',
    'input.gacode'
    )

# Writes profiles
iu.p2s(
    outpath = outpath,
    fgacode = fgacode,
    bulk_ion = 'D_ga',
    fast_ion = 'He3_ga',
    )