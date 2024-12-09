'''

Example script to turn TORIC E-field output
into SPIRAL input files

cjperks
Nov 6, 2024

'''

# Module
import sys, os

from transport_world.run_SPIRAL import input_utils as iu

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# File management
shot = 1140221012
toric_file = os.path.join(
    '/home/cjperks/work',
    #'2201_Pumpout/CMOD',
    #'shots/%i/TORIC'%(shot),
    #'toric_pnphi.ncdf'
    '2201_Pumpout/SPARC',
    'runs/PRD',
    'T0_6keV_W47ss62',
    'toric_pnphi.ncdf'
    )
outpath = os.path.join(
    '/home/cjperks/work',
    #'2201_Pumpout/CMOD',
    #'shots/%i/SPIRAL'%(shot),
    #'test_toric.cdf'
    '2201_Pumpout/SPARC',
    'runs/PRD',
    'T0_6keV_W47ss62/input_SPIRAL',
    'test_toric_2ant.cdf'
    )
#device = 'CMOD'
device = 'SPARC'

gfile = os.path.join(
    '/home/cjperks',
    'work/2201_Pumpout',
    'SPARC/runs',
    'PRD',
    'input.geq'
    )
#gfile = None
afile = None

# Writes SPIRAL input file
iu.t2s(
    toric_file = toric_file,
    outpath = outpath,
    #Prf_abs = 0.5*0.3,
    Prf_abs = 25,
    device = device,
    gfile = gfile,
    afile = afile,
    )