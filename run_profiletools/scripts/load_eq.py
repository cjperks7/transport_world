'''

Example script to load gfile using eqtools

'''

# Modules
import sys, os
from transport_world.run_profiletools import eqtools3 as eq

# File managment
gfile = os.path.join(
    '/home/cjperks',
    'work/2201_Pumpout',
    'SPARC/runs',
    'PRD',
    'input.geq'
    )
afile = None
machine = 'SPARC'

dedr = eq._get_eq(
    gfile = gfile,
    afile = afile,
    machine = machine
    )