'''

Script to read itm_waves.ascii files

cjperks
July 26, 2024

'''

# Modules
from transport_world.run_ASCOT.input_utils import _read_waves as rw
import os

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

file = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/ASCOT_template',
    'aug_ex.ascii'
    )

dout = rw.read_waves(
    file=file
    )
