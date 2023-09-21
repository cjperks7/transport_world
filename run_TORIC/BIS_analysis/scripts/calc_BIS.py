'''
Runs break-in-slope analysis

'''

# Modules
from transport_world.run_TORIC import BIS_analysis as BIS

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

ddata = BIS.get_BIS(
    shot = 1140221013,
    t1 = [0.75,1.0],
    t2 = [1.0,1.5],
    plt_all = True,
    )