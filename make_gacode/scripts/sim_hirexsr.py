'''

Script to overpolot HIREXSR spectra against
Aurora simulations

cjperks

8/7/23

'''

from transport_world.make_gacode import read_tokamak as rTok

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# User-controls
shot = 1140221013
tht = 0

# Obtain HIREXSR data
ddata = rTok.get_hirexsr(
    shot=shot,
    tht=tht,
    quants = ['int', 'spectra', 'moments'],
    plt_all = True,
    plt_ch = 7,
    )