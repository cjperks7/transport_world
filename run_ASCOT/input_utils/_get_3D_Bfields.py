'''

Script to use 3D vacuum B-field data recieved from
Steve Wolfe

NOTES: 
    1) C-Mod consistented of 20 copper coils alternating
        in a "short"/"tall" pattern induced an n =10, 20
        toroidal ripple
    2) phi =0 is mid-plane of a "tall" TF coil
    3) Steve's simulation extends in phi between two 
        "short" TF coils phi = [-18,18] degrees

cjperks
OCt 3, 2024

'''

# Module
import os, sys
from itertools import islice
import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

plt.rcParams.update({'font.size': 16})

# Path to data file
filename = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD',
    'ripple',
    'ripple_cmod.dat'
    )

# Rows with data to skip rest
first_row = 270
last_row = 8369

# Column quantities
quants = ['R', 'phi', 'Z', 'Br', 'Bphi', 'Bz']

# Init
ddata = {}
ddata['raw'] = {}

for qq in quants:
    ddata['raw'][qq] = []

# Opens data file
with open(filename, 'r') as file:
    # Reads data rows skipping unneccassary parts
    for line in islice(file, first_row -1, last_row):
        data = line.split()

        if len(data) <= 1:
            continue

        # Grabs data
        for ii, qq in enumerate(quants):
            ddata['raw'][qq].append(float(data[ii]))

# Formats data in arrays
for qq in quants:
    ddata['raw'][qq] = np.asarray(ddata['raw'][qq])

# 1D vectors for volume
RR_1d = np.unique(ddata['raw']['R']) # dim(nR,); [m]
ZZ_1d = np.unique(ddata['raw']['Z']) # dim(nR,); [m]
phi_1d = np.unique(ddata['raw']['phi']) # dim(nphi,); [deg]






ddata['grid'] = {}
ddata['grid']['R'], ddata['grid']['Z'], ddata['grid']['phi'] = np.meshgrid(
    RR_1d, ZZ_1d, phi_1d,
    indexing = 'ij'
    )

# Reshapes the B-field data into its components
ddata['Bfield'] = {}

bbs = ['Bphi', 'Bz', 'Br']
for bb in bbs:
    ddata['Bfield'][bb] = np.zeros((
        len(RR_1d), len(ZZ_1d), len(phi_1d)
        ))

    cnt = 0
    for pp in np.arange(len(phi_1d)):
        for zz in np.arange(len(ZZ_1d)):
            for rr in np.arange(len(RR_1d)):
                ddata['Bfield'][bb][rr,zz,pp] = ddata['raw'][bb][cnt]
                cnt += 1


# Extends the space to Z>0
ddata['grid']['Z'] = np.concatenate((
    ddata['grid']['Z'],
    np.flip(
        ddata['grid']['Z'][:,:-1,:] *-1,
        axis = 1
        )
        ),
    axis = 1
    )
ddata['grid']['R'] = np.concatenate((
    ddata['grid']['R'],
    ddata['grid']['R'][:,:-1,:]
        ),
    axis = 1
    )
ddata['grid']['phi'] = np.concatenate((
    ddata['grid']['phi'],
    ddata['grid']['phi'][:,:-1,:]
        ),
    axis = 1
    )

ddata['Bfield']['Bphi'] = np.concatenate((
    ddata['Bfield']['Bphi'],
    np.flip(
        ddata['Bfield']['Bphi'][:,:-1,:],
        axis = 1
        )
        ),
    axis = 1
    )
ddata['Bfield']['Bz'] = np.concatenate((
    ddata['Bfield']['Bz'],
    np.flip(
        ddata['Bfield']['Bz'][:,:-1,:],
        axis = 1
        )
        ),
    axis = 1
    )
ddata['Bfield']['Br'] = np.concatenate((
    ddata['Bfield']['Br'],
    np.flip(
        ddata['Bfield']['Br'][:,:-1,:] *-1,
        axis = 1
        )
        ),
    axis = 1
    )




def _plot(ddata=None, gs=None, ind=None):
    # Loop over B-fields
    for ii, bb in enumerate(ddata['Bfield'].keys()):
        ax = fig.add_subplot(gs[ii])
        ax.clear()

        pcm = ax.pcolormesh(
            ddata['grid']['R'][..., ind],
            ddata['grid']['Z'][..., ind],
            ddata['Bfield'][bb][..., ind],
            )
        cbar = plt.colorbar(pcm, ax=ax)

        pcm.set_clim(
            np.min(ddata['Bfield'][bb].flatten()),
            np.max(ddata['Bfield'][bb].flatten())
            )

        ax.set_xlabel('R [m]')
        ax.set_title(bb +' [T]')
        ax.set_aspect('equal')

        if ii == 0:
            ax.set_ylabel('Z [m]')




phi_init = 0

fig = plt.figure(figsize=(18,8))
gs = gridspec.GridSpec(1, 3, figure=fig)


# Adjust layout to make space for the slider
plt.subplots_adjust(bottom=0.2)

# Create a slider axis
slider_ax = plt.axes([0.2, 0.1, 0.6, 0.03])  # [left, bottom, width, height]
phi_slider = Slider(
    slider_ax, 
    'phi [deg]', 
    -18, 18, 
    valinit=phi_init)

def update(val):
    ind = np.argmin(abs(
        ddata['grid']['phi'][0,0,:] - phi_slider.val
        ))

    _plot(ddata=ddata,gs=gs,ind=ind)

    fig.canvas.draw_idle()

phi_slider.on_changed(update)

# Initial plot
_plot(ddata=ddata, gs=gs, ind=phi_init)



#### ---- Plots Bt v. phi at fixed (R,Z) ---- ####

# Controls
Rtarg = 0.877 # [m]
Ztarg = 0.0 # [m]

indR = np.argmin(abs(
    ddata['grid']['R'][:,0,0] - Rtarg
    ))
indZ = np.argmin(abs(
    ddata['grid']['Z'][0,:,0] - Ztarg
    ))

pp = ddata['grid']['phi'][indR,indZ,:]
bb = ddata['Bfield']['Bphi'][indR,indZ,:]
dd = (max(bb)-min(bb))/(max(bb)+min(bb))

figa, axa = plt.subplots(figsize = (10,6))

axa.plot(
    ddata['grid']['phi'][indR,indZ,:],
    ddata['Bfield']['Bphi'][indR,indZ,:],
    'o-'
    )

axa.set_xlabel('Toroidal angle [deg]')
axa.set_ylabel(r'$B_t$ [T]')
axa.grid('on')
axa.set_title(
    r'(R,Z) = (%0.3f, %0.3f) m; $\delta$ = %0.3f%%'%(
        ddata['grid']['R'][indR,0,0], ddata['grid']['Z'][0,indZ,0],
        dd*1e2
        )
    )
