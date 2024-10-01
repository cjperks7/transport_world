'''

Script to convert magnetic equilibirum data (e.g. eqdsk)
into ASCOT input files

cjperks
Sep 10, 2024

NOTE: The sign conventions assumed in ASCOT are
1) Looking down on the tokamak, CCW is positive (Bt, Ip)
2) In the poloidal plane, up in the outboard midplane is postive
    So a negative current would have a positive Bp, Psi (per RHR)

NOTE: All the C-Mod shots I'm interested in for the pumpout
project have CCW Bt and Ip

'''

# Modules
import numpy as np
from omfit_classes import omfit_eqdsk
import os
import matplotlib.pyplot as plt

__all__ = [
    'find_Xpt',
    'write_magn'
]

################################################
#
#           Main
#
################################################

# Makes plot to help find X-point (R,Z) coordinate
def find_Xpt(
    fgfile = None,
    Xpt_RZ = [0.55334, -0.3939],
    ):

    # Reads data
    ddata = _get_eqdsk(fgfile)
    
    # Init
    rhop = np.sqrt(ddata['AuxQuantities']['PSIRZ_NORM']) # dim(R, Z)
    R_1d = ddata['AuxQuantities']['R']
    Z_1d = ddata['AuxQuantities']['Z']

    # Plots
    fig, ax = plt.subplots()

    ax.contour(
        R_1d, Z_1d, rhop,
        levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        )

    ax.plot(
        Xpt_RZ[0], Xpt_RZ[1], 'r*'
        )

# Manages writing input files
def write_magn(
    ascot_path = None,
    fgfile = None,
    file_type = 'eqdsk',
    B_type = '2d',
    Xpt_RZ = [0.55334, -0.3939]
    ):

    # Loads data
    if file_type == 'eqdsk':
        ddata = _get_eqdsk(fgfile)
    else:
        print("NOT YET IMPLEMENTED!!!")

    # Writes input data
    if B_type == '2d':
        _write_bkg_2d(
            ddata = ddata,
            ascot_path =ascot_path,
            )
        _write_header_2d(
            ddata = ddata,
            ascot_path = ascot_path,
            Xpt_RZ = Xpt_RZ,
            )
    else:
        print("NOT YET IMPLEMENTED")


################################################
#
#           Writing
#
################################################

# Write B-field data header
def _write_header_2d(
    ddata = None,
    ascot_path = None,
    Xpt_RZ = None, # [R, Z]
    ):

    # Init
    f = open(
        os.path.join(
            ascot_path,
            'input.magn_header'
            ),
        'w'
        )

    # Shot details that I'll ignore
    f.write(
        '0'.rjust(8, ' ') # shot umber
        + ('%0.6f'%(0)).rjust(10, ' ') # time in shot
        + '0'.rjust(3, ' ') # modifier flag
        + '\n'
        )

    # Device name
    f.write('CMOD\n')

    # Something, plasma current
    f.write(
        '0'.rjust(4, ' ') # something unimportant? FPPkat
        + ('%0.6f'%(ddata['CURRENT'])).rjust(18, ' ') # [A]
        + '\n'
        )

    # Number of "special points" (magnetic axis, X-point)
    f.write(
        '2'.rjust(4, ' ')
        + '\n'
        )

    # Poloidal flux at special points
    tmp = 2*np.pi*(
        ddata['fluxSurfaces']['geo']['psi']
        -ddata['fluxSurfaces']['geo']['psi'][-1]
        ) # NOTE should have tmp <0 per helicity convention
    f.write(
        ('%0.6f'%(tmp[0])).rjust(10, ' ')
        + ('%0.6f'%(tmp[-1])).rjust(10, ' ')
        + '\n'
        )

    # Radius at special points
    f.write(
        ('%0.6f'%(ddata['RMAXIS'])).rjust(10, ' ')
        + ('%0.6f'%(Xpt_RZ[0])).rjust(10, ' ')
        + '\n'
        )

    # Height at special points
    f.write(
        ('%0.6f'%(ddata['ZMAXIS'])).rjust(10, ' ')
        + ('%0.6f'%(Xpt_RZ[1])).rjust(10, ' ')
        + '\n'
        )

    # SSQ?
    f.write(
        ('%0.6f'%(ddata['RCENTR'])).rjust(10, ' ')
        + ('%0.6f'%(ddata['ZMID'])).rjust(10, ' ')
        + ('%0.6f'%(
            ddata['fluxSurfaces']['geo']['a'][-1]
            )).rjust(10, ' ') # horizontal minor radius
        + ('%0.6f'%(
            ddata['fluxSurfaces']['geo']['a'][-1]
            * ddata['fluxSurfaces']['geo']['kap'][-1]
            )).rjust(10, ' ') # vertical minor radius
        + '\n'
        )

    # Len of 1d variables
    f.write(
        ('%i'%(len(ddata['fluxSurfaces']['geo']['psi']))).rjust(4, ' ')
        + '\n'
        )

    # rho vector
    _write_block(
        f=f,
        data=np.sqrt(ddata['fluxSurfaces']['geo']['psin'])[:,None]
        )

    # 1d poloidal flux vector
    _write_block(
        f=f,
        data=tmp[:,None]
        )

    # 1d volume vector
    _write_block(
        f=f,
        data=ddata['fluxSurfaces']['geo']['vol'][:,None]
        )

    # 1d area vector
    _write_block(
        f=f,
        data=ddata['fluxSurfaces']['geo']['cxArea'][:,None]
        )

    # q profile
    _write_block(
        f=f,
        data=ddata['QPSI'][:,None]
        )

    # Closes file
    f.close()

# Write background B-field data
def _write_bkg_2d(
    ddata = None,
    ascot_path = None,
    ):

    # Init
    f = open(
        os.path.join(
            ascot_path,
            'input.magn_bkg'
            ),
        'w'
        )

    # Writes the header data
    f.write(
        '0'.rjust(18, ' ') # origin toroidal coordinate
        + '0'.rjust(4, ' ') # number of sectors
        + '1'.rjust(4, ' ') # nphi per sector
        + '0'.rjust(4, ' ') # number of coils
        + '1'.rjust(4, ' ') # zero at coil
        + '\n'
        )

    # Writes (R,Z) mesh limits
    f.write(
        ('%0.10f'%(min(ddata['AuxQuantities']['R']))).rjust(18, ' ')
        + ('%0.10f'%(max(ddata['AuxQuantities']['R']))).rjust(18, ' ')
        + ('%i'%(len(ddata['AuxQuantities']['R']))).rjust(4, ' ')
        + '\n'
        )
    f.write(
        ('%0.10f'%(min(ddata['AuxQuantities']['Z']))).rjust(18, ' ')
        + ('%0.10f'%(max(ddata['AuxQuantities']['Z']))).rjust(18, ' ')
        + ('%i'%(len(ddata['AuxQuantities']['Z']))).rjust(4, ' ')
        + '\n'
        )

    # Relevant when nsector >0
    f.write('0\n0\n')

    # Write psi block
    _write_block(
        f=f,
        data=2*np.pi*(
            ddata['AuxQuantities']['PSIRZ']
            -ddata['fluxSurfaces']['geo']['psi'][-1]
            )
        )

    # Writes Br block
    _write_block(
        f=f,
        data=ddata['AuxQuantities']['Br']*-1 # -1 to account for helicity convention
        )

    # Writes Bphi block
    _write_block(
        f=f,
        data=ddata['AuxQuantities']['Bt']
        )

    # Writes Bz block
    _write_block(
        f=f,
        data=ddata['AuxQuantities']['Bz']*-1 # -1 to account for helicity convention
        )

    # Closes file
    f.close()
    

################################################
#
#           Utilities
#
################################################

# Write data block
def _write_block(
    f=None,
    data=None, # dim(n1,n2)
    ):

    # Number of columns
    ncol = 4

    # Number of rows
    nrow = int(np.ceil(
        data.shape[0]*data.shape[1]
        /ncol
        ))

    # Init
    line = ''
    dd = data.flatten()
    cnt = 0

    for ii in np.arange(nrow):
        for jj in np.arange(ncol):
            if cnt >= len(dd):
                continue
            line += (
                '%0.10f'%(
                    dd[cnt]
                    )
                ).rjust(18, ' ')
            cnt += 1

        line += '\n'

    f.write(line)

# If data type is an eqdsk
def _get_eqdsk(
    file = None,
    ):

    # Output
    return omfit_eqdsk.OMFITgeqdsk(file)
