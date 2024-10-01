'''

_get_waves.py writes on RFOF E-field input file

cjperks
Oct 2nd, 2023

'''

# Modules
import numpy as np
import os
from transport_world.run_TORIC.output_utils import toric_tools
from datetime import date
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt

__all__ = [
    'edit_waves',
    #'write_waves'
    ]


########################################
#
#              Main
#
########################################

# Makes RFOF input file from skeleton
def edit_waves(
    toric_path = None,
    ascot_path = None,
    toric_ncdf = ['toric_test.ncdf'],
    device = 'CMOD',
    RF_abs = [4.0e3], # [W]
    nphi = [int(10)],
    zion = [16],
    zn = [18],
    amn = [40],
    Rmaj = 68.496/100, # [m]
    nspec = 3, # matlab indexing
    ):

    # Path to skeleton file
    infile = '/home/cjperks/work/2201_Pumpout/ASCOT_template/skeleton_waves_v2.ascii'

    # Path to write RFOF input file
    if ascot_path is None:
        ascot_path = toric_path
    outfile = os.path.join(
        ascot_path,
        'itm_waves.ascii'
        )

    # Reads TORIC output
    #toric = toric_tools.toric_analysis(
    #    toric_name='/home/cjperks/2201_Pumpout/CMOD/shots/1140221013/toric_test.ncdf',
    #    mode='ICRF',
    #    path='/home/cjperks/2201_Pumpout/CMOD/shots/1140221013/'
    #    )
    toric = {}
    for ii in np.arange(len(RF_abs)):
        toric[ii] = toric_tools.toric_analysis(
            toric_name=os.path.join(toric_path, toric_ncdf[ii]),
            mode='ICRF',
            path=toric_path,
            )

    # Antenna data
    if device == 'CMOD':
        dant = {
            'freq': int(80e6),
            'ntor': nphi,
            'Rmaj': Rmaj
            }
    elif device == 'SPARC':
        dant = {
            'freq': int(120e6),
            'ntor': nphi,
            'Rmaj': Rmaj
        }
    dspec = {
        'zion': zion,
        'zn': zn,
        'amn': amn,
        'nspec': nspec,
        }

    # Copies skeleton file
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        # Loop over skeleton file
        prev_line='x%y'
        for line in f_in:

            # Edits skeleton file if special line
            output_line, prev_line = _edit_line(
                line=line,
                prev_line=prev_line,
                toric=toric,
                dant=dant,
                dspec=dspec,
                RF_abs=RF_abs, # [W]
                )

            # Writes
            f_out.write(output_line)

    # Closes files
    f_in.close()
    f_out.close()

'''
num = 0
with open('/home/cjperks/2201_Pumpout/CMOD/shots/1140221013/input_ASCOT/itm_waves.ascii', 'r') as f:
    for line in f:
        print(num)
        print(line)
        num+=1
'''
########################################
#
#              Utilities
#
########################################

# Edits special lines
def _edit_line(
    line=None,
    prev_line=None,
    toric=None,
    dant=None,
    dspec=None,
    RF_abs=None, # [W]
    ):

    # List of arrays to edit
    to_edit = [
        #### --- Admin data --- ####
        ['datainfo', 'dataprovider\n'],
        ['datainfo', 'putdate\n'],
        ['datainfo', 'source\n'],
        #### --- Defining wave --- ###
        [' waves', 'coherentwave\n'],
        ['type','id\n'],
        ['type','flag\n'],
        ### --- Defining species --- ###
        ['compositions', 'nuclei\n'],
        ['compositions', 'ions\n'],
        ### --- Defining antenna --- ###
        ['global_param', 'frequency\n'],
        ['global_param', 'ntor\n'],
        ### --- Defining powers --- ###
        ['global_param', 'pow_ntor_i\n'],
        ### --- Defining grid --- ###
        ['grid_2d', 'grid_type\n'],
        ['grid_2d', 'r\n'],
        ['grid_2d', 'z\n'],
        ### --- Defining E-field --- ###
        ['local', 'e_plus\n'],
        ['local', 'e_plus_ph\n'],
        ['local', 'e_minus\n'],
        ['local', 'e_minus_ph\n'],
        ]

    # If not editting the data entry
    if not prev_line.split('%')[-2:] in to_edit:
        output_line = line

    ### ----- Defining admin data ----- ###
    elif prev_line.split('%')[-2:] == ['datainfo', 'dataprovider\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            + 'C. Perks'
            + '\n'
            )
    elif prev_line.split('%')[-2:] == ['datainfo', 'putdate\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            + date.today().strftime("%Y/%m/%d")
            + '\n'
            )
    elif prev_line.split('%')[-2:] == ['datainfo', 'source\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            + 'TORIC5'
            + '\n'
            )

    ### --- Defining wave --- ####
    elif prev_line.split('%')[-2:] == [' waves', 'coherentwave\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            )
    elif prev_line.split('%')[-2:] == ['type', 'id\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            + 'IC'
            + '\n'
            )
    elif prev_line.split('%')[-2:] == ['type', 'flag\n']:
        output_line = (
            '0'.rjust(12, ' ')
            + '\n'
            + '3'.rjust(12, ' ')
            + '\n'
            )

    ### --- Defining species --- ###
    elif prev_line.split('%')[-2:] == ['compositions', 'nuclei\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + str(len(dspec['amn'])).rjust(12, ' ')
            + '\n'
            )
        for ii in np.arange(len(dspec['amn'])):
            output_line += (
                ' waves%coherentwave%compositions%nuclei%zn\n'
                + '0'.rjust(12, ' ')
                + '\n'
                + "{:1.10F}".format(dspec['zn'][ii]).rjust(21, ' ')
                + '\n'
                + ' waves%coherentwave%compositions%nuclei%amn\n'
                + '0'.rjust(12, ' ')
                + '\n'
                + "{:1.10F}".format(dspec['amn'][ii]).rjust(21, ' ')
                + '\n'
                + ' waves%coherentwave%compositions%nuclei%label\n'
                + '-1'.rjust(12, ' ')
                + '\n'
                )

    elif prev_line.split('%')[-2:] == ['compositions', 'ions\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + str(len(dspec['amn'])).rjust(12, ' ')
            + '\n'
            )
        for ii in np.arange(len(dspec['amn'])):
            output_line += (
                ' waves%coherentwave%compositions%ions%nucindex\n'
                + '0'.rjust(12, ' ')
                + '\n'
                + str(int(ii+1)).rjust(12, ' ')
                + '\n'
                + ' waves%coherentwave%compositions%ions%zion\n'
                + '0'.rjust(12, ' ')
                + '\n'
                + "{:1.10F}".format(dspec['zion'][ii]).rjust(21, ' ')
                + '\n'
                + ' waves%coherentwave%compositions%ions%imp_flag\n'
                + '-1'.rjust(12, ' ')
                + '\n'
                + ' waves%coherentwave%compositions%ions%label\n'
                + '-1'.rjust(12, ' ')
                + '\n'
                )


    ### --- Defining antenna --- ###
    elif prev_line.split('%')[-2:] == ['global_param', 'frequency\n']:
        output_line = (
            '0'.rjust(12,' ')
            + '\n'
            + "{:1.10F}".format(dant['freq']).rjust(21, ' ')
            + '\n'
            )
    elif prev_line.split('%')[-2:] == ['global_param', 'ntor\n']:
        output_line = (
            '1'.rjust(12,' ')
            + '\n'
            + str(len(dant['ntor'])).rjust(12,' ')
            + '\n'
            )
        for ii in np.arange(len(dant['ntor'])):
            output_line += (
                "{:1.0F}".format(dant['ntor'][ii]).rjust(12, ' ')
                )
        output_line += '\n'

    ### --- Defining absorbed power --- ###
    elif prev_line.split('%')[-2:] == ['global_param', 'pow_ntor_i\n']:
        output_line = (
            '2'.rjust(12,' ')
            + '\n'
            + str(len(dant['ntor'])).rjust(12,' ') 
            + str(len(dspec['amn'])).rjust(12, ' ')
            + '\n'
            )
        for ii in np.arange(len(dspec['amn'])):
            for nn in np.arange(len(dant['ntor'])):
                if ii+1 == dspec['nspec']:
                    output_line += "{:1.10F}".format(RF_abs[nn]).rjust(21, ' ')
                else:
                    output_line += "{:1.10F}".format(0).rjust(21, ' ')
            output_line += '\n'

    ### --- Defining grid --- ###
    elif prev_line.split('%')[-2:] == ['grid_2d', 'grid_type\n']:
        output_line = (
            '0'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            )

    # [m]
    elif prev_line.split('%')[-2:] == ['grid_2d', 'r\n']:
        output_line = (
            '2'.rjust(12, ' ')
            +'\n'
            )

        output_line += _write_line(
            RR = toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
            ZZ = toric[0].cdf_hdl.variables['Zplasma'].data/100,
            data = toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
            data_type = 'R'
            )

    # [m]
    elif prev_line.split('%')[-2:] == ['grid_2d', 'z\n']:
        output_line = (
            '2'.rjust(12, ' ')
            +'\n'
            )

        output_line += _write_line(
            RR = toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
            ZZ = toric[0].cdf_hdl.variables['Zplasma'].data/100,
            data = toric[0].cdf_hdl.variables['Zplasma'].data/100,
            data_type = 'Z'
            )

    ### --- Defining E-field --- ###
    elif prev_line.split('%')[-2:] == ['local', 'e_plus\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +str(len(dant['ntor'])).rjust(12, ' ')
            )

        output_line += _write_line(
            RR = [
                toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
                toric[1].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj']
                ],
            ZZ = [
                toric[0].cdf_hdl.variables['Zplasma'].data/100,
                toric[1].cdf_hdl.variables['Zplasma'].data/100
                ],
            data = [
                np.sqrt(
                    toric[0].cdf_hdl.variables['Re2Eplus'].data**2
                    + toric[0].cdf_hdl.variables['Im2Eplus'].data**2
                    ),
                np.sqrt(
                    toric[1].cdf_hdl.variables['Re2Eplus'].data**2
                    + toric[1].cdf_hdl.variables['Im2Eplus'].data**2
                    )
                ]
            )

    elif prev_line.split('%')[-2:] == ['local', 'e_minus\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +str(len(dant['ntor'])).rjust(12, ' ')
            )

        output_line += _write_line(
            RR = [
                toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
                toric[1].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj']
                ],
            ZZ = [
                toric[0].cdf_hdl.variables['Zplasma'].data/100,
                toric[1].cdf_hdl.variables['Zplasma'].data/100
                ],
            data = [
                np.sqrt(
                    toric[0].cdf_hdl.variables['Re2Eminus'].data**2
                    + toric[0].cdf_hdl.variables['Im2Eminus'].data**2
                    ),
                np.sqrt(
                    toric[1].cdf_hdl.variables['Re2Eminus'].data**2
                    + toric[1].cdf_hdl.variables['Im2Eminus'].data**2
                    )
                ]
            )

    elif prev_line.split('%')[-2:] == ['local', 'e_plus_ph\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +str(len(dant['ntor'])).rjust(12, ' ')
            )

        output_line += _write_line(
            RR = [
                toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
                toric[1].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj']
                ],
            ZZ = [
                toric[0].cdf_hdl.variables['Zplasma'].data/100,
                toric[1].cdf_hdl.variables['Zplasma'].data/100
                ],
            data = [
                np.arctan(
                    toric[0].cdf_hdl.variables['Im2Eplus'].data
                    / toric[0].cdf_hdl.variables['Re2Eplus'].data
                    ),
                np.arctan(
                    toric[1].cdf_hdl.variables['Im2Eplus'].data
                    / toric[1].cdf_hdl.variables['Re2Eplus'].data
                    )
                ]
            )

    elif prev_line.split('%')[-2:] == ['local', 'e_minus_ph\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +str(len(dant['ntor'])).rjust(12, ' ')
            )

        output_line += _write_line(
            RR = [
                toric[0].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj'],
                toric[1].cdf_hdl.variables['Xplasma'].data/100 + dant['Rmaj']
                ],
            ZZ = [
                toric[0].cdf_hdl.variables['Zplasma'].data/100,
                toric[1].cdf_hdl.variables['Zplasma'].data/100
                ],
            data = [
                np.arctan(
                    toric[0].cdf_hdl.variables['Im2Eminus'].data
                    / toric[0].cdf_hdl.variables['Re2Eminus'].data
                    ),
                np.arctan(
                    toric[1].cdf_hdl.variables['Im2Eminus'].data
                    / toric[1].cdf_hdl.variables['Re2Eminus'].data
                    )
                ]
            )

    # Output
    return output_line, line

# Flattens 2d to one line string
def _write_line(
    # Data
    RR = None, # dim(nR,nZ)
    ZZ = None, # dim(nR, nZ)
    data=None, # dim(nR,nZ)
    data_type = None,
    # Interpolataion
    nR = None,
    nZ = None,
    # PLotting
    plt_all = False,
    ):

    # Prepares regular mesh to interpolate data on
    if isinstance(RR, list):
        if nR is None:
            nR = RR[0].shape[0]
        if nZ is None:
            nZ = RR[0].shape[1]

        rlim = [
            np.min((
                np.min(RR[0].flatten()),
                np.min(RR[1].flatten())
                )),
            np.max((
                np.max(RR[0].flatten()),
                np.max(RR[1].flatten())
                )),
            ]
        zlim = [
            np.min((
                np.min(ZZ[0].flatten()),
                np.min(ZZ[1].flatten())
                )),
            np.max((
                np.max(ZZ[0].flatten()),
                np.max(ZZ[1].flatten())
                )),
            ]

    else:
        if nR is None:
            nR = RR.shape[0]
        if nZ is None:
            nZ = RR.shape[1]

        rlim = [
            np.min(RR.flatten()),
            np.max(RR.flatten())
            ]
        zlim = [
            np.min(ZZ.flatten()),
            np.max(ZZ.flatten())
            ]

    ZZ2, RR2 = np.meshgrid(
        np.linspace(zlim[0], zlim[1], nZ),
        np.linspace(rlim[0], rlim[1], nR),
        )

    # Interpolates data
    data2 = {}
    if data_type == 'R':
        data2[0] = RR2
    elif data_type == 'Z':
        data2[0] = ZZ2
    else:
        data2[0] = LinearNDInterpolator(
            (RR[0].flatten(), ZZ[0].flatten()),
            data[0].flatten(),
            fill_value = 0.0
            )(
                (RR2.flatten(), ZZ2.flatten())
                ).reshape(RR2.shape)
        data2[1] = LinearNDInterpolator(
            (RR[1].flatten(), ZZ[1].flatten()),
            data[1].flatten(),
            fill_value = 0.0
            )(
                (RR2.flatten(), ZZ2.flatten())
                ).reshape(RR2.shape)

    # Plotting
    if plt_all:
        fig2, ax2 = plt.subplots()
        ax2.contourf(RR2,ZZ2, data2[0])
        ax2.set_xlabel('R [m]')
        ax2.set_ylabel('Z [m]')


    # Writes block
    out = (
        str(nR).rjust(12, ' ')
        + str(nZ).rjust(12, ' ')
        + '\n'
        )

    cnt = 0
    for yy in np.arange(nZ):
        for xx in np.arange(nR):
            for kk in data2.keys():
                out += (
                    "{:1.16E}".format(data2[kk][xx,yy]).rjust(26, ' ')
                    )
                cnt+=1
                if cnt%3 == 0:
                    out += '\n'

    if cnt%3 != 0:
        out += '\n'

    return out
