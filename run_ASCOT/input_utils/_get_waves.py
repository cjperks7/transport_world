'''

_get_waves.py writes on RFOF E-field input file

cjperks
Oct 2nd, 2023

'''

# Modules
import numpy as np
import os
from transport_world.run_TORIC.output_utils import toric_tools

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
    in_path = None,
    device = 'CMOD',
    RF_abs = 4.0e3, # [W]
    ):

    # Path to skeleton file
    infile = '/home/cjperks/2201_Pumpout/ASCOT_template/skeleton_waves.ascii'
    #infile = '/home/cjperks/ASCOT4/ascot4-trunk/rfof/input/itm_waves_from_write_cpo.ascii'

    # Path to write RFOF input file
    #outfile = '/home/cjperks/2201_Pumpout/CMOD/shots/1140221013/input_ASCOT/itm_waves.ascii'
    outfile = os.path.join(
        in_path,
        'input_ASCOT',
        'itm_waves.ascii'
        )

    # Reads TORIC output
    #toric = toric_tools.toric_analysis(
    #    toric_name='/home/cjperks/2201_Pumpout/CMOD/shots/1140221013/toric_test.ncdf',
    #    mode='ICRF',
    #    path='/home/cjperks/2201_Pumpout/CMOD/shots/1140221013/'
    #    )
    toric = toric_tools.toric_analysis(
        toric_name=os.path.join(in_path, 'toric_test.ncdf'),
        mode='ICRF',
        path=in_path,
        )

    # Antenna data
    if device == 'CMOD':
        dant = {
            'freq': int(80e6),
            'ntor': int(10)
            }
    elif device == 'SPARC':
        dant = 0

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
    RF_abs=None, # [W]
    ):

    # List of arrays to edit
    to_edit = [
        [' waves', 'coherentwave\n'],
        #['type','flag\n'],
        #['type','id\n'],
        #['wave_id', 'name\n'],
        #['wave_id','index\n'],
        #['compositions', 'amn\n'],
        #['compositions', 'zn\n'],
        #['compositions', 'zion\n'],
        #['compositions', 'nucindex\n'],
        ['global_param', 'frequency\n'],
        ['global_param', 'ntor\n'],
        ['global_param', 'pow_ntor_i\n'],
        ['grid_2d', 'grid_type\n'],
        ['grid_2d', 'r\n'],
        ['grid_2d', 'z\n'],
        ['local', 'e_plus\n'],
        ['local', 'e_plus_ph\n'],
        ['local', 'e_minus\n'],
        ['local', 'e_minus_ph\n'],
        ]

    # If not editting the data entry
    if not prev_line.split('%')[-2:] in to_edit:
        output_line = line

    # Defining absorbed power
    elif prev_line.split('%')[-2:] == ['global_param', 'pow_ntor_i\n']:
        output_line = (
            '1'.rjust(12,' ')
            + '\n'
            + '1'.rjust(12,' ')
            + '\n'
            + "{:1.4E}".format(RF_abs).rjust(16, ' ')
            + '\n'
            )

    # Defining wave
    elif prev_line.split('%')[-2:] == [' waves', 'coherentwave\n']:
        output_line = (
            '1'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            )

    # Antenna data
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
            + '1'.rjust(12,' ')
            + '\n'
            + "{:1.0F}".format(int(dant['ntor'])).rjust(12, ' ')
            + '\n'
            )

    # (R, Z) mesh
    elif prev_line.split('%')[-2:] == ['grid_2d', 'grid_type\n']:
        output_line = (
            '0'.rjust(12, ' ')
            + '\n'
            + '1'.rjust(12, ' ')
            + '\n'
            )

    elif prev_line.split('%')[-2:] == ['grid_2d', 'r\n']:
        output_line = (
            '2'.rjust(12, ' ')
            +'\n'
            )

        output_line += _write_line(
            data = toric.cdf_hdl.variables['Xplasma'].data
            )

    elif prev_line.split('%')[-2:] == ['grid_2d', 'z\n']:
        output_line = (
            '2'.rjust(12, ' ')
            +'\n'
            )

        output_line += _write_line(
            data = toric.cdf_hdl.variables['Zplasma'].data
            )

    # Input E-fields in polar form
    elif prev_line.split('%')[-2:] == ['local', 'e_plus\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +'1'.rjust(12, ' ')
            )

        output_line += _write_line(
            data = np.sqrt(
                toric.cdf_hdl.variables['Re2Eplus'].data**2
                + toric.cdf_hdl.variables['Im2Eplus'].data**2
                )
            )

    elif prev_line.split('%')[-2:] == ['local', 'e_minus\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +'1'.rjust(12, ' ')
            )

        output_line += _write_line(
            data = np.sqrt(
                toric.cdf_hdl.variables['Re2Eminus'].data**2
                + toric.cdf_hdl.variables['Im2Eminus'].data**2
                )
            )

    elif prev_line.split('%')[-2:] == ['local', 'e_plus_ph\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +'1'.rjust(12, ' ')
            )

        output_line += _write_line(
            data = np.arctan(
                toric.cdf_hdl.variables['Im2Eplus'].data
                / toric.cdf_hdl.variables['Re2Eplus'].data
                )
            )

    elif prev_line.split('%')[-2:] == ['local', 'e_minus_ph\n']:
        output_line = (
            '3'.rjust(12, ' ')
            +'\n'
            +'1'.rjust(12, ' ')
            )

        output_line += _write_line(
            data = np.arctan(
                toric.cdf_hdl.variables['Im2Eminus'].data
                / toric.cdf_hdl.variables['Re2Eminus'].data
                )
            )

    # Output
    return output_line, line

# Flattens 2d to one line string
def _write_line(
    data=None,
    ):

    out = (
        str(data.shape[0]).rjust(12, ' ')
        + str(data.shape[1]).rjust(12, ' ')
        + '\n'
        )

    for xx in np.arange(data.shape[0]):
        for yy in np.arange(data.shape[1]):
            out += (
                "{:1.16E}".format(data[xx,yy]).rjust(26, ' ')
                )

    out += '\n'

    return out


