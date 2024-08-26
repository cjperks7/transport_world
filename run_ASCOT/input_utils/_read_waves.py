'''

Function to read itm_waves.ascii files

cjperks
July 26, 2024

'''

# Modules
import os, sys
import numpy as np
import re

__all__ = [
    'read_waves'
    ]

#########################################
#
#             Main
#
#########################################

# Reads itm_waves.ascii file
def read_waves(
    file = None,
    ):

    # Init 
    dout = {}
    data = False

    # Reads file
    ff = open(file, 'r')

    #ddd = 0
    for line in ff.readlines():
        #ddd += 1
        #if ddd > 40:
        #    continue

        # Organizes data block
        if data:
            # Reads the rank
            if cnt == 0:
                rank = int(line.strip())
                cnt += 1

                # Adds data
                tmp = dout
                for ii, key in enumerate(keys):
                    tmp = tmp[key]
                    if ii == len(keys)-1:
                        tmp['rank'] = rank

                # If skipping this data
                if rank == -1:
                    data = False
                elif rank == 0:
                    cnt += 1
                    # Adds data
                    tmp = dout
                    for ii, key in enumerate(keys):
                        tmp = tmp[key]
                        if ii == len(keys)-1:
                            tmp['dim'] = np.r_[1]

            # Reads the dimensions
            elif cnt == 1:
                dims = np.asarray([int(xx) for xx in line.split(' ') if len(xx)>0])
                cnt += 1

                # Adds data
                tmp = dout
                for ii, key in enumerate(keys):
                    tmp = tmp[key]
                    if ii == len(keys)-1:
                        tmp['dim'] = dims

            # Reads the data block
            elif cnt == 2:
                # If at next header
                if line.split('%')[0].strip() == 'waves':
                    data = False

                    # Adds data
                    tmp = dout
                    for ii, key in enumerate(keys):
                        tmp = tmp[key]
                        if ii == len(keys)-1 and len(vals)>0:
                            tmp['data'] = np.asarray(vals).reshape(tmp['dim'])

                # Adds data
                else:
                    for xx in line.split(' '):
                        yy = xx.strip()
                        if len(yy) == 0:
                            continue

                        if _is_number(yy):
                            vals += [float(yy)]
                        else:
                            vals += [yy]


        # If reading a header
        if line.split('%')[0].strip() == 'waves' and not data:
            data = True
            keys = []
            cnt = 0
            vals = []

            # Organizes header data
            for ii in np.arange(1, len(line.split('%'))):
                keys.append(line.split('%')[ii].strip())
            _create_nested_dict(dout, keys)


    # Closes file
    ff.close()

    # Output
    return dout

        
###################################
#
#           Utilities
#
###################################

# Checks if a string is words or numbers
def _is_number(string):
    # This regex pattern matches integers, floating-point numbers, and scientific notation
    pattern = re.compile(r'^-?\d+(\.\d+)?([eE][-+]?\d+)?$')
    return bool(pattern.match(string))

def _create_nested_dict(d, keys):
    for key in keys:
        d = d.setdefault(key, {})