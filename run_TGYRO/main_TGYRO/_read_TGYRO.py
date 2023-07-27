'''

Reads TGYRO output

'''

# Modules
from portals.gacode_tools import TGYROmodule,PROFILESmodule
from portals.misc_tools   import IOtools

__all__ = [
    'read_TGYRO',
    ]

######################################################
#
#               Main
#
######################################################

def read_TGYRO(
    fgacode = None,
    folder = None,
    ):

    # Finds input.gacode file if blank
    if fgacode is None:
        fgacode = folder + 'input.gacode'

    # Loads data
    profiles    = PROFILESmodule.PROFILES_GACODE(fgacode)
    tgyro_out   = TGYROmodule.TGYROoutput(folder,profiles=profiles)

    # Plots
    tgyro_out.plot()