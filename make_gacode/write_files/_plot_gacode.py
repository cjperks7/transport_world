'''

_plot_gacode is a function to conveninetly 
plot gacode files

'''

# Modules
from portals.gacode_tools import PROFILESmodule

__all__ = [
    'plot_ga',
    ]

##################################################
#
#           Main
#
##################################################

def plot_ga(
    fgacode=None,
    ):

    # Creates a profiles class
    prof = PROFILESmodule.PROFILES_GACODE(fgacode)

    # Plots
    prof.plot()