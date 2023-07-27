'''

Function managing running TGYRO from an input.gacode file

'''

# Modules
from portals.gacode_tools import TGYROmodule,PROFILESmodule
from portals.misc_tools   import IOtools

__all__ = [
    'run_TGYRO',
    ]

######################################################
#
#               Main
#
######################################################

def run_TGYRO(
    # File management
    fgacode = None,             # 'xxx/input.gacode'
    folder = None,              # 'xxx/'
    subfolder = None,           # 'xxx
    # Engaging settings
    minutesJob = 45,
    # Simulation settings
    # http://gafusion.github.io/doc/tgyro/tgyro_list.html#
    iterations = 20,            # number of solver iterations
    vectorRange = [0.3,0.8,4],  # [min(rhot), max(rhot), num. points]
    PredictionSet = [0,1,0],    # evolve [Te, Ti, ne]
    TGLFsettings = 4,           # 1 -> SAT0, 2 -> SAT1, 3 -> SAT1geo, 4 -> SAT2, 5 -> SAT2em
    solver = {                  # numerics settings
        'step_jac':    1E-2,        # step length used for finite-differencing
        'step_max':    1E-2,        # maximum step of any Newton step
        'res_method':     2,        # formula for resdiual in root finding
        'tgyro_method':   6,        # variations of Newton method used in root finding
        'relax_param':   0.1        # contols shrinkage of relaxation paramter
        },
    physics_options = {         # physics settings
        'TargetType':2              # 1 -> fixed sources, 2 -> evolve e-i exchange, 3 -> evolve fusion, radiation, e-i exchange
        }, 
    ):

    # Finds input.gacode file if blank
    if fgacode is None:
        fgacode = folder + 'input.gacode'

    # Creates a PROFILES class
    prof = PROFILESmodule.PROFILES_GACODE(gacode_file)
    #prof.plot()

    # Creates a TGYRO class
    tgyro = TGYROmodule.TGYRO()
    tgyro.prep(folder,profilesclass_custom=prof)

    # Runs TGYRO
    tgyro.run( 
        subFolderTGYRO        = subfolder +'/',
        iterations            = iterations,
        vectorRange           = vectorRange,
        PredictionSet         = PredictionSet,
        TGLFsettings          = TGLFsettings,
        TGYRO_solver_options  = solver,
        TGYRO_physics_options = physics_options,
        minutesJob            = minutesJob,
        )

    # Reads results
    tgyro.read(label=subfolder)

    # Plots results
    tgyro.plotRun(labels=[subfolder])
