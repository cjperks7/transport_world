'''

create_gacode.py is an example script to run this gacode writer

cjperks
06/01/23

'''

from transport_world import make_gacode as mkga

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2


#######################################################
#
#            User controls
#
#######################################################

# Experiment of interest
shot = 1140221013
t0 = 1.2 # [s]
dt = 0.05

# Ion concentration modeling
H2D = 0.3
dmodel = {
    'option': 'concentration',
    'ions': {
        'D': {
            'con': 0.7,
        },
        'H': {
            'con': H2D
        },
        'Ar': {
            'con': 1e-3,
        },
        'Mo': {
            'con': 3e-4,
        },
    },
    'rotation':{
        'option': 'zero',
    },
}
dmodel['ions']['H']['con'] *= dmodel['ions']['D']['con']

# Paths to fitted profiles and gfile
path_input = '/home/cjperks/transp_runs/CMOD/'+str(shot)
path_kin = '/prof_'+str(shot)+'.npy'
#gfile_d = '/g'+str(shot)+'.00999_968'
path_gfile = '/g'+str(shot)+'.01000'


#######################################################
#
#            Get experimental profiles
#
#######################################################

# Loads fitted kinetic profiles
dout = mkga.get_fits(
    path_input = path_input,
    path_kin = path_kin,
    t0 = t0,
    dt = dt,
    plt_all = False,
    )

# Loads magnetic equilibrium
dout = mkga.get_geq(
    dout = dout,
    path_input = path_input,
    path_gfile = path_gfile,
    plt_all = False,
    )

#######################################################
#
#            Calculates profiles
#
#######################################################


# Facilitates ion density modeling
dout = mkga.get_ions(
    dout=dout,
    dmodel=dmodel,
    plt_all = True,
    )

# Calculates resultant power densities
dout = mkga.get_powers(
    dout=dout,
    verb=2,
    Prf_in = 1, # [MW]
    plt_all = True,
    )

#######################################################
#
#            Loads C-Mod data
#
#######################################################

# Loads CMOD measuremens
dout = mkga.get_cmod(
    dout=dout,
    shot=shot,
    quants=None, # None -> do all
    )

# Comparison plots against CMOD values
mkga.plt_cmod(dout=dout)

# Rescales input profiles to CMOD values
dout = mkga.rscl_cmod(
    dout=dout,
    dquants={
        'fus': 0,
        'rad': 1/8,#3,
        'ohm': 5,#0.5,
        'rf': 5,#0.5,
        'Zeff': 1,
        },
    replot = True,
    )


#######################################################
#
#            Writes GACODE file
#
#######################################################

# Write GACODE file
mkga.write_ga(dout=dout)