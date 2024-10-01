'''

Reads distributions calculated by ASCOT

cjperks
Sep 26, 2024

'''

# Modules
import sys, os
import scipy.constants as cnt
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.ticker import ScalarFormatter

formatter = ScalarFormatter()
formatter.set_powerlimits((-2, 3))  # Change the limits as needed

plt.rcParams.update({'font.size': 16})

from transport_world.run_ASCOT import output_utils as utz
from transport_world import run_AURORA as rA
from transport_world.run_TORIC.output_utils import toric_tools

sys.path.insert(0,'/home/cjperks/usr/python3modules/eqtools3')
import eqtools
sys.path.pop(0)

from omfit_classes import omfit_eqdsk

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# File management
shot = '1140221012'
file = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD/shots',
    shot,
    'ASCOT/output'
    )
ascot = 'ascot_30990897.h5'
#iascot = 'input.markerDist.h5'
iascot = 'thermal.h5'

in_path = os.path.join(
    '/home/cjperks',
    'work/2201_Pumpout',
    'CMOD/shots',
    shot,
    'profiles'
    )
fgacode = 'input_t1.gacode'
fgfile = 'g'+shot+'.01000'

'''


#### --- End states --- ####

# Time to end state
end_t = ff['endstate']['time'][:]

'''

# Reads output distribution data
ddata = utz._read_dist(
    h5_file = os.path.join(file, ascot)
    )
idata = utz._read_dist(
    h5_file = os.path.join(file, iascot)
    )



#######################################################
#
#            Impurity Modeling
#
#######################################################


dmodel = rA.def_dmodel(
    option='Ar_TGLF',
    shot=shot,
    fgacode=fgacode,
    )

dimp = rA.calc_imp_dens(
    dmodel = dmodel,
    plt_all = False,
    )


#######################################################
#
#            TORIC Modeling
#
#######################################################

# Plots eqdsk contours
edr = eqtools.EqdskReader(
    gfile=os.path.join(in_path, fgfile),
    afile=os.path.join(in_path, 'a'+fgfile[1:])
    )

filet = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD',
    'shots',
    shot,
    'TORIC'
    )

ncd = [
    '/toric_pnphi.ncdf',
    '/toric_nnphi.ncdf'
    ]
pwrs = np.r_[
    675,
    877
    ]/1e6 # [MW]
MW_abs = 500e3/1e6/2*0.3

pwr_tot = np.zeros(dimp['rhop_fm'].shape[0])

for kk, cd in enumerate(ncd):

    # Reads TORIC output
    toric = toric_tools.toric_analysis(
        toric_name=filet+ cd,
        mode='ICRF',
        path=filet+'/'
        )
    # toric.cdf_hdl.variables['PwIH'].units --> [W/cm3/MW_abs]
    # toric.cdf_hdl.variables['Pw_abscissa'].long_name --> r/a

    trhop = edr.roa2psinorm(
        toric.cdf_hdl.variables['Pw_abscissa'].data,
        0,
        sqrt = True,
        )

    pwr_tot += interp1d(
        trhop,
        (
            toric.cdf_hdl.variables['PwIH'].data[:,2]*1e6 # [W/m3/MW_abs]
            *MW_abs
            ),
        bounds_error=False,
        fill_value = 0.0
        )(dimp['rhop_fm']) # [W/m3]




################################################################
#
#           Plotting
#
################################################################

# Plots evolution of density v. time
utz._plot_rhoDist(
    ind_dist = 1,
    ddata = ddata,
    dext = {
        'rhop': dimp['rhop_fm'],
        'vals': dimp['nz_fm'][:,16,-1]*1e6,
        'label': 'RF off (AURORA+TGLF)'
        },
    ylabel = r'$n_{Ar16+}$ [$1/m^3$]'
    )

# Plots evolution of jxB torque v. time
utz._plot_rhoDist(
    ind_dist = 5,
    ddata = ddata,
    ylabel = r'$j\times B$ torque [$N/m^2$]'
    )

# Plots evolution of absorbed ICRF power v. time
utz._plot_rhoDist(
    ind_dist = 21,
    ddata = ddata,
    dext = {
        'rhop': dimp['rhop_fm'],
        'vals': pwr_tot,
        'label': 'TORIC'
        },
    ylabel = r'$Q_{RF, Ar16+}$ [$W/m^3$]'
    )



# Plots evolution of other v. time
utz._plot_rhoDist(
    ind_dist = 27,
    ddata = ddata,
    ylabel = r'Power deposition to Lumped [$W/m^3$]'
    )

#apwr = np.sum(
#    (
#        ddata['rhoDist']['dists'][21]['vals'][..., 0]
#        *ddata['rhoDist']['shellVolume']['vals'][:,None]
#        /ddata['rhoDist']['cents']['dim2']['dbin']
#        ),
#    axis = 0
#    ) # [W]

#tpwr_0d = np.sum(pwrs*1e6) # [W]

# Plots ion energy distribution data
rhop_min = 0.9
rhop_max = 1.0
df = utz._plot_rzDist(
    ddata=ddata,
    rhop_min = rhop_min,
    rhop_max = rhop_max
    )
idf = utz._plot_rzDist(
    ddata=idata,
    rhop_min = rhop_min,
    rhop_max = rhop_max
    )

xx = utz._plot_rzDist(
    idf = idf,
    df=df,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    gfile = os.path.join(in_path, fgfile)
    )






'''
from scipy.interpolate import interp1d

RF_off = interp1d(
    dimp['rhop_fm'],
    dimp['nz_fm'][:,16,-1]*1e6
    )(ddata['rhoDist']['cents']['dim1']['vals'])
RF_on = tmp['vals'][:,-1,0]/20e-3

dn = (
    np.sum((RF_off-RF_on)*ddata['rhoDist']['shellVolume']['vals'])
    /np.sum(RF_off*ddata['rhoDist']['shellVolume']['vals'])
    )
'''


fig, ax = plt.subplots()

hist = ax.hist(
    ddata['endstate']['endcond']['vals'],
    bins = np.arange(
        np.min(ddata['endstate']['endcond']['vals']),
        np.max(ddata['endstate']['endcond']['vals'])+2
        ),
    weights = 100/102400 *np.ones(ddata['endstate']['endcond']['vals'].size)
    )
ax.set_xticks(0.5*(hist[1][1:]+hist[1][:-1]))
ax.set_xticklabels([f'{int(xx)}' for xx in hist[1][:-1]])

ax.set_xlabel('End condition')
ax.set_ylabel('frequency [%]')
ax.set_title('endcond')
ax.grid('on')
ax.yaxis.set_major_formatter(formatter)


fig, ax = plt.subplots()

bins = np.linspace(0,1.2,13)
histi = ax.hist(
    ddata['inistate']['rho']['vals'],
    bins = bins,
    alpha = 0.7,
    label = 'init',
    weights = 100/102400 *np.ones(ddata['inistate']['rho']['vals'].size)
    )

hist = ax.hist(
    ddata['endstate']['rho']['vals'],
    bins = bins,
    alpha = 0.7,
    label = 'end',
    weights = 100/102400 *np.ones(ddata['endstate']['rho']['vals'].size)
    )

leg = ax.legend()

ax.set_title(
    r'$n_{init,core}$=%i; $n_{end,core}$=%i; lost = %0.2f %%'%(
        np.sum(histi[0][:-2]),
        np.sum(hist[0][:-2]),
        (1-np.sum(hist[0][:-2])/np.sum(histi[0][:-2]))*100
        )
    )

ax.set_xlabel(r'$\sqrt{\Psi_n}$')
ax.set_ylabel('frequency [%]')
ax.grid('on')
ax.set_title('rho')


key = 'energy'
log = False
fig, ax = plt.subplots()

histi = ax.hist(
    ddata['inistate'][key]['vals'],
    bins = np.logspace(
        np.log10(np.min(ddata['inistate'][key]['vals'])),
        np.log10(np.max(ddata['inistate'][key]['vals'])),
        31
        ),
    alpha = 0.7,
    label = 'init',
    log = log,
    weights = 100/102400 *np.ones(ddata['inistate'][key]['vals'].size)
    )

hist = ax.hist(
    ddata['endstate'][key]['vals'],
    bins = np.logspace(
        np.log10(np.min(ddata['endstate'][key]['vals'])),
        np.log10(np.max(ddata['endstate'][key]['vals'])),
        81
        ),
    alpha = 0.7,
    label = 'end',
    log = log,
    weights = 100/102400 *np.ones(ddata['endstate'][key]['vals'].size)
    )
leg = ax.legend()
ax.grid('on')
ax.set_xscale('log')

ax.set_xlabel('Energy [eV]')
ax.set_ylabel('frequency [%]')

ax.set_title(key)




key = 'time'
fig, ax = plt.subplots()

hist = ax.hist(
    ddata['endstate'][key]['vals'],
    bins = 31,
    alpha = 0.7,
    label = 'end',
    weights = 100/102400 *np.ones(ddata['endstate'][key]['vals'].size)
    )
leg = ax.legend()
ax.grid('on')

ax.set_xlabel('time [s]')
ax.set_ylabel('frequency [%]')
ax.set_title('time')
ax.set_yscale('log')


# Useless = [
#   'Rprt', 'charge', 'dH', 'fPoloidal', 'fToroidal', 'id',
#   'mhdPhase', 'npassing', 'ntrapped', 'origin', 'phiprt',
#   'species', 'wallTile', 'weight', 'zprt',
#   ]
# Interesting = [
#   'R', 'cpuTime', 'phi', 'pitch', 'vR', 'vphi', 'vz',
#   'z'!!!!, 
#   ]

key = 'pitch'
fig, ax = plt.subplots()

histi = ax.hist(
    ddata['inistate'][key]['vals'],
    bins = 30,
    alpha = 0.7,
    label = 'init',
    weights = 100/102400 *np.ones(ddata['inistate'][key]['vals'].size)
    )

hist = ax.hist(
    ddata['endstate'][key]['vals'],
    bins = 30,
    alpha = 0.7,
    label = 'end',
    weights = 100/102400 *np.ones(ddata['endstate'][key]['vals'].size)
    )

leg = ax.legend()
ax.grid('on')
ax.set_ylabel('frequency [%]')

ax.set_title(key)
ax.xaxis.set_major_formatter(formatter)