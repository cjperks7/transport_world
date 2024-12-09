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
#ascot = 'ascot_31256662.h5' # IEDF opt to thermal
#ascot = 'ascot_30990897.h5' # IEDF opt to RF hot tail

iascot = 'input.markerDist.h5'
#iascot = 'thermal.h5'
#zascot = 'ascot_31306336.h5'

ascot = 'ascot_33673825.h5'
zascot = 'ascot_33201509.h5' # zeroRF sim with flipped Bpol

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
zdata = utz._read_dist(
    h5_file = os.path.join(file, zascot)
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


#######################################################
#
#            EQDSK
#
#######################################################


# Plots eqdsk contours
edr = eqtools.EqdskReader(
    gfile=os.path.join(in_path, fgfile),
    afile=None
    )
edr_psiRZ = edr.getFluxGrid()[0]
edr_psiLCFS = edr.getFluxLCFS()[0]
edr_psi0 = edr.getFluxAxis()[0]
edr_rGrid = edr.getRGrid()
edr_zGrid = edr.getZGrid()
RR, ZZ = np.meshgrid(edr_rGrid, edr_zGrid)

edr_RLCFS = edr.getRLCFS()[0]
edr_ZLCFS = edr.getZLCFS()[0]

edr_rhop2D = np.sqrt(-1*(edr_psiRZ+edr_psi0)/(edr_psiLCFS-edr_psi0))

dedr = {}
dedr['RLCFS'] = edr_RLCFS
dedr['ZLCFS'] = edr_ZLCFS

dedr['rGrid'] = edr_rGrid
dedr['zGrid'] = edr_zGrid

dedr['rhop2D'] = edr_rhop2D


dedr['wall_R'] = edr.getMachineCrossSection()[0]
dedr['wall_Z'] = edr.getMachineCrossSection()[1]

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

# Plots ion energy distribution data in (pitch, E)-space
rhop_min = 0.7
rhop_max = 0.8
df = utz._plot_rzDist(
    ddata=ddata,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    plt_all = False
    )
idf = utz._plot_rzDist(
    ddata=idata,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    plt_all = False
    )

xx = utz._plot_rzDist(
    idf = idf,
    df=df,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    gfile = os.path.join(in_path, fgfile),
    plt_whole = False
    )


# Plots ion energy distribution data in (Vpara, Vperp)-space
rhop_min = 0.6
rhop_max = 0.7
dfv = utz._plot_rzVDist(
    ddata=ddata,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    plt_all = False
    )
idfv = utz._plot_rzVDist(
    ddata=idata,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    plt_all = False
    )

xx = utz._plot_rzVDist(
    idf = idfv,
    df=dfv,
    rhop_min = rhop_min,
    rhop_max = rhop_max,
    gfile = os.path.join(in_path, fgfile)
    )








geq = omfit_eqdsk.OMFITgeqdsk(os.path.join(in_path, fgfile))
rmin_g = geq['fluxSurfaces']['geo']['a']
rhop_g = np.sqrt(
    geq['fluxSurfaces']['geo']['psin']
    /geq['fluxSurfaces']['geo']['psin'][-1]
    )

rmin = interp1d(
    rhop_g, rmin_g
    )(0.65)

# Critical angle
cos_thetaC = np.sqrt(2*rmin/geq['RMAXIS']) # [v_||/v]
thetaC = np.arccos(cos_thetaC)


vpe_1d = np.linspace(0,3e5, 1001)
vpa_1d = np.linspace(-3e5, 3e5, 2001)
MM = 40 *cnt.m_u
Ti = 1e3 *cnt.e

vpe_2d, vpa_2d = np.meshgrid(vpe_1d,vpa_1d)

'''
ff = MM/Ti *vpe *np.exp(-MM*vpe**2/2/Ti) *3e16

gg, aa = plt.subplots()

aa.plot(
    vpe,ff
    )
aa.set_yscale('log')
'''

ff = (
    (MM/2/np.pi/Ti)**(3/2)
    *vpe_2d
    *np.exp(-MM*vpa_2d**2/2/Ti)
    *np.exp(-MM*vpe_2d**2/2/Ti)
    ) *3e16

energy = 0.5*MM*(vpa_2d**2+vpe_2d**2)/cnt.e/1e3 # [keV]
#mask = np.where(energy > 8.0)
#ff[mask] = np.nan

gg,aa = plt.subplots()

pcm = aa.pcolormesh(
    vpa_2d,
    vpe_2d,
    np.log10(ff)
    )
cbar = plt.colorbar(pcm, ax=aa)
pcm.set_clim(-2,6.5)
cbar.set_label(r'$log_{10}(f_{Ar16+})$ [$log_{10}(s^2/m^5)$]')
aa.grid('on')
aa.yaxis.set_major_formatter(formatter)
aa.xaxis.set_major_formatter(formatter)

aa.set_xlabel(r'$v_{||}$ [$m/s$]')
aa.set_ylabel(r'$v_{\perp}$ [$m/s$]')

con = aa.contour(
    vpa_2d, vpe_2d,
    energy,
    [8],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    colors= 'k'
    )
    
aa.clabel(con, inline=True, fontsize = 8, fmt = '8 keV')
aa.set_aspect('equal')

aa.plot(
    np.r_[0, 3e5/np.tan(thetaC)],
    np.r_[0,3e5],
    'k--'
    )
aa.plot(
    np.r_[0, -3e5/np.tan(thetaC)],
    np.r_[0,3e5],
    'k--'
    )
aa.set_title(r'theory; $T_i$ =1keV; $n_{Ar16+}$ =3e16 $m^{-3}$')
#aa.set_xlim(-1.5e5, 1.5e5)
#aa.set_ylim(0, 1.5e4)













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

'''
ind0 = np.where(ddata['endstate']['energy']['vals'] <= 1e2)
ind1 = np.where(ddata['endstate']['energy']['vals'] > 1e2)
bins = [50,100]

fig, ax = plt.subplots(nrows=1,ncols=2, figsize=(14,8))

hist = ax[0].hist2d(
    ddata['endstate']['R']['vals'][ind0[0]],
    ddata['endstate']['z']['vals'][ind0[0]],
    bins = bins,
    #cmap = 'viridis'
    )
cb = plt.colorbar(hist[3], ax=ax[0])
cb.set_label('Counts')
ax[0].set_title('End energy <= 100eV')

hist = ax[1].hist2d(
    ddata['endstate']['R']['vals'][ind1[0]],
    ddata['endstate']['z']['vals'][ind1[0]],
    bins = bins,
    #cmap = 'viridis'
    )
cb = plt.colorbar(hist[3], ax=ax[1])
cb.set_label('Counts')
ax[1].set_title('End energy > 100eV')

for ii in np.arange(2):
    con = ax[ii].contour(
        dedr['rGrid'],dedr['zGrid'],
        dedr['rhop2D'],
        [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
        linestyles='solid',
        linewidth=2,
        zorder=2,
        colors= 'k'
        )
        
    ax[ii].clabel(con, inline=True, fontsize = 8)
    ax[ii].grid('on')
    ax[ii].set_xlabel('R [m]')
    ax[ii].set_ylabel('Z [m]')

    ax[ii].set_ylim(-0.45, 0.45)
    ax[ii].set_xlim(0.44, 0.92)
    ax[ii].set_aspect('equal')
'''




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

key = 'R'
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








fig, ax = plt.subplots(nrows=1,ncols=2, figsize=(14,8))

bins = [50,100]

histi = ax[0].hist2d(
    ddata['inistate']['R']['vals'],
    ddata['inistate']['z']['vals'],
    bins = bins,
    #cmap = 'viridis'
    )
cb = plt.colorbar(histi[3], ax=ax[0])
cb.set_label('Counts')

hist = ax[1].hist2d(
    ddata['endstate']['R']['vals'],
    ddata['endstate']['z']['vals'],
    bins = bins,
    #cmap = 'viridis'
    )
cb = plt.colorbar(hist[3], ax=ax[1])
cb.set_label('Counts')

for ii in np.arange(2):
    ax[ii].plot(
        dedr['wall_R'], 
        dedr['wall_Z'],
        'k',
        linewidth=3
        )

    ax[ii].plot(
        dedr['RLCFS'], 
        dedr['ZLCFS'],
        'r',
        linewidth=3
        )

    con = ax[ii].contour(
        dedr['rGrid'],dedr['zGrid'],
        dedr['rhop2D'],
        [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
        linestyles='solid',
        linewidth=2,
        zorder=2,
        colors= 'k'
        )
        
    ax[ii].clabel(con, inline=True, fontsize = 8)

ax[0].set_ylabel('Z [m]')
ax[0].set_xlabel('R [m]')

ax[1].set_ylabel('Z [m]')
ax[1].set_xlabel('R [m]')

ax[0].set_title('init')
ax[1].set_title('end')

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')

Rmax = np.max((
    np.max(histi[1]),
    np.max(hist[1])
    ))
Rmin = np.min((
    np.min(histi[1]),
    np.min(hist[1])
    ))
Zmax = np.max((
    np.max(histi[2]),
    np.max(hist[2])
    ))
Zmin = np.min((
    np.min(histi[2]),
    np.min(hist[2])
    ))

ax[0].set_xlim(Rmin,Rmax)
ax[1].set_xlim(Rmin,Rmax)

ax[0].set_ylim(Zmin, Zmax)
ax[1].set_ylim(Zmin,Zmax)