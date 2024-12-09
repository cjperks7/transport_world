'''

Script to compare results from ACTION_profiles in SPIRAL
against input profile data

cjperks
Nov 21, 2024

'''

# Modules
import sys, os
from scipy.interpolate import griddata
from scipy.spatial import ConvexHull, Delaunay
import matplotlib.gridspec as gridspec
import scipy.constants as cnt
from omfit_classes import omfit_eqdsk, omfit_gapy
import matplotlib.ticker as ticker

from transport_world import run_AURORA as rA
#from transport_world.run_SPIRAL.output_utils import _read_prf as rp
from transport_world.run_TORIC.output_utils import toric_tools

# Enables automatic reloading of modules
#%reload_ext autoreload
#%autoreload 2

sys.path.insert(0,'/home/cjperks/usr/python3modules/eqtools3')
import eqtools
sys.path.pop(0)

# File management
shot = 1140221012
prf = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD/shots',
    '%i'%(shot),
    'SPIRAL',
    'test_%i_01_prf_results.cdf'%(shot)
    )

# Paths
shot = '1140221012'
in_path = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout',
    'CMOD/shots',
    shot,
    'profiles'
    )
fgacode = 'input_t1.gacode'
Ar_title = 'Ar16+'


inga = omfit_gapy.OMFITgacode(
    os.path.join(
        in_path,
        fgacode,
        )
    )

# CMOD geom
Rmaj = 68.491/100# [m]
Zmax = [-40/100,40/100] # [m]


import netCDF4

ff = netCDF4.Dataset(prf, 'r')

# Init
ddata = {}

# Gets mesh
ddata['r_1d'] = {}
ddata['r_1d']['data'] = ff.variables['profile_r_axis'][:] # dim(nr,)
ddata['r_1d']['units'] = 'm'

ddata['z_1d'] = {}
ddata['z_1d']['data'] = ff.variables['profile_z_axis'][:] # dim(nz,)
ddata['z_1d']['units'] = 'm'

# Calculate the bin edges
ddata['r_edges'] = {}
bin_width = ddata['r_1d']['data'][1] - ddata['r_1d']['data'][0]
ddata['r_edges']['data'] = np.concatenate((
    [ddata['r_1d']['data'][0] - bin_width / 2],
    ddata['r_1d']['data'] + bin_width / 2,
    ))
ddata['z_edges'] = {}
bin_width = ddata['z_1d']['data'][1] - ddata['z_1d']['data'][0]
ddata['z_edges']['data'] = np.concatenate((
    [ddata['z_1d']['data'][0] - bin_width / 2],
    ddata['z_1d']['data'] + bin_width / 2,
    ))


ddata['r_2d'] = {}
ddata['z_2d'] = {}
ddata['r_2d']['units'] = ddata['z_2d']['units'] = 'm'
ddata['z_2d']['data'], ddata['r_2d']['data'] = np.meshgrid(
    ddata['z_1d']['data'], ddata['r_1d']['data']
    )

ddata['phi_1d'] = {}
ddata['phi_1d']['data'] = ff.variables['profile_phi_axis'][:] # dim(nphi,)
ddata['phi_1d']['units'] = 'deg'

# Gets fast ion density
ddata['density'] = {}
ddata['density']['data'] = ff.variables['profile_fast-ion_density'][:] # dim(nphi, nz,nr)
ddata['density']['units'] = '1/m^3'

# Gets fast ion effective temperature
ddata['T_eff'] = {}
ddata['T_eff']['data'] = ff.variables['profile_fast-ion_effective_temperature'][:] # dim(nphi, nz,nr)
ddata['T_eff']['units'] = 'keV'

# Gets voxel volume
ddata['volume'] = {}
ddata['volume']['data'] = ff.variables['profile_volume'][:] # dim(nphi, nz,nr)
ddata['volume']['units'] = 'm^3'

# Gets mode exchange energy
ddata['mode_energy'] = {}
ddata['mode_energy']['data'] = ff.variables['profile_mode_energy'][:] # dim(nphi, nz,nr)
ddata['mode_energy']['units'] = 'W'

# Gets total energy density
ddata['energy_density'] = {}
ddata['energy_density']['data'] = ff.variables['profile_fast-ion_pressure'][:] # dim(nphi, nz,nr)
ddata['energy_density']['units'] = 'J/m^3'

# Simulated time
#ddata['dt'] = 3e-3 # [s]
ddata['dt'] = 1e-4 # [s]

# Gets heat flux
ddata['heat_flux'] = {}
#ddata['heat_flux']['data'] = np.zeros_like(ddata['mode_energy']['data'])
ddata['heat_flux']['data'] = ddata['energy_density']['data']/ddata['dt']
ddata['heat_flux']['units'] = 'W/m^3'

##################################################
#
#           EQDSK
#
##################################################


edr = eqtools.EqdskReader(
    gfile=in_path+'/g'+shot+'.01000',
    afile=in_path+'/a'+shot+'.01000'
    )

dedr = {}

dedr['rGrid_1d'] = edr.getRGrid()
dedr['zGrid_1d'] = edr.getZGrid()
dedr['wall_R'] = edr.getMachineCrossSection()[0]
dedr['wall_Z'] = edr.getMachineCrossSection()[1]

dedr['rGrid_2d'], dedr['zGrid_2d'] = np.meshgrid(dedr['rGrid_1d'], dedr['zGrid_1d'])

dedr['RLCFS'] = edr.getRLCFS()[0]
dedr['ZLCFS'] = edr.getZLCFS()[0]

edr_psiRZ = edr.getFluxGrid()[0]
edr_psiLCFS = edr.getFluxLCFS()[0]
edr_psi0 = edr.getFluxAxis()[0]
dedr['rhop2D'] = np.sqrt(-1*(edr_psiRZ+edr_psi0)/(edr_psiLCFS-edr_psi0))

# Prepare interpolating 1D rho grid to 2D (R,Z) grid
xx = np.linspace(0,1,edr_psiRZ.shape[1])
yy = np.linspace(0,1,edr_psiRZ.shape[0])
XX, YY = np.meshgrid(xx,yy)
points = np.vstack((XX.ravel(), YY.ravel())).T
values = dedr['rhop2D'].ravel()

# Exclude points outside the Convex hull of the LCFS
LCFS = np.c_[dedr['RLCFS'], dedr['ZLCFS']]
hull = ConvexHull(LCFS)
hull_delaunay = Delaunay(LCFS[hull.vertices])

# Limits
dedr['Rlim'] = [0.9*np.min(dedr['RLCFS']), 1.1*np.max(dedr['RLCFS'])]
dedr['Zlim'] = [1.1*np.min(dedr['ZLCFS']), 1.1*np.max(dedr['ZLCFS'])]




#### --- Interpolating SPIRAL (R,Z) data onto rhop --- ####

ddata['rhop2D'] = {}
ddata['rhop2D']['data'] = griddata(
    (dedr['rGrid_2d'].flatten(), dedr['zGrid_2d'].flatten()),
    dedr['rhop2D'].flatten(),
    (ddata['r_2d']['data'].T.flatten(), ddata['z_2d']['data'].T.flatten()),
    method = 'linear',
    fill_value = np.nan
    ).reshape(ddata['r_2d']['data'].T.shape).T




ddata['rhop1D'] = {}
ddata['rhop1D']['data'] = np.linspace(0,1,51)

ddata['density_1d'] = {}
ddata['density_1d']['data'] = np.zeros_like(ddata['rhop1D']['data'])
ddata['T_eff_1d'] = {}
ddata['T_eff_1d']['data'] = np.zeros_like(ddata['rhop1D']['data'])
ddata['heat_flux_1d'] = {}
ddata['heat_flux_1d']['data'] = np.zeros_like(ddata['rhop1D']['data'])
ddata['pwr_tot'] = {}
ddata['pwr_tot']['data'] = 0


vol_rhop = np.zeros_like(ddata['rhop1D']['data'])


# Loop over R
for rr in np.arange(ddata['density']['data'].shape[-1]):
    # Loop over Z
    for zz in np.arange(ddata['density']['data'].shape[-2]):

        # Index of rhop bin
        ind = np.where(
            ddata['rhop2D']['data'][rr,zz] >= ddata['rhop1D']['data']
            )[0][-1]

        # Voxel volume
        voxel = (
            np.pi *(
                ddata['r_edges']['data'][rr+1]**2
                - ddata['r_edges']['data'][rr]**2
                )
            * (
                ddata['z_edges']['data'][zz+1]
                -ddata['z_edges']['data'][zz]
                )
            ) # [m3]


        # Volume integrate within rhop bin
        vol_rhop[ind] += voxel
        ddata['density_1d']['data'][ind] += ddata['density']['data'][0,zz,rr]*voxel
        ddata['T_eff_1d']['data'][ind] += ddata['T_eff']['data'][0,zz,rr]*voxel

        #ddata['heat_flux']['data'][0,zz,rr] = ddata['mode_energy']['data'][0,zz,rr]/voxel # [W/m^3]
        ddata['heat_flux_1d']['data'][ind] += ddata['heat_flux']['data'][0,zz,rr]*voxel

        ddata['pwr_tot']['data'] += ddata['mode_energy']['data'][0,zz,rr]/1e6 # [MW]

ddata['density_1d']['data'] /= vol_rhop # [1/m^3]
ddata['T_eff_1d']['data'] /= vol_rhop # [keV]
ddata['heat_flux_1d']['data'] /= vol_rhop # [W/m^3]

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

# Interpolates 1D profile to 2D
nz_2d = griddata(
    dimp['rhop_fm'],
    dimp['nz_fm'][:, 16, -1],
    values,
    method = 'linear'
    ).reshape(edr_psiRZ.shape)

# Cuts outside LCFS
for ii in np.arange(nz_2d.shape[0]):
    for jj in np.arange(nz_2d.shape[1]):
        if not hull_delaunay.find_simplex(
            np.c_[dedr['rGrid_2d'][ii,jj], dedr['zGrid_2d'][ii,jj]]
            ) >= 0:

            nz_2d[ii,jj] = np.nan


#######################################################
#
#               TORIC
#
#######################################################

# Absorbed power
Pabs = 0.5*0.3 # [MW]

# File Management
filet = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/CMOD',
    'shots',
    shot,
    'TORIC'
    )

# Reads TORIC output
toric = toric_tools.toric_analysis(
    toric_name=filet+ '/toric_nnphi.ncdf',
    mode='ICRF',
    path=filet+'/'
    )

tsol = {}

tsol['RR']  = toric.cdf_hdl.variables['Xplasma'].data/100 + Rmaj
tsol['ZZ']  = toric.cdf_hdl.variables['Zplasma'].data/100

tsol['Re2Eplus'] = toric.cdf_hdl.variables['Re2Eplus'].data
tsol['Im2Eplus'] = toric.cdf_hdl.variables['Re2Eplus'].data

tsol['pwr_2d'] = (toric.cdf_hdl.variables['TdPwIH'].data)[:,:,-2] * Pabs *1e6 # [W/m^3]

tsol['pwr_1d'] = (toric.cdf_hdl.variables['PwIH'].data)[:,-2] * Pabs *1e6 # [W/m^3]
tsol['rhop'] = edr.roa2psinorm(
        toric.cdf_hdl.variables['Pw_abscissa'].data,
        0,
        sqrt = True,
        )

#######################################################
#
#            Plotting
#
#######################################################


######## ---- Argon density ---- ############

fig = plt.figure(figsize = (12,16))
gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 0.5])


ax = fig.add_subplot(gs[0,0])
col = ax.pcolormesh(
    dedr['rGrid_2d'],dedr['zGrid_2d'],
    nz_2d,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    #vmax = np.max(np.sum(asim['nz_fm'], axis = 1))
    vmax = np.nanmax(nz_2d.flatten())
    )
ax.set_title('input; Ar16+')
#fig.colorbar(col, ax=ax[0], label = r'$n_{Ar}$ [$1/cm^3$]')

ax.plot(dedr['RLCFS'], dedr['ZLCFS'],'r',linewidth=3)
con = ax.contour(
    dedr['rGrid_1d'],dedr['zGrid_1d'],
    dedr['rhop2D'],
    [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    #cmap = 'hsv'
    colors = 'w'
    )
ax.clabel(con, inline=True, fontsize = 8)

ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')

ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])

ax.set_aspect('equal')





ax = fig.add_subplot(gs[0,1])
col = ax.pcolormesh(
    ddata['r_2d']['data'], ddata['z_2d']['data'] ,
    ddata['density']['data'][0,:,:].T/1e6,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    #vmax = np.max(np.sum(asim['nz_fm'], axis = 1))
    vmax = np.nanmax(nz_2d.flatten())
    )
ax.set_title('SPIRAL; Ar16+')
#fig.colorbar(col, ax=ax[1], label = r'$n_{Ar}$ [$1/cm^3$]')

ax.plot(dedr['RLCFS'], dedr['ZLCFS'],'r',linewidth=3)
con = ax.contour(
    dedr['rGrid_1d'],dedr['zGrid_1d'],
    dedr['rhop2D'],
    [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    #cmap = 'hsv'
    colors = 'w'
    )
ax.clabel(con, inline=True, fontsize = 8)

#ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')

ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])

ax.set_aspect('equal')

ax = fig.add_subplot(gs[0,2])
ax.set_axis_off()
fig.colorbar(col, ax=ax, label = r'$n_{Ar}$ [$1/cm^3$]')
#plt.tight_layout()


ax = fig.add_subplot(gs[1,0])


ax.plot(
    ddata['rhop2D']['data'].flatten(),
    ddata['density']['data'].T.flatten()/1e6,
    '*',
    label = 'SPIRAL; raw',
    )


ax.step(
    ddata['rhop1D']['data'],
    ddata['density_1d']['data']/1e6,
    where='post',
    label = 'SPIRAL; binned',
    linewidth = 4,
    )

ax.plot(
    dimp['rhop_fm'],
    dimp['nz_fm'][:, 16, -1],
    label = 'input',
    linewidth = 4
    )


ax.grid('on')
ax.set_ylabel(r'$n_{Ar}$ [$1/cm^3$]')
ax.set_xlabel(r'$\sqrt{\Psi_n}$')
ax.set_xlim(0,1)

leg = ax.legend(labelcolor='linecolor')










######## ---- effective Temperature ---- ############

fig = plt.figure(figsize = (8,16))
gs = gridspec.GridSpec(2, 1)

scale = 3.2
vmax = 160/scale

ax = fig.add_subplot(gs[0,0])
col = ax.pcolormesh(
    ddata['r_2d']['data'], ddata['z_2d']['data'] ,
    ddata['T_eff']['data'][0,:,:].T/cnt.e/1e3/scale,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    #vmax = np.max(np.sum(asim['nz_fm'], axis = 1))
    vmax = vmax
    )
ax.set_title('SPIRAL; Ar16+')
fig.colorbar(col, ax=ax, label = r'$T_{eff}$ [$keV$]')

ax.plot(dedr['RLCFS'], dedr['ZLCFS'],'r',linewidth=3)
con = ax.contour(
    dedr['rGrid_1d'],dedr['zGrid_1d'],
    dedr['rhop2D'],
    [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    #cmap = 'hsv'
    colors = 'w'
    )
ax.clabel(con, inline=True, fontsize = 8)

ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')

ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])

ax.set_aspect('equal')




ax = fig.add_subplot(gs[1,0])


ax.plot(
    ddata['rhop2D']['data'].flatten(),
    ddata['T_eff']['data'].T.flatten()/cnt.e/1e3/scale,
    '*',
    label = 'SPIRAL; raw',
    )


ax.step(
    ddata['rhop1D']['data'],
    ddata['T_eff_1d']['data']/cnt.e/1e3/scale,
    where='post',
    label = 'SPIRAL; binned',
    linewidth = 4,
    )

ax.plot(
    np.sqrt(inga['polflux']/inga['polflux'][-1]),
    inga['Ti_1'],
    label = 'input',
    linewidth = 4
    )



ax.grid('on')
ax.set_ylabel(r'$T_{eff}$ [$keV$]')
ax.set_xlabel(r'$\sqrt{\Psi_n}$')
ax.set_xlim(0,1)
ax.set_ylim(1e-1,vmax)
ax.set_yscale('log')

leg = ax.legend(labelcolor='linecolor')





######## ---- ICRF Heat Flux ---- ############

fig = plt.figure(figsize = (14,16))
gs = gridspec.GridSpec(2, 4, width_ratios=[1, 0.5, 1, 0.5])


scale = ddata['pwr_tot']['data']/(Pabs*0.01)
#scale = 1

ax = fig.add_subplot(gs[0,0])
col = ax.pcolormesh(
    tsol['RR'], tsol['ZZ'],
    tsol['pwr_2d'],
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    #vmax = np.nanmax(nz_2d.flatten())
    )
ax.set_title('TORIC; Ar16+')
#fig.colorbar(col, ax=ax, label = r'$Q_{RF,Ar16+}$ [$W/m^3$]')

ax.plot(dedr['RLCFS'], dedr['ZLCFS'],'r',linewidth=3)
con = ax.contour(
    dedr['rGrid_1d'],dedr['zGrid_1d'],
    dedr['rhop2D'],
    [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    #cmap = 'hsv'
    colors = 'w'
    )
ax.clabel(con, inline=True, fontsize = 8)

ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')

ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])

ax.set_aspect('equal')


ax = fig.add_subplot(gs[0,1])
ax.set_axis_off()
fig.colorbar(col, ax=ax, label = r'$Q_{RF,Ar16+}$ [$W/m^3$]', location='left')
#ax.set_position([0.0, 0.3, 0.1, 0.4])






ax = fig.add_subplot(gs[0,2])
col = ax.pcolormesh(
    ddata['r_2d']['data'], ddata['z_2d']['data'] ,
    ddata['heat_flux']['data'][0,:,:].T/scale,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    #vmax = np.nanmax(nz_2d.flatten())
    )
ax.set_title('SPIRAL; Ar16+')
#fig.colorbar(col, ax=ax[1], label = r'$n_{Ar}$ [$1/cm^3$]')

ax.plot(dedr['RLCFS'], dedr['ZLCFS'],'r',linewidth=3)
con = ax.contour(
    dedr['rGrid_1d'],dedr['zGrid_1d'],
    dedr['rhop2D'],
    [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9, 1.0],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    #cmap = 'hsv'
    colors = 'w'
    )
ax.clabel(con, inline=True, fontsize = 8)

#ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')

ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])

ax.set_aspect('equal')

ax = fig.add_subplot(gs[0,3])
ax.set_axis_off()
cbar = fig.colorbar(col, ax=ax, label = r'$Q_{RF,Ar16+}$ [$W/m^3$]')

cbar.formatter = ticker.ScalarFormatter(useMathText=True)
cbar.formatter.set_scientific(True)
cbar.formatter.set_powerlimits((-1,1))
cbar.update_ticks()


ax = fig.add_subplot(gs[1,2])


ax.plot(
    ddata['rhop2D']['data'].flatten(),
    ddata['heat_flux']['data'].T.flatten()/scale,
    '*',
    label = 'SPIRAL; raw',
    )


ax.step(
    ddata['rhop1D']['data'],
    ddata['heat_flux_1d']['data']/scale,
    where='post',
    label = 'SPIRAL; binned',
    linewidth = 4,
    )

#ax.plot(
#    dimp['rhop_fm'],
#    dimp['nz_fm'][:, 16, -1],
#    label = 'input',
#    linewidth = 4
#    )


ax.grid('on')
ax.set_ylabel(r'$Q_{Rf,Ar16+}$ [$W/m^3$]')
ax.set_xlabel(r'$\sqrt{\Psi_n}$')
ax.set_xlim(0,1)
ax.set_ylim(-1e8/scale,1.5e9/scale)

leg = ax.legend(labelcolor='linecolor')

ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.get_major_formatter().set_scientific(True)
ax.yaxis.get_major_formatter().set_powerlimits((-1, 1))


ax = fig.add_subplot(gs[1,0])

ax.plot(
    tsol['rhop'],
    tsol['pwr_1d'],
    label = 'TORIC',
    linewidth = 4,
    color = 'g'
    )


ax.grid('on')
ax.set_ylabel(r'$Q_{Rf,Ar16+}$ [$W/m^3$]')
ax.set_xlabel(r'$\sqrt{\Psi_n}$')
ax.set_xlim(0,1)


leg = ax.legend(labelcolor='linecolor')

ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
ax.yaxis.get_major_formatter().set_scientific(True)
ax.yaxis.get_major_formatter().set_powerlimits((-1, 1))


