'''

Script to learn to process SPIRAL ensemble output

cjperks
Dec 4th, 2024

'''

# Modules
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import copy
import scipy.constants as cnt

from transport_world.run_SPIRAL import output_utils as utz
from transport_world.run_profiletools import eqtools3 as eq
from transport_world.make_gacode import exp_utils as utils

# File management
file_in = prf = os.path.join(
    '/home/cjperks/work',
    '2201_Pumpout/SPARC',
    'SPIRAL_He3'
    )
prf = os.path.join(
    file_in,
    #'test_sparc_prf_6717797.cdf'
    'test_sparc_prf_6722649.cdf'
    )
ens = os.path.join(
    file_in, 
    #'test_sparc_6716555.txt'
    'test_sparc_6721918.txt'
    )
nmark = 1e4

# Gets ACTION_profiles data
dprf = utz.read_prf(prf=prf)

# Gets ACTION_thermal data
dens = utz.read_ens(
    nmark=nmark,
    ens=ens,
    )

##################################################
#
#           EQDSK
#
##################################################

# Paths
in_path = os.path.join(
    '/home/cjperks',
    'tofu_sparc/background_plasma',
    'PRD_plasma/run1'
    )

dedr, edr = eq._get_eq(
    gfile = in_path+'/input.geq',
    afile = None,
    #afile = in_path+'/workAround.aeq',
    machine = 'SPARC'
    )

# SPARC geom
Rmaj = 189.553/100# [m]
Zmax = [-100/100,100/100] # [m]

Rout, Zout = utils.get_harm_res(
    path_input = in_path+'/',
    path_gfile = 'input.geq',
    freq = 120e6,
    order = 1,
    cs = 2,
    amu = 3,
    )

# Mask particles outside the LCFS
mask_st = np.where(
    dens['PSI_START'] <= abs(dedr['psiLCFS'])
    )
mask_end = np.where(
    dens['PSI_END'] <= abs(dedr['psiLCFS'])
    )



##################################################
#
#           Plotting
#
##################################################

##### ---- Marker histogram in (R,Z) ---- #####

# Makes a 2D histogram of marker positions in (R,Z)
nbins_st = 41
nbins_end = 81
(
    hist_RZ_start, 
    R_edges_st, 
    Z_edges_st
    ) = np.histogram2d(
        dens['R_START'][mask_st], 
        dens['Z_START'][mask_st], 
        #bins=[dprf['r_edges']['data'], dprf['z_edges']['data']]
        bins = [nbins_st, nbins_end]
        )
R_cents_st = 0.5*(R_edges_st[1:] + R_edges_st[:-1])
Z_cents_st = 0.5*(Z_edges_st[1:] + Z_edges_st[:-1])
Z_cents_st_2d, R_cents_st_2d = np.meshgrid(Z_cents_st, R_cents_st)
(
    hist_RZ_end, 
    R_edges_end, 
    Z_edges_end
    ) = np.histogram2d(
        dens['R_END'][mask_end], 
        dens['Z_END'][mask_end], 
        #bins=[dprf['r_edges']['data'], dprf['z_edges']['data']]
        #bins = [101, 201]
        bins = [R_edges_st, Z_edges_st]
        )
R_cents_end = 0.5*(R_edges_end[1:] + R_edges_end[:-1])
Z_cents_end = 0.5*(Z_edges_end[1:] + Z_edges_end[:-1])
Z_cents_end_2d, R_cents_end_2d = np.meshgrid(Z_cents_end, R_cents_end)

'''
# Cuts outside LCFS
hist_RZ_end = eq._trim_hull(
    dedr = dedr,
    RR_2d = dprf['r_2d']['data'],
    ZZ_2d = dprf['z_2d']['data'],
    data_2d = hist_RZ_end
    )
'''

hist_RZ_start[np.where(hist_RZ_start == 0)] = np.nan
hist_RZ_end[np.where(hist_RZ_end == 0)] = np.nan



fig = plt.figure(figsize = (12,16))
gs = gridspec.GridSpec(1, 2)

ax = fig.add_subplot(gs[0,0])
col = ax.pcolormesh(
    R_cents_st_2d, Z_cents_st_2d,
    hist_RZ_start,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    )
ax.set_title('Start locations')

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

ax.plot(
    (Rout)[((Zout <Zmax[1]) & (Zout>Zmax[0]))],
    Zout[((Zout <Zmax[1]) & (Zout>Zmax[0]))],
    '--', 
    linewidth=3,
    color = 'r'
    )

ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')
ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])
ax.set_aspect('equal')

fig.colorbar(col, ax=ax, label = 'counts')



ax = fig.add_subplot(gs[0,1])
col = ax.pcolormesh(
    R_cents_end_2d, Z_cents_end_2d,
    hist_RZ_end,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    #vmax = vmax
    )
ax.set_title('End locations')

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

ax.plot(
    (Rout)[((Zout <Zmax[1]) & (Zout>Zmax[0]))],
    Zout[((Zout <Zmax[1]) & (Zout>Zmax[0]))],
    '--', 
    linewidth=3,
    color = 'r'
    )

#ax.set_ylabel('Z [m]')
ax.set_xlabel('R [m]')
ax.set_xlim(dedr['Rlim'])
ax.set_ylim(dedr['Zlim'])
ax.set_aspect('equal')

fig.colorbar(col, ax=ax, label = 'counts')




'''
vol = np.zeros_like(hist_RZ_start)
# Loop over R
for rr in np.arange(dprf['density']['data'].shape[-1]):
    # Loop over Z
    for zz in np.arange(dprf['density']['data'].shape[-2]):

        # Voxel volume
        vol[rr,zz] = (
            np.pi *(
                dprf['r_edges']['data'][rr+1]**2
                - dprf['r_edges']['data'][rr]**2
                )
            * (
                dprf['z_edges']['data'][zz+1]
                -dprf['z_edges']['data'][zz]
                )
            ) # [m3]

tmp = hist_RZ_end/vol # [cnt/m3]
'''


##### ---- Marker histogram in (vpar,vperp) ---- #####

dens_fil = copy.deepcopy(dens)


MM = 3 *cnt.m_u
vpe_1d = np.linspace(0,5e6, 1001)
vpa_1d = np.linspace(-5e6, 5e6, 2001)
vpe_2d, vpa_2d = np.meshgrid(vpe_1d,vpa_1d)
energy = 0.5*MM*(vpa_2d**2+vpe_2d**2)/cnt.e/1e3 # [keV]

# Makes a 2D histogram of marker positions in (vpar,vperp)
nbin_vpar = 101
nbin_vperp = 101
(
    hist_VparVperp_start, 
    Vpar_edges_st, 
    Vperp_edges_st
    ) = np.histogram2d(
        dens_fil['VPAR_START'][mask_st], 
        dens_fil['VPERP_START'][mask_st], 
        bins=[nbin_vpar, nbin_vperp]
        )
Vpar_cents_st = 0.5*(Vpar_edges_st[1:] + Vpar_edges_st[:-1])
Vperp_cents_st = 0.5*(Vperp_edges_st[1:] + Vperp_edges_st[:-1])
Vperp_cents_st_2d, Vpar_cents_st_2d = np.meshgrid(Vperp_cents_st, Vpar_cents_st)
(
    hist_VparVperp_end, 
    Vpar_edges_end, 
    Vperp_edges_end
    ) = np.histogram2d(
        dens_fil['VPAR_END'][mask_end], 
        dens_fil['VPERP_END'][mask_end], 
        bins=[nbin_vpar, nbin_vperp],
        #bins = [Vpar_edges_st, Vperp_edges_st]
        range = [[-2e6, 2e6], [0, 1e7]]
        )
Vpar_cents_end = 0.5*(Vpar_edges_end[1:] + Vpar_edges_end[:-1])
Vperp_cents_end = 0.5*(Vperp_edges_end[1:] + Vperp_edges_end[:-1])
Vperp_cents_end_2d, Vpar_cents_end_2d = np.meshgrid(Vperp_cents_end, Vpar_cents_end)

hist_VparVperp_start[np.where(hist_VparVperp_start == 0)] = np.nan
hist_VparVperp_end[np.where(hist_VparVperp_end == 0)] = np.nan



fig = plt.figure(figsize = (12,16))
gs = gridspec.GridSpec(1, 2)

ax = fig.add_subplot(gs[0,0])
col = ax.pcolormesh(
    Vpar_cents_st_2d, Vperp_cents_st_2d,
    hist_VparVperp_start,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    )
con = ax.contour(
    vpa_2d, vpe_2d,
    energy,
    [1,2,8,20,50,100],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    colors= 'r'
    )
ax.clabel(con, inline=True, fontsize = 8)

ax.set_title('Start locations')

ax.set_ylabel(r'$v_\perp$ [m/s]')
ax.set_xlabel(r'$v_{||}$ [m/s]')
ax.set_xlim(np.min(Vpar_edges_st), np.max(Vpar_edges_st))
ax.set_ylim(np.min(Vperp_edges_st), np.max(Vperp_edges_st))
ax.set_aspect('equal')

fig.colorbar(col, ax=ax, label = 'counts')



ax = fig.add_subplot(gs[0,1])
col = ax.pcolormesh(
    Vpar_cents_end_2d, Vperp_cents_end_2d,
    hist_VparVperp_end,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    )
con = ax.contour(
    vpa_2d, vpe_2d,
    energy,
    [1,2,8,20,50,100],
    linestyles='solid',
    linewidth=2,
    zorder=2,
    colors= 'r'
    )
ax.clabel(con, inline=True, fontsize = 8)

ax.set_title('End locations')

#ax.set_ylabel(r'$v_\perp$ [m/s]')
ax.set_xlabel(r'$v_{||}$ [m/s]')
ax.set_xlim(np.min(Vpar_edges_end), np.max(Vpar_edges_end))
ax.set_ylim(np.min(Vperp_edges_end), np.max(Vperp_edges_end))
ax.set_aspect('equal')

fig.colorbar(col, ax=ax, label = 'counts')





##### ---- Marker histogram in (energy,pitch) ---- #####


# Makes a 2D histogram of marker positions in (energy,pitch)
nbin_energy = 101
nbin_pitch = 41
(
    hist_energyPitch_start, 
    energy_edges_st, 
    pitch_edges_st
    ) = np.histogram2d(
        dens_fil['ENERGY_START'][mask_st], 
        dens_fil['PITCH_START'][mask_st], 
        bins=[nbin_energy, nbin_pitch]
        )
energy_cents_st = 0.5*(energy_edges_st[1:] + energy_edges_st[:-1])
pitch_cents_st = 0.5*(pitch_edges_st[1:] + pitch_edges_st[:-1])
pitch_cents_st_2d, energy_cents_st_2d = np.meshgrid(pitch_cents_st, energy_cents_st)
(
    hist_energyPitch_end, 
    energy_edges_end, 
    pitch_edges_end
    ) = np.histogram2d(
        dens_fil['ENERGY_END'][mask_end], 
        dens_fil['PITCH_END'][mask_end], 
        bins=[nbin_energy, nbin_pitch],
        #bins = [energy_edges_st, pitch_edges_st]
        range = [[0, 1e2], [-1,1]]
        )
energy_cents_end = 0.5*(energy_edges_end[1:] + energy_edges_end[:-1])
pitch_cents_end = 0.5*(pitch_edges_end[1:] + pitch_edges_end[:-1])
pitch_cents_end_2d, energy_cents_end_2d = np.meshgrid(pitch_cents_end, energy_cents_end)


hist_energyPitch_start[np.where(hist_energyPitch_start == 0)] = np.nan
hist_energyPitch_end[np.where(hist_energyPitch_end == 0)] = np.nan



fig = plt.figure(figsize = (12,16))
gs = gridspec.GridSpec(1, 2)

ax = fig.add_subplot(gs[0,0])
col = ax.pcolormesh(
    pitch_cents_st_2d, energy_cents_st_2d,
    hist_energyPitch_start,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    )

ax.set_title('Start locations')

ax.set_ylabel(r'Energy [keV]')
ax.set_xlabel(r'$v_{||}/v$ []')
ax.set_ylim(np.min(energy_edges_st), np.max(energy_edges_st))
ax.set_xlim(np.min(pitch_edges_st), np.max(pitch_edges_st))
#ax.set_aspect('equal')

fig.colorbar(col, ax=ax, label = 'counts')



ax = fig.add_subplot(gs[0,1])
col = ax.pcolormesh(
    pitch_cents_end_2d, energy_cents_end_2d,
    hist_energyPitch_end,
    shading='auto',
    cmap = 'viridis',
    vmin = 0,
    )

ax.set_title('End locations')

#ax.set_ylabel(r'$v_\perp$ [m/s]')
ax.set_xlabel(r'$v_{||}/v$ []')
ax.set_ylim(np.min(energy_edges_end), np.max(energy_edges_end))
ax.set_xlim(np.min(pitch_edges_end), np.max(pitch_edges_end))
#ax.set_aspect('equal')

fig.colorbar(col, ax=ax, label = 'counts')










###### ---- distribution of energies --- #####

Ti = 11 # [keV]
enrg = np.logspace(np.log10(1e-1), np.log10(1e2), 1001)
IEDF = (
    2*np.sqrt(enrg/np.pi)
    * Ti**(-3/2)
    * np.exp(-enrg/Ti)
    ) # [1/keV]


nbins = 101
bin_edges = np.logspace(np.log10(1e-1), np.log10(1e2), 101)
(
    hist_E_st, bins_E_st
    ) = np.histogram(
        dens['ENERGY_START'][mask_st],
        #bins = nbins
        bins = bin_edges
        )
bin_edges = np.logspace(np.log10(1e0), np.log10(100e3), 101)
(
    hist_E_end, bins_E_end
    ) = np.histogram(
        dens['ENERGY_END'][mask_end],
        #bins = nbins,
        #range = [0, 4e3]
        bins = bin_edges
        )




fig = plt.figure(figsize = (12,16))
gs = gridspec.GridSpec(1, 2)

ax = fig.add_subplot(gs[0,0])
ax.step(
    bins_E_st[:-1],
    hist_E_st/np.diff(bins_E_st),
    where='post',
    label='SPIRAL'
    )

scale1 = 0.8*(np.max((hist_E_st/np.diff(bins_E_st)).flatten())/np.max(IEDF))
ax.plot(
    enrg,
    IEDF* scale1,
    'r-',
    label='Max; <Ti>=%i keV'%(Ti)
    )

ax.set_title('Start locations')
ax.set_xlabel('Energy [keV]')
ax.set_ylabel('counts/keV')
ax.grid('on')
ax.set_xscale('log')

leg = ax.legend(labelcolor='linecolor')
leg.set_draggable('on')

ax = fig.add_subplot(gs[0,1])
ax.step(
    bins_E_end[:-1],
    hist_E_end/np.diff(bins_E_end),
    where='post'
    )

scale2 = 0.8*(np.max((hist_E_end/np.diff(bins_E_end)).flatten())/np.max(IEDF))
ax.plot(
    enrg,
    IEDF* scale2,
    'r-'
    )

ax.set_title('End locations')
ax.set_xlabel('Energy [keV]')
ax.set_ylabel('counts/keV')
ax.grid('on')
ax.set_xscale('log')






###### ---- distribution of pitch --- #####

nbins = 101
(
    hist_p_st, bins_p_st
    ) = np.histogram(
        dens['PITCH_START'][mask_st],
        bins = nbins
        #bins = bin_edges
        )
(
    hist_p_end, bins_p_end
    ) = np.histogram(
        dens['PITCH_END'][mask_end],
        bins = nbins,
        )


fig = plt.figure(figsize = (12,16))
gs = gridspec.GridSpec(1, 2)

ax = fig.add_subplot(gs[0,0])
ax.step(
    bins_p_st[:-1],
    hist_p_st,
    where='post'
    )

ax.set_title('Start locations')
ax.set_xlabel(r'$v_{||}/v$ []')
ax.set_ylabel('counts')
ax.grid('on')


ax = fig.add_subplot(gs[0,1])
ax.step(
    bins_p_end[:-1],
    hist_p_end,
    where='post'
    )

ax.set_title('End locations')
ax.set_xlabel(r'$v_{||}/v$ []')
ax.set_ylabel('counts')
ax.grid('on')
