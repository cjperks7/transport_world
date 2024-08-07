'''

Loads NetCDF data files generated by profiletools/gpfit

cjperks
June 28, 2024

'''

# Module
from transport_world.run_profiletools.profiles2 import read_gpfit_out as rg
import sys, os
from transport_world.make_gacode import read_tokamak as rt
sys.path.insert(0,'/home/cjperks/usr/python3modules/eqtools3')
import eqtools
sys.path.pop(0)

# Enables automatic reloading of modules
%reload_ext autoreload
%autoreload 2

# Shot of interest
shot = 1140221013

tmin = 0.86
tmax = 0.95

####################################################
#
#               Fits
#
####################################################

# File management
fol = os.path.join(
    '/home/cjperks',
    'test'
    )
fil = {
    'Ti': 'test',
    'Te': 'test2',
    'ne': 'test3'
    }

# Loop over kinetic profiles
profs = {}
for prof in fil.keys():
    profs[prof] = rg.read(file=os.path.join(fol, fil[prof]))

'''
# Writes fits to a convenient asciit
rg.write2ascii(
    fol=fol,
    profs=profs
    )
'''
####################################################
#
#               Exp Data
#
####################################################

# Gets data
dexp = rt.profs_cmod(shot=shot)

# Gets eqdsk
edr = eqtools.EqdskReader(
    gfile=os.path.join(
        '/home/cjperks/work',
        '2201_Pumpout/CMOD/shots',
        str(shot),
        'g'+str(shot)+'.01000'
        ),
    afile=None
    )

####################################################
#
#               PLotting
#
####################################################

fig, ax = plt.subplots(1,2)

fig.suptitle((
    str(profs['ne']['shot'])
    + '; t_min =%0.4fs, t_max =%0.4fs'%(tmin,tmax)
    ))

ax[0].plot(
    profs['ne']['sqrt{psi_n}']['data'],
    profs['ne']['n_e, CTS']['data'],
    color = 'blue',
    label = 'gpfit'
    )
ax[0].fill_between(
    profs['ne']['sqrt{psi_n}']['data'],
    (
        profs['ne']['n_e, CTS']['data'] 
        - profs['ne']['err_n_e, CTS']['data']
        ),
    (
        profs['ne']['n_e, CTS']['data'] 
        + profs['ne']['err_n_e, CTS']['data']
        ),
    alpha = 0.6,
    color = 'blue'
    )
ax[0].fill_between(
    profs['ne']['sqrt{psi_n}']['data'],
    (
        profs['ne']['n_e, CTS']['data'] 
        - 2*profs['ne']['err_n_e, CTS']['data']
        ),
    (
        profs['ne']['n_e, CTS']['data'] 
        + 2*profs['ne']['err_n_e, CTS']['data']
        ),
    alpha = 0.3,
    color = 'blue'
    )


# Loop over diagnostics
#colormap = plt.get_cmap('tab10')
#colors = colormap(np.linspace(0,1,len(dexp['ne'].keys())))
colors = ['m', 'g', 'c', 'y', 'k', 'r']

for dd, diag in enumerate(dexp['ne'].keys()):
    if len(dexp['ne'][diag].keys()) < 1:
        continue

    # Time window of interest
    indt = np.where(
        (dexp['ne'][diag]['t_s'] >= tmin)
        & (dexp['ne'][diag]['t_s'] <= tmax)
        )[0]

    if diag == 'TCI':
        ax[0].plot(
            [0,1],
            [np.mean(dexp['ne'][diag]['val_1e20m3'][indt])/1e20] *2,
            '--',
            color=colors[dd],
            label = diag
            )


    else:
        # Convert R_mid axis to sqrt(psi_n)
        rhop_tmp = np.zeros(dexp['ne'][diag]['r_m'][:,indt].shape)
        for jj,ii in enumerate(indt):
            rhop_tmp[:,jj] = np.sqrt(edr.rmid2psinorm(
                list(dexp['ne'][diag]['r_m'][:,ii]),
                dexp['ne'][diag]['t_s'][ii]
                ))

        # Takes weighted average
        weights = 1/dexp['ne'][diag]['err_1e20m3'][:,indt]**2
        weights[np.isinf(weights)] = 0
        weights[np.isnan(weights)] = 0
        rhop_avg = (rhop_tmp*weights).sum(axis=1)/weights.sum(axis=1)
        rhop_uncert = np.average((rhop_tmp - rhop_avg[:,None])**2, weights=weights, axis=1)
        data_avg = (dexp['ne'][diag]['val_1e20m3'][:,indt] *weights).sum(axis=1)/weights.sum(axis=1)
        uncert_avg = np.sqrt(1/weights.sum(axis=1))

        ax[0].errorbar(
            rhop_avg,
            data_avg/1e20,
            yerr=uncert_avg/1e20,
            xerr=rhop_uncert,
            fmt='.',
            color = colors[dd],
            ecolor = colors[dd],
            label = diag
            )

ax[0].set_xlabel(r'$\sqrt{\Psi_n}$')
ax[0].set_ylabel(r'$n_e$ ['+profs['ne']['n_e, CTS']['units']+']')

ax[0].set_xlim(0,1)
ax[0].grid('on')

leg0 = ax[0].legend()#labelcolor='linecolor')
leg0.set_draggable('on')




ax[1].plot(
    profs['Te']['sqrt{psi_n}']['data'],
    profs['Te']['T_e, CTS']['data'],
    color = 'blue',
    label = r'$T_e,$ gpfit'
    )
ax[1].plot(
    profs['Ti']['sqrt{psi_n}']['data'],
    profs['Ti']['T_i']['data'],
    color = 'red',
    label = r'$T_i,$ gpfit'
    )
ax[1].fill_between(
    profs['Te']['sqrt{psi_n}']['data'],
    (
        profs['Te']['T_e, CTS']['data'] 
        - profs['Te']['err_T_e, CTS']['data']
        ),
    (
        profs['Te']['T_e, CTS']['data'] 
        + profs['Te']['err_T_e, CTS']['data']
        ),
    alpha = 0.6,
    color = 'blue'
    )
ax[1].fill_between(
    profs['Te']['sqrt{psi_n}']['data'],
    (
        profs['Te']['T_e, CTS']['data'] 
        - 2*profs['Te']['err_T_e, CTS']['data']
        ),
    (
        profs['Te']['T_e, CTS']['data'] 
        + 2*profs['Te']['err_T_e, CTS']['data']
        ),
    alpha = 0.3,
    color = 'blue'
    )
ax[1].fill_between(
    profs['Ti']['sqrt{psi_n}']['data'],
    (
        profs['Ti']['T_i']['data'] 
        - profs['Ti']['err_T_i']['data']
        ),
    (
        profs['Ti']['T_i']['data'] 
        + profs['Ti']['err_T_i']['data']
        ),
    alpha = 0.6,
    color = 'red'
    )
ax[1].fill_between(
    profs['Ti']['sqrt{psi_n}']['data'],
    (
        profs['Ti']['T_i']['data'] 
        - 2*profs['Ti']['err_T_i']['data']
        ),
    (
        profs['Ti']['T_i']['data'] 
        + 2*profs['Ti']['err_T_i']['data']
        ),
    alpha = 0.3,
    color = 'red'
    )



# Loop over diagnostics
isel = -1
for xx in ['Te', 'Ti']:
    for dd, diag in enumerate(dexp[xx].keys()):
        if len(dexp[xx][diag].keys()) < 1:
            continue
        if diag in ['GPC', 'GPC2', 'FRC']:
            continue
        isel += 1

        # Time window of interest
        indt = np.where(
            (dexp[xx][diag]['t_s'] >= tmin)
            & (dexp[xx][diag]['t_s'] <= tmax)
            )[0]

        if 'r_m' in dexp[xx][diag].keys():
            # Convert R_mid axis to sqrt(psi_n)
            rhop_tmp = np.zeros(dexp[xx][diag]['r_m'][:,indt].shape)
            for jj,ii in enumerate(indt):
                rhop_tmp[:,jj] = np.sqrt(edr.rmid2psinorm(
                    list(dexp[xx][diag]['r_m'][:,ii]),
                    dexp[xx][diag]['t_s'][ii]
                    ))
        elif 'psin' in dexp[xx][diag].keys():
            rhop_tmp = np.sqrt(dexp[xx][diag]['psin'][:,indt])

        if 'err_keV' not in dexp[xx][diag].keys():
            dexp[xx][diag]['err_keV'] = 0.1*dexp[xx][diag]['val_keV']

        # Takes weighted average
        weights = 1/dexp[xx][diag]['err_keV'][:,indt]**2
        weights[np.isinf(weights)] = 0
        rhop_avg = (rhop_tmp*weights).sum(axis=1)/weights.sum(axis=1)
        rhop_uncert = np.average((rhop_tmp - rhop_avg[:,None])**2, weights=weights, axis=1)
        data_avg = (dexp[xx][diag]['val_keV'][:,indt] *weights).sum(axis=1)/weights.sum(axis=1)
        uncert_avg = np.sqrt(1/weights.sum(axis=1))

        ax[1].errorbar(
            rhop_avg,
            data_avg,
            yerr=uncert_avg,
            xerr=rhop_uncert,
            fmt='.',
            color = colors[isel],
            ecolor = colors[isel],
            label = xx +', '+diag
            )


ax[1].set_xlabel(r'$\sqrt{\Psi_n}$')
ax[1].set_ylabel(r'$T_s$ ['+profs['Te']['T_e, CTS']['units']+']')

ax[1].set_xlim(0,1)
ax[1].set_ylim(0, 1.3*np.nanmax(profs['Te']['T_e, CTS']['data']))
ax[1].grid('on')

leg = ax[1].legend()#labelcolor='linecolor')
leg.set_draggable('on')