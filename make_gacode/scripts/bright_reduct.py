'''

Calculates the brightness reductions on HIREX data
accounting for changes in ne/Te

cjperks
Sept 21, 2023

'''

# Modules
from transport_world.make_gacode import read_tokamak as rTok
import aurora
from omfit_classes import omfit_gapy
from scipy.interpolate import interp1d
import os

# User-controls
shot = 1140221013
tht = 0

shift_lya1 = 0.1

path_input = '/home/cjperks/2201_Pumpout/CMOD/shots/'+str(shot)
fgacode = [
    'input_t0.8s.gacode',
    'input_t1.2s.gacode'
    ]

dtimes = {
    't1': [0.75, 0.85],
    't2': [1.15, 1.25],
    }

# Obtain HIREXSR data
ddata = rTok.get_hirexsr(
    shot=shot,
    tht=tht,
    quants = ['moments'],
    plt_all = False,
    )

fig3, ax3 = plt.subplots(2, len(ddata['mom'].keys()))
fig3.tight_layout(pad=1.0)

# Collects rate coeffiecnt data for impurity ion, (scd-> ionization rates, acd-> recombination rates, ccd->charge exchange rates) 
atom_data = aurora.atomic.get_atom_data('Ar',['scd','acd'])

# Loop over lines
for ll, line in enumerate(list(ddata['mom'].keys())):
    # Time window to average data
    ind1 = np.where(
        (ddata['mom'][line]['t_s'] >= dtimes['t1'][0])
        & (ddata['mom'][line]['t_s'] <= dtimes['t1'][-1])
        )[0]
    ind2 = np.where(
        (ddata['mom'][line]['t_s'] >= dtimes['t2'][0])
        & (ddata['mom'][line]['t_s'] <= dtimes['t2'][-1])
        )[0]

    # Counts data
    data1 = np.mean(ddata['mom'][line]['data'][0,ind1,:], axis=-2) # [], dim(nchan,)
    data2 = np.mean(ddata['mom'][line]['data'][0,ind2,:], axis=-2) # [], dim(nchan,)

    sqpsin1 = np.sqrt(np.mean(ddata['mom'][line]['psin'][ind1,:], axis=-2)) # [], dim(nchan,)
    sqpsin2 = np.sqrt(np.mean(ddata['mom'][line]['psin'][ind2,:], axis=-2) )# [], dim(nchan,)

    if line == 'LYA1':
        sqpsin1_tmp = np.zeros(len(sqpsin1))
        sqpsin2_tmp = np.zeros(len(sqpsin2))

        # Loop over channel
        for ii in np.arange(len(sqpsin1)):
            if ii == len(sqpsin1)-1:
                sqpsin1_tmp[ii] = abs(
                    sqpsin1[ii]
                    + shift_lya1*np.sign(sqpsin1[ii]-sqpsin1[ii-1])
                    )
                sqpsin2_tmp[ii] = abs(
                    sqpsin2[ii]
                    + shift_lya1*np.sign(sqpsin2[ii]-sqpsin2[ii-1])
                    )
            else:
                sqpsin1_tmp[ii] = abs(
                    sqpsin1[ii]
                    + shift_lya1*np.sign(sqpsin1[ii+1]-sqpsin1[ii])
                    )
                sqpsin2_tmp[ii] = abs(
                    sqpsin2[ii]
                    + shift_lya1*np.sign(sqpsin2[ii+1]-sqpsin2[ii])
                    )

        sqpsin1 = sqpsin1_tmp
        sqpsin2 = sqpsin2_tmp


    ax3[0,ll].plot(
        sqpsin1,
        1- (data2/data1),
        'r*'
        )

    ax3[0,ll].set_xlabel(r'$\sqrt{\Psi_n}$')
    ax3[0,ll].set_ylabel('Brightness reduction (t2/t1)')
    ax3[0,ll].set_title(line)
    ax3[0,ll].set_ylim(0.7,1)
    ax3[0,ll].set_xlim(0,1)
    ax3[0,ll].grid('on')

    # Loads input.gacode profiles
    ga1 = omfit_gapy.OMFITgacode(
        os.path.join(
            path_input,
            fgacode[0],
            )
        )
    ga2 = omfit_gapy.OMFITgacode(
        os.path.join(
            path_input,
            fgacode[1],
            )
        )

    # Calculates fractional abundances assuming ionization equilibrium 
    _, fz1 = aurora.atomic.get_frac_abundances(
        atom_data, 
        ga1['ne']*1e13, ga1['Te']*1e3, 
        plot=False
        ) # dim(rhop, charge)
    _, fz2 = aurora.atomic.get_frac_abundances(
        atom_data, 
        ga2['ne']*1e13, ga2['Te']*1e3, 
        plot=False
        ) # dim(rhop, charge)

    if line == 'X':
        cs = 16
        lam = 3.96581
        adf15 = '/home/cjperks/tofu_sparc/atomic_data/ADAS_PEC_files/fsciortino_share/atomdb/pec#ar16.dat'
    elif line == 'Z':
        cs = 16
        lam = 3.99417
        adf15 = '/home/cjperks/tofu_sparc/atomic_data/ADAS_PEC_files/fsciortino_share/atomdb/pec#ar16.dat'
    elif line == 'LYA1':
        cs = 17
        lam = 3.73114
        adf15 = '/home/cjperks/tofu_sparc/atomic_data/ADAS_PEC_files/fsciortino_share/atomdb/pec#ar17.dat'
    dlam = 0.001

    # Loads all lines
    trs = aurora.read_adf15(adf15)

    # Loads specific line
    sub_trs = trs.loc[
        (trs['lambda [A]'] > lam - dlam)
        & (trs['lambda [A]'] < lam + dlam)
        ]

    # Calculates PEC assuming just excitation
    sel_trs = sub_trs.loc[(sub_trs['type'] == 'excit')]

    pec1 = sel_trs['log10 PEC fun'].iloc[0].ev(np.log10(ga1['ne']*1e13), np.log10(ga1['Te']*1e3)) # dim(rhop,)
    pec2 = sel_trs['log10 PEC fun'].iloc[0].ev(np.log10(ga2['ne']*1e13), np.log10(ga2['Te']*1e3)) # dim(rhop,)

    # Interpolates profiles
    ne1 = interp1d(
        ga1['polflux']/ga1['polflux'][-1],
        ga1['ne']
        )(sqpsin1**2)
    ne2 = interp1d(
        ga2['polflux']/ga2['polflux'][-1],
        ga2['ne']
        )(sqpsin2**2)

    fz1 = interp1d(
        ga1['polflux']/ga1['polflux'][-1],
        fz1[:,cs]
        )(sqpsin1**2)
    fz2 = interp1d(
        ga2['polflux']/ga2['polflux'][-1],
        fz2[:,cs]
        )(sqpsin2**2)

    pec1 = interp1d(
        ga1['polflux']/ga1['polflux'][-1],
        pec1
        )(sqpsin1**2)
    pec2 = interp1d(
        ga2['polflux']/ga2['polflux'][-1],
        pec2
        )(sqpsin2**2)

    ax3[1,ll].plot(
        sqpsin1,
        1 - (data2/data1)* (ne1*fz1*pec1)/(ne2*fz2*pec2),
        'b*'
        )

    ax3[1,ll].set_xlabel(r'$\sqrt{\Psi_n}$')
    ax3[1,ll].set_ylabel('proxy Ar density reduction (t2/t1)')
    #ax3[0,ll].set_title(line)
    ax3[1,ll].set_ylim(0.8,1)
    ax3[1,ll].set_xlim(0,1)
    ax3[1,ll].grid('on')