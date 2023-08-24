'''

_calc_imp_dens is a function to calculate impurity density 
profiles using Aurora

'''

# Modules
import aurora
from scipy.interpolate import interp1d
from omfit_classes import omfit_eqdsk, omfit_gapy
import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})

__all__ = [
    'calc_imp_dens',
    ]

###########################################
#
#           Main
#
###########################################

def calc_imp_dens(
    dmodel = None,
    plt_all = None,
    plt_cs = None,
    ):

    # Reads in magnetic equilibrium and kinetic profiles
    geqdsk = omfit_eqdsk.OMFITgeqdsk(
        os.path.join(
            dmodel['paths']['in_path'],
            dmodel['paths']['geqdsk']
            )
        )
    inputgacode = omfit_gapy.OMFITgacode(
        os.path.join(
            dmodel['paths']['in_path'],
            dmodel['paths']['fgacode'],
            )
        )

    # Initializes default AURORA namelist
    nml = aurora.load_default_namelist()

    # Defining impurity of interest
    nml['imp'] = dmodel['imp']['sym']

    # Defining impurity edge source flux
    nml['source_type'] = dmodel['AURORA']['source']['type']
    nml['source_rate'] = dmodel['AURORA']['source']['rate'] # [particles/s]

    # Limiter surface (R,Z) contour
    zlim = geqdsk['fluxsurfaces']['info']['zlim']
    rlim = geqdsk['fluxsurfaces']['info']['rlim']

    # kinetic profiles
    kp = nml['kin_profs']
    kp['Te']['rhop'] = kp['ne']['rhop'] = kp['Ti']['rhop'] = np.sqrt(inputgacode['polflux']/inputgacode['polflux'][-1])    # rho_pol coordinates from GACODE
    kp['ne']['vals'] = inputgacode['ne']*1e13 # [cm^-3]; dim(rhop,); electron denisty
    kp['Te']['vals'] = inputgacode['Te']*1e3 # [eV]; dim(rhop,); electron temperature
    kp['Ti']['vals'] = inputgacode['Ti_1']*1e3 # [eV]; dim(rhop,); ion temperature
    nml['device'] = dmodel['AURORA']['device']

    # -------------------
    # SOL modeling
    # -------------------

    # SOL decay length
    kp['ne']['fun'] = 'interp'

    kp['Te']['decay'] = dmodel['AURORA']['SOL']['decay']
    kp['Ti']['decay'] = dmodel['AURORA']['SOL']['decay']
    kp['ne']['decay'] = dmodel['AURORA']['SOL']['decay']

    # -------------------
    # Edge modeling
    # -------------------

    # Source location wrt LCFS
    nml['source_cm_out_lcfs'] = dmodel['AURORA']['edge']['source_cm_out_lcfs']

    # Recyling
    nml['recycling_flag'] = dmodel['AURORA']['edge']['recycling_flag'] # recycling from divertor
    nml['wall_recycling'] = dmodel['AURORA']['edge']['wall_recycling'] # recycling from limiter

    # Connection lengths
    nml['clen_divertor'], nml['cln_limiter'] = aurora.grids_utils.estimate_clen(geqdsk)

    # Edge time scales
    #nml['tau_div_SOL_ms'] # enter core from divertor through SOL
    #nml['tau_pump_ms'] # remove particles from sim via pump, default very long
    #nml['tau_rcl_ret_ms'] # time scale for recylcling retention via wall

    # SOL mach number
    #nml['SOL_mach'] # used to compute parallel loss rate

    # -------------------
    # Time Domain controls
    # -------------------

    # Start and end time
    nml['timing']['times'] = dmodel['AURORA']['times']

    # Intializes Aurora simulation
    asim = aurora.aurora_sim(nml, geqdsk=geqdsk)

    # Main ion mass, [amu]
    asim.main_ion_A = dmodel['AURORA']['main_ion_A']

    # note that to be able to give charge-dependent Dz and Vz,
    # one has to give also time dependence (see documentation of aurora.core.run_aurora()),
    # so a dummy time dimension with only one element is introduced here
    times_DV = np.array([0])

    # Impurity density profile initial condition
    nz_init  = np.zeros((asim.rvol_grid.size, asim.Z_imp+1)) # dim(fm_rhop, Z); [1/cm3]

    # Obtains D, V profiles given desired modeling
    if 'FACIT' in dmodel['options']:
        DZ, VZ, asym, TSC = _get_DV(
            dmodel = dmodel,
            asim = asim,
            geqdsk = geqdsk,
            inputgacode = inputgacode,
            times_DV = times_DV,
            nz_init = nz_init,
            ) # [cm2/s, cm/s]; dim(fm_rhop, t, Z)
    else:
        DZ, VZ = _get_DV(
            dmodel = dmodel,
            asim = asim,
            geqdsk = geqdsk,
            inputgacode = inputgacode,
            times_DV = times_DV,
            nz_init = nz_init,
            ) # [cm2/s, cm/s]; dim(fm_rhop, t, Z)

    # run Aurora forward model and plot results
    out = asim.run_aurora(
        DZ['tot'], VZ['tot'], 
        times_DV=times_DV, 
        nz_init=nz_init, 
        plot=False
        )

    # Checks particle conservation
    cons = asim.check_conservation(plot=plt_all)
    if not plt_all:
        Ntot = cons["integ_source"][-1]
        dN = np.trapz(
            (cons["total"] / Ntot - cons["integ_source"] / Ntot) ** 2, asim.time_out
        )
        dN /= np.trapz((cons["integ_source"] / Ntot) ** 2, asim.time_out)
        print("Particle conservation error %.1f%%" % (np.sqrt(dN) * 100))

    # extract densities and particle numbers in each simulation reservoir
    (
        nz, N_wall, N_div, N_pump,
        N_ret, N_tsu, N_dsu, N_dsul, 
        rcld_rate, rclw_rate
        ) = out

    # extract flux surface volumes from geqdsk
    psin_ref = geqdsk["fluxSurfaces"]["geo"]["psin"]
    rhop_ref = np.sqrt(psin_ref)  # sqrt(norm. pol. flux)
    vol = interp1d(
        rhop_ref, 
        geqdsk["fluxSurfaces"]["geo"]["vol"],
        bounds_error=False, 
        fill_value=geqdsk["fluxSurfaces"]["geo"]["vol"][-1]
        )(asim.rhop_grid)

    # Rescale to target impurity concentration
    nz = _rscl_nz(
        dmodel=dmodel, 
        nz=nz,
        geqdsk = geqdsk,
        asim = asim,
        vol = vol,
        )

    # Calculates radiated power
    asim.rad = aurora.compute_rad(
        asim.imp, 
        nz.transpose(2,1,0), 
        asim.ne, 
        asim.Te, 
        prad_flag=True
        )
    Prad_tot = np.trapz(asim.rad['tot'][-1,:],vol) # [MW]
    print('P_rad_imp ={:.3f} MW'.format(Prad_tot))

    # Calculates Zeff contribution
    asim.calc_Zeff()
    Zeff_avg = np.trapz(np.sum(asim.delta_Zeff, axis = 1)[:,-1], vol)/vol[-1]
    print('dZ_eff ={:.3f}'.format(Zeff_avg))

    # Calculates the average concentration
    c_imp = np.sum(nz, axis=1)[:,-1]/asim.ne[0,:]
    c_imp_avg = np.trapz(c_imp, vol)/vol[-1]
    print('c_imp ={:.3e}'.format(c_imp_avg))

    # Plotting
    if plt_all:
        _plot(
            asim=asim,
            nz=nz,
            DZ=DZ,
            VZ=VZ,
            plt_cs=plt_cs,
            )

    # Formats the outputted density profile to conform with gacode
    nz_ga = interp1d(
        asim.rhop_grid,
        nz[:,:,-1],
        axis=0,
        )(kp['ne']['rhop']) # [1/cm3], dim(rhop,cs)

    # Output impurity density profile
    if 'FACIT' in dmodel['options']:
        return {
            'rhop_fm': asim.rhop_grid,
            'nz_fm': nz, # [1/cm3], dim(rhop, cs, t)
            'DZ_fm': DZ, # [cm^2/s], dim(rhop, cs)
            'VZ_fm': VZ, # [cm/s], dim(rhop, cs)
            'asym_fm': asym,
            'TSC_fm': TSC,
            'nz_ga': nz_ga, # [1/cm3], dim(rhop,cs,t)
            }
    else:
        return {
            'rhop_fm': asim.rhop_grid,
            'nz_fm': nz, # [1/cm3], dim(rhop, cs, t)
            'DZ_fm': DZ, # [cm^2/s], dim(rhop, cs)
            'VZ_fm': VZ, # [cm/s], dim(rhop, cs)
            'nz_ga': nz_ga, # [1/cm3], dim(rhop,cs,t)
            }


###########################################
#
#           Extra
#
###########################################

def _plot(
    asim=None,
    nz=None,
    DZ=None,
    VZ=None,
    plt_cs=None,
    ):

    # Plots impurity density
    aurora.slider_plot(
        asim.rvol_grid, 
        #asim.rhop_grid,
        asim.time_out, 
        nz.transpose(1,0,2),
        xlabel=r'$r_V$ [cm]', ylabel='time [s]', zlabel=r'$n_z$ [$cm^{-3}$]',
        labels=[str(i) for i in np.arange(0,nz.shape[1])],
        plot_sum=True, x_line=asim.rvol_lcfs
        )

    # Plots radiation
    aurora.slider_plot(
        asim.rhop_grid,
        asim.time_out, 
        asim.rad['line_rad'].transpose(1,2,0),
        xlabel=r'$\rho_p$', ylabel='time [s]', zlabel=r'Total radiation [$MW/m^3$]',
        labels=[str(i) for i in np.arange(0,nz.shape[1])],
        plot_sum=True, x_line = 1.0
        )

    # Plots delta Zeff
    aurora.slider_plot(
        asim.rhop_grid,
        asim.time_out, 
        asim.delta_Zeff.transpose(1,0,2),
        xlabel=r'$\rho_p$', ylabel='time [s]', zlabel=r'$\Delta$ $Z_{eff}$',
        labels=[str(i) for i in np.arange(0,nz.shape[1])],
        plot_sum=True, x_line = 1.0
        )

    # Plots diffusion, convection coefficients
    if plt_cs is None:
        plt_cs = [cs for cs in np.arange(0,nz.shape[1])]

    fig,ax = plt.subplots(3,1)

    for src in DZ.keys():
        if src == 'tot':
            fmt = '-'
        else:
            fmt = '--'

        for cs in plt_cs:
            ax[0].plot(
                asim.rhop_grid,
                DZ[src][:,0,cs]/1e4,
                fmt,
                label = src+' '+str(cs),
                )

            ax[1].plot(
                asim.rhop_grid,
                VZ[src][:,0,cs]/1e2,
                fmt,
                label = src+' '+str(cs),
                )

            ax[2].plot(
                asim.rhop_grid,
                VZ[src][:,0,cs]/DZ[src][:,0,cs]*1e2,
                fmt,
                label = src+' '+str(cs),
                )

    leg = ax[0].legend()
    leg.set_draggable('on')

    ax[0].set_xlabel(r'$\rho_p$')
    ax[0].set_ylabel(r'$D_Z$ [$m^2/s$]')
    ax[0].grid('on')
    #ax[0].set_ylim(0,3)

    ax[1].set_xlabel(r'$\rho_p$')
    ax[1].set_ylabel(r'$V_Z$ [$m/s$]')
    ax[1].grid('on')
    #ax[1].set_ylim(-3,3)

    ax[2].set_xlabel(r'$\rho_p$')
    ax[2].set_ylabel(r'$V_Z/D_Z$ [$1/m$]')
    ax[2].grid('on')
    ax[2].set_ylim(-10,5)


def _rscl_nz(
    dmodel=None,
    nz=None,
    geqdsk = None,
    asim = None,
    vol = None,
    ):

    # Rescale to desired density
    if dmodel['AURORA']['target']['c_imp'] is not np.nan:
        # Calculates the average concentration
        c_imp_tmp = np.sum(nz, axis=1)[:,-1]/asim.ne[0,:]
        c_imp_avg_tmp = np.trapz(c_imp_tmp, vol)/vol[-1]

        # Rescales impurity density
        nz *= dmodel['AURORA']['target']['c_imp']/c_imp_avg_tmp

    elif dmodel['AURORA']['target']['P_rad'] is not np.nan:
        # Calculates radiated power
        rad_tmp = aurora.compute_rad(
            asim.imp, 
            nz.transpose(2,1,0), 
            asim.ne, 
            asim.Te, 
            prad_flag=True
            )
        Prad_tot_tmp = np.trapz(rad_tmp['tot'][-1,:],vol) # [MW]

        # Rescales impurity density
        nz *= dmodel['AURORA']['target']['P_rad']/Prad_tot_tmp

    # Rescale density profiles to rho=1 value
    elif dmodel['AURORA']['target']['BC'] is not np.nan:
        ind = np.argmin(
            abs(
                asim.rhop_grid-dmodel['AURORA']['target']['BC']['rhop']
                )
            )

        nz *= dmodel['AURORA']['target']['BC']['val']/np.sum(nz[ind,:,-1])


    return nz

def _get_DV(
    dmodel = None,
    asim = None,
    geqdsk = None,
    inputgacode = None,
    times_DV = None,
    nz_init = None,
    ):

    # Initialize D,V profiles
    DZ = {}
    VZ = {}
    DZ['tot'] = np.zeros(
        (asim.rvol_grid.size, times_DV.size, asim.Z_imp+1)
        )  # [cm2/s], dim(fm_rhop, t, Z)
    VZ['tot'] = np.zeros(
        (asim.rvol_grid.size, times_DV.size, asim.Z_imp+1)
        )  # [cm/s], dim(fm_rhop, t, Z)

    # Loop over modeling options
    for opt in dmodel['options']:
        # If user wants to just use flat profiles
        if opt == 'Flat':
            DZ_tmp = dmodel['Flat']['D'] * np.ones(DZ['tot'].shape)
            VZ_tmp = dmodel['Flat']['V'] * np.ones(DZ['tot'].shape)

        # If user wants to use TGLF predictions
        elif opt == 'TGLF':
            # Imports run_TGLF module
            from transport_world.run_TGYRO import run_TGLF as rT

            # Obtains data
            dout = rT.calc_imp_turbDV(
                fgacode = os.path.join(
                    dmodel['paths']['in_path'],
                    dmodel['paths']['fgacode'],
                    ),
                folder = os.path.join(
                    dmodel['paths']['in_path'],
                    dmodel['paths']['out_folder'],
                    dmodel['paths']['name_sim']
                    ),
                subfolder = dmodel['paths']['name_sim'],
                cs = dmodel['imp']['cs'],
                amu = dmodel['imp']['amu'],
                restart = dmodel[opt]['restart'],
                TGLFsettings = dmodel[opt]['TGLFsettings'],
                )

            # Smooths TGLF output
            dout = rT.smooth_DV(
                ddata=dout,
                rho_accept = dmodel[opt]['rho_accept']
            )

            rhop_TGLF = interp1d(
                inputgacode['rho'],
                np.sqrt(inputgacode['polflux']/inputgacode['polflux'][-1])
                )(dout['rhot'])

            # Interpolates onto simulation grid
            DZ_tmp = interp1d(
                rhop_TGLF,
                dout['DZ'],
                bounds_error = False,
                fill_value = 0.1,
                )(asim.rhop_grid)*1e4 # [cm2/s], dim(fm_rhop,)
            VZ_tmp = interp1d(
                rhop_TGLF,
                dout['VZ'],
                bounds_error = False,
                fill_value = -0.1,
                )(asim.rhop_grid)*1e2 # [cm/s], dim(fm_rhop,)

            # Reshape
            DZ_tmp = np.repeat(
                np.repeat(
                    DZ_tmp[:,None],
                    times_DV.size,
                    axis=-1
                    )[:,:,None],
                asim.Z_imp+1,
                axis=-1
                ) # dim(fm_rhop, t, Z)
            VZ_tmp = np.repeat(
                np.repeat(
                    VZ_tmp[:,None],
                    times_DV.size,
                    axis=-1
                    )[:,:,None],
                asim.Z_imp+1,
                axis=-1
                ) # dim(fm_rhop, t, Z)

        # If user wants to use FACIT predictions
        elif opt == 'FACIT':
            # Imports run_FACIT module
            from transport_world.run_AURORA import run_FACIT as rF

            # Obtains data
            dout = rF.calc_imp_neoDV(
                asim = asim,
                geqdsk = geqdsk,
                inputgacode = inputgacode,
                ICRH_option = dmodel['FACIT']['ICRH_option'],
                rotation_model = dmodel['FACIT']['rotation_model'],
                times_DV = times_DV,
                nz_init = nz_init,
                )

            # Adds neoclassical transport
            DZ_tmp = dout['DZ']
            VZ_tmp = dout['VZ']
            asym = dout['asym']
            TSC = dout['TSC']

        # Stores
        DZ['tot'] += DZ_tmp
        VZ['tot'] += VZ_tmp
        DZ[opt] = DZ_tmp
        VZ[opt] = VZ_tmp

    if 'FACIT' in dmodel['options']:
        return DZ, VZ, asym, TSC
    else:
        return DZ, VZ