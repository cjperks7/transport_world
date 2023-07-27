'''

_calc_imp_dens is a function to calculate impurity density 
profiles using Aurora

'''

# Modules
import aurora
from scipy.interpolate import interp1d
from omfit_classes import omfit_eqdsk, omfit_gapy

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
    ):

    # Reads in magnetic equilibrium and kinetic profiles
    geqdsk = omfit_eqdsk.OMFITgeqdsk(
        os.path.join(
            dmodel['paths']['folder'],
            dmodel['paths']['geqdsk']
            )
        )
    inputgacode = omfit_gapy.OMFITgacode(
        os.path.join(
            dmodel['paths']['folder'],
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
    # Edge modeling
    # -------------------

    # Source location wrt LCFS
    nml['source_cm_out_lcfs'] = dmodel['AURORA']['source_cm_out_lcfs']

    # Recyling
    nml['recycling_flag'] = dmodel['AURORA']['recycling_flag'] # recycling from divertor
    nml['wall_recycling'] = dmodel['AURORA']['wall_recycling'] # recycling from limiter

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
        DZ, VZ, 
        times_DV=times_DV, 
        nz_init=nz_init, 
        plot=plt_all
        )

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


###########################################
#
#           Extra
#
###########################################

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
    DZ = np.zeros(
        (asim.rvol_grid.size, times_DV.size, asim.Z_imp+1)
        )  # [cm2/s], dim(fm_rhop, t, Z)
    VZ = np.zeros(
        (asim.rvol_grid.size, times_DV.size, asim.Z_imp+1)
        )  # [cm/s], dim(fm_rhop, t, Z)

    # Loop over modeling options
    for opt in dmodel['options']:
        # If user wants to just use flat profiles
        if opt == 'Flat':
            DZ += dmodel['Flat']['D'] * np.ones(DZ.size)
            VZ += dmodel['Flat']['V'] * np.ones(DZ.size)

        # If user wants to use TGLF predictions
        elif opt == 'TGLF':
            # Imports run_TGLF module
            from transport_world.run_TGYRO import run_TGLF as rT

            # Obtains data
            dout = rT.calc_imp_turbDV(
                fgacode = dmodel['paths']['fgacode'],
                folder = dmodel['paths']['folder'],
                subfolder = dmodel['paths']['subfolders'],
                cs = dmodel['imp']['cs'],
                amu = dmodel['imp']['amu'],
                restart = dmodel[opt]['restart'],
                TGLFsettings = dmodel[opt]['TGLFsettings'],
                )

            rhop_TGLF = interp1d(
                inputgacode['rho'],
                np.sqrt(inputgacode['polflux']/inputgacode['polflux'][-1])
                )(dout['rhot'])

            # Interpolates onto simulation grid
            DZ += interp1d(
                rhop_TGLF,
                dout['DZ'],
                bounds_error = false,
                fill_value = (dout['DZ'][0], dout['DZ'][-1])
                )(asim.rhop_grid)[:,None,None] # [cm2/s], dim(fm_rhop, t, Z)
            VZ += interp1d(
                rhop_TGLF,
                dout['VZ'],
                bounds_error = false,
                fill_value = (dout['VZ'][0], dout['VZ'][-1])
                )(asim.rhop_grid)[:,None,None] # [cm/s], dim(fm_rhop, t, Z)

        # If user wants to use FACIT predictions
        elif opt == 'FACIT':
            # Imports run_FACIT module
            from transport_world.run_AURORA import run_FACIT as rF

            # Obtains data
            dout = rF.calc_imp_neoDV(
                asim = asim,
                geqdsk = geqdsk,
                inputgacode = inputgacode,
                ICRH_option = None,
                rotation_model = dmodel['FACIT']['rotation_model']
                times_DV = times_DV,
                nz_init = nz_init,
                )

            # Adds neoclassical transport
            DZ += dout['DZ']
            VZ += dout['VZ']
            asym = dout['asym']
            TSC = dout['TSC']

    if 'FACIT' in dmodel['options']:
        return DZ, VZ, asym, TSC
    else:
        return DZ, VZ