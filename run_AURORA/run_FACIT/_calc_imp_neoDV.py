'''

_calc_imp_neoDV is a module that calculate neoclassical impurity
D, V profiles using FACIT within Aurora

'''

# Module
import aurora
import scipy.constants as cnt
from scipy.interpolate import griddata

__all__ = [
    'calc_imp_neoDV',
    ]

########################################
#
#           Main
#
########################################

def calc_imp_neoDV(
    asim = None,        # Aurora simulation object
    geqdsk = None,      # Magnetic equilibrium
    inputgacode = None, # kinetic profiles
    ICRH_option = None, # Minority heating physics
    rotation_model = None,
    times_DV = None,    # Time-domain for D,V profiles
    nz_init = None,     # Initial impurity density
    ):

    # index to the spearatrix
    idxsep = np.argmin(
        np.abs(1.0 - asim.rhop_grid)
    ) 

    # volume-average radial coord
    rr = asim.rvol_grid / 100  # dim(fm_rhop); [m]

    # Minor radius
    amin = rr[idxsep]  # [m]

    # normalized radial coordinate
    roa = rr[: idxsep + 1] / amin # dim(fm_rhop); 0 -> 1

    # magnetic field on axis
    B0 = np.abs(geqdsk["BCENTR"]) # [T]

    # major radius
    R0 = geqdsk["fluxSurfaces"]["R0"] # [m]

    # safety factor profile
    qmag = np.interp(
        roa, 
        geqdsk["RHOVN"], 
        geqdsk["QPSI"]
        )[: idxsep + 1]  # dim(fm_rhop)

    # sq. norm. pol. flux
    rhop = asim.rhop_grid[: idxsep + 1] # dim(fm_rhop); 0_.1
    rhop_ga = np.sqrt(inputgacode['polflux']/inputgacode['polflux'][-1])

    # main ion density profile
    ni_tmp = 0
    for key in inputgacode['IONS']:
        if inputgacode['IONS'][key][0] in ['T', 'D']:
            ni_tmp += inputgacode['ni_'+str(key)]*1e13 # dim(rhop); [1/cm3] 
    Ni = (
        np.interp(
            rhop, 
            rhop_ga, 
            ni_tmp
            ) * 1e6
    )  # in m**3 instead of cm**3 in FACIT
    gradNi = np.gradient(Ni, roa*amin) # dim(fm_rhop); [1/m4]

    # Ion temperature
    Ti = np.interp(
        rhop, 
        rhop_ga, 
        inputgacode['Ti_1']
        ) *1e3 # dim(fm_rhop); [eV]
    gradTi = np.gradient(Ti, roa*amin) # dim(fm_rhop): [eV/m]

    # Ratio of electron to ion temperature
    TeovTi = np.interp(
        rhop, 
        rhop_ga, 
        inputgacode['Te']/inputgacode['Ti_1']
        ) # dim(fm_rhop)  

    # Z effective
    Zeff = np.interp(
        rhop,
        rhop_ga,
        inputgacode['Zeff']
        ) # dim(fm_rhop,) 

    # Ion Mach number profile
    Machi = np.interp(
        rhop,
        rhop_ga,
        (
            inputgacode['vtor']
            /np.sqrt(
                2*cnt.e*inputgacode['Ti_1']*1e3
                /asim.main_ion_A/cnt.m_p
                )
            )
        ) # dim(fm_rhop,)

    # Rotation modeling controls
    # (0-> non rotating, 1-> PS model, 2 -> BP and PS model)
    # Recommend 0 for light impurities, 2 for anything bigger than Ar

    # use non-rotating limit from Fajardo PPCF (2022)
    if rotation_model == 0:

        # Dummy values for unused parameters
        BV = None
        FV = None
        dpsidx = None
        full_geom = False
        RV = None
        ZV = None
        nth = 20

    # use PS model from Maget PPCf (2020)
    elif rotation_model == 1:

        # Poloidal angle grid
        nth = 51
        theta = np.linspace(0, 2 * np.pi, nth)

        # R,Z flux-surface contours
        RV, ZV = aurora.rhoTheta2RZ(geqdsk, rhop, theta, coord_in="rhop", n_line=201)
        RV, ZV = RV.T, ZV.T # dim(fm_rhop, fm_theta); [m]

        # 2D toroidal magnetic field; [T]
        BV_g = geqdsk['AuxQuantities']['Bt'] # dim(g_R, g_Z)
        R_g = geqdsk['AuxQuantities']['R'] # dim(g_R)
        Z_g = geqdsk['AuxQuantities']['Z'] # dim(g_Z)
        RR_g, ZZ_g = np.meshgrid(R_g, Z_g) # dim(g_R, g_Z)
        BV = griddata((RR_g.flatten(), ZZ_g.flatten()), 
            BV_g.flatten(),
            (RV.flatten(), ZV.flatten()),
            method='linear').reshape(RV.shape) # dim(fm_rhop, fm_theta)

        # Poloidal current flux; [T*m]
        FV = np.interp(
            roa, 
            geqdsk["RHOVN"], 
            geqdsk["FPOL"]
            )[: idxsep + 1] # dim(fm_rhop)
        
        # Radial coord transform from psi to rho; [T*m^2]
        psi = np.interp(
            roa, 
            geqdsk["RHOVN"], 
            geqdsk['fluxsurfaces']['geo']['psi']
            )[: idxsep + 1]
        dpsidx = np.gradient(psi, roa) # dim(fm_rhop)

        # Utilize circular (False) or full (True) geometry
        full_goem = True

    # use BP and PS model from Fajardo (PPCF 2023)
    elif rotation_model == 2:

        # Poloidal angle grid
        nth = 51
        theta = np.linspace(0, 2 * np.pi, nth)

        # R,Z flux-surface contours
        RV, ZV = aurora.rhoTheta2RZ(geqdsk, rhop, theta, coord_in="rhop", n_line=201)
        RV, ZV = RV.T, ZV.T # dim(fm_rhop, fm_theta); [m]

        # Dummy values for unused parameters
        BV = None
        FV = None
        dpsidx = None
        full_geom = False

    # -------------------
    # ICRH controls
    # -------------------

    if ICRH_option is None:
        fH = 0.0
        bC = 1.0
        sigH = 1.0
        TperpTpar_axis = 1.0
        
    else:
        print('ICRH not yet implemented')

    # -------------------
    # Runs FACIT
    # -------------------

    # Initializes matrices to store transport coefficients
    D_z = np.zeros((asim.rvol_grid.size, times_DV.size, asim.Z_imp+1)) # dim(rm_rhop, t, Z); [cm^2/s]
    V_z = np.zeros(D_z.shape) # dim(fm_rhop, t, Z); # [cm/s]
    horz_asym = np.zeros(D_z.shape) # dim(fm_rhop, t, Z)
    vert_asym = np.zeros(D_z.shape) # dim(fm_rhop, t, Z)
    TSC = np.zeros(D_z.shape) # dim(fm_rhop, t, Z)
    nn = np.zeros((asim.rvol_grid.size, nth, times_DV.size, asim.Z_imp+1)) # dim(fm_rhop, fm_theta, t, Z)

    # Loop over times
    for jj, tt in enumerate(times_DV):

        # Loop over charge states
        for ii, zz in enumerate(range(asim.Z_imp + 1)):

            if zz != 0:
                # Charge state density profile
                Nz     = nz_init[:idxsep+1,ii]*1e6 # dim(fm_rhop); [1/m3]
                gradNz = np.gradient(Nz, roa*amin) # dim(fm_rhop); [1/m4]

                # Runs FACIT
                fct = aurora.FACIT(
                    roa,                # dim(fm_rhop); rho = r/a
                    zz,                 # scalar; impurity charge 
                    asim.A_imp,         # scalar; impurity mass
                    asim.main_ion_Z,    # scalar; main ion charge
                    asim.main_ion_A,    # scalar; main ion mass
                    Ti,                 # dim(fm_rhop); [eV]; ion temperature
                    Ni,                 # dim(fm_rhop); [1/m3]; main ion density
                    Nz,                 # dim(fm_rhop); [1/m3]; Initial charge state density profile
                    Machi,              # scalar or dim(fm_rhop); Mach number
                    Zeff,               # scalar or dim(fm_rhop); Z effective
                    gradTi,             # dim(fm_rhop); [eV/m]; ion temperature gradient
                    gradNi,             # dim(fm_rhop); [1/m4]; main ion density gradient
                    gradNz,             # dim(fm_rhop); [1/m4]; impurity charge state density gradient
                    amin/R0,            # scalar; inverse aspect raio 
                    B0,                 # scalar; [T]; central magnetic field
                    R0,                 # scalar; [m]; major radius
                    qmag,               # dim(fm_rhop); q-profile
                    rotation_model = rotation_model, # scalar; rotation model switch
                    Te_Ti = TeovTi,     # optional; dim(fm_rhop); Te/Ti
                    RV = RV,            # dim(fm_rhop, fm_theta); ['m]; R flux-surface contours
                    ZV = ZV,            # dim(fm_rhop, fm_theta); ['m']; Z flux-surface contours
                    BV = BV,            # dim(fm_rhop, fm_theta); ['T']; B-field flux-surface contours
                    FV = FV,            # dim(fm_rhop, fm_theta); ['T*m']; poloidal current flux-surface contours
                    dpsidx = dpsidx,    # dim(fm_rhop); ['T*m2']; radial coord transform from psi to roa
                    fsaout = True,      # bool; if True, return FSA transport coefficients
                    full_geom = full_geom, # bool; if True, uses full flux-surface geom, if False, uses circular limit
                    nat_asym = True,    # bool; if True, includes friction induced asymmetry
                    nth = nth,           # scalar; number of poloidal grid points
                    regulopt = [1e-2,0.5,1e-5,1e2], # list; convergence params for iterative soln
                    # ICRH controls
                    fH = fH,            # scalar; hydrogen minority fraction
                    bC = bC,            # scalar; Bres/B0, where Bres is the resonant magnetic field of the ICRH
                    sigH = sigH,          # scalar; std. dev. of radial exp decay of temp anisotropy from ICRH
                    TperpTpar_axis = TperpTpar_axis,# scalar; main ion temp anisotropy at mangetic axis
                    )

                # Stores diffusion coefficient
                D_z[:idxsep+1,jj,ii] = fct.Dz*100**2 # dim(fm_rhop, t, Z); convert to cm**2/s
                V_z[:idxsep+1,jj,ii] = fct.Vconv*100 # dim(fm_rhop, t, Z); convert to cm/s

                # Stores poloidal asymmetry
                # Modeled as nz/<nz> = 1 + horz_asym * cos(theta) + vert_asym * sin(theta)
                horz_asym[:idxsep+1,jj,ii] = fct.horiz_asym # dim(fm_rhop, t, Z)
                vert_asym[:idxsep+1,jj,ii] = fct.vert_asym # dim(fm_rhop, t, Z)

                # Stores temperature screening coefficient
                TSC[:idxsep+1,jj,ii] = fct.Hz/fct.Kz # dim(fm_rhop, t, Z)

                # Stores poloidal distribution, nz(r,theta)/<nz>(r)
                nn[:idxsep+1,:,jj,ii] = fct.nn # dim(fm_rhop, fm_theta, t, Z)

    return {
        'DZ': D_z,
        'VZ': V_z,
        'TSC': TSC,
        'asym': {
            'hor_asym': hor_asym,
            'vert_asym': vert_asym,
            'nn': nn,
            'R': RV,
            'Z': ZV,
            'theta': theta,
            'rhop': asim.rhop_grid,
            }
        }

