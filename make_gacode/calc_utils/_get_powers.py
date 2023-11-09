'''

get_powers.py is a module to facilitate calcualting
    1) RF
    2) Ohmic
    3) Fusion
    4) Radiation power


'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cnt
from omfit_classes import omfit_eqdsk
import aurora
from . import _get_sigv as _sv
import math
from scipy.integrate import cumtrapz

plt.rcParams.update({'font.size': 16})

__all__ = [
    'get_powers',
    ]


##########################################################
#
#                     Main
#
##########################################################

def get_powers(
    dout=None,
    verb=None,
    Prf_in = 1, # [Mw]
    plt_all = None,
    ):

    # Initialize dictionary
    dout['powers'] = {}

    # Calculates powers
    dout = _get_qsync(
        dout=dout,
        verb=verb,
        )
    dout = _get_qrad(
        dout=dout,
        verb=verb,
        )
    dout = _get_qohm(
        dout=dout,
        verb=verb,
        )
    dout = _get_qrf(
        dout=dout,
        verb=verb,
        Prf_in=Prf_in,
        )
    dout= _get_qfus(
        dout=dout,
        verb=verb,
        )

    # Not implemented powers
    dout = _get_qNA(
        dout=dout,
        )

    # Plotting
    if plt_all:
        _plot(dout=dout)

    return dout

##########################################################
#
#                     PLotting
#
##########################################################

def _plot(dout=None):

    skip = [
        'beame',
        'beami',
        'ione',
        'ioni',
        'ei',
        'cxi',
        'par_beam',
        'par_wall',
        'mom',
        ]

    fig, ax =plt.subplots()

    # Loop over components
    for pwr in dout['powers'].keys():
        if pwr in skip:
            continue
        
        ax.plot(
            dout['rhot'],
            dout['powers'][pwr]['prof_MW/m3'],
            label = pwr
        )

    ax.set_xlabel(r'$\rho_t$')
    ax.set_ylabel(r'$q$ [MW/m3]')
    ax.legend()
    ax.grid('on')


##########################################################
#
#                     Extra
#
##########################################################

def _get_qNA(
    dout=None,
    ):

    # Not implemented sources
    pwr_NA = [      # heat sources
        'ei',       # electron-ion exchange
        'beame',    # beam to e's
        'beami',    # beam to ions's
        'ione',     # recombination to e's
        'ioni',     # recombination to ion's
        'cxi',      # charge-exchange to ion's
        ]
    par_NA = [      # particle sources
        'par_beam', # beam
        'par_wall', # wall
        ]
    mom_NA = [      # momentum sources
        'mom',
        ]

    # Stores
    for pwr in pwr_NA:
        dout['powers'][pwr] = {}
        dout['powers'][pwr]['tot_MW'] = 0
        dout['powers'][pwr]['prof_MW/m3'] = 0*dout['rhot']

    for par in par_NA:
        dout['powers'][par] = {}
        dout['powers'][par]['tot_1/s'] = 0
        dout['powers'][par]['prof_1/m3/s'] = 0*dout['rhot']

    for mom in mom_NA:
        dout['powers'][mom] = {}
        dout['powers'][mom]['tot_Nm'] = 0
        dout['powers'][mom]['prof_N/m2'] = 0*dout['rhot']

    return dout

##############
# Function to obtain synchrotron power density, Trubnikov, JETP Lett. 16 (1972) 25
# ... INPUT:  Te_keV  - vector, [keV], Electron temperature
#             ne_m3   - vector, [1e19 m^-3], Electron density
#             R       - scalar, [m], major radius
#             a       - scalar, [m], minor radius
#             B       - vector, [T], Magnetic field
#
# ... OUTPUT: qsync   - vector, [MW/m^3], synchrotron power density
def _get_qsync(
    dout=None,
    verb=None,
    ):
    # Reflection coefficient
    refl = 0.8

    # Physical constants
    pi = 4.0*math.atan(1.0)
    e = 4.8032e-10 # [statcoul]
    k = 1.6022e-12 # [erg/eV]
    me = 9.1094e-28 # [g]
    c = 2.9979e10 # [cm/s]

    # Converts units of inputs
    Te = dout['Te_keV']*1000 # [eV]
    ne = dout['ne_19m3']*1e13 # [cm^-3]
    R = float(dout['rcentr_m'])*100 # [cm]
    a = dout['rmin_m'][-1]*100 # [cm]
    B = abs(float(dout['bcentr_T']))*1e4 # [G]

    # Calculates aspect ratio
    A = R/a

    # Calculates the electron plasma frequency
    wpe = np.sqrt(4*pi*ne*e**2/me)

    # Calculates the electron cyclotron frequency
    wce = e*np.abs(B)/(me*c)

    # Calculates g-factor
    g = k*Te/(me*c**2)

    # Calculates phi-factor
    phi = 60*g**1.5*np.sqrt((1.0-refl)*(1+1/A/np.sqrt(g))/(a*wpe**2/c/wce))

    # Calculates synchrotron power density, [MW/m^2]
    qsync = me/(3*pi*c)*g*(wpe*wce)**2*phi/1e7

    # Printing values
    dout['powers']['sync'] = {}
    if verb == 2:
        Psync_tot = cumtrapz(qsync, dout['vol_m3'], initial=0.) # [MW], synchrotron radiation power
        print('Psync= {:1.2f} MW\n'.format(Psync_tot[-1]))
        dout['powers']['sync']['tot_MW'] = Psync_tot[-1]

    # Returns synchrotron power density, [MW/m^3], dim(rhot,)
    dout['powers']['sync']['prof_MW/m3'] = qsync
    return dout


##############
# Function to obtain continuum and line power density
# ... INPUT:  Te_keV  - vector, [keV], Electron temperature
#             ne_m3   - vector, [1e19 m^-3], Electron density
#             R       - scalar, [m], major radius
#             a       - scalar, [m], minor radius
#             B       - vector, [T], Magnetic field
#
# ... OUTPUT: qsync   - vector, [MW/m^3], synchrotron power density
def _get_qrad(
    dout=None,
    verb=None,
    ):
    # Output data from gfile
    gfile = omfit_eqdsk.OMFITgeqdsk(dout['paths']['input']+dout['paths']['gfile'])

    # Initialize
    qline = np.zeros(len(dout['rhot']))
    qcont = np.zeros(len(dout['rhot']))

    for ion in list(filter(None,dout['ions']['name'].split(' '))):
        if ion == 'D' or ion == 'H' or ion == 'T':
            imp = 'H'
        elif ion == 'LUMPED':
            imp = 'B'
        else:
            imp = ion
        # get charge state distributions from ionization equilibrium for impurity ion, 
        # (scd-> ionization rates, acd-> recombination rates, ccd->charge exchange rates) 
        atom_data = aurora.atomic.get_atom_data(imp,['scd','acd']) 

        if 'nz_19m3' in list(dout['ions'][ion].keys()):
            # Calculates the line & continuum radiated power for impurity, He, T, D
            rad = aurora.radiation_model(
                imp,
                dout['rhop'],
                dout['ne_19m3']*1e13,
                dout['Te_keV']*1000, 
                gfile,
                n0_cm3=None, 
                nz_cm3=dout['ions'][ion]['nz_19m3']*1e13, 
                plot=False
                ) # rad in [W/m^3]

        else:
            # Calls AURORA function to calculate fractional abundances assuming ionization equilibrium
            _, fz = aurora.atomic.get_frac_abundances(
                atom_data, 
                dout['ne_19m3']*1e13, 
                dout['Te_keV']*1e3,
                rho=dout['rhop'], 
                plot=None, 
                ax=None
                )# plot=plot,ax=plt.gca())

            # Calculates the line & continuum radiated power for impurity, He, T, D
            rad = aurora.radiation_model(
                imp,
                dout['rhop'],
                dout['ne_19m3']*1e13,
                dout['Te_keV']*1000, 
                gfile,
                n0_cm3=None, 
                nz_cm3=dout['ions'][ion]['ni_tot_19m3'][:,None]*fz*1e13, 
                plot=False
                ) # rad in [W/m^3]

        # Total radiated power density, [MW/m^3]
        qline += rad['line_rad_dens'].sum(0)/1e6
        qcont += rad['cont_rad_dens'].sum(0)/1e6

    # Printing values
    dout['powers']['cont'] = {}
    dout['powers']['line'] = {}
    dout['powers']['rad'] = {}
    if verb == 2:
        Pcont_tot = cumtrapz(qcont, dout['vol_m3'], initial=0.) # [MW], continuum radiation power
        Pline_tot = cumtrapz(qline, dout['vol_m3'], initial=0.) # [MW], line radiation power
        print('Pcont= {:1.2f} MW\n'.format(Pcont_tot[-1]))
        print('Pline= {:1.2f} MW\n'.format(Pline_tot[-1]))
        dout['powers']['cont']['tot_MW'] = Pcont_tot[-1]
        dout['powers']['line']['tot_MW'] = Pline_tot[-1]
        dout['powers']['rad']['tot_MW'] = Pline_tot[-1] + Pcont_tot[-1]

    # Returns radiation power density, [MW/m^3], dim(rhot,)
    dout['powers']['cont']['prof_MW/m3'] = qcont
    dout['powers']['line']['prof_MW/m3'] = qline
    dout['powers']['rad']['prof_MW/m3'] = qline + qcont

    return dout


##############
# Function to obtain ohmic power density
# ... INPUT:  Te      - vector, [keV], Electron temperature
#             ne      - vector, [1e19 m^-3], Electron density
#             R       - scalar, [m], major radius
#             a       - scalar, [m], minor radius
#             Zeff    - vector, [], Effective nuclear charge
#             johm    - vector, [MA/m^2], ohmic current density
#             jbstor  - vector, [MA/m^2], toroidal bootstrap current density
#
# ... OUTPUT: qomh    - vector, [MW/m^3], ohmic power density
def _get_qohm(
    dout=None,
    verb=None,
    ):
    # Constants
    e = 1.6021e-19
    me = 9.11e-31
    eps0 = 8.85e-12

    # Values
    Te = dout['Te_keV'] # [keV]
    ne = dout['ne_19m3'] # [1e19 m^-3]
    R = float(dout['rcentr_m']) # [m]
    a = dout['rmin_m'][-1] # [m]
    Zeff = dout['Zeff']['prof']
    johm = dout['johm_MA/m2'] # [MA/m^2]
    jbstor = dout['jbstor_MA/m2'] # [MA/m^2]

    # Calculates neoclassical resistivity correction factor
    gamma_R = (1 - 1.95*np.sqrt(a/R) + 0.95*a/R)**-1

    # Calculates charge correction factor
    F_Z = (1+1.198*Zeff+0.222*Zeff**2) /(1+2.966*Zeff+0.753*Zeff**2)

    # Calculates Coulomb logarithm
    lnLamb = 24 - np.log(np.sqrt(ne*1e19/1e6)*(Te*1000)**-1)

    # Calculates the electron-electron collision frequency [1/s]
    nu_ee = 1/(3*np.sqrt(np.pi)) *ne*1e19 *(e**2/(4*np.pi*eps0))**2 *4*np.pi/ \
        (np.sqrt(me)*(Te*1000*e)**(3/2)) *lnLamb

    # Calculates the neoclassical parallel Spitzer resistivity [Ohm *m]
    eta_para = (me*Zeff*nu_ee*gamma_R*F_Z)/(ne*1e19*e**2)

    # Calculates the ohmic power density [MW/m^3]
    qohm = eta_para * (johm*1e6 + jbstor*1e6)**2 / 1e6
    qohm[-1] = qohm[-2]

    # Printing values
    dout['powers']['ohm'] = {}
    if verb == 2:
        Pohm_tot = cumtrapz(qohm, dout['vol_m3'], initial=0.) # [MW], ohmic power
        print('Pohm= {:1.2f} MW\n'.format(Pohm_tot[-1]))
        dout['powers']['ohm']['tot_MW'] = Pohm_tot[-1]

    # Returns ohmic power density, [MW/m^3], dim(rhot,)
    dout['powers']['ohm']['prof_MW/m3'] = qohm

    return dout


##############
# Function to obtain LH heating power density
# ... INPUT:  a       - scalar, [m], minor radius
#             rho     - vector, [], sqrt. norm. tor. flux coordinate grid
#             vol     - vector, [m^3], Volume of the plasma
#             Prf_in  - scalar, [MW], RF power injected into the plasma
#
# ... OUTPUT: qLH     - vector, [MW/m^3], LH heating power density
def _get_qrf(
    dout=None,
    Prf_in = None, 
    verb = None,
    rf_file = None,
    ):

    # If it is wished that an Gaussian power deposition is used
    if rf_file is None:
        # Sets the spread (\sigma) of the guassian
        # Spread in real space (full width)
        dr = 0.1 # [m]

        # Calculates \sigma in r/a-space
        sigma = dr/2/(dout['rmin_m'][-1]) #[norm]

        # Sets the center (\mu) of the gaussian in r/a-space
        mu = 0.7 # [norm]

        # Initializes shape function for RF power deposition
        qrfe = np.exp(-((dout['rhot']-mu)/(sigma))**2)

        #Normalizes the shape function
        qrfe = Prf_in*qrfe/(np.trapz(qrfe, dout['vol_m3'])) # [MW/m^3]
        qrfe[qrfe<1e-10] = 0

    else:
        # Reads netCDF file storing power deposition profile
        from scipy.interpolate import interp1d
        from scipy.io import netcdf
        data = netcdf.NetCDFFile(file, 'r')

        ##### NOTE: These sims were done on a grid rho = sqrt(norm tor) ######

        # Loads in time-dependent power deposition profile
        rfpwr = np.copy(data.variables['powrft'].data)

        # Takes the power at the last time slice of RF sim
        pow1 = rfpwr[-1,:]

        # Loads in the rho = sqrt(norm tor) grid used
        rfrho = np.copy(data.variables['rya'][:].data)

        # Appends the necessary boundaries values
        rfrho = np.insert(rfrho, 0, 0)
        rfrho = np.append(rfrho, [1])
        pow1 = np.insert(pow1, 0, pow1[0])
        pow1 = np.append(pow1, pow1[-1])

        # Interpolates the power deposition power into GACODE grid
        qrfe = interp1d(rfrho, pow1)(rho)
        qrfe[qrfe<1e-10] = 0

    # Printing values
    dout['powers']['rf'] = {}
    dout['powers']['rfe'] = {}
    dout['powers']['rfi'] = {}
    if verb == 2:
        Prfe_tot = cumtrapz(qrfe, dout['vol_m3'], initial=0.) # [MW], rf total power
        print('Prfe= {:1.2f} MW\n'.format(Prfe_tot[-1]))
        dout['powers']['rf']['tot_MW'] = Prfe_tot[-1]
        dout['powers']['rfe']['tot_MW'] = Prfe_tot[-1]
        dout['powers']['rfi']['tot_MW'] = 0

    # Returns rf power density, [MW/m^3], dim(rhot,)
    dout['powers']['rf']['prof_MW/m3'] = qrfe
    dout['powers']['rfe']['prof_MW/m3'] = qrfe
    dout['powers']['rfi']['prof_MW/m3'] = 0*qrfe
    return dout

##############
# Function to obtain fusion power density
# ... INPUT:  n_fuel  -- [array], [1e19 1/m^3] fuel ion density (assumes 50/50 DT)
#             ne      -- [array], [1e19 1/m^3] electron density
#             n_ash   -- [array], [1e19 1/m^3] ash ion density
#             n_imp   -- [array], [1e19 1/m^3] impurity ion density
#             z_imp   -- [scalar], [] impurity charge
#             mimp_amu-- [scalar], [amu] impurity mass 
#             Te      -- [array], [keV] electron temperature
#             Ti      -- [array], [keV] ion temperature
#             nexp    -- [scalar], [] number of grid points
#
# ... OUTPUT: qfus    -- [array], [MW/m^3] total fusion power density
#             qfusi   -- [array], [MW/m^3] alpha heating power density to ions
#             qfuse   -- [array], [MW/m^3] alpha heating power density to electrons
def _get_qfus(
    dout=None,
    verb=None,
    ):
    # Elementary charge, [C]
    e = 1.6021e-19

    # Energy from DT fusion, [MJ]
    Efus_DT = 17.59*e

    # fusion type
    if 'T' in dout['ions'].keys() and 'D' in dout['ions'].keys():
        reacts = ["DT"]
        n_fuel2 = (
            dout['ions']['D']['ni_tot_19m3'] * dout['ions']['T']['ni_tot_19m3']
            ) * 1e19**2 # [1/m^6]
    else:
        reacts = ["DD_n3He","DD_pT"]
        n_fuel2 = (
            dout['ions']['D']['ni_tot_19m3']/2
            )**2 * 1e19**2 # [1/m^6]

    # Initializes array to store data
    SigV = np.zeros(np.size(dout['Ti_keV']))
    for react in reacts:
        # Calculates the fusion reactivity, [m^3/s]
        SigV += _sv._get_sigmav(react, dout['Ti_keV'])

    # Calculates the alpha-ion heating fraction
    frac_ai = _sv._get_ai_frac(
        dout=dout,
        )
    #print(frac_ai)

    # Calculates fusion power density, [MW/m^3]
    qfus = n_fuel2*SigV*Efus_DT

    # Printing values
    dout['powers']['fus'] = {}
    dout['powers']['alpi'] = {}
    dout['powers']['alpe'] = {}
    if verb == 2:
        Pfus_tot = cumtrapz(qfus, dout['vol_m3'], initial=0.) # [MW], fusion power
        print('Pfus= {:1.2f} MW\n'.format(Pfus_tot[-1]))
        dout['powers']['fus']['tot_MW'] = Pfus_tot[-1]

    # Return total alpha heating, to the ion, and to the electrons [MW/m^3], dim(rhot,)
    dout['powers']['fus']['prof_MW/m3'] = qfus
    dout['powers']['alpi']['prof_MW/m3'] = frac_ai*qfus/5
    dout['powers']['alpe']['prof_MW/m3'] = (1-frac_ai)*qfus/5
    return dout