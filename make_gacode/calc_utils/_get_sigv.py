'''

_get_sigv.py is a module storing various functions helpful
when calcualting fusion power

'''

import numpy as np
import warnings
import scipy.constants as cnt

def _get_sigmav(reaction, Ti_keV):
    '''
    get_sigmav is a function that calculates the fusion reactivity <\sigma*v> 
    as a function of electron temperature for the following reactions:
        D(T,n)4He, D(D,n)3He, and D(D,p)T, Hale/Bosch 1992
    INPUTS:     reaction -- [char] Fusion reaction of interest [options: DT,
                                   DD_n3He, DD_pT]
                T_i_kev -- [array] Ion temperature in units keV
    OUTPUTS:    sigmav -- [array] Fusion reactivity in units m^3/s
    
    '''
    
    # If interested in D(T,n)4He reaction
    if reaction == "DT":
        c0 = 6.6610 # keV^1/3
        c1 = 643.41*1e-16 # cm^3/s
        c2 = 15.136*1e-3 # keV^-1
        c3 = 75.189*1e-3 # keV^-1
        c4 = 4.6064*1e-3 # keV^-2
        c5 = 13.500*1e-3 # keV^-2
        c6 = -0.10675*1e-3 # keV^-3
        c7 = 0.01366*1e-3 # keV^-3
        
    # If interested in D(D,n)3He reaction 
    elif reaction == "DD_n3He":
        c0 = 6.2696 # keV^1/3
        c1 = 3.5741*1e-16 # cm^3/s
        c2 = 5.8577*1e-3 # keV^-1
        c3 = 7.6822*1e-3 # keV^-1
        c4 = 0*1e-3 # keV^-2
        c5 = -0.002964*1e-3 # keV^-2
        c6 = 0*1e-3 # keV^-3
        c7 = 0*1e-3 # keV^-3
      
    # If interested in D(D,p)T reaction    
    elif reaction == "DD_pT":
        c0 = 6.2696 # keV^1/3
        c1 = 3.7212*1e-16 # cm^3/s
        c2 = 5.4127*1e-3 # keV^-1
        c3 = 1.9917*1e-3 # keV^-1
        c4 = 0*1e-3 # keV^-2
        c5 = 0.010506*1e-3 # keV^-2
        c6 = 0*1e-3 # keV^-3
        c7 = 0*1e-3 # keV^-3

    # Calculates fit forms for <\sigma*v> parameterization
    zeta = 1 - (c2*Ti_keV + c4*Ti_keV**2 + c6*Ti_keV**3) / \
        (1 + c3*Ti_keV + c5*Ti_keV**2 + c7*Ti_keV**3)
    xsi = c0/Ti_keV**(1/3)

    # Calculates the fusion reactivity, [m^3/s]
    return c1/1e6 * zeta**(-5/6) * xsi**2 * np.exp(-3*zeta**(1/3)*xsi)

def _get_ai_frac(
    dout=None,
    fast_ion='He4',
    E_fast=17.59*1.6021e-19/5, # [MJ] alpha energy from DT fusion
    run=True,
    ):
    #n_fuel, ne, n_ash, n_imp, z_imp, mimp_amu, Te, nexp):
    '''
    get_ai_frac is a function to compute a low-accuracy but fast approximation to the ion-alpha
    heating fraction
                 x
              1  /     dy
     F(x) = --- | -----------
             x  /  1+y^(3/2)
                0
    Here, F is the fraction of the alpha energy transferred to ions (at common temperature Ti)
    by collisions, and 

        x = E_alpha/E_crit
    Details are given in Stix, Plasma Phys. 14 (1972) 367, eqn 17. The function F is derived from
    Sivukhin's energy loss equation

    (OLD) INPUTS:  n_fuel  -- [array], [1e19 1/m^3] fuel ion density (assumes 50/50 DT)
             ne      -- [array], [1e19 1/m^3] electron density
             n_ash   -- [array], [1e19 1/m^3] ash ion density
             n_imp   -- [array], [1e19 1/m^3] impurity ion density
             z_imp   -- [scalar], [] impurity charge
             mimp_amu-- [scalar], [amu] impurity mass 
             Te      -- [array], [keV] electron temperature
             nexp    -- [scalar], [] number of grid points

    OUTPUTS: frac_ai -- [array], [] alpha-ion heating fraction
    '''

    # Initializes array to store ai-fraction
    frac_ai = np.zeros(len(dout['rhot']))

    # This only really applies to DT
    #if 'T' in dout['ions'].keys() and 'D' in dout['ions'].keys():
    if run:

        # Calculates alpha heating coefficients [Stix, Plasma Phys. 14 (1972) 367], particularly Eqns 15, 17
        c_a = np.zeros(len(dout['rhot']))
        for ion in dout['ions'].keys():
            if ion in ['nion', 'name', 'charge', 'mass', 'itype']:
                continue
            c_a += (
                (
                    dout['ions'][ion]['ni_tot_19m3']
                    /dout['ne_19m3']
                )
                *dout['ions'][ion]['Z']**2
                /
                (
                    dout['ions'][ion]['M']
                    /dout['ions'][fast_ion]['M']
                )
            )

        # Calculates critical energy, [MJ]
        E_crit = (
            (dout['Te_keV']*1e3*cnt.e/1e6)
            *(4*np.sqrt(float(dout['e']['mass'])/dout['ions'][fast_ion]['M'])
            /(3*np.sqrt(np.pi)*c_a)
            )**(-2/3)
        )

        # Calculates fraction to alpha fusion energy
        x = E_fast/E_crit

        # Loop over radial grid points
        for i in np.arange(len(frac_ai)):
            # Performs heating fraction integral
            if x[i] > 0.1:
                # Large-x asymptotic formula
                if x[i] > 4.0:
                    f = (2*np.pi/3)/np.sin(2*np.pi/3)-2.0/np.sqrt(x[i])+0.5/(x[i]*x[i])

                # Numerical integration
                else:
                    # Number of integration points
                    n = 12
                
                    # Energy grid length
                    dy = x[i]/(n-1)

                    # Intializes scalar to store integration
                    f = 0

                    # Loop over energy grid points
                    for j in np.arange(n):
                        yi = j*dy
                        # Boundary points
                        if j == 0 or j == n-1:
                            f += 0.5/(1.0+yi**1.5)
                        # Interior points    
                        else:
                            f += 1.0/(1.0+yi**1.5)
                
                    # Integration value
                    f = f*dy

                # Alpha-ion heating fraction
                frac_ai[i] = f/x[i]

            # Small-x asymptotic series
            else:
                frac_ai[i] = 1.0-0.4*x[i]**1.5

    # Return alpha-ion heating fraction
    return frac_ai
